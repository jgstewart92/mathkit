/* mathkit: decimation of a closed (or open) triangle mesh without new folds, crossings or drift. Header only, C23.
   mk_decimate uses quadric error costs (Garland and Heckbert 1997) with half-edge collapses only, so every kept
   vertex is one of the input's at its input coordinates (a caller can write the input's own text back, digit for
   digit), and the link condition, so a closed 2-manifold keeps its topology. Each collapse is refused when it would
     turn   a triangle more than acos(turn) from where it was or from the input's normals there;
     edge   make an edge longer than maxedge;
     angle  make a triangle's smallest angle under minangle, unless the fan already had a smaller one;
     fold   lengthen, among the edges it touches, the summed length of edges whose two faces differ by over 90 degrees,
            or by over `fold` degrees (per region when a region function is given), so every fold length can only fall;
     dist   leave an input vertex farther than dist from the result (each removed vertex is owned by a triangle of the
            fan that replaced it and re-measured at each collapse: exact, one side), or put a sample point of a new
            triangle farther than dist from the input surface (the input's faces in the cells around it: the other side);
     cross  make a new triangle cross a live one it shares no vertex with (an edge through the other's interior,
            mk_seg_tri_strict, for the pairs whose boxes meet, as mk_bvh_crossings counts them), so an input with no
            crossing pairs gives none;
     low    make a triangle under minalt high (twice its area over its longest edge), unless the fan already had one as
            low, so an input with none gives none.
   Each vertex queues its cheapest target that passes the cheap refusals (turn, edge, angle, fold, low), one entry a
   vertex, of equal costs the lower-numbered vertex first; distance and crossing are tried when it comes up, and on a
   refusal the next-best target is queued. Stops at `target` vertices or when nothing more is allowed.
   keep (optional, a byte per vertex): 0 free; 1 or 3 never removed; 2 a seam vertex, removed only along the seam (onto
   a 2 or 3, the held side's triangles staying in their planes). An open surface's rim vertex goes only along the rim.
   clean (1): a triangle caught in a crossing (before decimation) or under minalt high (after) is removed by collapsing
   one of its edges under the same refusals; one that still stands is tried once more under the link, keep, seam and
   crossing refusals alone (a needle only by an edge under 20 minalt).
   Single-threaded.
   Size: the marks are 32-bit stamps, about 48 per collapse on an icosphere; past 2^32 of them (about 9e7 collapses, an
   input near 1e8 vertices) a stale mark can read as current. */
#ifndef MK_DEC_H
#define MK_DEC_H
#include "mk_geom.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MK_DEC_REGIONS 8
#define MK_DEC_MAXN 256 /* a vertex with this many triangles or more is never removed */
typedef struct {
    double dist;     /* the two-sided distance bound; 0: none */
    double maxedge;  /* the longest new edge; 0: none */
    double minangle; /* degrees, the smallest new angle; 0: none */
    double fold;     /* degrees: the second fold threshold (over 90 always counts); 0: only 90 */
    double turn;     /* a triangle's normal keeps a cosine of at least this with its old one and the input's */
    double minalt;   /* the lowest new triangle and the clean-up's needle height; 0: none (no low refusal, no needles removed) */
    int target;      /* stop at this many vertices; 0: when nothing more is allowed */
    int clean;       /* 1: the crossing and needle clean-up */
    const unsigned char *keep;                 /* NULL, or a byte per vertex (0, 1, 2, 3: see the top) */
    int nregion;                               /* 1, or up to MK_DEC_REGIONS with region */
    int (*region)(const double *p, void *ctx); /* the region of a point (an edge's midpoint, a triangle's centroid),
                                                  0 .. nregion - 1; any other value counts as 0 */
    void *region_ctx;
} mk_dec_opts;
typedef struct {
    int nv, nt;      /* the result's vertex and triangle counts */
    double maxown;   /* the largest distance of a removed input vertex from the triangle that owns it */
    /* fold lengths by region, over 90 and over `fold` degrees, before and after; areas by region, before and after */
    double fold0[MK_DEC_REGIONS][2], fold1[MK_DEC_REGIONS][2];
    double area0[MK_DEC_REGIONS], area1[MK_DEC_REGIONS];
    long refused[9]; /* refusals by: link, turn, edge, angle, fold, owned distance, sample distance, crossing, low */
    /* clean-up: triangles removed and left before decimation, removed and left after; removals by the relaxed try */
    long clean_before, left_before, clean_after, left_after, relaxed;
    long rim0, rim1; /* an open surface's rim vertices before and after */
} mk_dec_report;

/* defaults for a mesh in which one millimetre is mm units: dist 0.3 mm, maxedge 4 mm, minangle 5, fold 80, turn 0.2,
   minalt 0.001 mm, clean on, one region */
static inline void mk_dec_defaults(mk_dec_opts *o, double mm)
{
    memset(o, 0, sizeof *o);
    o->dist = 0.3 * mm, o->maxedge = 4 * mm, o->minangle = 5, o->fold = 80, o->turn = 0.2, o->minalt = 0.001 * mm, o->clean = 1;
    o->nregion = 1;
}

typedef struct { double x, y, z; } mk_dv3_;
static inline mk_dv3_ mk_dv_(double x, double y, double z) { return (mk_dv3_){x, y, z}; }
static inline mk_dv3_ mk_dsub_(mk_dv3_ a, mk_dv3_ b) { return mk_dv_(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline mk_dv3_ mk_dadd_(mk_dv3_ a, mk_dv3_ b) { return mk_dv_(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline mk_dv3_ mk_dmul_(mk_dv3_ a, double s) { return mk_dv_(a.x * s, a.y * s, a.z * s); }
static inline double mk_ddot_(mk_dv3_ a, mk_dv3_ b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline mk_dv3_ mk_dcross_(mk_dv3_ a, mk_dv3_ b) { return mk_dv_(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
static inline double mk_dlen_(mk_dv3_ a) { return sqrt(mk_ddot_(a, a)); }
static inline mk_dv3_ mk_dnorm_(mk_dv3_ a) { const double l = mk_dlen_(a); return l > 0 ? mk_dmul_(a, 1 / l) : a; }
static inline mk_dv3_ mk_dtri_n_(mk_dv3_ a, mk_dv3_ b, mk_dv3_ c) { return mk_dcross_(mk_dsub_(b, a), mk_dsub_(c, a)); }
/* the larger and the smaller of two numbers inline (fmin and fmax are C library calls): mk_dec_max_ as fmax, NaN
   included; mk_dec_hi_ and mk_dec_lo_ for numbers never NaN */
static inline double mk_dec_max_(double a, double b) { return a > b || b != b ? a : b; }
static inline double mk_dec_hi_(double a, double b) { return b > a ? b : a; }
static inline double mk_dec_lo_(double a, double b) { return b < a ? b : a; }
static inline double mk_dmax_cos_(mk_dv3_ a, mk_dv3_ b, mk_dv3_ c)
{
    const mk_dv3_ e0 = mk_dsub_(b, a), e1 = mk_dsub_(c, b), e2 = mk_dsub_(a, c);
    const double l0 = mk_dlen_(e0), l1 = mk_dlen_(e1), l2 = mk_dlen_(e2);
    if (!(l0 > 0) || !(l1 > 0) || !(l2 > 0)) return 1;
    const double ca = -mk_ddot_(e0, e2) / (l0 * l2), cb = -mk_ddot_(e0, e1) / (l0 * l1), cc = -mk_ddot_(e1, e2) / (l1 * l2);
    return mk_dec_max_(ca, mk_dec_max_(cb, cc));
}
/* an edge pq through triangle abc's interior (mk_seg_tri_strict) */
static inline int mk_dec_segtri_(mk_dv3_ p, mk_dv3_ q, mk_dv3_ a, mk_dv3_ b, mk_dv3_ c)
{
    const double P[3] = {p.x, p.y, p.z}, Q[3] = {q.x, q.y, q.z}, A[3] = {a.x, a.y, a.z}, B[3] = {b.x, b.y, b.z}, C[3] = {c.x, c.y, c.z};
    return mk_seg_tri_strict(P, Q, A, B, C);
}

typedef struct { double c; int u; } mk_dec_ent_;
typedef struct {
    const double *P;
    int nv, nt, nreg, relax, oom, gh;
    int *F, *VA, *NVT, *CVT;
    size_t *VO, vn, vcap;
    unsigned *MARK, MARKS; /* unsigned, so the marks wrap rather than overflow */
    double *N0, *Q;
    unsigned char *deadv, *deadt, *bnd;
    const unsigned char *keep;
    const mk_dec_opts *o;
    double dist2, maxedge2, cos2, cmin, mindot, gc, gs, samp, ftol, relaxlen, minalt, maxown2, xpad, dnear;
    /* the input's own faces T0 in the same cells (side gs): bucket k's at ntr[nst[k]] .. ntr[nst[k + 1] - 1], ngh
       buckets */
    const int *T0;
    int *nst, *ntr, ngh;
    long ref[9], nrelax;
    double (*nc)[9];
    int *ntb;
    /* the crossing grid's bucket k: its entries as at the last rebuild, ce[2 i] (a triangle) and ce[2 i + 1] (its
       version) for cst[k] <= i < cst[k + 1], then a list of those entered since from ghead[k], through go[3 e] (the
       next), go[3 e + 1] (the triangle) and go[3 e + 2] (its version), gn of gcap taken */
    int *cst, *ce, *ghead, *go, *tver, gn, gcap;
    double glo[3];
    int *ohead, *onext, *pend_p, *pend_t, pend_n, pend_c, pend_u, pend_v;
    mk_dec_ent_ *heap;
    int *hpos, *htv, *hsk, nh;
    /* the set of triangles one crossing query has met (open addressing, a slot current when its stamp is hmark) */
    int *hkey, hcap, hn;
    unsigned *hst, hmark;
    /* the fan of fu - 1 as it stood at collapse fcol (fu 0: none); level fw 1: each triangle's normal, its length and the
       fan's largest corner cosine; 2: as well the unit normals, each spoke's other triangle, the triangle across each
       outer edge and its unit normal, and the fold lengths among all the fan's edges */
    unsigned ncol, fcol;
    int fu, fw;
    double foldc, flow, fl0[MK_DEC_MAXN], fb[MK_DEC_REGIONS][2];
    mk_dv3_ fn0[MK_DEC_MAXN], fun[MK_DEC_MAXN], fon[MK_DEC_MAXN], fn1[MK_DEC_MAXN];
    int fot[MK_DEC_MAXN], fs[MK_DEC_MAXN][2], fpr[MK_DEC_MAXN][2];
    /* the new triangles of the collapse under test: corners other than v, boxes and the union of those not too large
       for the grid, too large for the grid */
    int xg[MK_DEC_MAXN][2];
    double xlo[MK_DEC_MAXN][3], xhi[MK_DEC_MAXN][3], xulo[3], xuhi[3];
    unsigned char xbig[MK_DEC_MAXN];
} mk_dec_;

static inline mk_dv3_ mk_dec_pt_(const mk_dec_ *c, int i) { return mk_dv_(c->P[3 * (size_t)i], c->P[3 * (size_t)i + 1], c->P[3 * (size_t)i + 2]); }
static inline int *mk_dec_f_(const mk_dec_ *c, int t) { return c->F + 3 * (size_t)t; }
static inline int mk_dec_has_(const mk_dec_ *c, int t, int x) { const int *f = mk_dec_f_(c, t); return f[0] == x || f[1] == x || f[2] == x; }
static inline int mk_dec_has3_(const int *g, int x) { return g[0] == x || g[1] == x || g[2] == x; }
static inline void *mk_dec_grow_(mk_dec_ *c, void *p, size_t n, size_t sz)
{
    if (c->oom) return p;
    void *q = n <= SIZE_MAX / sz ? realloc(p, n * sz) : NULL;
    if (!q) c->oom = 1;
    return q ? q : p;
}
/* each vertex's triangles: NVT[v] of them at VA + VO[v], room for CVT[v] (0: none, the vertex gone), the blocks in
   one array of vcap ints of which vn are taken; a full block moves to the end at twice the room, and a full array is
   compacted to its live blocks and a quarter more */
static inline int *mk_dec_vt_(const mk_dec_ *c, int v) { return c->VA + c->VO[v]; }
static inline void mk_dec_vpack_(mk_dec_ *c, size_t more)
{
    size_t live = 0;
    for (int v = 0; v < c->nv; v++) live += (size_t)c->CVT[v];
    const size_t cap = live + live / 4 + more + 4096;
    int *A = cap <= SIZE_MAX / sizeof *A ? malloc(sizeof *A * cap) : NULL;
    if (!A) { c->oom = 1; return; }
    size_t n = 0;
    for (int v = 0; v < c->nv; v++) {
        if (c->CVT[v] && c->NVT[v]) memcpy(A + n, mk_dec_vt_(c, v), sizeof *A * (size_t)c->NVT[v]);
        c->VO[v] = n, n += (size_t)c->CVT[v];
    }
    free(c->VA);
    c->VA = A, c->vn = n, c->vcap = cap;
}
/* the blocks for the faces F: each vertex's count and two more */
static inline void mk_dec_fans_(mk_dec_ *c)
{
    c->VO = malloc(sizeof *c->VO * (size_t)c->nv);
    if (!c->VO) { c->oom = 1; return; }
    for (int v = 0; v < c->nv; v++) c->NVT[v] = 0;
    for (size_t q = 0; q < 3 * (size_t)c->nt; q++) c->NVT[c->F[q]]++;
    for (int v = 0; v < c->nv; v++) c->CVT[v] = c->NVT[v] + 2, c->NVT[v] = 0;
    mk_dec_vpack_(c, 0);
    for (int t = 0; t < c->nt && !c->oom; t++)
        for (int k = 0; k < 3; k++) {
            const int v = c->F[3 * (size_t)t + k];
            mk_dec_vt_(c, v)[c->NVT[v]++] = t;
        }
}
static inline void mk_dec_fan_add_(mk_dec_ *c, int v, int t)
{
    if (c->NVT[v] == c->CVT[v]) {
        const int nc = c->CVT[v] ? 2 * c->CVT[v] : 8;
        if (c->vn + (size_t)nc > c->vcap) mk_dec_vpack_(c, (size_t)nc);
        if (c->oom) return;
        memmove(c->VA + c->vn, mk_dec_vt_(c, v), sizeof *c->VA * (size_t)c->NVT[v]);
        c->VO[v] = c->vn, c->vn += (size_t)nc, c->CVT[v] = nc;
    }
    mk_dec_vt_(c, v)[c->NVT[v]++] = t;
}
static inline void mk_dec_fan_del_(mk_dec_ *c, int v, int t)
{
    for (int i = 0; i < c->NVT[v]; i++)
        if (mk_dec_vt_(c, v)[i] == t) { mk_dec_vt_(c, v)[i] = mk_dec_vt_(c, v)[--c->NVT[v]]; return; }
}
static inline int mk_dec_nbrs_(mk_dec_ *c, int u, int *nb)
{
    int n = 0;
    c->MARKS += 2; /* past mk_dec_cheap_'s MARKS + 1, so its common vertices do not read as seen */
    for (int i = 0; i < c->NVT[u]; i++)
        for (int k = 0; k < 3; k++) {
            const int x = mk_dec_f_(c, mk_dec_vt_(c, u)[i])[k];
            if (x != u && c->MARK[x] != c->MARKS && n < MK_DEC_MAXN) c->MARK[x] = c->MARKS, nb[n++] = x;
        }
    return n;
}
static inline int mk_dec_region_(const mk_dec_ *c, mk_dv3_ m)
{
    if (c->nreg == 1) return 0;
    const double p[3] = {m.x, m.y, m.z};
    const int r = c->o->region(p, c->o->region_ctx);
    return r >= 0 && r < c->nreg ? r : 0;
}
static inline mk_dv3_ mk_dec_unit_(const mk_dec_ *c, int t)
{
    const int *g = mk_dec_f_(c, t);
    return mk_dnorm_(mk_dtri_n_(mk_dec_pt_(c, g[0]), mk_dec_pt_(c, g[1]), mk_dec_pt_(c, g[2])));
}
/* edge ab between faces of unit normals n1, n2: its length into out[r][0] when they differ by over 90 degrees, into
   out[r][1] when by over fold */
static inline void mk_dec_fedge_(const mk_dec_ *c, int a, int b, mk_dv3_ n1, mk_dv3_ n2, double out[MK_DEC_REGIONS][2])
{
    const double d = mk_ddot_(n1, n2);
    if (d >= c->cos2 && d >= 0) return;
    const mk_dv3_ pa = mk_dec_pt_(c, a), pb = mk_dec_pt_(c, b);
    const double L = mk_dlen_(mk_dsub_(pa, pb));
    const int r = mk_dec_region_(c, mk_dmul_(mk_dadd_(pa, pb), 0.5));
    if (d < 0) out[r][0] += L;
    if (d < c->cos2) out[r][1] += L;
}
/* the triangle other than those holding x that holds a and b, or -1 */
static inline int mk_dec_across_(const mk_dec_ *c, int a, int b, int x)
{
    for (int q = 0; q < c->NVT[a]; q++) {
        const int t = mk_dec_vt_(c, a)[q];
        if (mk_dec_has_(c, t, b) && !mk_dec_has_(c, t, x)) return t;
    }
    return -1;
}
/* a triangle's longest edge, squared */
static inline double mk_dec_e2_(mk_dv3_ a, mk_dv3_ b, mk_dv3_ c)
{
    const mk_dv3_ e0 = mk_dsub_(b, a), e1 = mk_dsub_(c, b), e2 = mk_dsub_(a, c);
    return mk_dec_hi_(mk_ddot_(e0, e0), mk_dec_hi_(mk_ddot_(e1, e1), mk_ddot_(e2, e2)));
}
/* the fan cache of u to level w (see mk_dec_) */
static inline void mk_dec_fan_(mk_dec_ *c, int u, int w)
{
    if (c->fu != u + 1 || c->fcol != c->ncol) {
        double oldc = 0, low = INFINITY;
        for (int i = 0; i < c->NVT[u]; i++) {
            const int *f = mk_dec_f_(c, mk_dec_vt_(c, u)[i]);
            const mk_dv3_ p0 = mk_dec_pt_(c, f[0]), p1 = mk_dec_pt_(c, f[1]), p2 = mk_dec_pt_(c, f[2]);
            oldc = mk_dec_max_(oldc, mk_dmax_cos_(p0, p1, p2));
            c->fn0[i] = mk_dtri_n_(p0, p1, p2), c->fl0[i] = mk_dlen_(c->fn0[i]);
            const double e0 = mk_dec_e2_(p0, p1, p2), h2 = e0 > 0 ? c->fl0[i] * c->fl0[i] / e0 : 0;
            low = h2 < low ? h2 : low;
        }
        c->foldc = oldc, c->flow = low, c->fu = u + 1, c->fcol = c->ncol, c->fw = 1;
    }
    if (w < 2 || c->fw >= 2) return;
    const int n = c->NVT[u];
    for (int r = 0; r < c->nreg; r++) c->fb[r][0] = c->fb[r][1] = 0;
    for (int i = 0; i < n; i++) {
        const int *f = mk_dec_f_(c, mk_dec_vt_(c, u)[i]);
        int k = 0;
        for (int q = 0; q < 3; q++)
            if (f[q] != u && k < 2) c->fs[i][k++] = f[q];
        c->fun[i] = mk_dnorm_(c->fn0[i]);
        c->fot[i] = k == 2 ? mk_dec_across_(c, c->fs[i][0], c->fs[i][1], u) : -1;
        if (c->fot[i] >= 0) c->fon[i] = mk_dec_unit_(c, c->fot[i]);
    }
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 2; k++) {
            c->fpr[i][k] = -1;
            for (int j = 0; j < n && c->fpr[i][k] < 0; j++)
                if (j != i && mk_dec_has_(c, mk_dec_vt_(c, u)[j], c->fs[i][k])) c->fpr[i][k] = j;
        }
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < 2; k++)
            if (c->fpr[i][k] > i) mk_dec_fedge_(c, u, c->fs[i][k], c->fun[i], c->fun[c->fpr[i][k]], c->fb);
        if (c->fot[i] >= 0) mk_dec_fedge_(c, c->fs[i][0], c->fs[i][1], c->fun[i], c->fon[i], c->fb);
    }
    c->fw = 2;
}
/* the fold lengths among the edges of the triangles of fan(u) moved onto v (fn1: their normals), against those among
   the edges of fan(u) as it stands (mk_dec_fan_ level 2): every other edge keeps its faces. 1 when a length rises */
static inline int mk_dec_fold_up_(mk_dec_ *c, int u, int v)
{
    double a[MK_DEC_REGIONS][2];
    mk_dv3_ nm[MK_DEC_MAXN];
    const int n = c->NVT[u];
    for (int r = 0; r < c->nreg; r++) a[r][0] = a[r][1] = 0;
    for (int i = 0; i < n; i++)
        if (!mk_dec_has_(c, mk_dec_vt_(c, u)[i], v)) nm[i] = mk_dnorm_(c->fn1[i]);
    for (int i = 0; i < n; i++) {
        if (mk_dec_has_(c, mk_dec_vt_(c, u)[i], v)) continue;
        for (int k = 0; k < 2; k++) {
            const int j = c->fpr[i][k], x = c->fs[i][k];
            if (j < 0) continue;
            if (!mk_dec_has_(c, mk_dec_vt_(c, u)[j], v)) {
                if (j > i) mk_dec_fedge_(c, v, x, nm[i], nm[j], a);
                continue;
            }
            const int t3 = mk_dec_across_(c, v, x, u); /* j dies: after the collapse a triangle of v's lies across vx */
            if (t3 >= 0) mk_dec_fedge_(c, v, x, nm[i], mk_dec_unit_(c, t3), a);
        }
        if (c->fot[i] >= 0) mk_dec_fedge_(c, c->fs[i][0], c->fs[i][1], nm[i], c->fon[i], a);
    }
    for (int r = 0; r < c->nreg; r++)
        if (a[r][0] > c->fb[r][0] + c->ftol || a[r][1] > c->fb[r][1] + c->ftol) return 1;
    return 0;
}
/* the crossing grid: an insert-only hash of cells (side gs) whose entries carry the triangle's version. A triangle
   whose box covers more than MK_DEC_BIG cells is listed once, in cell (-1, -1, -1), which every query reads; such a new
   triangle reads every live one instead of its cells, so one large triangle among small ones does not fill memory with
   its box's cells. */
#define MK_DEC_BIG 4096
static inline int mk_dec_cells_(const mk_dec_ *c, const mk_dv3_ *d, int lo[3], int hi[3])
{
    const double mn[3] = {mk_dec_lo_(d[0].x, mk_dec_lo_(d[1].x, d[2].x)), mk_dec_lo_(d[0].y, mk_dec_lo_(d[1].y, d[2].y)), mk_dec_lo_(d[0].z, mk_dec_lo_(d[1].z, d[2].z))};
    const double mx[3] = {mk_dec_hi_(d[0].x, mk_dec_hi_(d[1].x, d[2].x)), mk_dec_hi_(d[0].y, mk_dec_hi_(d[1].y, d[2].y)), mk_dec_hi_(d[0].z, mk_dec_hi_(d[1].z, d[2].z))};
    for (int k = 0; k < 3; k++) lo[k] = (int)floor((mn[k] - c->glo[k]) / c->gs), hi[k] = (int)floor((mx[k] - c->glo[k]) / c->gs);
#if defined(MK_MUTATE) && MK_MUTATE == 227
    return 0; /* 227: every triangle entered in every cell of its box */
#endif
    return (double)(hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1) > MK_DEC_BIG;
}
/* cell (x, y, z)'s bucket of n, a power of 2 */
static inline size_t mk_dec_hkey_(int x, int y, int z, int n)
{
    return (size_t)((((uint64_t)(uint32_t)x * 73856093u) ^ ((uint64_t)(uint32_t)y * 19349663u) ^ ((uint64_t)(uint32_t)z * 83492791u)) & (uint64_t)(n - 1));
}
static inline size_t mk_dec_gkey_(const mk_dec_ *c, int x, int y, int z) { return mk_dec_hkey_(x, y, z, c->gh); }
/* t into the list of each bucket of its box's cells (one bucket when it is too large for the grid) */
static inline void mk_dec_gins_(mk_dec_ *c, int t)
{
    const int *f = mk_dec_f_(c, t);
    const mk_dv3_ d[3] = {mk_dec_pt_(c, f[0]), mk_dec_pt_(c, f[1]), mk_dec_pt_(c, f[2])};
    int lo[3], hi[3];
    if (mk_dec_cells_(c, d, lo, hi)) lo[0] = lo[1] = lo[2] = hi[0] = hi[1] = hi[2] = -1;
    for (int x = lo[0]; x <= hi[0]; x++)
        for (int y = lo[1]; y <= hi[1]; y++)
            for (int z = lo[2]; z <= hi[2]; z++) {
                if (c->gn == c->gcap) {
                    if (c->gcap > INT_MAX / 6) { c->oom = 1; return; }
                    const int nc = c->gcap ? 2 * c->gcap : 1 << 16;
                    c->go = mk_dec_grow_(c, c->go, 3 * (size_t)nc, sizeof *c->go);
                    if (c->oom) return;
                    c->gcap = nc;
                }
                const size_t k = mk_dec_gkey_(c, x, y, z);
                int *g = c->go + 3 * (size_t)c->gn;
                g[0] = c->ghead[k], g[1] = t, g[2] = c->tver[t], c->ghead[k] = c->gn++;
            }
}
/* the buckets of t's cells: pass 0 counts each into cst[k + 1], pass 1 enters t at cst[k] and moves it on */
static inline void mk_dec_gcsr_(mk_dec_ *c, int t, int pass)
{
    const int *f = mk_dec_f_(c, t);
    const mk_dv3_ d[3] = {mk_dec_pt_(c, f[0]), mk_dec_pt_(c, f[1]), mk_dec_pt_(c, f[2])};
    int lo[3], hi[3];
    if (mk_dec_cells_(c, d, lo, hi)) lo[0] = lo[1] = lo[2] = hi[0] = hi[1] = hi[2] = -1;
    for (int x = lo[0]; x <= hi[0]; x++)
        for (int y = lo[1]; y <= hi[1]; y++)
            for (int z = lo[2]; z <= hi[2]; z++) {
                const size_t k = mk_dec_gkey_(c, x, y, z);
                if (!pass) c->cst[k + 1]++;
                else {
                    const size_t i = (size_t)c->cst[k]++;
                    c->ce[2 * i] = t, c->ce[2 * i + 1] = c->tver[t];
                }
            }
}
/* the grid anew from the live triangles, each bucket's entries side by side; room for half as many again entered
   later (past that it is rebuilt) */
static inline void mk_dec_grebuild_(mk_dec_ *c)
{
    if (!c->gh) {
        int gh = 1024;
        while (gh < c->nt / 2 && gh < (1 << 26)) gh *= 2;
        c->cst = malloc(sizeof *c->cst * ((size_t)gh + 1)), c->ghead = malloc(sizeof *c->ghead * (size_t)gh), c->gh = gh;
        if (!c->cst || !c->ghead) { c->oom = 1; return; }
    }
    for (int k = 0; k <= c->gh; k++) c->cst[k] = 0;
    for (int t = 0; t < c->nt; t++)
        if (!c->deadt[t]) mk_dec_gcsr_(c, t, 0);
    size_t need = 0;
    for (int k = 1; k <= c->gh; k++) {
        need += (size_t)c->cst[k];
        if (need > INT_MAX / 2) { c->oom = 1; return; }
        c->cst[k] = (int)need;
    }
    const size_t later = need / 2 + 65536;
    free(c->ce), c->ce = malloc(sizeof *c->ce * 2 * (need ? need : 1));
    if (later > (size_t)c->gcap || later < (size_t)c->gcap / 2) free(c->go), c->go = malloc(sizeof *c->go * 3 * later), c->gcap = (int)later;
    if (!c->ce || !c->go) { c->oom = 1, c->gcap = c->gn = 0; return; }
    for (int t = 0; t < c->nt; t++)
        if (!c->deadt[t]) mk_dec_gcsr_(c, t, 1);
    for (int k = c->gh; k > 0; k--) c->cst[k] = c->cst[k - 1];
    c->cst[0] = 0;
    for (int k = 0; k < c->gh; k++) c->ghead[k] = -1;
    c->gn = 0;
}
/* the next current triangle of bucket k, or -1: *i walks the rebuilt entries (start it at cst[k]), *e the list (start
   it at ghead[k]) */
static inline int mk_dec_gwalk_(const mk_dec_ *c, size_t k, int *i, int *e)
{
    for (; *i < c->cst[k + 1]; ++*i) {
        const int t = c->ce[2 * (size_t)*i];
        if (c->ce[2 * (size_t)*i + 1] == c->tver[t]) { ++*i; return t; }
    }
    while (*e >= 0) {
        const int *g = c->go + 3 * (size_t)*e;
        *e = g[0];
        if (g[2] == c->tver[g[1]]) return g[1];
    }
    return -1;
}
/* a new query of the met set */
static inline void mk_dec_hclear_(mk_dec_ *c)
{
    c->hn = 0;
    if (!++c->hmark) {
        memset(c->hst, 0, sizeof *c->hst * (size_t)c->hcap);
        c->hmark = 1;
    }
}
/* 1 if the query has met t already, else 0 (and t is marked met); out of memory counts as met */
static inline int mk_dec_hmet_(mk_dec_ *c, int t)
{
    if (2 * (c->hn + 1) > c->hcap) {
        const int nc = 2 * c->hcap;
        int *k = malloc(sizeof *k * (size_t)nc);
        unsigned *s = calloc((size_t)nc, sizeof *s);
        if (!k || !s) { free(k), free(s), c->oom = 1; return 1; }
        for (int i = 0; i < c->hcap; i++)
            if (c->hst[i] == c->hmark) {
                size_t j = (uint32_t)c->hkey[i] * 2654435761u & (uint32_t)(nc - 1);
                while (s[j]) j = (j + 1) & (size_t)(nc - 1);
                s[j] = 1, k[j] = c->hkey[i];
            }
        free(c->hkey), free(c->hst);
        c->hkey = k, c->hst = s, c->hcap = nc, c->hmark = 1;
    }
    size_t i = (uint32_t)t * 2654435761u & (uint32_t)(c->hcap - 1);
    while (c->hst[i] == c->hmark) {
        if (c->hkey[i] == t) return 1;
        i = (i + 1) & (size_t)(c->hcap - 1);
    }
    c->hst[i] = c->hmark, c->hkey[i] = t, c->hn++;
    return 0;
}
/* the box of triangle t, widened by xpad */
static inline void mk_dec_tbox_(const mk_dec_ *c, const int *f, double lo[3], double hi[3])
{
    for (int k = 0; k < 3; k++) {
        const double a = c->P[3 * (size_t)f[0] + k], b = c->P[3 * (size_t)f[1] + k], d = c->P[3 * (size_t)f[2] + k];
        lo[k] = mk_dec_lo_(a, mk_dec_lo_(b, d)) - c->xpad, hi[k] = mk_dec_hi_(a, mk_dec_hi_(b, d)) + c->xpad;
    }
}
static inline int mk_dec_meet_(const double *alo, const double *ahi, const double *blo, const double *bhi)
{
    return alo[0] <= bhi[0] && blo[0] <= ahi[0] && alo[1] <= bhi[1] && blo[1] <= ahi[1] && alo[2] <= bhi[2] && blo[2] <= ahi[2];
}
/* triangles d and f (corners) cross: an edge of one through the other's interior */
static inline int mk_dec_x6_(const mk_dec_ *c, const mk_dv3_ *d, const int *f)
{
#if defined(MK_MUTATE) && MK_MUTATE == 163
    (void)c, (void)d, (void)f;
    return 0; /* 163: the crossing test left out */
#else
    const mk_dv3_ tt[3] = {mk_dec_pt_(c, f[0]), mk_dec_pt_(c, f[1]), mk_dec_pt_(c, f[2])};
    for (int k = 0; k < 3; k++)
        if (mk_dec_segtri_(d[k], d[(k + 1) % 3], tt[0], tt[1], tt[2]) || mk_dec_segtri_(tt[k], tt[(k + 1) % 3], d[0], d[1], d[2])) return 1;
    return 0;
#endif
}
/* new triangle d (its corners g, its box lo hi) against live triangle t sharing none of its corners (nor skip): 1 when
   they cross. met: t is a grid candidate, tried once a query */
/* the box of the triangle of corners f, widened by xpad, meets box lo hi: mk_dec_tbox_ and mk_dec_meet_, an axis at a
   time */
static inline int mk_dec_tmeet_(const mk_dec_ *c, const int *f, const double *lo, const double *hi)
{
    const double *a = c->P + 3 * (size_t)f[0], *b = c->P + 3 * (size_t)f[1], *d = c->P + 3 * (size_t)f[2];
    for (int k = 0; k < 3; k++)
        if (mk_dec_hi_(a[k], mk_dec_hi_(b[k], d[k])) + c->xpad < lo[k] || mk_dec_lo_(a[k], mk_dec_lo_(b[k], d[k])) - c->xpad > hi[k]) return 0;
    return 1;
}
static inline int mk_dec_xone_(mk_dec_ *c, const mk_dv3_ *d, const int *g, int skip, const double *lo, const double *hi, int t, int met)
{
    if (c->deadt[t] || !mk_dec_tmeet_(c, mk_dec_f_(c, t), lo, hi) || (met && mk_dec_hmet_(c, t))) return 0;
    if (mk_dec_has_(c, t, skip) || mk_dec_has_(c, t, g[0]) || mk_dec_has_(c, t, g[1]) || mk_dec_has_(c, t, g[2])) return 0;
    return mk_dec_x6_(c, d, mk_dec_f_(c, t));
}
/* new triangle d (its corners g) crosses a live triangle sharing none of its corners (nor skip) */
static inline int mk_dec_crosses_(mk_dec_ *c, const mk_dv3_ *d, const int *g, int skip)
{
    int lo[3], hi[3];
    const double bl[3] = {mk_dec_lo_(d[0].x, mk_dec_lo_(d[1].x, d[2].x)), mk_dec_lo_(d[0].y, mk_dec_lo_(d[1].y, d[2].y)), mk_dec_lo_(d[0].z, mk_dec_lo_(d[1].z, d[2].z))};
    const double bh[3] = {mk_dec_hi_(d[0].x, mk_dec_hi_(d[1].x, d[2].x)), mk_dec_hi_(d[0].y, mk_dec_hi_(d[1].y, d[2].y)), mk_dec_hi_(d[0].z, mk_dec_hi_(d[1].z, d[2].z))};
    if (mk_dec_cells_(c, d, lo, hi)) {
        for (int t = 0; t < c->nt; t++)
            if (mk_dec_xone_(c, d, g, skip, bl, bh, t, 0)) return 1;
        return 0;
    }
    mk_dec_hclear_(c);
    for (int x = lo[0] - 1; x <= hi[0]; x++) /* x = lo[0] - 1: the bucket of the large triangles */
        for (int y = x < lo[0] ? -1 : lo[1]; y <= (x < lo[0] ? -1 : hi[1]); y++)
            for (int z = x < lo[0] ? -1 : lo[2]; z <= (x < lo[0] ? -1 : hi[2]); z++) {
                const size_t k = x < lo[0] ? mk_dec_gkey_(c, -1, -1, -1) : mk_dec_gkey_(c, x, y, z);
                for (int i = c->cst[k], e = c->ghead[k], t; (t = mk_dec_gwalk_(c, k, &i, &e)) >= 0;)
                    if (mk_dec_xone_(c, d, g, skip, bl, bh, t, 1)) return 1;
            }
    return 0;
}
/* candidate t (live, current) against the new triangles j < m of the collapse u -> v that are not too large for the
   grid: 1 when one crosses */
static inline int mk_dec_xcand_(mk_dec_ *c, int m, int u, int v, int t)
{
    if (c->deadt[t]) return 0;
    const int *f = mk_dec_f_(c, t);
    if (!mk_dec_tmeet_(c, f, c->xulo, c->xuhi) || mk_dec_has3_(f, u) || mk_dec_has3_(f, v) || mk_dec_hmet_(c, t)) return 0;
    double tl[3], th[3];
    mk_dec_tbox_(c, f, tl, th);
    for (int j = 0; j < m; j++) {
        if (c->xbig[j] || mk_dec_has3_(f, c->xg[j][0]) || mk_dec_has3_(f, c->xg[j][1]) || !mk_dec_meet_(tl, th, c->xlo[j], c->xhi[j])) continue;
        const mk_dv3_ d[3] = {mk_dv_(c->nc[j][0], c->nc[j][1], c->nc[j][2]), mk_dv_(c->nc[j][3], c->nc[j][4], c->nc[j][5]), mk_dv_(c->nc[j][6], c->nc[j][7], c->nc[j][8])};
        if (mk_dec_x6_(c, d, f)) return 1;
    }
    return 0;
}
/* the new triangles of the collapse u -> v (their corners in nc, the old triangles in ntb, m of them) against every
   live triangle sharing none of a new one's corners nor u, each candidate met once: 1 when a pair crosses */
static inline int mk_dec_xfan_(mk_dec_ *c, int u, int v, int m)
{
    int lo[3] = {INT_MAX, INT_MAX, INT_MAX}, hi[3] = {INT_MIN, INT_MIN, INT_MIN}, small = 0;
    for (int k = 0; k < 3; k++) c->xulo[k] = INFINITY, c->xuhi[k] = -INFINITY;
    for (int j = 0; j < m; j++) {
        const int *f = mk_dec_f_(c, c->ntb[j]);
        int g[3], q = 0;
        for (int k = 0; k < 3; k++) g[k] = f[k] == u ? v : f[k];
        for (int k = 0; k < 3; k++)
            if (g[k] != v) c->xg[j][q++ & 1] = g[k];
        const mk_dv3_ d[3] = {mk_dv_(c->nc[j][0], c->nc[j][1], c->nc[j][2]), mk_dv_(c->nc[j][3], c->nc[j][4], c->nc[j][5]), mk_dv_(c->nc[j][6], c->nc[j][7], c->nc[j][8])};
        for (int k = 0; k < 3; k++) {
            c->xlo[j][k] = mk_dec_lo_(c->nc[j][k], mk_dec_lo_(c->nc[j][3 + k], c->nc[j][6 + k]));
            c->xhi[j][k] = mk_dec_hi_(c->nc[j][k], mk_dec_hi_(c->nc[j][3 + k], c->nc[j][6 + k]));
        }
        int l[3], h[3];
        c->xbig[j] = (unsigned char)mk_dec_cells_(c, d, l, h);
        if (c->xbig[j]) { /* reads every live triangle */
            if (mk_dec_crosses_(c, d, g, u)) return 1;
            continue;
        }
        small = 1;
        for (int k = 0; k < 3; k++) {
            lo[k] = l[k] < lo[k] ? l[k] : lo[k], hi[k] = h[k] > hi[k] ? h[k] : hi[k];
            c->xulo[k] = mk_dec_lo_(c->xulo[k], c->xlo[j][k]), c->xuhi[k] = mk_dec_hi_(c->xuhi[k], c->xhi[j][k]);
        }
    }
    if (!small) return 0;
    mk_dec_hclear_(c);
    for (int x = lo[0] - 1; x <= hi[0]; x++) /* x = lo[0] - 1: the bucket of the large triangles */
        for (int y = x < lo[0] ? -1 : lo[1]; y <= (x < lo[0] ? -1 : hi[1]); y++)
            for (int z = x < lo[0] ? -1 : lo[2]; z <= (x < lo[0] ? -1 : hi[2]); z++) {
                const size_t k = x < lo[0] ? mk_dec_gkey_(c, -1, -1, -1) : mk_dec_gkey_(c, x, y, z);
                for (int i = c->cst[k], e = c->ghead[k], t; (t = mk_dec_gwalk_(c, k, &i, &e)) >= 0;)
                    if (mk_dec_xcand_(c, m, u, v, t)) return 1;
            }
    return 0;
}
/* the cheap refusals */
static inline int mk_dec_cheap_(mk_dec_ *c, int u, int v)
{
    const unsigned char *K = c->keep;
    (void)K;
    if (c->deadv[u] || c->deadv[v] || u == v || c->NVT[u] >= MK_DEC_MAXN || c->NVT[v] >= MK_DEC_MAXN) return 0;
#if !defined(MK_MUTATE) || MK_MUTATE != 165
    if (K && (K[u] & 1)) return 0;
    /* a seam vertex: onto a seam vertex or a cut point, one held-side triangle dying, the others in their planes */
    if (K && K[u] == 2) {
        if (K[v] < 2) return 0;
        int die = 0;
        for (int i = 0; i < c->NVT[u]; i++) {
            const int t = mk_dec_vt_(c, u)[i], *f = mk_dec_f_(c, t);
            if (!K[f[0]] || !K[f[1]] || !K[f[2]]) continue;
            if (mk_dec_has_(c, t, v)) { die++; continue; }
            mk_dv3_ p0[3], p1[3];
            for (int k = 0; k < 3; k++) p0[k] = mk_dec_pt_(c, f[k]), p1[k] = f[k] == u ? mk_dec_pt_(c, v) : p0[k];
            const mk_dv3_ n0 = mk_dtri_n_(p0[0], p0[1], p0[2]), n1 = mk_dtri_n_(p1[0], p1[1], p1[2]);
            const double l0 = mk_dlen_(n0), l1 = mk_dlen_(n1);
            if (!(l1 > 0) || mk_ddot_(n0, n1) < (1 - 1e-10) * l0 * l1) return 0;
        }
        if (die != 1) return 0;
    }
#endif
    int opp[2] = {-1, -1}, no = 0;
    for (int i = 0; i < c->NVT[u]; i++) {
        const int t = mk_dec_vt_(c, u)[i], *f = mk_dec_f_(c, t);
        if (!mk_dec_has_(c, t, v)) continue;
        if (no < 2) opp[no] = f[0] != u && f[0] != v ? f[0] : f[1] != u && f[1] != v ? f[1] : f[2];
        no++;
    }
    const int need = c->bnd[u] ? 1 : 2; /* a rim vertex goes only along the rim, onto a rim neighbour over a rim edge */
    if ((c->bnd[u] && !c->bnd[v]) || no != need || (need == 2 && opp[0] == opp[1])) return 0;
    c->MARKS += 2;
    for (int i = 0; i < c->NVT[u]; i++)
        for (int k = 0; k < 3; k++) c->MARK[mk_dec_f_(c, mk_dec_vt_(c, u)[i])[k]] = c->MARKS;
    int common = 0;
    for (int i = 0; i < c->NVT[v]; i++)
        for (int k = 0; k < 3; k++) {
            const int x = mk_dec_f_(c, mk_dec_vt_(c, v)[i])[k];
            if (x != u && x != v && c->MARK[x] == c->MARKS) c->MARK[x] = c->MARKS + 1, common++;
        }
#if !defined(MK_MUTATE) || MK_MUTATE != 164
    if (common != need) { c->ref[0]++; return 0; }
    int tu = 0, tv = 0;
    if (need == 2)
        for (int i = 0; i < c->NVT[u]; i++) {
            const int t = mk_dec_vt_(c, u)[i];
            tu |= mk_dec_has_(c, t, opp[0]) && mk_dec_has_(c, t, opp[1]) && !mk_dec_has_(c, t, v);
        }
    if (need == 2)
        for (int i = 0; i < c->NVT[v]; i++) {
            const int t = mk_dec_vt_(c, v)[i];
            tv |= mk_dec_has_(c, t, opp[0]) && mk_dec_has_(c, t, opp[1]) && !mk_dec_has_(c, t, u);
        }
    if (tu && tv) { c->ref[0]++; return 0; }
#else
    (void)common;
#endif
    mk_dec_fan_(c, u, 1);
    const double oldc = c->foldc;
    const mk_dv3_ pv = mk_dec_pt_(c, v);
    for (int i = 0; i < c->NVT[u]; i++) {
        const int t = mk_dec_vt_(c, u)[i], *f = mk_dec_f_(c, t);
        if (mk_dec_has_(c, t, v)) continue;
        mk_dv3_ p1[3];
        for (int k = 0; k < 3; k++) p1[k] = f[k] == u ? pv : mk_dec_pt_(c, f[k]);
        const mk_dv3_ n0 = c->fn0[i], n1 = mk_dtri_n_(p1[0], p1[1], p1[2]);
        const double l0 = c->fl0[i], l1 = mk_dlen_(n1);
        c->fn1[i] = n1;
        if (!(l1 > 1e-8 * l0) || (!c->relax && mk_ddot_(n0, n1) < c->mindot * l0 * l1)) { c->ref[1]++; return 0; }
        int g[3];
        for (int k = 0; k < 3; k++) g[k] = f[k] == u ? v : f[k];
        const mk_dv3_ n0i = mk_dnorm_(mk_dadd_(mk_dadd_(mk_dv_(c->N0[3 * (size_t)g[0]], c->N0[3 * (size_t)g[0] + 1], c->N0[3 * (size_t)g[0] + 2]),
                                                       mk_dv_(c->N0[3 * (size_t)g[1]], c->N0[3 * (size_t)g[1] + 1], c->N0[3 * (size_t)g[1] + 2])),
                                               mk_dv_(c->N0[3 * (size_t)g[2]], c->N0[3 * (size_t)g[2] + 1], c->N0[3 * (size_t)g[2] + 2])));
        if (!c->relax && mk_ddot_(n1, n0i) < c->mindot * l1) { c->ref[1]++; return 0; }
        for (int k = 0; k < 3; k++)
            if (g[k] != v && c->maxedge2 > 0) {
                const mk_dv3_ e = mk_dsub_(mk_dec_pt_(c, g[k]), pv);
                if (!c->relax && mk_ddot_(e, e) > c->maxedge2) { c->ref[2]++; return 0; }
            }
        const double mc = mk_dmax_cos_(p1[0], p1[1], p1[2]);
        if (!c->relax && mc > c->cmin && mc > oldc + 1e-12) { c->ref[3]++; return 0; }
#if !defined(MK_MUTATE) || MK_MUTATE != 320
        if (c->minalt > 0 && !c->relax) { /* low: under minalt high over its longest edge and lower than the fan's lowest */
            const double e1 = mk_dec_e2_(p1[0], p1[1], p1[2]);
            if (l1 * l1 < c->minalt * c->minalt * e1 && l1 * l1 < c->flow * e1) { c->ref[8]++; return 0; }
        }
#endif
    }
    if (c->relax) return 1;
#if !defined(MK_MUTATE) || MK_MUTATE != 161
    mk_dec_fan_(c, u, 2);
    if (mk_dec_fold_up_(c, u, v)) { c->ref[4]++; return 0; }
#endif
    return 1;
}
/* the input faces into the cells of their boxes (one bucket for those too large), side by side by bucket */
static inline void mk_dec_nbuild_(mk_dec_ *c)
{
    int gh = 1024;
    while (gh < c->nt / 2 && gh < (1 << 26)) gh *= 2;
    c->ngh = gh, c->nst = calloc((size_t)gh + 1, sizeof *c->nst);
    if (!c->nst) { c->oom = 1; return; }
    size_t need = 0;
    for (int pass = 0; pass < 2; pass++) {
        for (int t = 0; t < c->nt; t++) {
            const int *f = c->T0 + 3 * (size_t)t;
            const mk_dv3_ d[3] = {mk_dec_pt_(c, f[0]), mk_dec_pt_(c, f[1]), mk_dec_pt_(c, f[2])};
            int lo[3], hi[3];
            if (mk_dec_cells_(c, d, lo, hi)) lo[0] = lo[1] = lo[2] = hi[0] = hi[1] = hi[2] = -1;
            for (int x = lo[0]; x <= hi[0]; x++)
                for (int y = lo[1]; y <= hi[1]; y++)
                    for (int z = lo[2]; z <= hi[2]; z++) {
                        const size_t k = mk_dec_hkey_(x, y, z, gh);
                        if (!pass) c->nst[k + 1]++;
                        else c->ntr[c->nst[k]++] = t;
                    }
        }
        if (pass) break;
        for (int k = 1; k <= gh; k++) {
            need += (size_t)c->nst[k];
            if (need > INT_MAX) { c->oom = 1; return; }
            c->nst[k] = (int)need;
        }
        c->ntr = malloc(sizeof *c->ntr * (need ? need : 1));
        if (!c->ntr) { c->oom = 1; return; }
    }
    for (int k = gh; k > 0; k--) c->nst[k] = c->nst[k - 1];
    c->nst[0] = 0;
}
/* 1 when an input face lies within the bound of p (mk_closest_pt_tri, the face's corners in order): the faces of the
   cells within the bound and its rounding of p, each past a box test widened by that rounding */
static inline int mk_dec_near_(const mk_dec_ *c, const double *p)
{
    const double r = sqrt(c->dist2) + c->xpad, bs = mk_near_slack_(c->dist2, c->xpad);
    const mk_dv3_ q[3] = {mk_dv_(p[0] - r, p[1] - r, p[2] - r), mk_dv_(p[0] + r, p[1] + r, p[2] + r), mk_dv_(p[0], p[1], p[2])};
    int lo[3], hi[3];
    const int all = mk_dec_cells_(c, q, lo, hi);
    for (int x = all ? 0 : lo[0] - 1; x <= (all ? 0 : hi[0]); x++)
        for (int y = all || x < lo[0] ? -1 : lo[1]; y <= (all || x < lo[0] ? -1 : hi[1]); y++)
            for (int z = all || x < lo[0] ? -1 : lo[2]; z <= (all || x < lo[0] ? -1 : hi[2]); z++) {
                const int bx = all || x < lo[0] ? -1 : x;
                const size_t k = mk_dec_hkey_(bx, y, z, c->ngh);
                const int i0 = all ? 0 : c->nst[k], i1 = all ? c->nt : c->nst[k + 1];
                for (int i = i0; i < i1; i++) {
                    const int *f = c->T0 + 3 * (size_t)(all ? i : c->ntr[i]);
                    const double *a = c->P + 3 * (size_t)f[0], *b = c->P + 3 * (size_t)f[1], *d = c->P + 3 * (size_t)f[2];
                    double e2 = 0;
                    for (int j = 0; j < 3; j++) {
                        const double l = mk_dec_lo_(a[j], mk_dec_lo_(b[j], d[j])), h = mk_dec_hi_(a[j], mk_dec_hi_(b[j], d[j]));
                        const double g = p[j] < l ? l - p[j] : p[j] > h ? p[j] - h : 0;
                        e2 += g * g;
                    }
                    double w[3], o[3];
                    if (e2 <= bs && mk_closest_pt_tri(p, a, b, d, w, o) < c->dist2) return 1;
                }
            }
    return 0;
}
/* distance and crossing; a pass leaves the owned vertices' new triangles pending for mk_dec_collapse_ */
static inline int mk_dec_full_(mk_dec_ *c, int u, int v)
{
    int m = 0;
    const mk_dv3_ pv = mk_dec_pt_(c, v);
    c->pend_u = c->pend_v = -1, c->pend_n = 0;
    for (int i = 0; i < c->NVT[u]; i++) {
        const int t = mk_dec_vt_(c, u)[i], *f = mk_dec_f_(c, t);
        if (mk_dec_has_(c, t, v)) continue;
        c->ntb[m] = t;
        for (int k = 0; k < 3; k++) {
            const mk_dv3_ p = f[k] == u ? pv : mk_dec_pt_(c, f[k]);
            c->nc[m][3 * k] = p.x, c->nc[m][3 * k + 1] = p.y, c->nc[m][3 * k + 2] = p.z;
        }
        m++;
    }
    double worst = 0;
#if defined(MK_MUTATE) && MK_MUTATE == 229
#define MK_DEC_OWN_EARLY_ if (worst > c->maxown2) c->maxown2 = worst;
#else
#define MK_DEC_OWN_EARLY_
#endif
#if !defined(MK_MUTATE) || MK_MUTATE != 162
    if (c->dist2 > 0) {
        for (int i = -1; i < c->NVT[u]; i++)
            for (int p = i < 0 ? u : c->ohead[mk_dec_vt_(c, u)[i]]; p >= 0; p = i < 0 ? -1 : c->onext[p]) {
                double bd = INFINITY;
                int bt = -1;
                for (int j = 0; j < m; j++) {
                    double w[3], q[3];
                    const double d2 = mk_closest_pt_tri(c->P + 3 * (size_t)p, c->nc[j], c->nc[j] + 3, c->nc[j] + 6, w, q);
                    if (d2 < bd) bd = d2, bt = c->ntb[j];
                }
#if !defined(MK_MUTATE) || MK_MUTATE != 201
                /* when every triangle of u's holds v (a rim ear), no new triangle can own what u owns: in the
                   clean-up's relaxed try the nearest of v's triangles that stays owns it */
                for (int i2 = 0; !m && c->relax && i2 < c->NVT[v]; i2++) {
                    const int t2 = mk_dec_vt_(c, v)[i2], *g2 = mk_dec_f_(c, t2);
                    if (mk_dec_has_(c, t2, u)) continue;
                    double w[3], q[3];
                    const double d2 = mk_closest_pt_tri(c->P + 3 * (size_t)p, c->P + 3 * (size_t)g2[0], c->P + 3 * (size_t)g2[1], c->P + 3 * (size_t)g2[2], w, q);
                    if (d2 < bd) bd = d2, bt = t2;
                }
                if (bt < 0) { c->ref[5]++; return 0; }
#endif
                if (!(bd <= c->dist2) && !c->relax) { c->ref[5]++; return 0; }
                if (bd > worst) worst = bd;
                if (c->pend_n == c->pend_c) {
                    const int nc = c->pend_c ? 2 * c->pend_c : 4096;
                    c->pend_p = mk_dec_grow_(c, c->pend_p, (size_t)nc, sizeof *c->pend_p);
                    c->pend_t = mk_dec_grow_(c, c->pend_t, (size_t)nc, sizeof *c->pend_t);
                    if (c->oom) return 0;
                    c->pend_c = nc;
                }
                c->pend_p[c->pend_n] = p, c->pend_t[c->pend_n++] = bt;
            }
        for (int j = 0; j < m && !c->relax; j++) {
            double lm = 0;
            for (int k = 0; k < 3; k++) {
                double e = 0;
                for (int q = 0; q < 3; q++) {
                    const double d = c->nc[j][3 * k + q] - c->nc[j][3 * ((k + 1) % 3) + q];
                    e += d * d;
                }
                lm = mk_dec_hi_(lm, e);
            }
            const double nn = ceil(sqrt(lm) / c->samp);
            const int n = !(nn >= 2) ? 2 : nn > 4 ? 4 : (int)nn;
            for (int a = 0; a <= n; a++)
                for (int b = 0; a + b <= n; b++) {
                    const int cc = n - a - b;
                    if (a == n || b == n || cc == n) continue;
                    double p[3];
                    for (int k = 0; k < 3; k++) p[k] = (a * c->nc[j][k] + b * c->nc[j][3 + k] + cc * c->nc[j][6 + k]) / n;
                    double e2 = INFINITY; /* a corner is an input point: near enough, it settles the test */
                    for (int q = 0; q < 9; q += 3) {
                        const double dx = p[0] - c->nc[j][q], dy = p[1] - c->nc[j][q + 1], dz = p[2] - c->nc[j][q + 2];
                        e2 = mk_dec_lo_(e2, dx * dx + dy * dy + dz * dz);
                    }
                    if (e2 < c->dnear) continue;
                    if (!mk_dec_near_(c, p)) { c->ref[6]++; return 0; }
                }
        }
    }
#endif
    MK_DEC_OWN_EARLY_
    if (mk_dec_xfan_(c, u, v, m)) { c->ref[7]++; return 0; }
#if !(defined(MK_MUTATE) && MK_MUTATE == 229) /* 229: counted before the crossing test can refuse the collapse */
    if (worst > c->maxown2) c->maxown2 = worst; /* only a collapse that happens counts */
#endif
    c->pend_u = u, c->pend_v = v;
    return 1;
}
static inline double mk_dec_qeval_(const double *q, mk_dv3_ p)
{
    return q[0] * p.x * p.x + 2 * q[1] * p.x * p.y + 2 * q[2] * p.x * p.z + 2 * q[3] * p.x + q[4] * p.y * p.y + 2 * q[5] * p.y * p.z + 2 * q[6] * p.y +
           q[7] * p.z * p.z + 2 * q[8] * p.z + q[9];
}
static inline double mk_dec_cost_(const mk_dec_ *c, int u, int v)
{
    double q[10];
    for (int k = 0; k < 10; k++) q[k] = c->Q[10 * (size_t)u + k] + c->Q[10 * (size_t)v + k];
    const mk_dv3_ e = mk_dsub_(mk_dec_pt_(c, u), mk_dec_pt_(c, v));
    const double l2 = mk_ddot_(e, e);
    return mk_dec_qeval_(q, mk_dec_pt_(c, v)) + 1e-3 * l2 * l2;
}
/* the queue: a binary heap of one entry a vertex, the cheapest first and of equal costs the lower vertex, with each
   vertex's place in it (hpos, -1 none), its target and its skip count */
static inline int mk_dec_hless_(mk_dec_ent_ a, mk_dec_ent_ b) { return a.c < b.c || (a.c == b.c && a.u < b.u); }
static inline void mk_dec_hput_(mk_dec_ *c, int i, mk_dec_ent_ e) { c->heap[i] = e, c->hpos[e.u] = i; }
static inline void mk_dec_hmove_(mk_dec_ *c, int i, mk_dec_ent_ e)
{
    if (i && mk_dec_hless_(e, c->heap[(i - 1) / 2])) {
        for (; i && mk_dec_hless_(e, c->heap[(i - 1) / 2]); i = (i - 1) / 2) mk_dec_hput_(c, i, c->heap[(i - 1) / 2]);
    } else
        for (;;) {
            int k = 2 * i + 1;
            if (k >= c->nh) break;
            if (k + 1 < c->nh && mk_dec_hless_(c->heap[k + 1], c->heap[k])) k++;
            if (!mk_dec_hless_(c->heap[k], e)) break;
            mk_dec_hput_(c, i, c->heap[k]), i = k;
        }
    mk_dec_hput_(c, i, e);
}
/* u queued at cost cs onto v, skipping skip; mk_dec_hout_: u out of the queue */
static inline void mk_dec_hset_(mk_dec_ *c, int u, double cs, int v, int skip)
{
    c->htv[u] = v, c->hsk[u] = skip;
    mk_dec_hmove_(c, c->hpos[u] >= 0 ? c->hpos[u] : c->nh++, (mk_dec_ent_){cs, u});
}
static inline void mk_dec_hout_(mk_dec_ *c, int u)
{
    const int i = c->hpos[u];
    if (i < 0) return;
    c->hpos[u] = -1;
    const mk_dec_ent_ e = c->heap[--c->nh];
    if (i < c->nh) mk_dec_hmove_(c, i, e);
}
/* queue u's (skip+1)-th cheapest target that passes the cheap refusals, or take u out */
static inline void mk_dec_rate_(mk_dec_ *c, int u, int skip)
{
    if (c->deadv[u]) { mk_dec_hout_(c, u); return; }
    int nb[MK_DEC_MAXN];
    double cs[MK_DEC_MAXN];
    const int n = mk_dec_nbrs_(c, u, nb);
    for (int i = 0; i < n; i++) cs[i] = mk_dec_cost_(c, u, nb[i]);
    for (int got = 0;;) {
        int b = -1;
        for (int i = 0; i < n; i++)
            if (nb[i] >= 0 && (b < 0 || cs[i] < cs[b])) b = i;
        if (b < 0) { mk_dec_hout_(c, u); return; }
        const int v = nb[b];
        nb[b] = -1;
        if (!mk_dec_cheap_(c, u, v)) continue;
        if (got++ == skip) { mk_dec_hset_(c, u, cs[b], v, skip); return; }
    }
}
static inline void mk_dec_collapse_(mk_dec_ *c, int u, int v)
{
    int fan[MK_DEC_MAXN];
    const int n = c->NVT[u] < MK_DEC_MAXN ? c->NVT[u] : MK_DEC_MAXN;
    memcpy(fan, mk_dec_vt_(c, u), sizeof *fan * (size_t)n);
    c->ncol++;
    for (int i = 0; i < n; i++) {
        const int t = fan[i];
        int *f = mk_dec_f_(c, t);
        c->ohead[t] = -1;
        if (mk_dec_has_(c, t, v)) {
            c->deadt[t] = 1;
            for (int k = 0; k < 3; k++)
                if (f[k] != u) mk_dec_fan_del_(c, f[k], t);
        } else {
            for (int k = 0; k < 3; k++)
                if (f[k] == u) f[k] = v;
            mk_dec_fan_add_(c, v, t);
        }
    }
    if (c->pend_u == u && c->pend_v == v)
        for (int i = 0; i < c->pend_n; i++) c->onext[c->pend_p[i]] = c->ohead[c->pend_t[i]], c->ohead[c->pend_t[i]] = c->pend_p[i];
    c->pend_u = c->pend_v = -1;
    c->NVT[u] = c->CVT[u] = 0, c->deadv[u] = 1;
    for (int k = 0; k < 10; k++) c->Q[10 * (size_t)v + k] += c->Q[10 * (size_t)u + k];
    for (int i = 0; i < c->NVT[v]; i++) c->tver[mk_dec_vt_(c, v)[i]]++, mk_dec_gins_(c, mk_dec_vt_(c, v)[i]);
}
/* the whole mesh's fold lengths and areas by region */
static inline void mk_dec_total_(mk_dec_ *c, double fl[MK_DEC_REGIONS][2], double ar[MK_DEC_REGIONS])
{
    for (int r = 0; r < MK_DEC_REGIONS; r++) fl[r][0] = fl[r][1] = ar[r] = 0;
    for (int t = 0; t < c->nt; t++) {
        if (c->deadt[t]) continue;
        const int *f = mk_dec_f_(c, t);
        const mk_dv3_ a = mk_dec_pt_(c, f[0]), b = mk_dec_pt_(c, f[1]), cc = mk_dec_pt_(c, f[2]), n = mk_dtri_n_(a, b, cc), un = mk_dnorm_(n);
        ar[mk_dec_region_(c, mk_dmul_(mk_dadd_(mk_dadd_(a, b), cc), 1.0 / 3))] += 0.5 * mk_dlen_(n);
        for (int k = 0; k < 3; k++) {
            const int x = f[k], y = f[(k + 1) % 3];
            if (x > y) continue;
            for (int q = 0; q < c->NVT[x]; q++) {
                const int t2 = mk_dec_vt_(c, x)[q];
                if (t2 == t || !mk_dec_has_(c, t2, y)) continue;
                const int *g = mk_dec_f_(c, t2);
                const double d = mk_ddot_(un, mk_dnorm_(mk_dtri_n_(mk_dec_pt_(c, g[0]), mk_dec_pt_(c, g[1]), mk_dec_pt_(c, g[2]))));
                const mk_dv3_ px = mk_dec_pt_(c, x), py = mk_dec_pt_(c, y);
                const double L = mk_dlen_(mk_dsub_(px, py));
                const int r = mk_dec_region_(c, mk_dmul_(mk_dadd_(px, py), 0.5));
                if (d < 0) fl[r][0] += L;
                if (d < c->cos2) fl[r][1] += L;
                break;
            }
        }
    }
}
/* triangle d is under minalt high over its longest edge (a needle for the clean-up) */
static inline int mk_dec_needle_(const mk_dv3_ *d, double minalt)
{
    const double a2 = mk_dlen_(mk_dtri_n_(d[0], d[1], d[2]));
    const double l = mk_dec_hi_(mk_dlen_(mk_dsub_(d[0], d[1])), mk_dec_hi_(mk_dlen_(mk_dsub_(d[1], d[2])), mk_dlen_(mk_dsub_(d[2], d[0]))));
    return !(l > 0) || a2 / l < minalt;
}
static inline int mk_dec_bad_(mk_dec_ *c, int t, double minalt)
{
    const int *f = mk_dec_f_(c, t);
    const mk_dv3_ d[3] = {mk_dec_pt_(c, f[0]), mk_dec_pt_(c, f[1]), mk_dec_pt_(c, f[2])};
    if (minalt > 0 && mk_dec_needle_(d, minalt)) return 1;
    return mk_dec_crosses_(c, d, f, -1);
}
/* 1 when some live triangle is under minalt high (minalt > 0) or crosses another, each pair tried once, from its
   lower-numbered triangle: with 0 a clean-up round would find nothing */
static inline int mk_dec_anybad_(mk_dec_ *c, double minalt)
{
    for (int t = 0; t < c->nt && !c->oom; t++) {
        if (c->deadt[t]) continue;
        const int *f = mk_dec_f_(c, t);
        const mk_dv3_ d[3] = {mk_dec_pt_(c, f[0]), mk_dec_pt_(c, f[1]), mk_dec_pt_(c, f[2])};
        if (minalt > 0 && mk_dec_needle_(d, minalt)) return 1;
        int lo[3], hi[3];
        const double bl[3] = {mk_dec_lo_(d[0].x, mk_dec_lo_(d[1].x, d[2].x)), mk_dec_lo_(d[0].y, mk_dec_lo_(d[1].y, d[2].y)), mk_dec_lo_(d[0].z, mk_dec_lo_(d[1].z, d[2].z))};
        const double bh[3] = {mk_dec_hi_(d[0].x, mk_dec_hi_(d[1].x, d[2].x)), mk_dec_hi_(d[0].y, mk_dec_hi_(d[1].y, d[2].y)), mk_dec_hi_(d[0].z, mk_dec_hi_(d[1].z, d[2].z))};
        if (mk_dec_cells_(c, d, lo, hi)) {
            for (int t2 = t + 1; t2 < c->nt; t2++)
                if (mk_dec_xone_(c, d, f, -1, bl, bh, t2, 0)) return 1;
            continue;
        }
        mk_dec_hclear_(c);
        for (int x = lo[0] - 1; x <= hi[0]; x++) /* x = lo[0] - 1: the bucket of the large triangles */
            for (int y = x < lo[0] ? -1 : lo[1]; y <= (x < lo[0] ? -1 : hi[1]); y++)
                for (int z = x < lo[0] ? -1 : lo[2]; z <= (x < lo[0] ? -1 : hi[2]); z++) {
                    const size_t k = x < lo[0] ? mk_dec_gkey_(c, -1, -1, -1) : mk_dec_gkey_(c, x, y, z);
                    for (int i = c->cst[k], e = c->ghead[k], t2; (t2 = mk_dec_gwalk_(c, k, &i, &e)) >= 0;)
                        if (t2 > t && mk_dec_xone_(c, d, f, -1, bl, bh, t2, 1)) return 1;
                }
    }
    return c->oom;
}
static inline int mk_dec_kill_(mk_dec_ *c, int t, double maxlen)
{
    const int *f = mk_dec_f_(c, t);
    int ord[3] = {0, 1, 2};
    double len[3];
    for (int k = 0; k < 3; k++) len[k] = mk_dlen_(mk_dsub_(mk_dec_pt_(c, f[k]), mk_dec_pt_(c, f[(k + 1) % 3])));
    for (int i = 0; i < 3; i++)
        for (int j = i + 1; j < 3; j++)
            if (len[ord[j]] < len[ord[i]]) { const int x = ord[i]; ord[i] = ord[j], ord[j] = x; }
    for (int i = 0; i < 3; i++) {
        if (len[ord[i]] > maxlen) break;
        const int a = f[ord[i]], b = f[(ord[i] + 1) % 3];
        for (int dir = 0; dir < 2; dir++) {
            const int u = dir ? b : a, v = dir ? a : b;
            if (mk_dec_cheap_(c, u, v) && mk_dec_full_(c, u, v)) { mk_dec_collapse_(c, u, v); return 1; }
        }
    }
    return 0;
}
static inline long mk_dec_cleanup_(mk_dec_ *c, double minalt, long *alive, long *left)
{
    long fixed = 0;
    for (int round = 0; round < 6 && !c->oom; round++) {
        long n = 0, got = 0;
        if (!mk_dec_anybad_(c, minalt)) { *left = 0; break; }
        for (int t = 0; t < c->nt; t++) {
            if (c->deadt[t] || !mk_dec_bad_(c, t, minalt)) continue;
            n++;
            if (mk_dec_kill_(c, t, INFINITY)) got++, (*alive)--;
            else { /* a crossing triangle by any edge, a needle only by an edge under 20 minalt */
                c->relax = 1;
                if (mk_dec_kill_(c, t, mk_dec_bad_(c, t, 0) ? INFINITY : c->relaxlen)) got++, (*alive)--, c->nrelax++;
                c->relax = 0;
            }
        }
        fixed += got, *left = n - got;
        if (!got) break;
    }
    return fixed;
}
static inline void mk_dec_free_(mk_dec_ *c)
{
    free(c->VA), free(c->VO), free(c->F), free(c->NVT), free(c->CVT), free(c->MARK), free(c->N0), free(c->Q);
    free(c->deadv), free(c->deadt), free(c->bnd), free(c->nc), free(c->ntb), free(c->ghead), free(c->cst);
    free(c->ce), free(c->go), free(c->tver), free(c->hkey), free(c->hst), free(c->ohead), free(c->onext), free(c->pend_p), free(c->pend_t);
    free(c->heap), free(c->hpos), free(c->htv), free(c->hsk);
    free(c->nst), free(c->ntr);
}

/* Decimates (V[3 nv], T[3 nt]) under o (NULL: mk_dec_defaults with mm = 1e-3, metres). Out, in the caller's arrays:
   alive[nv] (1 kept, 0 removed) and the result's *nto faces in To (capacity 3 nt ints) in the input's vertex numbers;
   *rep if rep is not NULL. Returns 0, or -1 for a bad mesh (a size under 1 or over INT_MAX / 4, an index out of range,
   a face with a repeated corner, a coordinate not finite, or coordinates past 2^30 grid cells), bad options (a negative
   or NaN bound, turn outside [-1, 1], minangle or fold outside [0, 180), nregion outside 1 .. MK_DEC_REGIONS or above 1
   without region, target < 0), or no memory; alive and To are then untouched and *rep zeroed. mk_decimate_map does the
   same and, when face is not NULL, writes face[*nto] (capacity nt): the input face each result face is, its corners in
   the same order, a removed corner replaced by the vertex it went onto. A face's group or any other data carries over
   by it. */
static inline int mk_decimate_map(const double *V, int nv, const int *T, int nt, const mk_dec_opts *o, unsigned char *alive, int *To, int *nto,
                                  int *face, mk_dec_report *rep)
{
    mk_dec_opts od;
    if (rep) memset(rep, 0, sizeof *rep);
    if (!o) mk_dec_defaults(&od, 1e-3), o = &od;
    if (!V || !T || !alive || !To || !nto || nv < 1 || nt < 1 || nv > INT_MAX / 4 || nt > INT_MAX / 4) return -1;
    if (!(o->dist >= 0) || !(o->maxedge >= 0) || !(o->minalt >= 0) || !(o->minangle >= 0 && o->minangle < 180) || !(o->fold >= 0 && o->fold < 180) ||
        !(o->turn >= -1 && o->turn <= 1) || !isfinite(o->dist) || !isfinite(o->maxedge) || !isfinite(o->minalt) || o->target < 0 || o->nregion < 1 ||
        o->nregion > MK_DEC_REGIONS || (o->nregion > 1 && !o->region))
        return -1;
    double lo[3] = {INFINITY, INFINITY, INFINITY}, hi[3] = {-INFINITY, -INFINITY, -INFINITY};
    for (size_t q = 0; q < 3 * (size_t)nv; q++) {
        if (!isfinite(V[q])) return -1;
        lo[q % 3] = fmin(lo[q % 3], V[q]), hi[q % 3] = fmax(hi[q % 3], V[q]);
    }
    for (int t = 0; t < nt; t++) {
        const int *f = T + 3 * (size_t)t;
        for (int k = 0; k < 3; k++)
            if (f[k] < 0 || f[k] >= nv) return -1;
        if (f[0] == f[1] || f[1] == f[2] || f[2] == f[0]) return -1;
    }
    double el = 0;
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) {
            const double *a = V + 3 * (size_t)T[3 * (size_t)t + k], *b = V + 3 * (size_t)T[3 * (size_t)t + (k + 1) % 3];
            el += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
        }
    el /= 3.0 * nt; /* the mean edge */
    mk_dec_ c = {0};
    c.P = V, c.nv = nv, c.nt = nt, c.o = o, c.keep = o->keep, c.nreg = o->nregion, c.pend_u = c.pend_v = -1;
#if defined(MK_MUTATE) && MK_MUTATE == 228
    c.dist2 = o->dist * o->dist, c.mindot = o->turn; /* 228: a small bound squares to 0, which means none */
#else
    /* a positive bound stays positive when squared */
    c.dist2 = o->dist > 0 ? fmax(o->dist * o->dist, 0x1p-1074) : 0, c.mindot = o->turn;
#endif
    c.maxedge2 = o->maxedge > 0 ? fmax(o->maxedge * o->maxedge, 0x1p-1074) : 0;
    c.cos2 = o->fold > 0 ? cos(o->fold * 3.14159265358979323846 / 180) : 0, c.cmin = o->minangle > 0 ? cos(o->minangle * 3.14159265358979323846 / 180) : 2;
    c.gc = fmin(o->maxedge > 0 ? 0.5 * o->maxedge : INFINITY, 4 * el), c.samp = o->dist / 0.3, c.ftol = 1e-9 * el, c.minalt = o->minalt, c.relaxlen = 20 * o->minalt;
    if (!(c.gc > 0)) c.gc = 1;
#if defined(MK_MUTATE) && MK_MUTATE == 227
    c.gs = c.gc;
#else
    c.gs = fmax(c.gc, 0.5 * el); /* at least half the mean edge: an input triangle covers a few cells */
#endif
    double big = 0;
    for (int k = 0; k < 3; k++) {
        c.glo[k] = lo[k], big = fmax(big, fmax(fabs(lo[k]), fabs(hi[k])));
        if (!((hi[k] - lo[k]) / c.gc < 0x1p30)) return -1;
    }
    c.xpad = 0x1p-40 * big; /* boxes meet across rounding, so the box test never drops a pair the edge test would count */
    /* a sample this near a corner is inside the bound by more than any rounding of the distance (mk_geom's slack is 2^-49 big) */
    c.dnear = sqrt(c.dist2) > c.xpad ? (sqrt(c.dist2) - c.xpad) * (sqrt(c.dist2) - c.xpad) : 0;
    int ret = -1;
    c.F = malloc(sizeof *c.F * 3 * (size_t)nt);
    c.NVT = calloc((size_t)nv, sizeof *c.NVT), c.CVT = calloc((size_t)nv, sizeof *c.CVT);
    c.MARK = calloc((size_t)nv, sizeof *c.MARK), c.N0 = calloc(3 * (size_t)nv, sizeof *c.N0);
    c.Q = calloc(10 * (size_t)nv, sizeof *c.Q), c.deadv = calloc((size_t)nv, 1), c.deadt = calloc((size_t)nt, 1), c.bnd = calloc((size_t)nv, 1);
    c.nc = malloc(sizeof *c.nc * MK_DEC_MAXN), c.ntb = malloc(sizeof *c.ntb * MK_DEC_MAXN), c.tver = calloc((size_t)nt, sizeof *c.tver);
    c.hcap = 1 << 12, c.hkey = malloc(sizeof *c.hkey * (size_t)c.hcap), c.hst = calloc((size_t)c.hcap, sizeof *c.hst);
    c.ohead = malloc(sizeof *c.ohead * (size_t)nt), c.onext = malloc(sizeof *c.onext * (size_t)nv);
    c.heap = malloc(sizeof *c.heap * (size_t)nv), c.hpos = malloc(sizeof *c.hpos * (size_t)nv);
    c.htv = malloc(sizeof *c.htv * (size_t)nv), c.hsk = malloc(sizeof *c.hsk * (size_t)nv);
    if (!c.F || !c.NVT || !c.CVT || !c.MARK || !c.N0 || !c.Q || !c.deadv || !c.deadt || !c.bnd || !c.nc || !c.ntb || !c.tver ||
        !c.hkey || !c.hst || !c.ohead || !c.onext || !c.heap || !c.hpos || !c.htv || !c.hsk)
        goto out;
    for (int i = 0; i < nv; i++) c.hpos[i] = -1;
    memcpy(c.F, T, sizeof *c.F * 3 * (size_t)nt);
    mk_dec_fans_(&c);
    if (c.oom) goto out;
    c.T0 = T; /* the input's own faces: T is not written */
    for (int t = 0; t < nt; t++) c.ohead[t] = -1;
    for (int i = 0; i < nv; i++) c.onext[i] = -1;
    for (int t = 0; t < nt; t++) {
        const int *f = mk_dec_f_(&c, t);
        const mk_dv3_ nn = mk_dtri_n_(mk_dec_pt_(&c, f[0]), mk_dec_pt_(&c, f[1]), mk_dec_pt_(&c, f[2]));
        const double la = mk_dlen_(nn);
        if (!(la > 0)) continue;
        const mk_dv3_ un = mk_dmul_(nn, 1 / la);
        const double d = -mk_ddot_(un, mk_dec_pt_(&c, f[0])), w = 0.5 * la, n[4] = {un.x, un.y, un.z, d};
        for (int q = 0; q < 3; q++) {
            double *Qv = c.Q + 10 * (size_t)f[q];
            int k = 0;
            for (int i = 0; i < 4; i++)
                for (int j = i; j < 4; j++) Qv[k++] += w * n[i] * n[j];
            c.N0[3 * (size_t)f[q]] += nn.x, c.N0[3 * (size_t)f[q] + 1] += nn.y, c.N0[3 * (size_t)f[q] + 2] += nn.z;
        }
    }
    if (c.oom) goto out;
    for (int i = 0; i < nv; i++) {
        const mk_dv3_ n = mk_dnorm_(mk_dv_(c.N0[3 * (size_t)i], c.N0[3 * (size_t)i + 1], c.N0[3 * (size_t)i + 2]));
        c.N0[3 * (size_t)i] = n.x, c.N0[3 * (size_t)i + 1] = n.y, c.N0[3 * (size_t)i + 2] = n.z;
    }
    long nrim = 0;
    for (int u = 0; u < nv; u++) {
        for (int i = 0; i < c.NVT[u] && !c.bnd[u]; i++)
            for (int k = 0; k < 3 && !c.bnd[u]; k++) {
                const int x = mk_dec_f_(&c, mk_dec_vt_(&c, u)[i])[k];
                if (x == u) continue;
                int n = 0;
                for (int j = 0; j < c.NVT[u]; j++) n += mk_dec_has_(&c, mk_dec_vt_(&c, u)[j], x);
                if (n == 1) c.bnd[u] = 1;
            }
        nrim += c.bnd[u];
    }
    if (c.dist2 > 0) mk_dec_nbuild_(&c);
    mk_dec_grebuild_(&c);
    if (c.oom) goto out;
    mk_dec_report r = {0};
    r.rim0 = nrim;
    mk_dec_total_(&c, r.fold0, r.area0);
    long nalive = nv;
    if (o->clean) r.clean_before = mk_dec_cleanup_(&c, 0, &nalive, &r.left_before);
    for (int u = 0; u < nv && !c.oom; u++) mk_dec_rate_(&c, u, 0);
    while (nalive > o->target && !c.oom && c.nh) {
        const int u = c.heap[0].u, v = c.htv[u], k = c.hsk[u];
        mk_dec_hout_(&c, u);
        if (c.deadv[u] || c.deadv[v]) continue;
        if (!mk_dec_cheap_(&c, u, v)) { mk_dec_rate_(&c, u, 0); continue; }
        if (!mk_dec_full_(&c, u, v)) { mk_dec_rate_(&c, u, k + 1); continue; }
        int nb[MK_DEC_MAXN];
        mk_dec_collapse_(&c, u, v), nalive--;
        if (c.gn > c.gcap / 2 + c.gcap / 4) mk_dec_grebuild_(&c);
        const int m = mk_dec_nbrs_(&c, v, nb);
        mk_dec_rate_(&c, v, 0);
        for (int i = 0; i < m; i++) mk_dec_rate_(&c, nb[i], 0);
    }
    if (o->clean && !c.oom) r.clean_after = mk_dec_cleanup_(&c, o->minalt, &nalive, &r.left_after);
    if (c.oom) goto out;
    mk_dec_total_(&c, r.fold1, r.area1);
    int m = 0, kv = 0;
    for (int i = 0; i < nv; i++) alive[i] = !c.deadv[i], kv += !c.deadv[i], r.rim1 += c.bnd[i] && !c.deadv[i];
    for (int t = 0; t < nt; t++)
        if (!c.deadt[t]) {
            if (face) face[m] = t;
#if defined(MK_MUTATE) && MK_MUTATE == 321
            if (face && m) face[m - 1] = t; /* 321: each face mapped to the next one's input face */
#endif
            memcpy(To + 3 * (size_t)m++, mk_dec_f_(&c, t), 3 * sizeof *To);
        }
    *nto = m;
    r.nv = kv, r.nt = m, r.maxown = sqrt(c.maxown2), r.relaxed = c.nrelax;
    memcpy(r.refused, c.ref, sizeof r.refused);
    if (rep) *rep = r;
    ret = 0;
out:
    mk_dec_free_(&c);
    return ret;
}
static inline int mk_decimate(const double *V, int nv, const int *T, int nt, const mk_dec_opts *o, unsigned char *alive, int *To, int *nto,
                              mk_dec_report *rep)
{
    return mk_decimate_map(V, nv, T, nt, o, alive, To, nto, NULL, rep);
}
#endif
