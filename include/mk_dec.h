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
            triangle farther than dist from the input surface (through mk_geom's tree: the other side);
     cross  make a new triangle cross a live one it shares no vertex with (an edge through the other's interior,
            mk_seg_tri_strict), so an input with no crossing pairs gives none.
   Each vertex queues its cheapest target that passes the cheap refusals (turn, edge, angle, fold); distance and
   crossing are tried when it comes up, and on a refusal the next-best target is queued. Stops at `target` vertices or
   when nothing more is allowed. keep (optional, a byte per vertex): 0 free; 1 or 3 never removed; 2 a seam vertex,
   removed only along the seam (onto a 2 or 3, the held side's triangles staying in their planes). An open surface's
   rim vertex goes only along the rim. clean (1): a triangle caught in a crossing (before decimation) or under minalt
   high (after) is removed by collapsing one of its edges under the same refusals; one that still stands is tried once
   more under the link, keep, seam and crossing refusals alone (a needle only by an edge under 20 minalt).
   Experimental: its results may change between releases. Single-threaded.
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
    double minalt;   /* the clean-up's needle height; 0: no needles removed */
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
    long refused[8]; /* refusals by: link, turn, edge, angle, fold, owned distance, sample distance, crossing */
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
static inline double mk_dmax_cos_(mk_dv3_ a, mk_dv3_ b, mk_dv3_ c)
{
    const mk_dv3_ e0 = mk_dsub_(b, a), e1 = mk_dsub_(c, b), e2 = mk_dsub_(a, c);
    const double l0 = mk_dlen_(e0), l1 = mk_dlen_(e1), l2 = mk_dlen_(e2);
    if (!(l0 > 0) || !(l1 > 0) || !(l2 > 0)) return 1;
    const double ca = -mk_ddot_(e0, e2) / (l0 * l2), cb = -mk_ddot_(e0, e1) / (l0 * l1), cc = -mk_ddot_(e1, e2) / (l1 * l2);
    return fmax(ca, fmax(cb, cc));
}
/* an edge pq through triangle abc's interior (mk_seg_tri_strict) */
static inline int mk_dec_segtri_(mk_dv3_ p, mk_dv3_ q, mk_dv3_ a, mk_dv3_ b, mk_dv3_ c)
{
    const double P[3] = {p.x, p.y, p.z}, Q[3] = {q.x, q.y, q.z}, A[3] = {a.x, a.y, a.z}, B[3] = {b.x, b.y, b.z}, C[3] = {c.x, c.y, c.z};
    return mk_seg_tri_strict(P, Q, A, B, C);
}

typedef struct { double c; int u, v, s, k; } mk_dec_ent_;
typedef struct {
    const double *P;
    int nv, nt, nreg, relax, oom, gh;
    int *F, *OT, **VT, *NVT, *CVT, *STAMP;
    unsigned *MARK, MARKS, *TMARK, TMARKS; /* unsigned, so the marks wrap rather than overflow */
    double *N0, *Q;
    unsigned char *deadv, *deadt, *bnd;
    const unsigned char *keep;
    const mk_dec_opts *o;
    double dist2, maxedge2, cos2, cmin, mindot, gc, gs, samp, ftol, relaxlen, minalt, maxown2;
    mk_bvh ob;
    long ref[8], nrelax;
    int (*lt)[3];
    mk_dv3_ *ln;
    double (*nc)[9];
    int *ntb;
    int *ghead, *gnext, *gtri, *gver, *tver, gn, gcap;
    unsigned *seen, seenmark;
    double glo[3];
    int *ohead, *onext, *pend_p, *pend_t, pend_n, pend_c, pend_u, pend_v;
    mk_dec_ent_ *heap;
    size_t nh, ch;
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
static inline void mk_dec_fan_add_(mk_dec_ *c, int v, int t)
{
    if (c->NVT[v] == c->CVT[v]) {
        const int nc = c->CVT[v] ? 2 * c->CVT[v] : 8;
        int *q = mk_dec_grow_(c, c->VT[v], (size_t)nc, sizeof *q);
        if (c->oom) return;
        c->VT[v] = q, c->CVT[v] = nc;
    }
    c->VT[v][c->NVT[v]++] = t;
}
static inline void mk_dec_fan_del_(mk_dec_ *c, int v, int t)
{
    for (int i = 0; i < c->NVT[v]; i++)
        if (c->VT[v][i] == t) { c->VT[v][i] = c->VT[v][--c->NVT[v]]; return; }
}
static inline int mk_dec_nbrs_(mk_dec_ *c, int u, int *nb)
{
    int n = 0;
    c->MARKS += 2; /* past mk_dec_cheap_'s MARKS + 1, so its common vertices do not read as seen */
    for (int i = 0; i < c->NVT[u]; i++)
        for (int k = 0; k < 3; k++) {
            const int x = mk_dec_f_(c, c->VT[u][i])[k];
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
/* the fold lengths among the edges of fan(u) and fan(v): out[r][0] over 90 degrees, out[r][1] over fold; after: u moved
   onto v */
static inline void mk_dec_foldsum_(mk_dec_ *c, int u, int v, int after, double out[MK_DEC_REGIONS][2])
{
    int n = 0;
    c->TMARKS++;
    for (int r = 0; r < c->nreg; r++) out[r][0] = out[r][1] = 0;
    for (int pass = 0; pass < 2; pass++) {
        const int w = pass ? v : u;
        for (int i = 0; i < c->NVT[w]; i++) {
            const int t = c->VT[w][i];
            if (c->TMARK[t] == c->TMARKS) continue;
            c->TMARK[t] = c->TMARKS;
            if (after && mk_dec_has_(c, t, u) && mk_dec_has_(c, t, v)) continue;
            for (int k = 0; k < 3; k++) c->lt[n][k] = after && mk_dec_f_(c, t)[k] == u ? v : mk_dec_f_(c, t)[k];
            c->ln[n] = mk_dnorm_(mk_dtri_n_(mk_dec_pt_(c, c->lt[n][0]), mk_dec_pt_(c, c->lt[n][1]), mk_dec_pt_(c, c->lt[n][2])));
            n++;
        }
    }
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++) {
            const int a = c->lt[i][k], b = c->lt[i][(k + 1) % 3];
            int j = -1;
            mk_dv3_ nj;
            for (int q = 0; q < n; q++)
                if (q != i && mk_dec_has3_(c->lt[q], a) && mk_dec_has3_(c->lt[q], b)) { j = q; break; }
            if (j >= 0) {
                if (j < i) continue;
                nj = c->ln[j];
            } else {
                int t2 = -1;
                for (int q = 0; q < c->NVT[a]; q++) {
                    const int t = c->VT[a][q];
                    if (c->TMARK[t] != c->TMARKS && mk_dec_has_(c, t, b)) { t2 = t; break; }
                }
                if (t2 < 0) continue;
                const int *g = mk_dec_f_(c, t2);
                nj = mk_dnorm_(mk_dtri_n_(mk_dec_pt_(c, g[0]), mk_dec_pt_(c, g[1]), mk_dec_pt_(c, g[2])));
            }
            const double d = mk_ddot_(c->ln[i], nj);
            if (d >= c->cos2 && d >= 0) continue;
            const mk_dv3_ pa = mk_dec_pt_(c, a), pb = mk_dec_pt_(c, b);
            const double L = mk_dlen_(mk_dsub_(pa, pb));
            const int r = mk_dec_region_(c, mk_dmul_(mk_dadd_(pa, pb), 0.5));
            if (d < 0) out[r][0] += L;
            if (d < c->cos2) out[r][1] += L;
        }
}
/* the crossing grid: an insert-only hash of cells (side gs) whose entries carry the triangle's version. A triangle
   whose box covers more than MK_DEC_BIG cells is listed once, in cell (-1, -1, -1), which every query reads; such a new
   triangle reads every live one instead of its cells, so one large triangle among small ones does not fill memory with
   its box's cells. */
#define MK_DEC_BIG 4096
static inline int mk_dec_cells_(const mk_dec_ *c, const mk_dv3_ *d, int lo[3], int hi[3])
{
    const double mn[3] = {fmin(d[0].x, fmin(d[1].x, d[2].x)), fmin(d[0].y, fmin(d[1].y, d[2].y)), fmin(d[0].z, fmin(d[1].z, d[2].z))};
    const double mx[3] = {fmax(d[0].x, fmax(d[1].x, d[2].x)), fmax(d[0].y, fmax(d[1].y, d[2].y)), fmax(d[0].z, fmax(d[1].z, d[2].z))};
    for (int k = 0; k < 3; k++) lo[k] = (int)floor((mn[k] - c->glo[k]) / c->gs), hi[k] = (int)floor((mx[k] - c->glo[k]) / c->gs);
#if defined(MK_MUTATE) && MK_MUTATE == 227
    return 0; /* 227: every triangle entered in every cell of its box */
#endif
    return (double)(hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1) > MK_DEC_BIG;
}
static inline size_t mk_dec_gkey_(const mk_dec_ *c, int x, int y, int z)
{
    return (size_t)((((uint64_t)(uint32_t)x * 73856093u) ^ ((uint64_t)(uint32_t)y * 19349663u) ^ ((uint64_t)(uint32_t)z * 83492791u)) &
                    (uint64_t)(c->gh - 1));
}
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
                    if (c->gcap > INT_MAX / 2) { c->oom = 1; return; }
                    const int nc = c->gcap ? 2 * c->gcap : 1 << 16;
                    c->gnext = mk_dec_grow_(c, c->gnext, (size_t)nc, sizeof *c->gnext);
                    c->gtri = mk_dec_grow_(c, c->gtri, (size_t)nc, sizeof *c->gtri);
                    c->gver = mk_dec_grow_(c, c->gver, (size_t)nc, sizeof *c->gver);
                    if (c->oom) return;
                    c->gcap = nc;
                }
                const size_t k = mk_dec_gkey_(c, x, y, z);
                c->gtri[c->gn] = t, c->gver[c->gn] = c->tver[t], c->gnext[c->gn] = c->ghead[k], c->ghead[k] = c->gn++;
            }
}
static inline void mk_dec_grebuild_(mk_dec_ *c)
{
    for (int i = 0; i < c->gh; i++) c->ghead[i] = -1;
    c->gn = 0;
    for (int t = 0; t < c->nt; t++)
        if (!c->deadt[t]) mk_dec_gins_(c, t);
}
/* new triangle d (its corners g) against live triangle t sharing none of its corners (nor skip): 1 when they cross */
static inline int mk_dec_xone_(mk_dec_ *c, const mk_dv3_ *d, const int *g, int skip, int t)
{
    if (c->deadt[t] || c->seen[t] == c->seenmark) return 0;
    c->seen[t] = c->seenmark;
    if (mk_dec_has_(c, t, skip) || mk_dec_has_(c, t, g[0]) || mk_dec_has_(c, t, g[1]) || mk_dec_has_(c, t, g[2])) return 0;
    const int *f = mk_dec_f_(c, t);
    const mk_dv3_ tt[3] = {mk_dec_pt_(c, f[0]), mk_dec_pt_(c, f[1]), mk_dec_pt_(c, f[2])};
#if defined(MK_MUTATE) && MK_MUTATE == 163
    return 0; /* 163: the crossing test left out */
#endif
    for (int k = 0; k < 3; k++)
        if (mk_dec_segtri_(d[k], d[(k + 1) % 3], tt[0], tt[1], tt[2]) || mk_dec_segtri_(tt[k], tt[(k + 1) % 3], d[0], d[1], d[2])) return 1;
    return 0;
}
/* new triangle d (its corners g) crosses a live triangle sharing none of its corners (nor skip) */
static inline int mk_dec_crosses_(mk_dec_ *c, const mk_dv3_ *d, const int *g, int skip)
{
    int lo[3], hi[3];
    c->seenmark++;
    if (mk_dec_cells_(c, d, lo, hi)) {
        for (int t = 0; t < c->nt; t++)
            if (mk_dec_xone_(c, d, g, skip, t)) return 1;
        return 0;
    }
    for (int e = c->ghead[mk_dec_gkey_(c, -1, -1, -1)]; e >= 0; e = c->gnext[e])
        if (c->gver[e] == c->tver[c->gtri[e]] && mk_dec_xone_(c, d, g, skip, c->gtri[e])) return 1;
    for (int x = lo[0]; x <= hi[0]; x++)
        for (int y = lo[1]; y <= hi[1]; y++)
            for (int z = lo[2]; z <= hi[2]; z++)
                for (int e = c->ghead[mk_dec_gkey_(c, x, y, z)]; e >= 0; e = c->gnext[e])
                    if (c->gver[e] == c->tver[c->gtri[e]] && mk_dec_xone_(c, d, g, skip, c->gtri[e])) return 1;
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
            const int t = c->VT[u][i], *f = mk_dec_f_(c, t);
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
        const int t = c->VT[u][i], *f = mk_dec_f_(c, t);
        if (!mk_dec_has_(c, t, v)) continue;
        if (no < 2) opp[no] = f[0] != u && f[0] != v ? f[0] : f[1] != u && f[1] != v ? f[1] : f[2];
        no++;
    }
    const int need = c->bnd[u] ? 1 : 2; /* a rim vertex goes only along the rim, onto a rim neighbour over a rim edge */
    if ((c->bnd[u] && !c->bnd[v]) || no != need || (need == 2 && opp[0] == opp[1])) return 0;
    c->MARKS += 2;
    for (int i = 0; i < c->NVT[u]; i++)
        for (int k = 0; k < 3; k++) c->MARK[mk_dec_f_(c, c->VT[u][i])[k]] = c->MARKS;
    int common = 0;
    for (int i = 0; i < c->NVT[v]; i++)
        for (int k = 0; k < 3; k++) {
            const int x = mk_dec_f_(c, c->VT[v][i])[k];
            if (x != u && x != v && c->MARK[x] == c->MARKS) c->MARK[x] = c->MARKS + 1, common++;
        }
#if !defined(MK_MUTATE) || MK_MUTATE != 164
    if (common != need) { c->ref[0]++; return 0; }
    int tu = 0, tv = 0;
    if (need == 2)
        for (int i = 0; i < c->NVT[u]; i++) tu |= mk_dec_has_(c, c->VT[u][i], opp[0]) && mk_dec_has_(c, c->VT[u][i], opp[1]) && !mk_dec_has_(c, c->VT[u][i], v);
    if (need == 2)
        for (int i = 0; i < c->NVT[v]; i++) tv |= mk_dec_has_(c, c->VT[v][i], opp[0]) && mk_dec_has_(c, c->VT[v][i], opp[1]) && !mk_dec_has_(c, c->VT[v][i], u);
    if (tu && tv) { c->ref[0]++; return 0; }
#else
    (void)common;
#endif
    double oldc = 0;
    for (int i = 0; i < c->NVT[u]; i++) {
        const int *q = mk_dec_f_(c, c->VT[u][i]);
        oldc = fmax(oldc, mk_dmax_cos_(mk_dec_pt_(c, q[0]), mk_dec_pt_(c, q[1]), mk_dec_pt_(c, q[2])));
    }
    const mk_dv3_ pv = mk_dec_pt_(c, v);
    for (int i = 0; i < c->NVT[u]; i++) {
        const int t = c->VT[u][i], *f = mk_dec_f_(c, t);
        if (mk_dec_has_(c, t, v)) continue;
        mk_dv3_ p0[3], p1[3];
        for (int k = 0; k < 3; k++) p0[k] = mk_dec_pt_(c, f[k]), p1[k] = f[k] == u ? pv : p0[k];
        const mk_dv3_ n0 = mk_dtri_n_(p0[0], p0[1], p0[2]), n1 = mk_dtri_n_(p1[0], p1[1], p1[2]);
        const double l0 = mk_dlen_(n0), l1 = mk_dlen_(n1);
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
    }
    if (c->relax) return 1;
#if !defined(MK_MUTATE) || MK_MUTATE != 161
    double b[MK_DEC_REGIONS][2], a[MK_DEC_REGIONS][2];
    mk_dec_foldsum_(c, u, v, 0, b), mk_dec_foldsum_(c, u, v, 1, a);
    for (int r = 0; r < c->nreg; r++)
        if (a[r][0] > b[r][0] + c->ftol || a[r][1] > b[r][1] + c->ftol) { c->ref[4]++; return 0; }
#endif
    return 1;
}
/* distance and crossing; a pass leaves the owned vertices' new triangles pending for mk_dec_collapse_ */
static inline int mk_dec_full_(mk_dec_ *c, int u, int v)
{
    int m = 0;
    const mk_dv3_ pv = mk_dec_pt_(c, v);
    c->pend_u = c->pend_v = -1, c->pend_n = 0;
    for (int i = 0; i < c->NVT[u]; i++) {
        const int t = c->VT[u][i], *f = mk_dec_f_(c, t);
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
            for (int p = i < 0 ? u : c->ohead[c->VT[u][i]]; p >= 0; p = i < 0 ? -1 : c->onext[p]) {
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
                    const int t2 = c->VT[v][i2], *g2 = mk_dec_f_(c, t2);
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
        for (int j = 0; j < m; j++) {
            double lm = 0;
            for (int k = 0; k < 3; k++) {
                double e = 0;
                for (int q = 0; q < 3; q++) {
                    const double d = c->nc[j][3 * k + q] - c->nc[j][3 * ((k + 1) % 3) + q];
                    e += d * d;
                }
                lm = fmax(lm, e);
            }
            const double nn = ceil(sqrt(lm) / c->samp);
            const int n = !(nn >= 2) ? 2 : nn > 4 ? 4 : (int)nn;
            for (int a = 0; a <= n; a++)
                for (int b = 0; a + b <= n; b++) {
                    const int cc = n - a - b;
                    if (a == n || b == n || cc == n) continue;
                    double p[3], q[3];
                    int tr = -1;
                    for (int k = 0; k < 3; k++) p[k] = (a * c->nc[j][k] + b * c->nc[j][3 + k] + cc * c->nc[j][6 + k]) / n;
                    mk_bvh_nearest(&c->ob, p, c->dist2, &tr, q);
                    if (tr < 0 && !c->relax) { c->ref[6]++; return 0; }
                }
        }
    }
#endif
    MK_DEC_OWN_EARLY_
    for (int j = 0; j < m; j++) {
        const int t = c->ntb[j], *f = mk_dec_f_(c, t);
        mk_dv3_ d[3];
        int g[3];
        for (int k = 0; k < 3; k++) d[k] = mk_dv_(c->nc[j][3 * k], c->nc[j][3 * k + 1], c->nc[j][3 * k + 2]), g[k] = f[k] == u ? v : f[k];
        if (mk_dec_crosses_(c, d, g, u)) { c->ref[7]++; return 0; }
    }
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
static inline void mk_dec_push_(mk_dec_ *c, mk_dec_ent_ e)
{
    if (c->nh == c->ch) {
        const size_t nc = c->ch ? 2 * c->ch : 1 << 12;
        mk_dec_ent_ *h = mk_dec_grow_(c, c->heap, nc, sizeof *h);
        if (c->oom) return;
        c->heap = h, c->ch = nc;
    }
    size_t i = c->nh++;
    for (; i && c->heap[(i - 1) / 2].c > e.c; i = (i - 1) / 2) c->heap[i] = c->heap[(i - 1) / 2];
    c->heap[i] = e;
}
static inline int mk_dec_pop_(mk_dec_ *c, mk_dec_ent_ *out)
{
    if (!c->nh) return 0;
    *out = c->heap[0];
    const mk_dec_ent_ e = c->heap[--c->nh];
    size_t i = 0;
    for (;;) {
        size_t k = 2 * i + 1;
        if (k >= c->nh) break;
        if (k + 1 < c->nh && c->heap[k + 1].c < c->heap[k].c) k++;
        if (c->heap[k].c >= e.c) break;
        c->heap[i] = c->heap[k], i = k;
    }
    if (c->nh) c->heap[i] = e;
    return 1;
}
/* queue u's (skip+1)-th cheapest target that passes the cheap refusals */
static inline void mk_dec_rate_(mk_dec_ *c, int u, int skip)
{
    c->STAMP[u]++;
    if (c->deadv[u]) return;
    int nb[MK_DEC_MAXN];
    double cs[MK_DEC_MAXN];
    const int n = mk_dec_nbrs_(c, u, nb);
    for (int i = 0; i < n; i++) cs[i] = mk_dec_cost_(c, u, nb[i]);
    for (int got = 0;;) {
        int b = -1;
        for (int i = 0; i < n; i++)
            if (nb[i] >= 0 && (b < 0 || cs[i] < cs[b])) b = i;
        if (b < 0) return;
        const int v = nb[b];
        nb[b] = -1;
        if (!mk_dec_cheap_(c, u, v)) continue;
        if (got++ == skip) { mk_dec_push_(c, (mk_dec_ent_){cs[b], u, v, c->STAMP[u], skip}); return; }
    }
}
static inline void mk_dec_collapse_(mk_dec_ *c, int u, int v)
{
    int fan[MK_DEC_MAXN];
    const int n = c->NVT[u] < MK_DEC_MAXN ? c->NVT[u] : MK_DEC_MAXN;
    memcpy(fan, c->VT[u], sizeof *fan * (size_t)n);
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
    c->NVT[u] = 0, c->deadv[u] = 1;
    for (int k = 0; k < 10; k++) c->Q[10 * (size_t)v + k] += c->Q[10 * (size_t)u + k];
    for (int i = 0; i < c->NVT[v]; i++) c->tver[c->VT[v][i]]++, mk_dec_gins_(c, c->VT[v][i]);
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
                const int t2 = c->VT[x][q];
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
static inline int mk_dec_bad_(mk_dec_ *c, int t, double minalt)
{
    const int *f = mk_dec_f_(c, t);
    const mk_dv3_ d[3] = {mk_dec_pt_(c, f[0]), mk_dec_pt_(c, f[1]), mk_dec_pt_(c, f[2])};
    if (minalt > 0) {
        const double a2 = mk_dlen_(mk_dtri_n_(d[0], d[1], d[2])), l = fmax(mk_dlen_(mk_dsub_(d[0], d[1])), fmax(mk_dlen_(mk_dsub_(d[1], d[2])), mk_dlen_(mk_dsub_(d[2], d[0]))));
        if (!(l > 0) || a2 / l < minalt) return 1;
    }
    return mk_dec_crosses_(c, d, f, -1);
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
    if (c->VT)
        for (int i = 0; i < c->nv; i++) free(c->VT[i]);
    free(c->VT), free(c->F), free(c->OT), free(c->NVT), free(c->CVT), free(c->STAMP), free(c->MARK), free(c->TMARK), free(c->N0), free(c->Q);
    free(c->deadv), free(c->deadt), free(c->bnd), free(c->lt), free(c->ln), free(c->nc), free(c->ntb), free(c->ghead), free(c->gnext);
    free(c->gtri), free(c->gver), free(c->tver), free(c->seen), free(c->ohead), free(c->onext), free(c->pend_p), free(c->pend_t), free(c->heap);
    mk_bvh_free(&c->ob);
}

/* Decimates (V[3 nv], T[3 nt]) under o (NULL: mk_dec_defaults with mm = 1e-3, metres). Out, in the caller's arrays:
   alive[nv] (1 kept, 0 removed) and the result's *nto faces in To (capacity 3 nt ints) in the input's vertex numbers;
   *rep if rep is not NULL. Returns 0, or -1 for a bad mesh (a size under 1 or over INT_MAX / 4, an index out of range,
   a face with a repeated corner, a coordinate not finite, or coordinates past 2^30 grid cells), bad options (a negative
   or NaN bound, turn outside [-1, 1], minangle or fold outside [0, 180), nregion outside 1 .. MK_DEC_REGIONS or above 1
   without region, target < 0), or no memory; alive and To are then untouched and *rep zeroed. */
static inline int mk_decimate(const double *V, int nv, const int *T, int nt, const mk_dec_opts *o, unsigned char *alive, int *To, int *nto,
                              mk_dec_report *rep)
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
    for (int k = 0; k < 3; k++) {
        c.glo[k] = lo[k];
        if (!((hi[k] - lo[k]) / c.gc < 0x1p30)) return -1;
    }
    int ret = -1;
    c.F = malloc(sizeof *c.F * 3 * (size_t)nt), c.OT = malloc(sizeof *c.OT * 3 * (size_t)nt), c.VT = calloc((size_t)nv, sizeof *c.VT);
    c.NVT = calloc((size_t)nv, sizeof *c.NVT), c.CVT = calloc((size_t)nv, sizeof *c.CVT), c.STAMP = calloc((size_t)nv, sizeof *c.STAMP);
    c.MARK = calloc((size_t)nv, sizeof *c.MARK), c.TMARK = calloc((size_t)nt, sizeof *c.TMARK), c.N0 = calloc(3 * (size_t)nv, sizeof *c.N0);
    c.Q = calloc(10 * (size_t)nv, sizeof *c.Q), c.deadv = calloc((size_t)nv, 1), c.deadt = calloc((size_t)nt, 1), c.bnd = calloc((size_t)nv, 1);
    c.lt = malloc(sizeof *c.lt * 4 * MK_DEC_MAXN), c.ln = malloc(sizeof *c.ln * 4 * MK_DEC_MAXN), c.nc = malloc(sizeof *c.nc * MK_DEC_MAXN);
    c.ntb = malloc(sizeof *c.ntb * MK_DEC_MAXN), c.tver = calloc((size_t)nt, sizeof *c.tver), c.seen = calloc((size_t)nt, sizeof *c.seen);
    c.ohead = malloc(sizeof *c.ohead * (size_t)nt), c.onext = malloc(sizeof *c.onext * (size_t)nv);
    c.gh = 1024;
    while (c.gh < 2 * nt && c.gh < (1 << 26)) c.gh *= 2;
    c.ghead = malloc(sizeof *c.ghead * (size_t)c.gh);
    if (!c.F || !c.OT || !c.VT || !c.NVT || !c.CVT || !c.STAMP || !c.MARK || !c.TMARK || !c.N0 || !c.Q || !c.deadv || !c.deadt || !c.bnd || !c.lt ||
        !c.ln || !c.nc || !c.ntb || !c.tver || !c.seen || !c.ohead || !c.onext || !c.ghead)
        goto out;
    memcpy(c.F, T, sizeof *c.F * 3 * (size_t)nt), memcpy(c.OT, T, sizeof *c.OT * 3 * (size_t)nt);
    if (c.dist2 > 0 && mk_bvh_build(&c.ob, V, c.OT, nt)) goto out;
    for (int t = 0; t < nt; t++) c.ohead[t] = -1;
    for (int i = 0; i < nv; i++) c.onext[i] = -1;
    for (int t = 0; t < nt; t++) {
        const int *f = mk_dec_f_(&c, t);
        for (int q = 0; q < 3; q++) mk_dec_fan_add_(&c, f[q], t);
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
                const int x = mk_dec_f_(&c, c.VT[u][i])[k];
                if (x == u) continue;
                int n = 0;
                for (int j = 0; j < c.NVT[u]; j++) n += mk_dec_has_(&c, c.VT[u][j], x);
                if (n == 1) c.bnd[u] = 1;
            }
        nrim += c.bnd[u];
    }
    mk_dec_grebuild_(&c);
    if (c.oom) goto out;
    mk_dec_report r = {0};
    r.rim0 = nrim;
    mk_dec_total_(&c, r.fold0, r.area0);
    long nalive = nv, done = 0;
    if (o->clean) r.clean_before = mk_dec_cleanup_(&c, 0, &nalive, &r.left_before);
    for (int u = 0; u < nv && !c.oom; u++) mk_dec_rate_(&c, u, 0);
    mk_dec_ent_ e;
    while (nalive > o->target && !c.oom && mk_dec_pop_(&c, &e)) {
        if (e.s != c.STAMP[e.u] || c.deadv[e.u] || c.deadv[e.v]) continue;
        if (!mk_dec_cheap_(&c, e.u, e.v)) { mk_dec_rate_(&c, e.u, 0); continue; }
        if (!mk_dec_full_(&c, e.u, e.v)) { mk_dec_rate_(&c, e.u, e.k + 1); continue; }
        int nb[MK_DEC_MAXN];
        mk_dec_collapse_(&c, e.u, e.v), nalive--;
        if (++done % 200000 == 0) mk_dec_grebuild_(&c);
        const int m = mk_dec_nbrs_(&c, e.v, nb);
        mk_dec_rate_(&c, e.v, 0);
        for (int i = 0; i < m; i++) mk_dec_rate_(&c, nb[i], 0);
    }
    if (o->clean && !c.oom) r.clean_after = mk_dec_cleanup_(&c, o->minalt, &nalive, &r.left_after);
    if (c.oom) goto out;
    mk_dec_total_(&c, r.fold1, r.area1);
    int m = 0, kv = 0;
    for (int i = 0; i < nv; i++) alive[i] = !c.deadv[i], kv += !c.deadv[i], r.rim1 += c.bnd[i] && !c.deadv[i];
    for (int t = 0; t < nt; t++)
        if (!c.deadt[t]) memcpy(To + 3 * (size_t)m++, mk_dec_f_(&c, t), 3 * sizeof *To);
    *nto = m;
    r.nv = kv, r.nt = m, r.maxown = sqrt(c.maxown2), r.relaxed = c.nrelax;
    memcpy(r.refused, c.ref, sizeof r.refused);
    if (rep) *rep = r;
    ret = 0;
out:
    mk_dec_free_(&c);
    return ret;
}
#endif
