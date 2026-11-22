/* mathkit: a closed mesh changed only near a seam, the rest of it kept triangle for triangle. Header only, C23.
   A field's surface, by marching tetrahedra on a box of nodes, is cut at a seam function's zero and stitched into the
   mesh cut at the same zero: its rim is put into the mesh's rim and welded. The field is the caller's, on mk_sdf's grid
   layout (a grid from mk_sdf_mesh fits, and so does any field built from one; mk_sdf_mesh_box grids only the seam's
   neighbourhood).
   mk_mtet: marching tetrahedra, Freudenthal's six tetrahedra per cube round its diagonal, the same in every cube, so
   neighbours' faces match and the surface is a closed manifold inside the box whatever the field, open only on the
   box's faces (surface nets, mk_snet, can touch themselves where two sheets share a cell); a vertex on each node pair
   of unlike sign at the linear crossing, held between 0.1 and 0.9 of the pair, so no triangle has |2A| under 0.01 h^2
   (no needles from a crossing at a node); wound outward, the field negative inside. On a distance field the surface
   sinks by about h^2 H / 6 on average (H the mean curvature): a ball's volume is low by (h/R)^2 / 2, a tube's by
   (h/r)^2 / 6 (measured for R/h 5 to 40). A part with no solid or gap thinner than 2 h keeps its components and Euler
   characteristic; a thin edge, sheet or narrow gap anywhere on it can add specks or handles. On 46 anatomical meshes
   at h 0.5 mm, each on its own mk_sdf_mesh field: all closed and manifold, |2A| at least 0.01 h^2; the 24 deeper than
   6 mm within -0.22 % to +0.04 % of their volume (15 within 0.1 %), 17 of them with components and Euler
   characteristic kept (the liver gained 22 specks and voids, 0.09 ml in all); the 12 thinner than 2 mm fell apart (-20
   to -61 % of their volume; a membrane into 924 pieces), which this does not cover.
   A graft: mk_graft_init(S, psi) takes the untouched closed mesh S and the seam function psi (negative where the patch
   goes, zero at the seam): psi at S's vertices (nudged off zero by `nudge`, so the zero line misses them), linear over
   S's triangles. mk_graft_refine can first split S's long edges near the seam. mk_graft_patch(field) is the field's
   surface kept where psi < 0, cut on its edges at psi = 0 (bisected on the patch's psi: S's linear psi at the closest
   point of S within `near` of it, psi itself farther off). The caller may thin the patch (mk_decimate on an open
   surface keeps its rim on its rim), or stitch it as it is (t_graft's dent, unthinned, on S of 5 and 3 subdivisions:
   closed with no crossing pair).
   mk_graft_stitch(patch) cuts S at psi = 0 (S kept where psi > 0, its crossing triangles cut on their edges, a vertex
   just off the seam put on its neighbours' side by mk_graft_init), moves each patch rim vertex onto S's cut rim
   (no more than `snap`, or it refuses; one whose closest point lies within 1 % of an S rim edge's length of its end
   goes onto that S vertex, the extra move not counted against snap or in rep->snap) and the patch's inside near the
   rim with it (by (1 - r/R)^2 of the rim's move, R three times the move, so the first rings do not fold over the rim),
   pairs the rims as loops (each patch rim loop on one of S's, each of S's met by one, else counted), merges the patch
   rim vertices whose places on S's rim fall out of order into a rim neighbour, puts each rim's vertices into the
   other's rim edges (the rim triangles fanned, S's in their own planes), welds, and keeps the piece of most volume,
   either sign (closed specks the cut left off are dropped and counted; pieces wound against S are reported). It
   returns 0 only for one closed, consistently wound piece with S's orientation, no vertex pinched, no new triangle of
   zero area and no new crossing pair; every other outcome has a code and a reason (mk_graft_stitch's comment).
   out->keep marks for mk_decimate: 1 S's own vertices, 3 its cut points, 2 the patch's rim, 0 the rest.
   Single-threaded. */
#ifndef MK_GRAFT_H
#define MK_GRAFT_H
#include "mk_geom.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* a mesh handed back: V[3 nv], T[3 nt], keep[nv] (may be NULL), all malloc'd; free with mk_gmesh_free */
typedef struct { double *V; int *T; int nv, nt; unsigned char *keep; } mk_gmesh;
static inline void mk_gmesh_free(mk_gmesh *m) { free(m->V), free(m->T), free(m->keep), memset(m, 0, sizeof *m); }

/* a growing surface with rim edges (a, b, its triangle) */
typedef struct { double *V; int nv, cv, *T, nt, ct, *E, ne, ce, *map, nkeep, oom; } mk_gcut_;
static inline void mk_gcut_free_(mk_gcut_ *c) { free(c->V), free(c->T), free(c->E), free(c->map), memset(c, 0, sizeof *c); }
static inline int mk_gcut_v_(mk_gcut_ *c, const double *p)
{
    if (c->nv == c->cv) {
        const int n = c->cv ? (c->cv > INT_MAX / 2 ? -1 : 2 * c->cv) : 1024;
        double *q = n > 0 ? realloc(c->V, 3 * (size_t)n * sizeof *q) : NULL;
        if (!q) { c->oom = 1; return 0; }
        c->V = q, c->cv = n;
    }
    memcpy(c->V + 3 * (size_t)c->nv, p, 3 * sizeof *p);
    return c->nv++;
}
static inline int mk_gcut_t_(mk_gcut_ *c, int a, int b, int d)
{
    if (c->nt == c->ct) {
        const int n = c->ct ? (c->ct > INT_MAX / 2 ? -1 : 2 * c->ct) : 1024;
        int *q = n > 0 ? realloc(c->T, 3 * (size_t)n * sizeof *q) : NULL;
        if (!q) { c->oom = 1; return 0; }
        c->T = q, c->ct = n;
    }
    c->T[3 * (size_t)c->nt] = a, c->T[3 * (size_t)c->nt + 1] = b, c->T[3 * (size_t)c->nt + 2] = d;
    return c->nt++;
}
static inline void mk_gcut_e_(mk_gcut_ *c, int a, int b, int t)
{
    if (c->ne == c->ce) {
        const int n = c->ce ? (c->ce > INT_MAX / 2 ? -1 : 2 * c->ce) : 1024;
        int *q = n > 0 ? realloc(c->E, 3 * (size_t)n * sizeof *q) : NULL;
        if (!q) { c->oom = 1; return; }
        c->E = q, c->ce = n;
    }
    c->E[3 * (size_t)c->ne] = a, c->E[3 * (size_t)c->ne + 1] = b, c->E[3 * (size_t)c->ne + 2] = t, c->ne++;
}
/* an open-addressing map from a pair of ints to a vertex */
typedef struct { uint64_t *k; int *v; size_t n; } mk_gmap_;
static inline int mk_gmap_init_(mk_gmap_ *m, size_t want)
{
    m->n = 1024;
    while (m->n < 2 * want + 16 && m->n < ((size_t)1 << 40)) m->n *= 2;
    m->k = malloc(m->n * sizeof *m->k), m->v = malloc(m->n * sizeof *m->v);
#if defined(MK_MUTATE) && MK_MUTATE == 233
    if (!m->k || !m->v) { free(m->k), free(m->v); return -1; } /* 233: the pointers left set */
#else
    /* cleared, because the callers free them again on failure */
    if (!m->k || !m->v) { free(m->k), free(m->v), m->k = NULL, m->v = NULL; return -1; }
#endif
    for (size_t i = 0; i < m->n; i++) m->k[i] = UINT64_MAX;
    return 0;
}
static inline size_t mk_gmap_slot_(const mk_gmap_ *m, uint64_t key)
{
    size_t h = (size_t)((key * 0x9E3779B97F4A7C15ull) >> 20) & (m->n - 1);
    while (m->k[h] != UINT64_MAX && m->k[h] != key) h = (h + 1) & (m->n - 1);
    return h;
}

#if defined(MK_MUTATE) && MK_MUTATE == 232
#define MK_GRAFT_BISECT_(t_, ta_, tb_) /* 232: no bisection, so a psi of +inf gives NaN cut points */
#else
#define MK_GRAFT_BISECT_(t_, ta_, tb_) if (!(t_ >= ta_ && t_ <= tb_)) t_ = 0.5 * (ta_ + tb_);
#endif
#if defined(MK_MUTATE) && MK_MUTATE == 367
#define MK_MT_HOLD_(t_) /* 367: a vertex anywhere on its node pair, at a node when the field is 0 there */
#elif defined(MK_MUTATE) && MK_MUTATE == 368
#define MK_MT_HOLD_(t_) t_ = 0.5; /* 368: every vertex at its node pair's middle */
#else
#define MK_MT_HOLD_(t_) t_ = !(t_ >= 0.1) ? 0.1 : t_ > 0.9 ? 0.9 : t_;
#endif
#if defined(MK_MUTATE) && MK_MUTATE == 363
#define MK_GRAFT_HOLD_(t_) /* 363: a cut point anywhere on its edge, so a vertex just off the seam leaves slivers */
#else
#define MK_GRAFT_HOLD_(t_) t_ = t_ < (g ? 0.1 : 0.01) ? (g ? 0.1 : 0.01) : t_ > (g ? 0.9 : 0.99) ? (g ? 0.9 : 0.99) : t_; /* NaN left as it is (232) */
#endif
/* The field F's surface by marching tetrahedra: nodes lo + h (i, j, k), index (k n1 + j) n0 + i (mk_sdf's layout),
   negative inside. Out: out->V, out->T malloc'd (closed and outward inside the box; keep NULL). Returns 0, or -1 for a
   size under 2, more nodes than fit, h not positive and finite, a non-finite lo or far corner lo + (n - 1) h, or no
   memory (out then empty). A NaN node counts as outside. */
static inline int mk_mtet(const float *F, const int n[3], const double lo[3], double h, mk_gmesh *out)
{
    memset(out, 0, sizeof *out);
#if defined(MK_MUTATE) && MK_MUTATE == 230
    /* 230: a finite lo's sum can overflow */
    if (!F || n[0] < 2 || n[1] < 2 || n[2] < 2 || !(h > 0) || !isfinite(h) || !isfinite(lo[0] + lo[1] + lo[2])) return -1;
#else
    if (!F || n[0] < 2 || n[1] < 2 || n[2] < 2 || !(h > 0) || !isfinite(h) || !isfinite(lo[0]) || !isfinite(lo[1]) || !isfinite(lo[2])) return -1;
#endif
    if ((size_t)n[0] > (size_t)INT_MAX / (size_t)n[1] / (size_t)n[2] || (size_t)n[0] * n[1] * n[2] > SIZE_MAX / 64) return -1;
#if !(defined(MK_MUTATE) && MK_MUTATE == 231) /* 231: the far corner not looked at */
    for (int a = 0; a < 3; a++)
        if (!isfinite(lo[a] + (n[a] - 1) * h)) return -1; /* the far corner: every vertex lies between lo and it */
#endif
    static const int PERM[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}, PSGN[6] = {1, -1, -1, 1, 1, -1};
    const long nx = n[0], nxy = (long)n[0] * n[1];
    long ncube = 0;
    mk_gcut_ c = {0};
    mk_gmap_ M = {0};
#define MK_MT_IN_(q) (F[q] < 0)
    for (int pass = 0; pass < 2; pass++) {
        if (pass && mk_gmap_init_(&M, 8 * (size_t)ncube + 512)) goto fail;
#if defined(MK_MUTATE) && MK_MUTATE == 373
        for (int k = 1; k + 1 < n[2]; k++) /* 373: the first layer of cubes left out: the surface open inside the box */
#else
        for (int k = 0; k + 1 < n[2]; k++)
#endif
            for (int j = 0; j + 1 < n[1]; j++)
                for (int i = 0; i + 1 < n[0]; i++) {
                    const long q0 = (long)k * nxy + (long)j * nx + i;
                    int m8 = 0;
                    for (int v = 0; v < 8; v++) m8 += MK_MT_IN_(q0 + (v & 1) + (v >> 1 & 1) * nx + (v >> 2) * nxy);
                    if (m8 == 0 || m8 == 8) continue;
                    if (!pass) { ncube++; continue; }
                    for (int p = 0; p < 6; p++) {
                        long q[4];
                        int d[3] = {0, 0, 0}, in[4], m = 0;
                        q[0] = q0, d[PERM[p][0]] = 1, q[1] = q0 + d[0] + d[1] * nx + d[2] * nxy;
                        d[PERM[p][1]] = 1, q[2] = q0 + d[0] + d[1] * nx + d[2] * nxy, q[3] = q0 + 1 + nx + nxy;
                        for (int v = 0; v < 4; v++) in[v] = MK_MT_IN_(q[v]), m += in[v];
                        if (m == 0 || m == 4) continue;
                        int ord[4], no = 0, inv = 0;
                        for (int v = 0; v < 4; v++)
                            if (in[v]) ord[no++] = v;
                        for (int v = 0; v < 4; v++)
                            if (!in[v]) ord[no++] = v;
                        for (int x = 0; x < 4; x++)
                            for (int y = x + 1; y < 4; y++)
                                if (!in[x] && in[y]) inv++;
                        const int s = PSGN[p] * (inv & 1 ? -1 : 1);
#if defined(MK_MUTATE) && MK_MUTATE == 181
                        const int sw = -s; /* 181: every triangle wound inward */
#else
                        const int sw = s;
#endif
                        const long A = q[ord[0]], B = q[ord[1]], C = q[ord[2]], D = q[ord[3]];
                        int vv[4];
#define MK_MT_V_(x, y, slot)                                                                                                                      \
    do {                                                                                                                                          \
        long a_ = (x), b_ = (y);                                                                                                                  \
        if (a_ > b_) { const long t_ = a_; a_ = b_, b_ = t_; }                                                                                    \
        const uint64_t key_ = (uint64_t)a_ << 32 | (uint64_t)b_;                                                                                  \
        const size_t hs_ = mk_gmap_slot_(&M, key_);                                                                                               \
        if (M.k[hs_] == key_) slot = M.v[hs_];                                                                                                    \
        else {                                                                                                                                    \
            const double fa_ = F[a_], fb_ = F[b_];                                                                                                \
            double t_ = fa_ / (fa_ - fb_);                                                                                                        \
            MK_MT_HOLD_(t_)                                                                                                                       \
            const long ia_[3] = {a_ % nx, a_ / nx % n[1], a_ / nxy}, ib_[3] = {b_ % nx, b_ / nx % n[1], b_ / nxy};                                \
            double p_[3];                                                                                                                         \
            for (int k_ = 0; k_ < 3; k_++) p_[k_] = lo[k_] + ((double)ia_[k_] + t_ * (double)(ib_[k_] - ia_[k_])) * h;                           \
            slot = mk_gcut_v_(&c, p_);                                                                                                            \
            M.k[hs_] = key_, M.v[hs_] = slot;                                                                                                     \
        }                                                                                                                                         \
    } while (0)
                        if (m == 1) {
                            MK_MT_V_(A, B, vv[0]); MK_MT_V_(A, C, vv[1]); MK_MT_V_(A, D, vv[2]);
                            if (sw > 0) mk_gcut_t_(&c, vv[0], vv[1], vv[2]);
                            else mk_gcut_t_(&c, vv[0], vv[2], vv[1]);
                        } else if (m == 3) {
                            MK_MT_V_(D, A, vv[0]); MK_MT_V_(D, B, vv[1]); MK_MT_V_(D, C, vv[2]);
                            if (sw < 0) mk_gcut_t_(&c, vv[0], vv[2], vv[1]);
                            else mk_gcut_t_(&c, vv[0], vv[1], vv[2]);
                        } else {
                            MK_MT_V_(A, C, vv[0]); MK_MT_V_(A, D, vv[1]); MK_MT_V_(B, D, vv[2]); MK_MT_V_(B, C, vv[3]);
                            if (sw > 0) mk_gcut_t_(&c, vv[0], vv[1], vv[2]), mk_gcut_t_(&c, vv[0], vv[2], vv[3]);
                            else mk_gcut_t_(&c, vv[0], vv[2], vv[1]), mk_gcut_t_(&c, vv[0], vv[3], vv[2]);
                        }
#undef MK_MT_V_
                        if (c.oom) goto fail;
                    }
                }
    }
#undef MK_MT_IN_
    free(M.k), free(M.v);
    out->V = c.V, out->T = c.T, out->nv = c.nv, out->nt = c.nt;
    return 0;
fail:
    free(M.k), free(M.v), mk_gcut_free_(&c);
    return -1;
}

typedef double (*mk_graft_psi)(const double *p, void *ctx);
typedef struct {
    const double *V;
    const int *T;
    int nv, nt;
    mk_bvh bv;
    double *psi; /* at S's vertices, nudged off zero */
    mk_graft_psi fn;
    void *ctx;
    double near; /* the patch takes S's linear psi within this of S */
} mk_graft;
typedef struct {
    int patch_nv, patch_nt, rim_s, rim_p; /* the patch kept; the rim edges of S's cut and of the patch */
    /* patch rim vertices put into S's rim edges, and landed on S's vertices; patch rim triangles fanned, dropped (both
       ends on one point), and split at their centre (two rim edges) */
    long inserted, welded, fanned, dropped, split;
    double snap;                         /* the largest move of a patch rim vertex onto S's rim */
    int pieces;                          /* closed pieces before the specks were dropped */
    double dropped_vol, vol0, vol1;      /* their volume; S's volume before and the result's after */
    long dup, lone;                      /* directed edges used twice; without a reverse (both 0 if closed manifold) */
    int loops_s, loops_p, unpaired;      /* rim loops of S's cut and of the patch; unpaired loops (-1: not loops) */
    long merged;                         /* patch rim vertices merged into a neighbour, being out of order on S's rim */
    int inverted;                        /* pieces wound against S (volume of the other sign), the kept one included */
    double inverted_vol;                 /* their volume, summed (of the other sign) */
    int ndropped;                       /* pieces dropped */
    double dropped_abs;                  /* their volumes' sizes summed (dropped_vol, signed, can cancel) */
    long zero_area;                      /* new triangles of zero area (twice the area under 1e-6 e^2, e S's mean edge: on
                                            a body mesh in metres with edges of 1 mm or more, at least as strict as meshcheck's 1e-12 m^2) */
    long cross0, cross1;                 /* crossing pairs with a changed triangle: S's (one not kept whole), the result's (one
                                            not S's own) */
    long pinched;                        /* vertices whose triangles make more than one fan */
    int reason;                          /* why it failed: the return code when over 0, an MK_GRAFT_E_ code with -1; 0 for 0 */
    unsigned faults;                     /* every fault found, a bit (1u << code) each: rc is only the first by precedence */
} mk_graft_report;
/* mk_graft_stitch's return codes over 0, in order of precedence (each with out given) */
#define MK_GRAFT_INVERTED 2 /* the kept piece wound against S, or pieces wound against S over 1 % of its volume */
#define MK_GRAFT_OPEN 1     /* not a closed manifold: an edge used twice or without its reverse, or a pinched vertex */
#define MK_GRAFT_SPLIT 3    /* pieces dropped holding over 1 % of S's volume: the cut split the mesh */
#define MK_GRAFT_ZERO 4     /* a new triangle of zero area */
#define MK_GRAFT_CROSS 5    /* more crossing pairs with a changed triangle than S had (rep->cross0, cross1) */
/* rep->reason with -1 (out empty) */
#define MK_GRAFT_E_ARGS 10   /* a bad argument */
#define MK_GRAFT_E_MEM 11    /* no memory */
#define MK_GRAFT_E_RIM 12    /* S's cut rim is not closed loops, or S's cut or the patch has no rim */
#define MK_GRAFT_E_SNAP 13   /* a patch rim vertex farther than snap from S's rim (rep->snap: its distance) */
#define MK_GRAFT_E_MERGE 14  /* a rim vertex out of order farther than snap from the neighbour it merges into (rep->snap) */
#define MK_GRAFT_E_FOLLOW 15 /* a patch rim edge over more than 64 of S's rim edges, or its triangle lost: the rims part */
#define MK_GRAFT_E_EMPTY 16  /* nothing of volume left */

static inline void mk_graft_free(mk_graft *g) { mk_bvh_free(&g->bv), free(g->psi), memset(g, 0, sizeof *g); }
typedef struct { long long k; int m; } mk_graft_ekey_;
static inline int mk_graft_ekcmp_(const void *a, const void *b)
{
    const long long x = ((const mk_graft_ekey_ *)a)->k, y = ((const mk_graft_ekey_ *)b)->k;
    return (x > y) - (x < y);
}
/* S refined where the seam will cut it, before mk_graft_init: every edge longer than lmax whose ends' psi reach below
   band and above -band is split at its midpoint, in both its triangles (1, 2 or 3 split edges a triangle by the usual
   templates, two split edges on the quad's shorter diagonal, so no T-junction), in passes until none is left or 12 have
   run. The surface is unchanged (midpoints of flat triangles' edges); only the seam's rim gets as fine as a patch's,
   which a coarse S's long rim edges are not. out: the refined mesh (V, T malloc'd, keep NULL; S's vertices first, as
   they were). Returns the passes run, or -1 for a bad mesh, a NULL fn, lmax or band not positive and finite, or no
   memory (out then empty). */
static inline int mk_graft_refine(const double *V, int nv, const int *T, int nt, mk_graft_psi fn, void *ctx, double lmax, double band, mk_gmesh *out)
{
    memset(out, 0, sizeof *out);
    if (!V || !T || !fn || nv < 3 || nt < 1 || nv > INT_MAX / 4 || nt > INT_MAX / 8 || !(lmax > 0) || !isfinite(lmax) || !(band > 0) || !isfinite(band)) return -1;
    for (size_t q = 0; q < 3 * (size_t)nt; q++)
        if (T[q] < 0 || T[q] >= nv) return -1;
#if !(defined(MK_MUTATE) && MK_MUTATE == 234) /* 234: coordinates not looked at */
    for (size_t q = 0; q < 3 * (size_t)nv; q++)
        if (!isfinite(V[q])) return -1; /* a bad mesh, as for mk_graft_init: NaN or inf would reach the midpoints */
#endif
    int mv = nv, mt = nt, pass = 0, *mT = malloc(3 * (size_t)nt * sizeof *mT);
    double *mV = malloc(3 * (size_t)nv * sizeof *mV), *ps = malloc((size_t)nv * sizeof *ps);
    mk_graft_ekey_ *E = NULL;
    if (!mT || !mV || !ps) goto fail;
    memcpy(mV, V, 3 * (size_t)nv * sizeof *mV), memcpy(mT, T, 3 * (size_t)nt * sizeof *mT);
    for (int i = 0; i < nv; i++) ps[i] = fn(V + 3 * (size_t)i, ctx);
    for (; pass < 12; pass++) {
        if (mt > INT_MAX / 16 || mv > INT_MAX / 4) goto fail;
        E = malloc(3 * (size_t)mt * sizeof *E);
        if (!E) goto fail;
        size_t ne = 0;
        for (int t = 0; t < mt; t++)
            for (int k = 0; k < 3; k++) {
                int a = mT[3 * (size_t)t + k], b = mT[3 * (size_t)t + (k + 1) % 3];
                if (a > b) { const int x = a; a = b, b = x; }
                const double *p = mV + 3 * (size_t)a, *q = mV + 3 * (size_t)b;
                const double L = sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]));
                if (L > lmax && fmin(ps[a], ps[b]) < band && fmax(ps[a], ps[b]) > -band) E[ne].k = (long long)a * mv + b, E[ne].m = -1, ne++;
            }
        if (!ne) break;
        qsort(E, ne, sizeof *E, mk_graft_ekcmp_);
        size_t u = 0;
        for (size_t i = 0; i < ne; i++)
            if (!u || E[i].k != E[u - 1].k) E[u++] = E[i];
        const int nv0 = mv;
        double *V2 = realloc(mV, 3 * ((size_t)mv + u) * sizeof *mV), *p2 = V2 ? realloc(ps, ((size_t)mv + u) * sizeof *ps) : NULL;
        if (V2) mV = V2;
        if (p2) ps = p2;
        int *T2 = V2 && p2 ? malloc(3 * 4 * (size_t)mt * sizeof *T2) : NULL;
        if (!T2) goto fail;
        for (size_t i = 0; i < u; i++) {
            const int a = (int)(E[i].k / nv0), b = (int)(E[i].k % nv0);
            for (int k = 0; k < 3; k++) mV[3 * (size_t)mv + k] = 0.5 * (mV[3 * (size_t)a + k] + mV[3 * (size_t)b + k]);
            ps[mv] = fn(mV + 3 * (size_t)mv, ctx), E[i].m = mv++;
        }
        int nt2 = 0;
#define MK_GR_TRI_(x, y, z) (T2[3 * (size_t)nt2] = (x), T2[3 * (size_t)nt2 + 1] = (y), T2[3 * (size_t)nt2 + 2] = (z), nt2++)
        for (int t = 0; t < mt; t++) {
            int v[3], mid[3], ns = 0;
            for (int k = 0; k < 3; k++) v[k] = mT[3 * (size_t)t + k];
            for (int k = 0; k < 3; k++) {
                int a = v[k], b = v[(k + 1) % 3];
                if (a > b) { const int x = a; a = b, b = x; }
                const mk_graft_ekey_ key = {(long long)a * nv0 + b, 0}, *f = bsearch(&key, E, u, sizeof *E, mk_graft_ekcmp_);
                mid[k] = f ? f->m : -1, ns += f != NULL;
            }
#if defined(MK_MUTATE) && MK_MUTATE == 209
            if (ns == 1) ns = 0; /* 209: a triangle with one split edge left whole: a T-junction */
#endif
            if (!ns) MK_GR_TRI_(v[0], v[1], v[2]);
            else if (ns == 3) MK_GR_TRI_(v[0], mid[0], mid[2]), MK_GR_TRI_(v[1], mid[1], mid[0]), MK_GR_TRI_(v[2], mid[2], mid[1]), MK_GR_TRI_(mid[0], mid[1], mid[2]);
            else if (ns == 1) {
                const int k = mid[0] >= 0 ? 0 : mid[1] >= 0 ? 1 : 2; /* edge k: v[k] -> v[k + 1] */
                MK_GR_TRI_(v[k], mid[k], v[(k + 2) % 3]), MK_GR_TRI_(mid[k], v[(k + 1) % 3], v[(k + 2) % 3]);
            } else {
                const int k = mid[0] < 0 ? 0 : mid[1] < 0 ? 1 : 2; /* edge k whole; mb on b-c, mc on c-a */
                const int a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3], mb = mid[(k + 1) % 3], mc = mid[(k + 2) % 3];
                MK_GR_TRI_(mc, mb, c);
                const double *pa = mV + 3 * (size_t)a, *pb = mV + 3 * (size_t)b, *pmb = mV + 3 * (size_t)mb, *pmc = mV + 3 * (size_t)mc;
                const double d1 = (pa[0] - pmb[0]) * (pa[0] - pmb[0]) + (pa[1] - pmb[1]) * (pa[1] - pmb[1]) + (pa[2] - pmb[2]) * (pa[2] - pmb[2]);
                const double d2 = (pb[0] - pmc[0]) * (pb[0] - pmc[0]) + (pb[1] - pmc[1]) * (pb[1] - pmc[1]) + (pb[2] - pmc[2]) * (pb[2] - pmc[2]);
                if (d1 <= d2) MK_GR_TRI_(a, b, mb), MK_GR_TRI_(a, mb, mc);
                else MK_GR_TRI_(a, b, mc), MK_GR_TRI_(b, mb, mc);
            }
        }
#undef MK_GR_TRI_
        free(mT), free(E), E = NULL;
        mT = T2, mt = nt2;
    }
    free(ps), free(E);
    out->V = mV, out->T = mT, out->nv = mv, out->nt = mt;
    return pass;
fail:
    free(mT), free(mV), free(ps), free(E);
    return -1;
}
/* Sets up g for S (closed, V[3 nv], T[3 nt]; kept by pointer, so it must outlive g) and the seam function fn(p, ctx); a
   psi not finite at a vertex reads as 1, one under nudge as nudge, and one just off the seam goes to its neighbours'
   side (below). Free g with mk_graft_free. Returns 0, or -1 for a bad mesh, a NULL fn, near or
   nudge negative or not finite, or no memory. */
static inline int mk_graft_init(mk_graft *g, const double *V, int nv, const int *T, int nt, mk_graft_psi fn, void *ctx, double near, double nudge)
{
    memset(g, 0, sizeof *g);
    if (!V || !T || !fn || nv < 3 || nt < 1 || nv > INT_MAX / 4 || nt > INT_MAX / 4 || !(near >= 0) || !isfinite(near) || !(nudge >= 0) || !isfinite(nudge)) return -1;
    for (size_t q = 0; q < 3 * (size_t)nv; q++)
        if (!isfinite(V[q])) return -1;
    for (size_t q = 0; q < 3 * (size_t)nt; q++)
        if (T[q] < 0 || T[q] >= nv) return -1;
    g->V = V, g->T = T, g->nv = nv, g->nt = nt, g->fn = fn, g->ctx = ctx, g->near = near;
    g->psi = malloc((size_t)nv * sizeof *g->psi);
    if (!g->psi || mk_bvh_build(&g->bv, V, T, nt)) { free(g->psi), g->psi = NULL; return -1; }
    for (int i = 0; i < nv; i++) {
        double s = fn(V + 3 * (size_t)i, ctx);
        if (!isfinite(s)) s = 1;
        if (fabs(s) < nudge) s = nudge;
        g->psi[i] = s;
    }
#if !(defined(MK_MUTATE) && MK_MUTATE == 372) /* 372: a vertex just off the seam left on its side */
    /* a vertex just off the seam (its psi under 1 % of the way to a neighbour's across it) goes to the side most of its
       neighbours are on: S's cut circles such a vertex at 1 % of its edges, and on the other side that detour is an
       inlet round most of it, which the patch's fans fold over (a muscle cut against a thin membrane: a vertex 16-30 um
       inside it, 4 crossing pairs); on their side the detour spans its few neighbours across, under half a turn */
    int *side = calloc((size_t)nv, sizeof *side);
    unsigned char *nr = calloc((size_t)nv, 1);
    if (!side || !nr) { free(side), free(nr), mk_graft_free(g); return -1; }
    for (size_t q = 0; q < 3 * (size_t)nt; q++) {
        const int a = T[q], b = T[q - q % 3 + (q % 3 + 1) % 3];
        const double pa = g->psi[a], pb = g->psi[b];
        side[a] += pb > 0 ? 1 : -1, side[b] += pa > 0 ? 1 : -1;
        if ((pa > 0) != (pb > 0)) nr[a] |= fabs(pa) < 0.0101 * fabs(pb), nr[b] |= fabs(pb) < 0.0101 * fabs(pa);
    }
    for (int i = 0; i < nv; i++)
        if (nr[i] && side[i] && (side[i] > 0) != (g->psi[i] > 0)) g->psi[i] = -g->psi[i];
    free(side), free(nr);
#endif
    return 0;
}
/* the patch's psi: S's linear one at the closest point of S within near, psi itself farther off */
static inline double mk_graft_psi_p_(const mk_graft *g, const double *p)
{
    int tr = -1;
    double q[3], w[3];
    mk_bvh_nearest(&g->bv, p, g->near * g->near, &tr, q);
    if (tr < 0) return g->fn(p, g->ctx);
    const int *f = g->T + 3 * (size_t)tr;
    mk_closest_pt_tri(p, g->V + 3 * (size_t)f[0], g->V + 3 * (size_t)f[1], g->V + 3 * (size_t)f[2], w, q);
    return w[0] * g->psi[f[0]] + w[1] * g->psi[f[1]] + w[2] * g->psi[f[2]];
}
/* a surface cut at psi = 0, the side sgn psi > 0 kept: the kept vertices (map), then the cut points; rim edges a -> b
   with their triangle. g NULL (S): the cut points linear in psi, held between 1 % and 99 % of their edge; else (the
   patch) bisected on the patch's psi (8 steps), held between 10 % and 90 %, as mk_mtet's vertices are. So no piece of
   a cut triangle has under 1e-4 of its area (S's) or 1e-2 (the patch's): a vertex just off the seam left slivers of
   1e-12 m^2 on body meshes, and at 1 % a patch of 0.5 mm cells still left pieces 3-6 um across. The seam moves by
   that share of an edge at most, well inside snap. */
static inline int mk_graft_clip_(const double *V, int nv, const int *T, int nt, const double *psi, double sgn, const mk_graft *g, mk_gcut_ *c)
{
    memset(c, 0, sizeof *c);
    c->map = malloc((size_t)(nv ? nv : 1) * sizeof *c->map);
    unsigned char *use = calloc((size_t)(nv ? nv : 1), 1);
    mk_gmap_ M = {0};
    long nx = 0;
    if (!c->map || !use) goto fail;
    for (int t = 0; t < nt; t++) {
        int k = 0;
        for (int q = 0; q < 3; q++) k += sgn * psi[T[3 * (size_t)t + q]] > 0;
        if (k)
            for (int q = 0; q < 3; q++)
                if (sgn * psi[T[3 * (size_t)t + q]] > 0) use[T[3 * (size_t)t + q]] = 1;
        nx += k == 1 || k == 2;
    }
    for (int i = 0; i < nv; i++) c->map[i] = use[i] ? mk_gcut_v_(c, V + 3 * (size_t)i) : -1;
    c->nkeep = c->nv;
    if (c->oom || mk_gmap_init_(&M, 4 * (size_t)nx + 16)) goto fail;
#define MK_CL_X_(i_, j_, out_)                                                                                                                    \
    do {                                                                                                                                          \
        const int ii_ = (i_), jj_ = (j_);                                                                                                         \
        const uint64_t key_ = ii_ < jj_ ? (uint64_t)ii_ << 32 | (uint32_t)jj_ : (uint64_t)jj_ << 32 | (uint32_t)ii_;                              \
        const size_t hs_ = mk_gmap_slot_(&M, key_);                                                                                               \
        if (M.k[hs_] == key_) out_ = M.v[hs_];                                                                                                    \
        else {                                                                                                                                    \
            const double *a_ = V + 3 * (size_t)ii_, *b_ = V + 3 * (size_t)jj_;                                                                    \
            double ta_ = 0, tb_ = 1, fa_ = psi[ii_], fb_ = psi[jj_], p_[3], t_ = fa_ / (fa_ - fb_);                                               \
            for (int it_ = 0; g && it_ < 8; it_++) {                                                                                              \
                MK_GRAFT_BISECT_(t_, ta_, tb_) /* a psi of +-inf or NaN: bisect */                                                                \
                for (int k_ = 0; k_ < 3; k_++) p_[k_] = a_[k_] + t_ * (b_[k_] - a_[k_]);                                                          \
                const double f_ = mk_graft_psi_p_(g, p_);                                                                                         \
                if (f_ == 0) break;                                                                                                               \
                if ((f_ > 0) == (fa_ > 0)) ta_ = t_, fa_ = f_;                                                                                    \
                else tb_ = t_, fb_ = f_;                                                                                                          \
                t_ = it_ % 2 ? 0.5 * (ta_ + tb_) : ta_ + (tb_ - ta_) * fa_ / (fa_ - fb_);                                                         \
            }                                                                                                                                     \
            MK_GRAFT_HOLD_(t_)                                                                                                                    \
            for (int k_ = 0; k_ < 3; k_++) p_[k_] = a_[k_] + t_ * (b_[k_] - a_[k_]);                                                              \
            out_ = mk_gcut_v_(c, p_);                                                                                                             \
            M.k[hs_] = key_, M.v[hs_] = out_;                                                                                                     \
        }                                                                                                                                         \
    } while (0)
    for (int t = 0; t < nt && !c->oom; t++) {
        const int *f = T + 3 * (size_t)t;
        int in[3], k = 0;
        for (int q = 0; q < 3; q++) in[q] = sgn * psi[f[q]] > 0, k += in[q];
        if (!k) continue;
        if (k == 3) { mk_gcut_t_(c, c->map[f[0]], c->map[f[1]], c->map[f[2]]); continue; }
        int r = 0;
        if (k == 1) {
            while (!in[r]) r++;
        } else {
            while (in[r]) r++;
            r = (r + 1) % 3;
        }
        const int v0 = f[r], v1 = f[(r + 1) % 3], v2 = f[(r + 2) % 3];
        int p, q2;
        if (k == 1) {
            MK_CL_X_(v0, v1, p); MK_CL_X_(v2, v0, q2);
            mk_gcut_e_(c, p, q2, mk_gcut_t_(c, c->map[v0], p, q2));
        } else {
            MK_CL_X_(v1, v2, p); MK_CL_X_(v2, v0, q2);
            mk_gcut_t_(c, c->map[v0], c->map[v1], p);
            mk_gcut_e_(c, p, q2, mk_gcut_t_(c, c->map[v0], p, q2));
        }
    }
#undef MK_CL_X_
    free(M.k), free(M.v), free(use);
    if (c->oom) { mk_gcut_free_(c); return -1; }
    return 0;
fail:
    free(M.k), free(M.v), free(use), mk_gcut_free_(c);
    return -1;
}
/* The field's surface (mk_mtet) kept where psi < 0, an open patch: out->keep 1 on its rim; free with mk_gmesh_free.
   Returns 0, or -1 as mk_mtet, for a g not set up, or for no memory. */
static inline int mk_graft_patch(const mk_graft *g, const float *F, const int n[3], const double lo[3], double h, mk_gmesh *out)
{
    memset(out, 0, sizeof *out);
    mk_gmesh mt;
    if (!g || !g->psi || mk_mtet(F, n, lo, h, &mt)) return -1;
    double *ps = malloc((size_t)(mt.nv ? mt.nv : 1) * sizeof *ps);
    if (!ps) { mk_gmesh_free(&mt); return -1; }
    for (int i = 0; i < mt.nv; i++) {
        ps[i] = mk_graft_psi_p_(g, mt.V + 3 * (size_t)i);
        if (!(fabs(ps[i]) >= 1e-300)) ps[i] = -1e-300; /* a node on the seam goes to the patch */
    }
    mk_gcut_ c;
    const int r = mk_graft_clip_(mt.V, mt.nv, mt.T, mt.nt, ps, -1, g, &c);
    free(ps), mk_gmesh_free(&mt);
    if (r) return -1;
    out->keep = calloc((size_t)(c.nv ? c.nv : 1), 1);
    if (!out->keep) { mk_gcut_free_(&c); return -1; }
    for (int e = 0; e < c.ne; e++) out->keep[c.E[3 * (size_t)e]] = out->keep[c.E[3 * (size_t)e + 1]] = 1;
    out->V = c.V, out->T = c.T, out->nv = c.nv, out->nt = c.nt;
    free(c.E), free(c.map);
    return 0;
}
static inline int mk_graft_cmpu_(const void *a, const void *b)
{
    const uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}
typedef struct { uint64_t k; int c; } mk_graft_de_; /* a directed edge a -> b (a << 32 | b) and its triangle's third corner */
static inline int mk_graft_decmp_(const void *a, const void *b)
{
    const uint64_t x = ((const mk_graft_de_ *)a)->k, y = ((const mk_graft_de_ *)b)->k;
    return (x > y) - (x < y);
}
static inline int mk_graft_cmp3_(const void *a, const void *b)
{
    const double *x = a, *y = b;
    return x[0] != y[0] ? (x[0] > y[0]) - (x[0] < y[0]) : (x[1] > y[1]) - (x[1] < y[1]);
}
static inline int mk_graft_find_(int *pr, int x)
{
    while (pr[x] != x) x = pr[x] = pr[pr[x]];
    return x;
}
/* The patch's inside near its rim follows the rim's move (with the rim moved alone, an inside vertex 10 um from a rim
   moved 0.2 mm onto S's ends past it, and its triangle folds through S's rim triangles). V: the patch (nv vertices,
   T nt triangles, rim[] 1 on its rim, D[3 p] the rim's moves, nx[p], pv[p] its rim neighbours, -1 or -2 when none or
   several). An inside vertex x whose nearest rim vertex p by edges is within R_max moves by d (1 - r / R)^2, with d
   the move at x's closest point q on p's two rim edges (linear along them), r = |x - q| and R = 3 |d|: a field whose
   gradient stays under 2/3 of the rim's, so the first rings keep their order. Returns 0, or -1 for no memory. */
static inline int mk_graft_follow_(double *V, int nv, const int *T, int nt, const unsigned char *rim, const double *D, const int *nx, const int *pv)
{
    double rmax = 0;
    for (int p = 0; p < nv; p++)
        if (rim[p]) rmax = fmax(rmax, 3 * sqrt(D[3 * (size_t)p] * D[3 * (size_t)p] + D[3 * (size_t)p + 1] * D[3 * (size_t)p + 1] + D[3 * (size_t)p + 2] * D[3 * (size_t)p + 2]));
    if (!(rmax > 0)) return 0;
    int *off = calloc((size_t)nv + 1, sizeof *off), *adj = malloc(6 * (size_t)(nt ? nt : 1) * sizeof *adj), *src = malloc((size_t)nv * sizeof *src);
    double *dist = malloc((size_t)nv * sizeof *dist), *hk = malloc((6 * (size_t)nt + (size_t)nv + 1) * sizeof *hk);
    int *hv = malloc((6 * (size_t)nt + (size_t)nv + 1) * sizeof *hv), ret = -1;
    unsigned char *seen = calloc((size_t)nv + 1, 1); /* each vertex spread from once: at most 6 nt + nv pushes */
    if (!off || !adj || !src || !dist || !hk || !hv || !seen) goto done;
    for (size_t q = 0; q < 3 * (size_t)nt; q++) off[T[q] + 1] += 2;
    for (int v = 0; v < nv; v++) off[v + 1] += off[v];
    for (int v = 0; v < nv; v++) src[v] = off[v]; /* fill pointers */
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) {
            const int a = T[3 * (size_t)t + k], b = T[3 * (size_t)t + (k + 1) % 3];
            adj[src[a]++] = b, adj[src[b]++] = a;
        }
    /* nearest rim vertex by edge length (Dijkstra from the whole rim, a binary heap with stale entries) */
    size_t hn = 0;
#define MK_GF_PUSH_(key, val)                                                                                                                     \
    do {                                                                                                                                          \
        size_t i_ = hn++;                                                                                                                         \
        while (i_ && hk[(i_ - 1) / 2] > (key)) hk[i_] = hk[(i_ - 1) / 2], hv[i_] = hv[(i_ - 1) / 2], i_ = (i_ - 1) / 2;                           \
        hk[i_] = (key), hv[i_] = (val);                                                                                                           \
    } while (0)
    for (int v = 0; v < nv; v++) {
        dist[v] = rim[v] ? 0 : INFINITY, src[v] = rim[v] ? v : -1;
        if (rim[v]) MK_GF_PUSH_(0.0, v);
    }
    while (hn) {
        const double dk = hk[0];
        const int v = hv[0];
        const double lk = hk[--hn];
        const int lv = hv[hn];
        size_t i = 0;
        for (;;) { /* the last entry sifted down from the top */
            size_t c = 2 * i + 1;
            if (c >= hn) break;
            if (c + 1 < hn && hk[c + 1] < hk[c]) c++;
            if (!(hk[c] < lk)) break;
            hk[i] = hk[c], hv[i] = hv[c], i = c;
        }
        if (hn) hk[i] = lk, hv[i] = lv;
        if (seen[v] || dk > dist[v] || dk > rmax) continue;
        seen[v] = 1;
        for (int q = off[v]; q < off[v + 1]; q++) {
            const int w = adj[q];
            const double *a = V + 3 * (size_t)v, *b = V + 3 * (size_t)w;
            const double nd = dk + sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
            if (nd < dist[w] && nd <= rmax) { dist[w] = nd, src[w] = src[v]; MK_GF_PUSH_(nd, w); }
        }
    }
#undef MK_GF_PUSH_
    for (int x = 0; x < nv; x++) {
        if (rim[x] || src[x] < 0) continue;
        const int p = src[x], ends[2] = {pv[p], nx[p]};
        const double *X = V + 3 * (size_t)x, *P = V + 3 * (size_t)p;
        double d[3], br = 0;
        for (int k = 0; k < 3; k++) d[k] = D[3 * (size_t)p + k], br += (X[k] - P[k]) * (X[k] - P[k]);
        for (int e = 0; e < 2; e++) { /* the closer of p's two rim edges */
            const int o = ends[e];
            if (o < 0) continue;
            const double *O = V + 3 * (size_t)o;
            double ab[3], l2 = 0, s = 0, r2 = 0;
            for (int k = 0; k < 3; k++) ab[k] = O[k] - P[k], l2 += ab[k] * ab[k], s += (X[k] - P[k]) * ab[k];
            s = l2 > 0 ? fmin(1, fmax(0, s / l2)) : 0;
            for (int k = 0; k < 3; k++) r2 += (X[k] - P[k] - s * ab[k]) * (X[k] - P[k] - s * ab[k]);
            if (r2 < br)
                for (int k = 0; k < 3; k++) br = r2, d[k] = (1 - s) * D[3 * (size_t)p + k] + s * D[3 * (size_t)o + k];
        }
        const double R = 3 * sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]), r = sqrt(br);
        if (!(r < R)) continue;
        const double f = (1 - r / R) * (1 - r / R);
        for (int k = 0; k < 3; k++) V[3 * (size_t)x + k] += f * d[k];
    }
    ret = 0;
done:
    free(off), free(adj), free(src), free(dist), free(hk), free(hv), free(seen);
    return ret;
}
/* S cut at psi = 0 with the patch (PV[3 pnv], PT[3 pnt]: open, its rim on the seam) stitched in. out: the result, keep
   as above; free with mk_gmesh_free. Returns 0 only for one closed, consistently wound piece wound as S is, no vertex
   pinched, no new triangle of zero area and no more crossing pairs with a changed triangle than S had there (none
   where S was clean). Else, with out still given and every fault found a bit in rep->faults, the first of: 2
   (MK_GRAFT_INVERTED) the kept piece wound against S, or pieces wound against S holding over 1 % of its volume
   (rep->inverted, inverted_vol: part of the mesh turned inside out by the stitch); 1 (MK_GRAFT_OPEN) not a closed
   manifold (rep->dup, lone, pinched); 3 (MK_GRAFT_SPLIT) pieces dropped holding over 1 % of S's volume, a cut that split
   the mesh (rep->ndropped, dropped_abs); 4 (MK_GRAFT_ZERO) a new triangle of zero area (rep->zero_area);
   5 (MK_GRAFT_CROSS) new crossing pairs (rep->cross0, cross1). Or -1, out empty, rep->reason saying why: MK_GRAFT_E_ARGS
   a bad patch, or snap negative or not finite; E_MEM no memory; E_RIM S's cut rim not closed loops, or no rim; E_SNAP a
   patch rim vertex farther than snap from S's rim (rep->snap its distance); E_MERGE one merged farther than snap from
   its rim neighbour; E_FOLLOW a rim that does not follow S's; E_EMPTY nothing left. rep->reason is the return code
   when that is over 0. rep may be NULL; it is filled on -1 too, as far as the stitch got. rep->loops_s, loops_p and
   unpaired say how the rims paired up (unpaired -1: not simple loops). */
static inline int mk_graft_stitch(const mk_graft *g, const double *PV, int pnv, const int *PT, int pnt, double snap, mk_gmesh *out, mk_graft_report *rep)
{
    mk_graft_report R = {0};
    memset(out, 0, sizeof *out);
    R.reason = MK_GRAFT_E_ARGS;
    if (rep) *rep = R;
    if (!g || !g->psi || !PV || !PT || pnv < 3 || pnt < 1 || pnv > INT_MAX / 4 || pnt > INT_MAX / 8 || !(snap >= 0) || !isfinite(snap)) return -1;
    for (size_t q = 0; q < 3 * (size_t)pnv; q++)
        if (!isfinite(PV[q])) return -1;
    for (size_t q = 0; q < 3 * (size_t)pnt; q++)
        if (PT[q] < 0 || PT[q] >= pnv) return -1;
    R.reason = 0;
#define MK_GR_FAIL_(code) do { R.reason = (code); goto out; } while (0)
    int ret = -1;
    mk_gcut_ CS = {0}, CP = {0};
    uint64_t *de = NULL;
    mk_graft_de_ *ed = NULL;
    int *startS = NULL, *nextS = NULL, *pe = NULL, *weld = NULL, *newid = NULL, *OT = NULL, *pr = NULL, *fin = NULL, *rt = NULL, *sl = NULL, *so = NULL, *slen = NULL;
    int *smet = NULL, *startP = NULL, *lv = NULL, *ltl = NULL, *lpre = NULL, *lst = NULL;
    double *pt = NULL, (*ins)[3] = NULL, *OV = NULL, *pv = NULL, *pq = NULL, *pd = NULL, *lw = NULL;
    unsigned char *dead = NULL, *usep = NULL, *usef = NULL, *lkeep = NULL, *pseen = NULL, *mark = NULL;
    mk_bvh OB = {0};
    int *NT = NULL, *buf = NULL; /* the new triangles; a box query's triangles */
    double emean = 0; /* S's mean edge: the zero-area scale (meshcheck's 1e-12 m^2 is 1e-6 of a 1 mm edge squared) */
    {
        double v = 0, es = 0;
        for (int t = 0; t < g->nt; t++) {
            const double *a = g->V + 3 * (size_t)g->T[3 * (size_t)t], *b = g->V + 3 * (size_t)g->T[3 * (size_t)t + 1], *c = g->V + 3 * (size_t)g->T[3 * (size_t)t + 2];
            v += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
            es += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])) + sqrt((b[0] - c[0]) * (b[0] - c[0]) + (b[1] - c[1]) * (b[1] - c[1]) + (b[2] - c[2]) * (b[2] - c[2])) +
                  sqrt((c[0] - a[0]) * (c[0] - a[0]) + (c[1] - a[1]) * (c[1] - a[1]) + (c[2] - a[2]) * (c[2] - a[2]));
        }
        R.vol0 = v, emean = es / (3.0 * g->nt);
    }
    if (mk_graft_clip_(g->V, g->nv, g->T, g->nt, g->psi, 1, NULL, &CS)) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    /* the patch with its rim edges: a directed edge without its reverse; a triangle with two rim edges split at its
       centre */
    for (int i = 0; i < pnv; i++) mk_gcut_v_(&CP, PV + 3 * (size_t)i);
    de = malloc(3 * (size_t)pnt * sizeof *de);
    if (!de || CP.oom) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    for (int t = 0; t < pnt; t++) {
        mk_gcut_t_(&CP, PT[3 * (size_t)t], PT[3 * (size_t)t + 1], PT[3 * (size_t)t + 2]);
        for (int q = 0; q < 3; q++) de[3 * (size_t)t + q] = (uint64_t)PT[3 * (size_t)t + q] << 32 | (uint32_t)PT[3 * (size_t)t + (q + 1) % 3];
    }
    qsort(de, 3 * (size_t)pnt, sizeof *de, mk_graft_cmpu_);
    for (int t = 0; t < pnt && !CP.oom; t++) {
        int rm[3], nr = 0;
        const int f[3] = {PT[3 * (size_t)t], PT[3 * (size_t)t + 1], PT[3 * (size_t)t + 2]};
        for (int q = 0; q < 3; q++) {
            const uint64_t rv = (uint64_t)f[(q + 1) % 3] << 32 | (uint32_t)f[q];
            rm[q] = !bsearch(&rv, de, 3 * (size_t)pnt, sizeof *de, mk_graft_cmpu_), nr += rm[q];
        }
        if (nr <= 1) {
            for (int q = 0; q < 3; q++)
                if (rm[q]) mk_gcut_e_(&CP, f[q], f[(q + 1) % 3], t);
            continue;
        }
        double gc[3];
        for (int k = 0; k < 3; k++) gc[k] = (PV[3 * (size_t)f[0] + k] + PV[3 * (size_t)f[1] + k] + PV[3 * (size_t)f[2] + k]) / 3;
        const int gv = mk_gcut_v_(&CP, gc);
        R.split++;
        for (int q = 0; q < 3; q++) {
            int tt;
            if (q == 0) CP.T[3 * (size_t)t] = f[0], CP.T[3 * (size_t)t + 1] = f[1], CP.T[3 * (size_t)t + 2] = gv, tt = t;
            else tt = mk_gcut_t_(&CP, f[q], f[(q + 1) % 3], gv);
            if (rm[q]) mk_gcut_e_(&CP, f[q], f[(q + 1) % 3], tt);
        }
    }
    if (CP.oom || CS.oom) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    R.rim_s = CS.ne, R.rim_p = CP.ne, R.patch_nv = CP.nv, R.patch_nt = CP.nt;
    /* S's rim: next[e] = the rim edge that starts where e ends */
    startS = malloc((size_t)(CS.nv ? CS.nv : 1) * sizeof *startS), nextS = malloc((size_t)(CS.ne ? CS.ne : 1) * sizeof *nextS);
    if (!startS || !nextS) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    for (int i = 0; i < CS.nv; i++) startS[i] = -1;
    long bad = 0;
    for (int e = 0; e < CS.ne; e++) {
        if (startS[CS.E[3 * (size_t)e]] >= 0) bad++;
        startS[CS.E[3 * (size_t)e]] = e;
    }
    for (int e = 0; e < CS.ne; e++)
        if ((nextS[e] = startS[CS.E[3 * (size_t)e + 1]]) < 0) bad++;
    if (bad || !CS.ne || !CP.ne) MK_GR_FAIL_(MK_GRAFT_E_RIM);
    /* each patch rim vertex onto S's rim: its edge, where along it, S's vertex when it lands on one */
    pe = malloc((size_t)CP.nv * sizeof *pe), weld = malloc((size_t)CP.nv * sizeof *weld), pt = malloc((size_t)CP.nv * sizeof *pt);
    if (!pe || !weld || !pt) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    pq = malloc(3 * (size_t)CP.nv * sizeof *pq), pd = malloc((size_t)CP.nv * sizeof *pd), rt = malloc((size_t)CP.nv * sizeof *rt);
    if (!pq || !pd || !rt) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    for (int i = 0; i < CP.nv; i++) pe[i] = -1, weld[i] = -1, rt[i] = i;
    double snap2 = 0;
    for (int e2 = 0; e2 < CP.ne; e2++) {
        const int p = CP.E[3 * (size_t)e2];
        if (pe[p] >= 0) continue;
        const double *x = CP.V + 3 * (size_t)p;
        double bd = INFINITY, bq[3] = {0, 0, 0}, bt = 0, l2s = 0;
        int be = -1;
        for (int e = 0; e < CS.ne; e++) {
            const double *a = CS.V + 3 * (size_t)CS.E[3 * (size_t)e], *c = CS.V + 3 * (size_t)CS.E[3 * (size_t)e + 1];
            double ab[3], l2 = 0, t = 0, q[3], d2 = 0;
            for (int k = 0; k < 3; k++) ab[k] = c[k] - a[k], l2 += ab[k] * ab[k], t += (x[k] - a[k]) * ab[k];
            t = l2 > 0 ? t / l2 : 0, t = t < 0 ? 0 : t > 1 ? 1 : t;
            for (int k = 0; k < 3; k++) q[k] = a[k] + t * ab[k], d2 += (q[k] - x[k]) * (q[k] - x[k]);
            if (d2 < bd) bd = d2, be = e, bt = t, l2s = l2, memcpy(bq, q, sizeof q);
        }
        if (bd > snap2) snap2 = bd;
#if defined(MK_MUTATE) && MK_MUTATE == 366
        if (be < 0 || !(bd <= snap * snap)) goto out; /* 366: the snap exit with no reason, rep->snap 0 */
#else
        if (be < 0 || !(bd <= snap * snap)) { R.snap = sqrt(snap2); MK_GR_FAIL_(MK_GRAFT_E_SNAP); }
#endif
        const double *a = CS.V + 3 * (size_t)CS.E[3 * (size_t)be], *c = CS.V + 3 * (size_t)CS.E[3 * (size_t)be + 1];
        double da = 0, dc = 0;
        for (int k = 0; k < 3; k++) da += (bq[k] - a[k]) * (bq[k] - a[k]), dc += (bq[k] - c[k]) * (bq[k] - c[k]);
        const double tie = 1e-4 * l2s; /* on S's vertex: within 1 % of the edge's length */
        pd[p] = bd;
        if (da <= tie) pe[p] = be, pt[p] = 0, weld[p] = CS.E[3 * (size_t)be], R.welded++, memcpy(pq + 3 * (size_t)p, a, sizeof bq);
        else if (dc <= tie) pe[p] = nextS[be], pt[p] = 0, weld[p] = CS.E[3 * (size_t)be + 1], R.welded++, memcpy(pq + 3 * (size_t)p, c, sizeof bq);
        else pe[p] = be, pt[p] = bt, R.inserted++, memcpy(pq + 3 * (size_t)p, bq, sizeof bq);
    }
    R.snap = sqrt(snap2);
    /* the rims as loops. S's rim edges are numbered along each loop, so a place on it is (loop, edge number + t). Each
       patch rim loop should lie on one of S's and each of S's be met by one; the loops that do not (rep->unpaired: a
       patch loop with no counterpart, as a hole in the patch by the seam or an island of it, one on two of S's loops, a
       loop of S's met by none or by two) are counted and placed vertex by vertex, without the ordering below: their
       pieces are mostly specks the cut leaves, but they can break the result, which return 1 then says. Along a patch
       loop alone on its loop of S's the places must fall strictly: the longest such run, each place measured behind a
       start (the vertex nearest S's rim among those in order with both rim neighbours), is kept, and each vertex out of
       it (a dense rim's wiggle stepping back on S's coarser edges, which would put an edge in twice) is merged into the
       kept rim neighbour before or after it, moving no farther than snap (else -1); a merged vertex's rim triangle is
       dropped. Rims that are not simple loops are placed vertex by vertex too. */
    {
        sl = malloc((size_t)CS.ne * sizeof *sl), so = malloc((size_t)CS.ne * sizeof *so), slen = malloc((size_t)CS.ne * sizeof *slen);
        smet = calloc((size_t)CS.ne, sizeof *smet), startP = malloc((size_t)CP.nv * sizeof *startP), lv = malloc((size_t)CP.ne * sizeof *lv);
        lw = malloc((size_t)CP.ne * sizeof *lw), ltl = malloc((size_t)CP.ne * sizeof *ltl), lpre = malloc((size_t)CP.ne * sizeof *lpre);
        lkeep = malloc((size_t)CP.ne), pseen = calloc((size_t)CP.ne, 1), lst = malloc(((size_t)CP.ne + 1) * sizeof *lst);
        if (!sl || !so || !slen || !smet || !startP || !lv || !lw || !ltl || !lpre || !lkeep || !pseen || !lst) MK_GR_FAIL_(MK_GRAFT_E_MEM);
        int simple = 1, nlp = 0;
        for (int e = 0; e < CS.ne; e++) sl[e] = -1;
        for (int e = 0; e < CS.ne && simple; e++) {
            if (sl[e] >= 0) continue;
            int k = 0, cur = e;
            do {
                if (sl[cur] >= 0) { simple = 0; break; } /* S's rim is not simple loops */
                sl[cur] = R.loops_s, so[cur] = k++, cur = nextS[cur];
            } while (cur != e);
            slen[R.loops_s++] = k;
        }
        for (int i = 0; i < CP.nv; i++) startP[i] = -1;
        for (int e2 = 0; e2 < CP.ne && simple; e2++) {
            if (startP[CP.E[3 * (size_t)e2]] >= 0) simple = 0; /* two rim edges from one patch vertex */
            startP[CP.E[3 * (size_t)e2]] = e2;
        }
#define MK_GR_POS_(v) (so[pe[v]] + pt[v])
/* MK_GR_BACK_: how far b lies behind a along a loop of n edges, in [0, n) */
#define MK_GR_BACK_(a, b, n) (MK_GR_POS_(a) - MK_GR_POS_(b) < 0 ? MK_GR_POS_(a) - MK_GR_POS_(b) + (n) : MK_GR_POS_(a) - MK_GR_POS_(b))
        /* the patch's rim loops, one after another in lv (loop q from lst[q]) */
        for (int e0 = 0, k = 0; e0 < CP.ne && simple; e0++) {
            if (pseen[e0]) continue;
            int cur = e0;
            lst[nlp] = k;
            do {
                if (cur < 0 || pseen[cur]) { simple = 0; break; } /* the patch's rim is not simple loops */
                pseen[cur] = 1, lv[k++] = CP.E[3 * (size_t)cur], cur = startP[CP.E[3 * (size_t)cur + 1]];
            } while (cur != e0);
            lst[++nlp] = k;
        }
        if (!simple) R.loops_s = nlp = 0, R.unpaired = -1;
        R.loops_p = nlp;
        for (int q = 0; q < nlp; q++) {
            const int *v = lv + lst[q], k = lst[q + 1] - lst[q], L = sl[pe[v[0]]];
            int one = 1;
            for (int i = 1; i < k; i++) one &= sl[pe[v[i]]] == L;
            if (one) smet[L]++; /* on one of S's loops */
            else R.unpaired++;
        }
        for (int L = 0; L < R.loops_s; L++) R.unpaired += smet[L] == 0 ? 1 : smet[L] - 1; /* met by none, or by 2+ */
#if defined(MK_MUTATE) && MK_MUTATE == 205
        R.unpaired = 0; /* 205: the loops not counted */
#endif
        for (int q = 0; q < nlp; q++) {
            const int *v = lv + lst[q], k = lst[q + 1] - lst[q], L = sl[pe[v[0]]];
            int one = 1;
            for (int i = 1; i < k; i++) one &= sl[pe[v[i]]] == L;
            if (!one || smet[L] != 1) continue; /* not alone on one loop of S's: no ordering */
            const double n = slen[L];
            int s = -1; /* the start: nearest S's rim among vertices in order with both rim neighbours, else nearest */
            for (int i = 0; i < k; i++) {
                const int a = v[(i + k - 1) % k], b = v[i], c = v[(i + 1) % k];
                if (MK_GR_BACK_(a, b, n) < 0.5 * n && MK_GR_BACK_(b, c, n) < 0.5 * n && (s < 0 || pd[b] < pd[v[s]])) s = i;
            }
            for (int i = 0; s < 0 && i < k; i++)
                if (s < 0 || pd[v[i]] < pd[v[s]]) s = i;
            /* each place as how far behind the start it lies, in (-n, 0]: the longest strictly falling run goes once
               round */
            int nl = 0;
            lw[0] = 0, lkeep[0] = 1;
            for (int j = 1; j < k; j++) {
                lw[j] = -MK_GR_BACK_(v[s], v[(s + j) % k], n), lkeep[j] = 0, lpre[j] = -1;
#if defined(MK_MUTATE) && MK_MUTATE == 202
                lkeep[j] = 1; /* 202: every vertex kept, in or out of order */
                continue;
#endif
                if (!(lw[j] < 0 && lw[j] > -n)) continue;
                int lo = 0, hi = nl;
                while (lo < hi) {
                    const int m = (lo + hi) / 2;
                    if (lw[ltl[m]] > lw[j]) lo = m + 1;
                    else hi = m;
                }
                lpre[j] = lo ? ltl[lo - 1] : 0, ltl[lo] = j, nl += lo == nl;
            }
            for (int j = nl ? ltl[nl - 1] : 0; j > 0; j = lpre[j]) lkeep[j] = 1;
            /* each run out of order between kept a (before) and b (after): its first part into a, the rest into b, by
               which is nearer */
            for (int j = 0; j < k;) {
                int j2 = j + 1;
                while (j2 < k && !lkeep[j2]) j2++;
                const int a = v[(s + j) % k], b = v[(s + j2) % k];
                const double *fa = weld[a] >= 0 ? CS.V + 3 * (size_t)weld[a] : pq + 3 * (size_t)a, *fb = weld[b] >= 0 ? CS.V + 3 * (size_t)weld[b] : pq + 3 * (size_t)b;
                int na = 0;
                for (int i = j + 1; i < j2; i++) {
                    const double *x = CP.V + 3 * (size_t)v[(s + i) % k];
                    double da = 0, db = 0;
                    for (int z = 0; z < 3; z++) da += (x[z] - fa[z]) * (x[z] - fa[z]), db += (x[z] - fb[z]) * (x[z] - fb[z]);
                    na += da <= db;
                }
                for (int i = j + 1; i < j2; i++) {
                    const int w = v[(s + i) % k], r = i - j - 1 < na ? a : b;
                    const double *x = CP.V + 3 * (size_t)w, *f = r == a ? fa : fb;
                    double d2 = 0;
                    for (int z = 0; z < 3; z++) d2 += (x[z] - f[z]) * (x[z] - f[z]);
                    if (!(d2 <= snap * snap)) { R.snap = fmax(R.snap, sqrt(d2)); MK_GR_FAIL_(MK_GRAFT_E_MERGE); }
                    if (weld[w] >= 0) R.welded--;
                    else R.inserted--;
                    rt[w] = r, pe[w] = pe[r], pt[w] = pt[r], weld[w] = -1, R.merged++;
                }
                j = j2;
            }
        }
#undef MK_GR_POS_
#undef MK_GR_BACK_
#if !(defined(MK_MUTATE) && MK_MUTATE == 369) /* 369: the rim moved alone */
        {
            double *D = malloc(3 * (size_t)CP.nv * sizeof *D);
            int *nx = malloc((size_t)CP.nv * sizeof *nx), *pv2 = malloc((size_t)CP.nv * sizeof *pv2);
            unsigned char *rim = calloc((size_t)CP.nv, 1);
            int rf = !D || !nx || !pv2 || !rim;
            for (int p = 0; !rf && p < CP.nv; p++) {
                nx[p] = pv2[p] = -1;
                if (pe[p] < 0) continue;
                const int q = rt[p];
                const double *f = weld[q] >= 0 ? CS.V + 3 * (size_t)weld[q] : pq + 3 * (size_t)q;
                for (int k = 0; k < 3; k++) D[3 * (size_t)p + k] = f[k] - CP.V[3 * (size_t)p + k];
                rim[p] = 1;
            }
            for (int e2 = 0; !rf && e2 < CP.ne; e2++) {
                const int a = CP.E[3 * (size_t)e2], b = CP.E[3 * (size_t)e2 + 1];
                nx[a] = nx[a] == -1 ? b : -2, pv2[b] = pv2[b] == -1 ? a : -2;
            }
            rf = rf || mk_graft_follow_(CP.V, CP.nv, CP.T, CP.nt, rim, D, nx, pv2);
            free(D), free(nx), free(pv2), free(rim);
            if (rf) MK_GR_FAIL_(MK_GRAFT_E_MEM);
        }
#endif
        for (int p = 0; p < CP.nv; p++)
            if (pe[p] >= 0 && rt[p] == p && weld[p] < 0) memcpy(CP.V + 3 * (size_t)p, pq + 3 * (size_t)p, 3 * sizeof *pq);
    }
#define MK_GR_PID_(p) (weld[rt[p]] >= 0 ? weld[rt[p]] : CS.nv + rt[p])
    /* S's rim triangles fanned over the patch vertices put into their rim edge */
    {
        ins = malloc((size_t)(R.inserted ? R.inserted : 1) * sizeof *ins);
        if (!ins) MK_GR_FAIL_(MK_GRAFT_E_MEM);
        long m = 0;
        for (int p = 0; p < CP.nv; p++)
            if (pe[p] >= 0 && rt[p] == p && weld[p] < 0) ins[m][0] = pe[p], ins[m][1] = pt[p], ins[m][2] = p, m++;
        qsort(ins, (size_t)m, sizeof *ins, mk_graft_cmp3_);
        for (long i = 0; i < m;) {
            long j = i;
            while (j < m && ins[j][0] == ins[i][0]) j++;
            const int e = (int)ins[i][0], a = CS.E[3 * (size_t)e], c = CS.E[3 * (size_t)e + 1], t = CS.E[3 * (size_t)e + 2];
            int *f = CS.T + 3 * (size_t)t, r = 0;
            while (r < 3 && !(f[r] == a && f[(r + 1) % 3] == c)) r++;
            if (r == 3) MK_GR_FAIL_(MK_GRAFT_E_FOLLOW);
            const int x = f[(r + 2) % 3];
            int prev = a;
            for (long q = i; q < j; q++) {
                const int pid = CS.nv + (int)ins[q][2];
                if (q == i) f[0] = prev, f[1] = pid, f[2] = x;
                else mk_gcut_t_(&CS, prev, pid, x), f = CS.T + 3 * (size_t)t;
                prev = pid;
            }
            mk_gcut_t_(&CS, prev, c, x);
            i = j;
        }
        if (CS.oom) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    }
    /* the patch's rim triangles fanned over S's rim vertices between their two ends (S's vertex i carried as -1 - i) */
    const int ncpt = CP.nt;
    dead = calloc((size_t)ncpt + 1, 1);
    if (!dead) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    for (int e2 = 0; e2 < CP.ne; e2++) {
        const int p = CP.E[3 * (size_t)e2], p2 = CP.E[3 * (size_t)e2 + 1], t = CP.E[3 * (size_t)e2 + 2];
        if (pe[p] < 0 || pe[p2] < 0) MK_GR_FAIL_(MK_GRAFT_E_FOLLOW);
        int *f = CP.T + 3 * (size_t)t, r = 0;
        while (r < 3 && !(f[r] == p && f[(r + 1) % 3] == p2)) r++;
        if (r == 3) MK_GR_FAIL_(MK_GRAFT_E_FOLLOW);
        const int x = f[(r + 2) % 3];
        if (MK_GR_PID_(p) == MK_GR_PID_(p2)) { dead[t] = 1, R.dropped++; continue; }
        int mid[64], nm = 0, cur = pe[p2];
        if (pe[p2] != pe[p]) {
            for (;;) {
                if (nm == 64) MK_GR_FAIL_(MK_GRAFT_E_FOLLOW);
                mid[nm++] = CS.E[3 * (size_t)cur + 1], cur = nextS[cur];
                if (cur == pe[p]) break;
            }
            if (pt[p] == 0 && weld[rt[p]] >= 0) nm--;
        }
        if (!nm) continue;
        R.fanned++;
        int prev = p, first = 1;
        for (int q = nm - 1; q >= -1; q--) {
            const int nxt = q >= 0 ? -1 - mid[q] : p2;
            if (first) f[0] = prev, f[1] = nxt, f[2] = x, first = 0;
            else mk_gcut_t_(&CP, prev, nxt, x), f = CP.T + 3 * (size_t)t;
            prev = nxt;
        }
        if (CP.oom) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    }
    /* the pieces: triangles with final numbers; the piece of most volume kept */
    const int npv = CP.nv;
    newid = malloc((size_t)npv * sizeof *newid), usep = calloc((size_t)npv + 1, 1);
    if (!newid || !usep) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    for (int t = 0; t < CP.nt; t++) {
        if (t < ncpt && dead[t]) continue;
        for (int q = 0; q < 3; q++) {
            const int v = CP.T[3 * (size_t)t + q];
            if (v >= 0 && weld[rt[v]] < 0) usep[rt[v]] = 1;
        }
    }
    for (int t = 0; t < CS.nt; t++)
        for (int q = 0; q < 3; q++) {
            const int v = CS.T[3 * (size_t)t + q];
            if (v >= CS.nv) usep[v - CS.nv] = 1;
        }
    int no = CS.nv;
    for (int p = 0; p < npv; p++) newid[p] = usep[p] ? no++ : -1;
    OV = malloc(3 * (size_t)no * sizeof *OV), OT = malloc(3 * ((size_t)CS.nt + CP.nt) * sizeof *OT);
    if (!OV || !OT) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    memcpy(OV, CS.V, 3 * (size_t)CS.nv * sizeof *OV);
    for (int p = 0; p < npv; p++)
        if (usep[p]) memcpy(OV + 3 * (size_t)newid[p], CP.V + 3 * (size_t)p, 3 * sizeof *OV);
    long nto = 0;
#define MK_GR_EMIT_(a, b, c)                                                                                                                      \
    do {                                                                                                                                          \
        const int a_ = (a), b_ = (b), c_ = (c);                                                                                                   \
        if (a_ != b_ && b_ != c_ && c_ != a_) OT[3 * nto] = a_, OT[3 * nto + 1] = b_, OT[3 * nto + 2] = c_, nto++;                                \
    } while (0)
#define MK_GR_SID_(v) ((v) >= CS.nv ? newid[(v) - CS.nv] : (v))
#define MK_GR_QID_(v) ((v) < 0 ? -1 - (v) : weld[rt[v]] >= 0 ? weld[rt[v]] : newid[rt[v]])
    for (int t = 0; t < CS.nt; t++) MK_GR_EMIT_(MK_GR_SID_(CS.T[3 * (size_t)t]), MK_GR_SID_(CS.T[3 * (size_t)t + 1]), MK_GR_SID_(CS.T[3 * (size_t)t + 2]));
    for (int t = 0; t < CP.nt; t++) {
        if (t < ncpt && dead[t]) continue;
#if defined(MK_MUTATE) && MK_MUTATE == 182
        if (t >= ncpt) continue; /* 182: the patch's fans left out */
#endif
        MK_GR_EMIT_(MK_GR_QID_(CP.T[3 * (size_t)t]), MK_GR_QID_(CP.T[3 * (size_t)t + 1]), MK_GR_QID_(CP.T[3 * (size_t)t + 2]));
    }
    pr = malloc((size_t)no * sizeof *pr), pv = calloc((size_t)no, sizeof *pv);
    if (!pr || !pv) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    for (int i = 0; i < no; i++) pr[i] = i;
    for (long t = 0; t < nto; t++) {
        pr[mk_graft_find_(pr, OT[3 * t + 1])] = mk_graft_find_(pr, OT[3 * t]);
        pr[mk_graft_find_(pr, OT[3 * t + 2])] = mk_graft_find_(pr, OT[3 * t]);
    }
    for (long t = 0; t < nto; t++) {
        const double *a = OV + 3 * (size_t)OT[3 * t], *b = OV + 3 * (size_t)OT[3 * t + 1], *c = OV + 3 * (size_t)OT[3 * t + 2];
        pv[mk_graft_find_(pr, OT[3 * t])] += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
    }
    /* the piece of most volume whatever its sign (with the most positive, a result come out inverted or broken, its
       volume negative, would lose to a closed speck and pass); a piece wound against S is counted, and 2 returned when
       the kept one is */
    int big = -1;
    const double sg = R.vol0 < 0 ? -1 : 1;
    for (int i = 0; i < no; i++)
        if (pr[i] == i && pv[i] != 0) {
            R.pieces++, R.inverted += sg * pv[i] < 0, R.inverted_vol += sg * pv[i] < 0 ? pv[i] : 0;
#if defined(MK_MUTATE) && MK_MUTATE == 203
            if (big < 0 || pv[i] > pv[big]) big = i; /* 203: the most positive volume */
#else
            if (big < 0 || fabs(pv[i]) > fabs(pv[big])) big = i;
#endif
        }
    for (int i = 0; i < no; i++)
        if (pr[i] == i && pv[i] != 0 && i != big) R.dropped_vol += pv[i], R.dropped_abs += fabs(pv[i]), R.ndropped++;
    if (big < 0) MK_GR_FAIL_(MK_GRAFT_E_EMPTY);
    fin = malloc((size_t)no * sizeof *fin), usef = calloc((size_t)no + 1, 1);
    if (!fin || !usef) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    for (long t = 0; t < nto; t++)
        if (mk_graft_find_(pr, OT[3 * t]) == big) usef[OT[3 * t]] = usef[OT[3 * t + 1]] = usef[OT[3 * t + 2]] = 1;
    int nfin = 0;
    for (int i = 0; i < no; i++) fin[i] = usef[i] ? nfin++ : -1;
    out->V = malloc(3 * (size_t)(nfin ? nfin : 1) * sizeof *out->V), out->keep = malloc((size_t)(nfin ? nfin : 1)), out->T = malloc(3 * (size_t)(nto ? nto : 1) * sizeof *out->T);
    ed = malloc(3 * (size_t)(nto ? nto : 1) * sizeof *ed);
    if (!out->V || !out->keep || !out->T || !ed) MK_GR_FAIL_(MK_GRAFT_E_MEM);
    for (int i = 0; i < no; i++) {
        if (!usef[i]) continue;
        memcpy(out->V + 3 * (size_t)fin[i], OV + 3 * (size_t)i, 3 * sizeof *OV);
        out->keep[fin[i]] = i < CS.nkeep ? 1 : i < CS.nv ? 3 : 0;
    }
    for (int p = 0; p < npv; p++)
        if (usep[p] && usef[newid[p]] && pe[p] >= 0) out->keep[fin[newid[p]]] = 2;
    size_t ned = 0;
    for (long t = 0; t < nto; t++) {
        if (mk_graft_find_(pr, OT[3 * t]) != big) continue;
        const int a = fin[OT[3 * t]], b = fin[OT[3 * t + 1]], c = fin[OT[3 * t + 2]];
        int *o = out->T + 3 * (size_t)out->nt++;
        o[0] = a, o[1] = b, o[2] = c;
        ed[ned].k = (uint64_t)a << 32 | (uint32_t)b, ed[ned++].c = c, ed[ned].k = (uint64_t)b << 32 | (uint32_t)c, ed[ned++].c = a;
        ed[ned].k = (uint64_t)c << 32 | (uint32_t)a, ed[ned++].c = b;
    }
    out->nv = nfin, R.vol1 = pv[big];
    qsort(ed, ned, sizeof *ed, mk_graft_decmp_);
    for (size_t i = 0; i < ned; i++) {
        if (i && ed[i].k == ed[i - 1].k) R.dup++;
        const mk_graft_de_ rv = {(ed[i].k & 0xffffffffull) << 32 | ed[i].k >> 32, 0};
        if (!bsearch(&rv, ed, ned, sizeof *ed, mk_graft_decmp_)) R.lone++;
    }
#if !(defined(MK_MUTATE) && MK_MUTATE == 365) /* 365: a pinched vertex not looked for */
    /* a pinched vertex: closed and wound edge by edge, its triangles make more than one fan (two pieces touching at a
       point, a mirrored pair left by an ear). The fan from v's first edge v -> b: b, then the third corner c of v -> b,
       then that of v -> c, ...; it must take in all of v's edges */
    for (size_t i = 0; !R.dup && !R.lone && i < ned;) {
        const uint64_t v = ed[i].k >> 32;
        size_t j = i;
        while (j < ned && ed[j].k >> 32 == v) j++;
        uint64_t b = ed[i].k & 0xffffffffull;
        size_t n = 0;
        do {
            const mk_graft_de_ key = {v << 32 | b, 0}, *f = bsearch(&key, ed + i, j - i, sizeof *ed, mk_graft_decmp_);
            if (!f) break;
            b = (uint64_t)(uint32_t)f->c, n++;
        } while (b != (ed[i].k & 0xffffffffull) && n <= j - i);
        R.pinched += n != j - i;
        i = j;
    }
#endif
    /* the triangles that are not S's own (a corner not S's, where slivers and crossings arise): none of zero area,
       and no more crossing pairs with one of them than S had with its triangles not kept whole. The result's pairs: those
       of two new triangles on a tree of the new ones alone, then each new one against S's own through S's tree (a
       triangle of S's in the result: kept whole, in the kept piece), so the whole result needs no tree of its own (a
       946,000-triangle diaphragm: 4 s a stitch with one) */
    {
        mark = malloc((size_t)g->nt + (size_t)out->nt + 1);
        if (!mark) MK_GR_FAIL_(MK_GRAFT_E_MEM);
        unsigned char *mo = mark + g->nt;
        const double z2 = 1e-12 * emean * emean * emean * emean; /* |2A|^2 under (1e-6 e^2)^2 */
        for (int t = 0; t < g->nt; t++) {
            const int *f = g->T + 3 * (size_t)t;
            mark[t] = !(g->psi[f[0]] > 0 && g->psi[f[1]] > 0 && g->psi[f[2]] > 0);
        }
        R.cross0 = mk_bvh_crossings(&g->bv, mark);
        if (R.cross0 < 0) MK_GR_FAIL_(MK_GRAFT_E_MEM);
        long nn = 0;
        for (int t = 0; t < out->nt; t++) {
            const int *f = out->T + 3 * (size_t)t;
            mo[t] = out->keep[f[0]] != 1 || out->keep[f[1]] != 1 || out->keep[f[2]] != 1;
            if (!mo[t]) continue;
            nn++;
            const double *a = out->V + 3 * (size_t)f[0], *b = out->V + 3 * (size_t)f[1], *c = out->V + 3 * (size_t)f[2];
            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            const double x = u[1] * w[2] - u[2] * w[1], y = u[2] * w[0] - u[0] * w[2], z = u[0] * w[1] - u[1] * w[0];
            R.zero_area += !(x * x + y * y + z * z >= z2);
        }
        NT = malloc(3 * (size_t)(nn ? nn : 1) * sizeof *NT), buf = malloc(64 * sizeof *buf);
        if (!NT || !buf) MK_GR_FAIL_(MK_GRAFT_E_MEM);
        nn = 0;
        for (int t = 0; t < out->nt; t++)
            if (mo[t]) memcpy(NT + 3 * (size_t)nn++, out->T + 3 * (size_t)t, 3 * sizeof *NT);
        if (mk_bvh_build(&OB, out->V, NT, (int)nn) || (R.cross1 = mk_bvh_crossings(&OB, NULL)) < 0) MK_GR_FAIL_(MK_GRAFT_E_MEM);
        int cap = 64;
        for (long u = 0; u < nn; u++) {
            mk_aabb bx;
            mk_tri_box_(out->V, NT, (int)u, &bx);
            int n = mk_bvh_box(&g->bv, &bx, buf, cap);
            if (n > cap) {
                int *b2 = realloc(buf, (size_t)n * sizeof *b2);
                if (!b2) MK_GR_FAIL_(MK_GRAFT_E_MEM);
                buf = b2, cap = n, n = mk_bvh_box(&g->bv, &bx, buf, cap);
            }
            for (int i = 0; i < n; i++) {
                const int s = buf[i], *f = g->T + 3 * (size_t)s;
                if (mark[s]) continue; /* not kept whole */
                int Y[3], in = 1;
                for (int k = 0; k < 3; k++) Y[k] = fin[CS.map[f[k]]], in &= Y[k] >= 0;
                if (in) R.cross1 += mk_tri_cross(out->V, NT + 3 * (size_t)u, Y);
            }
        }
    }
    /* the faults, each a bit; rc the first by precedence. 2: the kept piece wound against S, or pieces wound against S
       holding over 1 % of its volume: the stitch turned part of the mesh inside out (where a stitch broke the result, 1
       to 70 %; the rounding specks of a clean stitch, under 0.1 %); it comes first, since an inverted piece is often
       open too, and rc 1 would hide it. 3: pieces dropped holding over 1 % of S's volume (a cut that split a muscle
       dropped 22 % of it) */
#if defined(MK_MUTATE) && MK_MUTATE == 204
    const int inv = 0; /* 204: a kept piece wound against S passes */
#elif defined(MK_MUTATE) && MK_MUTATE == 210
    const int inv = sg * R.vol1 < 0; /* 210: only the kept piece judged */
#else
    const int inv = sg * R.vol1 < 0 || fabs(R.inverted_vol) > 1e-2 * fabs(R.vol1);
#endif
    R.faults = (unsigned)inv << MK_GRAFT_INVERTED | (unsigned)(R.dup || R.lone || R.pinched) << MK_GRAFT_OPEN;
#if !(defined(MK_MUTATE) && MK_MUTATE == 361) /* 361: a split not judged */
    R.faults |= (unsigned)(R.dropped_abs > 1e-2 * fabs(R.vol0)) << MK_GRAFT_SPLIT;
#endif
#if !(defined(MK_MUTATE) && MK_MUTATE == 362) /* 362: zero area not judged */
    R.faults |= (unsigned)(R.zero_area > 0) << MK_GRAFT_ZERO;
#endif
#if !(defined(MK_MUTATE) && MK_MUTATE == 364) /* 364: crossings not judged */
    R.faults |= (unsigned)(R.cross1 > R.cross0) << MK_GRAFT_CROSS;
#endif
    {
#if defined(MK_MUTATE) && MK_MUTATE == 360
        static const int ORD[5] = {MK_GRAFT_OPEN, MK_GRAFT_INVERTED, MK_GRAFT_SPLIT, MK_GRAFT_ZERO, MK_GRAFT_CROSS}; /* 360: rc 1 first, hiding 2 */
#else
        static const int ORD[5] = {MK_GRAFT_INVERTED, MK_GRAFT_OPEN, MK_GRAFT_SPLIT, MK_GRAFT_ZERO, MK_GRAFT_CROSS};
#endif
        ret = 0;
        for (int k = 0; k < 5 && !ret; k++)
            if (R.faults >> ORD[k] & 1) ret = ORD[k];
        R.reason = ret;
    }
#undef MK_GR_EMIT_
#undef MK_GR_SID_
#undef MK_GR_QID_
#undef MK_GR_PID_
out:
#undef MK_GR_FAIL_
    if (ret < 0) mk_gmesh_free(out);
    if (rep) *rep = R;
    mk_gcut_free_(&CS), mk_gcut_free_(&CP), mk_bvh_free(&OB), free(mark), free(NT), free(buf);
    free(de), free(ed), free(startS), free(nextS), free(pe), free(weld), free(newid), free(OT), free(pr), free(fin), free(pt), free(ins), free(OV), free(pv);
    free(dead), free(usep), free(usef), free(rt), free(sl), free(so), free(slen), free(smet), free(startP), free(lv), free(ltl), free(lpre);
    free(pq), free(pd), free(lw), free(lkeep), free(pseen), free(lst);
    return ret;
}
#endif
