/* mathkit: a closed mesh changed only near a seam, the rest of it kept triangle for triangle. mesher-fable's bonegraft (2026-11-15):
   a field's surface by marching tetrahedra on a box of nodes, cut at a seam function's zero, stitched into the mesh cut at the same
   zero, its rim put into the mesh's rim and welded. The bone-specific field (bones grown, opened, blurred and blended into the
   muscle's own distance) stays with its owner; the field is the caller's, on mk_sdf's grid layout (a grid from mk_sdf_mesh fits,
   and so does any field built from one; mk_sdf_mesh_box grids only the seam's neighbourhood, 1.5.2).
   mk_mtet: marching tetrahedra, Freudenthal's six tetrahedra a cube round its diagonal, the same in every cube, so neighbours'
   faces match and the surface is a closed manifold inside the box whatever the field (surface nets, mk_snet, can touch themselves
   where two sheets share a cell); a vertex on each node pair of unlike sign at the linear crossing held between 0.1 and 0.9 of the
   edge (no needles from a crossing at a node); wound outward, the field negative inside.
   A graft: mk_graft_init(S, psi) takes the untouched closed mesh and the seam function psi (negative where the patch goes; zero at
   the seam), psi at S's vertices (nudged off zero by `nudge`, so the zero line misses them), linear over S's triangles.
   mk_graft_patch(field) is the field's surface kept where psi < 0, cut on its edges at psi = 0 (bisected on the patch's psi: S's
   linear psi at the closest point of S within `near` of it, psi itself farther off). The caller may thin the patch (mk_decimate on
   an open surface keeps its rim on its rim). Unthinned it still stitches closed (1.5.2), but its small rim triangles follow the
   field's ripple about S's plane and can cross S's rim triangles (t_graft's dent, by meshcheck: 4 crossing pairs, 34 on a coarse
   S; thinned, none): thin it where crossings matter. mk_graft_stitch(patch) cuts S at psi = 0 (S kept where psi > 0,
   its crossing triangles cut on their edges), moves each patch rim vertex onto S's cut rim (no more than `snap`, or it refuses;
   one whose closest point lies within 1 % of an S rim edge's length of its end goes on to that S vertex, uncounted in snap and
   rep->snap),
   pairs the rims as loops (each patch rim loop on one of S's, each of S's met by one, else counted), merges the patch rim vertices
   whose places on S's rim fall out of order into a rim neighbour (1.5.2), puts each rim's vertices into the other's rim edges (the
   rim triangles fanned, S's in their own planes), welds, and keeps the piece of most volume, either sign (closed specks the cut
   left off are dropped and counted; pieces wound against S are reported, 1.5.2). out->keep marks for mk_decimate: 1 S's own
   vertices, 3 its cut points, 2 the patch's rim, 0 the rest. Header only, C23.
   Experimental (1.6.0): not frozen with the numerical core; its results may still change between releases. */
#ifndef MK_GRAFT_H
#define MK_GRAFT_H
#include "mk_geom.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct { double *V; int *T; int nv, nt; unsigned char *keep; } mk_gmesh; /* malloc'd; keep may be NULL */
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
    if (!m->k || !m->v) { free(m->k), free(m->v); return -1; } /* 233: the pointers left set (1.5.3) */
#else
    if (!m->k || !m->v) { free(m->k), free(m->v), m->k = NULL, m->v = NULL; return -1; } /* the callers free them again on failure */
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
#define MK_GRAFT_BISECT_(t_, ta_, tb_) /* 232: no bisection, so a psi of +inf gives NaN cut points (1.5.3) */
#else
#define MK_GRAFT_BISECT_(t_, ta_, tb_) if (!(t_ >= ta_ && t_ <= tb_)) t_ = 0.5 * (ta_ + tb_);
#endif
/* the field F's surface by marching tetrahedra: nodes lo + h (i, j, k), index (k n1 + j) n0 + i (mk_sdf's layout), negative inside.
   Out: V, T malloc'd (closed and outward inside the box). 0, or -1 for a size under 2, more nodes than fit, h not positive and
   finite, a non-finite lo or far corner lo + (n - 1) h, or no memory (out then empty). A NaN node counts as outside. */
static inline int mk_mtet(const float *F, const int n[3], const double lo[3], double h, mk_gmesh *out)
{
    memset(out, 0, sizeof *out);
#if defined(MK_MUTATE) && MK_MUTATE == 230
    if (!F || n[0] < 2 || n[1] < 2 || n[2] < 2 || !(h > 0) || !isfinite(h) || !isfinite(lo[0] + lo[1] + lo[2])) return -1; /* 230: a finite lo's sum can overflow (1.5.3) */
#else
    if (!F || n[0] < 2 || n[1] < 2 || n[2] < 2 || !(h > 0) || !isfinite(h) || !isfinite(lo[0]) || !isfinite(lo[1]) || !isfinite(lo[2])) return -1;
#endif
    if ((size_t)n[0] > (size_t)INT_MAX / (size_t)n[1] / (size_t)n[2] || (size_t)n[0] * n[1] * n[2] > SIZE_MAX / 64) return -1;
#if !(defined(MK_MUTATE) && MK_MUTATE == 231) /* 231: the far corner not looked at (1.5.3) */
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
        for (int k = 0; k + 1 < n[2]; k++)
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
                        const int sw = -s; /* every triangle wound inward */
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
            t_ = !(t_ >= 0.1) ? 0.1 : t_ > 0.9 ? 0.9 : t_;                                                                                        \
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
    long inserted, welded, fanned, dropped, split; /* patch rim vertices put into S's rim edges, landed on its vertices; rim triangles */
    double snap;                         /* the largest move of a patch rim vertex onto S's rim */
    int pieces;                          /* closed pieces before the specks were dropped */
    double dropped_vol, vol0, vol1;      /* their volume; S's volume before and the result's after */
    long dup, lone;                      /* directed edges used twice; without their reverse (both 0 for a closed manifold) */
    int loops_s, loops_p, unpaired;      /* 1.5.2: the rim loops of S's cut and of the patch; loops without a partner (-1: not loops) */
    long merged;                         /* 1.5.2: patch rim vertices merged into a rim neighbour, their place on S's rim out of order */
    int inverted;                        /* 1.5.2: pieces wound against S (their volume of the other sign), the kept one included */
    double inverted_vol;                 /* 1.5.2: their volume, summed (of the other sign) */
} mk_graft_report;

static inline void mk_graft_free(mk_graft *g) { mk_bvh_free(&g->bv), free(g->psi), memset(g, 0, sizeof *g); }
typedef struct { long long k; int m; } mk_graft_ekey_;
static inline int mk_graft_ekcmp_(const void *a, const void *b)
{
    const long long x = ((const mk_graft_ekey_ *)a)->k, y = ((const mk_graft_ekey_ *)b)->k;
    return (x > y) - (x < y);
}
/* 1.5.2 (the maintainer's race-c1 H3 roomsplit refine): S refined where the seam will cut it, before mk_graft_init: every edge longer
   than lmax whose ends' psi reach below band and above -band is split at its midpoint, in both its triangles (1, 2 or 3 split
   edges a triangle by the usual templates, two split edges on the quad's shorter diagonal, so no T-junction), passes until none
   is left or 12 have run. The surface is unchanged (midpoints of flat triangles' edges); only the seam's rim gets as fine as a
   patch's, which a coarse S's long rim edges are not (the maintainer's H3 on 157 muscles, refining as a retry: 200 of 330 rows took,
   184 without). out: the refined mesh (V, T malloc'd, keep NULL; S's vertices first, as they were). The passes run, or -1 for
   a bad mesh, a NULL fn, lmax or band not positive and finite, or no memory (out then empty). */
static inline int mk_graft_refine(const double *V, int nv, const int *T, int nt, mk_graft_psi fn, void *ctx, double lmax, double band, mk_gmesh *out)
{
    memset(out, 0, sizeof *out);
    if (!V || !T || !fn || nv < 3 || nt < 1 || nv > INT_MAX / 4 || nt > INT_MAX / 8 || !(lmax > 0) || !isfinite(lmax) || !(band > 0) || !isfinite(band)) return -1;
    for (size_t q = 0; q < 3 * (size_t)nt; q++)
        if (T[q] < 0 || T[q] >= nv) return -1;
#if !(defined(MK_MUTATE) && MK_MUTATE == 234) /* 234: coordinates not looked at (1.5.3) */
    for (size_t q = 0; q < 3 * (size_t)nv; q++)
        if (!isfinite(V[q])) return -1; /* a bad mesh, as mk_graft_init's (the verification work R-1: NaN or inf passed into the midpoints) */
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
            if (ns == 1) ns = 0; /* a triangle with one split edge left whole: a T-junction (the 1-split template left out) */
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
/* S (closed, V[3 nv], T[3 nt]; kept by pointer, so it must outlive g) and the seam function. 0, or -1 for a bad mesh, a NULL fn,
   near or nudge negative or not finite, or no memory. */
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
/* a surface cut at psi = 0, the side sgn psi > 0 kept: the kept vertices (map), then the cut points; rim edges a -> b with their
   triangle. g NULL: the cut points linear in psi; else bisected on the patch's psi (8 steps). */
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
/* the field's surface (mk_mtet) kept where psi < 0, an open patch: out->keep 1 on its rim. 0, or -1 as mk_mtet, or no memory. */
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
/* S cut at psi = 0 with the patch (PV[3 pnv], PT[3 pnt]: open, its rim on the seam) stitched in. out: the result, keep as above.
   0; -1 for a bad patch, no memory, S's cut rim not closed loops, a patch rim vertex farther than snap from S's rim (or merged
   farther than snap from its rim neighbour), a rim that does not follow S's, or snap negative or not finite; 1 when the result
   is not a closed manifold (rep->dup, rep->lone; a vertex pinched between two fans, as from a pinched patch, is not looked
   for: the verification work G-6), 2 when it is one but wound against S (its volume of the other sign), or when
   pieces wound against S hold over 1 % of its volume (rep->inverted, inverted_vol: part of the body turned inside out by the
   stitch); out is still given with 1 and 2. rep is filled on -1 too, as far as the stitch got;
   rep->loops_s, loops_p and unpaired say how the rims paired up (unpaired -1: not simple loops). */
static inline int mk_graft_stitch(const mk_graft *g, const double *PV, int pnv, const int *PT, int pnt, double snap, mk_gmesh *out, mk_graft_report *rep)
{
    mk_graft_report R = {0};
    memset(out, 0, sizeof *out);
    if (rep) *rep = R;
    if (!g || !g->psi || !PV || !PT || pnv < 3 || pnt < 1 || pnv > INT_MAX / 4 || pnt > INT_MAX / 8 || !(snap >= 0) || !isfinite(snap)) return -1;
    for (size_t q = 0; q < 3 * (size_t)pnv; q++)
        if (!isfinite(PV[q])) return -1;
    for (size_t q = 0; q < 3 * (size_t)pnt; q++)
        if (PT[q] < 0 || PT[q] >= pnv) return -1;
    int ret = -1;
    mk_gcut_ CS = {0}, CP = {0};
    uint64_t *de = NULL, *ed = NULL;
    int *startS = NULL, *nextS = NULL, *pe = NULL, *weld = NULL, *newid = NULL, *OT = NULL, *pr = NULL, *fin = NULL, *rt = NULL, *sl = NULL, *so = NULL, *slen = NULL;
    int *smet = NULL, *startP = NULL, *lv = NULL, *ltl = NULL, *lpre = NULL, *lst = NULL;
    double *pt = NULL, (*ins)[3] = NULL, *OV = NULL, *pv = NULL, *pq = NULL, *pd = NULL, *lw = NULL;
    unsigned char *dead = NULL, *usep = NULL, *usef = NULL, *lkeep = NULL, *pseen = NULL;
    {
        double v = 0;
        for (int t = 0; t < g->nt; t++) {
            const double *a = g->V + 3 * (size_t)g->T[3 * (size_t)t], *b = g->V + 3 * (size_t)g->T[3 * (size_t)t + 1], *c = g->V + 3 * (size_t)g->T[3 * (size_t)t + 2];
            v += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
        }
        R.vol0 = v;
    }
    if (mk_graft_clip_(g->V, g->nv, g->T, g->nt, g->psi, 1, NULL, &CS)) goto out;
    /* the patch with its rim edges: a directed edge without its reverse; a triangle with two rim edges split at its centre */
    for (int i = 0; i < pnv; i++) mk_gcut_v_(&CP, PV + 3 * (size_t)i);
    de = malloc(3 * (size_t)pnt * sizeof *de);
    if (!de || CP.oom) goto out;
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
    if (CP.oom || CS.oom) goto out;
    R.rim_s = CS.ne, R.rim_p = CP.ne, R.patch_nv = CP.nv, R.patch_nt = CP.nt;
    /* S's rim: next[e] = the rim edge that starts where e ends */
    startS = malloc((size_t)(CS.nv ? CS.nv : 1) * sizeof *startS), nextS = malloc((size_t)(CS.ne ? CS.ne : 1) * sizeof *nextS);
    if (!startS || !nextS) goto out;
    for (int i = 0; i < CS.nv; i++) startS[i] = -1;
    long bad = 0;
    for (int e = 0; e < CS.ne; e++) {
        if (startS[CS.E[3 * (size_t)e]] >= 0) bad++;
        startS[CS.E[3 * (size_t)e]] = e;
    }
    for (int e = 0; e < CS.ne; e++)
        if ((nextS[e] = startS[CS.E[3 * (size_t)e + 1]]) < 0) bad++;
    if (bad || !CS.ne || !CP.ne) goto out;
    /* each patch rim vertex onto S's rim: its edge, where along it, S's vertex when it lands on one */
    pe = malloc((size_t)CP.nv * sizeof *pe), weld = malloc((size_t)CP.nv * sizeof *weld), pt = malloc((size_t)CP.nv * sizeof *pt);
    if (!pe || !weld || !pt) goto out;
    pq = malloc(3 * (size_t)CP.nv * sizeof *pq), pd = malloc((size_t)CP.nv * sizeof *pd), rt = malloc((size_t)CP.nv * sizeof *rt);
    if (!pq || !pd || !rt) goto out;
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
        if (be < 0 || !(bd <= snap * snap)) goto out;
        if (bd > snap2) snap2 = bd;
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
    /* 1.5.2: the rims as loops. S's rim edges numbered along each loop, so a place on it is (loop, edge number + t). Each patch rim
       loop should lie on one of S's and each of S's be met by one; the loops that do not (rep->unpaired: a patch loop with no
       counterpart, as a hole in the patch by the seam or an island of it, one on two of S's loops, a loop of S's met by none or by
       two) are counted and placed as 1.5.0 placed them: their pieces are mostly specks the cut leaves, but they can break the body
       (the maintainer's supraspinatus r against the humerus: S 3 loops, the patch 5), which 1 then says. Along a patch loop alone on its
       loop of S's the places must fall strictly: the longest such run, each place measured behind a start (the vertex nearest S's
       rim among those in order with both rim neighbours), is kept, and each vertex out of it (a dense rim's wiggle stepping back on
       S's coarser edges: 1.5.0 then put an edge in twice) is merged into the kept rim neighbour before or after it, moving no
       farther than snap (else -1); a merged vertex's rim triangle is dropped. Rims that are not simple loops: as 1.5.0. */
    {
        sl = malloc((size_t)CS.ne * sizeof *sl), so = malloc((size_t)CS.ne * sizeof *so), slen = malloc((size_t)CS.ne * sizeof *slen);
        smet = calloc((size_t)CS.ne, sizeof *smet), startP = malloc((size_t)CP.nv * sizeof *startP), lv = malloc((size_t)CP.ne * sizeof *lv);
        lw = malloc((size_t)CP.ne * sizeof *lw), ltl = malloc((size_t)CP.ne * sizeof *ltl), lpre = malloc((size_t)CP.ne * sizeof *lpre);
        lkeep = malloc((size_t)CP.ne), pseen = calloc((size_t)CP.ne, 1), lst = malloc(((size_t)CP.ne + 1) * sizeof *lst);
        if (!sl || !so || !slen || !smet || !startP || !lv || !lw || !ltl || !lpre || !lkeep || !pseen || !lst) goto out;
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
#define MK_GR_BACK_(a, b, n) (MK_GR_POS_(a) - MK_GR_POS_(b) < 0 ? MK_GR_POS_(a) - MK_GR_POS_(b) + (n) : MK_GR_POS_(a) - MK_GR_POS_(b)) /* b behind a, in [0, n) */
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
        for (int L = 0; L < R.loops_s; L++) R.unpaired += smet[L] == 0 ? 1 : smet[L] - 1; /* met by none, or by two or more */
#if defined(MK_MUTATE) && MK_MUTATE == 205
        R.unpaired = 0; /* the loops not counted (1.5.0) */
#endif
        for (int q = 0; q < nlp; q++) {
            const int *v = lv + lst[q], k = lst[q + 1] - lst[q], L = sl[pe[v[0]]];
            int one = 1;
            for (int i = 1; i < k; i++) one &= sl[pe[v[i]]] == L;
            if (!one || smet[L] != 1) continue; /* not alone on one loop of S's: as 1.5.0 */
            const double n = slen[L];
            int s = -1; /* the start: nearest S's rim among the vertices in order with both rim neighbours, else nearest */
            for (int i = 0; i < k; i++) {
                const int a = v[(i + k - 1) % k], b = v[i], c = v[(i + 1) % k];
                if (MK_GR_BACK_(a, b, n) < 0.5 * n && MK_GR_BACK_(b, c, n) < 0.5 * n && (s < 0 || pd[b] < pd[v[s]])) s = i;
            }
            for (int i = 0; s < 0 && i < k; i++)
                if (s < 0 || pd[v[i]] < pd[v[s]]) s = i;
            /* each place as how far behind the start it lies, in (-n, 0]: the longest strictly falling run goes once round */
            int nl = 0;
            lw[0] = 0, lkeep[0] = 1;
            for (int j = 1; j < k; j++) {
                lw[j] = -MK_GR_BACK_(v[s], v[(s + j) % k], n), lkeep[j] = 0, lpre[j] = -1;
#if defined(MK_MUTATE) && MK_MUTATE == 202
                lkeep[j] = 1; /* every vertex kept, in or out of order (1.5.0) */
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
            /* each run out of order between kept a (before) and b (after): its first part into a, the rest into b, by which is nearer */
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
                    if (!(d2 <= snap * snap)) goto out;
                    if (weld[w] >= 0) R.welded--;
                    else R.inserted--;
                    rt[w] = r, pe[w] = pe[r], pt[w] = pt[r], weld[w] = -1, R.merged++;
                }
                j = j2;
            }
        }
#undef MK_GR_POS_
#undef MK_GR_BACK_
        for (int p = 0; p < CP.nv; p++)
            if (pe[p] >= 0 && rt[p] == p && weld[p] < 0) memcpy(CP.V + 3 * (size_t)p, pq + 3 * (size_t)p, 3 * sizeof *pq);
    }
#define MK_GR_PID_(p) (weld[rt[p]] >= 0 ? weld[rt[p]] : CS.nv + rt[p])
    /* S's rim triangles fanned over the patch vertices put into their rim edge */
    {
        ins = malloc((size_t)(R.inserted ? R.inserted : 1) * sizeof *ins);
        if (!ins) goto out;
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
            if (r == 3) goto out;
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
        if (CS.oom) goto out;
    }
    /* the patch's rim triangles fanned over S's rim vertices between their two ends (S's vertices carried as -1 - id) */
    const int ncpt = CP.nt;
    dead = calloc((size_t)ncpt + 1, 1);
    if (!dead) goto out;
    for (int e2 = 0; e2 < CP.ne; e2++) {
        const int p = CP.E[3 * (size_t)e2], p2 = CP.E[3 * (size_t)e2 + 1], t = CP.E[3 * (size_t)e2 + 2];
        if (pe[p] < 0 || pe[p2] < 0) goto out;
        int *f = CP.T + 3 * (size_t)t, r = 0;
        while (r < 3 && !(f[r] == p && f[(r + 1) % 3] == p2)) r++;
        if (r == 3) goto out;
        const int x = f[(r + 2) % 3];
        if (MK_GR_PID_(p) == MK_GR_PID_(p2)) { dead[t] = 1, R.dropped++; continue; }
        int mid[64], nm = 0, cur = pe[p2];
        if (pe[p2] != pe[p]) {
            for (;;) {
                if (nm == 64) goto out;
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
        if (CP.oom) goto out;
    }
    /* the pieces: triangles with final numbers; the piece of most volume kept */
    const int npv = CP.nv;
    newid = malloc((size_t)npv * sizeof *newid), usep = calloc((size_t)npv + 1, 1);
    if (!newid || !usep) goto out;
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
    if (!OV || !OT) goto out;
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
        if (t >= ncpt) continue; /* the patch's fans left out */
#endif
        MK_GR_EMIT_(MK_GR_QID_(CP.T[3 * (size_t)t]), MK_GR_QID_(CP.T[3 * (size_t)t + 1]), MK_GR_QID_(CP.T[3 * (size_t)t + 2]));
    }
    pr = malloc((size_t)no * sizeof *pr), pv = calloc((size_t)no, sizeof *pv);
    if (!pr || !pv) goto out;
    for (int i = 0; i < no; i++) pr[i] = i;
    for (long t = 0; t < nto; t++) {
        pr[mk_graft_find_(pr, OT[3 * t + 1])] = mk_graft_find_(pr, OT[3 * t]);
        pr[mk_graft_find_(pr, OT[3 * t + 2])] = mk_graft_find_(pr, OT[3 * t]);
    }
    for (long t = 0; t < nto; t++) {
        const double *a = OV + 3 * (size_t)OT[3 * t], *b = OV + 3 * (size_t)OT[3 * t + 1], *c = OV + 3 * (size_t)OT[3 * t + 2];
        pv[mk_graft_find_(pr, OT[3 * t])] += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
    }
    /* 1.5.2: the piece of most volume whatever its sign (1.5.0 took the most positive, so a body come out inverted or broken, its
       volume negative, lost to a closed speck and passed: the maintainer's supraspinatus r became a 42-face speck); a piece wound against S
       counted, and 2 returned when the kept one is */
    int big = -1;
    const double sg = R.vol0 < 0 ? -1 : 1;
    for (int i = 0; i < no; i++)
        if (pr[i] == i && pv[i] != 0) {
            R.pieces++, R.inverted += sg * pv[i] < 0, R.inverted_vol += sg * pv[i] < 0 ? pv[i] : 0;
#if defined(MK_MUTATE) && MK_MUTATE == 203
            if (big < 0 || pv[i] > pv[big]) big = i; /* the most positive volume (1.5.0) */
#else
            if (big < 0 || fabs(pv[i]) > fabs(pv[big])) big = i;
#endif
        }
    for (int i = 0; i < no; i++)
        if (pr[i] == i && i != big) R.dropped_vol += pv[i];
    if (big < 0) goto out;
    fin = malloc((size_t)no * sizeof *fin), usef = calloc((size_t)no + 1, 1);
    if (!fin || !usef) goto out;
    for (long t = 0; t < nto; t++)
        if (mk_graft_find_(pr, OT[3 * t]) == big) usef[OT[3 * t]] = usef[OT[3 * t + 1]] = usef[OT[3 * t + 2]] = 1;
    int nfin = 0;
    for (int i = 0; i < no; i++) fin[i] = usef[i] ? nfin++ : -1;
    out->V = malloc(3 * (size_t)(nfin ? nfin : 1) * sizeof *out->V), out->keep = malloc((size_t)(nfin ? nfin : 1)), out->T = malloc(3 * (size_t)(nto ? nto : 1) * sizeof *out->T);
    ed = malloc(3 * (size_t)(nto ? nto : 1) * sizeof *ed);
    if (!out->V || !out->keep || !out->T || !ed) goto out;
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
        ed[ned++] = (uint64_t)a << 32 | (uint32_t)b, ed[ned++] = (uint64_t)b << 32 | (uint32_t)c, ed[ned++] = (uint64_t)c << 32 | (uint32_t)a;
    }
    out->nv = nfin, R.vol1 = pv[big];
    qsort(ed, ned, sizeof *ed, mk_graft_cmpu_);
    for (size_t i = 0; i < ned; i++) {
        if (i && ed[i] == ed[i - 1]) R.dup++;
        const uint64_t rv = (ed[i] & 0xffffffffull) << 32 | ed[i] >> 32;
        if (!bsearch(&rv, ed, ned, sizeof *ed, mk_graft_cmpu_)) R.lone++;
    }
#if defined(MK_MUTATE) && MK_MUTATE == 204
    ret = R.dup || R.lone ? 1 : 0; /* a kept piece wound against S passes (1.5.0) */
#elif defined(MK_MUTATE) && MK_MUTATE == 210
    ret = R.dup || R.lone ? 1 : sg * R.vol1 < 0 ? 2 : 0; /* only the kept piece judged */
#else
    /* 2: the kept piece wound against S, or pieces wound against S holding over 1 % of its volume: the stitch turned part of the
       body inside out (MEASURED on the maintainer's H3 muscles: 1 to 70 % where the body broke, under 0.1 % for the specks of rounding
       a clean stitch leaves) */
    ret = R.dup || R.lone ? 1 : sg * R.vol1 < 0 || fabs(R.inverted_vol) > 1e-2 * fabs(R.vol1) ? 2 : 0;
#endif
#undef MK_GR_EMIT_
#undef MK_GR_SID_
#undef MK_GR_QID_
#undef MK_GR_PID_
out:
    if (ret < 0) mk_gmesh_free(out);
    if (rep) *rep = R;
    mk_gcut_free_(&CS), mk_gcut_free_(&CP);
    free(de), free(ed), free(startS), free(nextS), free(pe), free(weld), free(newid), free(OT), free(pr), free(fin), free(pt), free(ins), free(OV), free(pv);
    free(dead), free(usep), free(usef), free(rt), free(sl), free(so), free(slen), free(smet), free(startP), free(lv), free(ltl), free(lpre);
    free(pq), free(pd), free(lw), free(lkeep), free(pseen), free(lst);
    return ret;
}
#endif
