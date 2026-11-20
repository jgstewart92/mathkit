/* mathkit: surface nets on a labelled voxel grid (Gibson 1998, "Constrained elastic surface nets"; s-forms, 2026-11-15, with
   the house's other copies merged in). One vertex in each dual cell (the cube between eight voxel centres) whose voxels do not
   all carry one label, shared by every label's surface there, so two labels that touch share their interface's vertices
   exactly (their triangles there are the same, wound oppositely). A vertex starts at the mean of the crossings on the cell's
   label-changing edges: the edge midpoints for labels, or for a scalar field (mk_snet_build_iso) the linear crossing of iso,
   held to [0.05, 0.95] of the edge (body-claude mesh.c: a sample exactly at iso would put the crossings of every edge round it
   at one point, giving zero-length edges). Then it relaxes by Taubin's shrink-free pair (lambda 0.5, mu -0.53) toward its six
   dual neighbours, held inside its own dual cell, which keeps the topology and the voxels' fidelity. A label's triangles: each
   voxel face between it and another label, as the quad of the four dual vertices round it, split on its shorter diagonal,
   wound anticlockwise seen from outside the label. No case table, nothing ambiguous.
   Grid: lab[(k ny + j) nx + i], voxel (i, j, k) spanning org + h [i, i+1] x [j, j+1] x [k, k+1]; outside the grid is label 0.
   For a scalar field phi at the voxel centres, label 1 is phi < iso (inside) and 0 the rest.
   Memory goes with the surface, not the grid (as experiment/cfd wrap.c keeps two layers): each vertex keeps its dual cell's
   index, ascending, and the faces are found through two rolling layers of the dual grid, so a 512^3 grid holding one ball
   needs no 4-byte-per-cell table. Results equal s-forms' full-table version bit for bit.
   Where a label's voxels in one dual cell (or the others) fall in two groups not joined by cube edges, its surface touches
   itself there (a non-manifold vertex or edge); mk_snet_split counts such cells (experiment/cfd wrap.c's test). Header only,
   C23. */
#ifndef MK_SNET_H
#define MK_SNET_H
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int nx, ny, nz;  /* the voxel grid; dual cells are (nx + 1) (ny + 1) (nz + 1), dual cell (i, j, k) between voxels i-1..i */
    double org[3], h;
    double *V;       /* 3 nv coordinates */
    int nv;
    int64_t *cell;   /* each vertex's dual cell, (k (ny + 1) + j) (nx + 1) + i, ascending */
    int *zs;         /* the vertices of dual layer k are zs[k] .. zs[k + 1] - 1 (nz + 2 entries) */
} mk_snet;

/* the label of voxel (i, j, k): lab's, or phi < iso; 0 outside the grid */
static inline int mk_snet_at_(const uint8_t *lab, const float *phi, float iso, int nx, int ny, int nz, int i, int j, int k)
{
    if (i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz) return 0;
    const size_t q = ((size_t)k * (size_t)ny + (size_t)j) * (size_t)nx + (size_t)i;
    return lab ? lab[q] : phi[q] < iso;
}
static inline void mk_snet_free(mk_snet *s) { free(s->V), free(s->cell), free(s->zs), memset(s, 0, sizeof *s); }
/* dual layer z's vertex ids into buf ((nx + 1) (ny + 1) ints, -1 where none) */
static inline void mk_snet_layer_(const mk_snet *s, int z, int *buf)
{
    const size_t lay = (size_t)(s->nx + 1) * (size_t)(s->ny + 1);
    for (size_t q = 0; q < lay; q++) buf[q] = -1;
    if (z < 0 || z > s->nz) return;
    for (int v = s->zs[z]; v < s->zs[z + 1]; v++) buf[s->cell[v] - (int64_t)z * (int64_t)lay] = v;
}
/* the vertex of dual cell (i, j, k), or -1 */
static inline int mk_snet_vertex(const mk_snet *s, int i, int j, int k)
{
    if (!s->zs || i < 0 || j < 0 || k < 0 || i > s->nx || j > s->ny || k > s->nz) return -1;
    const int64_t c = ((int64_t)k * (s->ny + 1) + j) * (s->nx + 1) + i;
    int lo = s->zs[k], hi = s->zs[k + 1] - 1;
    while (lo <= hi) {
        const int m = lo + (hi - lo) / 2;
        if (s->cell[m] == c) return m;
        if (s->cell[m] < c) lo = m + 1;
        else hi = m - 1;
    }
    return -1;
}

static inline int mk_snet_build_(const uint8_t *lab, const float *phi, float iso, int nx, int ny, int nz, const double org[3], double h, int iters, mk_snet *s)
{
    memset(s, 0, sizeof *s);
    if (nx < 1 || ny < 1 || nz < 1 || nx == INT_MAX || ny == INT_MAX || nz == INT_MAX || !(h > 0) || !isfinite(h) || iters < 0) return -1;
    for (int a = 0; a < 3; a++) if (!isfinite(org[a])) return -1;
    const int dx = nx + 1, dy = ny + 1, dz = nz + 1;
    const size_t lay = (size_t)dx * (size_t)dy;
    if (lay > (size_t)INT64_MAX / (size_t)dz || lay > SIZE_MAX / sizeof(int) / 3 || (size_t)nx * ny > SIZE_MAX / (size_t)nz) return -1;
    if (phi) {
        if (!isfinite(iso)) return -1;
        const size_t nn = (size_t)nx * ny * nz;
        for (size_t q = 0; q < nn; q++) if (!isfinite(phi[q])) return -1;
    }
    s->nx = nx, s->ny = ny, s->nz = nz, s->h = h;
    for (int a = 0; a < 3; a++) s->org[a] = org[a];
    size_t cap = 1 << 16;
    s->V = malloc(sizeof *s->V * 3 * cap), s->cell = malloc(sizeof *s->cell * cap), s->zs = malloc(sizeof *s->zs * ((size_t)dz + 1));
    int *L3 = malloc(sizeof *L3 * 3 * lay), *nb = NULL;
    double *W = NULL;
    if (!s->V || !s->cell || !s->zs || !L3) goto fail;
    size_t nv = 0;
    for (int k = 0; k < dz; k++) {
        s->zs[k] = (int)nv;
        for (int j = 0; j < dy; j++)
            for (int i = 0; i < dx; i++) {
                /* the eight voxels of dual cell (i, j, k): i-1..i, j-1..j, k-1..k */
                int L[8], mixed = 0;
                for (int c = 0; c < 8; c++) L[c] = mk_snet_at_(lab, phi, iso, nx, ny, nz, i - 1 + (c & 1), j - 1 + ((c >> 1) & 1), k - 1 + ((c >> 2) & 1)), mixed |= L[c] != L[0];
                if (!mixed) continue;
                /* the mean of the crossings on the edges whose two voxels differ, in units of h from the cell's low voxel centre */
                double p[3] = {0, 0, 0};
                int ne = 0;
                for (int c = 0; c < 8; c++)
                    for (int ax = 0; ax < 3; ax++) {
                        if (c & (1 << ax)) continue;
                        const int d = c | (1 << ax);
                        if (L[c] == L[d]) continue;
                        double t = 0.5;
                        if (phi) {
                            const int ic = i - 1 + (c & 1), jc = j - 1 + ((c >> 1) & 1), kc = k - 1 + ((c >> 2) & 1);
                            const int id = ic + (ax == 0), jd = jc + (ax == 1), kd = kc + (ax == 2);
                            if (ic >= 0 && jc >= 0 && kc >= 0 && id < nx && jd < ny && kd < nz) { /* both samples in the grid */
                                const double fa = phi[((size_t)kc * ny + jc) * nx + ic], fb = phi[((size_t)kd * ny + jd) * nx + id];
#if defined(MK_MUTATE) && MK_MUTATE == 124
                                (void)fa, (void)fb;
#else
                                t = (iso - fa) / (fb - fa);
                                t = t < 0.05 ? 0.05 : t > 0.95 ? 0.95 : t;
#endif
                            }
                        }
                        const double m[3] = {(c & 1) + (ax == 0 ? t : 0), ((c >> 1) & 1) + (ax == 1 ? t : 0), ((c >> 2) & 1) + (ax == 2 ? t : 0)};
                        for (int a = 0; a < 3; a++) p[a] += m[a];
                        ne++;
                    }
                if (nv == (size_t)INT_MAX) goto fail;
                if (nv == cap) {
                    cap *= 2;
                    double *v2 = realloc(s->V, sizeof *s->V * 3 * cap);
                    if (!v2) goto fail;
                    s->V = v2;
                    int64_t *c2 = realloc(s->cell, sizeof *s->cell * cap);
                    if (!c2) goto fail;
                    s->cell = c2;
                }
                const int ci[3] = {i, j, k};
                for (int a = 0; a < 3; a++) s->V[3 * nv + a] = org[a] + h * (ci[a] - 0.5 + p[a] / ne);
                s->cell[nv++] = ((int64_t)k * dy + j) * dx + i;
            }
    }
    s->zs[dz] = (int)nv, s->nv = (int)nv;
    /* relax (Taubin 1995: a step of lambda 0.5 toward the six dual neighbours' mean, then one of mu -0.53 away, so a convex part
       does not shrink as under plain averaging), each vertex held inside its own dual cell; one relaxation is both steps. The
       neighbours, +x -x +y -y +z -z, come from three rolling dual layers. */
    if (iters > 0 && nv) {
        nb = malloc(sizeof *nb * 6 * nv), W = malloc(sizeof *W * 3 * nv);
        if (!nb || !W) goto fail;
        int *Lm = L3, *L0 = L3 + lay, *Lp = L3 + 2 * lay;
        mk_snet_layer_(s, -1, Lm), mk_snet_layer_(s, 0, L0);
        for (int k = 0; k < dz; k++) {
            mk_snet_layer_(s, k + 1, Lp);
            for (int v = s->zs[k]; v < s->zs[k + 1]; v++) {
                const int64_t r = s->cell[v] - (int64_t)k * (int64_t)lay;
                const int i = (int)(r % dx), j = (int)(r / dx);
                nb[6 * (size_t)v] = i + 1 < dx ? L0[r + 1] : -1, nb[6 * (size_t)v + 1] = i > 0 ? L0[r - 1] : -1;
                nb[6 * (size_t)v + 2] = j + 1 < dy ? L0[r + dx] : -1, nb[6 * (size_t)v + 3] = j > 0 ? L0[r - dx] : -1;
                nb[6 * (size_t)v + 4] = Lp[r], nb[6 * (size_t)v + 5] = Lm[r];
            }
            int *t = Lm;
            Lm = L0, L0 = Lp, Lp = t;
        }
        for (int it = 0; it < 2 * iters; it++) {
            const double f = (it & 1) ? -0.53 : 0.5;
            for (size_t v = 0; v < nv; v++) {
                double q[3] = {0, 0, 0};
                int n = 0;
                for (int e = 0; e < 6; e++) {
                    const int u = nb[6 * v + e];
                    if (u < 0) continue;
                    for (int t = 0; t < 3; t++) q[t] += s->V[3 * (size_t)u + t];
                    n++;
                }
                const int64_t c = s->cell[v];
                const int ci[3] = {(int)(c % dx), (int)(c / dx % dy), (int)(c / ((int64_t)dx * dy))};
                (void)ci; /* unused under MK_MUTATE 121; not [[maybe_unused]], which clightgen refuses (the verification work N-B) */
                for (int t = 0; t < 3; t++) {
                    const double x0 = s->V[3 * v + t];
                    double x = n ? x0 + f * (q[t] / n - x0) : x0;
#if !defined(MK_MUTATE) || MK_MUTATE != 121
                    const double lo = org[t] + h * (ci[t] - 0.5), hi = lo + h; /* the dual cell */
                    x = x < lo ? lo : x > hi ? hi : x;
#endif
                    W[3 * v + t] = x;
                }
            }
            memcpy(s->V, W, sizeof *W * 3 * nv);
        }
    }
    free(nb), free(W), free(L3);
    return 0;
fail:
    free(nb), free(W), free(L3), mk_snet_free(s);
    return -1;
}
/* the shared vertices of a labelled grid, relaxed `iters` times (0 leaves them at the crossing means). 0, or -1 for a bad
   size, h or org, iters < 0, or no memory (*s is then empty). */
static inline int mk_snet_build(const uint8_t *lab, int nx, int ny, int nz, const double org[3], double h, int iters, mk_snet *s)
{
    if (!lab) { memset(s, 0, sizeof *s); return -1; }
    return mk_snet_build_(lab, NULL, 0, nx, ny, nz, org, h, iters, s);
}
/* the same for a scalar field phi at the voxel centres, inside where phi < iso; -1 also for a non-finite phi or iso */
static inline int mk_snet_build_iso(const float *phi, int nx, int ny, int nz, const double org[3], double h, float iso, int iters, mk_snet *s)
{
    if (!phi) { memset(s, 0, sizeof *s); return -1; }
    return mk_snet_build_(NULL, phi, iso, nx, ny, nz, org, h, iters, s);
}

static inline int mk_snet_faces_(const mk_snet *s, const uint8_t *lab, const float *phi, float iso, int want, int **T)
{
    const int nx = s->nx, ny = s->ny, nz = s->nz, dx = nx + 1, dy = ny + 1;
    *T = NULL;
    if (!s->zs) return -1;
    const size_t lay = (size_t)dx * dy;
    size_t cap = 1 << 12, nt = 0;
    int *t = malloc(sizeof *t * 3 * cap), *LL = malloc(sizeof *LL * 2 * lay);
    if (!t || !LL) { free(t), free(LL); return -1; }
    int *Lk = LL, *Lk1 = LL + lay; /* dual layers k and k + 1 */
    mk_snet_layer_(s, -1, Lk), mk_snet_layer_(s, 0, Lk1);
    for (int k = -1; k < nz; k++) {
        if (k >= 0) {
            int *x = Lk;
            Lk = Lk1, Lk1 = x;
            mk_snet_layer_(s, k + 1, Lk1);
        }
        for (int j = -1; j < ny; j++)
            for (int i = -1; i < nx; i++) {
                const int l0 = mk_snet_at_(lab, phi, iso, nx, ny, nz, i, j, k);
#if defined(MK_MUTATE) && MK_MUTATE == 123
                for (int ax = 0; ax < 2; ax++) {
#else
                for (int ax = 0; ax < 3; ax++) {
#endif
                    const int i1 = i + (ax == 0), j1 = j + (ax == 1), k1 = k + (ax == 2);
                    if (i1 >= nx + 1 || j1 >= ny + 1 || k1 >= nz + 1) continue;
                    const int l1 = mk_snet_at_(lab, phi, iso, nx, ny, nz, i1, j1, k1);
                    if (l0 == l1 || (l0 != want && l1 != want)) continue;
                    /* the face between voxel (i, j, k) and its +ax neighbour: the four dual cells round the edge joining their
                       centres, which share the dual coordinate i1 (j1, k1) along ax and take both values in the other two axes */
                    const int u = (ax + 1) % 3, w = (ax + 2) % 3;
                    int q[4];
                    for (int c = 0; c < 4; c++) {
                        int d[3] = {i + 1, j + 1, k + 1};
                        d[ax] = (ax == 0 ? i1 : ax == 1 ? j1 : k1);
                        d[u] -= (c == 1 || c == 2) ? 0 : 1, d[w] -= (c >= 2) ? 0 : 1;
                        q[c] = (d[2] == k + 1 ? Lk1 : Lk)[(size_t)d[1] * dx + d[0]];
                        if (q[c] < 0) { free(t), free(LL); return -1; } /* 1.6.1: not s's grid (1.6.0 read V[-3]) */
                    }
                    /* q[0..3] go anticlockwise about +ax; outward from the label is +ax when l0 is the label, -ax otherwise */
                    int out = l0 == want;
#if defined(MK_MUTATE) && MK_MUTATE == 122
                    out = !out;
#endif
                    if (!out) { const int x = q[1]; q[1] = q[3], q[3] = x; }
                    const double *A = s->V + 3 * (size_t)q[0], *B = s->V + 3 * (size_t)q[1], *C = s->V + 3 * (size_t)q[2], *D = s->V + 3 * (size_t)q[3];
                    double d02 = 0, d13 = 0;
                    for (int a = 0; a < 3; a++) d02 += (A[a] - C[a]) * (A[a] - C[a]), d13 += (B[a] - D[a]) * (B[a] - D[a]);
                    int tri[6];
                    if (d02 <= d13) tri[0] = q[0], tri[1] = q[1], tri[2] = q[2], tri[3] = q[0], tri[4] = q[2], tri[5] = q[3];
                    else tri[0] = q[0], tri[1] = q[1], tri[2] = q[3], tri[3] = q[1], tri[4] = q[2], tri[5] = q[3];
                    if (nt + 2 > cap) {
                        if (cap > (size_t)INT_MAX / 2) { free(t), free(LL); return -1; }
                        int *t2 = realloc(t, sizeof *t * 3 * (cap *= 2));
                        if (!t2) { free(t), free(LL); return -1; }
                        t = t2;
                    }
                    memcpy(t + 3 * nt, tri, sizeof tri), nt += 2;
                }
            }
    }
    free(LL);
    *T = t;
    return (int)nt;
}
/* the triangles of label `want`'s boundary (indices into s->V), anticlockwise seen from outside the label. *T is allocated
   (free it); the triangle count, or -1 for no memory, more than INT_MAX triangles, or a lab other than the grid s was built from
   where a face needs a vertex s has not (1.6.1). */
static inline int mk_snet_faces(const mk_snet *s, const uint8_t *lab, int want, int **T) { return mk_snet_faces_(s, lab, NULL, 0, want, T); }
/* the inside's boundary (phi < iso) of a net from mk_snet_build_iso */
static inline int mk_snet_faces_iso(const mk_snet *s, const float *phi, float iso, int **T) { return mk_snet_faces_(s, NULL, phi, iso, 1, T); }

/* 1 if the corners of a 2 x 2 x 2 cell with in[c] set fall in two or more groups joined by cube edges */
static inline int mk_snet_split8_(const int in[8])
{
    int n = 0, seen[8] = {0}, st[8], sp = 0, cnt = 0;
    for (int c = 0; c < 8; c++) n += in[c];
    if (n < 2) return 0;
    for (int c = 0; c < 8; c++) if (in[c]) { st[sp++] = c, seen[c] = 1; break; }
    while (sp) {
        const int c = st[--sp];
        cnt++;
        for (int b = 0; b < 3; b++) {
            const int r = c ^ (1 << b);
            if (in[r] && !seen[r]) seen[r] = 1, st[sp++] = r;
        }
    }
    return cnt != n;
}
static inline int mk_snet_split_(const mk_snet *s, const uint8_t *lab, const float *phi, float iso, int want)
{
    const int dx = s->nx + 1, dy = s->ny + 1;
    int bad = 0;
    for (int v = 0; v < s->nv; v++) {
        const int64_t c = s->cell[v];
        const int i = (int)(c % dx), j = (int)(c / dx % dy), k = (int)(c / ((int64_t)dx * dy));
        int a[8], b[8];
        for (int q = 0; q < 8; q++) a[q] = mk_snet_at_(lab, phi, iso, s->nx, s->ny, s->nz, i - 1 + (q & 1), j - 1 + ((q >> 1) & 1), k - 1 + ((q >> 2) & 1)) == want, b[q] = !a[q];
#if defined(MK_MUTATE) && MK_MUTATE == 125
        (void)b;
#else
        bad += mk_snet_split8_(a) || mk_snet_split8_(b);
#endif
    }
    return bad;
}
/* the dual cells where label want's voxels, or the others, fall in two or more groups not joined by cube edges: where want's
   surface touches itself. 0 means every vertex of want's surface has one disc of triangles round it. */
static inline int mk_snet_split(const mk_snet *s, const uint8_t *lab, int want) { return mk_snet_split_(s, lab, NULL, 0, want); }
static inline int mk_snet_split_iso(const mk_snet *s, const float *phi, float iso) { return mk_snet_split_(s, NULL, phi, iso, 1); }
#endif
