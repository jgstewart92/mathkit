/* mathkit: a triangle mesh to a signed distance grid, negative inside (s-bone's meshsdf, Batty's SDFGen 2007, with the
   house's other copies merged in; cudakit's sd_ primitives use the same sign). Double-precision mesh in mk_geom's layout
   (V[3 nv], T[3 nt]); float distances on the nodes lo + h (i, j, k), index (k n1 + j) n0 + i. Stored as float:
   a distance beyond FLT_MAX (3.4e38) reads inf, and one is rounded to 24 bits (the verification work S-E).
   Distance: the exact distance to the nearest triangle (mk_geom's mk_closest_pt_tri, Ericson 5.1.5) on the nodes within two
   cells of each triangle, then the nearest triangle passed on in two rounds of eight sweeps, each node keeping the nearest
   of its seven sweep-neighbours' triangles. Far from the mesh this is the distance to a triangle that is nearest or nearly so
   (within about 0.05 h on the tests); MK_SDF_EXACT then makes every node exact through mk_geom's BVH, the sweeps' value
   bounding each search.
   Sign: the parity of a ray's crossings along the first axis through each row of nodes, the rows set off a fraction of h
   from the nodes so no ray runs through a node, and the crossings bucketed by row and sorted along it (body-claude's sdfgrid
   tested every triangle for every row). A closed mesh crosses every whole row an even number of times; a row that crosses
   an odd number (a ray through an edge, a hole, a non-manifold edge) cannot be trusted. It is cast again on two nearby rays
   (O(triangles) each); if both are odd too, its nodes take the sign of the generalized winding number (mk_geom's
   mk_winding_number: |w| > 1/2 is inside), which also covers a mesh that is not watertight, at O(triangles) for each node of
   such a row. A node nearer the mesh than the rays' offset (2.7e-4 h; the test is 4e-4 h) may lie on the other side of its
   ray: it gets a ray of its own, set off by 1/50 of its distance, so every point between is on its side (the winding number
   if that ray is odd). The counts are reported. MK_SDF_FILL (sdfgrid's FILL, experiment/cfd levelset.h's flood): outside
   nodes that no path of outside nodes joins to the grid's faces are made inside (pockets of a surface that folds through itself, or
   covers a region twice; it fills a real cavity too, so it is for solids without one). MK_SDF_LATTICE puts the nodes on the
   world lattice h Z^3 (metal/ref/hair hairgen.c), so grids of one mesh with other margins, or of a cut of it, share nodes.
   mk_sdf_mesh_box (1.5.2, the maintainer's race-c1 H3 sdfbox): the same grid on a box's nodes only, a local field on a large mesh; on
   the lattice its nodes are the whole grid's, and with MK_SDF_EXACT their values too, bit for bit.
   Threads (mk_par.h, with -fopenmp), each giving 1.4.0's grid bit for bit at any thread count: the band by tiles of 8 x 8
   rows, each tile taking its triangles in index order, so every node sees its triangles in the serial order; the sweeps as a
   wavefront of 4 x 4 tiles of rows (a row needs only the rows before it in the sweep, which are done); MK_SDF_EXACT, the odd rows'
   winding numbers and the near nodes node by node; the rows' parity row by row. The crossings' sort and MK_SDF_FILL's
   flood are serial. Header only, C23. */
#ifndef MK_SDF_H
#define MK_SDF_H
#if defined(MK_MUTATE) && MK_MUTATE == 145
#define MK_SDF_MUT145_ 1 /* the near nodes keep the parity's sign */
#else
#define MK_SDF_MUT145_ 0
#endif
#include "mk_geom.h"
#include "mk_par.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { MK_SDF_FILL = 1, MK_SDF_EXACT = 2, MK_SDF_LATTICE = 4 };
typedef struct {
    int n[3];        /* nodes along each axis */
    double lo[3], h; /* node (0, 0, 0), the spacing */
    double in[6];    /* the mesh's box: lo x y z, hi x y z */
    float *d;        /* n0 n1 n2 signed distances (malloc'd: mk_sdf_free, or take it and free it yourself) */
    long odd_rows;   /* rows odd on three nearby rays (a hole, a non-manifold edge): signed by the winding number */
    long reray;      /* rows odd on the first ray but even on a nearby one (a ray through an edge): signed by that one */
    long filled;     /* MK_SDF_FILL: enclosed outside nodes made inside */
    long near;       /* nodes nearer the mesh than 4e-4 h, signed by a ray of their own (or the winding number) */
} mk_sdf;

static inline void mk_sdf_free(mk_sdf *g) { free(g->d), memset(g, 0, sizeof *g); }
typedef struct { int64_t row; double x; } mk_sdf_cross_;
static inline int mk_sdf_cmp_(const void *a, const void *b)
{
    const mk_sdf_cross_ *x = a, *y = b;
    return x->row != y->row ? (x->row < y->row ? -1 : 1) : (x->x > y->x) - (x->x < y->x);
}
/* the crossings of the ray along x through (y, z) with every triangle, ascending in *xs (grown as needed); their count, or -1
   for no memory. O(triangles). */
typedef struct { double x; int s; } mk_sdf_x_; /* a crossing along a row: where, and the sign of the face normal's x */
static inline int mk_sdf_dcmp_(const void *a, const void *b) { const double x = ((const mk_sdf_x_ *)a)->x, y = ((const mk_sdf_x_ *)b)->x; return (x > y) - (x < y); }
static inline long mk_sdf_rowx_(const double *V, const int *T, int nt, double y, double z, mk_sdf_x_ **xs, size_t *cap)
{
    long m = 0;
    for (int f = 0; f < nt; f++) {
        const double *A = V + 3 * (size_t)T[3 * (size_t)f], *B = V + 3 * (size_t)T[3 * (size_t)f + 1], *C = V + 3 * (size_t)T[3 * (size_t)f + 2];
        if ((A[1] < y && B[1] < y && C[1] < y) || (A[1] > y && B[1] > y && C[1] > y) || (A[2] < z && B[2] < z && C[2] < z) || (A[2] > z && B[2] > z && C[2] > z)) continue;
        const double ay = A[1] - y, az = A[2] - z, by = B[1] - y, bz = B[2] - z, cy = C[1] - y, cz = C[2] - z;
        const double w0 = by * cz - bz * cy, w1 = cy * az - cz * ay, w2 = ay * bz - az * by;
        if (!((w0 >= 0 && w1 >= 0 && w2 >= 0) || (w0 <= 0 && w1 <= 0 && w2 <= 0))) continue;
        const double s = w0 + w1 + w2;
        if (s == 0) continue;
        if ((size_t)m == *cap) {
            const size_t c2 = *cap ? 2 * *cap : 64;
            mk_sdf_x_ *x2 = realloc(*xs, sizeof **xs * c2);
            if (!x2) return -1;
            *xs = x2, *cap = c2;
        }
        (*xs)[m].x = (w0 * A[0] + w1 * B[0] + w2 * C[0]) / s, (*xs)[m++].s = s > 0 ? 1 : -1;
    }
    if (m > 1) qsort(*xs, (size_t)m, sizeof **xs, mk_sdf_dcmp_);
    return m;
}
static inline double mk_sdf_tri_(const double *V, const int *T, int f, const double *p)
{
    const double *a = V + 3 * (size_t)T[3 * (size_t)f], *b = V + 3 * (size_t)T[3 * (size_t)f + 1], *c = V + 3 * (size_t)T[3 * (size_t)f + 2];
    double w[3], q[3];
#if defined(MK_MUTATE) && MK_MUTATE == 143
    return mk_closest_pt_tri(p, a, b, c, w, q); /* the squared distance */
#else
    return sqrt(mk_closest_pt_tri(p, a, b, c, w, q));
#endif
}

/* the grid's node places, for the block functions the threads run */
typedef struct { int n[3], lat; double lo[3], bi[3], h; } mk_sdf_lat_;
static inline double mk_sdf_p_(const mk_sdf_lat_ *L, int a, long i) { return L->lat ? (L->bi[a] + (double)i) * L->h : L->lo[a] + (double)i * L->h; }
static inline size_t mk_sdf_id_(const mk_sdf_lat_ *L, long i, long j, long k) { return (size_t)i + (size_t)L->n[0] * ((size_t)j + (size_t)L->n[1] * (size_t)k); }
/* x as an int clamped to [lo, hi] before the conversion (NaN: lo). the verification work S-2: in a box's grid (mk_sdf_mesh_box) a triangle can lie
   more than INT_MAX cells from node 0, and (int) of such a cell index is undefined (x86: INT_MIN), which lost its band and its
   crossings, so whole rows took the wrong sign. Within the grid it is the plain conversion, so every grid is as before. */
static inline int mk_sdf_ci_(double x, int lo, int hi)
{
#if defined(MK_MUTATE) && MK_MUTATE == 237
    if (!(fabs(x) < 2147483647.0)) return lo; /* 237: a cell past an int is lost, as the cast's INT_MIN lost it (1.5.3) */
#endif
    return !(x >= lo) ? lo : x > hi ? hi : (int)x;
}
/* triangle f's nodes within two cells: [i0, i1] on each axis, clamped to the grid (i0 > i1: it reaches no node) */
static inline void mk_sdf_trange_(const double *V, const int *T, int f, const mk_sdf_lat_ *L, int i0[3], int i1[3])
{
    for (int a = 0; a < 3; a++) {
        const double v0 = V[3 * (size_t)T[3 * (size_t)f] + a], v1 = V[3 * (size_t)T[3 * (size_t)f + 1] + a], v2 = V[3 * (size_t)T[3 * (size_t)f + 2] + a];
        const double mn = fmin(v0, fmin(v1, v2)), mx = fmax(v0, fmax(v1, v2));
        i0[a] = mk_sdf_ci_(floor((mn - L->lo[a]) / L->h) - 1, 0, L->n[a]), i1[a] = mk_sdf_ci_(ceil((mx - L->lo[a]) / L->h) + 1, -1, L->n[a] - 1);
    }
}
enum { MK_SDF_TILE = 8 }; /* the band's tiles: 8 x 8 rows (j, k), a size fixed by the kit, not by the thread count */
/* the band on the tile of rows j in [j0, j1), k in [k0, k1): its triangles list[0 .. n) in index order, each putting its exact
   distance on its nodes there */
static inline void mk_sdf_band_tile_(const double *V, const int *T, const int *list, long n, const mk_sdf_lat_ *L, int j0, int j1, int k0, int k1,
                                     float *phi, int *near)
{
    for (long e = 0; e < n; e++) {
        const int f = list[e];
        int i0[3], i1[3];
        mk_sdf_trange_(V, T, f, L, i0, i1);
        const int ja = i0[1] > j0 ? i0[1] : j0, jb = i1[1] < j1 - 1 ? i1[1] : j1 - 1, ka = i0[2] > k0 ? i0[2] : k0, kb = i1[2] < k1 - 1 ? i1[2] : k1 - 1;
        for (int k = ka; k <= kb; k++)
            for (int j = ja; j <= jb; j++)
                for (int i = i0[0]; i <= i1[0]; i++) {
                    const double p[3] = {mk_sdf_p_(L, 0, i), mk_sdf_p_(L, 1, j), mk_sdf_p_(L, 2, k)};
                    const double d = mk_sdf_tri_(V, T, f, p);
                    const size_t id = mk_sdf_id_(L, i, j, k);
                    if (d < phi[id]) phi[id] = (float)d, near[id] = f;
                }
    }
}
/* one sweep (di, dj, dk) along row (j, k): each node takes the nearest of its seven sweep-neighbours' triangles. The row reads
   the rows (j - dj, k), (j, k - dk) and (j - dj, k - dk), so rows on one diagonal of the sweep can run at once. */
static inline void mk_sdf_sweep_row_(const double *V, const int *T, const mk_sdf_lat_ *L, float *phi, int *near, int j, int k, int di, int dj, int dk)
{
    const int nx = L->n[0], ia = di > 0 ? 1 : nx - 2, ib = di > 0 ? nx : -1;
    for (int i = ia; i != ib && i >= 0 && i < nx; i += di) {
        const size_t id = mk_sdf_id_(L, i, j, k);
        const double p[3] = {mk_sdf_p_(L, 0, i), mk_sdf_p_(L, 1, j), mk_sdf_p_(L, 2, k)};
        const size_t nb[7] = {mk_sdf_id_(L, i - di, j, k), mk_sdf_id_(L, i, j - dj, k), mk_sdf_id_(L, i - di, j - dj, k), mk_sdf_id_(L, i, j, k - dk),
                              mk_sdf_id_(L, i - di, j, k - dk), mk_sdf_id_(L, i, j - dj, k - dk), mk_sdf_id_(L, i - di, j - dj, k - dk)};
        for (int q = 0; q < 7; q++) {
            const int f = near[nb[q]];
            if (f < 0 || f == near[id]) continue;
            const double d = mk_sdf_tri_(V, T, f, p);
            if (d < phi[id]) phi[id] = (float)d, near[id] = f;
        }
    }
}
enum { MK_SDF_STILE = 4 }; /* the sweeps' tiles: 4 x 4 rows, so a sweep waits on its diagonals of tiles, not of rows */
/* one sweep over the tile of rows r in [r0, r0 + 4), q in [q0, q0 + 4) counted in the sweep's order (row ja + r dj, plane
   ka + q dk; mj, mk rows in all): each row after the rows before it in the tile, the tiles before it being done */
static inline void mk_sdf_sweep_tile_(const double *V, const int *T, const mk_sdf_lat_ *L, float *phi, int *near, int r0, int q0, int mj, int mk, int ja, int ka,
                                      int di, int dj, int dk)
{
    for (int q = q0; q < q0 + MK_SDF_STILE && q < mk; q++)
        for (int r = r0; r < r0 + MK_SDF_STILE && r < mj; r++) mk_sdf_sweep_row_(V, T, L, phi, near, ja + r * dj, ka + q * dk, di, dj, dk);
}
/* MK_SDF_EXACT on nodes [n0, n1): the distance through the tree, the sweeps' value bounding the search */
static inline void mk_sdf_exact_blk_(const mk_bvh *B, const mk_sdf_lat_ *L, float *phi, long n0, long n1)
{
    for (long id = n0; id < n1; id++) {
        const long i = id % L->n[0], j = id / L->n[0] % L->n[1], k = id / ((long)L->n[0] * L->n[1]);
        const double p[3] = {mk_sdf_p_(L, 0, i), mk_sdf_p_(L, 1, j), mk_sdf_p_(L, 2, k)}, b = (double)phi[id] * (1 + 1e-6);
        double q[3];
        int tri;
        const double d2 = mk_bvh_nearest(B, p, b * b, &tri, q);
        if (tri >= 0 && sqrt(d2) < phi[id]) phi[id] = (float)sqrt(d2);
    }
}
/* one row's sign from its crossings cr[0 .. m), sorted along the row: 0 signed by their parity, 1 by a nearby ray (theirs were
   odd), 2 odd on that ray too (the winding number signs it later), -1 no memory. *xs is the thread's own buffer. */
static inline int mk_sdf_row_(const double *V, const int *T, int nt, const mk_sdf_lat_ *L, const mk_sdf_cross_ *cr, long m, double oy, double oz, float *phi,
                              mk_sdf_x_ **xs, size_t *capx)
{
    const int j = (int)(cr[0].row % L->n[1]), k = (int)(cr[0].row / L->n[1]);
#if defined(MK_MUTATE) && MK_MUTATE == 142
    const int odd = 0;
#else
    const int odd = m & 1;
#endif
    if (!odd) {
        long q = 0;
        for (int i = 0; i < L->n[0]; i++) {
            const double x = mk_sdf_p_(L, 0, i);
            while (q < m && cr[q].x < x) q++;
            if (q & 1) phi[mk_sdf_id_(L, i, j, k)] = -phi[mk_sdf_id_(L, i, j, k)];
        }
        return 0;
    }
    long r2 = -1; /* a ray through an edge, or a hole: two more rays nearby */
    for (int r = 1; r <= 2; r++) {
        r2 = mk_sdf_rowx_(V, T, nt, mk_sdf_p_(L, 1, j) + oy * (r == 1 ? 0.61 : 1.53), mk_sdf_p_(L, 2, k) + oz * (r == 1 ? 1.37 : 0.47), xs, capx);
        if (r2 < 0) return -1;
        if (!(r2 & 1)) break;
    }
    if (r2 & 1) return 2;
    long q = 0;
    for (int i = 0; i < L->n[0]; i++) {
        const double x = mk_sdf_p_(L, 0, i);
        while (q < r2 && (*xs)[q].x < x) q++;
        if (q & 1) phi[mk_sdf_id_(L, i, j, k)] = -phi[mk_sdf_id_(L, i, j, k)];
    }
    return 1;
}
/* the odd rows' nodes [n0, n1) of rows[] (row j + ny k each), nx a row: each node's own winding number (through a hole it changes
   smoothly, not by a step at a crossing, so it cannot be carried along the row) */
static inline void mk_sdf_wind_blk_(const double *V, const int *T, int nt, const mk_sdf_lat_ *L, const int64_t *rows, float *phi, long n0, long n1)
{
    for (long q = n0; q < n1; q++) {
        const int i = (int)(q % L->n[0]), j = (int)(rows[q / L->n[0]] % L->n[1]), k = (int)(rows[q / L->n[0]] / L->n[1]);
        const double p[3] = {mk_sdf_p_(L, 0, i), mk_sdf_p_(L, 1, j), mk_sdf_p_(L, 2, k)};
        if (fabs(mk_winding_number(V, T, nt, p)) > 0.5) phi[mk_sdf_id_(L, i, j, k)] = -phi[mk_sdf_id_(L, i, j, k)];
    }
}
/* a node nearer the mesh than the rays' offset: its sign from a ray set off by 1/50 of its distance (the winding number if that
   ray is odd). 0, or -1 for no memory. */
static inline int mk_sdf_near_(const double *V, const int *T, int nt, const mk_sdf_lat_ *L, size_t id, float *phi, mk_sdf_x_ **xs, size_t *capx)
{
    const int i = (int)(id % (size_t)L->n[0]), j = (int)(id / (size_t)L->n[0] % (size_t)L->n[1]), k = (int)(id / ((size_t)L->n[0] * (size_t)L->n[1]));
    const double p[3] = {mk_sdf_p_(L, 0, i), mk_sdf_p_(L, 1, j), mk_sdf_p_(L, 2, k)}, a = fabsf(phi[id]);
    const long m = mk_sdf_rowx_(V, T, nt, p[1] + 0.0123 * a, p[2] + 0.0157 * a, xs, capx);
    if (m < 0) return -1;
    int in;
    if (!(m & 1)) {
        long q = 0;
        while (q < m && (*xs)[q].x < p[0]) q++;
        in = q & 1;
    } else in = fabs(mk_winding_number(V, T, nt, p)) > 0.5;
    phi[id] = in ? -(float)a : (float)a;
    return 0;
}
static inline void mk_sdf_init_blk_(float *phi, int *near, float big, long n0, long n1) { for (long i = n0; i < n1; i++) phi[i] = big, near[i] = -1; }

/* the mesh checked (nv, nt >= 1, indices in range, coordinates finite) and its box into in[6]: 0, or -1 */
static inline int mk_sdf_check_(const double *V, int nv, const int *T, int nt, double in[6])
{
    if (nv < 1 || nt < 1) return -1;
    for (size_t q = 0; q < 3 * (size_t)nv; q++) if (!isfinite(V[q])) return -1;
    for (int a = 0; a < 3; a++) in[a] = INFINITY, in[3 + a] = -INFINITY;
    for (size_t q = 0; q < 3 * (size_t)nt; q++) {
        if (T[q] < 0 || T[q] >= nv) return -1;
        for (int a = 0; a < 3; a++) in[a] = fmin(in[a], V[3 * (size_t)T[q] + a]), in[3 + a] = fmax(in[3 + a], V[3 * (size_t)T[q] + a]);
    }
    return 0;
}
static inline int mk_sdf_grid_(const double *V, const int *T, int nt, int flags, size_t N, double extra, mk_sdf *g);
/* fill g for the mesh with cell h and margin beyond the mesh's box on every side. 0, or -1 for a bad mesh (nv or nt < 1, an
   index out of range, a coordinate not finite), h not finite and > 0, margin not finite and >= 0, a grid whose node count
   overflows (refused before anything is allocated), or no memory; *g is then empty. */
static inline int mk_sdf_mesh(const double *V, int nv, const int *T, int nt, double h, double margin, int flags, mk_sdf *g)
{
    memset(g, 0, sizeof *g);
    double in[6];
    if (!(h > 0) || !isfinite(h) || !(margin >= 0) || !isfinite(margin) || mk_sdf_check_(V, nv, T, nt, in)) return -1;
    const double *lo = in, *hi = in + 3;
    size_t N = 1;
    for (int a = 0; a < 3; a++) {
        /* MK_SDF_LATTICE (metal/ref/hair hairgen.c): node 0 on the world lattice h Z, so a cut or shifted mesh keeps its nodes */
        const double o = flags & MK_SDF_LATTICE ? floor((lo[a] - margin) / h) * h : lo[a] - margin;
        const double m = flags & MK_SDF_LATTICE ? ceil((hi[a] + margin - o) / h) + 1 : ceil((hi[a] - lo[a] + 2 * margin) / h) + 1;
        if (!(m >= 1 && m <= INT_MAX - 1)) { memset(g, 0, sizeof *g); return -1; }
        g->n[a] = (int)m, g->lo[a] = o, g->in[a] = lo[a], g->in[3 + a] = hi[a];
        if ((size_t)g->n[a] > SIZE_MAX / sizeof(double) / N) { memset(g, 0, sizeof *g); return -1; }
        N *= (size_t)g->n[a];
    }
    g->h = h;
    return mk_sdf_grid_(V, T, nt, flags, N, 0, g);
}
/* 1.5.2: the same on the nodes of a box (lo x y z, hi x y z) instead of the mesh's whole box: a local field on a large mesh (the maintainer's
   race-c1 H3 grids a few cm round one contact of a muscle whose own box is 400 M nodes at 0.5 mm; mk_graft's patch needs only the
   seam's neighbourhood). Node 0 at lo (MK_SDF_LATTICE: on the world lattice at or below lo), nodes to hi or just past it. Every
   triangle still counts: a row's sign from all its crossings, inside the box or not. A node no triangle comes within two cells
   of the box's nodes to reach takes the distance to the nearest one that does (an upper bound; MK_SDF_EXACT makes every node
   exact). MK_SDF_FILL is refused (its flood starts at the grid's faces, which a box may cut inside the mesh). -1 as mk_sdf_mesh,
   or for a box not finite or with hi below lo. */
static inline int mk_sdf_mesh_box(const double *V, int nv, const int *T, int nt, double h, const double box[6], int flags, mk_sdf *g)
{
    memset(g, 0, sizeof *g);
    double in[6];
    if (!(h > 0) || !isfinite(h) || !box || (flags & MK_SDF_FILL) || mk_sdf_check_(V, nv, T, nt, in)) return -1;
    size_t N = 1;
    double extra = 0;
    for (int a = 0; a < 3; a++) {
#if defined(MK_MUTATE) && MK_MUTATE == 236
        if (!isfinite(box[a]) || !isfinite(box[3 + a]) || !(box[3 + a] >= box[a])) return -1; /* 236: axes before it left set (1.5.3) */
#else
        if (!isfinite(box[a]) || !isfinite(box[3 + a]) || !(box[3 + a] >= box[a])) { memset(g, 0, sizeof *g); return -1; } /* the verification work S-1: axes before it were set */
#endif
        const double o = flags & MK_SDF_LATTICE ? floor(box[a] / h) * h : box[a], m = ceil((box[3 + a] - o) / h) + 1;
        if (!(m >= 1 && m <= INT_MAX - 1)) { memset(g, 0, sizeof *g); return -1; }
        g->n[a] = (int)m, g->lo[a] = o, g->in[a] = in[a], g->in[3 + a] = in[3 + a];
        if ((size_t)g->n[a] > SIZE_MAX / sizeof(double) / N) { memset(g, 0, sizeof *g); return -1; }
        N *= (size_t)g->n[a];
#if !defined(MK_MUTATE) || MK_MUTATE != 207
        extra += fmax(in[3 + a], o + (m - 1) * h) - fmin(in[a], o); /* the box and the mesh's together: no node is farther from it */
#endif
    }
    g->h = h;
    return mk_sdf_grid_(V, T, nt, flags, N, extra, g);
}
/* the grid's work once its nodes are set (g->n, g->lo, g->h, g->in): band, sweeps, MK_SDF_EXACT, sign, MK_SDF_FILL. extra widens
   BIG, a node's value until a triangle reaches it, past the farthest a node lies from the mesh (0: the grid holds the mesh's box). */
static inline int mk_sdf_grid_(const double *V, const int *T, int nt, int flags, size_t N, double extra, mk_sdf *g)
{
    const double h = g->h;
    double bi[3]; /* MK_SDF_LATTICE: node 0's lattice index, so a node is (bi + i) h exactly, the same in every grid */
    const int lat = flags & MK_SDF_LATTICE;
    for (int a = 0; a < 3; a++) bi[a] = lat ? round(g->lo[a] / h) : 0;
#define MK_SDF_P_(a, i) (lat ? (bi[a] + (double)(i)) * h : g->lo[a] + (i) * h)
    const int nx = g->n[0], ny = g->n[1], nz = g->n[2];
    const mk_sdf_lat_ L = {{nx, ny, nz}, lat, {g->lo[0], g->lo[1], g->lo[2]}, {bi[0], bi[1], bi[2]}, h};
    float *phi = malloc(sizeof *phi * N);
    int *near = malloc(sizeof *near * N);
    if (!phi || !near) { free(phi), free(near), memset(g, 0, sizeof *g); return -1; }
    const float BIG = (float)(((double)nx + ny + nz) * h + extra);
    MK_PAR_FOR_((long)N, n0, n1, mk_sdf_init_blk_(phi, near, BIG, n0, n1));
#define MK_SDF_ID_(i, j, k) ((size_t)(i) + (size_t)nx * ((size_t)(j) + (size_t)ny * (size_t)(k)))
    /* the band: exact distances within two cells of each triangle, by tiles of rows; each tile takes the triangles that reach it
       in index order, so a node sees its triangles in the serial order whatever the threads */
    {
#if defined(MK_MUTATE) && MK_MUTATE == 147
        const int nthr = MK_PAR_NTHREADS_(), sj = ny, sk = (nz + nthr - 1) / nthr; /* 147: tiles cut by the thread count, a triangle listed in its first tile only */
#else
        const int sj = MK_SDF_TILE, sk = MK_SDF_TILE;
#endif
        const int tj = (ny + sj - 1) / sj, tk = (nz + sk - 1) / sk;
        const long ntile = (long)tj * tk;
        long *start = calloc((size_t)ntile + 1, sizeof *start);
        if (!start) { free(phi), free(near), memset(g, 0, sizeof *g); return -1; }
        for (int pass = 0; pass < 2; pass++) { /* count, then list (start[t] runs through tile t's places) */
            if (pass) {
                for (long t = 0; t < ntile; t++) start[t + 1] += start[t];
            }
            int *list = pass ? malloc(sizeof *list * (size_t)(start[ntile] ? start[ntile] : 1)) : NULL;
            if (pass && !list) { free(start), free(phi), free(near), memset(g, 0, sizeof *g); return -1; }
            for (int f = 0; f < nt; f++) {
                int i0[3], i1[3];
                mk_sdf_trange_(V, T, f, &L, i0, i1);
#if defined(MK_MUTATE) && MK_MUTATE == 208
                if (i0[0] >= i1[0] || i0[1] >= i1[1] || i0[2] >= i1[2]) continue; /* a range one node thick taken for none */
#else
                if (i0[0] > i1[0] || i0[1] > i1[1] || i0[2] > i1[2]) continue; /* off the grid (a box): it reaches no node */
#endif
#if defined(MK_MUTATE) && MK_MUTATE == 147
                i1[2] = i0[2];
#endif
                for (int b = i0[2] / sk; b <= i1[2] / sk; b++)
                    for (int a = i0[1] / sj; a <= i1[1] / sj; a++)
                        if (pass) list[start[(long)b * tj + a]++] = f;
                        else start[(long)b * tj + a + 1]++;
            }
            if (!pass) continue;
            for (long t = ntile; t > 0; t--) start[t] = start[t - 1];
            start[0] = 0;
            MK_OMP_(omp parallel for schedule(dynamic, 1) if(ntile > 1))
            for (long t = 0; t < ntile; t++) {
                const int a = (int)(t % tj), b = (int)(t / tj);
                mk_sdf_band_tile_(V, T, list + start[t], start[t + 1] - start[t], &L, a * sj, (a + 1) * sj < ny ? (a + 1) * sj : ny, b * sk,
                                  (b + 1) * sk < nz ? (b + 1) * sk : nz, phi, near);
            }
            free(list);
        }
        free(start);
    }
    /* the sweeps, each as a wavefront: rows (j, k) on one diagonal of the sweep need only the diagonals before it */
#if !(defined(MK_MUTATE) && MK_MUTATE == 141)
    if (ny > 1 && nz > 1) {
        const int mj = ny - 1, mk = nz - 1; /* rows a sweep visits along j and along k */
#if defined(MK_MUTATE) && MK_MUTATE == 146
        const int nthr = MK_PAR_NTHREADS_(), cs = (mk + nthr - 1) / nthr; /* 146: planes in chunks by the thread count, the last chunk first */
#endif
        MK_OMP_(omp parallel)
        for (int pass = 0; pass < 2; pass++)
            for (int sw = 0; sw < 8; sw++) {
                const int di = sw & 1 ? -1 : 1, dj = sw & 2 ? -1 : 1, dk = sw & 4 ? -1 : 1, ja = dj > 0 ? 1 : ny - 2, ka = dk > 0 ? 1 : nz - 2;
#if defined(MK_MUTATE) && MK_MUTATE == 146
                MK_OMP_(omp single)
                for (int c = nthr - 1; c >= 0; c--)
                    for (int q = c * cs; q < (c + 1) * cs && q < mk; q++)
                        for (int r = 0; r < mj; r++) mk_sdf_sweep_row_(V, T, &L, phi, near, ja + r * dj, ka + q * dk, di, dj, dk);
#else
                const int tj = (mj + MK_SDF_STILE - 1) / MK_SDF_STILE, tk = (mk + MK_SDF_STILE - 1) / MK_SDF_STILE;
                for (int d = 0; d <= tj + tk - 2; d++) { /* the tiles (a, c) with a + c = d */
                    const int c0 = d - (tj - 1) > 0 ? d - (tj - 1) : 0, c1 = d < tk - 1 ? d : tk - 1;
                    MK_OMP_(omp for schedule(dynamic, 1))
                    for (int c = c0; c <= c1; c++) mk_sdf_sweep_tile_(V, T, &L, phi, near, (d - c) * MK_SDF_STILE, c * MK_SDF_STILE, mj, mk, ja, ka, di, dj, dk);
                }
#endif
            }
    }
#endif
    free(near);
    /* MK_SDF_EXACT: every node's distance made exact through mk_geom's tree, the sweeps' value bounding the search */
    if (flags & MK_SDF_EXACT) {
        mk_bvh B;
        if (mk_bvh_build(&B, V, T, nt)) { free(phi), memset(g, 0, sizeof *g); return -1; }
#if defined(MK_MUTATE) && MK_MUTATE == 148
        MK_PAR_CUT_((long)N, n0, n1, mk_sdf_exact_blk_(&B, &L, phi, n0, n1)); /* 148: nodes cut by the thread count, each piece's last node left out */
#else
        const long nb = mk_par_nblk((long)N);
        MK_OMP_(omp parallel for schedule(dynamic, 1) if(nb > 1))
        for (long b = 0; b < nb; b++) {
            long n0, n1;
            mk_par_range((long)N, b, &n0, &n1), mk_sdf_exact_blk_(&B, &L, phi, n0, n1);
        }
#endif
        mk_bvh_free(&B);
    }
    /* the sign: crossings of a ray along the first axis through each row (j, k), set off from the row's nodes by fractions of
       h; each triangle adds its crossing to the rows its (y, z) box covers */
    size_t ncr = 0, capc = 1 << 16;
    mk_sdf_cross_ *cr = malloc(sizeof *cr * capc);
    const double oy = 1.234567e-4 * h, oz = 2.345678e-4 * h;
    for (int f = 0; f < nt && cr; f++) {
        const double *A = V + 3 * (size_t)T[3 * (size_t)f], *B = V + 3 * (size_t)T[3 * (size_t)f + 1], *C = V + 3 * (size_t)T[3 * (size_t)f + 2];
        const double ymn = fmin(A[1], fmin(B[1], C[1])), ymx = fmax(A[1], fmax(B[1], C[1])), zmn = fmin(A[2], fmin(B[2], C[2])), zmx = fmax(A[2], fmax(B[2], C[2]));
        /* clamped before the conversion (the verification work S-2: a triangle more than INT_MAX cells from a box's node 0) */
        const int j0 = mk_sdf_ci_(ceil((ymn - oy - g->lo[1]) / h), 0, ny), j1 = mk_sdf_ci_(floor((ymx - oy - g->lo[1]) / h), -1, ny - 1);
        const int k0 = mk_sdf_ci_(ceil((zmn - oz - g->lo[2]) / h), 0, nz), k1 = mk_sdf_ci_(floor((zmx - oz - g->lo[2]) / h), -1, nz - 1);
        for (int k = k0; k <= k1 && cr; k++)
            for (int j = j0; j <= j1; j++) {
                const double y = MK_SDF_P_(1, j) + oy, z = MK_SDF_P_(2, k) + oz;
                const double ay = A[1] - y, az = A[2] - z, by = B[1] - y, bz = B[2] - z, cy = C[1] - y, cz = C[2] - z;
                const double w0 = by * cz - bz * cy, w1 = cy * az - cz * ay, w2 = ay * bz - az * by;
                if (!((w0 >= 0 && w1 >= 0 && w2 >= 0) || (w0 <= 0 && w1 <= 0 && w2 <= 0))) continue;
                const double s = w0 + w1 + w2;
                if (s == 0) continue; /* the triangle seen edge-on */
                if (ncr == capc) {
                    mk_sdf_cross_ *c2 = realloc(cr, sizeof *cr * (capc *= 2));
                    if (!c2) { free(cr), cr = NULL; break; }
                    cr = c2;
                }
                cr[ncr].row = j + (int64_t)ny * k, cr[ncr].x = (w0 * A[0] + w1 * B[0] + w2 * C[0]) / s, ncr++;
            }
    }
    if (!cr) { free(phi), memset(g, 0, sizeof *g); return -1; }
    qsort(cr, ncr, sizeof *cr, mk_sdf_cmp_);
    /* each row with crossings, cr[seg[s] .. seg[s + 1]), signed on its own; the rows odd on every ray then node by node */
    long nseg = 0;
    for (size_t c = 0; c < ncr; c++) nseg += c == 0 || cr[c].row != cr[c - 1].row;
    size_t *seg = malloc(sizeof *seg * ((size_t)nseg + 1));
    signed char *st = malloc((size_t)nseg + 1);
    if (!seg || !st) { free(seg), free(st), free(cr), free(phi), memset(g, 0, sizeof *g); return -1; }
    for (size_t c = 0, s = 0; c < ncr; c++)
        if (c == 0 || cr[c].row != cr[c - 1].row) seg[s++] = c;
    seg[nseg] = ncr;
    MK_OMP_(omp parallel if(nseg > 1))
    {
        mk_sdf_x_ *xs = NULL; /* each thread's own crossings buffer */
        size_t capx = 0;
        MK_OMP_(omp for schedule(dynamic, 16))
        for (long s = 0; s < nseg; s++) st[s] = (signed char)mk_sdf_row_(V, T, nt, &L, cr + seg[s], (long)(seg[s + 1] - seg[s]), oy, oz, phi, &xs, &capx);
        free(xs);
    }
    int bad = 0;
    for (long s = 0; s < nseg; s++) bad |= st[s] < 0, g->reray += st[s] == 1, g->odd_rows += st[s] == 2;
    int64_t *rows = bad ? NULL : malloc(sizeof *rows * (size_t)(g->odd_rows ? g->odd_rows : 1));
    for (long s = 0, w = 0; rows && s < nseg; s++)
        if (st[s] == 2) rows[w++] = cr[seg[s]].row;
    free(seg), free(st), free(cr);
    if (!rows) { free(phi), memset(g, 0, sizeof *g); return -1; }
    const long nw = g->odd_rows * (long)nx; /* a node a block: each is O(triangles) */
    MK_OMP_(omp parallel for schedule(dynamic, 1) if(nw > 1))
    for (long q = 0; q < nw; q++) mk_sdf_wind_blk_(V, T, nt, &L, rows, phi, q, q + 1);
    free(rows);
    /* a node nearer the mesh than the rays' offset (4e-4 h; the offset is 2.7e-4 h) may lie on the other side of its ray: a
       ray set off by 1/50 of its distance decides it, every point within that distance being on its side; the winding number
       if that ray too is odd. Each such node on its own. */
    long ncand = 0;
    for (size_t id = 0; id < N; id++) ncand += fabsf(phi[id]) < 4e-4 * h && phi[id] != 0 && !(MK_SDF_MUT145_);
    if (ncand) {
        size_t *cand = malloc(sizeof *cand * (size_t)ncand);
        if (!cand) { free(phi), memset(g, 0, sizeof *g); return -1; }
        for (size_t id = 0, c = 0; id < N; id++)
            if (fabsf(phi[id]) < 4e-4 * h && phi[id] != 0) cand[c++] = id;
        MK_OMP_(omp parallel if(ncand > 1))
        {
            mk_sdf_x_ *xs = NULL;
            size_t capx = 0;
            MK_OMP_(omp for schedule(dynamic, 8))
            for (long c = 0; c < ncand; c++)
                if (mk_sdf_near_(V, T, nt, &L, cand[c], phi, &xs, &capx)) {
                    MK_OMP_(omp atomic write)
                    bad = 1;
                }
            free(xs);
        }
        free(cand);
        if (bad) { free(phi), memset(g, 0, sizeof *g); return -1; }
        g->near = ncand;
    }
    /* MK_SDF_FILL: outside nodes not joined to the grid's faces through outside nodes are inside */
    if (flags & MK_SDF_FILL) {
        unsigned char *seen = calloc(N, 1);
        size_t *qu = malloc(sizeof *qu * N), qh = 0, qt = 0;
        (void)qh;
        if (!seen || !qu) { free(seen), free(qu), free(phi), memset(g, 0, sizeof *g); return -1; }
        for (int k = 0; k < nz; k++)
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++)
                    if (i == 0 || j == 0 || k == 0 || i == nx - 1 || j == ny - 1 || k == nz - 1) {
                        const size_t id = MK_SDF_ID_(i, j, k);
                        if (phi[id] > 0 && !seen[id]) seen[id] = 1, qu[qt++] = id;
                    }
#if !(defined(MK_MUTATE) && MK_MUTATE == 144)
        while (qh < qt) {
            const size_t id = qu[qh++];
            const int i = (int)(id % (size_t)nx), j = (int)(id / (size_t)nx % (size_t)ny), k = (int)(id / ((size_t)nx * (size_t)ny));
            const int D[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
            for (int e = 0; e < 6; e++) {
                const int a = i + D[e][0], b = j + D[e][1], c = k + D[e][2];
                if (a < 0 || b < 0 || c < 0 || a >= nx || b >= ny || c >= nz) continue;
                const size_t nb = MK_SDF_ID_(a, b, c);
                if (!seen[nb] && phi[nb] > 0) seen[nb] = 1, qu[qt++] = nb;
            }
        }
#endif
        for (size_t id = 0; id < N; id++) if (phi[id] > 0 && !seen[id]) phi[id] = -phi[id], g->filled++;
        free(seen), free(qu);
    }
#undef MK_SDF_ID_
#undef MK_SDF_P_
    g->d = phi;
    return 0;
}
/* the distance at p by trilinear interpolation, clamped to the grid's box; NaN at a NaN coordinate (1.6.0: 1.5.6 read node 0) */
static inline double mk_sdf_sample(const mk_sdf *g, const double *p)
{
    if (!g->d || isnan(p[0]) || isnan(p[1]) || isnan(p[2])) return NAN;
    int i0[3];
    double t[3];
    for (int a = 0; a < 3; a++) {
        double f = (p[a] - g->lo[a]) / g->h;
        if (!(f > 0)) f = 0;
        if (f > g->n[a] - 1) f = g->n[a] - 1;
        int i = (int)floor(f);
        if (i > g->n[a] - 2) i = g->n[a] - 2;
        if (i < 0) i = 0;
        i0[a] = i, t[a] = g->n[a] > 1 ? f - i : 0;
    }
    double v = 0;
    for (int q = 0; q < 8; q++) {
        int c[3];
        double w = 1;
        for (int a = 0; a < 3; a++) {
            const int s = (q >> a) & 1;
            c[a] = i0[a] + s < g->n[a] ? i0[a] + s : i0[a];
            w *= s ? t[a] : 1 - t[a];
        }
        v += w * g->d[(size_t)c[0] + (size_t)g->n[0] * ((size_t)c[1] + (size_t)g->n[1] * (size_t)c[2])];
    }
    return v;
}
#endif
