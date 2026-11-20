/* mathkit: geometry on triangle meshes, double precision. Closest point on a triangle (Ericson's region method), watertight
   ray-triangle intersection (Woop, Benthin and Wald 2013: no ray slips between two triangles sharing an edge while the ray's sheared frame resolves
   them, the verification work F11), boxes, and a
   triangle bounding-volume tree (binned SAH build) with first-hit, any-hit, hit-count, nearest-point and point-in-mesh
   queries; the generalized winding number for a cross-check. Meshes are flat arrays: V[3 v + k], T[3 t + j]. Ties go to the
   lowest triangle index, so results do not depend on the tree. Crossings (1.5.2): an edge strictly through a triangle
   (meshcheck's test, mk_seg_tri_strict), two faces crossing (mk_tri_cross), and a mesh's crossing face pairs through the tree,
   all or those of marked faces (mk_bvh_crossings). Batches of points (mk_bvh_nearest_many, mk_bvh_inside_many,
   mk_winding_number_many) run on OpenMP threads when the caller builds with -fopenmp, each answer the single call's. Header
   only, C23. */
#ifndef MK_GEOM_H
#define MK_GEOM_H
#include "mk_par.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* the closest point's edge candidates for mk_closest_pt_tri: on a sliver, a repeated corner or a weight below 0 or NaN, each
   edge's own closest point replaces q when nearer (exact arithmetic would never prefer them; rounding there may). A candidate
   with a negative or NaN weight is dropped. Kept out of line: it is rare, and inlined it took mk_closest_pt_tri past gcc's inline
   limit, which cost mk_sdf_mesh 30 to 47 % at -O3 (1.3.0, MEASURED). */
#if defined(__GNUC__)
#define MK_COLD_ __attribute__((noinline, cold)) /* GNU syntax, not a C23 [[...]] attribute (the verification work N-B) */
#else
#define MK_COLD_
#endif
MK_COLD_ static double mk_closest_edges_(const double *p, const double *a, const double *b, const double *c, double d, double *w, double *q)
{
    int bad = !(w[0] >= 0 && w[1] >= 0 && w[2] >= 0);
    const double *P[3] = {a, b, c};
    for (int e = 0; e < 3; e++) {
        const double *u = P[e], *v = P[(e + 1) % 3];
        double l = 0, s = 0;
        for (int k = 0; k < 3; k++) l += (v[k] - u[k]) * (v[k] - u[k]), s += (v[k] - u[k]) * (p[k] - u[k]);
        double t = l > 0 ? fmin(1, fmax(0, s / l)) : 0, we[3] = {0, 0, 0}, qe[3], de = 0;
        we[e] = 1 - t, we[(e + 1) % 3] = t;
        for (int k = 0; k < 3; k++) qe[k] = we[0] * a[k] + we[1] * b[k] + we[2] * c[k], de += (qe[k] - p[k]) * (qe[k] - p[k]);
        if (bad || de < d) {
            bad = 0, d = de;
            for (int k = 0; k < 3; k++) w[k] = we[k], q[k] = qe[k];
        }
    }
    return d;
}
/* the closest point q on triangle abc to p, its barycentric weights w (q = w0 a + w1 b + w2 c); the squared distance */
static inline double mk_closest_pt_tri(const double *p, const double *a, const double *b, const double *c, double *w, double *q)
{
    double ab[3], ac[3], ap[3], bp[3], cp[3];
    for (int k = 0; k < 3; k++) ab[k] = b[k] - a[k], ac[k] = c[k] - a[k], ap[k] = p[k] - a[k], bp[k] = p[k] - b[k], cp[k] = p[k] - c[k];
#define MK_DOT3_(u, v) (u[0] * v[0] + u[1] * v[1] + u[2] * v[2])
    double d1 = MK_DOT3_(ab, ap), d2 = MK_DOT3_(ac, ap), d3 = MK_DOT3_(ab, bp), d4 = MK_DOT3_(ac, bp), d5 = MK_DOT3_(ab, cp), d6 = MK_DOT3_(ac, cp);
#undef MK_DOT3_
    double vc = d1 * d4 - d3 * d2, vb = d5 * d2 - d1 * d6, va = d3 * d6 - d5 * d4;
    int face = 0;
    if (d1 <= 0 && d2 <= 0) w[0] = 1, w[1] = w[2] = 0;
#if defined(MK_MUTATE) && MK_MUTATE == 91
    else if (d3 >= 0 && d4 >= d3) w[1] = 1, w[0] = w[2] = 0;
#else
    else if (d3 >= 0 && d4 <= d3) w[1] = 1, w[0] = w[2] = 0;
#endif
    else if (vc <= 0 && d1 >= 0 && d3 <= 0) { double v = d1 / (d1 - d3); w[0] = 1 - v, w[1] = v, w[2] = 0, face = 1; }
    else if (d6 >= 0 && d5 <= d6) w[2] = 1, w[0] = w[1] = 0;
    else if (vb <= 0 && d2 >= 0 && d6 <= 0) { double v = d2 / (d2 - d6); w[0] = 1 - v, w[2] = v, w[1] = 0, face = 1; }
    else if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) { double v = (d4 - d3) / ((d4 - d3) + (d5 - d6)); w[1] = 1 - v, w[2] = v, w[0] = 0, face = 1; }
    else { double den = 1 / (va + vb + vc), v = vb * den, ww = vc * den; w[0] = 1 - v - ww, w[1] = v, w[2] = ww, face = 1; }
    double d = 0;
    for (int k = 0; k < 3; k++) q[k] = w[0] * a[k] + w[1] * b[k] + w[2] * c[k], d += (q[k] - p[k]) * (q[k] - p[k]);
    if (face) { /* an edge or face region: the edges too, only on a sliver (va + vb + vc is |ab x ac|^2, here under
                   2^-20 |ab|^2 |ac|^2), a repeated corner or a weight below 0 or NaN (1.3.0: elsewhere the edges never win
                   but by an ulp) */
        const double lab = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2], lac = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
        if (!(w[0] >= 0 && w[1] >= 0 && w[2] >= 0) || !(va + vb + vc > 0x1p-20 * lab * lac)) d = mk_closest_edges_(p, a, b, c, d, w, q);
    }
    return d;
}

/* a ray prepared for the watertight test: the axis of largest direction becomes z, then a shear. A zero d gives NaN shears, and
   then no query hits anything (MEASURED for first_hit, any_hit and count_hits, 1.6.0). */
typedef struct { double o[3], d[3], inv[3], Sx, Sy, Sz; int kx, ky, kz; } mk_ray;
static inline mk_ray mk_ray_make(const double *o, const double *d)
{
    mk_ray r;
    for (int k = 0; k < 3; k++) r.o[k] = o[k], r.d[k] = d[k], r.inv[k] = 1 / d[k];
    int kz = fabs(d[0]) > fabs(d[1]) ? (fabs(d[0]) > fabs(d[2]) ? 0 : 2) : (fabs(d[1]) > fabs(d[2]) ? 1 : 2);
    int kx = (kz + 1) % 3, ky = (kx + 1) % 3;
    if (d[kz] < 0) { int t = kx; kx = ky, ky = t; }
    r.kx = kx, r.ky = ky, r.kz = kz, r.Sx = d[kx] / d[kz], r.Sy = d[ky] / d[kz], r.Sz = 1 / d[kz];
    return r;
}
/* a * b - c * d with relative error at most 2u (Kahan): its sign and its zero are exact */
static inline double mk_dop_(double a, double b, double c, double d)
{
    double w = c * d, e = fma(-c, d, w), f = fma(a, b, -w);
    return f + e;
}
/* mk_ray_tri, and with out its facing (1.5.3): 1 where the ray leaves through the face (runs along (b - a) x (c - a)), -1 where
   it enters. The sheared frame keeps orientation, and there det = U + V + W is minus that normal's component along the ray. */
static inline int mk_ray_tri_s_(const mk_ray *r, const double *a, const double *b, const double *c, double tmin, double tmax, double *t, double *bary, int *out)
{
    double A[3], B[3], C[3];
    for (int k = 0; k < 3; k++) A[k] = a[k] - r->o[k], B[k] = b[k] - r->o[k], C[k] = c[k] - r->o[k];
    double Ax = A[r->kx] - r->Sx * A[r->kz], Ay = A[r->ky] - r->Sy * A[r->kz];
    double Bx = B[r->kx] - r->Sx * B[r->kz], By = B[r->ky] - r->Sy * B[r->kz];
    double Cx = C[r->kx] - r->Sx * C[r->kz], Cy = C[r->ky] - r->Sy * C[r->kz];
    double U = mk_dop_(Cx, By, Cy, Bx), V = mk_dop_(Ax, Cy, Ay, Cx), W = mk_dop_(Bx, Ay, By, Ax);
#if defined(MK_MUTATE) && MK_MUTATE == 92
    if (U < 0 || V < 0 || W < 0) return 0;
#else
    if ((U < 0 || V < 0 || W < 0) && (U > 0 || V > 0 || W > 0)) return 0;
#endif
    double det = U + V + W;
    if (det == 0) return 0;
    double T = U * r->Sz * A[r->kz] + V * r->Sz * B[r->kz] + W * r->Sz * C[r->kz], tt = T / det;
    if (!(tt > tmin && tt < tmax)) return 0;
    *t = tt;
    if (bary) bary[0] = U / det, bary[1] = V / det, bary[2] = W / det;
#if defined(MK_MUTATE) && MK_MUTATE == 214
    if (out) *out = det > 0 ? 1 : -1; /* 214: the facing turned */
#else
    if (out) *out = det < 0 ? 1 : -1;
#endif
    return 1;
}
/* hit of r with triangle abc at a parameter in (tmin, tmax): 1 and *t, bary[3] (weights of a, b, c); both faces count */
static inline int mk_ray_tri(const mk_ray *r, const double *a, const double *b, const double *c, double tmin, double tmax, double *t, double *bary)
{
    return mk_ray_tri_s_(r, a, b, c, tmin, tmax, t, bary, NULL);
}

typedef struct { double lo[3], hi[3]; } mk_aabb;
static inline void mk_aabb_empty(mk_aabb *b) { for (int k = 0; k < 3; k++) b->lo[k] = INFINITY, b->hi[k] = -INFINITY; }
static inline void mk_aabb_add(mk_aabb *b, const double *p) { for (int k = 0; k < 3; k++) b->lo[k] = fmin(b->lo[k], p[k]), b->hi[k] = fmax(b->hi[k], p[k]); }
static inline void mk_aabb_join(mk_aabb *b, const mk_aabb *c) { for (int k = 0; k < 3; k++) b->lo[k] = fmin(b->lo[k], c->lo[k]), b->hi[k] = fmax(b->hi[k], c->hi[k]); }
static inline double mk_aabb_area(const mk_aabb *b)
{
    double x = b->hi[0] - b->lo[0], y = b->hi[1] - b->lo[1], z = b->hi[2] - b->lo[2];
    return x < 0 ? 0 : 2 * (x * y + y * z + z * x);
}
/* squared distance from p to the box (0 inside) */
static inline double mk_aabb_dist2(const mk_aabb *b, const double *p)
{
    double d = 0;
    for (int k = 0; k < 3; k++) {
        double e = p[k] < b->lo[k] ? b->lo[k] - p[k] : p[k] > b->hi[k] ? p[k] - b->hi[k] : 0;
        d += e * e;
    }
    return d;
}
/* the slab test; *tenter the entry parameter. The box is padded by a few ulps of its and the ray origin's coordinates, so a hit on its face is not lost. */
static inline int mk_aabb_ray(const mk_aabb *b, const mk_ray *r, double tmin, double tmax, double *tenter)
{
    double t0 = tmin, t1 = tmax, m = 1;
#if defined(MK_MUTATE) && MK_MUTATE == 235
    for (int k = 0; k < 3; k++) m = fmax(m, fmax(fabs(b->lo[k]), fabs(b->hi[k]))); /* 235: the box's own coordinates only (1.5.3) */
#else
    for (int k = 0; k < 3; k++) m = fmax(m, fmax(fmax(fabs(b->lo[k]), fabs(b->hi[k])), fabs(r->o[k]))); /* the ray origin too: (lo - o) rounds on its scale (the verification work H19 G1) */
#endif
    for (int k = 0; k < 3; k++) {
        double pad = 1e-12 * m;
        double a = (b->lo[k] - pad - r->o[k]) * r->inv[k], c = (b->hi[k] + pad - r->o[k]) * r->inv[k];
        t0 = fmax(t0, fmin(a, c)), t1 = fmin(t1, fmax(a, c));
    }
    *tenter = t0;
    return t0 <= t1;
}

typedef struct { mk_aabb box; int first, count; } mk_bvh_node; /* count > 0: a leaf of tri[first .. first+count); else children first, first+1 */
typedef struct { const double *V; const int *T; int nt, nnode; int *tri; mk_bvh_node *node; } mk_bvh;

static inline void mk_tri_box_(const double *V, const int *T, int t, mk_aabb *b)
{
    mk_aabb_empty(b);
    for (int j = 0; j < 3; j++) mk_aabb_add(b, V + 3 * (size_t)T[3 * (size_t)t + j]);
}
#define MK_BVH_SAH_DEPTH 200
static inline int mk_bvh_split_(mk_bvh *B, const double *cen, int node, int first, int count, int depth)
{
    enum { NB = 16, LEAF = 4 };
    mk_bvh_node *nd = &B->node[node];
    mk_aabb_empty(&nd->box);
    mk_aabb cb;
    mk_aabb_empty(&cb);
    for (int i = first; i < first + count; i++) {
        mk_aabb tb;
        mk_tri_box_(B->V, B->T, B->tri[i], &tb), mk_aabb_join(&nd->box, &tb), mk_aabb_add(&cb, cen + 3 * (size_t)B->tri[i]);
    }
    nd->first = first, nd->count = count;
    if (count <= LEAF) return 0;
    int ax = 0;
    for (int k = 1; k < 3; k++) if (cb.hi[k] - cb.lo[k] > cb.hi[ax] - cb.lo[ax]) ax = k;
    double ext = cb.hi[ax] - cb.lo[ax];
    int mid = first + count / 2;
    if (depth >= MK_BVH_SAH_DEPTH) ext = 0;
    if (ext > 0) {
        mk_aabb bb[NB];
        int bn[NB] = {0};
        for (int i = 0; i < NB; i++) mk_aabb_empty(&bb[i]);
        for (int i = first; i < first + count; i++) {
            int t = B->tri[i];
            double f = NB * (cen[3 * (size_t)t + ax] - cb.lo[ax]) / ext;
            int k = f >= NB ? NB - 1 : f > 0 ? (int)f : 0;
            mk_aabb tb;
            mk_tri_box_(B->V, B->T, t, &tb), mk_aabb_join(&bb[k], &tb), bn[k]++;
        }
        double best = INFINITY;
        int bs = -1;
        for (int s = 1; s < NB; s++) {
            mk_aabb l, r;
            int nl = 0, nr = 0;
            mk_aabb_empty(&l), mk_aabb_empty(&r);
            for (int i = 0; i < s; i++) mk_aabb_join(&l, &bb[i]), nl += bn[i];
            for (int i = s; i < NB; i++) mk_aabb_join(&r, &bb[i]), nr += bn[i];
            if (!nl || !nr) continue;
            double c = nl * mk_aabb_area(&l) + nr * mk_aabb_area(&r);
            if (c < best) best = c, bs = s;
        }
        if (bs < 0) ext = 0;
        else {
            if (best >= count * mk_aabb_area(&nd->box) && count <= 16) return 0; /* splitting costs more than testing all */
            int i = first, j = first + count - 1;
            while (i <= j) {
                int t = B->tri[i];
                double f = NB * (cen[3 * (size_t)t + ax] - cb.lo[ax]) / ext;
                int k = f >= NB ? NB - 1 : f > 0 ? (int)f : 0;
                if (k < bs) i++;
                else { B->tri[i] = B->tri[j], B->tri[j] = t, j--; }
            }
            mid = i;
        }
    }
    if (ext <= 0) mid = first + count / 2; /* centroids coincide: split by count */
    int c = B->nnode;
    B->nnode += 2;
    nd = &B->node[node];
    nd->first = c, nd->count = 0;
    if (mk_bvh_split_(B, cen, c, first, mid - first, depth + 1)) return -1;
    return mk_bvh_split_(B, cen, c + 1, mid, first + count - mid, depth + 1);
}
static inline void mk_bvh_free(mk_bvh *B)
{
    free(B->tri), free(B->node);
    B->tri = NULL, B->node = NULL, B->nnode = 0, B->nt = 0; /* 1.6.1: nt too, so a query after a failed build or a free finds no tree */
}
/* the tree keeps pointers to V and T, which must outlive it. 0 or -1. */
static inline int mk_bvh_build(mk_bvh *B, const double *V, const int *T, int nt)
{
    B->V = V, B->T = T, B->nt = nt, B->nnode = 1, B->tri = NULL, B->node = NULL;
    if (nt < 0 || nt > (__INT_MAX__ - 1) / 2) { B->nt = 0, B->nnode = 0; return -1; }
    B->tri = malloc(sizeof(int) * (size_t)(nt ? nt : 1));
    B->node = malloc(sizeof(mk_bvh_node) * (2 * (size_t)nt + 1));
    double *cen = malloc(sizeof(double) * 3 * (size_t)(nt ? nt : 1));
    if (!B->tri || !B->node || !cen) { free(cen), mk_bvh_free(B); return -1; }
    for (int t = 0; t < nt; t++) {
        B->tri[t] = t;
        for (int k = 0; k < 3; k++)
            cen[3 * (size_t)t + k] = (V[3 * (size_t)T[3 * (size_t)t] + k] + V[3 * (size_t)T[3 * (size_t)t + 1] + k] + V[3 * (size_t)T[3 * (size_t)t + 2] + k]) / 3;
    }
    int r = nt ? mk_bvh_split_(B, cen, 0, 0, nt, 0) : (mk_aabb_empty(&B->node[0].box), B->node[0].first = 0, B->node[0].count = 0, 0);
    free(cen);
    return r;
}
#define MK_BVH_STACK 256
static inline const double *mk_bvh_v_(const mk_bvh *B, int t, int j) { return B->V + 3 * (size_t)B->T[3 * (size_t)t + j]; }

/* the nearest hit in (tmin, tmax): its triangle (or -1), *t, bary */
static inline int mk_bvh_first_hit(const mk_bvh *B, const double *o, const double *d, double tmin, double tmax, double *t, double *bary)
{
    mk_ray r = mk_ray_make(o, d);
    int st[MK_BVH_STACK], sp = 0, best = -1;
    double bt = tmax, bb[3] = {0, 0, 0};
    if (!B->nt) return -1;
    st[sp++] = 0;
    while (sp) {
        const mk_bvh_node *nd = &B->node[st[--sp]];
        double te;
        if (!mk_aabb_ray(&nd->box, &r, tmin, bt, &te)) continue;
        if (nd->count) {
            for (int i = nd->first; i < nd->first + nd->count; i++) {
                int tr = B->tri[i];
                double th, b3[3];
                /* (tmin, bt]: a hit at exactly bt counts so the lower index wins a tie */
                if (mk_ray_tri(&r, mk_bvh_v_(B, tr, 0), mk_bvh_v_(B, tr, 1), mk_bvh_v_(B, tr, 2), tmin, best < 0 ? tmax : nextafter(bt, INFINITY), &th, b3) &&
                    (th < bt || (th == bt && (best < 0 || tr < best))))
                    bt = th, best = tr, bb[0] = b3[0], bb[1] = b3[1], bb[2] = b3[2];
            }
        } else if (sp + 2 <= MK_BVH_STACK) {
            double t0, t1;
            int h0 = mk_aabb_ray(&B->node[nd->first].box, &r, tmin, bt, &t0), h1 = mk_aabb_ray(&B->node[nd->first + 1].box, &r, tmin, bt, &t1);
            if (h0 && h1) { /* nearer child on top */
                if (t0 <= t1) st[sp++] = nd->first + 1, st[sp++] = nd->first;
                else st[sp++] = nd->first, st[sp++] = nd->first + 1;
            } else if (h0) st[sp++] = nd->first;
            else if (h1) st[sp++] = nd->first + 1;
        }
    }
    if (best >= 0) { *t = bt; if (bary) bary[0] = bb[0], bary[1] = bb[1], bary[2] = bb[2]; }
    return best;
}
/* is anything hit in (tmin, tmax)? (shadow rays) */
static inline int mk_bvh_any_hit(const mk_bvh *B, const double *o, const double *d, double tmin, double tmax)
{
    mk_ray r = mk_ray_make(o, d);
    int st[MK_BVH_STACK], sp = 0;
    if (!B->nt) return 0;
    st[sp++] = 0;
    while (sp) {
        const mk_bvh_node *nd = &B->node[st[--sp]];
        double te, th;
        if (!mk_aabb_ray(&nd->box, &r, tmin, tmax, &te)) continue;
        if (nd->count) {
            for (int i = nd->first; i < nd->first + nd->count; i++) {
                int tr = B->tri[i];
                if (mk_ray_tri(&r, mk_bvh_v_(B, tr, 0), mk_bvh_v_(B, tr, 1), mk_bvh_v_(B, tr, 2), tmin, tmax, &th, NULL)) return 1;
            }
        } else if (sp + 2 <= MK_BVH_STACK) st[sp++] = nd->first, st[sp++] = nd->first + 1;
    }
    return 0;
}
/* every hit with t > 0; *edge = 1 when one of them has a zero weight (on an edge or a vertex). With wind (1.5.3), *wind = the
   hits' facings summed (mk_ray_tri_s_): on a closed mesh (every edge run as often one way as the other) and a ray with no edge
   hit, the mesh's winding number at o, exactly, wherever the hits' t are not rounding-close to 0 */
static inline int mk_bvh_count_wind_(const mk_bvh *B, const double *o, const double *d, int *edge, int *wind)
{
    mk_ray r = mk_ray_make(o, d);
    int st[MK_BVH_STACK], sp = 0, n = 0, w = 0;
    *edge = 0;
    if (wind) *wind = 0;
    if (!B->nt) return 0;
    st[sp++] = 0;
    while (sp) {
        const mk_bvh_node *nd = &B->node[st[--sp]];
        double te, th, b3[3];
        if (!mk_aabb_ray(&nd->box, &r, 0, INFINITY, &te)) continue;
        if (nd->count) {
            for (int i = nd->first; i < nd->first + nd->count; i++) {
                int tr = B->tri[i], ex;
                if (mk_ray_tri_s_(&r, mk_bvh_v_(B, tr, 0), mk_bvh_v_(B, tr, 1), mk_bvh_v_(B, tr, 2), 0, INFINITY, &th, b3, &ex)) {
                    n++, w += ex;
                    if (b3[0] == 0 || b3[1] == 0 || b3[2] == 0) *edge = 1;
                }
            }
        } else if (sp + 2 <= MK_BVH_STACK) st[sp++] = nd->first, st[sp++] = nd->first + 1;
    }
    if (wind) *wind = w;
    return n;
}
/* every hit with t > 0 (1.6.0: the one traversal, mk_bvh_count_wind_'s) */
static inline int mk_bvh_count_hits(const mk_bvh *B, const double *o, const double *d)
{
    int edge;
    return mk_bvh_count_wind_(B, o, d, &edge, NULL);
}
/* the triangles whose boxes meet box q (closed: a shared face or corner counts): their count, the first cap of them in out[]
   (out may be NULL with cap 0, to count) */
static inline int mk_bvh_box(const mk_bvh *B, const mk_aabb *q, int *out, int cap)
{
    int st[MK_BVH_STACK], sp = 0, n = 0;
    if (B->nt) st[sp++] = 0;
    while (sp) {
        const mk_bvh_node *nd = &B->node[st[--sp]];
        int meet = 1;
        for (int k = 0; k < 3; k++) meet &= nd->box.lo[k] <= q->hi[k] && q->lo[k] <= nd->box.hi[k];
        if (!meet) continue;
        if (nd->count) {
            for (int i = nd->first; i < nd->first + nd->count; i++) {
                mk_aabb tb;
                mk_tri_box_(B->V, B->T, B->tri[i], &tb);
                int m = 1;
                for (int k = 0; k < 3; k++) m &= tb.lo[k] <= q->hi[k] && q->lo[k] <= tb.hi[k];
#if defined(MK_MUTATE) && MK_MUTATE == 95
                m &= tb.hi[0] > q->lo[0]; /* a box that only touches is lost */
#endif
                if (m) {
                    if (n < cap) out[n] = B->tri[i];
                    n++;
                }
            }
        } else if (sp + 2 <= MK_BVH_STACK) st[sp++] = nd->first, st[sp++] = nd->first + 1;
    }
    return n;
}
/* an edge pq strictly through triangle abc's interior (Moller-Trumbore, every weight and t strictly inside, the determinant judged
   against the lengths): meshcheck's crossing test. 1.5.2: mk_dec.h's, moved here with its arithmetic (mk_decimate unchanged). */
static inline int mk_seg_tri_strict(const double *p, const double *q, const double *a, const double *b, const double *c)
{
    double d[3], e1[3], e2[3], s[3];
    for (int k = 0; k < 3; k++) d[k] = q[k] - p[k], e1[k] = b[k] - a[k], e2[k] = c[k] - a[k], s[k] = p[k] - a[k];
    const double h[3] = {d[1] * e2[2] - d[2] * e2[1], d[2] * e2[0] - d[0] * e2[2], d[0] * e2[1] - d[1] * e2[0]};
    const double det = e1[0] * h[0] + e1[1] * h[1] + e1[2] * h[2];
    if (!(fabs(det) > 1e-15 * sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) * sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]) *
                          sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2])))
        return 0;
    const double f = 1 / det, u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
    if (u <= 0 || u >= 1) return 0;
    const double qq[3] = {s[1] * e1[2] - s[2] * e1[1], s[2] * e1[0] - s[0] * e1[2], s[0] * e1[1] - s[1] * e1[0]};
    const double v = f * (d[0] * qq[0] + d[1] * qq[1] + d[2] * qq[2]);
    if (v <= 0 || u + v >= 1) return 0;
    const double t = f * (e2[0] * qq[0] + e2[1] * qq[1] + e2[2] * qq[2]);
    return t > 0 && t < 1;
}
/* faces X and Y of one mesh cross: they share no corner, and an edge of one goes strictly through the other (1.5.2) */
static inline int mk_tri_cross(const double *V, const int *X, const int *Y)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (X[i] == Y[j]) return 0;
    for (int k = 0; k < 3; k++)
        if (mk_seg_tri_strict(V + 3 * (size_t)X[k], V + 3 * (size_t)X[(k + 1) % 3], V + 3 * (size_t)Y[0], V + 3 * (size_t)Y[1], V + 3 * (size_t)Y[2]) ||
            mk_seg_tri_strict(V + 3 * (size_t)Y[k], V + 3 * (size_t)Y[(k + 1) % 3], V + 3 * (size_t)X[0], V + 3 * (size_t)X[1], V + 3 * (size_t)X[2]))
            return 1;
    return 0;
}
/* the crossing face pairs of B's mesh (mk_tri_cross), each counted once; with mark (a byte a face, not NULL) only the pairs with a
   marked face, against every face near it: a local count, near a seam or a changed region (1.5.2, the maintainer's race-c1 H3). -1: no
   memory. */
static inline long mk_bvh_crossings(const mk_bvh *B, const unsigned char *mark)
{
    int cap = 64, *buf = malloc(sizeof *buf * (size_t)cap);
    long np = 0;
    if (!buf) return -1;
    for (int t = 0; t < B->nt; t++) {
        if (mark && !mark[t]) continue;
        mk_aabb bx;
        mk_tri_box_(B->V, B->T, t, &bx);
        int n = mk_bvh_box(B, &bx, buf, cap);
        if (n > cap) {
            int *b2 = realloc(buf, sizeof *b2 * (size_t)n);
            if (!b2) { free(buf); return -1; }
            buf = b2, cap = n, n = mk_bvh_box(B, &bx, buf, cap);
        }
        for (int i = 0; i < n; i++) {
            const int u = buf[i];
#if defined(MK_MUTATE) && MK_MUTATE == 206
            if (u <= t) continue; /* a marked face's pairs with a lower unmarked face lost */
#else
            if (u == t || (u < t && (!mark || mark[u]))) continue; /* counted from u */
#endif
            np += mk_tri_cross(B->V, B->T + 3 * (size_t)t, B->T + 3 * (size_t)u);
        }
    }
    free(buf);
    return np;
}
/* bd widened by the rounding a triangle's computed squared distance may carry (16 u of the coordinates in play, here those of
   the root's box and p, which bound every node's), so ties are not pruned. e once a query, the bound only when bd changes. */
static inline double mk_near_slack_(double bd, double e) { return bd + 2 * sqrt(bd) * e + e * e; }
/* the closest point on the mesh to p within sqrt(maxd2): returns its squared distance (maxd2 if none), *tri (-1 if none),
   q the point. mk_bvh_nearest_lab takes only the triangles t with lab[t] >= 0 and lab[t] != skip (lab NULL: all of them). */
static inline double mk_bvh_nearest_lab(const mk_bvh *B, const double *p, double maxd2, const int *lab, int skip, int *tri, double *q)
{
    int st[MK_BVH_STACK], sp = 0, best = -1;
    double bd = maxd2, e = 0, bs;
    if (B->nt) {
        double x = 0;
        for (int k = 0; k < 3; k++) x = fmax(x, fmax(fabs(p[k]), fmax(fabs(B->node[0].box.lo[k]), fabs(B->node[0].box.hi[k]))));
        e = 0x1p-49 * x, st[sp++] = 0;
    }
    bs = mk_near_slack_(bd, e);
    while (sp) {
        const mk_bvh_node *nd = &B->node[st[--sp]];
#if defined(MK_MUTATE) && MK_MUTATE == 93
        if ((void)bs, mk_aabb_dist2(&nd->box, p) > 0.25 * bd) continue;
#else
        if (mk_aabb_dist2(&nd->box, p) > bs) continue;
#endif
        if (nd->count) {
            for (int i = nd->first; i < nd->first + nd->count; i++) {
                int tr = B->tri[i];
                if (lab && (lab[tr] < 0 || lab[tr] == skip)) continue;
                double w[3], c[3], d2 = mk_closest_pt_tri(p, mk_bvh_v_(B, tr, 0), mk_bvh_v_(B, tr, 1), mk_bvh_v_(B, tr, 2), w, c);
                if (d2 < bd || (d2 == bd && best >= 0 && tr < best)) bd = d2, bs = mk_near_slack_(bd, e), best = tr, q[0] = c[0], q[1] = c[1], q[2] = c[2];
            }
        } else if (sp + 2 <= MK_BVH_STACK) {
            double d0 = mk_aabb_dist2(&B->node[nd->first].box, p), d1 = mk_aabb_dist2(&B->node[nd->first + 1].box, p);
            if (d0 <= d1) st[sp++] = nd->first + 1, st[sp++] = nd->first; /* nearer child on top */
            else st[sp++] = nd->first, st[sp++] = nd->first + 1;
        }
    }
    *tri = best;
    return bd;
}
static inline double mk_bvh_nearest(const mk_bvh *B, const double *p, double maxd2, int *tri, double *q)
{
    return mk_bvh_nearest_lab(B, p, maxd2, NULL, 0, tri, q);
}
static inline double mk_winding_number(const double *V, const int *T, int nt, const double *p);
/* inside a closed mesh: ray parity along three fixed, unaligned directions. A ray with a hit on an edge or a vertex may count
   one crossing twice, so it does not vote; the clean rays decide when they agree, otherwise the winding number does. */
static inline int mk_bvh_inside(const mk_bvh *B, const double *p)
{
    static const double dirs[3][3] = {{0.5773502691896258, 0.5773502691896258, 0.5773502691896258},
                                      {-0.3015113445777636, 0.9045340337332909, 0.3015113445777636},
                                      {0.2672612419124244, -0.5345224838248488, 0.8017837257372732}};
    int in = 0, out = 0;
    for (int k = 0; k < 3; k++) {
        int edge, c = mk_bvh_count_wind_(B, p, dirs[k], &edge, NULL);
        if (!edge) in += c & 1, out += !(c & 1);
    }
    if (in && !out) return 1;
    if (out && !in) return 0;
    return fabs(mk_winding_number(B->V, B->T, B->nt, p)) > 0.5; /* 1.6.1: either winding, as the parity and mk_sdf (1.6.0: outward only) */
}
/* generalized winding number of the mesh at p (sum of solid angles, Van Oosterom and Strackee; gauge's meshtop): 1 inside a
   closed outward surface, 0 outside, fractional for an open one. O(triangles). */
static inline double mk_winding_number(const double *V, const int *T, int nt, const double *p)
{
    double s = 0;
    for (int t = 0; t < nt; t++) {
        double A[3], Bv[3], C[3];
        for (int k = 0; k < 3; k++) A[k] = V[3 * (size_t)T[3 * (size_t)t] + k] - p[k], Bv[k] = V[3 * (size_t)T[3 * (size_t)t + 1] + k] - p[k], C[k] = V[3 * (size_t)T[3 * (size_t)t + 2] + k] - p[k];
        double la = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]), lb = sqrt(Bv[0] * Bv[0] + Bv[1] * Bv[1] + Bv[2] * Bv[2]), lc = sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);
        double num = A[0] * (Bv[1] * C[2] - Bv[2] * C[1]) - A[1] * (Bv[0] * C[2] - Bv[2] * C[0]) + A[2] * (Bv[0] * C[1] - Bv[1] * C[0]);
        double den = la * lb * lc + (A[0] * Bv[0] + A[1] * Bv[1] + A[2] * Bv[2]) * lc + (Bv[0] * C[0] + Bv[1] * C[1] + Bv[2] * C[2]) * la + (C[0] * A[0] + C[1] * A[1] + C[2] * A[2]) * lb;
        s += 2 * atan2(num, den);
    }
#if defined(MK_MUTATE) && MK_MUTATE == 94
    return -s / (4 * 3.14159265358979323846);
#else
    return s / (4 * 3.14159265358979323846);
#endif
}

/* Many queries at once, on OpenMP threads when the caller builds with -fopenmp (mk_par.h): blocks of MK_GEOM_QBLK points
   handed out as threads come free. Each answer is the one the single call gives for that point, bit for bit. */
enum { MK_GEOM_QBLK = 64 };
/* the serial forms on points [i0, i1) */
static inline void mk_bvh_nearest_blk_(const mk_bvh *B, const double *P, double maxd2, double *d2, int *tri, double *q, long i0, long i1)
{
    for (long i = i0; i < i1; i++) {
        double c[3];
        int t;
        const double d = mk_bvh_nearest(B, P + 3 * i, maxd2, &t, c);
        if (d2) d2[i] = d;
        if (tri) tri[i] = t;
        if (q && t >= 0) q[3 * i] = c[0], q[3 * i + 1] = c[1], q[3 * i + 2] = c[2];
    }
}
static inline void mk_bvh_inside_blk_(const mk_bvh *B, const double *P, int *in, long i0, long i1)
{
    for (long i = i0; i < i1; i++) in[i] = mk_bvh_inside(B, P + 3 * i);
}
static inline void mk_winding_blk_(const double *V, const int *T, int nt, const double *P, double *w, long i0, long i1)
{
    for (long i = i0; i < i1; i++) w[i] = mk_winding_number(V, T, nt, P + 3 * i);
}
/* for n points P[3 n]: d2[i] and tri[i] as mk_bvh_nearest gives them for P + 3 i, and q[3 i] the point when tri[i] >= 0 (left as it
   was otherwise, as the single call leaves it); any of d2, tri, q may be NULL */
static inline void mk_bvh_nearest_many(const mk_bvh *B, const double *P, long n, double maxd2, double *d2, int *tri, double *q)
{
#if defined(MK_MUTATE) && MK_MUTATE == 97
    MK_PAR_CUT_(n, i0, i1, mk_bvh_nearest_blk_(B, P, maxd2, d2, tri, q, i0, i1)); /* 97: points cut by the thread count, each piece's last point left out */
#else
    const long nb = n > 0 ? (n + MK_GEOM_QBLK - 1) / MK_GEOM_QBLK : 0;
    MK_OMP_(omp parallel for schedule(dynamic, 1) if(nb > 1))
    for (long b = 0; b < nb; b++) mk_bvh_nearest_blk_(B, P, maxd2, d2, tri, q, b * MK_GEOM_QBLK, (b + 1) * MK_GEOM_QBLK < n ? (b + 1) * MK_GEOM_QBLK : n);
#endif
}
/* in[i] = mk_bvh_inside(B, P + 3 i) */
static inline void mk_bvh_inside_many(const mk_bvh *B, const double *P, long n, int *in)
{
#if defined(MK_MUTATE) && MK_MUTATE == 98
    MK_PAR_CUT_(n, i0, i1, mk_bvh_inside_blk_(B, P, in, i0, i1)); /* 98: points cut by the thread count, each piece's last point left out */
#else
    const long nb = n > 0 ? (n + MK_GEOM_QBLK - 1) / MK_GEOM_QBLK : 0;
    MK_OMP_(omp parallel for schedule(dynamic, 1) if(nb > 1))
    for (long b = 0; b < nb; b++) mk_bvh_inside_blk_(B, P, in, b * MK_GEOM_QBLK, (b + 1) * MK_GEOM_QBLK < n ? (b + 1) * MK_GEOM_QBLK : n);
#endif
}
/* w[i] = mk_winding_number(V, T, nt, P + 3 i), O(triangles) each, so a point a block */
static inline void mk_winding_number_many(const double *V, const int *T, int nt, const double *P, long n, double *w)
{
#if defined(MK_MUTATE) && MK_MUTATE == 96
    for (long i = 0; i < n; i++) { /* 96: each point's sum split over the threads by OpenMP's reduction(+:), its order following their count */
        const double *p = P + 3 * i;
        double s = 0;
        MK_OMP_(omp parallel for schedule(static) reduction(+:s))
        for (int t = 0; t < nt; t++) {
            double A[3], Bv[3], C[3];
            for (int k = 0; k < 3; k++) A[k] = V[3 * (size_t)T[3 * (size_t)t] + k] - p[k], Bv[k] = V[3 * (size_t)T[3 * (size_t)t + 1] + k] - p[k], C[k] = V[3 * (size_t)T[3 * (size_t)t + 2] + k] - p[k];
            double la = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]), lb = sqrt(Bv[0] * Bv[0] + Bv[1] * Bv[1] + Bv[2] * Bv[2]), lc = sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);
            double num = A[0] * (Bv[1] * C[2] - Bv[2] * C[1]) - A[1] * (Bv[0] * C[2] - Bv[2] * C[0]) + A[2] * (Bv[0] * C[1] - Bv[1] * C[0]);
            double den = la * lb * lc + (A[0] * Bv[0] + A[1] * Bv[1] + A[2] * Bv[2]) * lc + (Bv[0] * C[0] + Bv[1] * C[1] + Bv[2] * C[2]) * la + (C[0] * A[0] + C[1] * A[1] + C[2] * A[2]) * lb;
            s += 2 * atan2(num, den);
        }
        w[i] = s / (4 * 3.14159265358979323846);
    }
#else
    MK_OMP_(omp parallel for schedule(dynamic, 1) if(n > 1))
    for (long i = 0; i < n; i++) mk_winding_blk_(V, T, nt, P, w, i, i + 1);
#endif
}
#endif
