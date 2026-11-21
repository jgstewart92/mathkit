/* mathkit: the volume two closed meshes share, and what each has that the other has not, by columns. Lines run along
   one axis through a lattice of cells of side h on the other two, set off by (0.5371, 0.5213) of a cell so no line runs
   along a lattice-aligned edge. Each mesh's crossings of a line are found one after another (mk_geom's watertight first
   hit; a hit at the same parameter as the last, a shared edge or vertex, counts once); a line is inside between an odd
   and the next even crossing. The shared length, and each mesh's own, are summed times h^2. A line with an odd count of
   either mesh (a tangent touch, a hole), or more than MK_OVL_CAP crossings, is left out and counted. The lattice spans
   the two meshes' common box, or both boxes with MK_OVL_ALL (needed for only_a, only_b and the volumes). Accuracy: a
   smooth solid's volume to O(h^2) relative, with an error of about the boundary length times h in a cross-section when
   a face lies along the line (a box aligned with the axis). Header only, C23.
   Experimental: its results may change between releases. */
#ifndef MK_OVL_H
#define MK_OVL_H
#include "mk_geom.h"
#include <math.h>
#include <stdlib.h>

#define MK_OVL_ALL 1
#define MK_OVL_CAP 4096
typedef struct {
    double both;           /* the volume inside both */
    /* each mesh's volume by the same columns (MK_OVL_ALL; else only on lines that cross B, each mesh's whole length
       along them) */
    double vol_a, vol_b;
    double only_a, only_b; /* inside one and not the other */
    long cols, odd;        /* lines that cross a mesh (without MK_OVL_ALL: that cross B); those left out */
} mk_ovl_report;

/* the crossings of the line o + t e_axis, t in [-len, len], of tree B, one at a time; n, or -1 past cap. The range is
   widened by 2^-48 len, past t's rounding, so a crossing at an end is found. */
static inline int mk_ovl_hits_(const mk_bvh *B, const double *o, int axis, double len, double *t, int cap)
{
    const double lim = len + len * 0x1p-48 + 0x1p-1060;
#if defined(MK_MUTATE) && MK_MUTATE == 226
    double d[3] = {0, 0, 0}, tm = 0, th, b[3]; /* planted fault: crossings from t = 0 only */
#else
    double d[3] = {0, 0, 0}, tm = -lim, th, b[3];
#endif
    int n = 0;
    d[axis] = 1;
    while (mk_bvh_first_hit(B, o, d, tm, lim, &th, b) >= 0) {
        if (n == cap) return -1;
        t[n++] = th;
        tm = nextafter(th, INFINITY);
    }
    return n;
}
/* the length inside both of two sorted crossing lists (inside between an odd and the next even one) */
static inline double mk_ovl_and_(const double *ta, int na, const double *tb, int nb)
{
    double s = 0;
    for (int p = 0; p + 1 < na; p += 2)
        for (int q = 0; q + 1 < nb; q += 2) {
#if defined(MK_MUTATE) && MK_MUTATE == 171
            /* planted fault: the union, not the overlap */
            const double lo = fmin(ta[p], tb[q]), hi = fmax(ta[p + 1], tb[q + 1]);
            if (hi > lo && ta[p + 1] > tb[q] && tb[q + 1] > ta[p]) s += hi - lo;
#else
            const double lo = fmax(ta[p], tb[q]), hi = fmin(ta[p + 1], tb[q + 1]);
            if (hi > lo) s += hi - lo;
#endif
        }
    return s;
}
static inline double mk_ovl_len_(const double *t, int n)
{
    double s = 0;
    for (int p = 0; p + 1 < n; p += 2) s += t[p + 1] - t[p];
    return s;
}
/* columns along axis (0, 1, 2) through trees A and B (closed meshes: mk_bvh_build's), cells of side h; flags MK_OVL_ALL
   or 0. 0, or -1 for h not positive and finite, a bad axis, an empty tree, more than 2^40 lines, or no memory; r is
   then zeroed. */
static inline int mk_ovl_bvh(const mk_bvh *A, const mk_bvh *B, double h, int axis, int flags, mk_ovl_report *r)
{
    mk_ovl_report z = {0};
    *r = z;
    if (!A || !B || !(h > 0) || !isfinite(h) || axis < 0 || axis > 2 || !A->nt || !B->nt) return -1;
    const mk_aabb *ba = &A->node[0].box, *bb = &B->node[0].box;
    const int u = (axis + 1) % 3, w = (axis + 2) % 3;
    double lo[3], hi[3], all_lo = fmin(ba->lo[axis], bb->lo[axis]), all_hi = fmax(ba->hi[axis], bb->hi[axis]);
    for (int k = 0; k < 3; k++) {
        lo[k] = flags & MK_OVL_ALL ? fmin(ba->lo[k], bb->lo[k]) : fmax(ba->lo[k], bb->lo[k]);
        hi[k] = flags & MK_OVL_ALL ? fmax(ba->hi[k], bb->hi[k]) : fmin(ba->hi[k], bb->hi[k]);
        if (!isfinite(lo[k]) || !isfinite(hi[k])) return -1;
    }
    if (!(lo[u] < hi[u]) || !(lo[w] < hi[w]) || !(lo[axis] < hi[axis])) return 0; /* no common box: nothing shared */
    const double i0 = floor(lo[u] / h) - 1, i1 = ceil(hi[u] / h) + 1, j0 = floor(lo[w] / h) - 1, j1 = ceil(hi[w] / h) + 1;
    if (!((i1 - i0 + 1) * (j1 - j0 + 1) <= 0x1p40)) return -1;
    /* the lines' origin: the point of the axis span nearest 0, so a crossing's t is as precise as its coordinate, and
       no span overflows */
#if defined(MK_MUTATE) && MK_MUTATE == 224
    /* planted fault: the origin below both meshes */
    const double pad = (all_hi - all_lo) * 1e-6 + h, za = all_lo - pad, len = all_hi + pad - za;
#else
    const double za = fmin(fmax(all_lo, 0.0), all_hi), len = fmax(all_hi - za, za - all_lo);
#endif
    double *ta = malloc(sizeof *ta * 2 * MK_OVL_CAP), *tb = ta ? ta + MK_OVL_CAP : NULL;
    if (!ta) return -1;
    double both = 0, va = 0, vb = 0, oa = 0, ob = 0;
    long cols = 0, odd = 0;
#if defined(MK_MUTATE) && MK_MUTATE == 222
    for (double i = i0; i <= i1; i++) /* planted fault: a double counter, which stops moving past 2^53 */
        for (long kj = 0; kj <= (long)(j1 - j0); kj++) {
            const double j = j0 + (double)kj;
#elif defined(MK_MUTATE) && MK_MUTATE == 223
    for (long ki = 0; ki <= (long)(i1 - i0); ki++)
        for (double j = j0; j <= j1; j++) { /* planted fault: the same on the other axis */
            const double i = i0 + (double)ki;
#else
    for (long ki = 0; ki <= (long)(i1 - i0); ki++)
        for (long kj = 0; kj <= (long)(j1 - j0); kj++) {
            const double i = i0 + (double)ki, j = j0 + (double)kj;
#endif
            double o[3];
#if defined(MK_MUTATE) && MK_MUTATE == 172
            o[u] = i * h, o[w] = j * h, o[axis] = za; /* planted fault: lines on the lattice itself */
#else
            o[u] = (i + 0.5371) * h, o[w] = (j + 0.5213) * h, o[axis] = za;
#endif
            const int nb = mk_ovl_hits_(B, o, axis, len, tb, MK_OVL_CAP);
            if (!nb && !(flags & MK_OVL_ALL)) continue;
            const int na = mk_ovl_hits_(A, o, axis, len, ta, MK_OVL_CAP);
            if (!na && !nb) continue;
            cols++;
#if defined(MK_MUTATE) && MK_MUTATE == 173
            if (na < 0 || nb < 0) { odd++; continue; } /* planted fault: an odd count used as it is */
#else
            if (na < 0 || nb < 0 || (na & 1) || (nb & 1)) { odd++; continue; }
#endif
            const double s = mk_ovl_and_(ta, na, tb, nb), la = mk_ovl_len_(ta, na), lb = mk_ovl_len_(tb, nb);
            both += s, va += la, vb += lb, oa += la - s, ob += lb - s;
        }
    free(ta);
#if defined(MK_MUTATE) && MK_MUTATE == 225
    const double c = h * h; /* planted fault: the cell's area first, which overflows or goes subnormal */
    r->both = both * c, r->vol_a = va * c, r->vol_b = vb * c, r->only_a = oa * c, r->only_b = ob * c, r->cols = cols, r->odd = odd;
#else
    r->both = both * h * h, r->vol_a = va * h * h, r->vol_b = vb * h * h, r->only_a = oa * h * h, r->only_b = ob * h * h, r->cols = cols, r->odd = odd;
#endif
    return 0;
}
/* the same from the meshes (VA[3 nva], TA[3 nta]; VB, TB): builds and frees the two trees. -1 also for a size under 1,
   an index out of range or a coordinate not finite. */
static inline int mk_overlap_volume(const double *VA, int nva, const int *TA, int nta, const double *VB, int nvb, const int *TB, int ntb, double h,
                                    int axis, int flags, mk_ovl_report *r)
{
    mk_ovl_report z = {0};
    *r = z;
    if (!VA || !TA || !VB || !TB || nva < 1 || nta < 1 || nvb < 1 || ntb < 1) return -1;
    for (int m = 0; m < 2; m++) {
        const double *V = m ? VB : VA;
        const int *T = m ? TB : TA, nv = m ? nvb : nva, nt = m ? ntb : nta;
        for (size_t q = 0; q < 3 * (size_t)nv; q++)
            if (!isfinite(V[q])) return -1;
        for (size_t q = 0; q < 3 * (size_t)nt; q++)
            if (T[q] < 0 || T[q] >= nv) return -1;
    }
    mk_bvh A, B;
    if (mk_bvh_build(&A, VA, TA, nta)) return -1;
    if (mk_bvh_build(&B, VB, TB, ntb)) { mk_bvh_free(&A); return -1; }
    const int ret = mk_ovl_bvh(&A, &B, h, axis, flags, r);
    mk_bvh_free(&A), mk_bvh_free(&B);
    return ret;
}
#endif
