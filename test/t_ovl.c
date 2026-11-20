/* mk_overlap_volume (mesher-fable's ovl, with stepcad's voldiff): two spheres' lens against its analytic volume, a sphere's own volume,
   a box pair whose faces lie on the lattice (exact), the three axes, disjoint and nested pairs, an open mesh's odd lines, refusals. */
#include "../include/mk_ovl.h"
#include "t.h"
#include "tmesh.h"

int main(void)
{
    const double R = 1, d = 1.2, h = 0.02, c0[3] = {0, 0, 0}, c1[3] = {d, 0, 0}, c2[3] = {5, 0, 0}, c3[3] = {0.1, 0.05, 0};
    tm_mesh A = {0}, B = {0}, D = {0}, S = {0};
    tm_sphere(&A, R, c0, 5), tm_sphere(&B, R, c1, 5), tm_sphere(&D, R, c2, 5), tm_sphere(&S, 0.4, c3, 4);
    mk_ovl_report r, r1, r2, ra;
    const double lens = T_PI * (4 * R + d) * (2 * R - d) * (2 * R - d) / 12, va = tm_volume(A.V, A.T, A.nt);
    int e = mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, h, 2, 0, &r);
    CHECK(e == 0 && fabs(r.both / lens - 1) < 2e-3 && r.odd == 0, "two unit spheres 1.2 apart: lens %.5f against %.5f analytic (%.3f %%), %ld lines, %ld odd", r.both, lens,
          100 * (r.both / lens - 1), r.cols, r.odd);
    e = mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, h, 2, MK_OVL_ALL, &ra);
    CHECK(e == 0 && fabs(ra.vol_a / va - 1) < 2e-3 && ra.both == r.both, "the sphere's own volume by columns %.5f against its mesh %.5f; the lens the same over both boxes", ra.vol_a, va);
    CHECK(fabs(ra.only_a + ra.both - ra.vol_a) <= 1e-12 * ra.vol_a && fabs(ra.only_b + ra.both - ra.vol_b) <= 1e-12 * ra.vol_b && fabs(ra.only_a - ra.only_b) < 1e-3 * ra.only_a,
          "only_a + both = vol_a, only_b + both = vol_b; the pair symmetric (only_a %.5f, only_b %.5f)", ra.only_a, ra.only_b);
    mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, h, 0, 0, &r1), mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, h, 1, 0, &r2);
    CHECK(fabs(r1.both / lens - 1) < 2e-3 && fabs(r2.both / lens - 1) < 2e-3, "along x and along y: %.5f, %.5f", r1.both, r2.both);
    e = mk_overlap_volume(A.V, A.nv, A.T, A.nt, D.V, D.nv, D.T, D.nt, h, 2, 0, &r);
    CHECK(e == 0 && r.both == 0 && r.cols == 0, "disjoint spheres: 0, no lines (%g, %ld)", r.both, r.cols);
    e = mk_overlap_volume(A.V, A.nv, A.T, A.nt, S.V, S.nv, S.T, S.nt, h, 2, 0, &r);
    const double vs = tm_volume(S.V, S.T, S.nt);
    CHECK(e == 0 && fabs(r.both / vs - 1) < 5e-3, "a sphere inside another: %.5f, its volume %.5f", r.both, vs);
    { /* boxes with faces on multiples of h: every line misses every edge, the count is exact */
        tm_mesh P = {0}, Q = {0};
        const double lo0[3] = {0, 0, 0}, s0[3] = {1, 1, 1}, lo1[3] = {0.5, 0.25, 0.25}, s1[3] = {1, 1, 1};
        tm_box(&P, lo0, s0, 4), tm_box(&Q, lo1, s1, 4);
        e = mk_overlap_volume(P.V, P.nv, P.T, P.nt, Q.V, Q.nv, Q.T, Q.nt, 0.01, 2, MK_OVL_ALL, &r);
        CHECK(e == 0 && fabs(r.both - 0.28125) < 1e-9 && fabs(r.vol_a - 1) < 1e-9 && r.odd == 0, "two unit boxes 0.5, 0.25, 0.25 apart: %.12f (0.28125), each %.12f, %ld odd", r.both,
              r.vol_a, r.odd);
        tm_free(&P), tm_free(&Q);
    }
    { /* an open mesh (a face less): its lines through the hole are left out and counted */
        tm_mesh O = {0};
        tm_sphere(&O, R, c1, 5);
        int best = 0;
        double bd = INFINITY;
        for (int t = 0; t < O.nt; t++) { /* a face of B inside A, facing up as much as along x */
            double g[3] = {0, 0, 0}, dd = 0;
            for (int k = 0; k < 3; k++) g[k] = (O.V[3 * O.T[3 * t] + k] + O.V[3 * O.T[3 * t + 1] + k] + O.V[3 * O.T[3 * t + 2] + k]) / 3;
            dd = (g[0] - 0.4) * (g[0] - 0.4) + g[1] * g[1] + (g[2] - 0.6) * (g[2] - 0.6);
            if (dd < bd) bd = dd, best = t;
        }
        (void)best;
        for (int t = 0; t < O.nt; t++) { /* every face within 0.08 of that point */
            double g[3];
            for (int k = 0; k < 3; k++) g[k] = (O.V[3 * O.T[3 * t] + k] + O.V[3 * O.T[3 * t + 1] + k] + O.V[3 * O.T[3 * t + 2] + k]) / 3;
            if ((g[0] - 0.4) * (g[0] - 0.4) + g[1] * g[1] + (g[2] - 0.6) * (g[2] - 0.6) < 0.0064) memcpy(O.T + 3 * t, O.T + 3 * (O.nt - 1), 3 * sizeof *O.T), O.nt--, t--;
        }
        e = mk_overlap_volume(A.V, A.nv, A.T, A.nt, O.V, O.nv, O.T, O.nt, h, 2, 0, &r);
        CHECK(e == 0 && r.odd > 0 && fabs(r.both / lens - 1) < 0.03, "one sphere with a hole 0.16 across inside the other: %ld odd lines left out and counted, the lens %.5f (%.2f %%)", r.odd, r.both, 100 * (r.both / lens - 1));
        tm_free(&O);
    }
    {
        int k = 0, n = 0;
        double bad[9] = {0, 0, 0, 1, 0, 0, 0, NAN, 0};
        const int T1[3] = {0, 1, 2}, Tb[3] = {0, 1, 3};
#define REF(call) (k++, n += (call) == -1 && r.both == 0)
        REF(mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, 0, 2, 0, &r));
        REF(mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, NAN, 2, 0, &r));
        REF(mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, INFINITY, 2, 0, &r));
        REF(mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, h, 3, 0, &r));
        REF(mk_overlap_volume(A.V, A.nv, A.T, A.nt, bad, 3, T1, 1, h, 2, 0, &r));
        REF(mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, 3, Tb, 1, h, 2, 0, &r));
        REF(mk_overlap_volume(A.V, 0, A.T, A.nt, B.V, B.nv, B.T, B.nt, h, 2, 0, &r));
        REF(mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, 0, h, 2, 0, &r));
        REF(mk_overlap_volume(A.V, A.nv, A.T, A.nt, B.V, B.nv, B.T, B.nt, 1e-12, 2, 0, &r));
#undef REF
        CHECK(n == k, "refusals: %d of %d (h 0, NaN, inf, axis 3, a NaN vertex, an index out of range, no vertices, no faces, 2^40 lines)", n, k);
    }
    tm_free(&A), tm_free(&B), tm_free(&D), tm_free(&S);
    return wrong;
}
