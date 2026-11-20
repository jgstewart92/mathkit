#define _DEFAULT_SOURCE /* mmap's MAP_ANONYMOUS and MAP_NORESERVE, for the I-C row */
#include <stdlib.h>
static size_t t_asked; /* what the header's last allocation asked for; t_refuse makes the next one fail */
static int t_refuse;
static void *t_malloc(size_t s)
{
    t_asked = s;
    if (t_refuse) { t_refuse = 0; return NULL; }
    return malloc(s);
}
#define malloc t_malloc
#include "../include/mk_interp.h"
#include "../include/mk_vec.h"
#undef malloc
#include <sys/mman.h>
#include "t.h"
#include <string.h>

/* the reference: aeromap.c's code (aeromap/aeromap.c, lines 306-340, 2026-11-13), copied verbatim but for names */
typedef struct { int first, interval, last; double t; } Locate;
static int sign(double x) { return (x > 0) - (x < 0); }
static double interior_slope(double d0, double d1, double h0, double h1)
{
    if (d0 * d1 <= 0) return 0;
    double w1 = 2 * h1 + h0, w2 = h1 + 2 * h0;
    return (w1 + w2) / (w1 / d0 + w2 / d1);
}
static double end_slope(double d0, double d1, double h0, double h1)
{
    double m = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
    if (sign(m) != sign(d0)) return 0;
    if (sign(d0) != sign(d1) && fabs(m) > 3 * fabs(d0)) return 3 * d0;
    return m;
}
static double pchip(const double *x, int n, const Locate *l, const double *y)
{
    int i = l->interval;
#define Y(j) y[(j) - l->first]
    double h = x[i + 1] - x[i], d = (Y(i + 1) - Y(i)) / h, m0, m1;
    if (n == 2) return Y(i) + l->t * (Y(i + 1) - Y(i));
    if (i > 0) m0 = interior_slope((Y(i) - Y(i - 1)) / (x[i] - x[i - 1]), d, x[i] - x[i - 1], h);
    else m0 = end_slope(d, (Y(2) - Y(1)) / (x[2] - x[1]), h, x[2] - x[1]);
    if (i + 2 < n) m1 = interior_slope(d, (Y(i + 2) - Y(i + 1)) / (x[i + 2] - x[i + 1]), h, x[i + 2] - x[i + 1]);
    else m1 = end_slope(d, (Y(i) - Y(i - 1)) / (x[i] - x[i - 1]), h, x[i] - x[i - 1]);
#undef Y
    double t = l->t, t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * y[i - l->first] + (t3 - 2 * t2 + t) * h * m0
         + (-2 * t3 + 3 * t2) * y[i + 1 - l->first] + (t3 - t2) * h * m1;
}
static double ref(const double *x, const double *y, int n, double q) /* am_eval's locate, one axis, y over all nodes */
{
    if (q < x[0]) q = x[0];
    if (q > x[n - 1]) q = x[n - 1];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) { int mid = (lo + hi) / 2; if (x[mid] <= q) lo = mid; else hi = mid; }
    Locate l = {0, lo, n - 1, (q - x[lo]) / (x[lo + 1] - x[lo])};
    return pchip(x, n, &l, y);
}

int main(void)
{
    long diff = 0, diff2 = 0, total = 0;
    for (int n = 2; n <= 12; n++)
        for (int rep = 0; rep < 50; rep++) {
            double x[12], y[12], m[12];
            x[0] = t_rand();
            for (int i = 1; i < n; i++) x[i] = x[i - 1] + 0.05 + fabs(t_rand());
            for (int i = 0; i < n; i++) y[i] = rep % 3 ? t_rand() : floor(3 * t_rand());
            mk_pchip_slopes(x, y, n, m);
            for (int k = 0; k < 40; k++) {
                double q = x[0] - 0.1 + (x[n - 1] - x[0] + 0.2) * (k + 0.5) / 40;
                double r = ref(x, y, n, q), a = mk_pchip_eval(x, y, n, q), b = mk_hermite_eval(x, y, m, n, q);
                diff += memcmp(&r, &a, sizeof r) != 0, diff2 += memcmp(&r, &b, sizeof r) != 0, total++;
            }
        }
    CHECK(diff == 0, "PCHIP equals aeromap.c's bit for bit: %ld of %ld queries differ (n 2..12, ties, outside the range)", diff, total);
    CHECK(diff2 == 0, "slope table plus Hermite equals it too: %ld of %ld differ", diff2, total);
    double xm[6] = {0, 1, 2, 3, 4, 5}, ym[6] = {0, 0.1, 0.1, 2, 2.1, 5}, mm[6], prev = -1;
    mk_pchip_slopes(xm, ym, 6, mm);
    int mono = 1;
    for (int k = 0; k <= 5000; k++) {
        double v = mk_hermite_eval(xm, ym, mm, 6, k * 0.001);
        if (v < prev - 1e-15) mono = 0;
        prev = v;
    }
    CHECK(mono, "monotone data with a flat step: the interpolant never turns down over 5001 points");
    enum { S = 41 };
    double xs[S], ys[S], ms[S], es = 0;
    for (int i = 0; i < S; i++) xs[i] = T_PI * i / (S - 1), ys[i] = sin(xs[i]);
    mk_spline_slopes(xs, ys, S, ms);
    for (int k = 0; k <= 1000; k++) {
        double q = T_PI * k / 1000;
        es = fmax(es, fabs(mk_hermite_eval(xs, ys, ms, S, q) - sin(q)));
    }
    double hq = 1e-4, c0 = (mk_hermite_eval(xs, ys, ms, S, 2 * hq) - 2 * mk_hermite_eval(xs, ys, ms, S, hq) + mk_hermite_eval(xs, ys, ms, S, 0)) / (hq * hq);
    CHECK(es < 2e-6 && fabs(c0) < 1e-2, "natural spline on sin over [0, pi], 41 nodes: error %.2g, curvature at 0 %.2g (0)", es, c0);
    double xl[3] = {0, 1, 3}, yl[3] = {1, 3, 7}, ml[3];
    mk_spline_slopes(xl, yl, 3, ml);
    CHECK(fabs(mk_hermite_eval(xl, yl, ml, 3, 2) - 5) < 1e-14 && fabs(mk_linear_eval(xl, yl, 3, 2) - 5) < 1e-15, "a straight line stays straight: spline and linear give 5 at x 2");
    { /* 1.0.5: the verification work's interp findings (Result 79) */
        double xb[5] = {NAN, NAN, 1, NAN, NAN}, yb[5] = {NAN, NAN, 2, NAN, NAN}, mb[5] = {7, 7, 7, 7, 7}, ms5[5] = {7, 7, 7, 7, 7};
        double *x1 = xb + 2, *y1 = yb + 2; /* one point; NaN on either side stands for memory past the table */
        double v0 = mk_linear_eval(x1, y1, 1, 0.5), v1 = mk_pchip_eval(x1, y1, 1, 3);
        mk_pchip_slopes(x1, y1, 1, mb + 2);
        int s1 = mk_spline_slopes(x1, y1, 1, ms5 + 2);
        double v2 = mk_hermite_eval(x1, y1, mb + 2, 1, 1), v3 = mk_linear_eval(x1, y1, 0, 1);
        int kept = mb[0] == 7 && mb[1] == 7 && mb[3] == 7 && mb[4] == 7 && ms5[0] == 7 && ms5[1] == 7 && ms5[3] == 7 && ms5[4] == 7;
        CHECK(v0 == 2 && v1 == 2 && v2 == 2 && mb[2] == 0 && s1 == 0 && ms5[2] == 0 && kept && isnan(v3),
              "I-A: a one-point table: linear, PCHIP and Hermite give y0 (%g %g %g), both slopes 0 (%g %g), nothing past it read or written; n 0 gives NaN",
              v0, v1, v2, mb[2], ms5[2]);
        double xs3[3] = {0, 1, 2}, ys3[3] = {0, 1e-170, 2e-170}, ms3[3];
        mk_pchip_slopes(xs3, ys3, 3, ms3);
        double vb = mk_pchip_eval(xs3, ys3, 3, 0.5), vh = mk_hermite_eval(xs3, ys3, ms3, 3, 0.5);
        CHECK(fabs(vb - 5e-171) < 1e-15 * 5e-171 && vh == vb && ms3[1] > 0 && isnan(mk_pchip_interior(NAN, 1, 1, 1)) && isnan(mk_pchip_interior(1, NAN, 1, 1)),
              "I-B: PCHIP of y = 1e-170 x: interior slope %.3g (not 0), value at 0.5 %.17g (5e-171); a NaN secant still gives a NaN slope", ms3[1], vb);
        double one3[3] = {0, 1, 2};
        t_refuse = 1, t_asked = 0;
        int rs = mk_spline_slopes(one3, one3, 1 << 30, one3);
        t_refuse = 0;
        CHECK(rs == -1 && t_asked == sizeof(double) * 4 * ((size_t)1 << 30),
              "I-D: spline slopes at n = 2^30 ask for %zu bytes (32 GiB), not a wrapped int product", t_asked);
#ifndef T_OLD
        int nl = (1 << 30) + 4;
        size_t bytes = sizeof(double) * (size_t)nl;
        double *xl = mmap(0, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
        if (xl == MAP_FAILED) CHECK(0, "I-C: mmap of an 8 GiB table (MAP_NORESERVE) failed: the row did not run");
        else {
            xl[nl - 1] = nl - 1, xl[nl - 2] = nl - 2; /* the rest read 0: the bisection climbs to the top */
            double tl;
            int il = mk_locate(xl, nl, nl - 1.5, &tl);
            CHECK(il == nl - 2 && tl == 0.5, "I-C: locate in a table of 2^30 + 4 (lo + hi past INT_MAX): interval %d (%d), t %g", il, nl - 2, tl);
            munmap(xl, bytes);
        }
#endif
    }
    { /* 1.6.0, the review's I1: the straight-line parts are mk_lerp (1.5.6: inf halfway, 0 at the last node) */
        const double x2[3] = {0, 1, 2}, yb[3] = {-1e308, 1e308, -1e308}, m0[3] = {0, 0, 0}, ye[2] = {1, 1e-17};
        const double a = mk_linear_eval(x2, yb, 2, 0.5), b = mk_linear_eval(x2, yb, 3, 1.5), c = mk_hermite_eval(x2, yb, m0, 2, 0.5), d = mk_pchip_eval(x2, yb, 2, 0.5);
        const double v = mk_v3_lerp(mk_v3_make(-1e308, 0, 0), mk_v3_make(1e308, 0, 0), 0.5).x;
        const double e = mk_linear_eval(x2, ye, 2, 1), f = mk_v3_lerp(mk_v3_make(1, 0, 0), mk_v3_make(1e-17, 0, 0), 1).x;
        CHECK(a == 0 && b == 0 && c == 0 && d == 0 && v == 0 && e == 1e-17 && f == 1e-17,
              "I1 (review): halfway between -1e308 and 1e308, linear, linear n 3, Hermite and PCHIP n 2, v3_lerp: %g %g %g %g %g (0); at the end of {1, 1e-17}: %g, %g", a, b, c, d,
              v, e, f);
    }
    return wrong;
}
