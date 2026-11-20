/* mathkit: 1-D interpolation on increasing nodes x[n]. Monotone cubic (PCHIP: Fritsch-Butland interior slopes, shape-
   preserving three-point ends, as scipy's PchipInterpolator) computed exactly as aeromap/aeromap.c does, so its results
   are bit-identical to that proved code, except where d0 * d1 underflows (aeromap.c gives a 0 slope there, 1.0.5 does
   not); natural cubic spline (C2); linear. n 1 gives y0 and n < 1 NaN (mk_locate itself needs n >= 2). Queries outside
   [x0, xn-1] are clamped to the end. The straight-line parts (mk_linear_eval, and n 2) are mk_vec.h's mk_lerp (1.6.0): the
   last node's value exactly, and no overflow between values of opposite sign near DBL_MAX. Header only, C23. */
#ifndef MK_INTERP_H
#define MK_INTERP_H
#include "mk_vec.h"
#include <math.h>
#include <stdlib.h>

static inline int mk_sign(double x) { return (x > 0) - (x < 0); }
static inline double mk_pchip_interior(double d0, double d1, double h0, double h1)
{
    /* signs, since d0 * d1 underflows for tiny secants; a NaN secant still gives NaN, as aeromap.c */
    if (mk_sign(d0) * mk_sign(d1) <= 0) return isnan(d0) || isnan(d1) ? NAN : 0;
    double w1 = 2 * h1 + h0, w2 = h1 + 2 * h0;
#if defined(MK_MUTATE) && MK_MUTATE == 71
    (void)w1, (void)w2;
    return 0.5 * (d0 + d1);
#else
    return (w1 + w2) / (w1 / d0 + w2 / d1);
#endif
}
static inline double mk_pchip_end(double d0, double d1, double h0, double h1)
{
    double m = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
    if (mk_sign(m) != mk_sign(d0)) return 0;
    if (mk_sign(d0) != mk_sign(d1) && fabs(m) > 3 * fabs(d0)) return 3 * d0;
    return m;
}
/* the interval holding q (clamped) and its fraction t */
static inline int mk_locate(const double *x, int n, double q, double *t)
{
    if (q < x[0]) q = x[0];
    if (q > x[n - 1]) q = x[n - 1];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) { int mid = lo + (hi - lo) / 2; if (x[mid] <= q) lo = mid; else hi = mid; }
    *t = (q - x[lo]) / (x[lo + 1] - x[lo]);
    return lo;
}
/* PCHIP slopes m[n] (n >= 2) */
static inline void mk_pchip_slopes(const double *x, const double *y, int n, double *m)
{
    if (n < 2) { if (n == 1) m[0] = 0; return; }
    if (n == 2) { m[0] = m[1] = (y[1] - y[0]) / (x[1] - x[0]); return; }
    for (int i = 1; i < n - 1; i++)
        m[i] = mk_pchip_interior((y[i] - y[i - 1]) / (x[i] - x[i - 1]), (y[i + 1] - y[i]) / (x[i + 1] - x[i]), x[i] - x[i - 1], x[i + 1] - x[i]);
    m[0] = mk_pchip_end((y[1] - y[0]) / (x[1] - x[0]), (y[2] - y[1]) / (x[2] - x[1]), x[1] - x[0], x[2] - x[1]);
    m[n - 1] = mk_pchip_end((y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]), (y[n - 2] - y[n - 3]) / (x[n - 2] - x[n - 3]),
                            x[n - 1] - x[n - 2], x[n - 2] - x[n - 3]);
}
/* cubic Hermite through (x, y) with slopes m; with n == 2 it is the straight line, as aeromap */
static inline double mk_hermite_eval(const double *x, const double *y, const double *m, int n, double q)
{
    if (n < 2) return n == 1 ? y[0] : NAN;
    double t;
    int i = mk_locate(x, n, q, &t);
    if (n == 2) return mk_lerp(y[i], y[i + 1], t);
    double h = x[i + 1] - x[i], t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * y[i] + (t3 - 2 * t2 + t) * h * m[i] + (-2 * t3 + 3 * t2) * y[i + 1] + (t3 - t2) * h * m[i + 1];
}
/* PCHIP at q without a slope table: only the two slopes it needs */
static inline double mk_pchip_eval(const double *x, const double *y, int n, double q)
{
    if (n < 2) return n == 1 ? y[0] : NAN;
    double t;
    int i = mk_locate(x, n, q, &t);
    if (n == 2) return mk_lerp(y[i], y[i + 1], t);
    double h = x[i + 1] - x[i], d = (y[i + 1] - y[i]) / h, m0, m1;
    if (i > 0) m0 = mk_pchip_interior((y[i] - y[i - 1]) / (x[i] - x[i - 1]), d, x[i] - x[i - 1], h);
    else m0 = mk_pchip_end(d, (y[2] - y[1]) / (x[2] - x[1]), h, x[2] - x[1]);
    if (i + 2 < n) m1 = mk_pchip_interior(d, (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]), h, x[i + 2] - x[i + 1]);
    else m1 = mk_pchip_end(d, (y[i] - y[i - 1]) / (x[i] - x[i - 1]), h, x[i] - x[i - 1]);
    double t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * y[i] + (t3 - 2 * t2 + t) * h * m0 + (-2 * t3 + 3 * t2) * y[i + 1] + (t3 - t2) * h * m1;
}
/* natural cubic spline slopes m[n] (second derivative 0 at both ends), for mk_hermite_eval. 0, or -1 for no memory. */
static inline int mk_spline_slopes(const double *x, const double *y, int n, double *m)
{
    if (n < 3) { mk_pchip_slopes(x, y, n, m); return 0; }
    double *w = malloc(sizeof(double) * 4 * (size_t)n); /* 4 * n overflows int for n >= 2^29 */
    if (!w) return -1;
    double *h = w, *M = h + n, *cp = M + n, *dp = cp + n;
    for (int i = 0; i < n - 1; i++) h[i] = x[i + 1] - x[i];
    cp[0] = 0, dp[0] = 0;
    for (int i = 1; i < n - 1; i++) {
        double a = h[i - 1], b = 2 * (h[i - 1] + h[i]), c = h[i];
        double r = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        double den = b - a * cp[i - 1];
        cp[i] = c / den, dp[i] = (r - a * dp[i - 1]) / den;
    }
    M[n - 1] = 0;
    for (int i = n - 2; i >= 1; i--) M[i] = dp[i] - cp[i] * M[i + 1];
#if defined(MK_MUTATE) && MK_MUTATE == 72
    M[0] = M[1];
#else
    M[0] = 0;
#endif
    for (int i = 0; i < n - 1; i++) m[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * M[i] + M[i + 1]) / 6;
    m[n - 1] = (y[n - 1] - y[n - 2]) / h[n - 2] + h[n - 2] * (M[n - 2] + 2 * M[n - 1]) / 6;
    free(w);
    return 0;
}
static inline double mk_linear_eval(const double *x, const double *y, int n, double q)
{
    if (n < 2) return n == 1 ? y[0] : NAN;
    double t;
    int i = mk_locate(x, n, q, &t);
    return mk_lerp(y[i], y[i + 1], t);
}
#endif
