/* mathkit: summary statistics. Median by selection (average of the two middle values for an even count); quantiles by
   linear interpolation between order statistics (Hyndman and Fan type 7: R's and numpy's default); Welford's running mean (shifted)
   variance (n - 1) and range. The array functions reorder their input and return NaN for an empty one; NaNs in the input
   are not allowed. Header only, C23. */
#ifndef MK_STATS_H
#define MK_STATS_H
#include "mk_vec.h"
#include <math.h>
#include <stdlib.h>

/* the k-th smallest (0-based) of a[n], in O(n) on average; a is partly reordered (Hoare's selection) */
static inline double mk_select(double *a, long n, long k)
{
    if (n <= 0 || k < 0 || k >= n) return NAN;
    long lo = 0, hi = n - 1;
    while (lo < hi) {
        double piv = a[lo + (hi - lo) / 2];
        long i = lo, j = hi;
        while (i <= j) {
            while (a[i] < piv) i++;
            while (a[j] > piv) j--;
            if (i <= j) { double t = a[i]; a[i] = a[j], a[j] = t; i++, j--; }
        }
        if (k <= j) hi = j;
        else if (k >= i) lo = i;
        else break;
    }
    return a[k];
}
static inline double mk_median(double *a, long n)
{
    if (n <= 0) return NAN;
    double hi = mk_select(a, n, n / 2);
#if defined(MK_MUTATE) && MK_MUTATE == 61
    return hi;
#endif
    if (n % 2) return hi;
    double lo = a[0]; /* the largest of the lower half, which selection left below n/2 */
    for (long i = 1; i < n / 2; i++) if (a[i] > lo) lo = a[i];
    const double m = 0.5 * (lo + hi);
    return isfinite(m) ? m : 0.5 * lo + 0.5 * hi; /* 1.6.1: halves first where the sum overflows (mk_lerp's idea; 1.6.0 gave inf) */
}

static inline int mk_cmp_double(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a > b) - (a < b);
}
/* q in [0, 1] of sorted s[n], type 7 */
static inline double mk_quantile_sorted(const double *s, long n, double q)
{
    if (n <= 0 || !(q >= 0 && q <= 1)) return NAN;
#if defined(MK_MUTATE) && MK_MUTATE == 62
    double h = n * q;
#else
    double h = (n - 1) * q;
#endif
    long i = (long)floor(h);
    if (i >= n - 1) return s[n - 1];
    double f = h - i;
    if (f == 0 || s[i] == s[i + 1]) return s[i]; /* an order statistic, or equal neighbours: exact, infinities included */
    if (!isinf(s[i]) != !isinf(s[i + 1])) return isinf(s[i]) ? s[i] : s[i + 1]; /* one infinite neighbour (isinf may return -1) */
    return mk_lerp(s[i], s[i + 1], f); /* 1.6.0: no overflow where s[i + 1] - s[i] would (the review's S1) */
}
/* sorts a */
static inline double mk_quantile(double *a, long n, double q)
{
    if (n <= 0) return NAN;
    qsort(a, (size_t)n, sizeof *a, mk_cmp_double);
    return mk_quantile_sorted(a, n, q);
}

/* Welford's running mean and variance, on data shifted by the first value (Chan, Golub and LeVeque 1983): the running
   mean of x - x0 stays near the spread, not near the mean, so a large mean costs no precision. mean is the true mean. An
   infinite value makes the mean NaN from the next value on (at once if it comes first; a plain mean reads inf), and the variance NaN. */
typedef struct { long n; double mean, m2, min, max, shift, smean; } mk_welford;
static inline void mk_welford_init(mk_welford *w) { *w = (mk_welford){0, 0, 0, INFINITY, -INFINITY, 0, 0}; }
static inline void mk_welford_add(mk_welford *w, double x)
{
#if !(defined(MK_MUTATE) && MK_MUTATE == 64)
    if (w->n == 0) w->shift = x;
#endif
    w->n++;
    double y = x - w->shift, d = y - w->smean;
    w->smean += d / w->n;
    w->m2 += d * (y - w->smean);
    w->mean = w->shift + w->smean;
    if (x < w->min) w->min = x;
    if (x > w->max) w->max = x;
}
/* a gets the statistics of a and b together (Chan's pairwise update, as cudakit's ck_welford_f64 merges threads) */
static inline void mk_welford_merge(mk_welford *a, const mk_welford *b)
{
    if (b->n == 0) return;
    if (a->n == 0) { *a = *b; return; }
    double n = (double)a->n + (double)b->n, d = ((b->shift - a->shift) + b->smean) - a->smean;
    a->smean += d * ((double)b->n / n);
#if defined(MK_MUTATE) && MK_MUTATE == 65
    a->m2 += b->m2;
#else
    a->m2 += b->m2 + d * d * ((double)a->n * (double)b->n / n);
#endif
    a->n += b->n;
    a->mean = a->shift + a->smean;
    a->min = fmin(a->min, b->min), a->max = fmax(a->max, b->max);
}
static inline double mk_welford_var(const mk_welford *w)
{
#if defined(MK_MUTATE) && MK_MUTATE == 63
    return w->n > 0 ? w->m2 / w->n : NAN;
#else
    return w->n > 1 ? w->m2 / (w->n - 1) : NAN;
#endif
}
static inline double mk_welford_sd(const mk_welford *w) { return sqrt(mk_welford_var(w)); }
#endif
