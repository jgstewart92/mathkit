/* mathkit: summary statistics. Median by selection (the mean of the two middle values for an even count); quantiles by
   linear interpolation between order statistics (Hyndman and Fan 1996 type 7, R's and numpy's default); Welford's
   running mean, variance (divisor n - 1) and range on data shifted by the first value, with a pairwise merge. The array
   functions reorder their input and return NaN for an empty one; NaNs in the input are not allowed.
   Header only, C23. */
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
    return hi; /* planted fault: the upper middle value for an even count */
#endif
    if (n % 2) return hi;
    double lo = a[0]; /* the largest of the lower half, which selection left below n/2 */
    for (long i = 1; i < n / 2; i++) if (a[i] > lo) lo = a[i];
    const double m = 0.5 * (lo + hi);
    return isfinite(m) ? m : 0.5 * lo + 0.5 * hi; /* halves first where the sum overflows */
}

static inline int mk_cmp_double(const void *x, const void *y)
{
    double a = *(const double *)x, b = *(const double *)y;
    return (a > b) - (a < b);
}
/* the q-quantile of sorted s[n], type 7; NaN for n < 1 or q outside [0, 1] */
static inline double mk_quantile_sorted(const double *s, long n, double q)
{
    if (n <= 0 || !(q >= 0 && q <= 1)) return NAN;
#if defined(MK_MUTATE) && MK_MUTATE == 62
    double h = n * q; /* planted fault: n q for (n - 1) q */
#else
    double h = (n - 1) * q;
#endif
    long i = (long)floor(h);
    if (i >= n - 1) return s[n - 1];
    double f = h - i;
    /* an order statistic, or equal neighbours: exact, infinities included */
    if (f == 0 || s[i] == s[i + 1]) return s[i];
    /* one infinite neighbour (isinf may return -1) */
    if (!isinf(s[i]) != !isinf(s[i + 1])) return isinf(s[i]) ? s[i] : s[i + 1];
    return mk_lerp(s[i], s[i + 1], f); /* no overflow where s[i + 1] - s[i] would */
}
/* sorts a in place, then mk_quantile_sorted */
static inline double mk_quantile(double *a, long n, double q)
{
    if (n <= 0) return NAN;
    qsort(a, (size_t)n, sizeof *a, mk_cmp_double);
    return mk_quantile_sorted(a, n, q);
}

/* Welford's running mean and variance (Welford 1962), on data shifted by the first value (Chan, Golub and LeVeque
   1983): the running mean of x - x0 stays near the spread, not near the mean, so a large mean costs no precision. mean
   is the true mean. An infinite value makes the mean NaN from the next value on (at once if it comes first; a plain
   mean reads inf), and the variance NaN. */
typedef struct { long n; double mean, m2, min, max, shift, smean; } mk_welford;
static inline void mk_welford_init(mk_welford *w) { *w = (mk_welford){0, 0, 0, INFINITY, -INFINITY, 0, 0}; }
static inline void mk_welford_add(mk_welford *w, double x)
{
#if !(defined(MK_MUTATE) && MK_MUTATE == 64) /* planted fault 64: no shift */
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
/* a gets the statistics of a and b together (Chan, Golub and LeVeque's pairwise update) */
static inline void mk_welford_merge(mk_welford *a, const mk_welford *b)
{
    if (b->n == 0) return;
    if (a->n == 0) { *a = *b; return; }
    double n = (double)a->n + (double)b->n, d = ((b->shift - a->shift) + b->smean) - a->smean;
    a->smean += d * ((double)b->n / n);
#if defined(MK_MUTATE) && MK_MUTATE == 65
    a->m2 += b->m2; /* planted fault: the term for the difference of means left out */
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
    return w->n > 0 ? w->m2 / w->n : NAN; /* planted fault: divisor n */
#else
    return w->n > 1 ? w->m2 / (w->n - 1) : NAN;
#endif
}
static inline double mk_welford_sd(const mk_welford *w) { return sqrt(mk_welford_var(w)); }
#endif
