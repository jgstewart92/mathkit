#include "../include/mk_stats.h"
#include "t.h"
#include <string.h>

int main(void)
{
    double a[7] = {9, 1, 8, 2, 7, 3, 5}, b[6] = {9, 1, 8, 2, 7, 3};
    CHECK(mk_median(a, 7) == 5, "median of 9 1 8 2 7 3 5 is 5");
    CHECK(mk_median(b, 6) == 5, "median of 9 1 8 2 7 3 is 5 (the two middle values averaged)");
    int bad = 0;
    for (int n = 1; n <= 60; n++)
        for (int rep = 0; rep < 20; rep++) {
            double x[60], s[60];
            for (int i = 0; i < n; i++) x[i] = s[i] = floor(t_rand() * 5); /* many ties */
            qsort(s, (size_t)n, sizeof *s, mk_cmp_double);
            double want = n % 2 ? s[n / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2]);
            if (mk_median(x, n) != want) bad++;
        }
    CHECK(bad == 0, "median by selection equals the sorted median, 1200 arrays of 1..60 with ties (%d differ)", bad);
    double c[10];
    for (int i = 0; i < 10; i++) c[i] = 10 - i;
    double q25 = mk_quantile(c, 10, 0.25), q90 = mk_quantile(c, 10, 0.9), q0 = mk_quantile(c, 10, 0), q1 = mk_quantile(c, 10, 1);
    CHECK(fabs(q25 - 3.25) < 1e-15 && fabs(q90 - 9.1) < 1e-14 && q0 == 1 && q1 == 10, "quantiles of 1..10 (type 7): 0.25 %.15g, 0.9 %.15g, 0 %g, 1 %g (3.25, 9.1, 1, 10)", q25, q90, q0, q1);
    mk_welford w;
    mk_welford_init(&w);
    double d[8] = {2, 4, 4, 4, 5, 5, 7, 9};
    for (int i = 0; i < 8; i++) mk_welford_add(&w, d[i]);
    CHECK(w.mean == 5 && fabs(mk_welford_var(&w) - 32.0 / 7) < 1e-14 && w.min == 2 && w.max == 9, "Welford on 2 4 4 4 5 5 7 9: mean %g, variance %.15g (32/7), range %g..%g", w.mean, mk_welford_var(&w), w.min, w.max);
    mk_welford_init(&w);
    for (int i = 0; i < 1000; i++) mk_welford_add(&w, 1e9 + (i % 2 ? 1 : -1));
    CHECK(fabs(mk_welford_var(&w) - 1000.0 / 999) < 1e-9, "1e9 +- 1, 1000 values: variance %.12f (1.001001...), no cancellation", mk_welford_var(&w));
    CHECK(isnan(mk_median(a, 0)) && isnan(mk_quantile(a, 0, 0.5)), "an empty array gives NaN");
    for (int pass = 0; pass < 2; pass++) {
        double off = pass ? 1e9 : 1e6;
        long nn = 1000000;
        double *xx = malloc(sizeof(double) * (size_t)nn);
        long double sm = 0, sv = 0;
        for (long i = 0; i < nn; i++) xx[i] = off + 1.7 * t_rand(), sm += xx[i];
        long double lm = sm / nn;
        for (long i = 0; i < nn; i++) sv += (xx[i] - lm) * (xx[i] - lm);
        double rv = (double)(sv / (nn - 1));
        mk_welford one, parts, piece;
        mk_welford_init(&one), mk_welford_init(&parts);
        for (long i = 0; i < nn; i++) mk_welford_add(&one, xx[i]);
        long cut[8] = {0, 3, 1000, 1001, 250000, 600000, 999990, nn};
        for (int c = 7; c > 0; c--) { /* pieces merged in reverse order and uneven sizes */
            mk_welford_init(&piece);
            for (long i = cut[c - 1]; i < cut[c]; i++) mk_welford_add(&piece, xx[i]);
            mk_welford_merge(&parts, &piece);
        }
        double ev = fabs(mk_welford_var(&one) - rv) / rv, em = fabs(one.mean - (double)lm) / (double)lm;
        double pv = fabs(mk_welford_var(&parts) - rv) / rv, pm = fabs(parts.mean - (double)lm) / (double)lm;
        CHECK(ev < 1e-12 && em < 1e-15, "Welford near a mean of %.0e, 1e6 values: variance %.2g and mean %.2g relative to a long-double two-pass", off, ev, em);
        CHECK(pv < 1e-12 && pm < 1e-15 && parts.n == nn && parts.min == one.min && parts.max == one.max, "seven uneven pieces merged (Chan): variance %.2g, mean %.2g relative, count and range equal", pv, pm);
        free(xx);
    }
    { /* 1.2.2: the verification work Result 84, S1 and S2 */
        double a3[3] = {1, 2, INFINITY}, a2[2] = {-INFINITY, 1}, a5[5] = {4, 3, 2, 1, 99};
        CHECK(mk_quantile_sorted(a3, 3, 0.5) == 2, "S1: quantile 0.5 of {1, 2, inf} is 2 (got %g)", mk_quantile_sorted(a3, 3, 0.5));
        CHECK(mk_quantile_sorted(a2, 2, 0) == -INFINITY, "S1: quantile 0 of {-inf, 1} is -inf (got %g)", mk_quantile_sorted(a2, 2, 0));
        CHECK(isnan(mk_select(a5, 4, 4)), "S2: select with k = n is refused (NaN), never a[n]");
    }
    { /* 1.6.0, the review's S1: neighbours whose difference overflows (1.5.6: inf) */
        const double big[2] = {-0x1.fffffffffffffp+1023, 0x1.fffffffffffffp+1023};
        const double q5 = mk_quantile_sorted(big, 2, 0.5), q25 = mk_quantile_sorted(big, 2, 0.25);
        CHECK(q5 == 0 && fabs(q25 + 0.5 * big[1]) <= 1e-15 * big[1], "S1 (review): quantiles of {-DBL_MAX, DBL_MAX}: 0.5 gives %g (0), 0.25 gives %.17g (-DBL_MAX / 2)", q5, q25);
    }
    { /* review 2 #5 (1.6.1): the even median's 0.5 (lo + hi) overflowed on finite data (1.6.0: inf) */
        double a[2] = {1e308, 1.7e308}, b[2] = {-1.7e308, -1e308}, c[4] = {0x1.fffffffffffffp1023, 0x1.fffffffffffffp1023, 1, 0x1.fffffffffffffp1023}, d[2] = {0x1p-1074, 0x1p-1074};
        const double ma = mk_median(a, 2), mb = mk_median(b, 2), mc = mk_median(c, 4), md = mk_median(d, 2);
        CHECK(ma == 1.35e308 && mb == -1.35e308 && mc == 0x1.fffffffffffffp1023 && md == 0x1p-1074, "review 2 #5: medians {1e308, 1.7e308} %g, {-1.7e308, -1e308} %g, {DBL_MAX x 3, 1} %g; two of the least subnormal %g",
              ma, mb, mc, md);
    }
    return wrong;
}
