#include "../include/mk_rng.h"
#ifdef MK_RIBBON_RNG_H   /* another project's RNG to compare with bit for bit: -DMK_RIBBON_RNG_H='"path/ribbon_rng.h"' */
#include MK_RIBBON_RNG_H
#endif
#include "t.h"

int main(void)
{
    /* the authors' reference: xoshiro256** from state {1, 2, 3, 4} */
    mk_rng g = {{1, 2, 3, 4}, 0, 0};
    uint64_t want[4] = {11520ULL, 0ULL, 1509978240ULL, 1215971899390074240ULL}, got[4];
    for (int i = 0; i < 4; i++) got[i] = mk_rng_next(&g);
    CHECK(got[0] == want[0] && got[1] == want[1] && got[2] == want[2] && got[3] == want[3],
          "xoshiro256** from {1,2,3,4}: 11520, 0, 1509978240, 1215971899390074240 (got %llu %llu %llu %llu)",
          (unsigned long long)got[0], (unsigned long long)got[1], (unsigned long long)got[2], (unsigned long long)got[3]);
    mk_rng a;
#ifdef MK_RIBBON_RNG_H
    RbRng b;
    mk_rng_seed(&a, 12345, 7), rb_seed(&b, 12345, 7);
    int same = 1;
    for (int i = 0; i < 100000; i++) same &= mk_rng_next(&a) == rb_next(&b);
    CHECK(same, "seeded (12345, stream 7): 100000 draws equal ribbon_rng.h's bit for bit");
#endif

    mk_rng_seed(&a, 1, 0);
    enum { N = 1000000 };
    double s1 = 0, s2 = 0, s4 = 0, u = 0, umin = 1;
    for (int i = 0; i < N; i++) {
        double z = mk_rng_gauss(&a);
        s1 += z, s2 += z * z, s4 += z * z * z * z;
        double v = mk_rng_u01(&a);
        u += v, umin = fmin(umin, v);
    }
    s1 /= N, s2 /= N, s4 /= N, u /= N;
    CHECK(fabs(s1) < 0.005 && fabs(s2 - 1) < 0.01 && fabs(s4 - 3) < 0.05, "Gaussian, 1e6 draws: mean %.4f, variance %.4f, 4th moment %.3f (0, 1, 3)", s1, s2, s4);
    CHECK(fabs(u - 0.5) < 0.002 && umin > 0, "uniform (0,1]: mean %.4f, smallest %.3g, never 0", u, umin);
    long cnt[3] = {0};
    for (int i = 0; i < 300000; i++) cnt[mk_rng_below(&a, 3)]++;
    double chi = 0;
    for (int k = 0; k < 3; k++) chi += (cnt[k] - 100000.0) * (cnt[k] - 100000.0) / 100000.0;
    CHECK(chi < 13.8, "below(3), 300000 draws: chi-squared %.2f on 2 degrees (limit 13.8, p 0.001)", chi);

    mk_sum s = {0, 0};
    double xs[4] = {1.0, 1e100, 1.0, -1e100};
    for (int i = 0; i < 4; i++) mk_sum_add(&s, xs[i]);
    CHECK(mk_sum_value(&s) == 2.0, "compensated sum of 1, 1e100, 1, -1e100 is 2 (got %g; plain sum 0)", mk_sum_value(&s));
    mk_sum t = {0, 0};
    for (int i = 0; i < 10000000; i++) mk_sum_add(&t, 0.1);
    CHECK(fabs(mk_sum_value(&t) - 1e6) < 1e-9, "1e7 additions of 0.1: %.12f (error %.2g)", mk_sum_value(&t), fabs(mk_sum_value(&t) - 1e6));
    double inf = INFINITY, sa[3] = {1, inf, 2}, sb[3] = {1, -inf, 5}, sc[3] = {inf, 1, -inf};
    mk_sum ua = {0, 0}, ub = {0, 0}, uc = {0, 0};
    for (int i = 0; i < 3; i++) mk_sum_add(&ua, sa[i]), mk_sum_add(&ub, sb[i]), mk_sum_add(&uc, sc[i]);
    double va = mk_sum_value(&ua), vb = mk_sum_value(&ub), vc = mk_sum_value(&uc);
    CHECK(va == inf && vb == -inf && isnan(vc), "non-finite addends give a plain sum's result: 1 + inf + 2 = %g, 1 - inf + 5 = %g, inf + 1 - inf = %g", va, vb, vc);
    { /* 1.2.2: the verification work Result 84, R1 and R2 */
        CHECK(mk_rotl64(0x8000000000000001ULL, 0) == 0x8000000000000001ULL && mk_rotl64(1, 64) == 1, "R1: rotl by 0 and by 64 is the identity");
#ifndef T_OLD
        mk_rng g0;
        mk_rng_seed(&g0, 7, 0);
        CHECK(mk_rng_below(&g0, 0) == 0, "R2: below(0) returns 0 instead of dividing by zero");
#endif
    }
    return wrong;
}
