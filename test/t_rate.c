/* mk_rate.h: each law exactly 1 (WLF: 0) at its reference and against its closed form in long double within an error
   bound taken from its own roundings; a Q10 over 10 K, a doubling over its width, the inverse with T and Tref swapped; the
   activation energy that doubles a rate from 298.15 to 308.15 K; the equivalent time of a history of constant pieces, and
   back; the WLF pair moved to another reference, its hold before the pole; Walther through its two points and at the
   temperature between whose log10 is their mean; every refusal. */
#include "../include/mk_rate.h"
#include <stdbool.h>
#include "t.h"

#define U 0x1p-53L /* half an ulp of 1 */
static long double relerr(double got, long double want) { return fabsl(((long double)got - want) / want); }
static double span(double lo, double hi) { return lo + (hi - lo) * (0.5 + 0.5 * t_rand()); }

static void q_rows(void)
{
    int one = 0, per = 0, n = 0;
    for (int i = 0; i < 2000; i++, n++) {
        const double T = span(200, 1500), q = span(0.05, 20), dT = span(-50, 50);
        one += mk_rate_q(T, T, q, dT) == 1;
        const double Tr = floor(T); /* a whole number and a whole width: T - Tref and the quotient are exact */
        per += mk_rate_q(Tr + 10, Tr, q, 10) == q;
    }
    CHECK(one == n, "mk_rate_q: exactly 1 at T = Tref (%d of %d)", one, n);
    CHECK(per == n, "mk_rate_q: a Q10 over 10 K gives q exactly (%d of %d)", per, n);
    CHECK(mk_rate_q(333.15, 303.15, 2, 30) == 2 && mk_rate_q(1234, 1219, 10, 15) == 10, "mk_rate_q: 2 per 30 K over 30 K is 2, 10 per 15 K over 15 K is 10");
    long double worst = 0, worstinv = 0;
    for (int i = 0; i < 200000; i++) {
        const double T = span(-60, 160), Tr = span(-20, 60), q = span(0.1, 10), dT = span(1, 40);
        const long double x = ((long double)T - Tr) / dT, want = powl(q, x), e = fabsl(x * logl(q));
        const double f = mk_rate_q(T, Tr, q, dT);
        worst = fmaxl(worst, relerr(f, want) / (U * (2 + 2 * e)));        /* T - Tref and the quotient round, then pow */
        worstinv = fmaxl(worstinv, fabsl((long double)f * mk_rate_q(Tr, T, q, dT) - 1) / (4 * U));
    }
    CHECK(worst <= 1, "mk_rate_q: 200000 cases against powl, the largest error %.2Lf of its bound 2^-53 (2 + 2 |ln f|)", worst);
    CHECK(worstinv <= 1.01L, "mk_rate_q: T and Tref swapped is the inverse, the product within %.2Lf of 2^-51 of 1", worstinv);
}

static void arrhenius_rows(void)
{
    int one = 0, n = 0;
    for (int i = 0; i < 2000; i++, n++) {
        const double T = span(1, 3000);
        one += mk_rate_arrhenius(T, T, span(-1e4, 1e5)) == 1;
    }
    CHECK(one == n, "mk_rate_arrhenius: exactly 1 at T = Tref (%d of %d)", one, n);
    long double worst = 0, worsta = 0;
    for (int i = 0; i < 200000; i++) {
        const double T = span(200, 1500), Tr = span(200, 800), th = span(-6000, 80000);
        const long double y = -(long double)th * (1 / (long double)T - 1 / (long double)Tr), want = expl(y);
        const double f = mk_rate_arrhenius(T, Tr, th);
        if (want > 1e-300L && want < 1e300L)
            worst = fmaxl(worst, relerr(f, want) / (U * (2 + 3 * fabsl(y) + fabsl(th) * (1 / (long double)T + 1 / (long double)Tr))));
        const long double a = expl(-(long double)th / T);
        if (a > 1e-300L && a < 1e300L) worsta = fmaxl(worsta, relerr(mk_rate_arrhenius(T, INFINITY, th), a) / (U * (2 + 3 * fabsl((long double)th / T))));
    }
    CHECK(worst <= 1, "mk_rate_arrhenius: 200000 cases against expl, the largest error %.2Lf of its bound", worst);
    CHECK(worsta <= 1, "mk_rate_arrhenius: Tref INFINITY is exp(-theta / T), the largest error %.2Lf of its bound", worsta);
    /* the energy that doubles a rate from 298.15 to 308.15 K, 52.9 kJ/mol */
    const double E = MK_RATE_GAS_R * log(2) * 298.15 * 308.15 / (308.15 - 298.15), f2 = mk_rate_arrhenius(308.15, 298.15, E / MK_RATE_GAS_R);
    CHECK(fabs(E - 52.948e3) < 5 && fabs(f2 - 2) < 1e-13, "mk_rate_arrhenius: %.3f kJ/mol doubles a rate from 298.15 to 308.15 K (2 %+.1e)", E * 1e-3, f2 - 2);
    CHECK(fabs(MK_RATE_GAS_R / (6.02214076e23 * 1.380649e-23) - 1) < 1e-15, "MK_RATE_GAS_R is Avogadro's number times Boltzmann's");
    CHECK(mk_rate_arrhenius(1e-300, 300, 9000) == 0 && mk_rate_arrhenius(INFINITY, INFINITY, 9000) == 1, "mk_rate_arrhenius: 0 toward 0 K, 1 with both infinite");
    /* a history of two constant pieces, 100 s at 300 K and 50 s at 350 K, at 320 K: the sum over 0.25 s steps against the
       two pieces' closed form, and the time at 350 K that does what the history did */
    const double th = 9000, Tr = 320;
    double teq = 0;
    for (int i = 0; i < 600; i++) teq += 0.25 * mk_rate_arrhenius(i < 400 ? 300 : 350, Tr, th);
    const long double f300 = expl(-9000.0L * (1 / 300.0L - 1 / 320.0L)), f350 = expl(-9000.0L * (1 / 350.0L - 1 / 320.0L));
    const long double want = 100 * f300 + 50 * f350;
    CHECK(relerr(teq, want) < 1e-13L, "equivalent time of 100 s at 300 K and 50 s at 350 K: %.6f s at 320 K (closed form %.6Lf)", teq, want);
    const double hot = teq / mk_rate_arrhenius(350, Tr, th), back = hot * mk_rate_arrhenius(350, Tr, th);
    CHECK(relerr(hot, 50 + 100 * expl(-9000.0L * (1 / 300.0L - 1 / 350.0L))) < 1e-13L && relerr(back, teq) < 1e-15L, "and back: %.6f s at 350 K does the same", hot);
    CHECK(relerr(100 * mk_rate_arrhenius(300, 350, th), 100 * expl(-9000.0L * (1 / 300.0L - 1 / 350.0L))) < 1e-13L, "100 s at 300 K is %.6f s at 350 K directly", 100 * mk_rate_arrhenius(300, 350, th));
}

static void wlf_rows(void)
{
    int zero = 0, n = 0;
    for (int i = 0; i < 2000; i++, n++) {
        const double Tr = span(150, 400);
        zero += mk_rate_wlf_log10(Tr, Tr, span(1, 30), span(20, 150)) == 0;
    }
    CHECK(zero == n, "mk_rate_wlf_log10: exactly 0 at Tref (%d of %d)", zero, n);
    long double worst = 0;
    for (int i = 0; i < 200000; i++) {
        const double Tr = span(150, 400), C1 = span(1, 30), C2 = span(20, 150), T = Tr + span(-0.95, 4) * C2;
        const long double d = (long double)T - Tr, want = -C1 * d / (C2 + d);
        /* d is off by U |d| from the subtraction; the law's slope turns that into U |d| C1 C2 / (C2 + d)^2 */
        const long double bound = U * (4 * fabsl(want) + (fabsl(d) + C2 + d) * C1 * C2 / ((C2 + d) * (C2 + d)));
        worst = fmaxl(worst, fabsl(mk_rate_wlf_log10(T, Tr, C1, C2) - want) / bound);
    }
    CHECK(worst <= 1, "mk_rate_wlf_log10: 200000 cases against long double, the largest error %.2Lf of its bound", worst);
    /* the universal pair at the glass, 17.44 and 51.6 K, carried 50 K up: C1 C2 / (C2 + 50) and C2 + 50, the pair 8.86, 101.6 */
    const double C1 = 17.44 * 51.6 / 101.6, C2 = 101.6;
    double off = 0;
    for (double T = 230; T < 480; T += 0.25)
        off = fmax(off, fabs(mk_rate_wlf_log10(T, 230, 17.44, 51.6) - mk_rate_wlf_log10(280, 230, 17.44, 51.6) - mk_rate_wlf_log10(T, 280, C1, C2)));
    CHECK(fabs(C1 - 8.86) < 0.005 && off < 1e-13, "mk_rate_wlf_log10: 17.44, 51.6 K at the glass is %.2f, 101.6 K 50 K above it, the laws within %.1e", C1, off);
    CHECK(fabs(mk_rate_wlf_log10(230 + 51.6, 230, 17.44, 51.6) + 17.44 / 2) < 1e-13, "mk_rate_wlf_log10: -C1 / 2 at C2 above the reference");
    bool held = true, down = true;
    double last = INFINITY;
    for (double T = -100; T < 400; T += 0.5) {
        const double v = mk_rate_wlf_log10(T, 230, 17.44, 51.6);
        held &= isfinite(v) && (T > 230 - 0.99 * 51.6 || fabs(v - 99 * 17.44) < 1e-10);
        down &= v <= last, last = v;
    }
    CHECK(held && down, "mk_rate_wlf_log10: 99 C1 at and below Tref - 0.99 C2 (the pole included), finite and never rising with T");
    CHECK(mk_rate_wlf_log10(100, 333.15, 0, 100) == 0 && mk_rate_wlf_log10(400, 333.15, 0, 100) == 0, "mk_rate_wlf_log10: C1 = 0 is no shift at any T");
    /* a time: 10 s at 330 K, 5 s at 310 K, in seconds at 320 K */
    const double teq = 10 / pow(10, mk_rate_wlf_log10(330, 320, 8.86, 101.6)) + 5 / pow(10, mk_rate_wlf_log10(310, 320, 8.86, 101.6));
    const long double want = 10 * powl(10, 8.86L * 10 / 111.6L) + 5 * powl(10, -8.86L * 10 / 91.6L);
    CHECK(relerr(teq, want) < 1e-13L, "equivalent time by a_T: %.6f s at 320 K (closed form %.6Lf)", teq, want);
}

static void walther_rows(void)
{
    /* an ISO VG 32 oil: 32 mm^2/s at 40 C, 5.4 at 100 C */
    double A = 0, B = 0;
    const int r = mk_rate_walther_fit(313.15, 32, 373.15, 5.4, &A, &B);
    CHECK(r == 0 && relerr(mk_rate_walther(313.15, A, B), 32) < 1e-13L && relerr(mk_rate_walther(373.15, A, B), 5.4) < 1e-13L, "mk_rate_walther: through its two points (A %.6f, B %.6f)", A, B);
    const long double z = (log10l(log10l(32.7L)) + log10l(log10l(6.1L))) / 2, tm = sqrtl(313.15L * 373.15L);
    const double mid = mk_rate_walther((double)tm, A, B);
    CHECK(relerr(mid, powl(10, powl(10, z)) - 0.7L) < 1e-13L, "mk_rate_walther: %.4f mm^2/s at %.3f K, where log10 T is the points' mean", mid, (double)tm);
    long double worst = 0;
    for (int i = 0; i < 100000; i++) {
        const double a = span(6, 11), b = span(2.2, 4.2), T = span(230, 460);
        const long double w = a - b * log10l(T), p = powl(10, w), v = powl(10, p), want = v - 0.7L;
        if (!(w > -1.5L && w < 0.9L) || !(want > 0.3L)) continue;
        const long double l10 = 2.302585092994045684L, ew = U * (fabsl(w) + 3 * fabsl(b * log10l(T))); /* the exponent's error */
        const long double bound = v / want * (p * l10 * (l10 * ew + 2 * U) + 2 * U) + 2 * U;
        worst = fmaxl(worst, relerr(mk_rate_walther(T, a, b), want) / bound);
    }
    CHECK(worst > 0 && worst <= 1, "mk_rate_walther: random lines against powl, the largest error %.2Lf of its bound", worst);
    CHECK(mk_rate_walther_fit(313.15, 7, 373.15, 7, &A, &B) == 0 && B == 0 && relerr(mk_rate_walther(250, A, B), 7) < 1e-14L, "mk_rate_walther: two equal viscosities, the same at any T");
    CHECK(mk_rate_walther_fit(313.15, 32, 373.15, 5.4, &A, &B) == 0 && mk_rate_walther(60, A, B) == INFINITY, "mk_rate_walther: INFINITY where it overflows");
}

static void refusals(void)
{
    const double N = NAN;
    const double q[][4] = {{N, 300, 2, 10}, {300, N, 2, 10}, {300, 300, N, 10}, {300, 300, 2, N}, {N, 300, 1, 10}, {300, N, 1, 10},
        {310, 300, 0, 10}, {310, 300, -2, 10}, {310, 300, 2, 0}, {300, 300, 2, 0}};
    int nq = 0, na = 0, nw = 0, nv = 0, nf = 0;
    for (size_t i = 0; i < sizeof q / sizeof *q; i++) nq += isnan(mk_rate_q(q[i][0], q[i][1], q[i][2], q[i][3])) != 0;
    CHECK(nq == 10, "mk_rate_q: NaN for a NaN (q = 1 too), q not above 0, dT 0 (%d of 10)", nq);
    const double a[][3] = {{N, 300, 9000}, {300, N, 9000}, {300, 310, N}, {0, 300, 9000}, {-0.0, 300, 9000}, {-5, 300, 9000}, {300, 0, 9000},
        {300, -1, 9000}, {-INFINITY, 300, 9000}, {N, 300, 0}};
    for (size_t i = 0; i < sizeof a / sizeof *a; i++) na += isnan(mk_rate_arrhenius(a[i][0], a[i][1], a[i][2])) != 0;
    CHECK(na == 10, "mk_rate_arrhenius: NaN for a NaN, T or Tref not above 0 (%d of 10)", na);
    const double w[][4] = {{N, 300, 8.86, 101.6}, {300, N, 8.86, 101.6}, {300, 310, N, 101.6}, {300, 310, 8.86, N}, {300, 310, 8.86, 0}, {300, 310, 8.86, -101.6},
        {310, 300, 8.86, -101.6}, {INFINITY, 300, 8.86, 101.6}};
    for (size_t i = 0; i < sizeof w / sizeof *w; i++) nw += isnan(mk_rate_wlf_log10(w[i][0], w[i][1], w[i][2], w[i][3])) != 0;
    CHECK(nw == 8, "mk_rate_wlf_log10: NaN for a NaN, C2 not above 0, T INFINITY (%d of 8)", nw);
    const double v[][3] = {{N, 8, 3}, {300, N, 3}, {300, 8, N}, {0, 8, 3}, {-300, 8, 3}};
    for (size_t i = 0; i < sizeof v / sizeof *v; i++) nv += isnan(mk_rate_walther(v[i][0], v[i][1], v[i][2])) != 0;
    CHECK(nv == 5, "mk_rate_walther: NaN for a NaN, T not above 0 (%d of 5)", nv);
    const double f[][4] = {{313, 0.3, 373, 5}, {313, 0.1, 373, 5}, {313, 32, 373, -1}, {313, 32, 313, 5}, {0, 32, 373, 5}, {313, 32, -373, 5},
        {N, 32, 373, 5}, {313, N, 373, 5}, {313, 32, N, 5}, {313, 32, 373, N}, {313, INFINITY, 373, 5}, {300, 5, 0, 8},
        {300, 5, -0.0, 8}, {300, 5, INFINITY, 8}, {INFINITY, 5, 300, 8}};
    for (size_t i = 0; i < sizeof f / sizeof *f; i++) {
        double A = 11, B = 22;
        nf += mk_rate_walther_fit(f[i][0], f[i][1], f[i][2], f[i][3], &A, &B) == -1 && A == 11 && B == 22;
    }
    CHECK(nf == 15, "mk_rate_walther_fit: -1 with A and B untouched for a point with no log10 log10, equal temperatures, a NaN, T 0 or INFINITY (%d of 15)", nf);
}

int main(void)
{
    q_rows(), arrhenius_rows(), wlf_rows(), walther_rows(), refusals();
    return wrong;
}
