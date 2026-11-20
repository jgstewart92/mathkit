#include "../include/mk_opt.h"
#include "t.h"

static double rosen(const double *p, void *ctx) { (void)ctx; return 100 * (p[1] - p[0] * p[0]) * (p[1] - p[0] * p[0]) + (1 - p[0]) * (1 - p[0]); }
static double quad4(const double *p, void *ctx)
{
    (void)ctx;
    double s = 0, c[4] = {1, -2, 3, 0.5};
    for (int i = 0; i < 4; i++) s += (i + 1) * (p[i] - c[i]) * (p[i] - c[i]);
    return s + 0.3 * (p[0] - c[0]) * (p[1] - c[1]);
}
static void rosen_r(const double *p, double *r, void *ctx) { (void)ctx; r[0] = 10 * (p[1] - p[0] * p[0]), r[1] = 1 - p[0]; }
static void rosen_j(const double *p, double *J, void *ctx) { (void)ctx; J[0] = -20 * p[0], J[1] = 10, J[2] = -1, J[3] = 0; }
enum { M = 40 };
static void exp_r(const double *p, double *r, void *ctx)
{
    (void)ctx;
    for (int i = 0; i < M; i++) {
        double x = i * 0.1;
        r[i] = p[0] * exp(p[1] * x) + p[2] - (2 * exp(-1.5 * x) + 0.5);
    }
}
static double parab(double x, void *ctx) { (void)ctx; return (x - 2) * (x - 2); } /* minimum 0: f resolves x to the tolerance */

static double rosen_fg(const double *x, double *g, void *ctx)
{
    int n = *(int *)ctx;
    double f = 0;
    for (int i = 0; i < n; i++) g[i] = 0;
    for (int i = 0; i + 1 < n; i++) {
        double a = x[i + 1] - x[i] * x[i], b = 1 - x[i];
        f += 100 * a * a + b * b;
        g[i] += -400 * x[i] * a - 2 * b, g[i + 1] += 200 * a;
    }
    return f;
}
enum { EN = 10 };
/* 1.0.2 regression helpers, from the verification work's lbfgs_demos.c and cmaes_demos.c (Result 71) */
static double fg_nan2(const double *x, double *g, void *c) { (void)x, (void)c; g[0] = NAN, g[1] = 0; return NAN; }
static double fg_scaled2(const double *x, double *g, void *c)
{
    double s = *(double *)c;
    g[0] = s * 2 * x[0], g[1] = s * 20 * x[1];
    return s * (x[0] * x[0] + 10 * x[1] * x[1]);
}
static double fg_quartic4(const double *x, double *g, void *c) /* the verification work's lbfgs_yy_underflow.c: f = |x|^4, n 4 */
{
    (void)c;
    long double q = 0;
    for (int j = 0; j < 4; j++) q += (long double)x[j] * x[j];
    for (int j = 0; j < 4; j++) g[j] = (double)(4 * q * x[j]);
    return (double)(q * q);
}
static long fg_calls;
static double fg_count1(const double *x, double *g, void *c) { (void)x, (void)c, (void)g; fg_calls++; return 1; }
static double sphn(const double *x, void *c) { int n = *(int *)c; double s = 0; for (int i = 0; i < n; i++) s += x[i] * x[i]; return s; }
static double sph_nan_start(const double *x, void *c) { return x[0] == 3 ? NAN : sphn(x, c); }
static int region_nan = 1;
static double sph_region(const double *x, void *c) { return x[0] > 1 ? (region_nan ? NAN : INFINITY) : sphn(x, c); }
static double Rm[EN * EN];
static double ellip(const double *x, void *ctx)
{
    (void)ctx;
    double f = 0;
    for (int i = 0; i < EN; i++) {
        double y = 0;
        for (int j = 0; j < EN; j++) y += Rm[i * EN + j] * x[j];
        f += pow(1e6, (double)i / (EN - 1)) * y * y;
    }
    return f;
}
static double rastrigin(const double *x, void *ctx)
{
    int n = *(int *)ctx;
    double f = 10.0 * n;
    for (int i = 0; i < n; i++) f += x[i] * x[i] - 10 * cos(2 * T_PI * x[i]);
    return f;
}
static double rosen2(const double *p, void *ctx) { (void)ctx; return 100 * (p[1] - p[0] * p[0]) * (p[1] - p[0] * p[0]) + (1 - p[0]) * (1 - p[0]); }
/* the verification work's findings 4 and 5: a kappa 1e7 linear fit, and Brown-Dennis and Penalty I from Moré, Garbow and Hillstrom */
enum { KM = 28 };
static double KA[KM * 2], KB[KM];
static void klin_r(const double *p, double *r, void *ctx) { (void)ctx; for (int i = 0; i < KM; i++) r[i] = KA[2 * i] * p[0] + KA[2 * i + 1] * p[1] - KB[i]; }
static void klin_j(const double *p, double *J, void *ctx) { (void)p, (void)ctx; memcpy(J, KA, sizeof KA); }
static void bden_r(const double *x, double *f, void *ctx)
{
    (void)ctx;
    for (int i = 1; i <= 20; i++) {
        double t = i / 5.0, a = x[0] + t * x[1] - exp(t), b = x[2] + x[3] * sin(t) - cos(t);
        f[i - 1] = a * a + b * b;
    }
}
static void bden_j(const double *x, double *J, void *ctx)
{
    (void)ctx;
    for (int i = 1; i <= 20; i++) {
        double t = i / 5.0, a = x[0] + t * x[1] - exp(t), b = x[2] + x[3] * sin(t) - cos(t), *row = J + 4 * (i - 1);
        row[0] = 2 * a, row[1] = 2 * a * t, row[2] = 2 * b, row[3] = 2 * b * sin(t);
    }
}
static int PN;
static void pen_r(const double *x, double *f, void *ctx)
{
    (void)ctx;
    double s = 0;
    for (int j = 0; j < PN; j++) f[j] = sqrt(1e-5) * (x[j] - 1), s += x[j] * x[j];
    f[PN] = s - 0.25;
}
static void pen_j(const double *x, double *J, void *ctx)
{
    (void)ctx;
    for (int i = 0; i <= PN; i++)
        for (int j = 0; j < PN; j++) J[i * PN + j] = i < PN ? (i == j ? sqrt(1e-5) : 0) : 2 * x[j];
}
int main(void)
{
    double p[2] = {-1.2, 1}, st[2] = {0.5, 0.5};
    long ev = 0;
    double f = mk_nm(rosen, NULL, 2, p, st, 10000, 1e-20, &ev);
    CHECK(f < 1e-12 && fabs(p[0] - 1) < 1e-5 && fabs(p[1] - 1) < 1e-5, "Nelder-Mead, Rosenbrock from (-1.2, 1): (%.7f, %.7f), f %.2g, %ld calls", p[0], p[1], f, ev);
    double q[4] = {0, 0, 0, 0}, s4[4] = {1, 1, 1, 1};
    double f4 = mk_nm(quad4, NULL, 4, q, s4, 20000, 1e-24, NULL);
    CHECK(f4 < 1e-10 && fabs(q[0] - 1) < 1e-4 && fabs(q[1] + 2) < 1e-4 && fabs(q[2] - 3) < 1e-4 && fabs(q[3] - 0.5) < 1e-4,
          "Nelder-Mead, a 4-D coupled quadratic: (%.5f, %.5f, %.5f, %.5f), f %.2g", q[0], q[1], q[2], q[3], f4);

    double r1[2] = {-1.2, 1}, r2[2] = {-1.2, 1};
    double c1 = mk_lm(rosen_r, rosen_j, NULL, 2, 2, r1, NULL), c2 = mk_lm(rosen_r, NULL, NULL, 2, 2, r2, NULL);
    CHECK(c1 < 1e-20 && fabs(r1[0] - 1) < 1e-10 && fabs(r1[1] - 1) < 1e-10, "Levenberg-Marquardt, Rosenbrock residuals, given Jacobian: (%.12f, %.12f), cost %.2g", r1[0], r1[1], c1);
    CHECK(c2 < 1e-20 && fabs(r2[0] - 1) < 1e-8 && fabs(r2[1] - 1) < 1e-8, "the same by finite differences: (%.12f, %.12f), cost %.2g", r2[0], r2[1], c2);
    double e[3] = {1, -1, 0};
    double ce = mk_lm(exp_r, NULL, NULL, M, 3, e, NULL);
    CHECK(ce < 1e-20 && fabs(e[0] - 2) < 1e-8 && fabs(e[1] + 1.5) < 1e-8 && fabs(e[2] - 0.5) < 1e-8, "Levenberg-Marquardt, a exp(b x) + c on 40 points: (%.10f, %.10f, %.10f), cost %.2g", e[0], e[1], e[2], ce);

    double fm, xm = mk_golden(parab, NULL, 0, 5, 1e-10, &fm);
    CHECK(fabs(xm - 2) < 1e-9 && fm < 1e-18, "golden section, (x - 2)^2 on [0, 5] to 1e-10: x %.12f, f %.2g", xm, fm);
    int n2 = 2, n10 = 10;
    double lx[2] = {-1.2, 1};
    mk_lbfgs_opts lo = {0};
    double lf = mk_lbfgs(rosen_fg, &n2, 2, lx, &lo);
    CHECK(lo.status >= 0 && lo.status <= 1 && lf < 1e-18 && fabs(lx[0] - 1) < 1e-8 && fabs(lx[1] - 1) < 1e-8 && lo.iters_used < 60,
          "L-BFGS, Rosenbrock from (-1.2, 1): (%.12f, %.12f), f %.2g, %d iterations, %ld calls", lx[0], lx[1], lf, lo.iters_used, lo.evals);
    double lx10[10];
    for (int i = 0; i < 10; i++) lx10[i] = i % 2 ? 1 : -1.2;
    mk_lbfgs_opts lo10 = {0};
    double lf10 = mk_lbfgs(rosen_fg, &n10, 10, lx10, &lo10);
    double le = 0;
    for (int i = 0; i < 10; i++) le = fmax(le, fabs(lx10[i] - 1));
    CHECK(lo10.status >= 0 && lo10.status <= 1 && lf10 < 1e-16 && le < 1e-7 && lo10.iters_used < 200, "L-BFGS, 10-D chained Rosenbrock: worst |x - 1| %.2g, %d iterations, %ld calls", le, lo10.iters_used, lo10.evals);

    double cx[2] = {-1.2, 1};
    mk_cmaes_opts co = {.sigma0 = 0.5, .maxevals = 20000, .ftarget = 1e-14, .seed = 1};
    double cf = mk_cmaes(rosen2, NULL, 2, cx, &co);
    CHECK(cf < 1e-14 && fabs(cx[0] - 1) < 1e-6, "CMA-ES, Rosenbrock: (%.8f, %.8f), f %.2g, %ld calls", cx[0], cx[1], cf, co.evals);
    for (int i = 0; i < EN * EN; i++) Rm[i] = i % (EN + 1) == 0;
    for (int g = 0; g < 60; g++) {
        int p = g % EN, q = (p + 1 + (g / EN) % (EN - 1)) % EN;
        double th = t_rand() * T_PI, c = cos(th), s = sin(th);
        for (int k = 0; k < EN; k++) {
            double a = Rm[k * EN + p], b = Rm[k * EN + q];
            Rm[k * EN + p] = c * a - s * b, Rm[k * EN + q] = s * a + c * b;
        }
    }
    double ex[EN];
    for (int i = 0; i < EN; i++) ex[i] = 1;
    mk_cmaes_opts eo = {.sigma0 = 1, .maxevals = 60000, .ftarget = 1e-10, .seed = 2};
    double ef = mk_cmaes(ellip, NULL, EN, ex, &eo);
    CHECK(eo.status == 0 && ef < 1e-10, "CMA-ES, a rotated 10-D ellipsoid, condition 1e6: f %.2g in %ld calls", ef, eo.evals);
    int n3 = 3;
    double bestf = INFINITY, start[3] = {3.3, -2.7, 4.1};
    long evtot = 0;
    int runs = 0;
    for (int lam = 0, k = 0; k < 7 && bestf > 1e-8; k++, lam = lam ? 2 * lam : 14, runs++) {
        double rx[3] = {start[0], start[1], start[2]};
        mk_cmaes_opts ro = {.sigma0 = 2.5, .lambda = lam, .maxevals = 20000, .ftarget = 1e-10, .seed = 10 + (uint64_t)k};
        double rf = mk_cmaes(rastrigin, &n3, 3, rx, &ro);
        evtot += ro.evals;
        if (rf < bestf) bestf = rf;
    }
    CHECK(bestf < 1e-8, "CMA-ES with doubling populations (IPOP), 3-D Rastrigin: global minimum reached, f %.2g after %d runs, %ld calls", bestf, runs, evtot);
    double one[3] = {start[0], start[1], start[2]};
    mk_cmaes_opts so = {.sigma0 = 0.3, .maxevals = 20000, .ftarget = 1e-10, .seed = 3};
    double sf = mk_cmaes(rastrigin, &n3, 3, one, &so);
    printf("info  one small-step run from the same start stops in a local basin: f %.4g at (%.3f, %.3f, %.3f)\n", sf, one[0], one[1], one[2]);
    /* ---- the verification work's fuzz findings, 2026-11-14 (FINDINGS.md 4, 5) and the contract note on mk_nm ---- */
    {
        long double u1[KM], u2[KM], d = 0, n1 = 0, n2 = 0;
        for (int i = 0; i < KM; i++) u1[i] = t_rand(), n1 += u1[i] * u1[i];
        for (int i = 0; i < KM; i++) u1[i] /= sqrtl(n1), u2[i] = t_rand(), d += u1[i] * u2[i];
        for (int i = 0; i < KM; i++) u2[i] -= d * u1[i], n2 += u2[i] * u2[i];
        for (int i = 0; i < KM; i++) u2[i] /= sqrtl(n2);
        long double s1 = 1 / sqrtl(1 + 1e-14L), s2 = 1e-7L * s1, th = 0.7L, c = cosl(th), s = sinl(th);
        for (int i = 0; i < KM; i++) { /* A = U diag(s1, s2) V^T, V a rotation by 0.7 */
            KA[2 * i] = (double)(u1[i] * s1 * c - u2[i] * s2 * s);
            KA[2 * i + 1] = (double)(u1[i] * s1 * s + u2[i] * s2 * c);
            KB[i] = t_rand();
        }
        long double ata[3] = {0, 0, 0}, atb[2] = {0, 0}; /* the exact minimum, long-double normal equations */
        for (int i = 0; i < KM; i++)
            ata[0] += (long double)KA[2 * i] * KA[2 * i], ata[1] += (long double)KA[2 * i] * KA[2 * i + 1], ata[2] += (long double)KA[2 * i + 1] * KA[2 * i + 1],
                atb[0] += (long double)KA[2 * i] * KB[i], atb[1] += (long double)KA[2 * i + 1] * KB[i];
        long double det = ata[0] * ata[2] - ata[1] * ata[1], x0 = (ata[2] * atb[0] - ata[1] * atb[1]) / det, x1 = (ata[0] * atb[1] - ata[1] * atb[0]) / det, cs = 0;
        for (int i = 0; i < KM; i++) { long double e = KA[2 * i] * x0 + KA[2 * i + 1] * x1 - KB[i]; cs += e * e; }
        double kp[2] = {0, 0}, kc = mk_lm(klin_r, klin_j, NULL, KM, 2, kp, NULL), xn = (double)sqrtl(x0 * x0 + x1 * x1);
        double bound = 1e-8 * (double)cs; /* a converged fit is at the minimum to rounding: the patched one is within 2e-11 */
        CHECK(kc - (double)cs <= bound, "finding 4: a kappa 1e7 linear fit (28 x 2, |x*| %.2g) reaches its minimum: excess %.3g against a bound of 1e-8 of the minimum, %.3g", xn, kc - (double)cs, bound);
    }
    mk_lm_opts mo = {.iters = 1000};
    double bx[4] = {25, 5, -5, -1}, bcj = mk_lm(bden_r, bden_j, NULL, 20, 4, bx, &mo);
    double bx2[4] = {25, 5, -5, -1}, bcd = mk_lm(bden_r, NULL, NULL, 20, 4, bx2, &mo);
    CHECK(fabs(bcj - 85822.2) < 0.1 && fabs(bcd - 85822.2) < 0.1, "finding 5: Brown-Dennis reaches its listed minimum 85822.2 (Jacobian %.6g, differences %.6g)", bcj, bcd);
    PN = 4;
    double px[10] = {1, 2, 3, 4}, pcd = mk_lm(pen_r, NULL, NULL, 5, 4, px, &mo);
    PN = 10;
    double px2[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, pcj = mk_lm(pen_r, pen_j, NULL, 11, 10, px2, &mo);
    CHECK(fabs(pcd - 2.24997e-5) < 1e-9 && fabs(pcj - 7.08765e-5) < 1e-9, "finding 5: Penalty I reaches its listed minima (n 4 by differences %.6g, 2.24997e-5; n 10 with J %.6g, 7.08765e-5)", pcd, pcj);
    long nev = 12345;
    double nst[1] = {1}, nps[1] = {0};
    double nv = mk_nm(rosen, NULL, 0, nps, nst, 10, 1e-9, &nev);
    CHECK(isnan(nv) && nev == 0, "contract: Nelder-Mead with n out of range returns NaN and sets the call count to 0 (got %ld)", nev);
    { /* 1.0.2: the verification work's L-BFGS and CMA-ES findings (Result 71) */
        mk_lbfgs_opts o = {0};
        double x[2] = {1, 1};
        mk_lbfgs(fg_nan2, NULL, 2, x, &o);
        CHECK(o.status != 0, "lbfgs: a NaN gradient is not reported as converged (status %d)", o.status);
        double s = 1e160, y[2] = {1, 1};
        mk_lbfgs_opts p = {0};
        mk_lbfgs(fg_scaled2, &s, 2, y, &p);
        CHECK((p.status == 0 || p.status == 1) && fabs(y[0]) < 1e-6 && fabs(y[1]) < 1e-6, "lbfgs: |g| near 1e161 still converges (status %d, x %g %g)", p.status, y[0], y[1]);
        mk_lbfgs_opts q = {.m = 1061067738};
        double z[1] = {1};
        fg_calls = 0;
        mk_lbfgs(fg_count1, NULL, -1086567296, z, &q);
        CHECK(q.status == -1 && fg_calls == 0, "lbfgs: a negative n is refused before any call (status %d, calls %ld)", q.status, fg_calls);
        int n3 = 3;
        double c[3] = {3, 3, 3};
        mk_cmaes_opts co = {.sigma0 = 1, .ftarget = 1e-10, .maxevals = 5000, .seed = 1};
        double fc = mk_cmaes(sph_nan_start, &n3, 3, c, &co);
        CHECK(fc <= 1e-10, "cmaes: a NaN start value is replaced (f %g)", fc);
        int n5 = 5;
        double r[2];
        long ev[2];
        for (region_nan = 1; region_nan >= 0; region_nan--) {
            double y5[5] = {0.9, 2, -2, 1, 2};
            mk_cmaes_opts po = {.sigma0 = 2, .ftarget = 1e-10, .maxevals = 20000, .seed = 101};
            r[region_nan] = mk_cmaes(sph_region, &n5, 5, y5, &po), ev[region_nan] = po.evals;
        }
        CHECK(r[1] == r[0] && ev[1] == ev[0], "cmaes: a NaN region ranks like +Inf, bit for bit (f %g/%g, evals %ld/%ld)", r[1], r[0], ev[1], ev[0]);
    }
#ifndef T_OLD
    { /* 1.0.3: the verification work Result 74: 3 m wrapped an int, malloc got 16 bytes and res() wrote past it */
        double p3[3] = {0, 0, 0};
        CHECK(mk_lm(rosen_r, NULL, NULL, 715827883, 3, p3, NULL) == -1, "mk_lm refuses m = 715827883, n = 3 before allocating or calling res");
    }
#endif
    { /* 1.0.4: the verification work finding H (Result 78): y.y underflowed once |y| < 1e-154, gamma went to Inf and the search failed */
        double x4[4] = {0.7, -1.3, 0.2, 2.1};
        mk_lbfgs_opts o4 = {.gtol = 1e-300, .ftol = 1e-300, .iters = 3000};
        double f4 = mk_lbfgs(fg_quartic4, NULL, 4, x4, &o4);
        CHECK(o4.status >= 0 && f4 < 1e-290, "lbfgs: |x|^4 to f below 1e-290 without y.y underflow (status %d, f %g)", o4.status, f4);
    }
    { /* 1.6.0, the review's C1: lambda 1 has no parent (mu 0, mueff 1/0); 1.5.6 ran a generation and gave up with no reason */
        int n2 = 2;
        double x1[2] = {1, 1};
        mk_cmaes_opts o1 = {.sigma0 = 0.5, .lambda = 1, .seed = 1};
        const double f1 = mk_cmaes(sphn, &n2, 2, x1, &o1);
        CHECK(isnan(f1) && o1.status == -1 && o1.evals == 0 && x1[0] == 1 && x1[1] == 1, "C1 (review): cmaes refuses lambda 1 before evaluating (f %g, status %d, %ld evaluations)", f1,
              o1.status, o1.evals);
    }
    { /* review 2 #2 (1.6.1): a start simplex whose values tie stopped Nelder-Mead at its start (1.6.0: x^2 from -0.5, step 1, gave 0.25) */
        int n1 = 1, n2 = 2;
        double p1[1] = {-0.5}, s1[1] = {1}, p2[2] = {-0.5, -0.5}, s2[2] = {1, 1};
        long e1, e2;
        const double f1 = mk_nm(sphn, &n1, 1, p1, s1, 1000, 1e-9, &e1), f2 = mk_nm(sphn, &n2, 2, p2, s2, 1000, 1e-9, &e2);
        CHECK(f1 < 1e-8 && f2 < 1e-8, "review 2 #2: Nelder-Mead from tied starts: x^2 from -0.5 reaches %.2g (%ld calls), x^2 + y^2 from (-0.5, -0.5) %.2g (%ld calls)", f1, e1, f2, e2);
    }
    { /* review 2 #2 again (1.6.2): a tie reached later still stopped 1.6.1 (x^2 from -3.5, step 1: {-0.5, 0.5}, 0.25 after 6 calls; 14 of
         the 20 half-integer starts -9.5 to 9.5 stopped at 0.25) */
        int n1 = 1, n2 = 2, away = 0;
        double p1[1] = {-3.5}, s1[1] = {1}, p2[2] = {-3.5, -3.5}, s2[2] = {1, 1};
        long e1, e2;
        const double f1 = mk_nm(sphn, &n1, 1, p1, s1, 1000, 1e-9, &e1), f2 = mk_nm(sphn, &n2, 2, p2, s2, 1000, 1e-9, &e2);
        for (double x0 = -9.5; x0 <= 9.5; x0 += 1) {
            double q[1] = {x0};
            away += mk_nm(sphn, &n1, 1, q, s1, 1000, 1e-9, NULL) > 1e-8;
        }
        CHECK(f1 < 1e-8 && f2 < 1e-8 && away == 0, "review 2 #2, later ties: x^2 from -3.5 reaches %.2g (%ld calls), x^2 + y^2 from (-3.5, -3.5) %.2g (%ld calls); %d of 20 half-integer starts stop away",
              f1, e1, f2, e2, away);
    }
    { /* review 2 A5 (1.6.2): a refusal left the counts of the run before in a reused opts (cmaes evals, lbfgs iters_used and evals) */
        int n2 = 2;
        double x1[2] = {1, 1}, x2[2] = {1, 2};
        mk_cmaes_opts o1 = {.sigma0 = 0.5, .seed = 1};
        mk_cmaes(sphn, &n2, 2, x1, &o1);
        const long ran = o1.evals;
        o1.lambda = 1;
        const double f1 = mk_cmaes(sphn, &n2, 2, x1, &o1);
        mk_lbfgs_opts o2 = {0};
        mk_lbfgs(fg_quartic4, NULL, 4, (double[4]){0.7, -1.3, 0.2, 2.1}, &o2);
        const long ran2 = o2.evals;
        const double f2 = mk_lbfgs(fg_quartic4, NULL, -1, x2, &o2);
        CHECK(ran > 0 && isnan(f1) && o1.status == -1 && o1.evals == 0 && ran2 > 0 && isnan(f2) && o2.status == -1 && o2.evals == 0 && o2.iters_used == 0,
              "review 2 A5: a refusal after a run counts 0: cmaes lambda 1 %ld evaluations (the run %ld), lbfgs n -1 %ld and %d iterations (the run %ld)", o1.evals, ran, o2.evals,
              o2.iters_used, ran2);
    }
    return wrong;
}
