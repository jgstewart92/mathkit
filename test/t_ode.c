#include "../include/mk_ode.h"
#include "t.h"
#include <limits.h>

static void osc(double t, const double *y, double *dy, void *ctx) { (void)t, (void)ctx; dy[0] = y[1], dy[1] = -y[0]; }
static void gaussf(double t, const double *y, double *dy, void *ctx) { (void)ctx; dy[0] = -2 * t * y[0]; }
static void acc(const double *x, double *a, void *ctx) { (void)ctx; a[0] = -x[0]; }

static double rk4_err(double h)
{
    double y[2] = {1, 0}, w[MK_RK4_WORK(2)];
    int n = (int)lround(1 / h);
    for (int i = 0; i < n; i++) mk_rk4_step(osc, NULL, 2, i * h, y, h, w);
    return fabs(y[0] - cos(1.0));
}
static double verlet_err(double h)
{
    double x = 1, v = 0, a = -1, w;
    int n = (int)lround(1 / h);
    for (int i = 0; i < n; i++) mk_verlet_step(acc, NULL, 1, &x, &v, &a, h, &w);
    return fabs(x - cos(1.0));
}

static void rober(double t, const double *y, double *dy, void *ctx)
{
    (void)t, (void)ctx;
    dy[0] = -0.04 * y[0] + 1e4 * y[1] * y[2];
    dy[1] = 0.04 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1];
    dy[2] = 3e7 * y[1] * y[1];
}
static void rober_j(double t, const double *y, double *J, void *ctx)
{
    (void)t, (void)ctx;
    J[0] = -0.04, J[1] = 1e4 * y[2], J[2] = 1e4 * y[1];
    J[3] = 0.04, J[4] = -1e4 * y[2] - 6e7 * y[1], J[5] = -1e4 * y[1];
    J[6] = 0, J[7] = 6e7 * y[1], J[8] = 0;
}
static const double MU = 1000;
static void vdp(double t, const double *y, double *dy, void *ctx) { (void)t, (void)ctx; dy[0] = y[1], dy[1] = MU * (1 - y[0] * y[0]) * y[1] - y[0]; }
static void lin(double t, const double *y, double *dy, void *ctx) { (void)ctx; dy[0] = -2 * (y[0] - cos(t)) - sin(t); }
static double bdf_err(int order, double h)
{
    double y = 1, yp = 1, t = 0, hp = 0;
    mk_newton_ctl c = {1e-12, 1e-14, 10, 0, 0, 0};
    int n = (int)lround(1 / h);
    for (int i = 0; i < n; i++) {
        mk_bdf_step(lin, NULL, NULL, 1, i == 0 ? 1 : order, t, &y, &yp, h, hp, &c);
        t += h, hp = h;
    }
    return fabs(y - cos(1.0));
}
static void one(double t, const double *y, double *dy, void *ctx) { (void)t, (void)y, (void)ctx; dy[0] = 1; }
static void decay(double t, const double *y, double *dy, void *ctx) { (void)t, (void)ctx; dy[0] = -y[0]; }
static long cos_calls;
static double cos_tmax;
static void cos40(double t, const double *y, double *dy, void *ctx)
{
    (void)y, (void)ctx;
    if (cos_calls++ < 7 && t > cos_tmax) cos_tmax = t; /* the first step's seven calls */
    dy[0] = cos(40 * t);
}
/* review 2 (1.6.1): an f that turns NaN part-way */
static void sqrt1t(double t, const double *y, double *dy, void *ctx) { (void)y, (void)ctx; dy[0] = sqrt(1 - t); }
static void msqrt(double t, const double *y, double *dy, void *ctx) { (void)t, (void)ctx; dy[0] = -sqrt(y[0]); }
static void halfj(double t, const double *y, double *J, void *ctx) { (void)t, (void)y, (void)ctx; J[0] = -0.5; }
static int late_calls;
static void late(double t, const double *y, double *dy, void *ctx) { (void)t, (void)ctx; dy[0] = late_calls++ ? NAN : -y[0] * y[0] * y[0]; }
static void latej(double t, const double *y, double *J, void *ctx) { (void)t, (void)ctx; J[0] = -3 * y[0] * y[0]; }
/* review 2 (1.6.2): two components, the NaN in the first; a step of a few ulps of t */
static void sqrt1t2(double t, const double *y, double *dy, void *ctx) { (void)y, (void)ctx; dy[0] = sqrt(1 - t), dy[1] = 1; }
static void msqrt2(double t, const double *y, double *dy, void *ctx) { (void)t, (void)ctx; dy[0] = -sqrt(y[0]), dy[1] = -y[1]; }
static void halfj2(double t, const double *y, double *J, void *ctx) { (void)t, (void)y, (void)ctx; J[0] = -0.5, J[1] = J[2] = 0, J[3] = -1; }
static void late2(double t, const double *y, double *dy, void *ctx) { (void)t, (void)ctx; dy[0] = late_calls++ ? NAN : -y[0], dy[1] = -y[1]; }
static void late2j(double t, const double *y, double *J, void *ctx) { (void)t, (void)y, (void)ctx; J[0] = -1, J[1] = J[2] = 0, J[3] = -1; }
static double lin_k;
static void dec05(double t, const double *y, double *dy, void *ctx) { (void)t, (void)ctx; dy[0] = lin_k * y[0]; }
static void dec05j(double t, const double *y, double *J, void *ctx) { (void)t, (void)y, (void)ctx; J[0] = lin_k; }
int main(void)
{
    double e1 = rk4_err(0.1), e2 = rk4_err(0.05), e3 = rk4_err(0.025);
    double o1 = log2(e1 / e2), o2 = log2(e2 / e3);
    CHECK(o1 > 3.8 && o1 < 4.3 && o2 > 3.8 && o2 < 4.3, "RK4 on x'' = -x to t 1: errors %.2g %.2g %.2g, order %.2f %.2f (4)", e1, e2, e3, o1, o2);

    double y[2] = {1, 0}, w[MK_DOPRI_WORK(2)], h = 0;
    mk_dopri_ctl c = {1e-10, 1e-12, 0, 0, 0, 0, 0};
    int r = mk_dopri(osc, NULL, 2, 0, y, 20 * T_PI, &h, w, &c);
    double eo = hypot(y[0] - 1, y[1]);
    CHECK(!r && eo < 1e-8, "Dormand-Prince, 10 periods at rtol 1e-10: error %.2g, %ld steps, %ld rejected, %ld calls", eo, c.accepted, c.rejected, c.evals);
    double errs[3], tols[3] = {1e-5, 1e-8, 1e-11};
    long steps[3];
    for (int k = 0; k < 3; k++) {
        double g = 1, hg = 0, wg[MK_DOPRI_WORK(1)];
        mk_dopri_ctl cg = {tols[k], tols[k] * 1e-3, 0, 0, 0, 0, 0};
        mk_dopri(gaussf, NULL, 1, 0, &g, 3, &hg, wg, &cg);
        errs[k] = fabs(g - exp(-9.0)) / exp(-9.0) * exp(-9.0) + fabs(g - exp(-9.0));
        steps[k] = cg.accepted;
    }
    CHECK(errs[0] > errs[1] && errs[1] > errs[2] && errs[2] < 1e-11, "y' = -2 t y to t 3, tolerance 1e-5, 1e-8, 1e-11: errors %.2g %.2g %.2g, steps %ld %ld %ld",
          errs[0], errs[1], errs[2], steps[0], steps[1], steps[2]);
    double yb[2] = {cos(5.0), -sin(5.0)}, hb = 0;
    mk_dopri_ctl cb = {1e-10, 1e-12, 0, 0, 0, 0, 0};
    mk_dopri(osc, NULL, 2, 5, yb, 0, &hb, w, &cb);
    CHECK(hypot(yb[0] - 1, yb[1]) < 1e-8, "backward in time, 5 to 0: error %.2g", hypot(yb[0] - 1, yb[1]));

    double v1 = verlet_err(0.1), v2 = verlet_err(0.05), ov = log2(v1 / v2);
    CHECK(ov > 1.8 && ov < 2.3, "velocity Verlet on x'' = -x: errors %.2g %.2g, order %.2f (2)", v1, v2, ov);
    double x = 1, v = 0, a = -1, wv, emax = 0;
    for (int i = 0; i < 100000; i++) {
        mk_verlet_step(acc, NULL, 1, &x, &v, &a, 0.1, &wv);
        emax = fmax(emax, fabs(0.5 * (x * x + v * v) - 0.5));
    }
    CHECK(emax < 0.003, "Verlet, 100000 steps of 0.1 (1600 periods): energy error stays under %.3g of 0.5", emax);
    /* Robertson's kinetics to t 40 (reference: Hairer and Wanner, Solving ODEs II) */
    const double ref[3] = {0.7158270687193941, 9.185534764557338e-06, 0.2841637457401608};
    double ry[3] = {1, 0, 0}, rh = 0;
    mk_stiff_ctl sc = {1e-7, 1e-12, 0, 1, 0, 0, 0, 0, 0, 200000};
    int rr = mk_ode23s(rober, rober_j, NULL, 3, 0, ry, 40, &rh, &sc);
    double re1 = fabs(ry[0] - ref[0]) / ref[0], re2 = fabs(ry[1] - ref[1]) / ref[1], re3 = fabs(ry[2] - ref[2]) / ref[2];
    CHECK(!rr && re1 < 1e-5 && re2 < 1e-4 && re3 < 1e-5 && fabs(ry[0] + ry[1] + ry[2] - 1) < 1e-13,
          "ode23s, Robertson to t 40: relative errors %.2g %.2g %.2g, mass %.2g off 1, %ld steps, %ld rejected", re1, re2, re3, ry[0] + ry[1] + ry[2] - 1, sc.accepted, sc.rejected);
    double dy3[3] = {1, 0, 0}, dh = 0, dw[MK_DOPRI_WORK(3)];
    mk_dopri_ctl dc = {1e-7, 1e-12, 0, 0, 0, 0, 0};
    mk_dopri(rober, NULL, 3, 0, dy3, 40, &dh, dw, &dc);
    CHECK(sc.accepted * 5 < dc.accepted, "Dormand-Prince on the same problem and tolerance needs %ld steps (%ld calls) against ode23s's %ld (%ld calls): %.0f times",
          dc.accepted, dc.evals, sc.accepted, sc.evals, (double)dc.accepted / sc.accepted);

    /* Van der Pol, mu 1000: the period against its asymptotic value (3 - 2 ln 2) mu + 3 a mu^(-1/3), a = 2.338107 */
    double vy[2] = {2, 0}, vh = 0, t = 0, cross[4];
    int nc = 0;
    mk_stiff_ctl vc = {1e-6, 1e-6, 0, 1, 0, 0, 0, 0, 0, 200000};
    while (t < 3000 && nc < 2) {
        double y0 = vy[0], t0 = t, t1 = t + 1;
        if (mk_ode23s(vdp, NULL, NULL, 2, t0, vy, t1, &vh, &vc)) break;
        if (y0 > 0 && vy[0] <= 0) cross[nc++] = t0 + y0 / (y0 - vy[0]);
        t = t1;
    }
    double per = nc == 2 ? cross[1] - cross[0] : 0, want = (3 - 2 * log(2.0)) * MU + 3 * 2.338107 * pow(MU, -1.0 / 3);
    CHECK(nc == 2 && fabs(per - want) < 1e-3 * want, "ode23s, Van der Pol mu 1000 (numerical Jacobian): period %.2f against %.2f, %ld steps", per, want, vc.accepted);
    double dv[2] = {2, 0}, dvh = 0, dvw[MK_DOPRI_WORK(2)];
    mk_dopri_ctl dvc = {1e-6, 1e-6, 0, 0, 0, 0, 0};
    mk_stiff_ctl vc2 = {1e-6, 1e-6, 0, 1, 0, 0, 0, 0, 0, 200000};
    double w2[2] = {2, 0}, w2h = 0;
    mk_ode23s(vdp, NULL, NULL, 2, 0, w2, 3000, &w2h, &vc2);
    mk_dopri(vdp, NULL, 2, 0, dv, 3000, &dvh, dvw, &dvc);
    CHECK(vc2.accepted * 20 < dvc.accepted && fabs(dv[0] - w2[0]) < 1e-2 * fabs(w2[0]), "to t 3000: Dormand-Prince %ld steps against ode23s %ld (%.0f times); both end at y1 %.4f and %.4f",
          dvc.accepted, vc2.accepted, (double)dvc.accepted / vc2.accepted, dv[0], w2[0]);

    double b1a = bdf_err(1, 0.01), b1b = bdf_err(1, 0.005), b2a = bdf_err(2, 0.01), b2b = bdf_err(2, 0.005);
    double q1 = log2(b1a / b1b), q2 = log2(b2a / b2b);
    CHECK(q1 > 0.9 && q1 < 1.1 && q2 > 1.8 && q2 < 2.2, "BDF1 and BDF2 on y' = -2 (y - cos t) - sin t: orders %.2f and %.2f (1, 2)", q1, q2);
    double by[3] = {1, 0, 0}, bp[3] = {1, 0, 0}, bt = 0, bh = 1e-6, bhp = 0;
    mk_newton_ctl nc2 = {1e-6, 1e-12, 10, 0, 0, 0};
    long bsteps = 0;
    int bfail = 0;
    while (bt < 40 && !bfail) {
        double h = fmin(bh, 40 - bt);
        if (mk_bdf_step(rober, rober_j, NULL, 3, bsteps ? 2 : 1, bt, by, bp, h, bhp, &nc2)) { bfail = 1; break; }
        bt += h, bhp = h, bh *= 1.05, bsteps++;
    }
    double be1 = fabs(by[0] - ref[0]) / ref[0], be3 = fabs(by[2] - ref[2]) / ref[2];
    CHECK(!bfail && be1 < 2e-3 && be3 < 5e-3, "BDF2, Robertson with steps growing 5 %% from 1e-6: %ld steps, relative errors %.2g %.2g, %ld Newton iterations", bsteps, be1, be3, nc2.iters);
    { /* 1.0.5: the verification work's ode findings (Result 79) and the NaN times left to us */
        int lbad = 0;
        double lworst = 0;
        for (int k = 2; k <= 50; k++) {
            double y = 0, h = 1.0 / k, w[MK_DOPRI_WORK(1)], y2 = 0, h2 = 1.0 / k;
            mk_dopri_ctl c = {.rtol = 1e-6, .atol = 1e-9, .hmax = 1.0 / k};
            mk_stiff_ctl s = {.rtol = 1e-6, .atol = 1e-9, .hmax = 1.0 / k, .autonomous = 1};
            int r = mk_dopri(one, NULL, 1, 0.0, &y, 1.0, &h, w, &c), r2 = mk_ode23s(one, NULL, NULL, 1, 0.0, &y2, 1.0, &h2, &s);
            if (r || r2) lbad++;
            lworst = fmax(lworst, fmax(fabs(y - 1), fabs(y2 - 1)));
        }
        CHECK(lbad == 0 && lworst < 1e-14, "O-A: y' = 1 over [0, 1] with hmax 1/k, k 2..50: %d of 49 fail in Dormand-Prince or ode23s, y(1) off by %.2g", lbad, lworst);
        double cy = 0, ch = NAN, cw[MK_DOPRI_WORK(1)];
        mk_dopri_ctl cc = {.rtol = 1e-6, .atol = 1e-9, .hmax = 0.01};
        cos_calls = 0, cos_tmax = 0;
        int cr = mk_dopri(cos40, NULL, 1, 0.0, &cy, 1.0, &ch, cw, &cc);
        CHECK(cr == 0 && cos_tmax <= 0.01 && fabs(cy - sin(40.0) / 40) < 1e-6,
              "O-D: Dormand-Prince with *h NaN and hmax 0.01: the first step reaches t %.3g (at most 0.01), y(1) %.9f (%.9f)", cos_tmax, cy, sin(40.0) / 40);
        double ny = 1, nh = 0, nw[MK_DOPRI_WORK(1)], tnan[5] = {0};
        mk_dopri_ctl nc = {.rtol = 1e-6, .atol = 1e-9};
        mk_stiff_ctl ns = {.rtol = 1e-6, .atol = 1e-9, .autonomous = 1};
        mk_newton_ctl nn = {1e-6, 1e-9, 10, 0, 0, 0};
        int o1 = mk_dopri(decay, NULL, 1, 0.0, &ny, NAN, &nh, nw, &nc), o2 = mk_dopri(decay, NULL, 1, NAN, &ny, 1.0, &nh, nw, &nc);
        int o3 = mk_ode23s(decay, NULL, NULL, 1, 0.0, &ny, NAN, &nh, &ns), o4 = mk_ode23s(decay, NULL, NULL, 1, NAN, &ny, 1.0, &nh, &ns);
        double yp = 1;
        int o5 = mk_bdf_step(decay, NULL, NULL, 1, 1, NAN, &ny, &yp, 0.1, 0, &nn), o6 = mk_bdf_step(decay, NULL, NULL, 1, 1, 0, &ny, &yp, NAN, 0, &nn);
        tnan[0] = (double)(nc.evals + ns.evals + nn.evals);
        CHECK(o1 == -1 && o2 == -1 && o3 == -1 && o4 == -1 && o5 == -1 && o6 == -1 && ny == 1 && tnan[0] == 0,
              "O-C: a NaN t1 or t0 is refused (-1) by Dormand-Prince (%d %d) and ode23s (%d %d), a NaN t or h by BDF (%d %d), before f is called (%g calls)",
              o1, o2, o3, o4, o5, o6, tnan[0]);
#ifndef T_OLD
        double by1[1] = {0}, bw1[16], bh1 = 0;
        mk_dopri_ctl bc = {.rtol = 1e-6, .atol = 1e-9};
        CHECK(mk_dopri(one, NULL, -1, 0, by1, 1, &bh1, bw1, &bc) == -1, "O-B: Dormand-Prince refuses n = -1 (1.0.4 called memcpy with a negative size)");
#endif
    }
    { /* the verification work O-E' (Result 86): maxsteps = LONG_MAX; 1.2.2 computed LONG_MAX + 1 only on its LONG_MAX-th step, so this row
         passes there too: the overflow is out of a test's reach */
        double y1[1] = {1}, h1 = 0, w1[MK_DOPRI_WORK(1)], y2[1] = {1}, h2 = 0;
        mk_dopri_ctl c1 = {1e-9, 1e-12, 0, 0, 0, 0, LONG_MAX};
        mk_stiff_ctl c2 = {1e-6, 1e-9, 0, 0, 0, 0, 0, 0, 0, LONG_MAX};
        const int r1 = mk_dopri(gaussf, NULL, 1, 0, y1, 1, &h1, w1, &c1), r2 = mk_ode23s(gaussf, NULL, NULL, 1, 0, y2, 1, &h2, &c2);
        CHECK(r1 == 0 && r2 == 0 && fabs(y1[0] - exp(-1)) < 1e-8 && fabs(y2[0] - exp(-1)) < 1e-4, "O-E': maxsteps LONG_MAX: both solvers run to the end (y %.9f, %.6f; e^-1 %.9f)", y1[0], y2[0], exp(-1));
    }
    { /* 1.6.0, the review's O1: ode23s integrates forward only; 1.5.6 answered t1 < t0 with 0 and y unchanged */
        double yb[1] = {1}, hb = 0, ye[1] = {1}, he = 0;
        mk_stiff_ctl cb = {.rtol = 1e-6, .atol = 1e-9}, ce = {.rtol = 1e-6, .atol = 1e-9};
        const int rb = mk_ode23s(gaussf, NULL, NULL, 1, 5, yb, 0, &hb, &cb), re = mk_ode23s(gaussf, NULL, NULL, 1, 5, ye, 5, &he, &ce);
        CHECK(rb == -1 && yb[0] == 1 && re == 0 && ye[0] == 1, "O1 (review): ode23s from t 5 back to 0 is refused (%d), y as it was; t1 = t0 is 0 (%d)", rb, re);
    }
    { /* review 2 #1 (1.6.1): the error norm kept a NaN out (fmax), so 1.6.0 took the NaN step: dopri and ode23s returned 0 with y NaN,
         bdf_step called Newton converged */
        double y[1] = {0}, h = 0, w[MK_DOPRI_WORK(1)], z[1] = {1}, hz = 0, b[1] = {1}, bp[1] = {1};
        mk_dopri_ctl c = {.rtol = 1e-8, .atol = 1e-10};
        mk_stiff_ctl s = {.rtol = 1e-6, .atol = 1e-9, .autonomous = 1};
        mk_newton_ctl nc = {.rtol = 1e-8, .atol = 1e-10};
        const int rd = mk_dopri(sqrt1t, NULL, 1, 0, y, 2, &h, w, &c), rs = mk_ode23s(msqrt, halfj, NULL, 1, 0, z, 3, &hz, &s);
        late_calls = 0;
        const int rb = mk_bdf_step(late, latej, NULL, 1, 1, 0, b, bp, 0.5, 0, &nc);
        CHECK(rd == -1 && fabs(y[0] - 2.0 / 3) < 1e-6 && rs == -1 && z[0] >= 0 && z[0] < 1e-6 && rb == -1 && b[0] == 1,
              "review 2 #1: f NaN past a point: dopri on sqrt(1 - t) to 2 is %d at y %.9f (2/3), ode23s on -sqrt(y) to 3 is %d at y %.3g; bdf_step with f NaN on its second call %d, y %g as it was",
              rd, y[0], rs, z[0], rb, b[0]);
    }
    { /* review 2 #8 (1.6.1): each step is (t + s) - t, the one t can take; 1.6.0 moved y by s where t's ulp rounded the step */
        double ya[1] = {0}, yb[1] = {0}, yc[1] = {1}, ha = 0, hb = 0, hc = 0, w[MK_DOPRI_WORK(1)], za[1] = {0}, hz = 0;
        mk_dopri_ctl ca = {.rtol = 1e-8, .atol = 1e-12}, cb = ca, cc = ca;
        mk_stiff_ctl s = {.rtol = 1e-6, .atol = 1e-12, .autonomous = 1};
        const double t0 = 1.7e9, t1 = t0 + 0.01, span = t1 - t0; /* epoch seconds: ulp 2.4e-7 */
        const int ra = mk_dopri(one, NULL, 1, t0, ya, t1, &ha, w, &ca), rz = mk_ode23s(one, NULL, NULL, 1, t0, za, t1, &hz, &s);
        const int rb = mk_dopri(one, NULL, 1, 1e17, yb, 1e17 + 1000, &hb, w, &cb), rc = mk_dopri(decay, NULL, 1, 1e17, yc, 1e17 + 1000, &hc, w, &cc);
        CHECK(ra == 0 && fabs(ya[0] - span) <= 1e-14 * span && rz == 0 && fabs(za[0] - span) <= 1e-14 * span && rb == 0 && fabs(yb[0] - 992) <= 1e-12 && rc == -1 && yc[0] == 1,
              "review 2 #8: y' = 1 from t 1.7e9 over 0.01: dopri %d, y off the span by %.2g of it, ode23s %d, %.2g; from 1e17 over 1000 (ulp 16, span 992): %d, y %.17g; y' = -y there refused (%d, y %g)",
              ra, fabs(ya[0] - span) / span, rz, fabs(za[0] - span) / span, rb, yb[0], rc, yc[0]);
    }
    { /* review 2 #1 again (1.6.2): 1.6.1's max kept a NaN only in the last component, so dopri took the step with y[0] NaN (0) */
        double y[2] = {0, 0}, h = 0, w[MK_DOPRI_WORK(2)], z[2] = {1, 1}, hz = 0, b[2] = {1, 1}, bp[2] = {1, 1};
        mk_dopri_ctl c = {.rtol = 1e-8, .atol = 1e-10};
        mk_stiff_ctl s = {.rtol = 1e-6, .atol = 1e-9, .autonomous = 1};
        mk_newton_ctl nc = {.rtol = 1e-8, .atol = 1e-10};
        const int rd = mk_dopri(sqrt1t2, NULL, 2, 0, y, 2, &h, w, &c), rs = mk_ode23s(msqrt2, halfj2, NULL, 2, 0, z, 3, &hz, &s);
        late_calls = 0;
        const int rb = mk_bdf_step(late2, late2j, NULL, 2, 1, 0, b, bp, 0.5, 0, &nc);
        CHECK(rd == -1 && fabs(y[0] - 2.0 / 3) < 1e-6 && fabs(y[1] - 1) < 1e-6 && rs == -1 && z[0] >= 0 && z[0] < 1e-6 && rb == -1 && b[0] == 1 && b[1] == 1,
              "review 2 #1, n 2 with the NaN in component 0: dopri on (sqrt(1 - t), 1) to 2 is %d at y (%.9f, %.9f); ode23s %d at y0 %.3g; bdf_step %d, y (%g, %g) as it was",
              rd, y[0], y[1], rs, z[0], rb, b[0], b[1]);
    }
    { /* review 2 A1 (1.6.2): after a rejection the snapped step rounded back up to the step rejected, so 1.6.1 retried it to maxsteps:
         y' = -0.05 y from t 1e17 (ulp 16) over 992. dopri at rtol 1e-2 stuck at 32 (1.6.1), now finishes; at rtol 1e-4 the step needed is
         under t's resolution: -1 at once. ode23s with W = 1 - 16 d J exactly 0 (its LU refused) from t0 1e17 + 16, whose t + 8 ties up to 16 */
        const double t0 = 1e17, t1 = t0 + 992, d = 1 / (2 + sqrt(2.0));
        double ya[1] = {1}, yb[1] = {1}, yc[1] = {1}, ha = 0, hb = 0, hc = 16, w[MK_DOPRI_WORK(1)];
        mk_dopri_ctl ca = {.rtol = 1e-2, .atol = 1e-12, .maxsteps = 100000}, cb = {.rtol = 1e-4, .atol = 1e-12, .maxsteps = 100000};
        mk_stiff_ctl s = {.rtol = 1e-3, .atol = 1e-12, .autonomous = 1, .maxsteps = 100000};
        lin_k = -0.05;
        const int ra = mk_dopri(dec05, NULL, 1, t0, ya, t1, &ha, w, &ca), rb = mk_dopri(dec05, NULL, 1, t0, yb, t1, &hb, w, &cb);
        lin_k = 1 / (16 * d);
        for (int k = 0; k < 64 && -16 * d * lin_k + 1 != 0; k++) lin_k = nextafter(lin_k, -16 * d * lin_k + 1 > 0 ? INFINITY : 0);
        const int rc = mk_ode23s(dec05, dec05j, NULL, 1, t0 + 16, yc, t1 + 16, &hc, &s);
        CHECK(ra == 0 && fabs(ya[0]) < 1e-11 && ca.rejected < 100 && rb == -1 && cb.rejected < 100 && -16 * d * lin_k + 1 == 0 && rc == -1 && s.rejected < 10,
              "review 2 A1: steps of a few ulps of t 1e17: dopri at rtol 1e-2 %d (y %.2g, e^-49.6 = 2.9e-22, %ld rejected), at 1e-4 %d after %ld rejections; ode23s with its LU refused at 16: %d after %ld",
              ra, ya[0], ca.rejected, rb, cb.rejected, rc, s.rejected);
    }
    return wrong;
}
