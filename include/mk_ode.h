/* mathkit: time steppers for y' = f(t, y). Classic RK4 (fixed step); Dormand-Prince 5(4) with adaptive steps, mixed
   absolute and relative tolerance and first-same-as-last reuse; velocity Verlet for x'' = a(x), symplectic, so energy
   errors stay bounded over long runs; for stiff problems, Shampine's ode23s (Rosenbrock, adaptive, L-stable) and
   BDF1/BDF2 steps solved by Newton. The explicit steppers take caller-owned work arrays sized by the macros; the stiff
   ones allocate per call. Header only, C23; the stiff part uses mk_linalg.h. */
#ifndef MK_ODE_H
#define MK_ODE_H
#include "mk_linalg.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef void (*mk_odefn)(double t, const double *y, double *dy, void *ctx);

#define MK_RK4_WORK(n) (5 * (n))
/* one step of h from t; y updated in place; work holds MK_RK4_WORK(n) doubles */
static inline void mk_rk4_step(mk_odefn f, void *ctx, int n, double t, double *y, double h, double *work)
{
    double *k1 = work, *k2 = k1 + n, *k3 = k2 + n, *k4 = k3 + n, *z = k4 + n;
    f(t, y, k1, ctx);
    for (int i = 0; i < n; i++) z[i] = y[i] + 0.5 * h * k1[i];
    f(t + 0.5 * h, z, k2, ctx);
    for (int i = 0; i < n; i++) z[i] = y[i] + 0.5 * h * k2[i];
    f(t + 0.5 * h, z, k3, ctx);
    for (int i = 0; i < n; i++) z[i] = y[i] + h * k3[i];
    f(t + h, z, k4, ctx);
#if defined(MK_MUTATE) && MK_MUTATE == 51
    for (int i = 0; i < n; i++) y[i] += h / 4 * (k1[i] + k2[i] + k3[i] + k4[i]); /* planted fault: equal weights */
#else
    for (int i = 0; i < n; i++) y[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
#endif
}

typedef struct {
    double rtol, atol; /* error per component <= atol + rtol * |y| */
    /* hmax 0: no limit; met to within half an ulp of t (a step is one t can take: hmax 9 at t 1e17 steps 16) */
    double hmax;
    /* the counts are added to by every call (a run in pieces counts as one): zero them to count one call */
    long accepted, rejected, evals;
    long maxsteps;     /* accepted plus rejected steps per call before giving up (-1); 0 = 10 million */
} mk_dopri_ctl;

#define MK_DOPRI_WORK(n) (9 * (n))
/* Integrate y from t0 to t1 (either direction). *h: the first step to try (0 picks one) and, after, the last one taken,
   to carry on from. work holds MK_DOPRI_WORK(n) doubles. Each step is one t can take, (t + s) - t, so y moves as far as
   t does, and a step after a rejection is smaller than the one rejected. 0, or -1 if the step fell below 1e-14 of the
   span or below t's resolution, f gave a NaN (in any component) the steps could not shrink past, maxsteps ran out, t0
   or t1 is not finite, or n < 0. */
static inline int mk_dopri(mk_odefn f, void *ctx, int n, double t0, double *y, double t1, double *h, double *work, mk_dopri_ctl *c)
{
    static const double c2 = 1. / 5, c3 = 3. / 10, c4 = 4. / 5, c5 = 8. / 9;
    static const double a21 = 1. / 5, a31 = 3. / 40, a32 = 9. / 40, a41 = 44. / 45, a42 = -56. / 15, a43 = 32. / 9,
                        a51 = 19372. / 6561, a52 = -25360. / 2187, a53 = 64448. / 6561, a54 = -212. / 729,
                        a61 = 9017. / 3168, a62 = -355. / 33, a63 = 46732. / 5247, a64 = 49. / 176, a65 = -5103. / 18656,
                        b1 = 35. / 384, b3 = 500. / 1113, b4 = 125. / 192, b5 = -2187. / 6784, b6 = 11. / 84;
#if defined(MK_MUTATE) && MK_MUTATE == 52
    static const double e1 = 0, e3 = 0, e4 = 0, e5 = 0, e6 = 0, e7 = 0; /* planted fault: no error estimate */
#else
    static const double e1 = 71. / 57600, e3 = -71. / 16695, e4 = 71. / 1920, e5 = -17253. / 339200, e6 = 22. / 525, e7 = -1. / 40;
#endif
    double *k1 = work, *k2 = k1 + n, *k3 = k2 + n, *k4 = k3 + n, *k5 = k4 + n, *k6 = k5 + n, *k7 = k6 + n, *z = k7 + n, *yn = z + n;
    double span = t1 - t0, dir = span >= 0 ? 1 : -1, t = t0;
    /* a NaN or infinite *h picks one too, of at least 16 ulps of t0 */
    double hh = *h != 0 && isfinite(*h) ? fabs(*h) : fmax(fabs(span) * 1e-3, 16 * (nextafter(fabs(t0), INFINITY) - fabs(t0)));
    if (c->hmax > 0 && hh > c->hmax) hh = c->hmax;
    if (!isfinite(span)) return -1; /* a NaN or infinite t0 or t1 */
    if (span == 0) return 0;
    if (n < 0) { *h = dir * hh; return -1; }
    f(t, y, k1, ctx), c->evals++;
    long cap = c->maxsteps > 0 ? c->maxsteps : 10000000, used = 0;
    double hrej = INFINITY; /* the step just rejected */
    while (dir * (t1 - t) > 0) {
        if (used >= cap) { /* compare, then count: used never passes cap */ *h = dir * hh; return -1; }
        used++;
        double rem = dir * (t1 - t), hs = rem - hh < 1e-12 * fabs(span) ? rem : hh, err; /* no sliver left before t1 */
        if (hs != rem) hs = dir * ((t + dir * hs) - t); /* the step t can take (0 below t's resolution) */
        /* smaller than the step rejected, which the snap can round up to */
        if (hs >= hrej) hs = dir * (nextafter(t + dir * hrej, t) - t);
        if (hs < 1e-14 * fabs(span)) { *h = dir * hh; return -1; }
        double s = dir * hs;
        for (int i = 0; i < n; i++) z[i] = y[i] + s * a21 * k1[i];
        f(t + c2 * s, z, k2, ctx);
        for (int i = 0; i < n; i++) z[i] = y[i] + s * (a31 * k1[i] + a32 * k2[i]);
        f(t + c3 * s, z, k3, ctx);
        for (int i = 0; i < n; i++) z[i] = y[i] + s * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
        f(t + c4 * s, z, k4, ctx);
        for (int i = 0; i < n; i++) z[i] = y[i] + s * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
        f(t + c5 * s, z, k5, ctx);
        for (int i = 0; i < n; i++) z[i] = y[i] + s * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
        f(t + s, z, k6, ctx);
        for (int i = 0; i < n; i++) yn[i] = y[i] + s * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
        f(t + s, yn, k7, ctx);
        c->evals += 6;
        err = 0;
        for (int i = 0; i < n; i++) {
            double ei = s * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
            double sc = c->atol + c->rtol * fmax(fabs(y[i]), fabs(yn[i]));
            const double q = fabs(ei) / sc;
            if (q > err || q != q) err = q; /* a NaN in any component stays NaN (fmax would drop it) */
        }
        if (!isfinite(err)) err = 1e10;
        hh = hs * fmin(5, fmax(0.2, 0.9 * pow(err + 1e-300, -0.2)));
        if (c->hmax > 0 && hh > c->hmax) hh = c->hmax;
        if (err > 1) { c->rejected++, hrej = hs; continue; }
        c->accepted++, hrej = INFINITY;
        t = dir * (t1 - t - s) <= 0 ? t1 : t + s;
        memcpy(y, yn, sizeof *y * (size_t)n);
        memcpy(k1, k7, sizeof *k1 * (size_t)n); /* first same as last */
    }
    *h = dir * hh;
    return 0;
}

/* velocity Verlet for x'' = acc(x): x and v of n, a the acceleration at x on entry and on return (call acc once before
   the first step). work holds n doubles. Second order and symplectic. */
typedef void (*mk_accfn)(const double *x, double *a, void *ctx);
static inline void mk_verlet_step(mk_accfn acc, void *ctx, int n, double *x, double *v, double *a, double h, double *work)
{
#if defined(MK_MUTATE) && MK_MUTATE == 53
    /* planted fault: explicit Euler, whose energy grows */
    for (int i = 0; i < n; i++) x[i] += h * v[i], v[i] += h * a[i];
    acc(x, a, ctx);
    (void)work;
#else
    for (int i = 0; i < n; i++) v[i] += 0.5 * h * a[i], x[i] += h * v[i];
    acc(x, work, ctx);
    for (int i = 0; i < n; i++) v[i] += 0.5 * h * work[i], a[i] = work[i];
#endif
}

/* ---- stiff problems ---------------------------------------------------------------------------------------------- */

/* the Jacobian J[i * n + j] = d f_i / d y_j at (t, y) */
typedef void (*mk_odejac)(double t, const double *y, double *J, void *ctx);

static inline void mk_numjac_(mk_odefn f, void *ctx, int n, double t, double *y, const double *fy, double *J, double *tmp)
{
    for (int j = 0; j < n; j++) {
        double keep = y[j], d = 1.4901161193847656e-8 * fmax(fabs(keep), 1e-5);
        y[j] = keep + d;
        f(t, y, tmp, ctx);
        y[j] = keep;
        for (int i = 0; i < n; i++) J[(size_t)i * n + j] = (tmp[i] - fy[i]) / d;
    }
}

typedef struct {
    double rtol, atol; /* error per component <= atol + rtol * |y| */
    double hmax;       /* 0: no limit; as mk_dopri_ctl's */
    int autonomous;    /* 1: f does not depend on t (skips df/dt) */
    long accepted, rejected, evals, jacs, lus; /* added to by every call, as mk_dopri_ctl's */
    long maxsteps;     /* steps per call before giving up (-1); 0 = 10 million */
} mk_stiff_ctl;

/* Shampine's ode23s (Shampine and Reichelt 1997, the MATLAB ODE suite): a Rosenbrock (W) method, second order with a
   third-order error estimate, L-stable, one Jacobian and one LU per step, the last stage reused as the next first.
   Chosen over ROS3P because L-stability damps the fastest modes completely (ROS3P's stability function tends to 0.73 at
   infinity, so very stiff transients ring), which is what Robertson's kinetics and a relaxation oscillator need; ROS3P
   (Lang and Verwer 2001) is the better pick for semi-discretized parabolic PDEs at tight tolerance. jac may be NULL
   (forward differences). Integrates forward from t0 to t1 >= t0; *h the first step to try (0, or anything but a finite
   positive step, picks one, as in mk_dopri) and, after, the last. Steps as mk_dopri's: the size t can take, smaller
   after a rejection, a NaN step refused. 0; -1 if the step underflowed, maxsteps ran out, n is outside 0..46340, t0 or
   t1 is not finite, or t1 < t0; -2 no memory. The floor is 1e-14 of the span, so a fast start on a long span needs a
   call of its own (Robertson's from 0 to 1e11 in one call is -1; 0 to 1, then 1 to 1e11, matches Hairer and Wanner
   1996); W's LU refuses a pivot under 1e-13 of W's largest entry, which on Robertson's turns away steps past about
   1e10. */
static inline int mk_ode23s(mk_odefn f, mk_odejac jac, void *ctx, int n, double t0, double *y, double t1, double *h, mk_stiff_ctl *c)
{
    const double d = 1 / (2 + sqrt(2.0)), e32 = 6 + sqrt(2.0);
    /* the sizes mk_lu refuses; a NaN or infinite t0 or t1; backward */
    if (n < 0 || n > 46340 || !isfinite(t1 - t0) || t1 < t0) return -1;
    size_t nn = (size_t)n * n;
    double *mem = malloc(sizeof(double) * (2 * nn + 11 * (size_t)n));
    int *piv = malloc(sizeof(int) * (size_t)n);
    if (!mem || !piv) { free(mem), free(piv); return -2; }
    double *J = mem, *W = J + nn, *F0 = W + nn, *F1 = F0 + n, *F2 = F1 + n, *T = F2 + n, *k1 = T + n, *k2 = k1 + n, *k3 = k2 + n;
    double *yn = k3 + n, *z = yn + n, *tmp = z + n, *ys = tmp + n;
    /* the first guess as mk_dopri's */
    double t = t0, span = t1 - t0, hh = *h > 0 && isfinite(*h) ? *h : fmax(span * 1e-6, 16 * (nextafter(fabs(t0), INFINITY) - fabs(t0)));
    int ret = 0, needj = 1;
    if (span <= 0) { free(mem), free(piv); return 0; }
    f(t, y, F0, ctx), c->evals++;
    long cap = c->maxsteps > 0 ? c->maxsteps : 10000000, used = 0;
    double hrej = INFINITY; /* as mk_dopri's */
    while (t < t1) {
        if (used >= cap) { ret = -1; break; }
        used++;
        if (needj) {
            memcpy(ys, y, sizeof *y * (size_t)n);
            if (jac) jac(t, y, J, ctx);
            else mk_numjac_(f, ctx, n, t, ys, F0, J, tmp), c->evals += n;
            if (c->autonomous) for (int i = 0; i < n; i++) T[i] = 0;
            else {
                double dt = 1.4901161193847656e-8 * fmax(fabs(t), 1);
                f(t + dt, y, tmp, ctx), c->evals++;
                for (int i = 0; i < n; i++) T[i] = (tmp[i] - F0[i]) / dt;
            }
            c->jacs++, needj = 0;
        }
        if (c->hmax > 0 && hh > c->hmax) hh = c->hmax;
        double rem = t1 - t, hs = rem - hh < 1e-12 * span ? rem : hh; /* no sliver left before t1 */
        if (hs != rem) hs = (t + hs) - t; /* the step t can take, as mk_dopri's */
        if (hs >= hrej) hs = nextafter(t + hrej, t) - t; /* smaller than the step rejected, as mk_dopri's */
        if (hs < 1e-14 * span) { ret = -1; break; }
        for (size_t i = 0; i < nn; i++) W[i] = -hs * d * J[i];
        for (int i = 0; i < n; i++) W[(size_t)i * n + i] += 1;
        c->lus++;
        if (mk_lu(W, n, piv)) { hh = 0.5 * hs, hrej = hs, c->rejected++; continue; }
        for (int i = 0; i < n; i++) k1[i] = F0[i] + hs * d * T[i];
        mk_lu_solve(W, n, piv, k1);
        for (int i = 0; i < n; i++) z[i] = y[i] + 0.5 * hs * k1[i];
        f(t + 0.5 * hs, z, F1, ctx);
        for (int i = 0; i < n; i++) k2[i] = F1[i] - k1[i];
        mk_lu_solve(W, n, piv, k2);
#if !(defined(MK_MUTATE) && MK_MUTATE == 54) /* planted fault 54: k1 not added back to k2 */
        for (int i = 0; i < n; i++) k2[i] += k1[i];
#endif
        for (int i = 0; i < n; i++) yn[i] = y[i] + hs * k2[i];
        f(t + hs, yn, F2, ctx);
        for (int i = 0; i < n; i++) k3[i] = F2[i] - e32 * (k2[i] - F1[i]) - 2 * (k1[i] - F0[i]) + hs * d * T[i];
        mk_lu_solve(W, n, piv, k3);
        c->evals += 2;
        double err = 0;
        for (int i = 0; i < n; i++) {
            double ei = hs / 6 * (k1[i] - 2 * k2[i] + k3[i]), sc = c->atol + c->rtol * fmax(fabs(y[i]), fabs(yn[i]));
            const double q = fabs(ei) / sc;
            if (q > err || q != q) err = q; /* a NaN in any component stays NaN (fmax would drop it) */
        }
        if (!isfinite(err)) err = 1e10;
        if (err > 1) {
            c->rejected++, hrej = hs;
            hh = hs * fmax(0.1, 0.8 * pow(err, -1.0 / 3));
            continue;
        }
        c->accepted++, hrej = INFINITY;
        t = t1 - t - hs <= 0 ? t1 : t + hs;
        memcpy(y, yn, sizeof *y * (size_t)n), memcpy(F0, F2, sizeof *F0 * (size_t)n);
        hh = hs * fmin(5, 0.8 * pow(fmax(err, 1e-10), -1.0 / 3));
        needj = 1;
    }
    *h = hh;
    free(mem), free(piv);
    return ret;
}

typedef struct {
    double rtol, atol; /* Newton stops when every |change| <= 1e-2 (atol + rtol |y|) */
    int maxit;         /* Newton iterations, default 10 */
    long evals, jacs, iters;
} mk_newton_ctl;

/* One step of BDF1 (backward Euler) or BDF2 with a variable step: order 2 uses yprev from hprev before, the
   coefficients for the step ratio w = h / hprev, alpha (1+w)^2/(1+2w), (w^2)/(1+2w), beta (1+w)/(1+2w). Solved by
   Newton with one Jacobian and one LU (simplified Newton; mk_linalg.h); the predictor is linear extrapolation. y
   advances to t + h and yprev gets the old y. 0, -1 if Newton did not converge (y is left as it was: halve h and retry)
   or the step was refused (n outside 0..46340, t + h not finite), -2 no memory. */
static inline int mk_bdf_step(mk_odefn f, mk_odejac jac, void *ctx, int n, int order, double t, double *y, double *yprev,
                              double h, double hprev, mk_newton_ctl *c)
{
    if (n < 0 || n > 46340 || !isfinite(t + h)) return -1; /* the sizes mk_lu refuses; a NaN or infinite t or h */
    size_t nn = (size_t)n * n;
    double *mem = malloc(sizeof(double) * (nn + 5 * (size_t)n));
    int *piv = malloc(sizeof(int) * (size_t)n);
    if (!mem || !piv) { free(mem), free(piv); return -2; }
    double *J = mem, *Y = J + nn, *r = Y + n, *fy = r + n, *rhs = fy + n, *tmp = rhs + n;
    double w = order == 2 && hprev > 0 ? h / hprev : 0, a1 = 1, a2 = 0, beta = 1;
    if (order == 2 && hprev > 0) {
        a1 = (1 + w) * (1 + w) / (1 + 2 * w), a2 = w * w / (1 + 2 * w), beta = (1 + w) / (1 + 2 * w);
#if defined(MK_MUTATE) && MK_MUTATE == 55
        beta = 1; /* planted fault: backward Euler's beta in BDF2 */
#endif
    }
    for (int i = 0; i < n; i++) {
        rhs[i] = a1 * y[i] - a2 * (order == 2 && hprev > 0 ? yprev[i] : 0);
        Y[i] = order == 2 && hprev > 0 ? y[i] + w * (y[i] - yprev[i]) : y[i];
    }
    int maxit = c->maxit > 0 ? c->maxit : 10, ok = 0;
    for (int it = 0; it < maxit; it++) {
        f(t + h, Y, fy, ctx), c->evals++;
        if (it == 0) {
            if (jac) jac(t + h, Y, J, ctx);
            else mk_numjac_(f, ctx, n, t + h, Y, fy, J, tmp), c->evals += n;
            c->jacs++;
            for (size_t i = 0; i < nn; i++) J[i] *= -beta * h;
            for (int i = 0; i < n; i++) J[(size_t)i * n + i] += 1;
            if (mk_lu(J, n, piv)) break;
        }
        for (int i = 0; i < n; i++) r[i] = -(Y[i] - rhs[i] - beta * h * fy[i]);
        mk_lu_solve(J, n, piv, r);
        double big = 0;
        for (int i = 0; i < n; i++) {
            Y[i] += r[i];
            const double q = fabs(r[i]) / (c->atol + c->rtol * fabs(Y[i]));
            if (q > big || q != q) big = q; /* a NaN in any component fails the step (fmax would drop it) */
        }
        c->iters++;
        if (!isfinite(big)) break;
        if (big <= 1e-2) { ok = 1; break; }
    }
    if (ok) for (int i = 0; i < n; i++) yprev[i] = y[i], y[i] = Y[i];
    free(mem), free(piv);
    return ok ? 0 : -1;
}
#endif
