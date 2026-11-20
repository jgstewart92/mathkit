/* mathkit: minimization. Nelder-Mead (no derivatives, up to MK_NM_MAX parameters), Levenberg-Marquardt for sums of squares
   (Jacobian given or by forward differences), golden-section search on an interval, L-BFGS with a strong-Wolfe line search
   (gradients given), CMA-ES (no derivatives; for rugged, rotated or ill-conditioned problems). Each finds a low place, not
   necessarily the lowest: start from several places where that matters. Header only, C23; uses mk_linalg.h and mk_rng.h. */
#ifndef MK_OPT_H
#define MK_OPT_H
#include "mk_linalg.h"
#include "mk_rng.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef double (*mk_fn)(const double *p, void *ctx);

/* Nelder-Mead with the standard coefficients (reflect 1, expand 2, contract 1/2, shrink 1/2; Lagarias et al. 1998).
   p: start and result; step: the first simplex's size along each parameter. Stops when the spread of values is below tol on
   two iterations running (1.6.2: a simplex whose values tie across the minimum, x^2 from -3.5 with step 1 at {-0.5, 0.5}, stopped
   there; 1.6.1 skipped only the start's tie) or after iters iterations. Returns the value at p; *evals (if given) counts calls of fn. */
#define MK_NM_MAX 16
static inline double mk_nm(mk_fn fn, void *ctx, int n, double *p, const double *step, int iters, double tol, long *evals)
{
    if (n < 1 || n > MK_NM_MAX) { if (evals) *evals = 0; return NAN; }
    double s[MK_NM_MAX + 1][MK_NM_MAX], v[MK_NM_MAX + 1], c[MK_NM_MAX], r[MK_NM_MAX], e[MK_NM_MAX];
    long ne = 0;
    int flat = 0;
    for (int i = 0; i <= n; i++) {
        memcpy(s[i], p, sizeof *p * (size_t)n);
        if (i) s[i][i - 1] += step[i - 1];
        v[i] = fn(s[i], ctx), ne++;
    }
    for (int it = 0; it < iters; it++) {
        int hi = 0, lo = 0, nh;
        for (int i = 1; i <= n; i++) { if (v[i] > v[hi]) hi = i; if (v[i] < v[lo]) lo = i; }
        nh = lo;
        for (int i = 0; i <= n; i++) if (i != hi && v[i] > v[nh]) nh = i; /* second worst */
        if (!(v[hi] - v[lo] >= tol)) { if (flat) break; flat = 1; } else flat = 0; /* 1.6.2: under tol twice running (a tie across the minimum moves on) */
#if defined(MK_MUTATE) && MK_MUTATE == 41
        for (int j = 0; j < n; j++) { c[j] = 0; for (int i = 0; i <= n; i++) c[j] += s[i][j]; c[j] /= n + 1; } /* the worst point kept in the centroid */
#else
        for (int j = 0; j < n; j++) { c[j] = 0; for (int i = 0; i <= n; i++) if (i != hi) c[j] += s[i][j]; c[j] /= n; }
#endif
        for (int j = 0; j < n; j++) r[j] = 2 * c[j] - s[hi][j];
        double vr = fn(r, ctx);
        ne++;
        if (vr < v[lo]) {
            for (int j = 0; j < n; j++) e[j] = 3 * c[j] - 2 * s[hi][j];
            double ve = fn(e, ctx);
            ne++;
            if (ve < vr) memcpy(s[hi], e, sizeof *e * (size_t)n), v[hi] = ve;
            else memcpy(s[hi], r, sizeof *r * (size_t)n), v[hi] = vr;
            continue;
        }
        if (vr < v[nh]) { memcpy(s[hi], r, sizeof *r * (size_t)n), v[hi] = vr; continue; }
        int outside = vr < v[hi];
        for (int j = 0; j < n; j++) e[j] = outside ? c[j] + 0.5 * (r[j] - c[j]) : c[j] + 0.5 * (s[hi][j] - c[j]);
        double vk = fn(e, ctx);
        ne++;
        if (vk < (outside ? vr : v[hi])) { memcpy(s[hi], e, sizeof *e * (size_t)n), v[hi] = vk; continue; }
        for (int i = 0; i <= n; i++) /* shrink toward the best */
            if (i != lo) {
                for (int j = 0; j < n; j++) s[i][j] = s[lo][j] + 0.5 * (s[i][j] - s[lo][j]);
                v[i] = fn(s[i], ctx), ne++;
            }
    }
    int lo = 0;
    for (int i = 1; i <= n; i++) if (v[i] < v[lo]) lo = i;
    memcpy(p, s[lo], sizeof *p * (size_t)n);
    if (evals) *evals = ne;
    return v[lo];
}

/* residuals r[m] at p[n]; the Jacobian J[i * n + j] = d r_i / d p_j at p */
typedef void (*mk_resfn)(const double *p, double *r, void *ctx);
typedef void (*mk_jacfn)(const double *p, double *J, void *ctx);

typedef struct {
    int iters;      /* outer iterations, default 200 */
    double ftol;    /* stop when an accepted step lowers the cost by less than ftol of it, default 1e-12 */
    double lambda0; /* first damping, default 1e-3 */
} mk_lm_opts;

/* Levenberg-Marquardt: minimize sum r_i^2 over p (n <= MK_LM_MAX). Marquardt's scaling, damping lambda * (diag(JTJ) + tiny),
   divided by 3 on success and multiplied by 4 on failure (10 when the damped system is not positive definite).
   jac may be NULL: forward differences with step sqrt(eps) * max(|p_j|, 1). Returns the final sum of squares, or -1 for
   no memory or a size out of range. */
#define MK_LM_MAX 64
static inline double mk_lm(mk_resfn res, mk_jacfn jac, void *ctx, int m, int n, double *p, const mk_lm_opts *o)
{
    if (n < 1 || n > MK_LM_MAX || m < 1 || m > INT_MAX / (3 + n)) return -1; /* keeps the block under 2^34 bytes (the verification work, Results 74, 77) */
    int iters = o && o->iters > 0 ? o->iters : 200;
    double ftol = o && o->ftol > 0 ? o->ftol : 1e-12, lambda = o && o->lambda0 > 0 ? o->lambda0 : 1e-3;
    double *r = malloc(sizeof(double) * (3 * (size_t)m + (size_t)m * n));
    if (!r) return -1;
    double *r2 = r + m, *rh = r2 + m, *J = rh + m;
    double JTJ[MK_LM_MAX * MK_LM_MAX], A[MK_LM_MAX * MK_LM_MAX], g[MK_LM_MAX], d[MK_LM_MAX], pt[MK_LM_MAX], dk[MK_LM_MAX] = {0};
    res(p, r, ctx);
    double cost = 0;
    for (int i = 0; i < m; i++) cost += r[i] * r[i];
    for (int it = 0; it < iters && cost > 0; it++) {
        if (jac) jac(p, J, ctx);
        else
            for (int j = 0; j < n; j++) {
                double h = 1.4901161193847656e-8 * fmax(fabs(p[j]), 1), keep = p[j];
                p[j] = keep + h;
                res(p, rh, ctx);
                p[j] = keep;
                for (int i = 0; i < m; i++) J[(size_t)i * n + j] = (rh[i] - r[i]) / h;
            }
        double dmax = 0;
        for (int a = 0; a < n; a++) {
            double s = 0;
            for (int i = 0; i < m; i++) s += J[(size_t)i * n + a] * r[i];
#if defined(MK_MUTATE) && MK_MUTATE == 42
            g[a] = s;
#else
            g[a] = -s;
#endif
            for (int b = 0; b <= a; b++) {
                double t = 0;
                for (int i = 0; i < m; i++) t += J[(size_t)i * n + a] * J[(size_t)i * n + b];
                JTJ[a * n + b] = JTJ[b * n + a] = t;
            }
            dmax = fmax(dmax, JTJ[a * n + a]);
        }
        int accepted = 0;
        for (int tries = 0; tries < 30 && !accepted; tries++) {
            memcpy(A, JTJ, sizeof(double) * (size_t)(n * n));
            for (int a = 0; a < n; a++) dk[a] = fmax(dk[a], JTJ[a * n + a]), A[a * n + a] += lambda * (dk[a] + 1e-12 * dmax + 1e-300);
            memcpy(d, g, sizeof(double) * (size_t)n);
            if (mk_chol(A, n)) { lambda *= 10; continue; }
            mk_chol_solve(A, n, d);
            for (int a = 0; a < n; a++) pt[a] = p[a] + d[a];
            res(pt, r2, ctx);
            double c2 = 0;
            for (int i = 0; i < m; i++) c2 += r2[i] * r2[i];
            if (c2 < cost) {
                double drop = cost - c2;
                memcpy(p, pt, sizeof(double) * (size_t)n), memcpy(r, r2, sizeof(double) * (size_t)m);
                cost = c2, lambda = fmax(lambda / 3, 1e-30), accepted = drop > ftol * (cost + drop) || lambda > 1e-9 ? 1 : 2; /* a small drop under heavy damping is not convergence */
            } else lambda *= 4;
        }
        if (accepted != 1) break;
    }
    free(r);
    return cost;
}

/* golden-section search for a minimum of a unimodal f on [a, b], until the bracket is narrower than tol (absolute);
   returns the minimizer, *fmin its value */
typedef double (*mk_fn1)(double x, void *ctx);
static inline double mk_golden(mk_fn1 f, void *ctx, double a, double b, double tol, double *fmin)
{
    const double g = 0.61803398874989484820;
    double x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = f(x1, ctx), f2 = f(x2, ctx);
    for (int it = 0; it < 400 && fabs(b - a) > tol; it++) {
#if defined(MK_MUTATE) && MK_MUTATE == 43
        if (f1 > f2) {
#else
        if (f1 < f2) {
#endif
            b = x2, x2 = x1, f2 = f1, x1 = b - g * (b - a), f1 = f(x1, ctx);
        } else {
            a = x1, x1 = x2, f1 = f2, x2 = a + g * (b - a), f2 = f(x2, ctx);
        }
    }
    double x = f1 < f2 ? x1 : x2;
    if (fmin) *fmin = f1 < f2 ? f1 : f2;
    return x;
}

/* f at x, and its gradient into g */
typedef double (*mk_fgfn)(const double *x, double *g, void *ctx);

typedef struct {
    int m;        /* pairs kept, default 8 */
    int iters;    /* default 1000 */
    double gtol;  /* stop when max |g_i| <= gtol, default 1e-10 */
    double ftol;  /* stop when a step lowers f by less than ftol max(|f|, 1), default 1e-16 */
    int iters_used, status; /* out: status 0 gradient small, 1 no progress, 2 iteration limit, -1 failure */
    long evals;
} mk_lbfgs_opts;

/* phi(alpha) = f(x + alpha d) and its slope, into xn and gn */
static inline double mk_ls_eval_(mk_fgfn fg, void *ctx, int n, const double *x, const double *d, double a, double *xn, double *gn, double *slope, long *ev)
{
    for (int i = 0; i < n; i++) xn[i] = x[i] + a * d[i];
    double f = fg(xn, gn, ctx), s = 0;
    for (int i = 0; i < n; i++) s += gn[i] * d[i];
    (*ev)++;
    *slope = s;
    return f;
}
/* a step satisfying the strong Wolfe conditions (Nocedal and Wright, Algorithms 3.5 and 3.6; c1 1e-4, c2 0.9), with
   safeguarded cubic interpolation in the zoom. 0 found (xn, gn, *fn hold it), -1 not. */
static inline int mk_wolfe_(mk_fgfn fg, void *ctx, int n, const double *x, double f0, double s0, const double *d, double a1,
                            double *xn, double *gn, double *fn, long *ev)
{
    const double c1 = 1e-4, c2 = 0.9;
    double ap = 0, fp = f0, sp = s0, a = a1, lo = 0, flo = 0, slo = 0, hi = 0, fhi = 0, shi = 0;
    int zoom = 0;
    for (int i = 0; i < 40 && !zoom; i++) {
        double s, f = mk_ls_eval_(fg, ctx, n, x, d, a, xn, gn, &s, ev);
        if (!isfinite(f) || f > f0 + c1 * a * s0 || (i > 0 && f >= fp)) { lo = ap, flo = fp, slo = sp, hi = a, fhi = f, shi = s, zoom = 1; break; }
        if (fabs(s) <= -c2 * s0) { *fn = f; return 0; }
        if (s >= 0) { lo = a, flo = f, slo = s, hi = ap, fhi = fp, shi = sp, zoom = 1; break; }
        ap = a, fp = f, sp = s, a *= 2;
    }
    if (!zoom) return -1;
    for (int i = 0; i < 60; i++) {
        double lft = fmin(lo, hi), rgt = fmax(lo, hi), w = rgt - lft, at;
        if (w <= 1e-16 * fmax(1, rgt)) break;
        /* the cubic through (lo, flo, slo) and (hi, fhi, shi), if finite and well inside */
        double d1 = slo + shi - 3 * (flo - fhi) / (lo - hi), q = d1 * d1 - slo * shi;
        at = 0.5 * (lo + hi);
        if (isfinite(fhi) && q >= 0) {
            double d2 = (hi > lo ? 1 : -1) * sqrt(q), c = hi - (hi - lo) * (shi + d2 - d1) / (shi - slo + 2 * d2);
            if (isfinite(c) && c > lft + 0.1 * w && c < rgt - 0.1 * w) at = c;
        }
        double s, f = mk_ls_eval_(fg, ctx, n, x, d, at, xn, gn, &s, ev);
        if (!isfinite(f) || f > f0 + c1 * at * s0 || f >= flo) hi = at, fhi = f, shi = s;
        else {
            if (fabs(s) <= -c2 * s0) { *fn = f; return 0; }
            if (s * (hi - lo) >= 0) hi = lo, fhi = flo, shi = slo;
            lo = at, flo = f, slo = s;
        }
    }
    if (lo > 0 && flo < f0) { /* the best sufficient-decrease point found, short of the curvature condition */
        double s;
        *fn = mk_ls_eval_(fg, ctx, n, x, d, lo, xn, gn, &s, ev);
        return 0;
    }
    return -1;
}

/* L-BFGS (Nocedal 1980, the two-loop recursion; the first scaling s.y / y.y). x: start and result. Returns f at x. */
static inline double mk_lbfgs(mk_fgfn fg, void *ctx, int n, double *x, mk_lbfgs_opts *o)
{
    int m = o && o->m > 0 ? o->m : 8, iters = o && o->iters > 0 ? o->iters : 1000;
    double gtol = o && o->gtol > 0 ? o->gtol : 1e-10, ftol = o && o->ftol > 0 ? o->ftol : 1e-16;
    size_t per = 2 * (size_t)m + 5; /* in size_t: 2 m + 5 overflows int from m = 2^30 - 2 */
    if (n < 0 || (n > 0 && per > (SIZE_MAX / sizeof(double) - 2 * (size_t)m) / (size_t)n)) { if (o) o->status = -1, o->iters_used = 0, o->evals = 0; return NAN; } /* 1.6.2: no count left from an earlier run */
    double *w = malloc(sizeof(double) * (per * (size_t)n + 2 * (size_t)m));
    if (!w) { if (o) o->status = -1, o->iters_used = 0, o->evals = 0; return NAN; }
    double *S = w, *Y = S + (size_t)m * n, *g = Y + (size_t)m * n, *d = g + n, *xn = d + n, *gn = xn + n, *q = gn + n;
    double *rho = q + n, *al = rho + m;
    long ev = 1;
    double f = fg(x, g, ctx);
    int used = 0, head = 0, status = 2, it;
    for (it = 0; it < iters; it++) {
        double gmax = 0;
        for (int i = 0; i < n; i++) gmax = fabs(g[i]) > gmax || isnan(g[i]) ? fabs(g[i]) : gmax; /* a NaN stays: not small */
        if (gmax <= gtol) { status = 0; break; }
        for (int i = 0; i < n; i++) q[i] = g[i];
        for (int j = 0; j < used; j++) {
            int s = (head - 1 - j + m) % m;
            double a = 0;
            for (int i = 0; i < n; i++) a += S[(size_t)s * n + i] * q[i];
            al[s] = a *= rho[s];
            for (int i = 0; i < n; i++) q[i] -= a * Y[(size_t)s * n + i];
        }
        double gam = 1;
        if (used) {
            int s = (head - 1 + m) % m;
            double sy = 0, yy = 0, ym = 0; /* y scaled by a power of two (exact): y.y underflows once |y| < 1e-154 */
            for (int i = 0; i < n; i++) ym = fmax(ym, fabs(Y[(size_t)s * n + i]));
            int e = ilogb(ym);
            for (int i = 0; i < n; i++) {
                double yi = ldexp(Y[(size_t)s * n + i], -e);
                sy += S[(size_t)s * n + i] * yi, yy += yi * yi;
            }
            gam = ldexp(sy / yy, -e);
        }
        for (int i = 0; i < n; i++) q[i] *= gam;
        for (int j = used - 1; j >= 0; j--) {
            int s = (head - 1 - j + m) % m;
            double b = 0;
            for (int i = 0; i < n; i++) b += Y[(size_t)s * n + i] * q[i];
            b *= rho[s];
#if defined(MK_MUTATE) && MK_MUTATE == 44
            for (int i = 0; i < n; i++) q[i] -= S[(size_t)s * n + i] * (al[s] - b);
#else
            for (int i = 0; i < n; i++) q[i] += S[(size_t)s * n + i] * (al[s] - b);
#endif
        }
        double s0 = 0;
        for (int i = 0; i < n; i++) d[i] = -q[i], s0 += d[i] * g[i];
        if (!(s0 < 0) || !used) { /* not a descent direction: drop the memory, go downhill; d = -g / |g|, a first step
                                      min(1, |g|) as before, but |g|^2 does not overflow for |g| up to the largest double */
            used = 0, s0 = 0;
            double gs = 0;
            for (int i = 0; i < n; i++) gs += (g[i] / gmax) * (g[i] / gmax);
            double gnorm = gmax * sqrt(gs);
            for (int i = 0; i < n; i++) d[i] = -(g[i] / gmax) / sqrt(gs), s0 += d[i] * g[i];
            if (!(s0 < 0)) s0 = -gnorm; /* an infinite or NaN g: the search fails at once (-1) */
        }
        double a1 = used ? 1 : fmin(1, -s0), fnew;
        if (mk_wolfe_(fg, ctx, n, x, f, s0, d, a1, xn, gn, &fnew, &ev)) {
            if (used) { used = 0; continue; } /* retry once from steepest descent */
            status = -1;
            break;
        }
        double sy = 0, ss = 0, yy = 0, sn = 0, sm = 0, ym = 0; /* the guard on s, y scaled by powers of two (exact), so ss yy can't underflow or overflow */
        for (int i = 0; i < n; i++) sm = fmax(sm, fabs(xn[i] - x[i])), ym = fmax(ym, fabs(gn[i] - g[i]));
        int es = sm > 0 ? ilogb(sm) : 0, ey = ym > 0 ? ilogb(ym) : 0;
        for (int i = 0; i < n; i++) {
            double si = xn[i] - x[i], yi = gn[i] - g[i], a = ldexp(si, -es), b = ldexp(yi, -ey);
            sy += si * yi, sn += a * b, ss += a * a, yy += b * b;
        }
        if (sn > 1e-12 * sqrt(ss * yy) && sy >= 0x1p-1022) { /* body-motion-support minimize.c's guard; a rejected pair leaves the kept ones alone */
            for (int i = 0; i < n; i++) S[(size_t)head * n + i] = xn[i] - x[i], Y[(size_t)head * n + i] = gn[i] - g[i];
            rho[head] = 1 / sy, head = (head + 1) % m, used = used < m ? used + 1 : m;
        }
        double drop = f - fnew;
        for (int i = 0; i < n; i++) x[i] = xn[i], g[i] = gn[i];
        f = fnew;
        if (drop <= ftol * fmax(fabs(f), 1)) {
            double gm = 0;
            for (int i = 0; i < n; i++) gm = fabs(g[i]) > gm || isnan(g[i]) ? fabs(g[i]) : gm;
            status = gm <= gtol ? 0 : 1;
            it++;
            break;
        }
    }
    if (o) o->iters_used = it, o->status = status, o->evals = ev;
    free(w);
    return f;
}

typedef struct {
    double sigma0;   /* first step size (about a quarter of the search range) */
    int lambda;      /* population, 0 = 4 + floor(3 ln n); 1 is refused (status -1) */
    long maxevals;   /* default 1000 n^2 */
    double ftarget;  /* stop at f <= ftarget; 0 = no target (to stop at f = 0, give 0x1p-1074) */
    double xtol;     /* stop when sigma times the largest axis is under xtol, default 1e-12 */
    uint64_t seed;
    const double *scale; /* per-coordinate scale (as glance-core's gaittrain.c), NULL = 1: the first covariance is diag(scale^2) */
    long evals;      /* out */
    int status;      /* out: 0 target, 1 xtol, 2 budget, -1 failure */
} mk_cmaes_opts;

/* CMA-ES, (mu/mu_w, lambda) with rank-one and rank-mu updates and cumulative step-size control (Hansen 2016, "The CMA
   Evolution Strategy: A Tutorial", its default constants). n <= MK_CMA_MAX. x: start (the first mean) and the best found.
   For a multimodal f, restart with a doubled lambda (IPOP) and keep the best. Returns f at x. */
#define MK_CMA_MAX 64
static inline double mk_cmaes(mk_fn fn, void *ctx, int n, double *x, mk_cmaes_opts *o)
{
    if (n < 1 || n > MK_CMA_MAX || !o || o->lambda == 1) { if (o) o->status = -1, o->evals = 0; return NAN; } /* 1.6.0: lambda 1 has no parent (mu 0); 1.6.2: evals 0, not an earlier run's */
    int lam = o->lambda > 0 ? o->lambda : 4 + (int)floor(3 * log(n)), mu = lam / 2;
    long maxev = o->maxevals > 0 ? o->maxevals : 1000L * n * n;
    double ftarget = o->ftarget != 0 ? o->ftarget : -INFINITY, xtol = o->xtol > 0 ? o->xtol : 1e-12;
    size_t nn = (size_t)n * n;
    double *mem = malloc(sizeof(double) * (4 * nn + (size_t)lam * (2 * n + 1) + 8 * (size_t)n + (size_t)mu + 1));
    int *idx = malloc(sizeof(int) * (size_t)lam);
    if (!mem || !idx) { free(mem), free(idx); o->status = -1, o->evals = 0; return NAN; }
    double *C = mem, *B = C + nn, *A = B + nn, *Vv = A + nn, *Z = Vv + nn, *Yk = Z + (size_t)lam * n, *fk = Yk + (size_t)lam * n;
    double *m = fk + lam, *mold = m + n, *ps = mold + n, *pc = ps + n, *D = pc + n, *yw = D + n, *t = yw + n, *best = t + n, *wt = best + n;
    double wsum = 0, w2 = 0;
    for (int i = 0; i < mu; i++) wt[i] = log((lam + 1) / 2.0) - log(i + 1.0), wsum += wt[i];
    for (int i = 0; i < mu; i++) wt[i] /= wsum, w2 += wt[i] * wt[i];
    double mueff = 1 / w2;
    double cs = (mueff + 2) / (n + mueff + 5), ds = 1 + 2 * fmax(0, sqrt((mueff - 1) / (n + 1)) - 1) + cs;
    double cc = (4 + mueff / n) / (n + 4 + 2 * mueff / n), c1 = 2 / ((n + 1.3) * (n + 1.3) + mueff);
    double cmu = fmin(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((n + 2) * (n + 2) + mueff));
#if defined(MK_MUTATE) && MK_MUTATE == 45
    c1 = cmu = 0;
#endif
    double chin = sqrt(n) * (1 - 1.0 / (4 * n) + 1.0 / (21.0 * n * n)), sigma = o->sigma0 > 0 ? o->sigma0 : 0.3;
    for (int i = 0; i < n; i++) m[i] = x[i], ps[i] = pc[i] = 0, D[i] = 1;
    for (size_t i = 0; i < nn; i++) C[i] = B[i] = 0;
    for (int i = 0; i < n; i++) {
        double s = o->scale && o->scale[i] > 0 ? o->scale[i] : 1;
        C[i * n + i] = s * s, B[i * n + i] = 1, D[i] = s;
    }
    mk_rng rng;
    mk_rng_seed(&rng, o->seed, 0x434d41ULL);
    double fbest = fn(m, ctx);
    long ev = 1;
    for (int i = 0; i < n; i++) best[i] = m[i];
    int status = 2;
    for (long gen = 0; ev < maxev; gen++) {
        if (fbest <= ftarget) { status = 0; break; }
        for (int k = 0; k < lam; k++) {
            double *z = Z + (size_t)k * n, *y = Yk + (size_t)k * n;
            for (int i = 0; i < n; i++) z[i] = mk_rng_gauss(&rng);
            for (int i = 0; i < n; i++) {
                double s = 0;
                for (int j = 0; j < n; j++) s += B[i * n + j] * D[j] * z[j];
                y[i] = s, t[i] = m[i] + sigma * s;
            }
            fk[k] = fn(t, ctx), ev++, idx[k] = k;
            if (fk[k] < fbest || (isnan(fbest) && !isnan(fk[k]))) { fbest = fk[k]; for (int i = 0; i < n; i++) best[i] = t[i]; } /* a NaN start is replaced */
        }
        for (int i = 1; i < lam; i++) /* rank by f, a NaN as the worst value */
            for (int j = i; j > 0 && (fk[idx[j]] < fk[idx[j - 1]] || (isnan(fk[idx[j - 1]]) && !isnan(fk[idx[j]]))); j--) { int s = idx[j]; idx[j] = idx[j - 1], idx[j - 1] = s; }
        for (int i = 0; i < n; i++) {
            mold[i] = m[i], yw[i] = 0;
            for (int k = 0; k < mu; k++) yw[i] += wt[k] * Yk[(size_t)idx[k] * n + i];
            m[i] += sigma * yw[i];
        }
        /* C^-1/2 yw = B D^-1 B^T yw */
        for (int j = 0; j < n; j++) {
            double s = 0;
            for (int i = 0; i < n; i++) s += B[i * n + j] * yw[i];
            t[j] = s / D[j];
        }
        double psn = 0, k1 = sqrt(cs * (2 - cs) * mueff);
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (int j = 0; j < n; j++) s += B[i * n + j] * t[j];
            ps[i] = (1 - cs) * ps[i] + k1 * s, psn += ps[i] * ps[i];
        }
        psn = sqrt(psn);
        int hs = psn / sqrt(1 - pow(1 - cs, 2.0 * (gen + 1))) < (1.4 + 2.0 / (n + 1)) * chin;
        double k2 = sqrt(cc * (2 - cc) * mueff);
        for (int i = 0; i < n; i++) pc[i] = (1 - cc) * pc[i] + (hs ? k2 * yw[i] : 0);
        double keep = 1 - c1 - cmu + (hs ? 0 : c1 * cc * (2 - cc));
        for (int i = 0; i < n; i++)
            for (int j = 0; j <= i; j++) {
                double r = 0;
                for (int k = 0; k < mu; k++) r += wt[k] * Yk[(size_t)idx[k] * n + i] * Yk[(size_t)idx[k] * n + j];
                C[i * n + j] = C[j * n + i] = keep * C[i * n + j] + c1 * pc[i] * pc[j] + cmu * r;
            }
        sigma *= exp(fmin(1, (cs / ds) * (psn / chin - 1)));
        memcpy(A, C, sizeof(double) * nn);
        if (mk_sym_jacobi(A, n, Vv, D) < 0) { status = -1; break; }
        memcpy(B, Vv, sizeof(double) * nn);
        double dmax = 0;
        for (int i = 0; i < n; i++) D[i] = sqrt(fmax(D[i], 1e-300)), dmax = fmax(dmax, D[i]);
        if (sigma * dmax < xtol) { status = 1; break; }
        if (!isfinite(sigma)) { status = -1; break; }
    }
    if (fbest <= ftarget) status = 0;
    for (int i = 0; i < n; i++) x[i] = best[i];
    o->evals = ev, o->status = status;
    free(mem), free(idx);
    return fbest;
}
#endif
