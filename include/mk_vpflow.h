/* mathkit: a volume-preserving deformation, the flow of an exactly divergence-free field (b-helper, 2026-11-15). The field is
   u = curl A, A's three components uniform cubic B-splines on a regular grid (control points at the nodes), so div u = 0 at
   every point, not only on the grid: the flow keeps every volume up to the RK4 step and the polygon of a moved surface. The
   coordinates are any (x0, x1, x2) whose volume element is constant: Cartesian, or (angle, r^2 / 2, length) about an axis,
   whose volume element is d(angle) d(r^2 / 2) ds. Any axis may be periodic (an angle); on the others, nodes past the grid
   count as 0, so the field fades out C2-smoothly over two cells and is 0 beyond. The time stepper is mk_ode.h's mk_rk4_step.
   Header only, C23.

   mk_vp: n[3] nodes, origin o[3] (node 0's coordinate), step h[3], per[3] (1: periodic, period n h), a[3] the components at
   the nodes, index (k n1 + j) n0 + i (NULL: a zero component). A grid is refused (NaN out) unless every n >= 1, every h is
   finite and > 0, every o finite, and n0 n1 n2 doubles fit in memory's address range.
   mk_vp_eval(P, x, A, J): A(x) and its Jacobian J[c][d] = dA_c / dx_d. mk_vp_vel(P, x, u): u = curl A.
   mk_vp_div(P, x, eps): div u by central differences of width eps (a check: 0 up to rounding and eps^2).
   mk_vp_advect(P, X, n, T, steps): each of n points (3 doubles each) carried for time T in `steps` RK4 steps; the largest
   move, or NaN (no point moved) for a refused grid, steps < 1 or T not finite.
   mk_vp_jacdet(P, x, T, steps, eps): the map's Jacobian determinant at x, central differences of width eps. */
#ifndef MK_VPFLOW_H
#define MK_VPFLOW_H
#include "mk_ode.h"
#include <math.h>
#include <stddef.h>
#include <limits.h>
#include <stdint.h>

typedef struct {
    int n[3];
    double o[3], h[3];
    int per[3];
    const double *a[3];
} mk_vp;

/* 0 if the grid can be used, -1 if not */
static inline int mk_vp_ok(const mk_vp *P)
{
    size_t m = 1;
    for (int ax = 0; ax < 3; ax++) {
        if (P->n[ax] < 1 || P->n[ax] > INT_MAX - 4 || !(P->h[ax] > 0) || !isfinite(P->h[ax]) || !isfinite(P->o[ax])) return -1;
        if ((size_t)P->n[ax] > SIZE_MAX / sizeof(double) / m) return -1;
        m *= (size_t)P->n[ax];
    }
    return 0;
}
/* uniform cubic B-spline weights and their derivatives at t in [0, 1): nodes i-1, i, i+1, i+2 */
static inline void mk_vp_bw_(double t, double w[4], double d[4])
{
    const double t2 = t * t, t3 = t2 * t, s = 1 - t;
    w[0] = s * s * s / 6, w[1] = (3 * t3 - 6 * t2 + 4) / 6, w[2] = (-3 * t3 + 3 * t2 + 3 * t + 1) / 6, w[3] = t3 / 6;
    d[0] = -s * s / 2, d[1] = (3 * t2 - 4 * t) / 2, d[2] = (-3 * t2 + 2 * t + 1) / 2, d[3] = t2 / 2;
}
/* a node index along an axis, or -1 off a non-periodic grid */
static inline int mk_vp_ix_(const mk_vp *P, int ax, int i)
{
    if (P->per[ax]) { const int n = P->n[ax]; i %= n; return i < 0 ? i + n : i; }
    return i < 0 || i >= P->n[ax] ? -1 : i;
}
static inline void mk_vp_eval(const mk_vp *P, const double *x, double A[3], double J[3][3])
{
    for (int c = 0; c < 3; c++) { A[c] = 0; for (int d = 0; d < 3; d++) J[c][d] = 0; }
    if (mk_vp_ok(P)) goto bad; /* a refused grid: NaN, as for mk_vp_advect */
    int i0[3];
    double w[3][4], dw[3][4];
    for (int ax = 0; ax < 3; ax++) {
        double g = (x[ax] - P->o[ax]) / P->h[ax];
        if (P->n[ax] < 1 || P->n[ax] > INT_MAX - 4 || !isfinite(g)) goto bad; /* a NaN or infinite x, o or h, a zero h, a bad axis */
        if (P->per[ax]) {
            if (fabs(g) > 0x1p30) g = fmod(g, P->n[ax]); /* whole periods out first, so the node index fits an int */
        } else if (g < -3 || g > P->n[ax] + 2) { /* two cells past the grid every node is missing: zero, unless x is not finite */
            if (!isfinite(x[0] + x[1] + x[2])) goto bad; /* one test, so the function stays small enough to inline */
            return;
        }
        const double f = floor(g);
        i0[ax] = (int)f, mk_vp_bw_(g - f, w[ax], dw[ax]);
        for (int k = 0; k < 4; k++) dw[ax][k] /= P->h[ax];
    }
    for (int kc = 0; kc < 4; kc++) {
        const int k = mk_vp_ix_(P, 2, i0[2] - 1 + kc);
        if (k < 0) continue;
        for (int jc = 0; jc < 4; jc++) {
            const int j = mk_vp_ix_(P, 1, i0[1] - 1 + jc);
            if (j < 0) continue;
            for (int ic = 0; ic < 4; ic++) {
                const int i = mk_vp_ix_(P, 0, i0[0] - 1 + ic);
                if (i < 0) continue;
                const size_t q = ((size_t)k * (size_t)P->n[1] + (size_t)j) * (size_t)P->n[0] + (size_t)i;
                const double b = w[0][ic] * w[1][jc] * w[2][kc], b0 = dw[0][ic] * w[1][jc] * w[2][kc], b1 = w[0][ic] * dw[1][jc] * w[2][kc],
                             b2 = w[0][ic] * w[1][jc] * dw[2][kc];
                for (int c = 0; c < 3; c++) {
                    if (!P->a[c]) continue;
                    const double v = P->a[c][q];
                    A[c] += b * v, J[c][0] += b0 * v, J[c][1] += b1 * v, J[c][2] += b2 * v;
                }
            }
        }
    }
    return;
bad:
    for (int c = 0; c < 3; c++) { A[c] = NAN; for (int d = 0; d < 3; d++) J[c][d] = NAN; }
}
static inline void mk_vp_vel(const mk_vp *P, const double *x, double u[3])
{
    double A[3], J[3][3];
    mk_vp_eval(P, x, A, J);
    u[0] = J[2][1] - J[1][2], u[1] = J[0][2] - J[2][0], u[2] = J[1][0] - J[0][1];
#if defined(MK_MUTATE) && MK_MUTATE == 131
    u[0] = J[2][1] + J[1][2]; /* not a curl, so not divergence-free */
#endif
}
static inline double mk_vp_div(const mk_vp *P, const double *x, double eps)
{
    double s = 0;
    for (int d = 0; d < 3; d++) {
        double xp[3] = {x[0], x[1], x[2]}, xm[3] = {x[0], x[1], x[2]}, up[3], um[3];
        xp[d] += eps, xm[d] -= eps;
        mk_vp_vel(P, xp, up), mk_vp_vel(P, xm, um);
        s += (up[d] - um[d]) / (2 * eps);
    }
    return s;
}
static inline void mk_vp_f_(double t, const double *y, double *dy, void *ctx) { (void)t, mk_vp_vel((const mk_vp *)ctx, y, dy); }
static inline double mk_vp_advect(const mk_vp *P, double *X, size_t n, double T, int steps)
{
    if (mk_vp_ok(P) || steps < 1 || !isfinite(T)) return NAN;
    const double dt = T / steps;
    double mx = 0, work[MK_RK4_WORK(3)];
    for (size_t p = 0; p < n; p++) {
        double *x = X + 3 * p;
        const double x0[3] = {x[0], x[1], x[2]};
        for (int s = 0; s < steps; s++) {
#if defined(MK_MUTATE) && MK_MUTATE == 132
            double u[3]; /* forward Euler, first order */
            mk_vp_vel(P, x, u);
            for (int d = 0; d < 3; d++) x[d] += dt * u[d];
            (void)work;
#else
            mk_rk4_step(mk_vp_f_, (void *)P, 3, s * dt, x, dt, work);
#endif
        }
        const double m = sqrt((x[0] - x0[0]) * (x[0] - x0[0]) + (x[1] - x0[1]) * (x[1] - x0[1]) + (x[2] - x0[2]) * (x[2] - x0[2]));
        if (m > mx || isnan(m)) mx = m;
    }
    return mx;
}
static inline double mk_vp_jacdet(const mk_vp *P, const double *x, double T, int steps, double eps)
{
    if (!(eps > 0)) return NAN;
    double M[3][3];
    for (int d = 0; d < 3; d++) {
        double p[3] = {x[0], x[1], x[2]}, m[3] = {x[0], x[1], x[2]};
        p[d] += eps, m[d] -= eps;
        if (isnan(mk_vp_advect(P, p, 1, T, steps)) || isnan(mk_vp_advect(P, m, 1, T, steps))) return NAN;
        for (int c = 0; c < 3; c++) M[c][d] = (p[c] - m[c]) / (2 * eps);
    }
    return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
           M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
}
#endif
