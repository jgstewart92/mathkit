/* mk_vpflow on known cases (b-helper's rows, 2026-11-15, plus the kit's): a ball (icosphere) carried by the flow of curl A,
   A = 0.02 (sin 2 pi y cos 2 pi z, sin 2 pi z cos 2 pi x, sin 2 pi x cos 2 pi y) at the nodes of a 33^3 grid on the unit
   cube: div u, the ball's volume, the map's Jacobian, RK4 against half the step; a rigid rotation from a quadratic potential
   (exact answer, fourth order); a uniform flow across a periodic axis's seam; the refusals. */
#include "../include/mk_vpflow.h"
#include "t.h"
#include <stdlib.h>
#include <string.h>

typedef struct { double *V; int nv, cap; int (*T)[3]; int nt; } msh;
static int addv(msh *m, double x, double y, double z)
{
    if (m->nv == m->cap) m->V = realloc(m->V, sizeof(double) * 3 * (size_t)(m->cap = m->cap ? 2 * m->cap : 1024));
    m->V[3 * m->nv] = x, m->V[3 * m->nv + 1] = y, m->V[3 * m->nv + 2] = z;
    return m->nv++;
}
static int mid(msh *m, int a, int b, int *key, int *val, int nk)
{
    const int lo = a < b ? a : b, hi = a < b ? b : a;
    int h = (int)(((unsigned)lo * 73856093u ^ (unsigned)hi * 19349663u) % (unsigned)nk);
    while (key[2 * h] >= 0) {
        if (key[2 * h] == lo && key[2 * h + 1] == hi) return val[h];
        h = (h + 1) % nk;
    }
    double x = (m->V[3 * a] + m->V[3 * b]) / 2, y = (m->V[3 * a + 1] + m->V[3 * b + 1]) / 2, z = (m->V[3 * a + 2] + m->V[3 * b + 2]) / 2;
    const double l = sqrt(x * x + y * y + z * z);
    key[2 * h] = lo, key[2 * h + 1] = hi;
    return val[h] = addv(m, x / l, y / l, z / l);
}
static void icosphere(msh *m, int sub)
{
    const double p = (1 + sqrt(5.0)) / 2, s = sqrt(1 + p * p);
    const double v[12][3] = {{-1, p, 0}, {1, p, 0}, {-1, -p, 0}, {1, -p, 0}, {0, -1, p}, {0, 1, p}, {0, -1, -p}, {0, 1, -p}, {p, 0, -1}, {p, 0, 1}, {-p, 0, -1}, {-p, 0, 1}};
    static const int f[20][3] = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                                 {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
    memset(m, 0, sizeof *m);
    for (int i = 0; i < 12; i++) addv(m, v[i][0] / s, v[i][1] / s, v[i][2] / s);
    m->T = malloc(sizeof *m->T * 20), m->nt = 20, memcpy(m->T, f, sizeof f);
    for (int r = 0; r < sub; r++) {
        const int nk = 8 * m->nt + 1;
        int *key = malloc(sizeof(int) * 2 * (size_t)nk), *val = malloc(sizeof(int) * (size_t)nk);
        for (int i = 0; i < 2 * nk; i++) key[i] = -1;
        int (*T2)[3] = malloc(sizeof *T2 * 4 * (size_t)m->nt);
        for (int t = 0; t < m->nt; t++) {
            const int a = m->T[t][0], b = m->T[t][1], c = m->T[t][2], ab = mid(m, a, b, key, val, nk), bc = mid(m, b, c, key, val, nk), ca = mid(m, c, a, key, val, nk);
            const int q[4][3] = {{a, ab, ca}, {b, bc, ab}, {c, ca, bc}, {ab, bc, ca}};
            memcpy(T2[4 * t], q, sizeof q);
        }
        free(m->T), m->T = T2, m->nt *= 4, free(key), free(val);
    }
}
static double vol(const double *V, int (*T)[3], int nt)
{
    double s = 0;
    for (int t = 0; t < nt; t++) {
        const double *a = V + 3 * T[t][0], *b = V + 3 * T[t][1], *c = V + 3 * T[t][2];
        s += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    return s / 6;
}
static double dist3(const double *a, const double *b) { return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])); }

int main(void)
{
    const int n = 33;
    const double h = 1.0 / (n - 1);
    double *a[3];
    for (int c = 0; c < 3; c++) a[c] = malloc(sizeof(double) * (size_t)n * n * n);
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                const double x = i * h, y = j * h, z = k * h;
                const size_t q = ((size_t)k * n + j) * n + i;
                a[0][q] = 0.02 * sin(2 * T_PI * y) * cos(2 * T_PI * z), a[1][q] = 0.02 * sin(2 * T_PI * z) * cos(2 * T_PI * x), a[2][q] = 0.02 * sin(2 * T_PI * x) * cos(2 * T_PI * y);
            }
    const mk_vp P = {{n, n, n}, {0, 0, 0}, {h, h, h}, {0, 0, 0}, {a[0], a[1], a[2]}};
    double dmax = 0, umax = 0;
    for (int r = 0; r < 200; r++) {
        const double x[3] = {0.5 + 0.35 * t_rand(), 0.5 + 0.35 * t_rand(), 0.5 + 0.35 * t_rand()};
        double u[3];
        mk_vp_vel(&P, x, u);
        umax = fmax(umax, sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2])), dmax = fmax(dmax, fabs(mk_vp_div(&P, x, 1e-5)));
    }
    CHECK(dmax < 1e-6 * umax / h, "div u at 200 random points: largest %.2e against |u| %.3f / h (bar 1e-6 of that)", dmax, umax);
    msh m;
    icosphere(&m, 5);
    for (int v = 0; v < m.nv; v++) for (int d = 0; d < 3; d++) m.V[3 * v + d] = 0.5 + 0.3 * m.V[3 * v + d];
    double *X = calloc(3 * (size_t)m.nv, sizeof(double));
    memcpy(X, m.V, sizeof(double) * 3 * (size_t)m.nv);
    const double v0 = vol(m.V, m.T, m.nt), mv = mk_vp_advect(&P, X, (size_t)m.nv, 1.0, 40), v1 = vol(X, m.T, m.nt);
    CHECK(mv > 0.1, "the flow moves the ball: largest vertex move %.3f (radius 0.3)", mv);
    CHECK(fabs(v1 - v0) / v0 < 2e-5, "ball volume %.6f -> %.6f (%.2e relative; bar 2e-5, the %d-face polyhedron's own error)", v0, v1, (v1 - v0) / v0, m.nt);
    double jmin = 1e9, jmax = -1e9;
    for (int r = 0; r < 60; r++) {
        const double x[3] = {0.5 + 0.25 * t_rand(), 0.5 + 0.25 * t_rand(), 0.5 + 0.25 * t_rand()};
        const double J = mk_vp_jacdet(&P, x, 1.0, 40, 1e-4);
        jmin = fmin(jmin, J), jmax = fmax(jmax, J);
    }
    CHECK(jmin > 1 - 1e-6 && jmax < 1 + 1e-6, "Jacobian determinant of the map at 60 points: %.9f .. %.9f (bar 1 +- 1e-6)", jmin, jmax);
    double y40[3] = {0.4, 0.55, 0.6}, y80[3] = {0.4, 0.55, 0.6};
    mk_vp_advect(&P, y40, 1, 1.0, 40), mk_vp_advect(&P, y80, 1, 1.0, 80);
    CHECK(dist3(y40, y80) < 1e-7, "a point after 40 steps against 80: %.2e apart (bar 1e-7, fourth order)", dist3(y40, y80));

    { /* a rigid rotation about x2 at rate w: A2 = -w (x0^2 + x1^2) / 2, whose B-spline coefficients are x^2 - h^2 / 3 (cubic
         B-splines reproduce quadratics exactly), on 17^3 nodes over [-1, 1]^3. A quarter turn has an exact answer. */
        const int nr = 17;
        const double hr = 2.0 / (nr - 1), w = 1.3;
        double *a2 = malloc(sizeof(double) * nr * nr * nr);
        for (int k = 0; k < nr; k++)
            for (int j = 0; j < nr; j++)
                for (int i = 0; i < nr; i++) {
                    const double x = -1 + i * hr, y = -1 + j * hr;
                    a2[((size_t)k * nr + j) * nr + i] = -w / 2 * (x * x - hr * hr / 3 + y * y - hr * hr / 3);
                }
        const mk_vp R = {{nr, nr, nr}, {-1, -1, -1}, {hr, hr, hr}, {0, 0, 0}, {NULL, NULL, a2}};
        const double T = T_PI / 2 / w, x0[3] = {0.31, -0.17, 0.05}, ex[3] = {0.17, 0.31, 0.05};
        double e[4];
        for (int s = 0; s < 4; s++) {
            double x[3] = {x0[0], x0[1], x0[2]};
            mk_vp_advect(&R, x, 1, T, 5 << s);
            e[s] = dist3(x, ex);
        }
        const double o1 = log2(e[0] / e[1]), o2 = log2(e[1] / e[2]), o3 = log2(e[2] / e[3]);
        CHECK(e[3] < 5e-8 && o1 > 3.8 && o2 > 3.8 && o3 > 3.8 && o1 < 4.2 && o2 < 4.2 && o3 < 4.2,
              "a quarter turn of a rigid rotation, 5, 10, 20, 40 steps: errors %.2g %.2g %.2g %.2g, orders %.2f %.2f %.2f (4)", e[0], e[1], e[2], e[3], o1, o2, o3);
        double u[3];
        const double xs[3] = {0.4, 0.2, -0.3};
        mk_vp_vel(&R, xs, u);
        CHECK(fabs(u[0] + w * 0.2) < 1e-13 && fabs(u[1] - w * 0.4) < 1e-13 && fabs(u[2]) < 1e-13,
              "the rotation's field is exact: u (%.15f %.15f %.15f), want (%.15f %.15f 0)", u[0], u[1], u[2], -w * 0.2, w * 0.4);
        free(a2);
    }
    { /* a periodic axis: x0 an angle of period 2 pi with 90 nodes; A2 depends on x0 only, so the field must agree at 0 and 2 pi */
        const int na = 90, nb = 9, nc = 9;
        double *b2 = malloc(sizeof(double) * na * nb * nc);
        for (int k = 0; k < nc; k++)
            for (int j = 0; j < nb; j++)
                for (int i = 0; i < na; i++) b2[((size_t)k * nb + j) * na + i] = 0.1 * sin(2 * T_PI * i / na) * j;
        const mk_vp Q = {{na, nb, nc}, {0, 0, 0}, {2 * T_PI / na, 0.1, 0.1}, {1, 0, 0}, {NULL, NULL, b2}};
        double u0[3], u1[3];
        const double x0[3] = {1e-9, 0.43, 0.41}, x1[3] = {2 * T_PI - 1e-9, 0.43, 0.41};
        mk_vp_vel(&Q, x0, u0), mk_vp_vel(&Q, x1, u1);
        CHECK(fabs(u0[0] - u1[0]) + fabs(u0[1] - u1[1]) < 1e-6 && fabs(u0[1]) > 1e-3, "periodic axis: the field at 0 and 2 pi agrees (%.6f %.6f against %.6f %.6f)", u0[0], u0[1], u1[0], u1[1]);
        /* a uniform flow along the periodic axis (A2 = c x1, linear, reproduced exactly) carries a point across the seam */
        for (int k = 0; k < nc; k++)
            for (int j = 0; j < nb; j++)
                for (int i = 0; i < na; i++) b2[((size_t)k * nb + j) * na + i] = 0.25 * (0.1 * j);
        double xp[3] = {2 * T_PI - 0.1, 0.43, 0.41}, up[3], uq[3];
        const double far[3] = {xp[0] + 2 * T_PI * 3, xp[1], xp[2]}, huge[3] = {1e300, 0.43, 0.41};
        mk_vp_vel(&Q, xp, up), mk_vp_vel(&Q, far, uq);
        mk_vp_advect(&Q, xp, 1, 2.0, 40);
        double uh[3];
        mk_vp_vel(&Q, huge, uh);
        CHECK(fabs(xp[0] - (2 * T_PI - 0.1 + 0.5)) < 1e-13 && fabs(xp[1] - 0.43) < 1e-15 && fabs(up[0] - uq[0]) < 1e-13 && isfinite(uh[0]) && fabs(uh[0] - 0.25) < 1e-12,
              "uniform flow 0.25 across the seam at 2 pi: moved %.15f (0.5); three periods on, u0 %.15f against %.15f; at x0 1e300 u0 %.15f",
              xp[0] - (2 * T_PI - 0.1), uq[0], up[0], uh[0]);
        free(b2);
    }
    { /* refusals: NaN out and no point moved */
        double y[3] = {0.4, 0.5, 0.6}, u[3];
        mk_vp B = P;
        B.n[1] = 0;
        const double r1 = mk_vp_advect(&B, y, 1, 1.0, 10);
        B = P, B.h[2] = 0;
        const double r2 = mk_vp_advect(&B, y, 1, 1.0, 10);
        const double r3 = mk_vp_advect(&P, y, 1, 1.0, 0), r4 = mk_vp_advect(&P, y, 1, NAN, 10);
        const double xn[3] = {0.4, NAN, 0.6};
        mk_vp_vel(&P, xn, u);
        B = P, B.n[0] = B.n[1] = B.n[2] = 1 << 21;
        CHECK(isnan(r1) && isnan(r2) && isnan(r3) && isnan(r4) && y[0] == 0.4 && y[1] == 0.5 && y[2] == 0.6 && isnan(u[0]) && mk_vp_ok(&B) == -1 && mk_vp_ok(&P) == 0,
              "refused with NaN and no point moved: an empty axis, h 0, steps 0, T NaN; x NaN gives u NaN; 2^63 nodes are refused");
    }
    for (int c = 0; c < 3; c++) free(a[c]);
    free(X), free(m.V), free(m.T);
    { /* the verification work's findings on 1.2.2 (Result 87; scratch/the verification work/mk-new/findings/demo.c) */
        mk_vp P2 = {{INT_MAX - 3, 1, 1}, {0, 0, 0}, {1, 1, 1}, {0, 0, 0}, {NULL, NULL, NULL}};
        double x2[3] = {(double)(INT_MAX - 3) + 2, 0.5, 0.5}, A2[3], J2[3][3];
        mk_vp_eval(&P2, x2, A2, J2);
        CHECK(mk_vp_ok(&P2) == -1 && isnan(A2[0]) && isnan(mk_vp_advect(&P2, x2, 1, 1.0, 1)), "V-A: INT_MAX - 3 nodes on an axis are refused (the node index would pass INT_MAX): NaN");
        double a8[8] = {1, 2, 3, 4, 5, 6, 7, 8};
        mk_vp P3 = {{2, 2, 2}, {0, 0, 0}, {-1, 1, 1}, {0, 0, 0}, {a8, a8, a8}}, P4 = {{2, 0, 2}, {0, 0, 0}, {1, 1, 1}, {0, 0, 0}, {a8, a8, a8}};
        double x3[3] = {-0.5, 0.5, 0.5}, x4[3] = {100, 0.5, 0.5}, A3[3], J3[3][3], u3[3], A4[3], J4[3][3];
        mk_vp_eval(&P3, x3, A3, J3), mk_vp_vel(&P3, x3, u3), mk_vp_eval(&P4, x4, A4, J4);
        CHECK(isnan(A3[0]) && isnan(u3[0]) && isnan(mk_vp_div(&P3, x3, 1e-3)) && isnan(A4[0]), "V-B: a grid mk_vp_ok refuses (h < 0; an empty axis behind an out-of-range one) gives NaN from eval, vel and div, not numbers");
    }
    return wrong;
}
