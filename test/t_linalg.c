/* mk_linalg.h: LU, Cholesky, least squares, the symmetric value problem and the tridiagonal solve on systems with known
   answers (hand-solved 3x3s, a built Q diag Q^T, an exact line fit, a discrete Poisson problem), residuals in long
   double, cases found by fuzzing (stored in test/data), size guards and non-finite input. */
#include "../include/mk_linalg.h"
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include "t.h"

static double resid(const double *A, int n, const double *x, const double *b)
{
    double r = 0;
    for (int i = 0; i < n; i++) {
        double s = -b[i];
        for (int j = 0; j < n; j++) s += A[i * n + j] * x[j];
        r = fmax(r, fabs(s));
    }
    return r;
}

/* a least-squares case from test/data: "op m n", then A row-major and b, as hex floats */
static int load_ls(const char *path, int *m, int *n, double **A, double **b)
{
    FILE *f = fopen(path, "r");
    int op;
    if (!f || fscanf(f, "%d %d %d", &op, m, n) != 3) { if (f) fclose(f); return -1; }
    *A = malloc(sizeof(double) * (size_t)(*m * *n)), *b = malloc(sizeof(double) * (size_t)*m);
    for (int i = 0; i < *m * *n; i++) if (fscanf(f, "%la", &(*A)[i]) != 1) { fclose(f); return -1; }
    for (int i = 0; i < *m; i++) if (fscanf(f, "%la", &(*b)[i]) != 1) { fclose(f); return -1; }
    fclose(f);
    return 0;
}
/* mk_lstsq on a file's problem: the returned residual against x's own (long double), and the largest |x| */
static void ls_case(const char *path, double *ret, double *own, double *xmax)
{
    int m, n;
    double *A, *b;
    *ret = *own = *xmax = NAN;
    if (load_ls(path, &m, &n, &A, &b)) return;
    double *a2 = malloc(sizeof(double) * (size_t)(m * n)), *b2 = malloc(sizeof(double) * (size_t)m), x[64];
    memcpy(a2, A, sizeof(double) * (size_t)(m * n)), memcpy(b2, b, sizeof(double) * (size_t)m);
    *ret = mk_lstsq(a2, m, n, b2, x);
    long double r = 0;
    *xmax = 0;
    for (int j = 0; j < n; j++) *xmax = fmax(*xmax, fabs(x[j]));
    for (int i = 0; i < m; i++) {
        long double s = -(long double)b[i];
        for (int j = 0; j < n; j++) s += (long double)A[i * n + j] * x[j];
        r += s * s;
    }
    *own = (double)sqrtl(r);
    free(A), free(b), free(a2), free(b2);
}
int main(void)
{
    /* a zero first pivot: solvable only with row exchanges */
    double A[9] = {0, 2, 1, 1, 1, 1, 2, 1, 3}, LU[9], b[3] = {5, 4, 7}, x[3];
    int piv[3];
    memcpy(LU, A, sizeof A);
    int r = mk_lu(LU, 3, piv);
    memcpy(x, b, sizeof b);
    if (!r) mk_lu_solve(LU, 3, piv, x);
    CHECK(!r && fabs(x[0] - 1) < 1e-14 && fabs(x[1] - 2) < 1e-14 && fabs(x[2] - 1) < 1e-14, "LU with a zero first pivot solves to (1, 2, 1)");
    CHECK(!r && fabs(mk_lu_det(LU, 3, piv) - (-3)) < 1e-13, "its determinant is -3 (%.17g)", r ? NAN : mk_lu_det(LU, 3, piv));
    double S[9] = {1, 2, 3, 2, 4, 6, 1, 0, 1};
    CHECK(mk_lu(S, 3, piv) == -1, "a singular matrix is refused");

    enum { N = 8 };
    double B[N * N], M[N * N], L[N * N], c[N], y[N];
    for (int i = 0; i < N * N; i++) B[i] = t_rand();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double s = i == j ? N : 0;
            for (int k = 0; k < N; k++) s += B[k * N + i] * B[k * N + j];
            M[i * N + j] = s;
        }
    for (int i = 0; i < N; i++) c[i] = t_rand();
    memcpy(L, M, sizeof M);
    int rc = mk_chol(L, N);
    memcpy(y, c, sizeof c);
    if (!rc) mk_chol_solve(L, N, y);
    CHECK(!rc && resid(M, N, y, c) < 1e-12, "Cholesky of B^T B + 8 I solves to residual %.2g", resid(M, N, y, c));
    double neg[4] = {1, 2, 2, 1};
    CHECK(mk_chol(neg, 2) == -1, "an indefinite matrix is refused by Cholesky");
    int pv[N];
    memcpy(L, B, sizeof B);
    memcpy(y, c, sizeof c);
    int rl = mk_lu(L, N, pv);
    if (!rl) mk_lu_solve(L, N, pv, y);
    CHECK(!rl && resid(B, N, y, c) < 1e-12, "LU of a random 8 x 8 solves to residual %.2g", resid(B, N, y, c));

    /* least squares: y = 2 + 300 x on 20 points (columns of very different size), exact; and a doubled column */
    enum { P = 20 };
    double D[P * 2], e[P], z[2], D3[P * 3], e3[P], z3[3];
    for (int i = 0; i < P; i++) {
        double t = i * 0.01;
        D[i * 2] = 1, D[i * 2 + 1] = t, e[i] = 2 + 300 * t;
        D3[i * 3] = 1, D3[i * 3 + 1] = t, D3[i * 3 + 2] = 2 * t, e3[i] = 2 + 300 * t;
    }
    double rn = mk_lstsq(D, P, 2, e, z);
    CHECK(fabs(z[0] - 2) < 1e-11 && fabs(z[1] - 300) < 1e-9 && rn < 1e-10, "least squares recovers 2 + 300 x (%.15g, %.15g, residual %.2g)", z[0], z[1], rn);
    double rn3 = mk_lstsq(D3, P, 3, e3, z3);
    double fit = 0;
    for (int i = 0; i < P; i++) fit = fmax(fit, fabs(z3[0] + (z3[1] + 2 * z3[2]) * i * 0.01 - (2 + 300 * i * 0.01)));
    CHECK(isfinite(z3[0] + z3[1] + z3[2]) && fit < 1e-9 && rn3 < 1e-9, "a repeated column still solves: fit error %.2g", fit);

    /* symmetric: Q diag(-2, 0.5, 1, 3, 7, 11) Q^T with Q from Givens rotations */
    enum { K = 6 };
    double lam[K] = {7, -2, 11, 0.5, 3, 1}, Q[K * K], As[K * K], V[K * K], val[K];
    for (int i = 0; i < K * K; i++) Q[i] = i % (K + 1) == 0;
    for (int g = 0; g < 40; g++) {
        int p = g % K, q = (p + 1 + (g / K) % (K - 1)) % K;
        if (p == q) continue;
        double th = t_rand() * T_PI, cs = cos(th), sn = sin(th);
        for (int k = 0; k < K; k++) {
            double a = Q[k * K + p], bq = Q[k * K + q];
            Q[k * K + p] = cs * a - sn * bq, Q[k * K + q] = sn * a + cs * bq;
        }
    }
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++) {
            double s = 0;
            for (int k = 0; k < K; k++) s += Q[i * K + k] * lam[k] * Q[j * K + k];
            As[i * K + j] = s;
        }
    double Ac[K * K];
    memcpy(Ac, As, sizeof As);
    int sw = mk_sym_jacobi(Ac, K, V, val);
    double want[K] = {-2, 0.5, 1, 3, 7, 11}, ev = 0, av = 0;
    for (int k = 0; k < K; k++) {
        ev = fmax(ev, fabs(val[k] - want[k]));
        for (int i = 0; i < K; i++) {
            double s = 0;
            for (int j = 0; j < K; j++) s += As[i * K + j] * V[j * K + k];
            av = fmax(av, fabs(s - val[k] * V[i * K + k]));
        }
    }
    CHECK(sw >= 0 && ev < 1e-12 && av < 1e-12, "symmetric 6 x 6: values -2 .. 11 ascending (worst %.2g), A v = value v (worst %.2g), %d sweeps", ev, av, sw);

    double A3[3][3] = {{0, 1, 2}, {3, 0, 1}, {1, 1, 0}}, b3[3] = {8, 6, 3}, x3[3];
    int r3 = mk_solve3(A3, b3, x3);
    CHECK(!r3 && fabs(x3[0] - 1) < 1e-14 && fabs(x3[1] - 2) < 1e-14 && fabs(x3[2] - 3) < 1e-14, "3x3 with a zero first pivot solves to (1, 2, 3)");
    CHECK(fabs(mk_det3(A3) - 7) < 1e-14, "3x3 determinant 7 (%.17g)", mk_det3(A3));
    double T3[3][3] = {{4, 1, 0.5}, {1, 3, 0.2}, {0.5, 0.2, 1}}, v3[3], w3[3][3], worst = 0;
    mk_sym3(T3, v3, w3);
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < 3; i++) worst = fmax(worst, fabs(T3[i][0] * w3[k][0] + T3[i][1] * w3[k][1] + T3[i][2] * w3[k][2] - v3[k] * w3[k][i]));
    double hand = w3[2][0] * (w3[0][1] * w3[1][2] - w3[0][2] * w3[1][1]) + w3[2][1] * (w3[0][2] * w3[1][0] - w3[0][0] * w3[1][2]) + w3[2][2] * (w3[0][0] * w3[1][1] - w3[0][1] * w3[1][0]);
    CHECK(worst < 1e-13 && v3[0] <= v3[1] && v3[1] <= v3[2] && fabs(hand - 1) < 1e-13, "symmetric 3x3: A v = value v (%.2g), ascending, right-handed (%.15g)", worst, hand);
    /* ---- cases found by fuzzing: rank near the span, scaling, values near overflow, infinite entries ---- */
    double lr, lo, lx;
    ls_case("test/data/lstsq_keeps_dependent.txt", &lr, &lo, &lx);
    CHECK(fabs(lr - lo) <= 1e-12 * lo && lx < 1e6, "least squares: a column 2.3e-14 from the span is dropped (41 x 37): residual %.16g, x's own %.16g, max |x| %.2g", lr, lo, lx);
    ls_case("test/data/lstsq_drops_independent.txt", &lr, &lo, &lx);
    CHECK(fabs(lr - lo) <= 1e-12 * lo && lr <= 2.9039116231485418 * (1 + 1e-12), "least squares: a column 1.1e-5 from the span is kept (34 x 24): residual %.16g, x's own %.16g (minimum 2.9039116231485418)", lr, lo);
    double scs[10] = {1e150, 1e154, 1e155, 1e200, 1e300, 1e-150, 1e-160, 1e-170, 1e-200, 1e-300}, wsc = 0;
    for (int t = 0; t < 10; t++) {
        double a[16] = {0}, bb[4] = {1, 1, 1, 1}, xx[4];
        for (int i = 0; i < 4; i++) a[i * 4 + i] = 1;
        a[0] = scs[t];
        double r = mk_lstsq(a, 4, 4, bb, xx);
        wsc = fmax(wsc, fmax(fabs(xx[0] * scs[t] - 1), r == 0 ? 0 : 1));
    }
    CHECK(wsc < 1e-15, "least squares: diag(c, 1, 1, 1) x = 1 for c 1e-300 .. 1e300: x0 = 1/c and residual 0 (worst %.2g)", wsc);
    double jb[4] = {0, 1e155, 1e155, 0}, jv[4], jw[2];
    int jr = mk_sym_jacobi(jb, 2, jv, jw);
    CHECK(jr >= 0 && fabs(jw[0] + 1e155) < 1e140 && fabs(jw[1] - 1e155) < 1e140, "sym_jacobi: [[0, 1e155], [1e155, 0]] has values -1e155, 1e155 (got %g, %g, return %d)", jw[0], jw[1], jr);
    double js[4] = {0, 1e-170, 1e-170, 0};
    jr = mk_sym_jacobi(js, 2, jv, jw);
    CHECK(jr >= 0 && fabs(jw[0] + 1e-170) < 1e-185 && fabs(jw[1] - 1e-170) < 1e-185, "sym_jacobi: [[0, 1e-170], [1e-170, 0]] has values -1e-170, 1e-170 (got %g, %g)", jw[0], jw[1]);
    double ji[4] = {1, INFINITY, INFINITY, 1}, j3[3][3] = {{1, 0, 0}, {0, INFINITY, 0}, {0, 0, 1}}, v3i[3], w3i[3][3];
    CHECK(mk_sym_jacobi(ji, 2, jv, jw) == -1 && mk_sym3(j3, v3i, w3i) == -1, "an infinite entry is refused (-1) by mk_sym_jacobi and mk_sym3");
    enum { DN = 200 };
    double *dg = calloc(DN * DN, sizeof(double));
    int dp[DN];
    double det1, det2;
    for (int i = 0; i < DN; i++) dg[i * DN + i] = i < DN / 2 ? 1e6 : 1e-6;
    mk_lu(dg, DN, dp), det1 = mk_lu_det(dg, DN, dp);
    for (int i = 0; i < DN * DN; i++) dg[i] = 0;
    for (int i = 0; i < DN; i++) dg[i * DN + i] = i < DN / 2 ? 1e-6 : 1e6;
    mk_lu(dg, DN, dp), det2 = mk_lu_det(dg, DN, dp);
    free(dg);
    CHECK(fabs(det1 - 1) < 1e-12 && fabs(det2 - 1) < 1e-12, "LU determinant: 100 pivots of 1e6 and 100 of 1e-6 give det 1 in either order (%.17g, %.17g)", det1, det2);
    double ci[16] = {0};
    for (int i = 0; i < 4; i++) ci[i * 4 + i] = 1;
    ci[5] = INFINITY;
    CHECK(mk_chol(ci, 4) == -1, "an infinite diagonal is refused by Cholesky");
    /* documented contracts: non-finite A refused, subnormal pivots accepted with the non-finite answer reported */
    double na[16] = {0}, nb2[4] = {1, 1, 1, 1}, nx[4], ia[16] = {0}, ib[4] = {1, 1, 1, 1};
    for (int i = 0; i < 4; i++) na[i * 4 + i] = ia[i * 4 + i] = 1;
    na[0] = NAN, ia[6] = INFINITY;
    double rn1 = mk_lstsq(na, 4, 4, nb2, nx), rn2 = mk_lstsq(ia, 4, 4, ib, nx);
    CHECK(rn1 == -1 && rn2 == -1, "mk_lstsq refuses NaN or Inf in A with -1 (got %g, %g)", rn1, rn2);
    double ta[4] = {1e-310, 0, 0, 1e-310}, tb[2] = {1, 1};
    int tp[2];
    int tl = mk_lu(ta, 2, tp), ts = tl ? 0 : mk_lu_solve(ta, 2, tp, tb);
    double ga[4] = {2, 1, 1, 3}, gb[2] = {3, 4};
    int gp[2];
    mk_lu(ga, 2, gp);
    int gs = mk_lu_solve(ga, 2, gp, gb);
    double ca[4] = {1e-310, 0, 0, 1e-310}, cb[2] = {1, 1};
    int cl = mk_chol(ca, 2), cs = cl ? 0 : mk_chol_solve(ca, 2, cb);
    CHECK(tl == 0 && ts == -1 && gs == 0 && cl == 0 && cs == -1, "1e-310 I passes LU's relative test, and the solves report the non-finite answer (lu_solve %d, chol_solve %d; a normal solve %d)", ts, cs, gs);
    { /* subnormal symmetric matrices converge (power-of-two prescale), and bad sizes are refused */
        double s = 2.1e-313, A[9] = {0, s, 0, s, 0, s, 0, s, 0}, V[9], w[3];
        int sw = mk_sym_jacobi(A, 3, V, w);
        double e = sqrt(2.0) * s;
        CHECK(sw >= 0 && fabs(w[0] + e) <= 1e-8 * e && fabs(w[1]) <= 1e-8 * e && fabs(w[2] - e) <= 1e-8 * e,
              "sym_jacobi on a subnormal 3x3 (entries 2.1e-313): -sqrt2 s, 0, sqrt2 s (sweeps %d)", sw);
        double B[17 * 17], W[17 * 17], v17[17];
        for (int i = 0; i < 17; i++)
            for (int j = 0; j <= i; j++) B[i * 17 + j] = B[j * 17 + i] = (t_rand() + (i == j)) * 3.0e-313;
        sw = mk_sym_jacobi(B, 17, W, v17);
        CHECK(sw >= 0, "sym_jacobi on a subnormal 17x17 converges (sweeps %d)", sw);
        double one[1] = {1}, vv[1], ww[1];
        int pv[1];
        CHECK(mk_lu(one, -1, pv) == -1, "mk_lu refuses n = -1");
        CHECK(mk_sym_jacobi(one, -1, vv, ww) == -1, "mk_sym_jacobi refuses n = -1");
        CHECK(mk_lu(one, 46341, pv) == -1, "mk_lu refuses n = 46341 (n * n would overflow an int) before reading a");
        CHECK(mk_sym_jacobi(one, 46341, vv, ww) == -1, "mk_sym_jacobi refuses n = 46341 before reading a");
    }
    { /* size guards: negative n, and n * n or m * n past INT_MAX */
        double one[1] = {4}, b1[1] = {1}, x1[1];
        CHECK(mk_chol(one, -1) == -1, "mk_chol refuses n = -1");
        CHECK(mk_chol_solve(one, -1, b1) == -1, "mk_chol_solve refuses n = -1");
        CHECK(mk_chol(one, 46341) == -1, "mk_chol refuses n = 46341 before reading a");
        CHECK(mk_lstsq(one, INT_MAX / 2 + 1, 2, b1, x1) == -1, "mk_lstsq refuses m * n past INT_MAX before reading a");
        (void)x1;
    }
    { /* mk_tridiag: an analytic solution, backward error on random systems, a plain elimination's bits, the refusals */
        enum { N9 = 9, NR = 1000 };
        double ta[N9], tb[N9], tc[N9], td[N9], te = 0;
        for (int i = 0; i < N9; i++) ta[i] = tc[i] = -1, tb[i] = 2, td[i] = 1;
        int tr = mk_tridiag(N9, ta, tb, tc, td, td);
        for (int i = 0; i < N9; i++) te = fmax(te, fabs(td[i] - (i + 1) * (N9 - i) / 2.0) / ((i + 1) * (N9 - i) / 2.0));
        CHECK(tr == 0 && te < 1e-14, "tridiag, -x[i-1] + 2 x[i] - x[i+1] = 1 on 9 nodes: x = i (10 - i) / 2 to %.2g relative", te);
        static double ra[NR], rb[NR], rc[NR], rd[NR], rb0[NR], rd0[NR], rx[NR];
        double bw = 0;
        int rr = 0;
        for (int rep = 0; rep < 20; rep++) {
            for (int i = 0; i < NR; i++) {
                ra[i] = i ? t_rand() : 0, rc[i] = i + 1 < NR ? t_rand() : 0;
                rb0[i] = rb[i] = (fabs(ra[i]) + fabs(rc[i]) + 0.01 + fabs(t_rand())) * (t_rand() < 0 ? -1 : 1);
                rd0[i] = rd[i] = t_rand() * 1e3;
            }
            rr |= mk_tridiag(NR, ra, rb, rc, rd, rx);
            for (int i = 0; i < NR; i++) { /* componentwise backward error */
                double s = rb0[i] * rx[i] - rd0[i], m = fabs(rb0[i] * rx[i]) + fabs(rd0[i]);
                if (i) s += ra[i] * rx[i - 1], m += fabs(ra[i] * rx[i - 1]);
                if (i + 1 < NR) s += rc[i] * rx[i + 1], m += fabs(rc[i] * rx[i + 1]);
                bw = fmax(bw, fabs(s) / m);
            }
        }
        CHECK(rr == 0 && bw < 4 * DBL_EPSILON, "tridiag, 20 random diagonally dominant systems of 1000: residual %.2g of the row's terms (under 4 eps)", bw);
        /* the Thomas algorithm without pivoting, written out: mk_tridiag must give its bits */
        long tdiff = 0, ttot = 0;
        for (int rep = 0; rep < 200; rep++) {
            int n = 1 + rep % 40;
            double a[41], b[41], cc[41], d[41], x[41], b2[41], d2[41], x2[41];
            for (int i = 0; i < n; i++) {
                a[i] = i ? -fabs(t_rand()) : 0, cc[i] = i + 1 < n ? -fabs(t_rand()) : 0;
                b2[i] = b[i] = fabs(a[i]) + fabs(cc[i]) + fabs(t_rand()) * 1e3 + 1e-3, d2[i] = d[i] = t_rand() * 300;
            }
            for (int i = 1; i < n; i++) {
                const double f = a[i] / b[i - 1];
                b[i] -= f * cc[i - 1], d[i] -= f * d[i - 1];
            }
            x[n - 1] = d[n - 1] / b[n - 1];
            for (int i = n - 2; i >= 0; i--) x[i] = (d[i] - cc[i] * x[i + 1]) / b[i];
            if (mk_tridiag(n, a, b2, cc, d2, x2)) tdiff += n;
            else for (int i = 0; i < n; i++) tdiff += memcmp(&x[i], &x2[i], sizeof x[i]) != 0;
            ttot += n;
        }
        CHECK(tdiff == 0, "tridiag equals a plain Thomas elimination bit for bit: %ld of %ld values differ", tdiff, ttot);
        double za[2] = {0, 1}, zb[2] = {1, 1}, zc[2] = {1, 0}, zd[2] = {1, 2};
        double pa[2] = {0, 1}, pb[2] = {1e-17, 1}, pc[2] = {1, 0}, pd[2] = {1, 2};
        double na[2] = {0, 1}, nb[2] = {4, 4}, nc[2] = {1, 0}, nd[2] = {1, NAN};
        double ob[1] = {4}, od[1] = {2};
        int z1 = mk_tridiag(2, za, zb, zc, zd, zd), z2 = mk_tridiag(2, pa, pb, pc, pd, pd), z3 = mk_tridiag(2, na, nb, nc, nd, nd);
        int z4 = mk_tridiag(0, za, zb, zc, zd, zd), z5 = mk_tridiag(1, za, ob, zc, od, od);
        CHECK(z1 == -1 && z2 == -1 && z3 == -1 && z4 == -1 && z5 == 0 && od[0] == 0.5,
              "tridiag refuses a zero pivot (%d), a pivot 1e-17 of its row that needs pivoting (%d), a NaN right side (%d), n 0 (%d); n 1 gives 2 / 4 (%d)",
              z1, z2, z3, z4, z5);
    }
    return wrong;
}
