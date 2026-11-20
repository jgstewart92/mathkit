#include <limits.h>
#include <stdlib.h>
static size_t t_asked; /* what the header's last allocation asked for; t_refuse makes the next one fail */
static int t_refuse;
static void *t_malloc(size_t s)
{
    t_asked = s;
    if (t_refuse) { t_refuse = 0; return NULL; }
    return malloc(s);
}
/* t_rev makes the header's next qsort leave equal keys in reverse input order, as C allows (glibc 2.44's does too once its buffer
   cannot be had) */
static int t_rev;
static void t_qsort(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *))
{
    if (t_rev) {
        unsigned char *c = b, tmp;
        for (size_t i = 0; i < n / 2; i++)
            for (size_t k = 0; k < s; k++) tmp = c[i * s + k], c[i * s + k] = c[(n - 1 - i) * s + k], c[(n - 1 - i) * s + k] = tmp;
        t_rev = 0;
    }
    qsort(b, n, s, cmp);
}
#define malloc t_malloc
#define qsort t_qsort
#include "../include/mk_sparse.h"
#undef qsort
#undef malloc
#include "../include/mk_linalg.h"
#include "t.h"

/* the 5-point Laplacian on an m x m grid, Dirichlet (4 on the diagonal, -1 to each neighbour) */
static int lap2(mk_csr *A, int m)
{
    int n = m * m, k = 0, *ri = malloc(sizeof(int) * 10 * n), *ci = ri + 5 * n;
    double *v = malloc(sizeof(double) * 5 * n);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < m; i++) {
            int r = j * m + i;
            ri[k] = r, ci[k] = r, v[k++] = 4;
            if (i > 0) ri[k] = r, ci[k] = r - 1, v[k++] = -1;
            if (i < m - 1) ri[k] = r, ci[k] = r + 1, v[k++] = -1;
            if (j > 0) ri[k] = r, ci[k] = r - m, v[k++] = -1;
            if (j < m - 1) ri[k] = r, ci[k] = r + m, v[k++] = -1;
        }
    int e = mk_csr_from_triplets(A, n, k, ri, ci, v);
    free(ri), free(v);
    return e;
}
static double maxerr(const double *a, const double *b, int n)
{
    double e = 0;
    for (int i = 0; i < n; i++) e = fmax(e, fabs(a[i] - b[i]));
    return e;
}

/* diag(a (1 + 0.1 i)), n x n, as the verification work's demo.c */
static void diagm(mk_csr *A, int n, double a)
{
    int ri[64], ci[64];
    double v[64];
    for (int i = 0; i < n; i++) ri[i] = ci[i] = i, v[i] = a * (1 + 0.1 * i);
    mk_csr_from_triplets(A, n, n, ri, ci, v);
}
/* the true |b - A x| / |b| in long double, both scaled by b's largest entry; -1 if x is not finite */
static double truerel(const mk_csr *A, const double *b, const double *x)
{
    long double big = 0, num = 0, den = 0;
    for (int i = 0; i < A->n; i++) {
        if (!isfinite(x[i])) return -1;
        big = fmaxl(big, fabsl(b[i]));
    }
    for (int i = 0; i < A->n; i++) {
        long double s = b[i];
        for (int e = A->rowptr[i]; e < A->rowptr[i + 1]; e++) s -= (long double)A->val[e] * x[A->col[e]];
        num += (s / big) * (s / big), den += (b[i] / big) * (b[i] / big);
    }
    return (double)sqrtl(num / den);
}
/* solve with b all bval from x = 0: the return and the true relative residual */
static int solve1(const mk_csr *A, double bval, int pc, int bi, double tol, double *tr)
{
    double b[64], x[64] = {0};
    int it;
    for (int i = 0; i < A->n; i++) b[i] = bval;
    int r = bi ? mk_bicgstab(A, b, x, pc, tol, 100, &it, NULL) : mk_cg(A, b, x, pc, tol, 100, &it, NULL);
    *tr = truerel(A, b, x);
    return r;
}
static void lapscaled(mk_csr *A, int m, double s)
{
    lap2(A, m);
    for (int e = 0; e < A->rowptr[A->n]; e++) A->val[e] *= s;
}

int main(void)
{
    enum { M = 30, N = M * M };
    mk_csr A;
    lap2(&A, M);
    static double xs[N], b[N], x[N];
    for (int i = 0; i < N; i++) xs[i] = t_rand();
    mk_csr_matvec(&A, xs, b);
    double kap = 1 / pow(tan(T_PI / (2 * (M + 1))), 2), tol = 1e-10;
    int bound = (int)ceil(0.5 * sqrt(kap) * log(2 * sqrt(kap) / tol)), it, it_j, it_ic;
    double rr;
    for (int i = 0; i < N; i++) x[i] = 0;
    int r = mk_cg(&A, b, x, MK_PC_NONE, tol, 10 * N, &it, &rr);
    CHECK(!r && maxerr(x, xs, N) < 1e-8, "CG, 30 x 30 Laplacian: converged, error %.2g, residual %.2g", maxerr(x, xs, N), rr);
    CHECK(it <= bound && it >= 0.5 * sqrt(kap), "CG iterations %d within theory: sqrt(kappa) %.1f, bound 1/2 sqrt(kappa) ln(2 sqrt(kappa)/tol) = %d", it, sqrt(kap), bound);
    for (int i = 0; i < N; i++) x[i] = 0;
    mk_cg(&A, b, x, MK_PC_JACOBI, tol, 10 * N, &it_j, NULL);
    for (int i = 0; i < N; i++) x[i] = 0;
    int ric = mk_cg(&A, b, x, MK_PC_IC0, tol, 10 * N, &it_ic, &rr);
    CHECK(!ric && maxerr(x, xs, N) < 1e-8 && it_ic < 0.6 * it, "IC(0): error %.2g in %d iterations against CG's %d (Jacobi %d: the diagonal is constant)", maxerr(x, xs, N), it_ic, it, it_j);
    /* one discrete mode as right-hand side: CG is exact after one step */
    for (int j = 0; j < M; j++)
        for (int i = 0; i < M; i++) xs[j * M + i] = sin(T_PI * (i + 1) / (M + 1)) * sin(2 * T_PI * (j + 1) / (M + 1));
    mk_csr_matvec(&A, xs, b);
    for (int i = 0; i < N; i++) x[i] = 0;
    mk_cg(&A, b, x, MK_PC_NONE, 1e-12, 10, &it, NULL);
    CHECK(it == 1 && maxerr(x, xs, N) < 1e-13, "a single mode of the Laplacian is solved in %d step (error %.2g)", it, maxerr(x, xs, N));

    /* ck_spd3's layout: a chain of 40 nodes, 3x3 blocks, strong coupling inside each node */
    enum { NB = 40 };
    int row[NB + 1], col[2 * NB];
    double low[9 * 2 * NB];
    int e = 0;
    for (int i = 0; i < NB; i++) {
        row[i] = e;
        if (i > 0) {
            col[e] = i - 1;
            for (int a = 0; a < 9; a++) low[9 * e + a] = a % 4 == 0 ? -1 : 0;
            e++;
        }
        col[e] = i;
        for (int a = 0; a < 3; a++)
            for (int c = 0; c < 3; c++) low[9 * e + 3 * a + c] = a == c ? 7 + 0.1 * i : 2.9 - 0.01 * (a + c);
        e++;
    }
    row[NB] = e;
    mk_csr K;
    mk_csr_from_bsr3_lower(&K, NB, row, col, low);
    static double Kd[3 * NB * 3 * NB], kb[3 * NB], kx[3 * NB], kl[3 * NB];
    int piv[3 * NB], n3 = 3 * NB;
    for (int i = 0; i < n3 * n3; i++) Kd[i] = 0;
    for (int i = 0; i < NB; i++)
        for (int q = row[i]; q < row[i + 1]; q++)
            for (int a = 0; a < 3; a++)
                for (int c = 0; c < 3; c++) Kd[(3 * i + a) * n3 + 3 * col[q] + c] = Kd[(3 * col[q] + c) * n3 + 3 * i + a] = low[9 * q + 3 * a + c];
    double asym = 0;
    for (int i = 0; i < n3; i++)
        for (int q = K.rowptr[i]; q < K.rowptr[i + 1]; q++) asym = fmax(asym, fabs(K.val[q] - Kd[i * n3 + K.col[q]]));
    CHECK(asym == 0 && K.rowptr[n3] == NB * 9 + (NB - 1) * 18, "block lower triangle expands to the full matrix exactly (%d entries)", K.rowptr[n3]);
    for (int i = 0; i < n3; i++) kb[i] = t_rand(), kx[i] = 0, kl[i] = kb[i];
    int ib, ij;
    int rb = mk_cg(&K, kb, kx, MK_PC_BJACOBI3, 1e-12, 1000, &ib, NULL);
    for (int i = 0; i < n3; i++) x[i] = 0;
    mk_cg(&K, kb, x, MK_PC_JACOBI, 1e-12, 1000, &ij, NULL);
    mk_lu(Kd, n3, piv), mk_lu_solve(Kd, n3, piv, kl);
    CHECK(!rb && maxerr(kx, kl, n3) < 1e-10 && ib < ij, "3x3 block-Jacobi PCG matches a dense LU (%.2g) in %d iterations, scalar Jacobi %d", maxerr(kx, kl, n3), ib, ij);

    /* periodic 1-D Laplacian: singular, constants in its null space; a mean-free right-hand side */
    enum { P = 200 };
    int pr[3 * P], pc[3 * P];
    double pv[3 * P], pb[P], px[P], pxs[P];
    for (int i = 0; i < P; i++) {
        pr[3 * i] = pr[3 * i + 1] = pr[3 * i + 2] = i;
        pc[3 * i] = i, pv[3 * i] = 2, pc[3 * i + 1] = (i + 1) % P, pv[3 * i + 1] = -1, pc[3 * i + 2] = (i + P - 1) % P, pv[3 * i + 2] = -1;
        pxs[i] = sin(2 * T_PI * 3 * i / P) + 0.3 * cos(2 * T_PI * 7 * i / P), px[i] = 0;
    }
    mk_csr Pm;
    mk_csr_from_triplets(&Pm, P, 3 * P, pr, pc, pv);
    mk_csr_matvec(&Pm, pxs, pb);
    int ip, rp = mk_cg(&Pm, pb, px, MK_PC_JACOBI | MK_PC_MEANFREE, 1e-12, 2000, &ip, NULL);
    CHECK(!rp && maxerr(px, pxs, P) < 1e-9, "periodic (singular) Laplacian with the mean kept out: error %.2g, %d iterations", maxerr(px, pxs, P), ip);

    /* nonsymmetric: convection-diffusion, upwind, 30 x 30 */
    int k = 0, *ri = malloc(sizeof(int) * 10 * N), *ci = ri + 5 * N;
    double *v = malloc(sizeof(double) * 5 * N), cx = 3, cy = -2;
    for (int j = 0; j < M; j++)
        for (int i = 0; i < M; i++) {
            int r0 = j * M + i;
            ri[k] = r0, ci[k] = r0, v[k++] = 4 + fabs(cx) + fabs(cy);
            if (i > 0) ri[k] = r0, ci[k] = r0 - 1, v[k++] = -1 - (cx > 0 ? cx : 0);
            if (i < M - 1) ri[k] = r0, ci[k] = r0 + 1, v[k++] = -1 + (cx < 0 ? cx : 0);
            if (j > 0) ri[k] = r0, ci[k] = r0 - M, v[k++] = -1 - (cy > 0 ? cy : 0);
            if (j < M - 1) ri[k] = r0, ci[k] = r0 + M, v[k++] = -1 + (cy < 0 ? cy : 0);
        }
    mk_csr C;
    mk_csr_from_triplets(&C, N, k, ri, ci, v);
    for (int i = 0; i < N; i++) xs[i] = t_rand(), x[i] = 0;
    mk_csr_matvec(&C, xs, b);
    int ic, rc = mk_bicgstab(&C, b, x, MK_PC_JACOBI, 1e-11, 2000, &ic, &rr);
    CHECK(!rc && maxerr(x, xs, N) < 1e-8, "BiCGSTAB, upwind convection-diffusion 30 x 30: error %.2g in %d iterations", maxerr(x, xs, N), ic);
    for (int i = 0; i < N; i++) x[i] = 0;
    CHECK(mk_cg(&C, b, x, MK_PC_NONE, 1e-11, 300, &ic, NULL) != 0 || maxerr(x, xs, N) > 1e-6, "CG on the same nonsymmetric matrix does not solve it (as it should not)");
    free(ri), free(v), mk_csr_free(&A), mk_csr_free(&K), mk_csr_free(&Pm), mk_csr_free(&C);
    { /* 1.0.6: the verification work's sparse findings (Result 80) */
        mk_csr D;
        double t1, t2, t3, t4;
        diagm(&D, 3, 1);
#ifndef T_OLD
        int f1 = mk_cg(&D, (double[3]){1, 1, 1}, (double[3]){0}, 4, 1e-10, 10, NULL, NULL);
        int f1b = mk_bicgstab(&D, (double[3]){1, 1, 1}, (double[3]){0}, 7, 1e-10, 10, NULL, NULL);
        CHECK(f1 == -1 && f1b == -1, "F1: preconditioner kinds outside 0..3 are refused (CG kind 4: %d, BiCGSTAB kind 7: %d)", f1, f1b);
#endif
        int a1 = solve1(&D, 1e-170, MK_PC_NONE, 0, 1e-10, &t1), a2 = solve1(&D, 1e-170, MK_PC_NONE, 1, 1e-10, &t2);
        CHECK(a1 == 0 && a2 == 0 && t1 >= 0 && t1 < 1e-9 && t2 >= 0 && t2 < 1e-9,
              "F2a: b = 1e-170 (|b|^2 underflows) is solved, not passed back as x0: CG %d, true residual %.2g; BiCGSTAB %d, %.2g", a1, t1, a2, t2);
        int b1 = solve1(&D, 1e170, MK_PC_NONE, 0, 1e-10, &t1), b2 = solve1(&D, 1e170, MK_PC_NONE, 1, 1e-10, &t2);
        CHECK(b1 == 0 && b2 == 0 && t1 >= 0 && t1 < 1e-9 && t2 >= 0 && t2 < 1e-9,
              "F2b: b = 1e170 (|b|^2 overflows) is solved: CG %d, true residual %.2g; BiCGSTAB %d, %.2g", b1, t1, b2, t2);
        mk_csr_free(&D);
        diagm(&D, 3, 1e-87);
        int c1 = solve1(&D, 1e-121, MK_PC_NONE, 0, 1e-10, &t1), c2 = solve1(&D, 1e-121, MK_PC_JACOBI, 0, 1e-10, &t2);
        CHECK(c1 == 0 && c2 == 0 && t1 >= 0 && t1 < 1e-9, "F2c: A = 1e-87 I, b = 1e-121 (p A p underflowed): CG %d, true residual %.2g, as Jacobi (%d)", c1, t1, c2);
        mk_csr_free(&D);
        { /* the scaling is by powers of two: b times 2^600 gives x times 2^600, bit for bit, in the same iterations */
            mk_csr L;
            lap2(&L, 12);
            int nl = L.n, same = 1, i1, i2, i3, i4;
            static double bl[144], bh[144], x1[144], x2[144], x3[144], x4[144];
            for (int i = 0; i < nl; i++) bl[i] = t_rand(), bh[i] = ldexp(bl[i], 600), x1[i] = x2[i] = x3[i] = x4[i] = 0;
            int q1 = mk_cg(&L, bl, x1, MK_PC_JACOBI, 1e-12, 500, &i1, NULL), q2 = mk_cg(&L, bh, x2, MK_PC_JACOBI, 1e-12, 500, &i2, NULL);
            int q3 = mk_bicgstab(&L, bl, x3, MK_PC_NONE, 1e-12, 500, &i3, NULL), q4 = mk_bicgstab(&L, bh, x4, MK_PC_NONE, 1e-12, 500, &i4, NULL);
            for (int i = 0; i < nl; i++) same &= ldexp(x1[i], 600) == x2[i] && ldexp(x3[i], 600) == x4[i];
            CHECK(q1 == 0 && q2 == 0 && q3 == 0 && q4 == 0 && same && i1 == i2 && i3 == i4,
                  "F2: b times 2^600 gives x times 2^600 bit for bit (CG %d and %d iterations, BiCGSTAB %d and %d)", i1, i2, i3, i4);
            mk_csr_free(&L);
        }
        diagm(&D, 3, 1e-101);
        int d1 = solve1(&D, 1, MK_PC_BJACOBI3, 0, 1e-10, &t1);
        mk_csr_free(&D);
        diagm(&D, 3, 1e110);
        int d2 = solve1(&D, 1, MK_PC_BJACOBI3, 0, 1e-10, &t2);
        mk_csr_free(&D);
        int sing = 0; /* a singular block, exact at every power-of-two scale, still gets a zero inverse */
        for (int k = 0; k < 3; k++) {
            int ri[9], ci[9];
            double v[9], u[3] = {1, 2, 3}, s = k == 0 ? 1 : k == 1 ? 0x1p-700 : 0x1p700;
            for (int a = 0; a < 3; a++) for (int c = 0; c < 3; c++) ri[3 * a + c] = a, ci[3 * a + c] = c, v[3 * a + c] = u[a] * u[c] * s;
            mk_csr S;
            mk_csr_from_triplets(&S, 3, 9, ri, ci, v);
            mk_pc_ P = {0};
            if (!mk_pc_build_(&P, &S, MK_PC_BJACOBI3)) for (int t = 0; t < 9; t++) sing += P.d[t] != 0;
            else sing += 100;
            mk_pc_free_(&P), mk_csr_free(&S);
        }
        CHECK(d1 == 0 && d2 == 0 && t1 >= 0 && t1 < 1e-9 && t2 >= 0 && t2 < 1e-9 && sing == 0,
              "F3: block-Jacobi on 1e-101 I and 1e110 I solves (%d %d, true residuals %.2g %.2g); a singular block at scales 1, 2^-700, 2^700 still gets a zero inverse (%d entries not 0)",
              d1, d2, t1, t2, sing);
        mk_csr E = {7, (int *)&sing, NULL, NULL};
        int e5 = mk_csr_from_triplets(&E, -1, 0, NULL, NULL, NULL);
        CHECK(e5 == -1 && E.n == 0 && !E.rowptr, "F5: a matrix of n -1 is refused (%d, A.n %d)", e5, E.n);
        if (!e5) mk_csr_free(&E);
#ifndef T_OLD
        int z0 = 0;
        double z1 = 1;
        int e6 = mk_csr_from_triplets(&E, 1, -(1L << 60), &z0, &z0, &z1);
        CHECK(e6 == -1, "F6: nnz = -2^60 is refused before allocating (1.0.5: malloc(0), then qsort of 2^64 - 2^60 entries)");
#endif
        mk_csr H = {1 << 29, NULL, NULL, NULL};
        double hb = 0, hx = 0;
        t_refuse = 1, t_asked = 0;
        int g1 = mk_cg(&H, &hb, &hx, MK_PC_NONE, 1e-8, 1, NULL, NULL);
        size_t s1 = t_asked;
        H.n = 1 << 28, t_refuse = 1, t_asked = 0;
        int g2 = mk_bicgstab(&H, &hb, &hx, MK_PC_NONE, 1e-8, 1, NULL, NULL);
        size_t s2 = t_asked;
        H.n = 715827885, t_refuse = 1, t_asked = 0; /* 3 n past INT_MAX */
        mk_pc_ HP = {0};
        int g3 = mk_pc_build_(&HP, &H, MK_PC_BJACOBI3);
        size_t s3 = t_asked;
        t_refuse = 0;
        CHECK(g1 == -1 && g2 == -1 && g3 == -1 && s1 == 32 * ((size_t)1 << 29) && s2 == 64 * ((size_t)1 << 28) && s3 == 24 * (size_t)715827885,
              "F7: work sizes formed in size_t: CG at n 2^29 asks %zu bytes, BiCGSTAB at 2^28 %zu, block-Jacobi at 715827885 %zu (16, 16 and 16 GiB)", s1, s2, s3);
        int brow[2] = {0, 1}, bcol[1] = {0};
        double blow[9] = {4, 0, 0, 0, 4, 0, 0, 0, 4};
        mk_csr B = {77, (int *)&sing, (int *)&sing, NULL};
        t_refuse = 1;
        int e8 = mk_csr_from_bsr3_lower(&B, 1, brow, bcol, blow);
        t_refuse = 0;
        CHECK(e8 == -1 && B.n == 0 && !B.rowptr && !B.col, "F8: bsr3 with no memory resets *A (%d, A.n %d)", e8, B.n);
#ifndef T_OLD
        int e8b = mk_csr_from_bsr3_lower(&B, -1, brow + 1, bcol, blow), e8c = mk_csr_from_bsr3_lower(&B, INT_MAX / 3 + 1, brow, bcol, blow);
        CHECK(e8b == -1 && e8c == -1, "F8: bsr3 refuses nb -1 and nb past INT_MAX / 3 before reading row[nb]");
#endif
        { /* F9: CG at tol 0 runs past attainable accuracy; it stops with a finite x at any scale of b */
            mk_csr L;
            lap2(&L, 12);
            int nl = L.n, its[3], rets[3];
            double trs[3];
            static double b9[144], x9[144];
            for (int k = 0; k < 3; k++) {
                double sb = k == 0 ? 1 : k == 1 ? 0x1p837 : 0x1p-900;
                for (int i = 0; i < nl; i++) b9[i] = sb * (1 + 0.3 * sin(i)), x9[i] = 0;
                rets[k] = mk_cg(&L, b9, x9, MK_PC_JACOBI, 0, 100000, &its[k], NULL);
                trs[k] = truerel(&L, b9, x9);
            }
            CHECK(rets[0] == 1 && rets[1] == 1 && rets[2] == 1 && trs[0] >= 0 && trs[0] < 1e-13 && trs[1] == trs[0] && trs[2] == trs[0] && its[1] == its[0] && its[2] == its[0],
                  "F9: CG with tol 0 ends (1) at attainable accuracy, true residual %.2g, finite x, at b 1, 2^837 and 2^-900 alike (%d %d %d iterations; returns %d %d %d)",
                  trs[0], its[0], its[1], its[2], rets[0], rets[1], rets[2]);
            mk_csr G;
            lapscaled(&G, 12, 1e290);
            int ng = G.n, i9, r9;
            for (int i = 0; i < ng; i++) b9[i] = t_rand(), x9[i] = 0;
            r9 = mk_cg(&G, b9, x9, MK_PC_JACOBI, 1e-10, 1000, &i9, NULL);
            t3 = truerel(&G, b9, x9);
            CHECK(r9 == 0 && t3 >= 0 && t3 < 1e-9, "F9: the underflow stop does not end a solve early: A near 1e290 with Jacobi (r z near 1e-290) converges, %d in %d iterations, true residual %.2g", r9, i9, t3);
            mk_csr_free(&L), mk_csr_free(&G);
        }
        diagm(&D, 3, 1e-300);
        int h1 = solve1(&D, 1e300, MK_PC_NONE, 0, 1e-10, &t1), h2 = solve1(&D, 1e300, MK_PC_NONE, 1, 1e-10, &t2), h3 = solve1(&D, 1e300, MK_PC_NONE, 1, NAN, &t4);
        mk_csr_free(&D);
        diagm(&D, 3, 1);
        int h4 = solve1(&D, 1, MK_PC_NONE, 1, NAN, &t3);
        mk_csr_free(&D);
        CHECK(h1 == -1 && h2 == -1 && h3 == -1 && (h4 < 0 || t3 >= 0),
              "F10: x past the largest double (1e300 / 1e-300) gives -1 from CG (%d) and BiCGSTAB (%d, and %d with tol NaN), never a status with x infinite; tol NaN on I gives %d with x finite",
              h1, h2, h3, h4);
        (void)t4;
    }
    { /* 1.2.1: the verification work F11 (Result 83): A near 2^-1015 is not rescaled, so a * p overflowed before the multiply by the b scale */
        enum { N = 40 };
        int ri[3 * N], ci[3 * N], k = 0;
        double v[3 * N], s = ldexp(1, -1015);
        for (int i = 0; i < N; i++) {
            ri[k] = i, ci[k] = i, v[k++] = 2.5 * s;
            if (i > 0) ri[k] = i, ci[k] = i - 1, v[k++] = -s;
            if (i + 1 < N) ri[k] = i, ci[k] = i + 1, v[k++] = -s;
        }
        mk_csr A;
        int ok = mk_csr_from_triplets(&A, N, k, ri, ci, v) == 0;
        double b[N], x[N], rr = -1;
        for (int i = 0; i < N; i++) b[i] = s, x[i] = -1000;
        int it = -1, ret = ok ? mk_cg(&A, b, x, MK_PC_NONE, 1e-13, 1000, &it, &rr) : -9;
        long double res = 0;
        for (int i = 0; i < N; i++) {
            long double q = 2.5L * x[i] - (i > 0 ? x[i - 1] : 0) - (i + 1 < N ? x[i + 1] : 0);
            res += (1 - q) * (1 - q);
        }
        CHECK((ret == 0 || (ret == 1 && rr > 1e-13)) && sqrtl(res / N) < 1e-12, "F11: CG on 2^-1015 tridiag(-1, 2.5, -1), x0 -1000: solves (ret %d, relres %.3g; 1.6.1: 1 where the true residual is past tol 1e-13; true relres %.3Lg)", ret, rr,
              sqrtl(res / N));
        double b2[N], x2[N];
        for (int i = 0; i < N; i++) b2[i] = 1e308, x2[i] = 0; /* A is s-scaled: the answer overflows; refused, never a wrong finite x */
        ret = ok ? mk_cg(&A, b2, x2, MK_PC_NONE, 1e-13, 1000, &it, &rr) : -9;
        int fin = 1;
        for (int i = 0; i < N; i++) fin &= isfinite(x2[i]) != 0;
        CHECK(ret == -1 || (ret >= 0 && fin), "b = 1e308 against A near 2^-1015: an unrepresentable answer is refused, not returned (ret %d)", ret);
        if (ok) mk_csr_free(&A);
    }
    { /* 1.6.0, the review's B1: a full block matrix (a block above the diagonal) is refused; 1.5.6 mirrored it and summed it twice */
        double fl[36];
        for (int k = 0; k < 36; k++) fl[k] = k % 9 % 4 == 0 ? (k / 9 % 3 == 0 ? 10 : 2) : 0; /* diagonal blocks 10 I, off-diagonal 2 I */
        const int frow[3] = {0, 2, 4}, fcol[4] = {0, 1, 0, 1};
        mk_csr F;
        const int rf = mk_csr_from_bsr3_lower(&F, 2, frow, fcol, fl);
        CHECK(rf == -1 && F.n == 0 && !F.rowptr, "B1 (review): bsr3 lower refuses a block above the diagonal (%d)", rf);
        if (!rf) mk_csr_free(&F);
    }
    { /* review 2 #10 (1.6.1): 0 and relres were the recurrence's residual; from x0 = 1e6 random the true one was 173 times larger
         (1.6.0: CG 0 with relres 8.9e-13, true 1.5e-10). Now the true residual is formed again and the solve restarts from x. */
        enum { N1 = 500 };
        static int ri[3 * N1], ci[3 * N1];
        static double v[3 * N1], xs[N1], b1[N1], x1[N1], x2[N1];
        long k = 0;
        for (int i = 0; i < N1; i++) {
            ri[k] = i, ci[k] = i, v[k++] = 2;
            if (i) ri[k] = i, ci[k] = i - 1, v[k++] = -1;
            if (i + 1 < N1) ri[k] = i, ci[k] = i + 1, v[k++] = -1;
            xs[i] = t_rand();
        }
        mk_csr L;
        mk_csr_from_triplets(&L, N1, k, ri, ci, v);
        mk_csr_matvec(&L, xs, b1);
        for (int i = 0; i < N1; i++) x1[i] = x2[i] = 1e6 * t_rand();
        int i1, i2;
        double q1, q2;
        const int r1 = mk_cg(&L, b1, x1, MK_PC_NONE, 1e-12, 100000, &i1, &q1), r2 = mk_bicgstab(&L, b1, x2, MK_PC_NONE, 1e-12, 100000, &i2, &q2);
        const double t1 = truerel(&L, b1, x1), t2 = truerel(&L, b1, x2);
        CHECK(r1 == 0 && t1 <= 1e-12 && fabs(q1 - t1) <= 1e-2 * t1 && r2 == 0 && t2 <= 1e-12 && fabs(q2 - t2) <= 1e-2 * t2,
              "review 2 #10: 1-D Laplacian of 500 from x0 1e6 random, tol 1e-12: CG %d in %d iterations, relres %.3g, true %.3g; BiCGSTAB %d in %d, %.3g, true %.3g", r1, i1, q1, t1, r2, i2, q2, t2);
        mk_csr_free(&L);
    }
    { /* review 2 #18 (1.6.1): a negative r z (Jacobi on a negative diagonal entry) stopped CG with 1, not the breakdown's -1; BiCGSTAB dropped
         MK_PC_MEANFREE and iterated without it */
        enum { N2 = 50 };
        int ri[3 * N2], ci[3 * N2];
        double v[3 * N2], b2[N2], x[N2];
        long k = 0;
        for (int i = 0; i < N2; i++) {
            ri[k] = i, ci[k] = i, v[k++] = i == 7 ? -2 : 2;
            if (i) ri[k] = i, ci[k] = i - 1, v[k++] = -1;
            if (i + 1 < N2) ri[k] = i, ci[k] = i + 1, v[k++] = -1;
            b2[i] = 1, x[i] = 0;
        }
        mk_csr D;
        mk_csr_from_triplets(&D, N2, k, ri, ci, v);
        int ic = -1, ib = -1;
        const int rc = mk_cg(&D, b2, x, MK_PC_JACOBI, 1e-10, 1000, &ic, NULL);
        for (int i = 0; i < N2; i++) x[i] = 0;
        const int rb = mk_bicgstab(&D, b2, x, MK_PC_NONE | MK_PC_MEANFREE, 1e-10, 1000, &ib, NULL);
        CHECK(rc == -1 && rb == -1 && ib == 0, "review 2 #18: CG with Jacobi on a negative diagonal entry %d (after %d iterations); BiCGSTAB given MK_PC_MEANFREE %d after %d", rc, ic, rb, ib);
        mk_csr_free(&D);
    }
    { /* review 2 A2 (1.6.2): with tol under the attainable accuracy each pass reached tol on the recurrence and restarted, to maxit
         (1.6.1: 1e6 iterations on t_sparse F9's 144 unknowns); now a restart that does not lower b - A x ends with 1 */
        mk_csr L;
        lap2(&L, 12);
        static double b[144], x1[144], x2[144];
        for (int i = 0; i < 144; i++) b[i] = 1 + 0.3 * sin(i), x1[i] = x2[i] = 0;
        int i1 = -1, i2 = -1;
        double q1 = -1, q2 = -1;
        const int r1 = mk_cg(&L, b, x1, MK_PC_JACOBI, 1e-15, 1000000, &i1, &q1), r2 = mk_bicgstab(&L, b, x2, MK_PC_NONE, 1e-15, 1000000, &i2, &q2);
        const double t1 = truerel(&L, b, x1), t2 = truerel(&L, b, x2);
        CHECK(r1 == 1 && i1 < 1000 && fabs(q1 - t1) <= 0.5 * t1 && t1 < 1e-14 && r2 == 1 && i2 < 1000 && fabs(q2 - t2) <= 0.5 * t2 && t2 < 1e-14, /* relres in double, true in long double */
              "review 2 A2: tol 1e-15 under the floor: CG %d after %d iterations (relres %.3g, true %.3g), BiCGSTAB %d after %d (%.3g, %.3g)", r1, i1, q1, t1, r2, i2, q2, t2);
        mk_csr_free(&L);
    }
    { /* review 2 N4 (1.6.2): precond bits past MK_PC_MEANFREE were ignored (32 ran as MK_PC_NONE, -16 as NONE | MEANFREE) */
        mk_csr D;
        diagm(&D, 3, 2);
        const double b[3] = {1, 2, 3};
        double x[3] = {0, 0, 0};
        int it = 0;
        const int c32 = mk_cg(&D, b, x, 32, 1e-12, 100, &it, NULL), cm = mk_cg(&D, b, x, -16, 1e-12, 100, &it, NULL);
        const int b32 = mk_bicgstab(&D, b, x, MK_PC_JACOBI | 64, 1e-12, 100, &it, NULL), c0 = mk_cg(&D, b, x, MK_PC_NONE, 1e-12, 100, &it, NULL);
        CHECK(c32 == -1 && cm == -1 && b32 == -1 && c0 == 0, "review 2 N4: unknown precond bits refused: CG 32 %d, -16 %d, BiCGSTAB JACOBI | 64 %d (MK_PC_NONE %d)", c32, cm, b32, c0);
        mk_csr_free(&D);
    }
    { /* review 2 #15 (1.6.2): repeated positions are summed in input order whatever order qsort leaves equal keys in; 1.6.1 took qsort's
         (1, 1e16, -1e16 sums to 0 in input order, 1 reversed) */
        const int ri[5] = {0, 1, 1, 0, 1}, ci[5] = {0, 1, 1, 0, 1};
        const double v[5] = {2, 1, 1e16, 3, -1e16};
        mk_csr A, B;
        const int ea = mk_csr_from_triplets(&A, 2, 5, ri, ci, v);
        t_rev = 1;
        const int eb = mk_csr_from_triplets(&B, 2, 5, ri, ci, v);
        CHECK(ea == 0 && eb == 0 && A.val[1] == 0 && B.val[1] == 0 && A.val[0] == 5 && B.val[0] == 5,
              "review 2 #15: (1, 1e16, -1e16) at one position sums to %g, and %g when qsort reverses equal keys (0 in input order)", A.val[1], B.val[1]);
        mk_csr_free(&A), mk_csr_free(&B);
    }
    return wrong;
}
