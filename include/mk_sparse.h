/* mathkit: sparse matrices and Krylov solvers. CSR (compressed rows, columns sorted, duplicates summed) built from
   triplets, or from a symmetric matrix's lower triangle in 3x3 block CSR. Conjugate gradients with no, Jacobi, 3x3
   block-Jacobi or incomplete-Cholesky IC(0) preconditioning; BiCGSTAB for nonsymmetric systems. Both solvers: x holds
   the first guess and gets the answer; tol is on |b - A x| / |b|; return 0 converged, 1 not within maxit, -1 breakdown,
   no memory, a bad argument or an answer that is not finite. The solvers work in units of a power of two near b's
   largest entry, so any scale of b gives the same iterations and x scaled exactly. A itself is not rescaled: entries
   past about 1e300 or below 1e-300 can still under- or overflow the dots. Threads (mk_par.h, with -fopenmp): the
   product, the dots, the vector updates and the Jacobi and block-Jacobi steps run on blocks of rows; IC(0) is serial. A
   dot over more than one block (n > 4096) adds the blocks' partials in block order, so it differs in rounding from one
   serial sum but not between thread counts. Header only, C23. */
#ifndef MK_SPARSE_H
#define MK_SPARSE_H
#include "mk_par.h"
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int n;
    int *rowptr, *col; /* rowptr[n + 1]; col sorted within each row */
    double *val;
} mk_csr;

static inline void mk_csr_free(mk_csr *A)
{
    free(A->rowptr), free(A->col), free(A->val);
    A->rowptr = A->col = NULL, A->val = NULL, A->n = 0;
}

/* k, the triplet's place in the input, until the sort is done */
typedef struct { int r, c; union { long k; double v; }; } mk_trip_;
static inline int mk_trip_cmp_(const void *x, const void *y)
{
    const mk_trip_ *a = x, *b = y;
    return a->r != b->r ? (a->r > b->r) - (a->r < b->r) : a->c != b->c ? (a->c > b->c) - (a->c < b->c) : (a->k > b->k) - (a->k < b->k);
}
/* n x n from nnz triplets (row, col, value); repeated positions are summed in input order with any C library (the sort
   breaks ties by input place, so qsort's order of equal keys does not matter). 0, or -1 for a negative size, an index
   out of range or no memory. At most INT_MAX distinct positions (rowptr is int; about 2.1e9 entries, near 100 GB while
   building): not checked. */
static inline int mk_csr_from_triplets(mk_csr *A, int n, long nnz, const int *ri, const int *ci, const double *v)
{
    A->n = 0, A->rowptr = A->col = NULL, A->val = NULL;
    if (n < 0 || nnz < 0) return -1;
    mk_trip_ *t = malloc(sizeof *t * (size_t)(nnz ? nnz : 1));
    if (!t) return -1;
    for (long k = 0; k < nnz; k++) {
        if (ri[k] < 0 || ri[k] >= n || ci[k] < 0 || ci[k] >= n) { free(t); return -1; }
        t[k] = (mk_trip_){ri[k], ci[k], {k}};
    }
    qsort(t, (size_t)nnz, sizeof *t, mk_trip_cmp_);
    for (long k = 0; k < nnz; k++) t[k].v = v[t[k].k];
    long m = 0;
    for (long k = 0; k < nnz; k++)
        if (m && t[m - 1].r == t[k].r && t[m - 1].c == t[k].c) t[m - 1].v += t[k].v;
        else t[m++] = t[k];
    A->rowptr = calloc((size_t)n + 1, sizeof *A->rowptr);
    A->col = malloc(sizeof *A->col * (size_t)(m ? m : 1));
    A->val = malloc(sizeof *A->val * (size_t)(m ? m : 1));
    if (!A->rowptr || !A->col || !A->val) { free(t); mk_csr_free(A); return -1; }
    for (long k = 0; k < m; k++) A->rowptr[t[k].r + 1]++, A->col[k] = t[k].c, A->val[k] = t[k].v;
    for (int i = 0; i < n; i++) A->rowptr[i + 1] += A->rowptr[i];
    A->n = n;
    free(t);
    return 0;
}
/* the full scalar matrix of a symmetric matrix of 3x3 blocks given by its lower triangle in block CSR: row has nb + 1
   entries, block row i holds blocks row[i] .. row[i + 1] - 1, block e sits in block column col[e] <= i, and lower[9 e]
   starts its 3x3 values, row-major (a diagonal block in full). Each off-diagonal block is also placed, transposed,
   above the diagonal. 3 nb unknowns. 0, or -1 for a block above the diagonal, a bad nb or no memory. */
static inline int mk_csr_from_bsr3_lower(mk_csr *A, int nb, const int *row, const int *col, const double *lower)
{
    A->n = 0, A->rowptr = A->col = NULL, A->val = NULL;
    if (nb < 0 || nb > INT_MAX / 3 || row[nb] < 0) return -1;
    long nnzb = row[nb], cap = 18 * nnzb, k = 0;
    int *ri = malloc(sizeof(int) * (size_t)(2 * cap + 1));
    double *v = malloc(sizeof(double) * (size_t)(cap + 1));
    if (!ri || !v) { free(ri), free(v); return -1; }
    int *ci = ri + cap;
    for (int i = 0; i < nb; i++)
        for (int e = row[i]; e < row[i + 1]; e++) {
            int j = col[e];
            if (j > i) { free(ri), free(v); return -1; }
            const double *b = lower + 9L * e;
            for (int a = 0; a < 3; a++)
                for (int c = 0; c < 3; c++) {
                    ri[k] = 3 * i + a, ci[k] = 3 * j + c, v[k++] = b[3 * a + c];
                    if (j != i) ri[k] = 3 * j + c, ci[k] = 3 * i + a, v[k++] = b[3 * a + c];
                }
        }
    int r = mk_csr_from_triplets(A, 3 * nb, k, ri, ci, v);
    free(ri), free(v);
    return r;
}
/* y = A x on rows [i0, i1): the serial form each thread runs */
static inline void mk_csr_matvec_rows_(const mk_csr *A, const double *x, double *y, long i0, long i1)
{
    for (long i = i0; i < i1; i++) {
        double s = 0;
        for (int e = A->rowptr[i]; e < A->rowptr[i + 1]; e++) s += A->val[e] * x[A->col[e]];
        y[i] = s;
    }
}
/* y = A x; x and y must not overlap */
static inline void mk_csr_matvec(const mk_csr *A, const double *x, double *y)
{
#if defined(MK_MUTATE) && MK_MUTATE == 108
    /* planted fault: rows cut by the thread count, each piece's last row left out */
    MK_PAR_CUT_(A->n, i0, i1, mk_csr_matvec_rows_(A, x, y, i0, i1));
#else
    MK_PAR_FOR_(A->n, i0, i1, mk_csr_matvec_rows_(A, x, y, i0, i1));
#endif
}

enum { MK_PC_NONE = 0, MK_PC_JACOBI = 1, MK_PC_BJACOBI3 = 2, MK_PC_IC0 = 3, MK_PC_MEANFREE = 16 };

/* IC(0): L on the lower pattern of A, A ~ L L^T. A breakdown (a pivot <= 0) retries with the diagonal raised by
   1e-3, 4e-3, ... up to 0.256 (Manteuffel 1980's shift). */
typedef struct { int n; int *rowptr, *col; double *val; double shift; } mk_ic0_;
static inline void mk_ic0_free_(mk_ic0_ *L) { free(L->rowptr), free(L->col), free(L->val); }
static inline int mk_ic0_build_(mk_ic0_ *L, const mk_csr *A)
{
    int n = A->n;
    long m = 0;
    L->rowptr = malloc(sizeof(int) * ((size_t)n + 1));
    for (int i = 0; i < n; i++)
        for (int e = A->rowptr[i]; e < A->rowptr[i + 1]; e++) m += A->col[e] <= i;
    L->col = malloc(sizeof(int) * (size_t)(m ? m : 1));
    L->val = malloc(sizeof(double) * (size_t)(m ? m : 1));
    if (!L->rowptr || !L->col || !L->val) { mk_ic0_free_(L); return -1; }
    L->n = n;
    long k = 0;
    L->rowptr[0] = 0;
    for (int i = 0; i < n; i++) {
        int hasd = 0;
        for (int e = A->rowptr[i]; e < A->rowptr[i + 1]; e++)
            if (A->col[e] <= i) L->col[k] = A->col[e], hasd |= A->col[e] == i, k++;
        if (!hasd || L->col[k - 1] != i) { mk_ic0_free_(L); return -1; } /* the diagonal must be stored */
        L->rowptr[i + 1] = (int)k;
    }
    for (double shift = 0; shift < 1; shift = shift ? 4 * shift : 1e-3) {
        k = 0;
        for (int i = 0; i < n; i++)
            for (int e = A->rowptr[i]; e < A->rowptr[i + 1]; e++)
                if (A->col[e] <= i) L->val[k++] = A->val[e] * (A->col[e] == i ? 1 + shift : 1);
        int ok = 1;
        for (int i = 0; i < n && ok; i++) {
            int r0 = L->rowptr[i], r1 = L->rowptr[i + 1];
            for (int e = r0; e < r1; e++) {
                int j = L->col[e];
                double s = L->val[e];
                int p = r0, q = L->rowptr[j], qe = j == i ? e : L->rowptr[j + 1] - 1; /* row j without its diagonal */
                while (p < e && q < qe) {
                    if (L->col[p] == L->col[q]) s -= L->val[p++] * L->val[q++];
                    else if (L->col[p] < L->col[q]) p++;
                    else q++;
                }
                if (j < i) L->val[e] = s / L->val[L->rowptr[j + 1] - 1];
                else {
                    if (!(s > 0)) { ok = 0; break; }
                    L->val[e] = sqrt(s);
                }
            }
        }
        if (ok) { L->shift = shift; return 0; }
    }
    mk_ic0_free_(L);
    return -1;
}
static inline void mk_ic0_apply_(const mk_ic0_ *L, const double *r, double *z)
{
    int n = L->n;
    for (int i = 0; i < n; i++) {
        double s = r[i];
        int r1 = L->rowptr[i + 1] - 1;
#if !(defined(MK_MUTATE) && MK_MUTATE == 103) /* planted fault 103: the forward solve without L's off-diagonal */
        for (int e = L->rowptr[i]; e < r1; e++) s -= L->val[e] * z[L->col[e]];
#endif
        z[i] = s / L->val[r1];
    }
    for (int i = n - 1; i >= 0; i--) {
        int r1 = L->rowptr[i + 1] - 1;
        z[i] /= L->val[r1];
        for (int e = L->rowptr[i]; e < r1; e++) z[L->col[e]] -= L->val[e] * z[i];
    }
}

typedef struct { int kind; double *d; mk_ic0_ ic; } mk_pc_;
static inline int mk_pc_build_(mk_pc_ *P, const mk_csr *A, int kind)
{
    int n = A->n;
    P->kind = kind, P->d = NULL;
    if (kind < MK_PC_NONE || kind > MK_PC_IC0) { P->kind = MK_PC_NONE; return -1; }
    if (kind == MK_PC_NONE) return 0;
    if (kind == MK_PC_IC0) return mk_ic0_build_(&P->ic, A);
    P->d = malloc(sizeof(double) * (kind == MK_PC_BJACOBI3 ? 3 * (size_t)n : (size_t)n));
    if (!P->d) return -1;
    if (kind == MK_PC_JACOBI) {
        for (int i = 0; i < n; i++) {
            double a = 0;
            for (int e = A->rowptr[i]; e < A->rowptr[i + 1]; e++) if (A->col[e] == i) a = A->val[e];
            P->d[i] = a != 0 ? 1 / a : 1;
        }
        return 0;
    }
    if (n % 3) { free(P->d); return -1; }
    for (int b = 0; b < n / 3; b++) { /* the inverse of each 3x3 diagonal block (0 where it is singular) */
        double m[9] = {0};
        for (int a = 0; a < 3; a++) {
            int i = 3 * b + a;
            for (int e = A->rowptr[i]; e < A->rowptr[i + 1]; e++)
                if (A->col[e] / 3 == b) m[3 * a + A->col[e] % 3] = A->val[e];
        }
        double mx = 0; /* invert 2^-e M (largest entry in [0.5, 1)): the 1e-300 test is then relative to |M|^3 */
        int ex = 0;
        for (int t = 0; t < 9; t++) mx = fmax(mx, fabs(m[t]));
#if !(defined(MK_MUTATE) && MK_MUTATE == 106)
        if (mx > 0 && isfinite(mx)) frexp(mx, &ex); /* planted fault 106 skips this, so the test is absolute */
#endif
        for (int t = 0; t < 9; t++) m[t] = ldexp(m[t], -ex);
        double a0 = m[0], b0 = m[1], c0 = m[2], d0 = m[3], e0 = m[4], f0 = m[5], g0 = m[6], h0 = m[7], k0 = m[8];
        double AA = e0 * k0 - f0 * h0, BB = -(d0 * k0 - f0 * g0), CC = d0 * h0 - e0 * g0, det = a0 * AA + b0 * BB + c0 * CC;
        double s = fabs(det) > 1e-300 ? 1.0 / det : 0.0, *o = P->d + 9L * b;
        o[0] = AA * s, o[1] = -(b0 * k0 - c0 * h0) * s, o[2] = (b0 * f0 - c0 * e0) * s;
        o[3] = BB * s, o[4] = (a0 * k0 - c0 * g0) * s, o[5] = -(a0 * f0 - c0 * d0) * s;
        o[6] = CC * s, o[7] = -(a0 * h0 - b0 * g0) * s, o[8] = (a0 * e0 - b0 * d0) * s;
        for (int t = 0; t < 9; t++) o[t] = ldexp(o[t], -ex);
    }
    return 0;
}
/* Jacobi on unknowns [i0, i1), or block-Jacobi on 3x3 blocks [i0, i1) */
static inline void mk_pc_apply_blk_(const mk_pc_ *P, const double *r, double *z, long i0, long i1)
{
    if (P->kind == MK_PC_JACOBI) for (long i = i0; i < i1; i++) z[i] = P->d[i] * r[i];
    else
        for (long b = i0; b < i1; b++) {
            const double *o = P->d + 9 * b, *x = r + 3 * b;
            for (int a = 0; a < 3; a++) z[3 * b + a] = o[3 * a] * x[0] + o[3 * a + 1] * x[1] + o[3 * a + 2] * x[2];
        }
}
static inline void mk_pc_apply_(const mk_pc_ *P, int n, const double *r, double *z)
{
    if (P->kind == MK_PC_NONE) memcpy(z, r, sizeof *z * (size_t)n);
    else if (P->kind == MK_PC_IC0) mk_ic0_apply_(&P->ic, r, z);
    else MK_PAR_FOR_(P->kind == MK_PC_JACOBI ? n : n / 3, i0, i1, mk_pc_apply_blk_(P, r, z, i0, i1));
}
static inline void mk_pc_free_(mk_pc_ *P)
{
    if (P->kind == MK_PC_IC0) mk_ic0_free_(&P->ic);
    free(P->d);
}
/* the serial pieces the threads run on [i0, i1) */
static inline double mk_dot_blk_(const double *a, const double *b, long i0, long i1)
{
    double s = 0;
    for (long i = i0; i < i1; i++) s += a[i] * b[i];
    return s;
}
static inline double mk_sum_blk_(const double *a, long i0, long i1)
{
    double s = 0;
    for (long i = i0; i < i1; i++) s += a[i];
    return s;
}
static inline double mk_absmax_blk_(const double *a, long i0, long i1)
{
    double m = 0;
    for (long i = i0; i < i1; i++) m = fmax(m, fabs(a[i]));
    return m;
}
static inline void mk_unit_blk_(const double *b, double *r, double *u, int e, long i0, long i1)
{
    for (long i = i0; i < i1; i++) r[i] = ldexp(r[i], -e), u[i] = ldexp(b[i], -e);
}
static inline void mk_shift_blk_(double *r, double m, long i0, long i1) { for (long i = i0; i < i1; i++) r[i] -= m; }
static inline void mk_sub_blk_(const double *b, const double *q, double *r, long i0, long i1) { for (long i = i0; i < i1; i++) r[i] = b[i] - q[i]; }
/* x += c p, r -= a q (CG's update) */
static inline void mk_cg_xr_blk_(double *x, double *r, const double *p, const double *q, double c, double a, long i0, long i1)
{
    for (long i = i0; i < i1; i++) x[i] += c * p[i], r[i] -= a * q[i];
}
/* p = z + beta p */
static inline void mk_xpby_blk_(const double *z, double beta, double *p, long i0, long i1) { for (long i = i0; i < i1; i++) p[i] = z[i] + beta * p[i]; }
/* CG's passes that also give a block's partial of a dot, so one pass over the block does both: q = A p with p.q; the
   update x += c p, r -= a q with r.r; Jacobi's z = P r with r.z. Each partial is the one mk_dot_blk_ gives on the
   block. */
static inline double mk_cg_aq_blk_(const mk_csr *A, const double *p, double *q, long i0, long i1)
{
    mk_csr_matvec_rows_(A, p, q, i0, i1);
    return mk_dot_blk_(p, q, i0, i1);
}
static inline double mk_cg_xrr_blk_(double *x, double *r, const double *p, const double *q, double c, double a, long i0, long i1)
{
    mk_cg_xr_blk_(x, r, p, q, c, a, i0, i1);
    return mk_dot_blk_(r, r, i0, i1);
}
static inline double mk_cg_pz_blk_(const mk_pc_ *P, const double *r, double *z, long i0, long i1)
{
    mk_pc_apply_blk_(P, r, z, i0, i1);
    return mk_dot_blk_(r, z, i0, i1);
}
/* the dot of a and b: blocks' partials added in block order (mk_par.h), the same at any thread count */
static inline double mk_dot_(const double *a, const double *b, int n)
{
    double s;
    MK_PAR_SUM_(n, i0, i1, mk_dot_blk_(a, b, i0, i1), s);
    return s;
}
/* the solvers work in units of 2^e, e the exponent of the largest |b_i| (0 for b = 0, at most 1023), so their dots
   neither underflow nor overflow at extreme scales of b: r is scaled in place, u gets the scaled b, the return is 2^e.
   The largest is taken block by block (a maximum does not depend on the order). */
static inline double mk_unit_(const double *b, double *r, int n, double *u)
{
    double big = 0, part[MK_PAR_NBLK];
    int e = 0;
    const long nb = mk_par_nblk(n);
    MK_OMP_(omp parallel for schedule(static) if(nb > 1))
    for (long k = 0; k < nb; k++) {
        long i0, i1;
        mk_par_range(n, k, &i0, &i1), part[k] = mk_absmax_blk_(b, i0, i1);
    }
    for (long k = 0; k < nb; k++) big = fmax(big, part[k]);
    if (big > 0 && isfinite(big)) frexp(big, &e);
    if (e > 1023) e = 1023;
    MK_PAR_FOR_(n, i0, i1, mk_unit_blk_(b, r, u, e, i0, i1));
#if defined(MK_MUTATE) && MK_MUTATE == 105
    for (int i = 0; i < n; i++) r[i] = ldexp(r[i], e), u[i] = b[i]; /* planted fault: unscaled */
    e = 0;
#endif
    return ldexp(1.0, e);
}
static inline void mk_meanfree_(double *r, int n)
{
    double m;
    MK_PAR_SUM_(n, i0, i1, mk_sum_blk_(r, i0, i1), m);
    m /= n;
    MK_PAR_FOR_(n, i0, i1, mk_shift_blk_(r, m, i0, i1));
}

/* (preconditioned) conjugate gradients for symmetric positive-definite A. precond: MK_PC_NONE, _JACOBI, _BJACOBI3 or
   _IC0, plus MK_PC_MEANFREE for a singular A whose null space is the constants (periodic or all-Neumann Poisson): the
   residual is kept mean-free; any other bit is refused (-1). 0 and 1 report the true residual b - A x in relres: when
   the recurrence's residual reaches tol, b - A x is formed again, CG restarts from x while that is still past tol, and
   stops with 1 when a restart did not lower it (tol under the attainable accuracy). *iters counts the iterations of all
   passes. On -1 relres is the last the iteration had. A negative r z (an indefinite preconditioner) is -1. */
static inline int mk_cg(const mk_csr *A, const double *b, double *x, int precond, double tol, int maxit, int *iters, double *relres)
{
    int n = A->n, meanfree = precond & MK_PC_MEANFREE, kind = precond & 15, it = 0, ret = 1, last = 0;
    if (iters) *iters = 0;
    if (relres) *relres = NAN;
    mk_pc_ P = {0};
    double *w = malloc(sizeof(double) * 4 * (size_t)n);
    /* bits past MK_PC_MEANFREE refused */
    if (!w || (precond & ~31) || mk_pc_build_(&P, A, kind)) { free(w); return -1; }
    double *r = w, *z = r + n, *p = z + n, *q = p + n, sc, bn, rn, rnp = 0;
    /* a pass: the true residual, then CG from x until the recurrence's reaches tol, maxit, a stall or a breakdown */
    for (int pass = 0;; pass++) {
        mk_csr_matvec(A, x, q);
        MK_PAR_FOR_(n, i0, i1, mk_sub_blk_(b, q, r, i0, i1));
        if (meanfree) mk_meanfree_(r, n);
        sc = mk_unit_(b, r, n, q), bn = sqrt(mk_dot_(q, q, n)), rn = sqrt(mk_dot_(r, r, n));
        if (bn == 0) bn = 1;
        if (rn / bn <= tol) { ret = 0; break; }
        /* a restart that did not lower b - A x ends it */
        if (last || it >= maxit || (pass && rn >= rnp)) { ret = 1; break; }
        rnp = rn;
        mk_pc_apply_(&P, n, r, z);
        if (meanfree) mk_meanfree_(z, n);
        memcpy(p, z, sizeof *p * (size_t)n);
        double rz = mk_dot_(r, z, n), rz0 = rz;
        const int it0 = it;
        ret = 1;
        while (ret == 1 && it < maxit) {
            double pq, rr, rz2;
            MK_PAR_SUM_(n, i0, i1, mk_cg_aq_blk_(A, p, q, i0, i1), pq); /* q = A p and p.q in one pass */
            if (!(pq > 0)) { ret = -1; break; } /* lost positive-definiteness, or NaN */
            double a = rz / pq;
            /* x += c p, r -= a q, then r.r; without MK_PC_MEANFREE the update and r.r in one pass */
            if (meanfree) {
                MK_PAR_FOR_(n, i0, i1, mk_cg_xr_blk_(x, r, p, q, a * sc, a, i0, i1));
                mk_meanfree_(r, n);
                rr = mk_dot_(r, r, n);
            } else MK_PAR_SUM_(n, i0, i1, mk_cg_xrr_blk_(x, r, p, q, a * sc, a, i0, i1), rr);
            rn = sqrt(rr);
            it++;
            if (rn / bn <= tol) { ret = 0; break; }
            const double *zz = z;
            /* z and r.z: with no preconditioner z = r, so r.z is r.r bit for bit; Jacobi forms z = P r and r.z in one
               pass */
            if (kind == MK_PC_NONE && !meanfree) zz = r, rz2 = rr;
            else if (kind == MK_PC_JACOBI && !meanfree) MK_PAR_SUM_(n, i0, i1, mk_cg_pz_blk_(&P, r, z, i0, i1), rz2);
            else {
                mk_pc_apply_(&P, n, r, z);
                if (meanfree) mk_meanfree_(z, n);
                rz2 = mk_dot_(r, z, n);
            }
            /* r z negative or NaN: a breakdown. r z zero, or underflowed after falling 2^-100 below its start (the
               residual about 1e-15 of its first: past attainable accuracy, as with tol 0): a stall. An underflow alone
               is not a stop: a preconditioner near 1e290 starts r z near 1e-290, and those solves still converge. */
            if (!(rz2 >= 0)) { ret = -1; break; }
            if (!(rz2 > 0) || (rz2 < DBL_MIN && rz2 < rz0 * 0x1p-100)) break;
#if defined(MK_MUTATE) && MK_MUTATE == 101
            double beta = 0; (void)rz2; /* planted fault: beta 0, steepest descent */
#else
            double beta = rz2 / rz;
#endif
            MK_PAR_FOR_(n, i0, i1, mk_xpby_blk_(zz, beta, p, i0, i1));
            rz = rz2;
        }
        /* a restart that broke down before its first step (r z or p A p underflowing at an extreme scale of A) leaves x
           as the last pass had it, and rn its true residual: not within tol (1), not a breakdown */
        if (ret < 0) {
            if (pass && it == it0) ret = 1;
            break;
        }
        last = ret; /* 1 (maxit or a stall): the next pass reads the true residual and stops */
    }
#if !(defined(MK_MUTATE) && MK_MUTATE == 104) /* planted fault 104: a non-finite x not refused */
    if (ret >= 0) for (int i = 0; i < n; i++) if (!isfinite(x[i])) { ret = -1; break; }
#endif
    if (iters) *iters = it;
    if (relres) *relres = rn / bn;
    mk_pc_free_(&P), free(w);
    return ret;
}

/* BiCGSTAB's vector steps on [i0, i1) */
static inline void mk_bi_start_blk_(const double *r, double *r0, double *p, double *v, long i0, long i1)
{
    for (long i = i0; i < i1; i++) r0[i] = r[i], p[i] = v[i] = 0;
}
static inline void mk_bi_p_blk_(const double *r, const double *v, double beta, double om, double *p, long i0, long i1)
{
    for (long i = i0; i < i1; i++) p[i] = r[i] + beta * (p[i] - om * v[i]);
}
static inline void mk_bi_s_blk_(const double *r, const double *v, double alpha, double *s, long i0, long i1)
{
    for (long i = i0; i < i1; i++) s[i] = r[i] - alpha * v[i];
}
static inline void mk_bi_x_blk_(double *x, const double *ph, double c, long i0, long i1) { for (long i = i0; i < i1; i++) x[i] += c * ph[i]; }
static inline void mk_bi_xr_blk_(double *x, double *r, const double *ph, const double *sh, const double *s, const double *t, double ca, double co, double om,
                                 long i0, long i1)
{
    for (long i = i0; i < i1; i++) x[i] += ca * ph[i] + co * sh[i], r[i] = s[i] - om * t[i];
}
/* BiCGSTAB (van der Vorst 1992), right-preconditioned (MK_PC_NONE, _JACOBI or _BJACOBI3) for a general A. As mk_cg: 0
   and 1 with the true residual, restarts from x while it is past tol and falling. MK_PC_IC0, MK_PC_MEANFREE and any
   other bit are refused (-1). */
static inline int mk_bicgstab(const mk_csr *A, const double *b, double *x, int precond, double tol, int maxit, int *iters, double *relres)
{
    int n = A->n, it = 0, ret = 1, last = 0;
    if (iters) *iters = 0;
    if (relres) *relres = NAN;
    mk_pc_ P = {0};
    if ((precond & 15) == MK_PC_IC0 || (precond & ~15)) return -1; /* IC(0), MK_PC_MEANFREE and any higher bit */
    double *w = malloc(sizeof(double) * 8 * (size_t)n);
    if (!w || mk_pc_build_(&P, A, precond & 15)) { free(w); return -1; }
    double *r = w, *r0 = r + n, *p = r0 + n, *v = p + n, *s = v + n, *t = s + n, *ph = t + n, *sh = ph + n, sc, bn, rn, rnp = 0;
    for (int pass = 0;; pass++) { /* a pass, as mk_cg's */
        mk_csr_matvec(A, x, v);
        MK_PAR_FOR_(n, i0, i1, mk_sub_blk_(b, v, r, i0, i1));
        sc = mk_unit_(b, r, n, s), bn = sqrt(mk_dot_(s, s, n)), rn = sqrt(mk_dot_(r, r, n));
        if (bn == 0) bn = 1;
        if (rn / bn <= tol) { ret = 0; break; }
        if (last || it >= maxit || (pass && rn >= rnp)) { ret = 1; break; } /* as mk_cg's */
        rnp = rn;
        double rho = 1, alpha = 1, om = 1;
        MK_PAR_FOR_(n, i0, i1, mk_bi_start_blk_(r, r0, p, v, i0, i1));
        const int it0 = it;
        ret = 1;
        while (ret == 1 && it < maxit) {
            double rho2 = mk_dot_(r0, r, n);
            if (rho2 == 0 || om == 0) { ret = -1; break; }
            double beta = rho2 / rho * alpha / om;
            rho = rho2;
            MK_PAR_FOR_(n, i0, i1, mk_bi_p_blk_(r, v, beta, om, p, i0, i1));
            mk_pc_apply_(&P, n, p, ph);
            mk_csr_matvec(A, ph, v);
            double r0v = mk_dot_(r0, v, n);
            if (r0v == 0) { ret = -1; break; }
            alpha = rho / r0v;
            MK_PAR_FOR_(n, i0, i1, mk_bi_s_blk_(r, v, alpha, s, i0, i1));
            it++;
            double sn = sqrt(mk_dot_(s, s, n));
            if (sn / bn <= tol) {
                MK_PAR_FOR_(n, i0, i1, mk_bi_x_blk_(x, ph, alpha * sc, i0, i1));
                rn = sn, ret = 0;
                break;
            }
            mk_pc_apply_(&P, n, s, sh);
            mk_csr_matvec(A, sh, t);
            double tt = mk_dot_(t, t, n);
            om = tt > 0 ? mk_dot_(t, s, n) / tt : 0;
#if defined(MK_MUTATE) && MK_MUTATE == 102
            /* planted fault: the omega step left out of x */
            for (int i = 0; i < n; i++) x[i] += alpha * ph[i], r[i] = s[i] - om * t[i];
#else
            MK_PAR_FOR_(n, i0, i1, mk_bi_xr_blk_(x, r, ph, sh, s, t, alpha * sc, om * sc, om, i0, i1));
#endif
            rn = sqrt(mk_dot_(r, r, n));
            if (!isfinite(rn)) { ret = -1; break; }
            if (rn / bn <= tol) ret = 0;
        }
        if (ret < 0) { /* as mk_cg's */
            if (pass && it == it0) ret = 1;
            break;
        }
        last = ret;
    }
#if !(defined(MK_MUTATE) && MK_MUTATE == 104)
    if (ret >= 0) for (int i = 0; i < n; i++) if (!isfinite(x[i])) { ret = -1; break; }
#endif
    if (iters) *iters = it;
    if (relres) *relres = rn / bn;
    mk_pc_free_(&P), free(w);
    return ret;
}
#endif
