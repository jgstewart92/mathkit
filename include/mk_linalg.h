/* mathkit: small dense linear algebra, row-major double a[i * n + j]. LU with partial pivoting, Cholesky, Householder
   least squares that survives rank deficiency, symmetric decomposition by Jacobi rotations (values ascending, vectors
   in columns), a tridiagonal (Thomas) solve and 3x3 fast paths. Every solver returns -1 rather than divide by a
   vanishing or NaN pivot. Header only, C23. */
#ifndef MK_LINALG_H
#define MK_LINALG_H
#include <limits.h>
#include <math.h>
#include <string.h>
_Static_assert(sizeof(long) == 8, "mk_linalg.h needs LP64: it counts m * n entries in a long");

/* LU in place (L's unit diagonal not stored) with row pivots piv[n]. -1 when a pivot is below 1e-13 of the largest
   entry, or NaN, or n is outside 0..46340. */
static inline int mk_lu(double *a, int n, int *piv)
{
    if (n < 0 || n > 46340) return -1; /* n * n must fit an int, and a negative n has no matrix */
    double scale = 0;
    for (int i = 0; i < n * n; i++) scale = fmax(scale, fabs(a[i]));
    for (int k = 0; k < n; k++) {
        int p = k;
#if !(defined(MK_MUTATE) && MK_MUTATE == 21) /* planted fault 21: no pivot search */
        for (int i = k + 1; i < n; i++)
            if (fabs(a[i * n + k]) > fabs(a[p * n + k])) p = i;
#endif
        piv[k] = p;
        if (!(fabs(a[p * n + k]) > 1e-13 * scale)) return -1;
        if (p != k)
            for (int j = 0; j < n; j++) { double t = a[k * n + j]; a[k * n + j] = a[p * n + j]; a[p * n + j] = t; }
        for (int i = k + 1; i < n; i++) {
            double f = a[i * n + k] /= a[k * n + k];
            for (int j = k + 1; j < n; j++) a[i * n + j] -= f * a[k * n + j];
        }
    }
    return 0;
}
/* 0, or -1 if the answer is not finite: LU's refusal is relative to the largest entry by design, so a well-conditioned
   but tiny matrix (1e-310 I) passes it and the answer can overflow; the solve reports that rather than returning it
   quietly */
static inline int mk_lu_solve(const double *lu, int n, const int *piv, double *b)
{
    for (int k = 0; k < n; k++)
        if (piv[k] != k) { double t = b[k]; b[k] = b[piv[k]]; b[piv[k]] = t; }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++) b[i] -= lu[i * n + j] * b[j];
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) b[i] -= lu[i * n + j] * b[j];
        b[i] /= lu[i * n + i];
    }
    int ok = 1;
    for (int i = 0; i < n; i++) ok &= isfinite(b[i]) != 0;
    return ok ? 0 : -1;
}
static inline double mk_lu_det(const double *lu, int n, const int *piv)
{
    double d = 1; /* mantissa and exponent kept apart, so a long product overflows or underflows only at the end */
    int e = 0, k2;
    for (int k = 0; k < n; k++) d *= frexp(piv[k] != k ? -lu[k * n + k] : lu[k * n + k], &k2), e += k2, d = frexp(d, &k2), e += k2;
    return ldexp(d, e);
}

/* Cholesky, lower triangle in place (the upper is left as it was). -1 unless positive definite. */
static inline int mk_chol(double *a, int n)
{
    if (n < 0 || n > 46340) return -1; /* j * n + k must fit an int */
    for (int j = 0; j < n; j++) {
        double d = a[j * n + j];
        for (int k = 0; k < j; k++) d -= a[j * n + k] * a[j * n + k];
        if (!(d > 0) || isinf(d)) return -1; /* NaN, not positive, or an infinite diagonal */
#if !(defined(MK_MUTATE) && MK_MUTATE == 22) /* planted fault 22: no square root */
        d = sqrt(d);
#endif
        a[j * n + j] = d;
        for (int i = j + 1; i < n; i++) {
            double s = a[i * n + j];
            for (int k = 0; k < j; k++) s -= a[i * n + k] * a[j * n + k];
            a[i * n + j] = s / d;
        }
    }
    return 0;
}
/* 0, or -1 if the answer is not finite (as mk_lu_solve) */
static inline int mk_chol_solve(const double *l, int n, double *b)
{
    if (n < 0 || n > 46340) return -1;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < i; k++) b[i] -= l[i * n + k] * b[k];
        b[i] /= l[i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int k = i + 1; k < n; k++) b[i] -= l[k * n + i] * b[k];
        b[i] /= l[i * n + i];
    }
    int ok = 1;
    for (int i = 0; i < n; i++) ok &= isfinite(b[i]) != 0;
    return ok ? 0 : -1;
}

/* Least squares min |a x - b|, a m x n (m >= n, n <= MK_LSTSQ_MAX), Householder QR with column scaling and pivoting. a
   and b are destroyed. A column whose pivot falls below 1e-10 of its scaled norm is dropped and its unknown set to 0,
   so a rank-deficient system still solves. Returns the residual norm, or -1 for a size out of range or a NaN or
   infinite entry in a. */
#define MK_LSTSQ_MAX 64
static inline double mk_lstsq(double *a, int m, int n, double *b, double *x)
{
    if (n < 1 || n > MK_LSTSQ_MAX || m < n || m > INT_MAX / n) return -1; /* i * n + j must fit an int */
    for (long i = 0; i < (long)m * n; i++)
        if (!isfinite(a[i])) return -1;
    double colmax[MK_LSTSQ_MAX], colnorm[MK_LSTSQ_MAX], diag[MK_LSTSQ_MAX];
    int perm[MK_LSTSQ_MAX], rk = 0;
    /* the column's norm scaled by its largest entry, so neither overflows nor underflows */
    for (int j = 0; j < n; j++) {
        double mx = 0, s = 0;
        for (int i = 0; i < m; i++) mx = fmax(mx, fabs(a[i * n + j]));
        if (!(mx > 0 && isfinite(mx))) mx = 1;
        for (int i = 0; i < m; i++) { double t = a[i * n + j] / mx; s += t * t; }
        colmax[j] = mx, colnorm[j] = s > 0 ? sqrt(s) : 1.0;
        for (int i = 0; i < m; i++) a[i * n + j] = a[i * n + j] / mx / colnorm[j];
        perm[j] = j;
    }
    /* Householder QR with column pivoting: each step takes the remaining column farthest from the span of those taken,
       and stops when that distance is below 1e-10 (columns have unit norm): the rest are dropped, unknowns 0 */
    for (int t = 0; t < n; t++) {
        int best = -1;
        double bn = 0;
        for (int q = t; q < n; q++) {
            double s = 0;
            for (int i = t; i < m; i++) s += a[i * n + perm[q]] * a[i * n + perm[q]];
            if (s > bn) bn = s, best = q;
        }
        if (best < 0 || !(sqrt(bn) > 1e-10)) break;
        int k = perm[best];
        perm[best] = perm[t], perm[t] = k;
        double norm = sqrt(bn), alpha = a[t * n + k] > 0 ? -norm : norm;
        double v0 = a[t * n + k] - alpha;
        a[t * n + k] = v0;
        double vnorm2 = v0 * v0;
        for (int i = t + 1; i < m; i++) vnorm2 += a[i * n + k] * a[i * n + k];
        diag[t] = alpha;
        for (int q = t + 1; q < n; q++) {
            int j = perm[q];
            double s = 0;
            for (int i = t; i < m; i++) s += a[i * n + k] * a[i * n + j];
            s *= 2 / vnorm2;
            for (int i = t; i < m; i++) a[i * n + j] -= s * a[i * n + k];
        }
        double s = 0;
        for (int i = t; i < m; i++) s += a[i * n + k] * b[i];
        s *= 2 / vnorm2;
        for (int i = t; i < m; i++) b[i] -= s * a[i * n + k];
        rk = t + 1;
    }
    for (int j = 0; j < n; j++) x[j] = 0;
    for (int t = rk - 1; t >= 0; t--) {
        double s = b[t];
        for (int q = t + 1; q < rk; q++) s -= a[t * n + perm[q]] * x[perm[q]];
        x[perm[t]] = s / diag[t];
    }
    double r = 0;
    for (int i = rk; i < m; i++) r += b[i] * b[i];
    (void)colmax, (void)colnorm;
#if !(defined(MK_MUTATE) && MK_MUTATE == 23) /* planted fault 23: the column scaling not undone */
    for (int j = 0; j < n; j++) x[j] = x[j] / colnorm[j] / colmax[j];
#endif
    return sqrt(r);
}

/* Symmetric a (n x n, destroyed): a = V diag(values) V^T, values ascending, vectors[i * n + k] the k-th column.
   Cyclic Jacobi with the stable tangent; stops when the off-diagonal is below 1e-13 of the whole (in norm). Returns the
   sweeps used, or -1 if 100 sweeps did not converge. */
static inline int mk_sym_jacobi(double *a, int n, double *vectors, double *values)
{
    if (n < 0 || n > 46340) return -1; /* n * n must fit an int */
    int sweep, done = 0;
    int e = 0, fin = 1; /* scale once by a power of two (exact) so the largest entry is near 1: subnormal and huge
                           entries otherwise round the rotations to nothing or past the largest double */
    double m0 = 0;
    for (int i = 0; i < n * n; i++) fin &= isfinite(a[i]) != 0, m0 = fmax(m0, fabs(a[i]));
    if (fin && m0 > 0) { e = ilogb(m0); for (int i = 0; i < n * n; i++) a[i] = ldexp(a[i], -e); }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) vectors[i * n + j] = i == j;
    for (sweep = 0; sweep < 100; sweep++) {
        double off = 0, tot = 0, mx = 0; /* sums of squares scaled by the largest entry: no overflow or underflow */
        int finite = 1;
        for (int i = 0; i < n * n; i++) finite &= isfinite(a[i]) != 0, mx = fmax(mx, fabs(a[i]));
        if (!finite) break;              /* a NaN or infinite entry: not converged (-1) */
        if (mx == 0) { done = 1; break; }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                double t = a[i * n + j] / mx;
                tot += t * t;
                if (i != j) off += t * t;
            }
        if (off <= 1e-26 * tot) { done = 1; break; }
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++) {
                double apq = a[p * n + q];
                if (apq == 0) continue;
                double theta = (a[q * n + q] - a[p * n + p]) / (2 * apq);
                double t = (theta >= 0 ? 1 : -1) / (fabs(theta) + sqrt(theta * theta + 1));
#if defined(MK_MUTATE) && MK_MUTATE == 24
                t = -t; /* planted fault: the rotation's sign */
#endif
                double c = 1 / sqrt(t * t + 1), s = t * c;
                for (int k = 0; k < n; k++) {
                    double akp = a[k * n + p], akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq, a[k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; k++) {
                    double apk = a[p * n + k], aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk, a[q * n + k] = s * apk + c * aqk;
                }
                for (int k = 0; k < n; k++) {
                    double vkp = vectors[k * n + p], vkq = vectors[k * n + q];
                    vectors[k * n + p] = c * vkp - s * vkq, vectors[k * n + q] = s * vkp + c * vkq;
                }
            }
    }
    for (int i = 0; i < n; i++) values[i] = ldexp(a[i * n + i], e);
    for (int i = 1; i < n; i++) /* insertion sort, columns follow */
        for (int j = i; j > 0 && values[j] < values[j - 1]; j--) {
            double t = values[j]; values[j] = values[j - 1], values[j - 1] = t;
            for (int k = 0; k < n; k++) { t = vectors[k * n + j]; vectors[k * n + j] = vectors[k * n + j - 1], vectors[k * n + j - 1] = t; }
        }
    return done ? sweep : -1;
}

/* Tridiagonal solve (Thomas): a[i] x[i-1] + b[i] x[i] + c[i] x[i+1] = d[i] for i = 0..n-1 (a[0] and c[n-1] are not
   read). No pivoting; b and d are overwritten (b with the pivots) and x may be d. 0, or -1 for n < 1, a pivot that is
   NaN, zero or under 2^-52 of its row's |a| + |b| + |c| (singular to rounding, or the system needs pivoting: use
   mk_lu), or an answer that is not finite. Diagonally dominant systems, as heat columns give, need no pivoting. */
static inline int mk_tridiag(int n, const double *a, double *b, const double *c, double *d, double *x)
{
    if (n < 1) return -1;
#if defined(MK_MUTATE) && MK_MUTATE == 27
#define MK_TRIDIAG_SMALL_(p, row) (!((p) != 0) || !isfinite(row)) /* planted fault: only an exact zero refused */
#else
#define MK_TRIDIAG_SMALL_(p, row) (!(fabs(p) > 0x1p-52 * (row)))
#endif
    if (MK_TRIDIAG_SMALL_(b[0], fabs(b[0]) + (n > 1 ? fabs(c[0]) : 0))) return -1;
    for (int i = 1; i < n; i++) {
        double row = fabs(a[i]) + fabs(b[i]) + (i + 1 < n ? fabs(c[i]) : 0), f = a[i] / b[i - 1];
#if defined(MK_MUTATE) && MK_MUTATE == 26
        d[i] -= f * d[i - 1]; /* planted fault: the pivot left as b[i] */
#else
        b[i] -= f * c[i - 1], d[i] -= f * d[i - 1];
#endif
        if (MK_TRIDIAG_SMALL_(b[i], row)) return -1;
    }
#undef MK_TRIDIAG_SMALL_
    x[n - 1] = d[n - 1] / b[n - 1];
    int ok = isfinite(x[n - 1]);
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
        ok &= isfinite(x[i]) != 0;
    }
    return ok ? 0 : -1;
}

/* 3x3: Gaussian elimination with partial pivoting; 0, or -1 for a zero pivot or a non-finite answer. A pivot small from
   rounding is not refused (mk_lu's 1e-13 rule is): [[1 2 3] [4 5 6] [7 8 9]] gives 0 with x near 9e15 for
   b (1, 0, 0). */
static inline int mk_solve3(const double A[3][3], const double b[3], double x[3])
{
    double m[3][4];
    for (int r = 0; r < 3; r++) { for (int c = 0; c < 3; c++) m[r][c] = A[r][c]; m[r][3] = b[r]; }
    for (int c = 0; c < 3; c++) {
        int p = c;
#if !(defined(MK_MUTATE) && MK_MUTATE == 25) /* planted fault 25: no pivot search */
        for (int r = c + 1; r < 3; r++) if (fabs(m[r][c]) > fabs(m[p][c])) p = r;
#endif
        if (!(fabs(m[p][c]) > 0)) return -1;
        for (int k = 0; k < 4; k++) { double t = m[c][k]; m[c][k] = m[p][k]; m[p][k] = t; }
        for (int r = c + 1; r < 3; r++) { double f = m[r][c] / m[c][c]; for (int k = c; k < 4; k++) m[r][k] -= f * m[c][k]; }
    }
    for (int r = 2; r >= 0; r--) { double v = m[r][3]; for (int c = r + 1; c < 3; c++) v -= m[r][c] * x[c]; x[r] = v / m[r][r]; }
    return isfinite(x[0]) && isfinite(x[1]) && isfinite(x[2]) ? 0 : -1;
}
static inline double mk_det3(const double A[3][3])
{
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
           A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}
/* symmetric 3x3 (A not changed): values ascending, vec[k] the k-th unit vector, right-handed */
static inline int mk_sym3(const double A[3][3], double val[3], double vec[3][3])
{
    double a[9], v[9];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) a[i * 3 + j] = A[i][j];
    int r = mk_sym_jacobi(a, 3, v, val);
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < 3; i++) vec[k][i] = v[i * 3 + k];
    double c[3] = {vec[0][1] * vec[1][2] - vec[0][2] * vec[1][1], vec[0][2] * vec[1][0] - vec[0][0] * vec[1][2],
                   vec[0][0] * vec[1][1] - vec[0][1] * vec[1][0]};
    if (c[0] * vec[2][0] + c[1] * vec[2][1] + c[2] * vec[2][2] < 0)
        for (int i = 0; i < 3; i++) vec[2][i] = -vec[2][i];
    return r;
}
#endif
