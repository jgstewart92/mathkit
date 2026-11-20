/* mathkit: numerical integration. Gauss-Legendre nodes and weights for any n (Newton on the Legendre recurrence, exact for
   polynomials to degree 2n - 1); adaptive Gauss-Kronrod 7-15 with an error estimate, which also takes integrable endpoint
   singularities (nodes never touch the ends); rules on the reference triangle (0,0) (1,0) (0,1) and tetrahedron
   (0,0,0) (1,0,0) (0,1,0) (0,0,1), weights summing to the area 1/2 and the volume 1/6. Header only, C23. */
#ifndef MK_QUAD_H
#define MK_QUAD_H
#include <math.h>
#include <stdlib.h>

#define MK_QUAD_PI 3.14159265358979323846

/* n nodes x[] (ascending) and weights w[] on [-1, 1]; 0, or -1 if Newton did not settle or n is outside 1..2^30. O(n^2): n 65536
   takes 23 s (MEASURED 1.6.1, the development machine); every n from 21 to 2000 and 2^11 to 2^16 settles, weights summing to 2 within 1.8e-15. */
static inline int mk_gauss_legendre(int n, double *x, double *w)
{
    if (n < 1 || n > (1 << 30)) return -1; /* 2 * k - 1 and (n + 1) / 2 stay inside an int */
    for (int i = 0; i < (n + 1) / 2; i++) {
        double z = cos(MK_QUAD_PI * (i + 0.75) / (n + 0.5)), dp = 0;
        int ok = 0;
        for (int it = 0; it < 100; it++) {
            double p0 = 1, p1 = z;
            for (int k = 2; k <= n; k++) {
                double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
                p0 = p1, p1 = p2;
            }
            if (n == 1) p0 = 1, p1 = z;
            dp = n * (z * p1 - p0) / (z * z - 1); /* P_n' from P_n and P_{n-1} */
            double dz = p1 / dp;
            z -= dz;
            if (fabs(dz) < 1e-16) { ok = 1; break; }
        }
        if (!ok && fabs(z) > 1e-15) return -1;
        double p0 = 1, p1 = z;
        for (int k = 2; k <= n; k++) {
            double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
            p0 = p1, p1 = p2;
        }
        dp = n * (z * p1 - p0) / (z * z - 1);
        x[n - 1 - i] = z, x[i] = -z;
#if defined(MK_MUTATE) && MK_MUTATE == 111
        w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * dp * dp) * (1 + 1e-6);
#else
        w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * dp * dp);
#endif
    }
    if (n % 2) x[n / 2] = 0;
    return 0;
}

/* Kronrod 15 and its embedded Gauss 7 (QUADPACK qk15): nodes on [0, 1), symmetric */
static const double mk_gk15_x_[8] = {0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                                     0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                                     0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
                                     0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
static const double mk_gk15_wk_[8] = {0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                                      0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                                      0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                                      0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
static const double mk_gk15_wg_[4] = {0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                                      0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

typedef double (*mk_qfn)(double x, void *ctx);

/* one Gauss-Kronrod 7-15 pass on [a, b]: the Kronrod value, *err = |Kronrod - Gauss| */
static inline double mk_gk15(mk_qfn f, void *ctx, double a, double b, double *err)
{
    double c = 0.5 * (a + b), h = 0.5 * (b - a), fc = f(c, ctx);
    double k = fc * mk_gk15_wk_[7], g = fc * mk_gk15_wg_[3];
    for (int j = 0; j < 7; j++) {
        double dx = h * mk_gk15_x_[j], s = f(c - dx, ctx) + f(c + dx, ctx);
        k += mk_gk15_wk_[j] * s;
#if defined(MK_MUTATE) && MK_MUTATE == 112
        if (j % 2 == 1) g += mk_gk15_wg_[j / 2] * s * 1.001;
#else
        if (j % 2 == 1) g += mk_gk15_wg_[j / 2] * s;
#endif
    }
    if (err) *err = fabs((k - g) * h);
    return k * h;
}

/* adaptive: bisect the interval with the largest error until the sum of errors is under max(abstol, reltol |I|) or
   maxint intervals are in use. Returns the integral; *err the estimate; *status 0 met, 1 interval budget spent, -1 no memory
   or maxint >= 2^29 (refused before f is called or anything is allocated). */
static inline double mk_quad_adapt(mk_qfn f, void *ctx, double a, double b, double abstol, double reltol, int maxint, double *err, int *status)
{
    if (err) *err = NAN;
    if (maxint >= (1 << 29)) { if (status) *status = -1; return NAN; } /* 4 * i + 3 stays inside an int */
    if (maxint < 1) maxint = 1;
    double *iv = malloc(sizeof(double) * 4 * (size_t)maxint); /* a, b, value, error */
    if (!iv) { if (status) *status = -1; return NAN; }
    int m = 1;
    iv[0] = a, iv[1] = b, iv[2] = mk_gk15(f, ctx, a, b, &iv[3]);
    double tot, et;
    int st = 1;
    for (;;) {
        int w = 0; /* the sums are re-formed each pass (the search is O(m) anyway): running updates drift by cancellation */
        tot = 0, et = 0;
        for (int i = 0; i < m; i++) { tot += iv[4 * i + 2], et += iv[4 * i + 3]; if (iv[4 * i + 3] > iv[4 * w + 3]) w = i; }
        if (et <= fmax(abstol, reltol * fabs(tot))) { st = 0; break; }
        if (m >= maxint) break;
        double l = iv[4 * w], r = iv[4 * w + 1], c = 0.5 * (l + r), e1, e2;
        double v1 = mk_gk15(f, ctx, l, c, &e1), v2 = mk_gk15(f, ctx, c, r, &e2);
        iv[4 * w] = l, iv[4 * w + 1] = c, iv[4 * w + 2] = v1, iv[4 * w + 3] = e1;
        iv[4 * m] = c, iv[4 * m + 1] = r, iv[4 * m + 2] = v2, iv[4 * m + 3] = e2;
        m++;
    }
    tot = 0, et = 0; /* m may have grown by one since the last sums */
    for (int i = 0; i < m; i++) tot += iv[4 * i + 2], et += iv[4 * i + 3];
    free(iv);
    if (err) *err = et;
    if (status) *status = st;
    return tot;
}

/* triangle and tetrahedron rules: point i at barycentric-free coordinates p[i][..], weight w[i]; the count returned.
   degree 1 (centroid), 2 (3 points, as body-claude's bfe.c), 5 (Radon's 7 points); 3 and 4 get the 7, past 5 the count is 0 (no rule). */
static inline int mk_tri_rule(int degree, double p[][2], double w[])
{
    if (degree <= 1) { p[0][0] = p[0][1] = 1.0 / 3, w[0] = 0.5; return 1; }
    if (degree == 2) {
        const double q[3][2] = {{1.0 / 6, 1.0 / 6}, {2.0 / 3, 1.0 / 6}, {1.0 / 6, 2.0 / 3}};
        for (int i = 0; i < 3; i++) p[i][0] = q[i][0], p[i][1] = q[i][1], w[i] = 1.0 / 6;
#if defined(MK_MUTATE) && MK_MUTATE == 113
        p[1][0] = 0.6;
#endif
        return 3;
    }
    if (degree > 5) return 0;
    double s = sqrt(15.0), a = (6 - s) / 21, b = (6 + s) / 21, wa = (155 - s) / 2400, wb = (155 + s) / 2400;
    p[0][0] = p[0][1] = 1.0 / 3, w[0] = 9.0 / 80;
    double pa[3][2] = {{a, a}, {1 - 2 * a, a}, {a, 1 - 2 * a}}, pb[3][2] = {{b, b}, {1 - 2 * b, b}, {b, 1 - 2 * b}};
    for (int i = 0; i < 3; i++) {
        p[1 + i][0] = pa[i][0], p[1 + i][1] = pa[i][1], w[1 + i] = wa;
        p[4 + i][0] = pb[i][0], p[4 + i][1] = pb[i][1], w[4 + i] = wb;
    }
    return 7;
}
/* degree 1 (centroid), 2 (4 points, as bfe.c and fable), 3 (Keast's 5 points: one negative weight); past 3 the count is 0 */
static inline int mk_tet_rule(int degree, double p[][3], double w[])
{
    if (degree <= 1) { p[0][0] = p[0][1] = p[0][2] = 0.25, w[0] = 1.0 / 6; return 1; }
    if (degree == 2) {
#if defined(MK_MUTATE) && MK_MUTATE == 114
        const double a = 0.58, b = 0.14;
#else
        const double a = 0.5854101966249685, b = 0.1381966011250105;
#endif
        for (int g = 0; g < 4; g++) {
            double L[4] = {b, b, b, b};
            L[g] = a;
            p[g][0] = L[1], p[g][1] = L[2], p[g][2] = L[3], w[g] = 1.0 / 24;
        }
        return 4;
    }
    if (degree > 3) return 0;
    p[0][0] = p[0][1] = p[0][2] = 0.25, w[0] = -4.0 / 30;
    for (int g = 0; g < 4; g++) {
        double L[4] = {1.0 / 6, 1.0 / 6, 1.0 / 6, 1.0 / 6};
        L[g] = 0.5;
        p[1 + g][0] = L[1], p[1 + g][1] = L[2], p[1 + g][2] = L[3], w[1 + g] = 9.0 / 120;
    }
    return 5;
}
#endif
