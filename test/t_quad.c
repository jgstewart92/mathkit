#include <stdlib.h>
static long t_mallocs; /* allocations the header asked for */
static void *t_malloc(size_t s) { t_mallocs++; return malloc(s); }
#define malloc t_malloc
#include "../include/mk_quad.h"
#undef malloc
#include "t.h"
#include <limits.h>

static double fact(int k) { double f = 1; for (int i = 2; i <= k; i++) f *= i; return f; }
static int pw;
static double mono(double x, void *c) { (void)c; return pow(x, pw); }
static double isq(double x, void *c) { (void)c; return 1 / sqrt(x); }
static double lg(double x, void *c) { (void)c; return log(x); }
static double osc(double x, void *c) { (void)c; return cos(30 * x); }

static long spike_calls;
static double spike(double x, void *c)
{
    const double *f = c;
    spike_calls++;
    return f[2] * exp(-f[0] * (x - f[1]) * (x - f[1])) + 1;
}

int main(void)
{
    double x[40], w[40];
    int bad = 0, sharp = 1;
    for (int n = 1; n <= 20; n++) {
        if (mk_gauss_legendre(n, x, w)) { bad++; continue; }
        for (int k = 0; k <= 2 * n; k++) {
            double s = 0, exact = k % 2 ? 0 : 2.0 / (k + 1);
            for (int i = 0; i < n; i++) s += w[i] * pow(x[i], k);
            if (k <= 2 * n - 1 && fabs(s - exact) > 1e-14 * 2) bad++;
            if (k == 2 * n && fabs(s - exact) < 1e-13 * exact) sharp = 0; /* the error there is about pi 4^-n */
        }
    }
    CHECK(bad == 0, "Gauss-Legendre n 1..20: exact for x^k, k <= 2n-1 (%d misses)", bad);
    CHECK(sharp, "and not exact at k = 2n (the degree is not overstated)");
    mk_gauss_legendre(5, x, w);
    double x5[5] = {-0.9061798459386640, -0.5384693101056831, 0, 0.5384693101056831, 0.9061798459386640};
    double w5[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891}, e5 = 0;
    for (int i = 0; i < 5; i++) e5 = fmax(e5, fmax(fabs(x[i] - x5[i]), fabs(w[i] - w5[i])));
    CHECK(e5 < 1e-15, "n 5 against the tabulated nodes and weights (worst %.2g)", e5);

    int kb = 0, eb = 0, en = 0;
    for (pw = 0; pw <= 23; pw++) {
        double err, v = mk_gk15(mono, NULL, 0, 1, &err), exact = 1.0 / (pw + 1);
        if (pw <= 22 && fabs(v - exact) > 2e-15) kb++;
        if (pw <= 13 && err > 1e-15) eb++;
        if (pw >= 16 && pw <= 22 && err < 1e-8) en++;
    }
    CHECK(kb == 0, "Kronrod 15: exact for x^k on [0, 1] to k 22 (%d misses)", kb);
    CHECK(eb == 0 && en == 0, "its error estimate is 0 to k 13 (Gauss 7 is exact) and not 0 for k 16..22 (%d, %d wrong)", eb, en);
    double err;
    int st;
    double v1 = mk_quad_adapt(isq, NULL, 0, 1, 1e-12, 0, 2000, &err, &st);
    CHECK(st == 0 && fabs(v1 - 2) < 1e-10 && err >= fabs(v1 - 2) * 0.1, "adaptive, 1/sqrt(x) on [0, 1] (singular end): %.14f, error %.2g, estimate %.2g", v1, fabs(v1 - 2), err);
    double v2 = mk_quad_adapt(lg, NULL, 0, 1, 1e-12, 0, 2000, &err, &st);
    CHECK(st == 0 && fabs(v2 + 1) < 1e-11, "adaptive, log x on [0, 1]: %.14f (-1)", v2);
    double v3 = mk_quad_adapt(osc, NULL, 0, 10, 1e-13, 0, 2000, &err, &st), want = sin(300.0) / 30;
    CHECK(st == 0 && fabs(v3 - want) < 1e-12, "adaptive, cos 30x on [0, 10]: error %.2g", fabs(v3 - want));

    double tp[7][2], tw[7], qp[5][3], qw[5];
    int tdeg[3] = {1, 2, 5}, tbad = 0, tsharp = 1;
    for (int r = 0; r < 3; r++) {
        int m = mk_tri_rule(tdeg[r], tp, tw);
        for (int a = 0; a <= tdeg[r] + 1; a++)
            for (int b = 0; a + b <= tdeg[r] + 1; b++) {
                double s = 0, ex = fact(a) * fact(b) / fact(a + b + 2);
                for (int i = 0; i < m; i++) s += tw[i] * pow(tp[i][0], a) * pow(tp[i][1], b);
                if (a + b <= tdeg[r] && fabs(s - ex) > 1e-15) tbad++;
                if (a + b == tdeg[r] + 1 && a == tdeg[r] + 1 && fabs(s - ex) < 1e-12 * ex) tsharp = 0;
            }
    }
    CHECK(tbad == 0 && tsharp, "triangle rules of degree 1, 2, 5: exact for x^a y^b to their degree, not past it (%d misses)", tbad);
    int qdeg[3] = {1, 2, 3}, qbad = 0, qsharp = 1;
    for (int r = 0; r < 3; r++) {
        int m = mk_tet_rule(qdeg[r], qp, qw);
        for (int a = 0; a <= qdeg[r] + 1; a++)
            for (int b = 0; a + b <= qdeg[r] + 1; b++)
                for (int c = 0; a + b + c <= qdeg[r] + 1; c++) {
                    double s = 0, ex = fact(a) * fact(b) * fact(c) / fact(a + b + c + 3);
                    for (int i = 0; i < m; i++) s += qw[i] * pow(qp[i][0], a) * pow(qp[i][1], b) * pow(qp[i][2], c);
                    if (a + b + c <= qdeg[r] && fabs(s - ex) > 1e-15) qbad++;
                    if (a == qdeg[r] + 1 && fabs(s - ex) < 1e-12 * ex) qsharp = 0;
                }
    }
    CHECK(qbad == 0 && qsharp, "tetrahedron rules of degree 1, 2, 3: exact for x^a y^b z^c to their degree, not past it (%d misses)", qbad);
    { /* 1.0.5: the verification work's quad finding Q-A (Result 79) and the size guards its proofs asked for */
        static const double qa[6][4] = { /* p, q, s, reltol: s exp(-p (x - q)^2) + 1 on [-1, 1], from a search on 1.0.4 */
            {0x1.1308842f326f9p+15, -0x1.839c35f9524cap-1, 0x1.04725e871e93p+27, 0x1.6849b86a12b9bp-47},
            {0x1.88b8f7e7806fdp+18, 0x1.9f8456d11faa8p-2, 0x1.9886d8d403ba3p+22, 0x1.5fd7fe1796495p-37},
            {0x1.74b00cc38f519p+14, -0x1.7c15b4abda456p-1, 0x1.3614f2e9ace36p+23, 0x1.19799812dea11p-40},
            {0x1.7d1b2efadb14dp+15, -0x1.25c86c23ad20cp-1, 0x1.ad080ee9d070cp+31, 0x1.c25c268497682p-44},
            {0x1.65486cf69359p+15, -0x1.77415fa0a0b4ep-1, 0x1.cb65123e11e96p+30, 0x1.c25c268497682p-44},
            {0x1.ae21ad510f538p+16, 0x1.2b20973b1f53p-1, 0x1.af2c88c40426ap+37, 0x1.6849b86a12b9bp-47}};
        int qbadst = 0;
        double worst = 0;
        for (int r = 0; r < 6; r++) {
            double e;
            int st;
            double I = mk_quad_adapt(spike, (void *)qa[r], -1, 1, 0, qa[r][3], 2000, &e, &st);
            if (st == 0 && e > qa[r][3] * fabs(I) * (1 + 1e-12)) qbadst++, worst = fmax(worst, e / (qa[r][3] * fabs(I)));
        }
        CHECK(qbadst == 0, "Q-A: status 0 means the error meets the limit: %d of 6 spikes say 0 above it (worst %.3g x the limit)", qbadst, worst);
        long m0 = t_mallocs;
        double e2 = 0;
        int st2 = 9;
        spike_calls = 0;
        static const double flat[3] = {1, 0, 0}; /* f = 1: 1.0.4 allocates and finishes in one pass */
        double I2 = mk_quad_adapt(spike, (void *)flat, -1, 1, 0, 1e-10, 1 << 29, &e2, &st2);
        CHECK(st2 == -1 && isnan(I2) && isnan(e2) && spike_calls == 0 && t_mallocs == m0,
              "mk_quad_adapt refuses maxint 2^29 (4 i + 3 past INT_MAX) before allocating or calling f: status %d, %ld calls, %ld allocations",
              st2, spike_calls, t_mallocs - m0);
#ifndef T_OLD
        double gx[1], gw[1];
        CHECK(mk_gauss_legendre((1 << 30) + 1, gx, gw) == -1 && mk_gauss_legendre(INT_MAX, gx, gw) == -1,
              "mk_gauss_legendre refuses n = 2^30 + 1 and INT_MAX (2 k - 1 and (n + 1) / 2 past an int) before writing x or w");
#endif
    }
    return wrong;
}
