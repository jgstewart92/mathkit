/* mk_fft.h: radix-2 1-D, 2-D and 3-D transforms and mk_fftm's any-length transform against a direct O(n^2) DFT;
   Parseval, round trips, the real-input transform against the complex one, and refused sizes. */
#include "../include/mk_fft.h"
#include "t.h"
#include <stdlib.h>
#include <string.h>

static void dft(const double complex *x, double complex *X, int n, int inv)
{
    for (int k = 0; k < n; k++) {
        double complex s = 0;
        for (int j = 0; j < n; j++) {
            double a = (inv ? 2 : -2) * T_PI * (double)((long)j * k % n) / n;
            s += x[j] * CMPLX(cos(a), sin(a));
        }
        X[k] = s;
    }
}

int main(void)
{
    enum { N = 64 };
    mk_fft p;
    CHECK(mk_fft_plan(&p, 48) == -1, "a size not a power of 2 is refused");
    mk_fft_plan(&p, N);
    double complex a[N];
    for (int j = 0; j < N; j++) a[j] = 3 * cos(2 * T_PI * 5 * j / N);
    mk_fft_c(&p, a, 0);
    double off = 0;
    for (int k = 0; k < N; k++) if (k != 5 && k != N - 5) off = fmax(off, cabs(a[k]));
    CHECK(cabs(a[5] - 96) < 1e-12 && cabs(a[N - 5] - 96) < 1e-12 && off < 1e-12, "3 cos(2 pi 5 j / 64): bins 5 and 59 hold 96, the rest under 1e-12 (%.2g)", off);

    int n = 1024;
    mk_fft q;
    mk_fft_plan(&q, n);
    double complex *x = malloc(sizeof *x * n), *y = malloc(sizeof *y * n);
    double ex = 0, eX = 0;
    for (int j = 0; j < n; j++) x[j] = CMPLX(t_rand(), t_rand()), ex += creal(x[j] * conj(x[j]));
    memcpy(y, x, sizeof *x * n);
    mk_fft_c(&q, y, 0);
    for (int j = 0; j < n; j++) eX += creal(y[j] * conj(y[j]));
    CHECK(fabs(eX / n - ex) < 1e-10 * ex, "Parseval, n 1024: sum |X|^2 / n = sum |x|^2 (rel %.2g)", fabs(eX / n - ex) / ex);
    mk_fft_c(&q, y, 1);
    double rt = 0;
    for (int j = 0; j < n; j++) rt = fmax(rt, cabs(y[j] / n - x[j]));
    CHECK(rt < 1e-13, "forward then inverse / n gives x back (worst %.2g)", rt);

    enum { M = 32 };
    mk_fft pm;
    mk_fft_plan(&pm, M);
    double complex u[M], U[M], V[M];
    double re[M], im[M];
    for (int j = 0; j < M; j++) u[j] = CMPLX(t_rand(), t_rand()), re[j] = creal(u[j]), im[j] = cimag(u[j]);
    dft(u, U, M, 0);
    memcpy(V, u, sizeof u);
    mk_fft_c(&pm, V, 0);
    mk_fft_ri(&pm, re, im, 0);
    double ed = 0, es = 0;
    for (int k = 0; k < M; k++) ed = fmax(ed, cabs(V[k] - U[k])), es = fmax(es, cabs(CMPLX(re[k], im[k]) - V[k]));
    CHECK(ed < 1e-12, "complex form against a direct DFT, n 32 (worst %.2g)", ed);
    CHECK(es == 0, "split form equals the complex form bit for bit (worst %.2g)", es);

    enum { N1 = 8, N2 = 4, N3 = 2 };
    double r2[N1 * N2], i2[N1 * N2];
    double complex g[N1 * N2], G[N1 * N2];
    for (int j = 0; j < N1 * N2; j++) g[j] = CMPLX(t_rand(), t_rand()), r2[j] = creal(g[j]), i2[j] = cimag(g[j]);
    for (int k2 = 0; k2 < N2; k2++)
        for (int k1 = 0; k1 < N1; k1++) {
            double complex s = 0;
            for (int j2 = 0; j2 < N2; j2++)
                for (int j1 = 0; j1 < N1; j1++) {
                    double ang = -2 * T_PI * ((double)j1 * k1 / N1 + (double)j2 * k2 / N2);
                    s += g[j2 * N1 + j1] * CMPLX(cos(ang), sin(ang));
                }
            G[k2 * N1 + k1] = s;
        }
    mk_fft2_ri(r2, i2, N1, N2, 0);
    double e2 = 0;
    for (int j = 0; j < N1 * N2; j++) e2 = fmax(e2, cabs(CMPLX(r2[j], i2[j]) - G[j]));
    CHECK(e2 < 1e-12, "2-D 8 x 4 against a direct DFT (worst %.2g)", e2);

    double complex h[N1 * N2 * N3], H[N1 * N2 * N3];
    for (int j = 0; j < N1 * N2 * N3; j++) h[j] = CMPLX(t_rand(), t_rand());
    for (int k3 = 0; k3 < N3; k3++)
        for (int k2 = 0; k2 < N2; k2++)
            for (int k1 = 0; k1 < N1; k1++) {
                double complex s = 0;
                for (int j3 = 0; j3 < N3; j3++)
                    for (int j2 = 0; j2 < N2; j2++)
                        for (int j1 = 0; j1 < N1; j1++) {
                            double ang = -2 * T_PI * ((double)j1 * k1 / N1 + (double)j2 * k2 / N2 + (double)j3 * k3 / N3);
                            s += h[(j3 * N2 + j2) * N1 + j1] * CMPLX(cos(ang), sin(ang));
                        }
                H[(k3 * N2 + k2) * N1 + k1] = s;
            }
    mk_fft3_c(h, N1, N2, N3, 0);
    double e3 = 0;
    for (int j = 0; j < N1 * N2 * N3; j++) e3 = fmax(e3, cabs(h[j] - H[j]));
    CHECK(e3 < 1e-12, "3-D 8 x 4 x 2 against a direct DFT (worst %.2g)", e3);

    double xr[N];
    double complex full[N], R[N / 2 + 1];
    for (int j = 0; j < N; j++) xr[j] = t_rand(), full[j] = xr[j];
    mk_fft_c(&p, full, 0);
    mk_fft half;
    mk_fft_plan(&half, N / 2);
    mk_rfft(&half, xr, R);
    double er = 0;
    for (int k = 0; k <= N / 2; k++) er = fmax(er, cabs(R[k] - full[k]));
    CHECK(er < 1e-12, "real-input transform, n 64: bins 0..32 match the complex transform (worst %.2g)", er);
    mk_fft_free(&p), mk_fft_free(&q), mk_fft_free(&pm), mk_fft_free(&half), free(x), free(y);
    /* a size that is not a power of 2 is refused before the data is touched */
    double complex g3[48], g30[48];
    for (int i = 0; i < 48; i++) g3[i] = g30[i] = CMPLX(i % 5, 0);
    int r3 = mk_fft3_c(g3, 4, 4, 3, 0), ch = 0;
    for (int i = 0; i < 48; i++) ch += memcmp(&g3[i], &g30[i], sizeof g3[i]) != 0;
    CHECK(r3 == -1 && ch == 0, "mk_fft3_c(4, 4, 3) returns -1 with the data untouched (%d of 48 values changed)", ch);
    { /* n above 2^29 would overflow mk_fft_ri's int doubling: refused before any allocation */
        mk_fft big;
        CHECK(mk_fft_plan(&big, 1 << 30) == -1 && big.n == 0, "mk_fft_plan refuses n = 2^30");
        CHECK(mk_fft_plan(&big, 1 << 20) == 0, "mk_fft_plan still accepts n = 2^20");
        mk_fft_free(&big);
    }
    { /* mk_fftm at any length (factors 4, 2, 3, 5, 7 and the prime 97) against the direct DFT, both directions */
        static const int ns[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 16, 25, 49, 60, 97, 180, 210, 288, 1080, 2160};
        enum { NN = sizeof ns / sizeof ns[0] };
        double worst[2] = {0, 0};
        int wn[2] = {0, 0}, okp = 1;
        for (int t = 0; t < NN; t++) {
            const int m = ns[t];
            mk_fftm pl;
            double complex *xa = malloc(sizeof *xa * (size_t)m), *xb = malloc(sizeof *xb * (size_t)m), *Xr = malloc(sizeof *Xr * (size_t)m),
                           *w = malloc(sizeof *w * (size_t)m);
            if (mk_fftm_plan(&pl, m)) okp = 0;
            for (int inv = 0; inv < 2; inv++) {
                double big = 0, e = 0;
                for (int j = 0; j < m; j++) xa[j] = xb[j] = CMPLX(t_rand(), t_rand());
                dft(xa, Xr, m, inv);
                if (pl.n == m) mk_fftm_c(&pl, xb, w, inv);
                for (int k = 0; k < m; k++) big = fmax(big, cabs(Xr[k])), e = fmax(e, cabs(xb[k] - Xr[k]));
                if (!(e / big <= worst[inv])) worst[inv] = e / big, wn[inv] = m;
            }
            mk_fftm_free(&pl), free(xa), free(xb), free(Xr), free(w);
        }
        CHECK(okp && worst[0] < 1e-13, "mk_fftm forward, 21 lengths 1 to 2160, against a direct DFT (worst %.2g of the largest bin, at n %d)", worst[0], wn[0]);
        CHECK(okp && worst[1] < 1e-13, "mk_fftm inverse, the same lengths, against a direct inverse DFT (worst %.2g, at n %d)", worst[1], wn[1]);
        mk_fftm pl;
        enum { NB = 3840 };
        double complex *xa = malloc(sizeof *xa * NB), *xb = malloc(sizeof *xb * NB), *w = malloc(sizeof *w * NB);
        for (int j = 0; j < NB; j++) xa[j] = xb[j] = CMPLX(t_rand(), t_rand());
        double rt3 = INFINITY;
        if (!mk_fftm_plan(&pl, NB)) {
            mk_fftm_c(&pl, xb, w, 0), mk_fftm_c(&pl, xb, w, 1), rt3 = 0;
            for (int j = 0; j < NB; j++) rt3 = fmax(rt3, cabs(xb[j] / NB - xa[j]));
        }
        CHECK(rt3 < 1e-14, "mk_fftm n 3840: forward then inverse / n gives x back (worst %.2g)", rt3);
        mk_fftm_free(&pl);
        mk_fft p2;
        double e2p = INFINITY;
        if (!mk_fftm_plan(&pl, 1024) && !mk_fft_plan(&p2, 1024)) {
            for (int j = 0; j < 1024; j++) xa[j] = xb[j] = CMPLX(t_rand(), t_rand());
            mk_fftm_c(&pl, xa, w, 0), mk_fft_c(&p2, xb, 0), e2p = 0;
            for (int j = 0; j < 1024; j++) e2p = fmax(e2p, cabs(xa[j] - xb[j]));
            mk_fft_free(&p2);
        }
        CHECK(e2p < 1e-12, "mk_fftm n 1024 equals mk_fft_c to %.2g (another order of the same sums)", e2p);
        mk_fftm_free(&pl), free(xa), free(xb), free(w);
        mk_fftm bad;
        CHECK(mk_fftm_plan(&bad, 0) == -1 && mk_fftm_plan(&bad, -3) == -1 && mk_fftm_plan(&bad, (1 << 29) + 1) == -1 && bad.n == 0,
              "mk_fftm_plan refuses n 0, -3 and 2^29 + 1");
    }
    return wrong;
}
