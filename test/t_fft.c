/* mk_fft.h: radix-2 1-D, 2-D and 3-D transforms against a direct O(n^2) DFT; Parseval, round trips, the real-input
   transform against the complex one, and refused sizes; mk_fftm against a long-double DFT within its stated bound at
   every n to 520 and at chosen lengths to 2^20, its stored twiddles and chirps correctly rounded, and its direct sum
   for a Bluestein length without memory. */
#include "../include/mk_fft.h"
#include "t.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>

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

/* mk_fftm's checks: a long-double reference, its angles reduced in integers so that even tiny sines keep their digits */
static const long double T_PIL = 3.141592653589793238462643383279502884L;
static long double sinpi_ld(long long j, long long n) /* sin(pi j / n) */
{
    long double sg = 1;
    j %= 2 * n;
    if (j < 0) j += 2 * n;
    if (j >= n) j -= n, sg = -1;
    if (2 * j > n) j = n - j;
    return sg * sinl(T_PIL * (long double)j / (long double)n);
}
static long double cospi_ld(long long j, long long n) { return sinpi_ld(n + 2 * j, 2 * n); }
/* X at the bins ks[0 .. nk), forward (inv 0, e^{-i}) or inverse, in long double */
static void dft_ld(const double complex *x, long n, const long *ks, long nk, int inv, long double *Xr, long double *Xi)
{
    long double *c = malloc(sizeof *c * (size_t)n), *s = malloc(sizeof *s * (size_t)n);
    for (long m = 0; m < n; m++) c[m] = cospi_ld(2 * m, n), s[m] = (inv ? 1 : -1) * sinpi_ld(2 * m, n);
    for (long q = 0; q < nk; q++) {
        long double tr = 0, ti = 0;
        for (long j = 0, idx = 0; j < n; j++, idx = idx + ks[q] < n ? idx + ks[q] : idx + ks[q] - n) {
            const long double xr = creal(x[j]), xi = cimag(x[j]);
            tr += xr * c[idx] - xi * s[idx], ti += xr * s[idx] + xi * c[idx];
        }
        Xr[q] = tr, Xi[q] = ti;
    }
    free(c), free(s);
}
/* mk_fftm at n against dft_ld at every bin (n <= 8192) or at 64 (0, 1, n - 1 and 61 drawn), each direction: r[d][0] the
   largest |error| over 2 eps (log2 n + 4) times the largest |X|, r[d][1] the rms error over eps (log2 n + 4) times the
   rms |X|; with direct, through the no-memory fallback. 0, or -1 (no plan). */
static int fftm_bound(int n, int direct, double r[2][2])
{
    mk_fftm p;
    if (mk_fftm_plan(&p, n)) return -1;
    const long nk = n <= 8192 ? n : 64;
    double complex *x = malloc(sizeof *x * (size_t)n), *y = malloc(sizeof *y * (size_t)n), *w = malloc(sizeof *w * (size_t)n);
    long *ks = malloc(sizeof *ks * (size_t)nk);
    long double *Xr = malloc(sizeof *Xr * (size_t)nk), *Xi = malloc(sizeof *Xi * (size_t)nk);
    for (long q = 0; q < nk; q++) ks[q] = n <= 8192 ? q : q == 0 ? 0 : q == 1 ? 1 : q == 2 ? n - 1 : (long)((t_rand() + 1) / 2 * n) % n;
    for (int j = 0; j < n; j++) x[j] = CMPLX(t_rand(), t_rand());
    const double b = DBL_EPSILON * (log2((double)n) + 4);
    for (int d = 0; d < 2; d++) {
        memcpy(y, x, sizeof *x * (size_t)n);
        if (direct) mk_fftm_x_(&p, y, w, NULL, d);
        else mk_fftm_c(&p, y, w, d);
        dft_ld(x, n, ks, nk, d, Xr, Xi);
        long double big = 0, s2 = 0, eb = 0, e2 = 0;
        for (long q = 0; q < nk; q++) {
            const long double er = creal(y[ks[q]]) - Xr[q], ei = cimag(y[ks[q]]) - Xi[q], a2 = Xr[q] * Xr[q] + Xi[q] * Xi[q];
            big = fmaxl(big, sqrtl(a2)), s2 += a2, eb = fmaxl(eb, sqrtl(er * er + ei * ei)), e2 += er * er + ei * ei;
        }
        r[d][0] = big > 0 ? (double)(eb / big) / (2 * b) : (double)eb, r[d][1] = s2 > 0 ? (double)sqrtl(e2 / s2) / b : (double)sqrtl(e2);
    }
    mk_fftm_free(&p), free(x), free(y), free(w), free(ks), free(Xr), free(Xi);
    return 0;
}
/* how far a stored root lies from e^{-i pi j / n} (j over n of a half turn), in ulps of the exact parts */
static double root_ulps(double complex v, long long j, long long n)
{
    const long double ex[2] = {cospi_ld(j, n), -sinpi_ld(j, n)};
    const double got[2] = {creal(v), cimag(v)};
    double worst = 0;
    for (int k = 0; k < 2; k++) {
        int e;
        if (ex[k] == 0) { worst = got[k] == 0 ? worst : INFINITY; continue; }
        frexpl(fabsl(ex[k]), &e);
        worst = fmax(worst, (double)(fabsl(got[k] - ex[k]) / ldexpl(1, e - 53)));
    }
    return worst;
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
    { /* mk_fftm's bound at every length to 520 (radices 2 to 7 written out, the other primes to 61 by their pairs, the
         rest by Bluestein) and at chosen lengths to 2^20, both directions */
        static const int big[] = {1031, 1080, 2062, 2160, 3840, 4099, 4320, 8191, 44100, 65536, 65537, 1048573, 1048576};
        enum { NBIG = sizeof big / sizeof big[0] };
        for (int set = 0; set < 2; set++) {
            double wr[2][2] = {{0, 0}, {0, 0}};
            int wn[2][2] = {{0, 0}, {0, 0}}, okp = 1;
            for (int t = 0; t < (set ? NBIG : 520); t++) {
                double r[2][2];
                const int n = set ? big[t] : t + 1;
                if (fftm_bound(n, 0, r)) { okp = 0; continue; }
                for (int d = 0; d < 2; d++)
                    for (int k = 0; k < 2; k++)
                        if (!(r[d][k] <= wr[d][k])) wr[d][k] = r[d][k], wn[d][k] = n;
            }
            for (int d = 0; d < 2; d++)
                CHECK(okp && wr[d][0] <= 1 && wr[d][1] <= 1,
                      "mk_fftm %s, %s, against a long-double DFT: max |error| <= 2 eps (log2 n + 4) max |X| and rms <= eps (log2 n + 4) rms |X| (worst shares of the bounds %.3f at n %d, %.3f at n %d)",
                      d ? "inverse" : "forward", set ? "13 chosen lengths from 1031 to 2^20 (64 bins sampled above 8192)" : "every n from 1 to 520", wr[d][0], wn[d][0], wr[d][1], wn[d][1]);
        }
    }
    { /* the plans' stored roots, against a reference reduced its own way */
        static const int tn[] = {236, 1031, 1080, 2062, 2187, 2401, 3125, 4096, 4099, 65537};
        double worst = 0;
        long cnt = 0;
        int okp = 1;
        for (int i = 0; i < (int)(sizeof tn / sizeof tn[0]); i++) {
            mk_fftm p;
            if (mk_fftm_plan(&p, tn[i])) { okp = 0; continue; }
            for (long long st = 0, s = 1; st < p.nf; st++) {
                const long long r = p.f[st], m = p.L / s / r;
                const double complex *W = p.tw + p.off[st];
                for (long long j = 0; j < m; j++)
                    for (long long t = 1; t < r; t++) worst = fmax(worst, root_ulps(W[j * (r - 1) + t - 1], 2 * j * t * s, p.L)), cnt++;
                if (r > 7)
                    for (long long k = 0; k < r; k++) worst = fmax(worst, root_ulps(W[m * (r - 1) + k], 2 * k, r)), cnt++;
                s *= r;
            }
            if (p.w)
                for (long long k = 0; k < p.n; k++) worst = fmax(worst, root_ulps(p.w[k], k * k % (2LL * p.n), p.n)), cnt++;
            mk_fftm_free(&p);
        }
        CHECK(okp && worst <= 0.501, "mk_fftm's twiddles and chirps rounded right: within %.4f ulp of each exact part (%ld values in 10 plans; 0.5 is correct rounding)", worst, cnt);
    }
    { /* a Bluestein length without its scratch: the direct sum */
        static const int dn[] = {67, 97, 1031};
        double wr = 0;
        int wn = 0, okp = 1;
        for (int i = 0; i < 3; i++) {
            double r[2][2];
            if (fftm_bound(dn[i], 1, r)) { okp = 0; continue; }
            for (int d = 0; d < 2; d++)
                for (int k = 0; k < 2; k++)
                    if (!(r[d][k] <= wr)) wr = r[d][k], wn = dn[i];
        }
        CHECK(okp && wr <= 1, "mk_fftm without its scratch (the direct sum), n 67, 97 and 1031 both ways: within the same bounds (worst share %.3f at n %d)", wr, wn);
    }
    {
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
