/* mathkit: radix-2 FFT, planned, in place. Forward is e^{-i}, inverse e^{+i}; neither is normalized (divide by n after an
   inverse). Complex (double complex) and split (re[], im[]) forms share one plan. 2-D and 3-D run along each axis in turn.
   A real-input transform of length n gives n/2+1 bins through one complex transform of n/2. Any other length (1.5.6):
   mk_fftm, mixed radix, at the end of this file. Header only, C23.
   Experimental (1.6.0): mk_fftm alone, not frozen with the rest of this file (the radix-2 transforms, mk_rfft). */
#ifndef MK_FFT_H
#define MK_FFT_H
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MK_FFT_PI 3.14159265358979323846

typedef struct {
    int n;
    int *rev;          /* bit-reversed index */
    double *c, *s;     /* cos and sin of 2 pi k / n, k < n/2 */
} mk_fft;

static inline void mk_fft_free(mk_fft *p)
{
    free(p->rev), free(p->c);
    p->rev = NULL, p->c = p->s = NULL, p->n = 0;
}
/* n a power of 2, n >= 1; 0 or -1 */
static inline int mk_fft_plan(mk_fft *p, int n)
{
    p->rev = NULL, p->c = p->s = NULL, p->n = 0;
    if (n < 1 || n > (1 << 29) || (n & (n - 1))) return -1; /* mk_fft_ri doubles an int up to 2n: n <= 2^29 keeps it in range */
    int lg = 0;
    while ((1 << lg) < n) lg++;
    int h = n > 1 ? n / 2 : 1;
    p->rev = malloc(sizeof *p->rev * (size_t)n);
    p->c = malloc(sizeof *p->c * (size_t)(2 * h));
    if (!p->rev || !p->c) { mk_fft_free(p); return -1; }
    p->s = p->c + h;
    p->n = n;
    for (int i = 0; i < n; i++) {
        int r = 0;
        for (int b = 0; b < lg; b++)
            if (i >> b & 1) r |= 1 << (lg - 1 - b);
        p->rev[i] = r;
    }
    for (int k = 0; k < n / 2; k++) p->c[k] = cos(2 * MK_FFT_PI * k / n), p->s[k] = sin(2 * MK_FFT_PI * k / n);
    return 0;
}

static inline void mk_fft_ri(const mk_fft *p, double *re, double *im, int inverse)
{
    int n = p->n;
    for (int i = 0; i < n; i++) {
        int j = p->rev[i];
#if defined(MK_MUTATE) && MK_MUTATE == 12
        j = i;
#endif
        if (i < j) {
            double t = re[i];
            re[i] = re[j], re[j] = t, t = im[i], im[i] = im[j], im[j] = t;
        }
    }
    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1, step = n / len;
        for (int i = 0; i < n; i += len)
            for (int k = 0; k < half; k++) {
#if defined(MK_MUTATE) && MK_MUTATE == 11
                double wr = p->c[k * step], wi = -p->s[k * step];
                (void)inverse;
#else
                double wr = p->c[k * step], wi = inverse ? p->s[k * step] : -p->s[k * step];
#endif
                int a = i + k, b = a + half;
                double xr = re[b] * wr - im[b] * wi, xi = re[b] * wi + im[b] * wr;
                re[b] = re[a] - xr, im[b] = im[a] - xi, re[a] += xr, im[a] += xi;
            }
    }
}

static inline void mk_fft_c(const mk_fft *p, double complex *a, int inverse)
{
    /* the same arithmetic as the split form, on interleaved storage */
    int n = p->n;
    for (int i = 0; i < n; i++) {
        int j = p->rev[i];
#if defined(MK_MUTATE) && MK_MUTATE == 12
        j = i;
#endif
        if (i < j) { double complex t = a[i]; a[i] = a[j], a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1, step = n / len;
        for (int i = 0; i < n; i += len)
            for (int k = 0; k < half; k++) {
#if defined(MK_MUTATE) && MK_MUTATE == 11
                double wr = p->c[k * step], wi = -p->s[k * step];
                (void)inverse;
#else
                double wr = p->c[k * step], wi = inverse ? p->s[k * step] : -p->s[k * step];
#endif
                int ia = i + k, ib = ia + half;
                double br = creal(a[ib]), bi = cimag(a[ib]);
                double xr = br * wr - bi * wi, xi = br * wi + bi * wr;
                double ar = creal(a[ia]), ai = cimag(a[ia]);
                a[ib] = CMPLX(ar - xr, ai - xi);
                a[ia] = CMPLX(ar + xr, ai + xi);
            }
    }
}

/* 2-D on split storage, x fastest: re[j * n1 + i]. Columns go eight at a time through a gathered buffer (metal's wind.c).
   0 or -1 (a size not a power of 2, or no memory). */
static inline int mk_fft2_ri(double *re, double *im, int n1, int n2, int inverse)
{
    enum { B = 8 };
    mk_fft p1, p2;
    if (mk_fft_plan(&p1, n1)) return -1;
    if (mk_fft_plan(&p2, n2)) { mk_fft_free(&p1); return -1; }
    double *cr = malloc(sizeof(double) * (size_t)n2 * B * 2), *ci = cr ? cr + (size_t)n2 * B : NULL;
    if (!cr) { mk_fft_free(&p1), mk_fft_free(&p2); return -1; }
    for (int j = 0; j < n2; j++) mk_fft_ri(&p1, re + (size_t)j * n1, im + (size_t)j * n1, inverse);
    for (int i0 = 0; i0 < n1; i0 += B) {
        int nb = n1 - i0 < B ? n1 - i0 : B;
        for (int j = 0; j < n2; j++)
            for (int b = 0; b < nb; b++)
                cr[(size_t)b * n2 + j] = re[(size_t)j * n1 + i0 + b], ci[(size_t)b * n2 + j] = im[(size_t)j * n1 + i0 + b];
        for (int b = 0; b < nb; b++) mk_fft_ri(&p2, cr + (size_t)b * n2, ci + (size_t)b * n2, inverse);
        for (int j = 0; j < n2; j++)
            for (int b = 0; b < nb; b++)
                re[(size_t)j * n1 + i0 + b] = cr[(size_t)b * n2 + j], im[(size_t)j * n1 + i0 + b] = ci[(size_t)b * n2 + j];
    }
    free(cr), mk_fft_free(&p1), mk_fft_free(&p2);
    return 0;
}

/* 3-D on complex storage, x fastest: a[(k * n2 + j) * n1 + i]. 0 or -1. */
static inline int mk_fft3_c(double complex *a, int n1, int n2, int n3, int inverse)
{
    int nn[3] = {n1, n2, n3}, mx = n1 > n2 ? n1 : n2;
    mx = mx > n3 ? mx : n3;
    for (int ax = 0; ax < 3; ax++) if (nn[ax] < 1 || (nn[ax] & (nn[ax] - 1))) return -1; /* refuse before touching the data */
    size_t stride[3] = {1, (size_t)n1, (size_t)n1 * n2};
    double complex *line = malloc(sizeof *line * (size_t)mx);
    if (!line) return -1;
    for (int ax = 0; ax < 3; ax++) {
        mk_fft p;
        if (mk_fft_plan(&p, nn[ax])) { free(line); return -1; }
        int o1 = (ax + 1) % 3, o2 = (ax + 2) % 3;
        for (int u = 0; u < nn[o1]; u++)
            for (int v = 0; v < nn[o2]; v++) {
                size_t base = (size_t)u * stride[o1] + (size_t)v * stride[o2];
                for (int t = 0; t < nn[ax]; t++) line[t] = a[base + (size_t)t * stride[ax]];
                mk_fft_c(&p, line, inverse);
                for (int t = 0; t < nn[ax]; t++) a[base + (size_t)t * stride[ax]] = line[t];
            }
        mk_fft_free(&p);
    }
    free(line);
    return 0;
}

/* Real input x[n] (n a power of 2, n >= 2) to bins X[0..n/2], forward, unnormalized; the rest are conj(X[n-k]).
   half is a plan for n/2. */
static inline void mk_rfft_packed_(const mk_fft *half, double complex *X);
static inline void mk_rfft(const mk_fft *half, const double *x, double complex *X)
{
    for (int k = 0; k < half->n; k++) X[k] = CMPLX(x[2 * k], x[2 * k + 1]);
    mk_rfft_packed_(half, X);
}
/* mk_rfft on X[0..n/2) already holding the pairs (x[2k], x[2k+1]) (mk_stft_frame packs its frame straight into X) */
static inline void mk_rfft_packed_(const mk_fft *half, double complex *X)
{
    int m = half->n, n = 2 * m;
    mk_fft_c(half, X, 0);
    double complex z0 = X[0];
    X[m] = CMPLX(creal(z0) - cimag(z0), 0);
    X[0] = CMPLX(creal(z0) + cimag(z0), 0);
    for (int k = 1; k <= m / 2; k++) {
        double complex a = X[k], b = conj(X[m - k]);
        double complex e = (a + b) / 2, o = (a - b) / 2;
        double ang = -2 * MK_FFT_PI * k / n;
        double complex w = CMPLX(cos(ang), sin(ang));
        double complex xk = e - I * w * o;
        double complex wm = -conj(w); /* e^(-2 pi i (m - k) / n) = -conj(e^(-2 pi i k / n)) */
        double complex e2 = conj(e), o2 = -conj(o);
        X[m - k] = e2 - I * wm * o2;
        X[k] = xk;
    }
}

/* Any length (1.5.6): a mixed-radix complex FFT, self-sorting (Stockham, decimation in frequency), for lengths such as 1080,
   1920 and 3840 that the radix-2 plan refuses. The same signs as mk_fft (forward e^{-i}, neither direction normalized).
   Radices 4, 2 and 3 by hand, any other prime factor by its direct sum, so a length with a large prime factor is slow but
   right. mk_fftm_plan: n from 1 to 2^29; 0, or -1 (a bad n, no memory). mk_fftm_c: a[n] in place, work[n] scratch (a
   caller's, so one plan serves many threads). */
typedef struct {
    int n, nf, f[32]; /* the radices, in the order of the stages */
    double *c, *s;    /* cos and sin of 2 pi k / n, k < n */
} mk_fftm;

static inline void mk_fftm_free(mk_fftm *p)
{
    free(p->c);
    p->c = p->s = NULL, p->n = p->nf = 0;
}
static inline int mk_fftm_plan(mk_fftm *p, int n)
{
    p->c = p->s = NULL, p->n = p->nf = 0;
    if (n < 1 || n > (1 << 29)) return -1;
    int m = n, nf = 0;
    while (m % 4 == 0) p->f[nf++] = 4, m /= 4;
    while (m % 2 == 0) p->f[nf++] = 2, m /= 2;
    for (int q = 3; m > 1; q += 2) {
        if ((long)q * q > m) q = m; /* what is left is prime */
        while (m % q == 0) p->f[nf++] = q, m /= q;
    }
    p->c = malloc(sizeof *p->c * (size_t)n * 2);
    if (!p->c) return -1;
    p->s = p->c + n;
    for (int k = 0; k < n; k++) p->c[k] = cos(2 * MK_FFT_PI * k / n), p->s[k] = sin(2 * MK_FFT_PI * k / n);
    p->n = n, p->nf = nf;
    return 0;
}
static inline void mk_fftm_c(const mk_fftm *p, double complex *a, double complex *work, int inverse)
{
    const int n = p->n;
    const double sg = inverse ? 1 : -1, h3 = 0.86602540378443864676; /* the exponent's sign; sin(2 pi / 3) */
    double complex *x = a, *y = work;
    for (int st = 0, l = n, s = 1; st < p->nf; st++) {
        const int r = p->f[st], m = l / r; /* l: the length of the transforms this stage splits; s: their stride */
        for (int j = 0; j < m; j++) {
            double wr[4], wi[4]; /* the twiddles e^{sg 2 pi i j t / l}, t < 4; radix 5 and up read the table as they go */
            for (int t = 0; t < (r < 4 ? r : 4); t++) {
                const long k = (long)j * t * s;
                wr[t] = p->c[k], wi[t] = sg * p->s[k];
            }
            const double complex *in = x + (size_t)s * (size_t)j; /* the stage reads in[q + s m u], u < r ... */
            double complex *out = y + (size_t)s * (size_t)r * (size_t)j; /* ... and writes out[q + s t], t < r */
            const size_t sm = (size_t)s * (size_t)m;
            if (r == 4) {
                for (int q = 0; q < s; q++) {
                    const double complex a0 = in[q], a1 = in[q + sm], a2 = in[q + 2 * sm], a3 = in[q + 3 * sm];
                    const double pr = creal(a0) + creal(a2), pi = cimag(a0) + cimag(a2), mr = creal(a0) - creal(a2), mi = cimag(a0) - cimag(a2);
                    const double qr = creal(a1) + creal(a3), qi = cimag(a1) + cimag(a3), dr = creal(a1) - creal(a3), di = cimag(a1) - cimag(a3);
#if defined(MK_MUTATE) && MK_MUTATE == 238
                    const double e = -1; /* 238: the quarter turn always forward's */
#else
                    const double e = sg;
#endif
                    const double b1r = mr - e * di, b1i = mi + e * dr, b2r = pr - qr, b2i = pi - qi, b3r = mr + e * di, b3i = mi - e * dr;
                    out[q] = CMPLX(pr + qr, pi + qi);
                    out[q + s] = CMPLX(b1r * wr[1] - b1i * wi[1], b1r * wi[1] + b1i * wr[1]);
                    out[q + 2 * (size_t)s] = CMPLX(b2r * wr[2] - b2i * wi[2], b2r * wi[2] + b2i * wr[2]);
                    out[q + 3 * (size_t)s] = CMPLX(b3r * wr[3] - b3i * wi[3], b3r * wi[3] + b3i * wr[3]);
                }
            } else if (r == 2) {
                for (int q = 0; q < s; q++) {
                    const double complex a0 = in[q], a1 = in[q + sm];
                    const double br = creal(a0) - creal(a1), bi = cimag(a0) - cimag(a1);
                    out[q] = a0 + a1;
                    out[q + s] = CMPLX(br * wr[1] - bi * wi[1], br * wi[1] + bi * wr[1]);
                }
            } else if (r == 3) {
                for (int q = 0; q < s; q++) {
                    const double complex a0 = in[q], a1 = in[q + sm], a2 = in[q + 2 * sm];
                    const double tr = creal(a1) + creal(a2), ti = cimag(a1) + cimag(a2);
#if defined(MK_MUTATE) && MK_MUTATE == 239
                    const double dr = sg * h3 * (creal(a2) - creal(a1)), di = sg * h3 * (cimag(a2) - cimag(a1)); /* 239: the difference turned */
#else
                    const double dr = sg * h3 * (creal(a1) - creal(a2)), di = sg * h3 * (cimag(a1) - cimag(a2));
#endif
                    const double cr = creal(a0) - tr / 2, ci = cimag(a0) - ti / 2;
                    const double b1r = cr - di, b1i = ci + dr, b2r = cr + di, b2i = ci - dr;
                    out[q] = CMPLX(creal(a0) + tr, cimag(a0) + ti);
                    out[q + s] = CMPLX(b1r * wr[1] - b1i * wi[1], b1r * wi[1] + b1i * wr[1]);
                    out[q + 2 * (size_t)s] = CMPLX(b2r * wr[2] - b2i * wi[2], b2r * wi[2] + b2i * wr[2]);
                }
            } else {
                const long nr = n / r;
                for (int q = 0; q < s; q++)
                    for (int t = 0; t < r; t++) {
                        double br = 0, bi = 0;
                        for (int u = 0; u < r; u++) {
                            const long k = (long)u * t % r * nr;
                            const double cr = p->c[k], ci = sg * p->s[k], ar = creal(in[q + (size_t)u * sm]), ai = cimag(in[q + (size_t)u * sm]);
                            br += ar * cr - ai * ci, bi += ar * ci + ai * cr;
                        }
                        const long k = (long)j * t * s;
                        const double tr = p->c[k], ti = sg * p->s[k];
                        out[q + (size_t)t * s] = CMPLX(br * tr - bi * ti, br * ti + bi * tr);
                    }
            }
        }
        double complex *z = x;
        x = y, y = z, l = m, s *= r;
    }
    if (x != a) memcpy(a, x, sizeof *a * (size_t)n);
}
#endif
