/* mathkit: spectral analysis on mk_fft.h. Windows (rectangular, Hann, Hamming, 4-term Blackman-Harris, Kaiser),
   periodic (DFT-even) as spectral analysis wants; STFT frames; periodograms and Welch averages, one-sided, as a density
   (units^2 / Hz: white noise of variance s^2 reads 2 s^2 / fs) or as a power spectrum (a sine of amplitude A centred on
   a bin reads A^2/2 there); cross-spectra (conj(X) Y, scipy's convention), transfer function and coherence. nfft a
   power of 2, >= 4. Header only. */
#ifndef MK_SPEC_H
#define MK_SPEC_H
#include "mk_fft.h"
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum { MK_WIN_RECT = 0, MK_WIN_HANN = 1, MK_WIN_HAMMING = 2, MK_WIN_BLACKMAN_HARRIS = 3, MK_WIN_KAISER = 4 };
enum { MK_SPEC_DENSITY = 0, MK_SPEC_POWER = 1 };

/* modified Bessel function of the first kind, order 0, by its power series (all terms positive: no cancellation) */
static inline double mk_bessel_i0(double x)
{
    double q = 0.25 * x * x, term = 1, sum = 1;
    for (int k = 1; k < 500; k++) {
        term *= q / ((double)k * k);
        sum += term;
#if defined(MK_MUTATE) && MK_MUTATE == 82
        if (k >= 4) break; /* planted fault: the series cut short */
#endif
        if (term < 1e-17 * sum) break;
    }
    return sum;
}
/* w[n], periodic; beta is Kaiser's shape (8.6 is near Blackman's) and unused otherwise */
static inline void mk_window(int type, int n, double beta, double *w)
{
    const double tp = 2 * MK_FFT_PI;
    for (int i = 0; i < n; i++) {
        double c = cos(tp * i / n);
        switch (type) {
        case MK_WIN_HANN: w[i] = 0.5 - 0.5 * c; break;
        case MK_WIN_HAMMING: w[i] = 0.54 - 0.46 * c; break;
        case MK_WIN_BLACKMAN_HARRIS:
#if defined(MK_MUTATE) && MK_MUTATE == 83 /* planted fault: a wrong coefficient */
            w[i] = 0.35875 - 0.48 * c + 0.14128 * cos(2 * tp * i / n) - 0.01168 * cos(3 * tp * i / n);
#else
            w[i] = 0.35875 - 0.48829 * c + 0.14128 * cos(2 * tp * i / n) - 0.01168 * cos(3 * tp * i / n);
#endif
            break;
        case MK_WIN_KAISER: {
            double r = 2.0 * i / n - 1;
            w[i] = mk_bessel_i0(beta * sqrt(fmax(0, 1 - r * r))) / mk_bessel_i0(beta);
            break;
        }
        default: w[i] = 1;
        }
    }
}

typedef struct {
    int nfft;
    double fs, s1, s2; /* sum of w, sum of w^2 */
    double *w;
    mk_fft half;
} mk_spec;

static inline void mk_spec_free(mk_spec *s)
{
    free(s->w), mk_fft_free(&s->half);
    s->w = NULL;
}
/* 0, or -1 for a size not a power of 2 (or under 4), a window not one of MK_WIN_*, fs not positive and finite, or no
   memory. mk_spec_free releases it. */
static inline int mk_spec_init(mk_spec *s, int nfft, int window, double beta, double fs)
{
    s->w = NULL, s->s1 = s->s2 = 0, s->nfft = 0, s->fs = 1;
    s->half.rev = NULL, s->half.c = s->half.s = NULL, s->half.n = 0;
    if (nfft < 4 || (nfft & (nfft - 1)) || window < MK_WIN_RECT || window > MK_WIN_KAISER || !(fs > 0 && fs < INFINITY) || mk_fft_plan(&s->half, nfft / 2)) return -1;
    s->w = malloc(sizeof(double) * (size_t)nfft);
    if (!s->w) { mk_fft_free(&s->half); return -1; }
    s->nfft = nfft, s->fs = fs;
    mk_window(window, nfft, beta, s->w);
    s->s1 = s->s2 = 0;
    for (int i = 0; i < nfft; i++) s->s1 += s->w[i], s->s2 += s->w[i] * s->w[i];
    return 0;
}
/* the windowed frame x[start .. start + nfft) to X[0 .. nfft/2]; samples outside [0, len) count as 0; detrend 1 removes
   the frame's mean first (over the samples inside). The plan is only read, so one plan serves any number of threads. */
static inline void mk_stft_frame(const mk_spec *s, const double *x, long len, long start, int detrend, double complex *X)
{
    int n = s->nfft;
    double mean = 0;
    long cnt = 0;
    if (detrend)
        for (int i = 0; i < n; i++) if (start + i >= 0 && start + i < len) mean += x[start + i], cnt++;
    mean = cnt ? mean / cnt : 0;
    for (int i = 0; i < n; i += 2) {
        const long j = start + i;
        X[i / 2] = CMPLX(j >= 0 && j < len ? (x[j] - mean) * s->w[i] : 0, j + 1 >= 0 && j + 1 < len ? (x[j + 1] - mean) * s->w[i + 1] : 0);
    }
    mk_rfft_packed_(&s->half, X);
}
/* |X|^2 (or a cross product) to the one-sided scale: density or power */
static inline double mk_spec_scale_(const mk_spec *s, int k, int scaling)
{
    double one = (k == 0 || k == s->nfft / 2) ? 1 : 2;
#if defined(MK_MUTATE) && MK_MUTATE == 81
    one = 1; /* planted fault: the one-sided spectrum not doubled */
#endif
    return scaling == MK_SPEC_POWER ? one / (s->s1 * s->s1) : one / (s->fs * s->s2);
}
/* Welch: segments of nfft stepping nfft - noverlap (nfft if that is not positive) through x[len], the last partial one
   dropped, each detrended to its mean if detrend, periodograms averaged into p[0 .. nfft/2]. Returns the segments used
   (0: x shorter than one), or -1 for no memory. */
static inline long mk_welch(const mk_spec *s, const double *x, long len, int noverlap, int detrend, int scaling, double *p)
{
    int n = s->nfft, nb = n / 2 + 1, step = n - noverlap > 0 ? n - noverlap : n;
    double complex *X = malloc(sizeof *X * (size_t)nb);
    if (!X) return -1;
    for (int k = 0; k < nb; k++) p[k] = 0;
    long seg = 0;
    for (long st = 0; st + n <= len; st += step, seg++) {
        mk_stft_frame(s, x, len, st, detrend, X);
        for (int k = 0; k < nb; k++) p[k] += creal(X[k]) * creal(X[k]) + cimag(X[k]) * cimag(X[k]);
    }
    for (int k = 0; k < nb; k++) p[k] = seg ? p[k] / seg * mk_spec_scale_(s, k, scaling) : 0;
    free(X);
    return seg;
}
/* Welch cross-spectral density of x and y: pxy = <conj(X) Y>, and pxx, pyy (any may be NULL), density scaling; segments
   and return as mk_welch */
static inline long mk_csd(const mk_spec *s, const double *x, const double *y, long len, int noverlap, int detrend,
                          double complex *pxy, double *pxx, double *pyy)
{
    int n = s->nfft, nb = n / 2 + 1, step = n - noverlap > 0 ? n - noverlap : n;
    double complex *X = malloc(sizeof *X * 2 * (size_t)nb), *Y = X ? X + nb : NULL;
    if (!X) return -1;
    for (int k = 0; k < nb; k++) {
        if (pxy) pxy[k] = 0;
        if (pxx) pxx[k] = 0;
        if (pyy) pyy[k] = 0;
    }
    long seg = 0;
    for (long st = 0; st + n <= len; st += step, seg++) {
        mk_stft_frame(s, x, len, st, detrend, X), mk_stft_frame(s, y, len, st, detrend, Y);
        for (int k = 0; k < nb; k++) {
#if defined(MK_MUTATE) && MK_MUTATE == 84
            if (pxy) pxy[k] += X[k] * Y[k]; /* planted fault: no conjugate */
#else
            if (pxy) pxy[k] += conj(X[k]) * Y[k];
#endif
            if (pxx) pxx[k] += creal(X[k] * conj(X[k]));
            if (pyy) pyy[k] += creal(Y[k] * conj(Y[k]));
        }
    }
    for (int k = 0; k < nb; k++) {
        double c = seg ? mk_spec_scale_(s, k, MK_SPEC_DENSITY) / seg : 0;
        if (pxy) pxy[k] *= c;
        if (pxx) pxx[k] *= c;
        if (pyy) pyy[k] *= c;
    }
    free(X);
    return seg;
}
/* magnitude-squared coherence |Pxy|^2 / (Pxx Pyy), 0 where a spectrum is 0; and the transfer function H = Pxy / Pxx */
static inline void mk_coherence(const double complex *pxy, const double *pxx, const double *pyy, int nb, double *coh)
{
    for (int k = 0; k < nb; k++) {
        double d = pxx[k] * pyy[k];
        coh[k] = d > 0 ? creal(pxy[k] * conj(pxy[k])) / d : 0;
    }
}
static inline void mk_transfer(const double complex *pxy, const double *pxx, int nb, double complex *H)
{
    for (int k = 0; k < nb; k++) H[k] = pxx[k] > 0 ? pxy[k] / pxx[k] : 0;
}
static inline double mk_spec_freq(const mk_spec *s, int k) { return k * s->fs / s->nfft; }
#endif
