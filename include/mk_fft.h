/* mathkit: radix-2 FFT, planned, in place. Forward is e^{-i}, inverse e^{+i}; neither is normalized (divide by n after
   an inverse). Complex (double complex) and split (re[], im[]) forms share one plan. 2-D and 3-D run along each axis in
   turn. A real-input transform of length n gives n/2+1 bins through one complex transform of n/2. Any other length:
   mk_fftm, mixed radix, at the end of this file. Header only, C23. */
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
    /* mk_fft_ri doubles an int up to 2n: n <= 2^29 keeps it in range */
    if (n < 1 || n > (1 << 29) || (n & (n - 1))) return -1;
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
        j = i; /* planted fault: no bit reversal */
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
                double wr = p->c[k * step], wi = -p->s[k * step]; /* planted fault: the inverse runs forward */
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
        j = i; /* planted fault: no bit reversal */
#endif
        if (i < j) { double complex t = a[i]; a[i] = a[j], a[j] = t; }
    }
    for (int len = 2; len <= n; len <<= 1) {
        int half = len >> 1, step = n / len;
        for (int i = 0; i < n; i += len)
            for (int k = 0; k < half; k++) {
#if defined(MK_MUTATE) && MK_MUTATE == 11
                double wr = p->c[k * step], wi = -p->s[k * step]; /* planted fault: the inverse runs forward */
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

/* 2-D on split storage, x fastest: re[j * n1 + i]. Columns go eight at a time through a gathered buffer.
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
    /* refuse before touching the data */
    for (int ax = 0; ax < 3; ax++) if (nn[ax] < 1 || (nn[ax] & (nn[ax] - 1))) return -1;
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

/* Any length: mk_fftm, a complex FFT of any length n from 1 to 2^29 with no padding: the transform of
   the n values themselves, for lengths such as 1080, 1920 and 3840 that the radix-2 plan refuses. The same signs as
   mk_fft (forward e^{-i}, neither direction normalized).
   Method: a length whose prime factors are small goes through self-sorting passes (Stockham, decimation in frequency)
   of radices 4, 2, 3, 5 and 7 written out, and any other small prime p by its pairs (outputs t and p - t from the same
   two sums). A length with a large prime factor goes by Bluestein's chirp z: a cyclic convolution of length M >= 2n - 1
   whose only factors are 2, 3, 5 and 7, two transforms of length M and the chirp's transform kept in the plan. The plan
   picks the pairs or the chirp, and M, by a fixed operation count, so a length always gets the same plan and the same
   bits. Twiddles: each e^{-2 pi i k / n} reduced in integers to an angle in [0, pi/4], whose cosine and sine are rounded
   from long double; with a 64-bit long double significand (x86-64) every twiddle is within half an ulp plus 2^-62
   relative of the exact value.
   Accuracy, against the exact transform, eps = DBL_EPSILON: max |error| <= 2 eps (log2 n + 4) max |X| and the rms
   error <= eps (log2 n + 4) rms |X|, the same bound in both directions (test/t_fft.c checks it at every n to 520 and at
   chosen lengths to 2^20).
   Speed, one thread, -O2, against mk_fft_c at N, the next power of 2 >= n: no slower at powers of 2; within 1.5 times
   at lengths 2^a 3^b 5^c 7^d; within 4 times at any n to 2^20.
   mk_fftm_plan(p, n): 0, or -1 (n out of range, no memory; p is then empty). The plan holds about 16 n bytes, a
   Bluestein length 16 (n + 2 M). mk_fftm_free releases it.
   mk_fftm_c(p, a, work, inverse): a[n] in place, work[n] the caller's scratch, so one plan serves many threads. A
   Bluestein length also takes 32 M bytes for the call (malloc, freed before it returns); without that memory it falls
   back to the direct sum, right but slow (n^2 operations), with an error that grows as sqrt n. */
typedef struct {
    int n, nf, f[32];   /* the length; the radices of the passes over L, in their order */
    long long L;        /* the passes' length: n, or Bluestein's M */
    size_t off[33];     /* pass st's twiddles start at tw + off[st]; a pass of another prime r keeps its r roots after them */
    double complex *tw; /* e^{-2 pi i j t / l}, j < l / r, 0 < t < r, pass by pass; one allocation with w and bh */
    double complex *w;  /* Bluestein: the chirp e^{-i pi k^2 / n}, k < n; NULL otherwise */
    double complex *bh; /* Bluestein: the length-M transform of the chirp's conjugate laid out at k and M - k, over M */
} mk_fftm;
enum { MK_FFTM_GMAX_ = 61 }; /* the largest prime taken by its pairs; above it, Bluestein */

static inline void mk_fftm_free(mk_fftm *p)
{
    free(p->tw);
    memset(p, 0, sizeof *p);
}
/* e^{-2 pi i k / n}, 0 <= k < n: the angle reduced in integers to o eighths of a turn and a rest in [0, pi/4] */
static inline double complex mk_fftm_cis_(long long k, long long n)
{
    const long long k8 = 8 * k, o = k8 / n, r = k8 - o * n, a = o & 1 ? n - r : r;
    const long double ph = 0.78539816339744830961566084581987572L * (long double)a / (long double)n;
    const double c = (double)cosl(ph), s = (double)(o & 1 ? -sinl(ph) : sinl(ph));
    switch (((o + 1) / 2) & 3) { /* the angle is (o + 1) / 2 quarter turns, plus ph for an even o, minus it for odd */
    case 0: return CMPLX(c, -s);
    case 1: return CMPLX(-s, -c);
    case 2: return CMPLX(-c, s);
    default: return CMPLX(s, c);
    }
}
/* t[k] = e^{-2 pi i k / L}, k < L: up to an eighth, a quarter or half of them by mk_fftm_cis_, the rest by the exact
   symmetries of the circle */
static inline void mk_fftm_roots_(double complex *t, long long L)
{
    for (long long k = 0; k < L; k++) {
#if defined(MK_MUTATE) && MK_MUTATE == 305
        t[k] = CMPLX(cos(2 * MK_FFT_PI * k / L), -sin(2 * MK_FFT_PI * k / L)); /* planted fault: the twiddles of a double angle */
        continue;
#endif
        if (L % 4 == 0 && k >= L / 4) t[k] = CMPLX(cimag(t[k - L / 4]), -creal(t[k - L / 4])); /* a quarter turn on */
        else if (L % 2 == 0 && k >= L / 2) t[k] = CMPLX(-creal(t[k - L / 2]), -cimag(t[k - L / 2]));
        else if (L % 4 == 0 && 8 * k > L) t[k] = CMPLX(-cimag(t[L / 4 - k]), -creal(t[L / 4 - k])); /* mirrored in pi/4 */
        else if (L % 2 == 0 && 4 * k > L) t[k] = CMPLX(-creal(t[L / 2 - k]), cimag(t[L / 2 - k])); /* mirrored in pi/2 */
        else if (2 * k > L) t[k] = conj(t[L - k]);
        else t[k] = mk_fftm_cis_(k, L);
    }
}
/* the radices of the passes over L: 4s, a 2, then the odd primes upward; their count, or -1 */
static inline int mk_fftm_radices_(long long L, int *f)
{
    int nf = 0;
    long long m = L;
    for (long long q = 4; m > 1; q = q == 4 ? 2 : q == 2 ? 3 : q + 2) {
        if (q > 4 && q * q > m) q = m; /* what is left is prime */
        while (m % q == 0) {
            if (nf == 32) return -1;
            f[nf++] = (int)q, m /= q;
            if (q == 2) break;
        }
    }
    return nf;
}
/* the plan's count of a pass of radix r, per point, in passes of radix 4 (weights fitted to one-thread timings at -O2
   on an x86-64 desktop; they decide plans, never results' accuracy) */
static inline double mk_fftm_rcost_(int r)
{
    switch (r) {
    case 2: return 0.9;
    case 3: return 1.0;
    case 4: return 1.0;
    case 5: return 1.4;
    case 7: return 1.66;
    default: return 2.5 + 0.27 * r;
    }
}
static inline double mk_fftm_count_(long long L)
{
    int f[32], nf = mk_fftm_radices_(L, f);
    double c = 0;
    for (int i = 0; i < nf; i++) c += mk_fftm_rcost_(f[i]);
    return c * (double)L;
}
static inline double complex mk_fftm_mul_(double complex a, double wr, double wi)
{
    return CMPLX(creal(a) * wr - cimag(a) * wi, creal(a) * wi + cimag(a) * wr);
}
/* One pass of radix r: the r-point transforms of x[q + s (j + m u)], u < r, each output t times e^{sg 2 pi i j t / (m r)},
   to y[q + s (t + r j)]. W: the pass's twiddles (forward); sg -1 forward, +1 inverse. The last pass (m = 1) may run in
   place (y = x): each transform reads all its inputs before it writes. */
static inline void mk_fftm_p2_(const double complex *x, double complex *y, size_t m, size_t s, const double complex *W, double sg)
{
    const size_t sm = s * m;
    for (size_t j = 0; j < m; j++) {
        const double complex *in = x + s * j;
        double complex *out = y + 2 * s * j;
        const double w1r = creal(W[j]), w1i = -sg * cimag(W[j]);
        for (size_t q = 0; q < s; q++) {
            const double complex a0 = in[q], a1 = in[q + sm];
            const double complex b1 = CMPLX(creal(a0) - creal(a1), cimag(a0) - cimag(a1));
            out[q] = CMPLX(creal(a0) + creal(a1), cimag(a0) + cimag(a1));
            out[q + s] = j ? mk_fftm_mul_(b1, w1r, w1i) : b1;
        }
    }
}
static inline void mk_fftm_p3_(const double complex *x, double complex *y, size_t m, size_t s, const double complex *W, double sg)
{
    const size_t sm = s * m;
    const double h3 = 0.866025403784438646763723170752936; /* sin(2 pi / 3) */
    for (size_t j = 0; j < m; j++) {
        const double complex *in = x + s * j, *wt = W + 2 * j;
        double complex *out = y + 3 * s * j;
        const double w1r = creal(wt[0]), w1i = -sg * cimag(wt[0]), w2r = creal(wt[1]), w2i = -sg * cimag(wt[1]);
        for (size_t q = 0; q < s; q++) {
            const double complex a0 = in[q], a1 = in[q + sm], a2 = in[q + 2 * sm];
            const double tr = creal(a1) + creal(a2), ti = cimag(a1) + cimag(a2);
#if defined(MK_MUTATE) && MK_MUTATE == 239
            /* planted fault: the difference's sign turned */
            const double dr = sg * h3 * (creal(a2) - creal(a1)), di = sg * h3 * (cimag(a2) - cimag(a1));
#else
            const double dr = sg * h3 * (creal(a1) - creal(a2)), di = sg * h3 * (cimag(a1) - cimag(a2));
#endif
            const double cr = creal(a0) - tr / 2, ci = cimag(a0) - ti / 2;
            const double complex b1 = CMPLX(cr - di, ci + dr), b2 = CMPLX(cr + di, ci - dr);
            out[q] = CMPLX(creal(a0) + tr, cimag(a0) + ti);
            out[q + s] = j ? mk_fftm_mul_(b1, w1r, w1i) : b1;
            out[q + 2 * s] = j ? mk_fftm_mul_(b2, w2r, w2i) : b2;
        }
    }
}
static inline void mk_fftm_p4_(const double complex *x, double complex *y, size_t m, size_t s, const double complex *W, double sg)
{
    const size_t sm = s * m;
#if defined(MK_MUTATE) && MK_MUTATE == 238
    const double e = -1; /* planted fault: the quarter turn always the forward one */
#else
    const double e = sg;
#endif
    for (size_t j = 0; j < m; j++) {
        const double complex *in = x + s * j, *wt = W + 3 * j;
        double complex *out = y + 4 * s * j;
        const double w1r = creal(wt[0]), w1i = -sg * cimag(wt[0]), w2r = creal(wt[1]), w2i = -sg * cimag(wt[1]);
        const double w3r = creal(wt[2]), w3i = -sg * cimag(wt[2]);
        for (size_t q = 0; q < s; q++) {
            const double complex a0 = in[q], a1 = in[q + sm], a2 = in[q + 2 * sm], a3 = in[q + 3 * sm];
            const double pr = creal(a0) + creal(a2), pi = cimag(a0) + cimag(a2), mr = creal(a0) - creal(a2), mi = cimag(a0) - cimag(a2);
            const double qr = creal(a1) + creal(a3), qi = cimag(a1) + cimag(a3), dr = creal(a1) - creal(a3), di = cimag(a1) - cimag(a3);
            const double complex b1 = CMPLX(mr - e * di, mi + e * dr), b2 = CMPLX(pr - qr, pi - qi), b3 = CMPLX(mr + e * di, mi - e * dr);
            out[q] = CMPLX(pr + qr, pi + qi);
            if (j) {
                out[q + s] = mk_fftm_mul_(b1, w1r, w1i);
                out[q + 2 * s] = mk_fftm_mul_(b2, w2r, w2i);
                out[q + 3 * s] = mk_fftm_mul_(b3, w3r, w3i);
            } else out[q + s] = b1, out[q + 2 * s] = b2, out[q + 3 * s] = b3;
        }
    }
}
static inline void mk_fftm_p5_(const double complex *x, double complex *y, size_t m, size_t s, const double complex *W, double sg)
{
    const size_t sm = s * m;
    const double c1 = 0.309016994374947424102293417182819, c2 = -0.809016994374947424102293417182819; /* cos(2 pi k / 5) */
    const double t1 = sg * 0.951056516295153572116439333379382, t2 = sg * 0.587785252292473129168705954639073; /* sg sin(2 pi k / 5) */
    for (size_t j = 0; j < m; j++) {
        const double complex *in = x + s * j, *wt = W + 4 * j;
        double complex *out = y + 5 * s * j;
        const double w1r = creal(wt[0]), w1i = -sg * cimag(wt[0]), w2r = creal(wt[1]), w2i = -sg * cimag(wt[1]);
        const double w3r = creal(wt[2]), w3i = -sg * cimag(wt[2]), w4r = creal(wt[3]), w4i = -sg * cimag(wt[3]);
        for (size_t q = 0; q < s; q++) {
            const double complex a0 = in[q], a1 = in[q + sm], a2 = in[q + 2 * sm], a3 = in[q + 3 * sm], a4 = in[q + 4 * sm];
            const double s1r = creal(a1) + creal(a4), s1i = cimag(a1) + cimag(a4), s2r = creal(a2) + creal(a3), s2i = cimag(a2) + cimag(a3);
            const double d1r = creal(a1) - creal(a4), d1i = cimag(a1) - cimag(a4), d2r = creal(a2) - creal(a3), d2i = cimag(a2) - cimag(a3);
            const double A1r = creal(a0) + c1 * s1r + c2 * s2r, A1i = cimag(a0) + c1 * s1i + c2 * s2i;
            const double A2r = creal(a0) + c2 * s1r + c1 * s2r, A2i = cimag(a0) + c2 * s1i + c1 * s2i;
            const double B1r = t1 * d1r + t2 * d2r, B1i = t1 * d1i + t2 * d2i;
#if defined(MK_MUTATE) && MK_MUTATE == 300
            const double B2r = t2 * d1r + t1 * d2r, B2i = t2 * d1i + t1 * d2i; /* planted fault: a sign of the second pair */
#else
            const double B2r = t2 * d1r - t1 * d2r, B2i = t2 * d1i - t1 * d2i;
#endif
            const double complex b1 = CMPLX(A1r - B1i, A1i + B1r), b2 = CMPLX(A2r - B2i, A2i + B2r), b3 = CMPLX(A2r + B2i, A2i - B2r), b4 = CMPLX(A1r + B1i, A1i - B1r);
            out[q] = CMPLX(creal(a0) + s1r + s2r, cimag(a0) + s1i + s2i);
            if (j) {
                out[q + s] = mk_fftm_mul_(b1, w1r, w1i), out[q + 2 * s] = mk_fftm_mul_(b2, w2r, w2i);
                out[q + 3 * s] = mk_fftm_mul_(b3, w3r, w3i), out[q + 4 * s] = mk_fftm_mul_(b4, w4r, w4i);
            } else out[q + s] = b1, out[q + 2 * s] = b2, out[q + 3 * s] = b3, out[q + 4 * s] = b4;
        }
    }
}
static inline void mk_fftm_p7_(const double complex *x, double complex *y, size_t m, size_t s, const double complex *W, double sg)
{
    const size_t sm = s * m;
    const double c1 = 0.623489801858733530525004884004240, c2 = -0.222520933956314404288902564496795,
                 c3 = -0.900968867902419126236102319507445; /* cos(2 pi k / 7) */
    const double t1 = sg * 0.781831482468029808708444526674057, t2 = sg * 0.974927912181823607018131682993931,
                 t3 = sg * 0.433883739117558120475768332848358; /* sg sin(2 pi k / 7) */
    for (size_t j = 0; j < m; j++) {
        const double complex *in = x + s * j, *wt = W + 6 * j;
        double complex *out = y + 7 * s * j;
        const double w1r = creal(wt[0]), w1i = -sg * cimag(wt[0]), w2r = creal(wt[1]), w2i = -sg * cimag(wt[1]);
        const double w3r = creal(wt[2]), w3i = -sg * cimag(wt[2]), w4r = creal(wt[3]), w4i = -sg * cimag(wt[3]);
        const double w5r = creal(wt[4]), w5i = -sg * cimag(wt[4]), w6r = creal(wt[5]), w6i = -sg * cimag(wt[5]);
        for (size_t q = 0; q < s; q++) {
            const double complex a0 = in[q], a1 = in[q + sm], a2 = in[q + 2 * sm], a3 = in[q + 3 * sm], a4 = in[q + 4 * sm], a5 = in[q + 5 * sm], a6 = in[q + 6 * sm];
            const double s1r = creal(a1) + creal(a6), s1i = cimag(a1) + cimag(a6), d1r = creal(a1) - creal(a6), d1i = cimag(a1) - cimag(a6);
            const double s2r = creal(a2) + creal(a5), s2i = cimag(a2) + cimag(a5), d2r = creal(a2) - creal(a5), d2i = cimag(a2) - cimag(a5);
            const double s3r = creal(a3) + creal(a4), s3i = cimag(a3) + cimag(a4), d3r = creal(a3) - creal(a4), d3i = cimag(a3) - cimag(a4);
            const double A1r = creal(a0) + c1 * s1r + c2 * s2r + c3 * s3r, A1i = cimag(a0) + c1 * s1i + c2 * s2i + c3 * s3i;
#if defined(MK_MUTATE) && MK_MUTATE == 301
            const double A2r = creal(a0) + c2 * s1r + c1 * s2r + c3 * s3r, A2i = cimag(a0) + c2 * s1i + c1 * s2i + c3 * s3i; /* planted fault: two cosines swapped */
#else
            const double A2r = creal(a0) + c2 * s1r + c3 * s2r + c1 * s3r, A2i = cimag(a0) + c2 * s1i + c3 * s2i + c1 * s3i;
#endif
            const double A3r = creal(a0) + c3 * s1r + c1 * s2r + c2 * s3r, A3i = cimag(a0) + c3 * s1i + c1 * s2i + c2 * s3i;
            const double B1r = t1 * d1r + t2 * d2r + t3 * d3r, B1i = t1 * d1i + t2 * d2i + t3 * d3i;
            const double B2r = t2 * d1r - t3 * d2r - t1 * d3r, B2i = t2 * d1i - t3 * d2i - t1 * d3i;
            const double B3r = t3 * d1r - t1 * d2r + t2 * d3r, B3i = t3 * d1i - t1 * d2i + t2 * d3i;
            const double complex b1 = CMPLX(A1r - B1i, A1i + B1r), b2 = CMPLX(A2r - B2i, A2i + B2r), b3 = CMPLX(A3r - B3i, A3i + B3r);
            const double complex b4 = CMPLX(A3r + B3i, A3i - B3r), b5 = CMPLX(A2r + B2i, A2i - B2r), b6 = CMPLX(A1r + B1i, A1i - B1r);
            out[q] = CMPLX(creal(a0) + s1r + s2r + s3r, cimag(a0) + s1i + s2i + s3i);
            if (j) {
                out[q + s] = mk_fftm_mul_(b1, w1r, w1i), out[q + 2 * s] = mk_fftm_mul_(b2, w2r, w2i), out[q + 3 * s] = mk_fftm_mul_(b3, w3r, w3i);
                out[q + 4 * s] = mk_fftm_mul_(b4, w4r, w4i), out[q + 5 * s] = mk_fftm_mul_(b5, w5r, w5i), out[q + 6 * s] = mk_fftm_mul_(b6, w6r, w6i);
            } else out[q + s] = b1, out[q + 2 * s] = b2, out[q + 3 * s] = b3, out[q + 4 * s] = b4, out[q + 5 * s] = b5, out[q + 6 * s] = b6;
        }
    }
}
/* any other odd prime r <= MK_FFTM_GMAX_, by its pairs: y_t, y_{r - t} = a0 + sum S_u cos(2 pi u t / r) +- i sg sum D_u
   sin(2 pi u t / r), S_u and D_u the sum and difference of a_u and a_{r - u} */
static inline void mk_fftm_pg_(const double complex *x, double complex *y, size_t m, size_t s, const double complex *W, int r, double sg)
{
    const size_t sm = s * m;
    const double complex *R = W + m * (size_t)(r - 1); /* the roots e^{-2 pi i k / r}, k < r */
    const int h = (r - 1) / 2;
    double cs[MK_FFTM_GMAX_], sn[MK_FFTM_GMAX_], Sr[MK_FFTM_GMAX_ / 2 + 1], Si[MK_FFTM_GMAX_ / 2 + 1], Dr[MK_FFTM_GMAX_ / 2 + 1], Di[MK_FFTM_GMAX_ / 2 + 1];
    for (int k = 0; k < r; k++) cs[k] = creal(R[k]), sn[k] = -sg * cimag(R[k]); /* sn: sg sin(2 pi k / r) */
    for (size_t j = 0; j < m; j++) {
        const double complex *in = x + s * j, *wt = W + (size_t)(r - 1) * j;
        double complex *out = y + (size_t)r * s * j;
        for (size_t q = 0; q < s; q++) {
            const double complex a0 = in[q];
            double y0r = creal(a0), y0i = cimag(a0);
            for (int u = 1; u <= h; u++) {
                const double complex au = in[q + (size_t)u * sm], av = in[q + (size_t)(r - u) * sm];
                Sr[u] = creal(au) + creal(av), Si[u] = cimag(au) + cimag(av), Dr[u] = creal(au) - creal(av), Di[u] = cimag(au) - cimag(av);
                y0r += Sr[u], y0i += Si[u];
            }
            for (int t = 1; t <= h; t += 2) { /* two outputs a sweep, t and v = t + 1 (t again when h is odd): eight sums */
                const int v = t < h ? t + 1 : t;
                double Ar = creal(a0), Ai = cimag(a0), Br = 0, Bi = 0, Er = creal(a0), Ei = cimag(a0), Fr = 0, Fi = 0;
                for (int u = 1, k = t, l = v; u <= h; u++, k = k + t < r ? k + t : k + t - r, l = l + v < r ? l + v : l + v - r) {
#if defined(MK_MUTATE) && MK_MUTATE == 302
                    const int kk = k + 1 < r ? k + 1 : 0, ll = l + 1 < r ? l + 1 : 0; /* planted fault: the roots read one place on */
#else
                    const int kk = k, ll = l;
#endif
                    Ar += Sr[u] * cs[kk], Ai += Si[u] * cs[kk], Br += Dr[u] * sn[kk], Bi += Di[u] * sn[kk];
                    Er += Sr[u] * cs[ll], Ei += Si[u] * cs[ll], Fr += Dr[u] * sn[ll], Fi += Di[u] * sn[ll];
                }
                const double complex bt = CMPLX(Ar - Bi, Ai + Br), bu = CMPLX(Ar + Bi, Ai - Br), bv = CMPLX(Er - Fi, Ei + Fr), bw = CMPLX(Er + Fi, Ei - Fr);
                out[q + (size_t)t * s] = j ? mk_fftm_mul_(bt, creal(wt[t - 1]), -sg * cimag(wt[t - 1])) : bt;
                out[q + (size_t)(r - t) * s] = j ? mk_fftm_mul_(bu, creal(wt[r - t - 1]), -sg * cimag(wt[r - t - 1])) : bu;
                if (v != t) {
                    out[q + (size_t)v * s] = j ? mk_fftm_mul_(bv, creal(wt[v - 1]), -sg * cimag(wt[v - 1])) : bv;
                    out[q + (size_t)(r - v) * s] = j ? mk_fftm_mul_(bw, creal(wt[r - v - 1]), -sg * cimag(wt[r - v - 1])) : bw;
                }
            }
            out[q] = CMPLX(y0r, y0i);
        }
    }
}
/* the passes over L in place on a[L], w[L] the scratch */
static inline void mk_fftm_run_(const mk_fftm *p, double complex *a, double complex *w, int inverse)
{
    const double sg = inverse ? 1 : -1;
    double complex *x = a, *y = w;
    size_t l = (size_t)p->L, s = 1;
    for (int st = 0; st < p->nf; st++) {
        const int r = p->f[st];
        const size_t m = l / (size_t)r;
        const double complex *W = p->tw + p->off[st];
        double complex *out = st == p->nf - 1 && p->nf % 2 ? x : y; /* an odd count of passes: the last one (m = 1) in place */
        switch (r) {
        case 2: mk_fftm_p2_(x, out, m, s, W, sg); break;
        case 3: mk_fftm_p3_(x, out, m, s, W, sg); break;
        case 4: mk_fftm_p4_(x, out, m, s, W, sg); break;
        case 5: mk_fftm_p5_(x, out, m, s, W, sg); break;
        case 7: mk_fftm_p7_(x, out, m, s, W, sg); break;
        default: mk_fftm_pg_(x, out, m, s, W, r, sg); break;
        }
        if (out == y) y = x, x = out;
        l = m, s *= (size_t)r;
    }
}
/* the passes' length for n (n itself, or Bluestein's M at the lowest count from 2n - 1 to the next power of 2) and,
   in *cost, the count of a transform: what the plan picks by, and what a caller can compare lengths by */
static inline long long mk_fftm_len_(int n, double *cost)
{
    int f[32], nf = mk_fftm_radices_(n, f), big = 1;
    for (int i = 0; i < nf; i++) big = f[i] > big ? f[i] : big;
    const double cd = big > MK_FFTM_GMAX_ || nf < 0 ? HUGE_VAL : mk_fftm_count_(n);
    long long L = n;
    double c0 = cd;
    if (big > 7) { /* the pairs, or Bluestein */
#if defined(MK_MUTATE) && MK_MUTATE == 309
        const long long lo = n; /* planted fault: a convolution only as long as the transform */
#else
        const long long lo = 2LL * n - 1;
#endif
        long long hi = 1, best = 0;
        double cb = HUGE_VAL;
        while (hi < lo) hi *= 2;
        for (long long o7 = 1; o7 <= hi; o7 *= 7)
            for (long long o5 = o7; o5 <= hi; o5 *= 5)
                for (long long o3 = o5; o3 <= hi; o3 *= 3) {
                    long long M = o3;
                    while (M < lo) M *= 2;
                    /* two transforms of M; the chirp in and out, the product and the zeros: 3 n + M at 0.8 each (fitted) */
                    const double c = M > hi ? HUGE_VAL : 2 * mk_fftm_count_(M) + 0.8 * (3.0 * n + (double)M);
                    if (c < cb || (c == cb && M < best)) cb = c, best = M;
                }
        if (cb < cd) L = best, c0 = cb;
    }
    if (cost) *cost = c0;
    return L;
}
static inline int mk_fftm_plan(mk_fftm *p, int n)
{
    memset(p, 0, sizeof *p);
    if (n < 1 || n > (1 << 29)) return -1;
    int f[32];
    const long long L = mk_fftm_len_(n, NULL);
    const int nf = mk_fftm_radices_(L, f);
    if (nf < 0) return -1;
    size_t nt = 0;
    for (long long st = 0, l = L; st < nf; st++) nt += (size_t)(l / f[st]) * (size_t)(f[st] - 1) + (f[st] > 7 ? (size_t)f[st] : 0), l /= f[st];
    const size_t nb = L != n ? (size_t)n + (size_t)L : 0, nr = L != n ? 2 * (size_t)L : (size_t)L; /* Bluestein: the 2n roots, then B and its scratch */
    double complex *tw = malloc(sizeof *tw * (nt + nb + 1)), *t = malloc(sizeof *t * nr);
    if (!tw || !t) { free(tw), free(t); return -1; }
    p->tw = tw, p->n = n, p->nf = nf, p->L = L;
    memcpy(p->f, f, sizeof f);
    if (L != n) { /* the chirp from the 2n-th roots: e^{-i pi k^2 / n} = e^{-2 pi i (k^2 mod 2n) / 2n} */
        p->w = tw + nt, p->bh = p->w + n;
        mk_fftm_roots_(t, 2LL * n);
        for (long long k = 0; k < n; k++) {
#if defined(MK_MUTATE) && MK_MUTATE == 303
            p->w[k] = t[k * k % n]; /* planted fault: the square taken mod n */
#else
            p->w[k] = t[k * k % (2LL * n)];
#endif
        }
    }
    mk_fftm_roots_(t, L);
    size_t o = 0, sl = 1;
    for (int st = 0; st < nf; st++) {
        const size_t r = (size_t)f[st], m = (size_t)L / sl / r;
        p->off[st] = o;
        for (size_t j = 0; j < m; j++)
            for (size_t u = 1; u < r; u++) tw[o++] = t[j * u * sl];
        if (r > 7)
            for (size_t k = 0; k < r; k++) tw[o++] = t[k * ((size_t)L / r)];
        sl *= r;
    }
    p->off[nf] = o;
    if (L != n) { /* the chirp's conjugate at k and M - k, transformed and divided by M */
        double complex *B = t, *S = t + L;
        for (long long k = 0; k < L; k++) B[k] = 0;
        for (long long k = 0; k < n; k++) B[k] = conj(p->w[k]), B[k ? L - k : 0] = conj(p->w[k]);
        mk_fftm_run_(p, B, S, 0);
        for (long long k = 0; k < L; k++) {
#if defined(MK_MUTATE) && MK_MUTATE == 307
            p->bh[k] = B[k]; /* planted fault: not divided by M */
#else
            p->bh[k] = CMPLX(creal(B[k]) / (double)L, cimag(B[k]) / (double)L);
#endif
        }
    }
    free(t);
    return 0;
}
/* the scratch mk_fftm_x_ takes (double complex values): 2 M for a Bluestein length, else 0 */
static inline size_t mk_fftm_scratch_(const mk_fftm *p) { return p->w ? 2 * (size_t)p->L : 0; }
/* mk_fftm_c with the caller's scratch: scr[mk_fftm_scratch_(p)], or NULL for the direct sum */
static inline void mk_fftm_x_(const mk_fftm *p, double complex *a, double complex *work, double complex *scr, int inverse)
{
    if (!p->w) { mk_fftm_run_(p, a, work, inverse); return; }
    const size_t n = (size_t)p->n, M = (size_t)p->L;
#if defined(MK_MUTATE) && MK_MUTATE == 304
    const double cj = 1; /* planted fault: the inverse without its conjugates */
#else
    const double cj = inverse ? -1 : 1; /* the inverse is the conjugate of the forward transform of the conjugate */
#endif
    if (!scr) { /* no memory: the convolution as its direct sum, into work */
        for (size_t j = 0; j < n; j++) a[j] = mk_fftm_mul_(CMPLX(creal(a[j]), cj * cimag(a[j])), creal(p->w[j]), cimag(p->w[j]));
        for (size_t k = 0; k < n; k++) {
            double sr = 0, si = 0;
            for (size_t j = 0; j < n; j++) {
                const double complex c = p->w[k > j ? k - j : j - k];
                sr += creal(a[j]) * creal(c) + cimag(a[j]) * cimag(c), si += cimag(a[j]) * creal(c) - creal(a[j]) * cimag(c);
            }
#if defined(MK_MUTATE) && MK_MUTATE == 306
            work[k] = CMPLX(sr, cj * si); /* planted fault: the direct sum without the chirp out */
#else
            const double complex z = mk_fftm_mul_(CMPLX(sr, si), creal(p->w[k]), cimag(p->w[k]));
            work[k] = CMPLX(creal(z), cj * cimag(z));
#endif
        }
        memcpy(a, work, sizeof *a * n);
        return;
    }
    double complex *A = scr, *S = scr + M;
    for (size_t j = 0; j < n; j++) A[j] = mk_fftm_mul_(CMPLX(creal(a[j]), cj * cimag(a[j])), creal(p->w[j]), cimag(p->w[j]));
    for (size_t j = n; j < M; j++) A[j] = 0;
    mk_fftm_run_(p, A, S, 0);
    for (size_t k = 0; k < M; k++) A[k] = mk_fftm_mul_(A[k], creal(p->bh[k]), cimag(p->bh[k]));
    mk_fftm_run_(p, A, S, 1);
    for (size_t k = 0; k < n; k++) {
        const double complex z = mk_fftm_mul_(A[k], creal(p->w[k]), cimag(p->w[k]));
        a[k] = CMPLX(creal(z), cj * cimag(z));
    }
}
static inline void mk_fftm_c(const mk_fftm *p, double complex *a, double complex *work, int inverse)
{
    if (p->n <= 4 && p->n != 3) { /* n 1, 2 and 4: the same sums as their one pass, without its calls */
        const double e = inverse ? 1 : -1;
        if (p->n == 2) {
            const double complex a0 = a[0], a1 = a[1];
            a[0] = CMPLX(creal(a0) + creal(a1), cimag(a0) + cimag(a1)), a[1] = CMPLX(creal(a0) - creal(a1), cimag(a0) - cimag(a1));
        } else if (p->n == 4) {
            const double pr = creal(a[0]) + creal(a[2]), pi = cimag(a[0]) + cimag(a[2]), mr = creal(a[0]) - creal(a[2]), mi = cimag(a[0]) - cimag(a[2]);
            const double qr = creal(a[1]) + creal(a[3]), qi = cimag(a[1]) + cimag(a[3]), dr = creal(a[1]) - creal(a[3]), di = cimag(a[1]) - cimag(a[3]);
            a[0] = CMPLX(pr + qr, pi + qi), a[1] = CMPLX(mr - e * di, mi + e * dr), a[2] = CMPLX(pr - qr, pi - qi), a[3] = CMPLX(mr + e * di, mi - e * dr);
        }
        return;
    }
    if (!p->w) { mk_fftm_run_(p, a, work, inverse); return; }
    double complex *scr = malloc(sizeof *scr * mk_fftm_scratch_(p));
    mk_fftm_x_(p, a, work, scr, inverse);
    free(scr);
}
#endif
