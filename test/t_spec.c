/* mk_spec.h: Parseval and amplitude scaling with every window, Welch's density of white noise, the coherence and
   transfer function of a known linear system, Blackman-Harris's published sidelobe, STFT zero padding, the refusals. */
#include "../include/mk_spec.h"
#include "../include/mk_rng.h"
#include "t.h"

int main(void)
{
    CHECK(fabs(mk_bessel_i0(5) - 27.239871823604442) < 1e-12 && fabs(mk_bessel_i0(1) - 1.2660658777520082) < 1e-15,
          "I0(5) = %.15g, I0(1) = %.16g (27.239871823604442, 1.2660658777520082)", mk_bessel_i0(5), mk_bessel_i0(1));
    enum { NF = 256 };
    static double x[NF], p[NF / 2 + 1];
    mk_rng g;
    mk_rng_seed(&g, 5, 0);
    for (int i = 0; i < NF; i++) x[i] = mk_rng_gauss(&g) + 0.3;
    int pbad = 0, abad = 0;
    double worst = 0, aw = 0;
    for (int wt = 0; wt < 5; wt++) {
        mk_spec s;
        mk_spec_init(&s, NF, wt, 8.6, 1000);
        mk_welch(&s, x, NF, 0, 0, MK_SPEC_DENSITY, p);
        double sum = 0, ms = 0;
        for (int k = 0; k <= NF / 2; k++) sum += p[k] * s.fs / NF;
        for (int i = 0; i < NF; i++) ms += x[i] * x[i] * s.w[i] * s.w[i];
        ms /= s.s2;
        worst = fmax(worst, fabs(sum - ms) / ms);
        if (fabs(sum - ms) > 1e-12 * ms) pbad++;
        /* a sine of amplitude 3 centred on bin 40 */
        double y[NF];
        for (int i = 0; i < NF; i++) y[i] = 3 * sin(2 * T_PI * 40 * i / NF + 0.7);
        mk_welch(&s, y, NF, 0, 0, MK_SPEC_POWER, p);
        double a = sqrt(2 * p[40]);
        aw = fmax(aw, fabs(a - 3));
        if (fabs(a - 3) > 1e-5) abad++; /* Kaiser has no exact spectral zeros: its negative-frequency image leaks about 1e-6 in */
        mk_spec_free(&s);
    }
    CHECK(pbad == 0, "Parseval, one segment, every window: density summed over frequency = the windowed mean square (worst %.2g)", worst);
    CHECK(abad == 0, "power scaling: a sine's amplitude 3 recovered from its bin with every window, to 1e-5 (worst %.2g, Kaiser's leakage)", aw);

    long n = 1L << 18;
    double *w = malloc(sizeof(double) * (size_t)n), *z = malloc(sizeof(double) * (size_t)n);
    for (long i = 0; i < n; i++) w[i] = 2 * mk_rng_gauss(&g);
    mk_spec s;
    mk_spec_init(&s, NF, MK_WIN_HANN, 0, 1000);
    long seg = mk_welch(&s, w, n, NF / 2, 1, MK_SPEC_DENSITY, p);
    double lvl = 0;
    for (int k = 1; k < NF / 2; k++) lvl += p[k];
    lvl /= NF / 2 - 1;
    CHECK(fabs(lvl - 0.008) < 0.01 * 0.008, "Welch, Hann, half overlap, %ld segments: white noise of variance 4 at 1 kHz reads %.6f /Hz (2 s^2 / fs = 0.008)", seg, lvl);

    /* y = x + 0.5 x[-1] + noise of sd 0.5: coherence |H|^2 / (|H|^2 + 1/16) (x has variance 4), H = 1 + 0.5 e^-iw */
    for (long i = 0; i < n; i++) z[i] = w[i] + (i ? 0.5 * w[i - 1] : 0) + 0.5 * mk_rng_gauss(&g);
    static double complex pxy[NF / 2 + 1], H[NF / 2 + 1];
    static double pxx[NF / 2 + 1], pyy[NF / 2 + 1], coh[NF / 2 + 1];
    mk_csd(&s, w, z, n, NF / 2, 1, pxy, pxx, pyy);
    mk_coherence(pxy, pxx, pyy, NF / 2 + 1, coh);
    mk_transfer(pxy, pxx, NF / 2 + 1, H);
    double ce = 0, he = 0;
    for (int k = 1; k < NF / 2; k++) {
        double om = 2 * T_PI * k / NF, h2 = 1.25 + cos(om), want = 4 * h2 / (4 * h2 + 0.25);
        double complex hw = 1 + 0.5 * CMPLX(cos(om), -sin(om));
        ce = fmax(ce, fabs(coh[k] - want)), he = fmax(he, cabs(H[k] - hw));
    }
    CHECK(ce < 0.03 && he < 0.05, "coherence of a known linear system with noise: worst error %.3f; transfer function worst error %.3f", ce, he);

    /* Blackman-Harris's highest sidelobe: about -92 dB (Harris 1978) */
    enum { WN = 64, Z = 8192 };
    double bw[WN];
    mk_window(MK_WIN_BLACKMAN_HARRIS, WN, 0, bw);
    static double re[Z], im[Z];
    for (int i = 0; i < Z; i++) re[i] = i < WN ? bw[i] : 0, im[i] = 0;
    mk_fft pz;
    mk_fft_plan(&pz, Z);
    mk_fft_ri(&pz, re, im, 0);
    double dc = hypot(re[0], im[0]), side = 0;
    for (int k = (int)(4.2 * Z / WN); k <= Z / 2; k++) side = fmax(side, hypot(re[k], im[k]));
    double db = 20 * log10(side / dc);
    CHECK(db < -90, "Blackman-Harris highest sidelobe %.1f dB (published -92)", db);
    /* an STFT frame hanging off the end equals the same frame with the missing samples set to 0 */
    static double complex F1[NF / 2 + 1], F2[NF / 2 + 1];
    double pad[NF];
    for (int i = 0; i < NF; i++) pad[i] = i < 100 ? w[n - 100 + i] : 0;
    mk_stft_frame(&s, w, n, n - 100, 0, F1);
    mk_stft_frame(&s, pad, NF, 0, 0, F2);
    double fe = 0;
    for (int k = 0; k <= NF / 2; k++) fe = fmax(fe, cabs(F1[k] - F2[k]));
    CHECK(fe == 0, "an STFT frame past the end of the signal is zero-padded (difference %.2g)", fe);
    free(w), free(z), mk_fft_free(&pz), mk_spec_free(&s);
    { /* an unknown window, or an fs that is not positive and finite, is refused (-1) */
        mk_spec b;
        const int r7 = mk_spec_init(&b, 64, 7, 0, 1000);
        mk_spec_free(&b);
        const int rm = mk_spec_init(&b, 64, -1, 0, 1000);
        mk_spec_free(&b);
        const int r0 = mk_spec_init(&b, 64, MK_WIN_HANN, 0, 0);
        mk_spec_free(&b);
        const int rn = mk_spec_init(&b, 64, MK_WIN_HANN, 0, NAN);
        mk_spec_free(&b);
        const int ri = mk_spec_init(&b, 64, MK_WIN_HANN, 0, INFINITY);
        mk_spec_free(&b);
        const int rk = mk_spec_init(&b, 64, MK_WIN_KAISER, 8.6, 48000);
        mk_spec_free(&b);
        CHECK(r7 == -1 && rm == -1 && r0 == -1 && rn == -1 && ri == -1 && rk == 0, "spec init refuses window 7 and -1, fs 0, NaN and inf (%d %d %d %d %d); Kaiser at 48 kHz is 0 (%d)", r7,
              rm, r0, rn, ri, rk);
    }
    return wrong;
}
