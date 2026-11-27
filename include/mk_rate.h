/* mathkit: how a rate, a relaxation time or a viscosity follows temperature. Each law takes its reference temperature as
   an argument and is 1 there (the WLF law, a log10, is 0), so two callers that disagree about it say so in their calls.

   mk_rate_q(T, Tref, q, dT)           q^((T - Tref) / dT): a factor q per dT of warming. A Q10 is q with dT 10.
   mk_rate_arrhenius(T, Tref, theta)   exp(-theta (1/T - 1/Tref)), theta = E / R in K (Arrhenius 1889). theta < 0 for what
                                       falls with warming, such as a liquid's viscosity (Andrade 1930). Tref INFINITY gives
                                       exp(-theta / T), the factor of a rate constant A exp(-E / (R T)). With a boiling
                                       point as Tref and theta = latent heat / R it is the vapour pressure over the pressure
                                       there (Clausius and Clapeyron, the latent heat constant).
   mk_rate_wlf_log10(T, Tref, C1, C2)  log10 of the shift factor a_T = tau(T) / tau(Tref) of a polymer above its glass
                                       temperature, -C1 (T - Tref) / (C2 + T - Tref) (Williams, Landel and Ferry 1955).
                                       A log10, so that the ratio of two shifts is one pow of a difference.
   mk_rate_walther(T, A, B)            an oil's kinematic viscosity in mm^2/s from log10 log10(nu + 0.7) = A - B log10 T
                                       (Walther 1928; ASTM D341, whose 0.7 serves from 2 mm^2/s up);
                                       mk_rate_walther_fit takes A and B from two measured points.

   T and Tref are in K for Arrhenius and Walther; mk_rate_q and the WLF law use only their difference, so either scale
   serves. A factor and its inverse: mk_rate_q and mk_rate_arrhenius with T and Tref swapped give the inverse, to rounding.

   A history of temperatures: dt at T does what dt times the factor does at Tref, so the equivalent time at Tref is the sum
   of factor(T) dt over the history's pieces, and the time at T that equals t at Tref is t / factor(T). a_T is a factor on
   a time, so for it the two change places: the equivalent time is the sum of dt / a_T.

   A NaN argument, T or Tref not above 0 in a law in K, q not above 0, dT 0, C2 not above 0 or T INFINITY in the WLF law
   gives NaN. Each call costs one exp or pow (the WLF law none, Walther two pow and a log10) and has no state. Header
   only, C23. */
#ifndef MK_RATE_H
#define MK_RATE_H
#include <math.h>

#define MK_RATE_GAS_R 8.31446261815324 /* J / (mol K), exact: theta = E / MK_RATE_GAS_R for E in J/mol */

/* the factor q^((T - Tref) / dT) */
static inline double mk_rate_q(double T, double Tref, double q, double dT)
{
#if defined(MK_MUTATE) && MK_MUTATE == 400
    const double x = (T - Tref) / (dT != 0 ? 10 : dT); /* planted fault: always per 10 */
#else
    const double x = (T - Tref) / dT;
#endif
#if defined(MK_MUTATE) && MK_MUTATE == 401
    return q > 0 && dT != 0 ? pow(q, x) : NAN; /* planted fault: pow(1, NaN) is 1 */
#else
    return q > 0 && dT != 0 && !isnan(x) ? pow(q, x) : NAN;
#endif
}

/* the factor exp(-theta (1/T - 1/Tref)); T, Tref and theta in K, Tref may be INFINITY */
static inline double mk_rate_arrhenius(double T, double Tref, double theta)
{
#if defined(MK_MUTATE) && MK_MUTATE == 402
    return T > 0 && Tref > 0 ? exp(-theta * (1 / T)) : NAN; /* planted fault: the reference dropped */
#elif defined(MK_MUTATE) && MK_MUTATE == 403
    return exp(-theta * (1 / T - 1 / Tref)); /* planted fault: no refusal at and below 0 K */
#else
    if (!(T > 0 && Tref > 0)) return NAN;
    const double d = 1 / T - 1 / Tref;
#if !(defined(MK_MUTATE) && MK_MUTATE == 408) /* planted fault: a reciprocal's overflow left */
    /* a reciprocal overflows under 2^-1024 K: both temperatures 2^600 times as large, and the power given back last */
    if (!isfinite(d)) return exp(ldexp(-theta * (1 / (T * 0x1p600) - 1 / (Tref * 0x1p600)), 600));
#endif
    return exp(-theta * d);
#endif
}

/* log10 a_T; below Tref - 0.99 C2 (the law's pole is Tref - C2) it is held at its value there, 99 C1 */
static inline double mk_rate_wlf_log10(double T, double Tref, double C1, double C2)
{
    double d = T - Tref;
#if !(defined(MK_MUTATE) && MK_MUTATE == 404) /* planted fault: not held before the pole */
    if (d <= -0.99 * C2) d = -0.99 * C2;
#endif
    return C2 > 0 ? -C1 * d / (C2 + d) : NAN;
}

/* A and B through nu1 at T1 and nu2 at T2 (mm^2/s, K); -1 with A and B untouched for a point that has no log10 log10
   (nu not above 0.3, T not above 0 or not finite), equal temperatures or a NaN, else 0 */
static inline int mk_rate_walther_fit(double T1, double nu1, double T2, double nu2, double *A, double *B)
{
    const double z1 = log10(log10(nu1 + 0.7)), z2 = log10(log10(nu2 + 0.7)), l1 = log10(T1), l2 = log10(T2);
    const double b = (z1 - z2) / (l2 - l1), a = z1 + b * l1;
#if !(defined(MK_MUTATE) && MK_MUTATE == 405) /* planted fault: no refusal */
#if defined(MK_MUTATE) && MK_MUTATE == 407 /* planted fault: T2 of 0 or INFINITY accepted */
    if (!isfinite(a) || !isfinite(b)) return -1;
#else
    if (!isfinite(a) || !isfinite(b) || !isfinite(l2)) return -1;
#endif
#endif
    *A = a, *B = b;
    return 0;
}
/* the viscosity at T, mm^2/s; INFINITY where it overflows */
static inline double mk_rate_walther(double T, double A, double B)
{
#if defined(MK_MUTATE) && MK_MUTATE == 406
    return T > 0 ? pow(10, pow(10, A - B * log10(T))) : NAN; /* planted fault: the 0.7 not taken off */
#else
    return T > 0 ? pow(10, pow(10, A - B * log10(T))) - 0.7 : NAN;
#endif
}
#endif
