/* mathkit: random numbers and compensated sums. mk_rng is xoshiro256** as published (Blackman and Vigna 2021), its
   state seeded by splitmix64 from a seed and a stream number, so separate streams are independent and reproducible.
   Uniform doubles in (0,1] (safe for log) and [0,1), unbiased integers below n, Gaussian by Box-Muller with the second
   value kept. mk_sum is Neumaier's compensated sum. Header only, C23. */
#ifndef MK_RNG_H
#define MK_RNG_H
#include <math.h>
#include <stdint.h>

#define MK_RNG_PI 3.14159265358979323846

/* splitmix64: a step of a 64-bit counter through a mixer; also a good hash of an integer */
static inline uint64_t mk_splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline uint64_t mk_mix64(uint64_t x) { return mk_splitmix64(&x); }

typedef struct {
    uint64_t s[4];
    double spare;   /* Box-Muller's second value */
    int has_spare;
} mk_rng;

static inline void mk_rng_seed(mk_rng *g, uint64_t seed, uint64_t stream)
{
    uint64_t x = seed * 0x2545f4914f6cdd1dULL + stream * 0x9e3779b97f4a7c15ULL + 1;
    for (int i = 0; i < 4; i++) g->s[i] = mk_splitmix64(&x);
    if (!(g->s[0] | g->s[1] | g->s[2] | g->s[3])) g->s[0] = 1; /* the one state xoshiro never leaves */
    g->has_spare = 0, g->spare = 0;
}
static inline uint64_t mk_rotl64(uint64_t v, int k) { return (v << (k & 63)) | (v >> (-k & 63)); }
static inline uint64_t mk_rng_next(mk_rng *g)
{
    uint64_t *s = g->s, r = mk_rotl64(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0], s[3] ^= s[1], s[1] ^= s[2], s[0] ^= s[3], s[2] ^= t;
#if defined(MK_MUTATE) && MK_MUTATE == 31
    s[3] = mk_rotl64(s[3], 44); /* planted fault: the wrong rotation */
#else
    s[3] = mk_rotl64(s[3], 45);
#endif
    return r;
}
static inline double mk_rng_u01(mk_rng *g) { return ((double)(mk_rng_next(g) >> 11) + 1.0) * 0x1.0p-53; } /* (0,1] */
static inline double mk_rng_u53(mk_rng *g) { return (double)(mk_rng_next(g) >> 11) * 0x1.0p-53; }         /* [0,1) */
/* uniform integer in [0, n), n >= 1, without modulo bias */
static inline uint64_t mk_rng_below(mk_rng *g, uint64_t n)
{
    if (!n) return 0; /* n = 0 is outside the contract: 0, not a division by zero */
    uint64_t lim = -n % n; /* 2^64 mod n */
    for (;;) {
        uint64_t r = mk_rng_next(g);
        if (r >= lim) return r % n;
    }
}
static inline double mk_rng_gauss(mk_rng *g)
{
    if (g->has_spare) { g->has_spare = 0; return g->spare; }
    double u = mk_rng_u01(g), v = mk_rng_u53(g);
#if defined(MK_MUTATE) && MK_MUTATE == 32
    double r = sqrt(-log(u)); /* planted fault: the radius without its factor 2 */
#else
    double r = sqrt(-2 * log(u));
#endif
    g->spare = r * sin(2 * MK_RNG_PI * v), g->has_spare = 1;
    return r * cos(2 * MK_RNG_PI * v);
}

/* Neumaier's sum (Neumaier 1974): Kahan's, also right when an addend is larger than the running total. Once the sum is
   not finite the compensation stops, so the value is what a plain sum gives: inf, -inf, or NaN (inf + -inf, or a NaN
   addend). */
typedef struct { double s, c; } mk_sum;
static inline void mk_sum_add(mk_sum *a, double x)
{
    double t = a->s + x;
#if !(defined(MK_MUTATE) && MK_MUTATE == 33) /* planted fault 33: no compensation */
#if defined(MK_MUTATE) && MK_MUTATE == 34
    if (1) { /* planted fault: compensation past an infinite sum, giving NaN */
#else
    if (isfinite(t)) { /* else (s - t) + x is inf - inf: NaN */
#endif
        if (fabs(a->s) >= fabs(x)) a->c += (a->s - t) + x;
        else a->c += (x - t) + a->s;
    }
#endif
    a->s = t;
}
static inline double mk_sum_value(const mk_sum *a) { return a->s + a->c; }
#endif
