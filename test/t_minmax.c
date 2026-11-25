/* mk_minmax.h: mk_fmax, mk_fmin, mk_fmaxf and mk_fminf bit for bit against the C library's own functions, on every pair
   of special values (signed zeros, subnormals, the extremes, infinities, NaNs with payloads and signs) and on 4 million
   random pairs (bit patterns, equal pairs, opposite pairs); and MK_MINMAX_LIBM turning fmax and the others in this unit
   into the inline calls, with <math.h> included after the header. */
#define MK_MINMAX_LIBM
#include "../include/mk_minmax.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "t.h"

static double (*volatile const ref_fmax)(double, double) = (fmax);
static double (*volatile const ref_fmin)(double, double) = (fmin);
static float (*volatile const ref_fmaxf)(float, float) = (fmaxf);
static float (*volatile const ref_fminf)(float, float) = (fminf);

static uint64_t bits64(double x) { uint64_t u; memcpy(&u, &x, 8); return u; }
static uint32_t bits32(float x) { uint32_t u; memcpy(&u, &x, 4); return u; }
static double from64(uint64_t u) { double x; memcpy(&x, &u, 8); return x; }
static float from32(uint32_t u) { float x; memcpy(&x, &u, 4); return x; }
static uint64_t next64(void)
{
    t_state ^= t_state << 13, t_state ^= t_state >> 7, t_state ^= t_state << 17;
    return t_state;
}

enum { NS = 22, NR = 1 << 22 };
static const uint64_t special64[NS] = {
    0x0000000000000000, 0x8000000000000000, 0x0000000000000001, 0x8000000000000001, 0x0010000000000000,
    0x8010000000000000, 0x3ff0000000000000, 0xbff0000000000000, 0x7fefffffffffffff, 0xffefffffffffffff,
    0x7ff0000000000000, 0xfff0000000000000, 0x7ff8000000000000, 0xfff8000000000000, 0x7ff8000000000123,
    0xfff8000000000456, 0x7ff4000000000001, 0x3ff0000000000001, 0x000fffffffffffff, 0x800fffffffffffff,
    0x4000000000000000, 0xc000000000000000};
static const uint32_t special32[NS] = {
    0x00000000, 0x80000000, 0x00000001, 0x80000001, 0x00800000, 0x80800000, 0x3f800000, 0xbf800000,
    0x7f7fffff, 0xff7fffff, 0x7f800000, 0xff800000, 0x7fc00000, 0xffc00000, 0x7fc00123, 0xffc00456,
    0x7fa00001, 0x3f800001, 0x007fffff, 0x807fffff, 0x40000000, 0xc0000000};

/* the pairs a test runs: every special pair, then random bit patterns, each also against itself and its negation */
static long check64(double (*f)(double, double), double (*volatile const ref)(double, double), long *nan_pairs)
{
    long diff = 0, nan = 0;
    for (int i = 0; i < NS; i++)
        for (int j = 0; j < NS; j++) {
            const double x = from64(special64[i]), y = from64(special64[j]);
            diff += bits64(f(x, y)) != bits64(ref(x, y)), nan += isnan(x) || isnan(y);
        }
    for (long k = 0; k < NR; k++) {
        const double x = from64(next64()), y = k % 3 == 0 ? x : k % 3 == 1 ? -x : from64(next64());
        diff += bits64(f(x, y)) != bits64(ref(x, y)), nan += isnan(x) || isnan(y);
    }
    *nan_pairs = nan;
    return diff;
}
static long check32(float (*f)(float, float), float (*volatile const ref)(float, float), long *nan_pairs)
{
    long diff = 0, nan = 0;
    for (int i = 0; i < NS; i++)
        for (int j = 0; j < NS; j++) {
            const float x = from32(special32[i]), y = from32(special32[j]);
            diff += bits32(f(x, y)) != bits32(ref(x, y)), nan += isnan(x) || isnan(y);
        }
    for (long k = 0; k < NR; k++) {
        const float x = from32((uint32_t)next64()), y = k % 3 == 0 ? x : k % 3 == 1 ? -x : from32((uint32_t)next64());
        diff += bits32(f(x, y)) != bits32(ref(x, y)), nan += isnan(x) || isnan(y);
    }
    *nan_pairs = nan;
    return diff;
}

int main(void)
{
    const long total = (long)NS * NS + NR;
    long nan, d;
    d = check64(mk_fmax, ref_fmax, &nan);
    CHECK(d == 0, "mk_fmax equals the C library's fmax bit for bit: %ld of %ld pairs differ (%ld with a NaN)", d, total, nan);
    d = check64(mk_fmin, ref_fmin, &nan);
    CHECK(d == 0, "mk_fmin equals fmin: %ld of %ld pairs differ (%ld with a NaN)", d, total, nan);
    d = check32(mk_fmaxf, ref_fmaxf, &nan);
    CHECK(d == 0, "mk_fmaxf equals fmaxf: %ld of %ld pairs differ (%ld with a NaN)", d, total, nan);
    d = check32(mk_fminf, ref_fminf, &nan);
    CHECK(d == 0, "mk_fminf equals fminf: %ld of %ld pairs differ (%ld with a NaN)", d, total, nan);
    const double pz = 0.0, nz = -0.0;
    CHECK(bits64(mk_fmax(nz, pz)) == bits64(pz) && bits64(mk_fmax(pz, nz)) == bits64(pz) &&
              bits64(mk_fmin(nz, pz)) == bits64(nz) && bits64(mk_fmin(pz, nz)) == bits64(nz) && mk_fmax(NAN, 2) == 2 &&
              mk_fmin(3, NAN) == 3 && isnan(mk_fmax(NAN, NAN)),
          "-0 counts below +0 either way round, and a NaN gives the other argument");
#ifdef fmax
    const int macro = 1;
#else
    const int macro = 0;
#endif
    volatile double a = -0.0, b = 0.0;
    CHECK(macro && bits64(fmax(a, b)) == bits64(ref_fmax(a, b)) && bits32(fminf((float)a, (float)b)) == bits32(ref_fminf((float)a, (float)b)),
          "with MK_MINMAX_LIBM, fmax and the others are the inline calls, <math.h> included after (macro %d)", macro);
    return wrong;
}
