/* mk_minmax.h: fmax, fmin, fmaxf and fminf inline, with glibc's results bit for bit. The C library's are calls the
   compiler doesn't inline, which costs real time in code that clamps a lot.

   Ordered arguments are decided here. Unequal: the larger (fmax) or the smaller (fmin). Equal, +0 and -0 included:
   fmax gives y when x's sign bit is set and x otherwise, fmin gives x when x's sign bit is set and y otherwise, so -0
   counts below +0 and the result is glibc's. A pair with a NaN goes to the C library's own function, through a pointer
   the compiler cannot see into, so a NaN result has the library's bits (asked as builtins, gcc folds NaN cases by its own
   rule, which differs from glibc's).

   mk_fmax, mk_fmin, mk_fmaxf and mk_fminf are the calls. Defined before this header (or passed as -DMK_MINMAX_LIBM with
   -include mk_minmax.h), MK_MINMAX_LIBM also makes fmax, fmin, fmaxf and fminf in the including unit mean these. The header
   doesn't include <math.h>, so it can be -included ahead of a unit's own feature macros; <math.h> included later is
   fine. With another C library a NaN still goes to that library; its ordered results may differ from these only where it
   orders signed zeros differently. */
#ifndef MK_MINMAX_H
#define MK_MINMAX_H

double fmax(double, double);
double fmin(double, double);
float fmaxf(float, float);
float fminf(float, float);
static double (*volatile const mk_fmax_lib_)(double, double) = fmax;
static double (*volatile const mk_fmin_lib_)(double, double) = fmin;
static float (*volatile const mk_fmaxf_lib_)(float, float) = fmaxf;
static float (*volatile const mk_fminf_lib_)(float, float) = fminf;

static inline double mk_fmax(double x, double y)
{
    if (__builtin_isunordered(x, y)) return mk_fmax_lib_(x, y);
#if defined(MK_MUTATE) && MK_MUTATE == 380
    if (x == y) return x; /* planted fault: equal arguments give x, so fmax(-0, +0) is -0 */
#else
    if (x == y) return __builtin_signbit(x) ? y : x;
#endif
    return x > y ? x : y;
}
static inline double mk_fmin(double x, double y)
{
    if (__builtin_isunordered(x, y)) return mk_fmin_lib_(x, y);
    if (x == y) return __builtin_signbit(x) ? x : y;
    return x < y ? x : y;
}
static inline float mk_fmaxf(float x, float y)
{
    if (__builtin_isunordered(x, y)) return mk_fmaxf_lib_(x, y);
    if (x == y) return __builtin_signbit(x) ? y : x;
    return x > y ? x : y;
}
static inline float mk_fminf(float x, float y)
{
#if defined(MK_MUTATE) && MK_MUTATE == 381
    if (__builtin_isunordered(x, y)) return x < y ? x : y; /* planted fault: a NaN decided here, not by the library */
#else
    if (__builtin_isunordered(x, y)) return mk_fminf_lib_(x, y);
#endif
    if (x == y) return __builtin_signbit(x) ? x : y;
    return x < y ? x : y;
}

#ifdef MK_MINMAX_LIBM
#define fmax(x, y) mk_fmax(x, y)
#define fmin(x, y) mk_fmin(x, y)
#define fmaxf(x, y) mk_fmaxf(x, y)
#define fminf(x, y) mk_fminf(x, y)
#endif
#endif
