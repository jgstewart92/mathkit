/* mathkit: how the kit threads. A caller that builds with -fopenmp gets OpenMP threads in the kit's hot loops; one that
   does not gets no threads at all (the pragmas compile out) and the same bits. Work is cut into blocks whose length
   depends on the item count only, never on the thread count. Each block is a serial function; a sum gives each block a
   partial and adds the partials in block order. So a result is the same at any thread count and without OpenMP, and the
   block function is the serial reference form of the loop. Inside a caller's own parallel region the kit's loops run on
   that one thread (OpenMP's default of one active level). The underscore-suffixed macros are internal to the kit.
   Header only, C23. */
#ifndef MK_PAR_H
#define MK_PAR_H
#ifdef _OPENMP
#define MK_OMP_(...) _Pragma(#__VA_ARGS__)
#else
#define MK_OMP_(...)
#endif
#if defined(MK_MUTATE) && defined(_OPENMP)
#include <omp.h>
#define MK_PAR_NTHREADS_() omp_get_max_threads() /* planted faults only: a result made to depend on the thread count */
#else
#define MK_PAR_NTHREADS_() 1
#endif
#ifdef MK_MUTATE
/* planted faults only: n items cut into as many pieces as there are threads, each piece's last item left out */
#define MK_PAR_CUT_(n, i0, i1, stmt)                                                                                    \
    do {                                                                                                                 \
        const long mk_n_ = (n), mk_t_ = MK_PAR_NTHREADS_(), mk_c_ = (mk_n_ + mk_t_ - 1) / mk_t_;                         \
        for (long i0 = 0, i1; i0 < mk_n_; i0 += mk_c_) {                                                                \
            i1 = i0 + mk_c_ - 1 < mk_n_ ? i0 + mk_c_ - 1 : mk_n_;                                                         \
            stmt;                                                                                                        \
        }                                                                                                                \
    } while (0)
#endif

enum { MK_PAR_BLK = 4096, MK_PAR_NBLK = 256 };
/* the block length for n items: MK_PAR_BLK, or longer so there are at most MK_PAR_NBLK blocks; a function of n alone.
   Ceilings are n / k + (n % k != 0), which cannot overflow within k of LONG_MAX as (n + k - 1) / k would. */
static inline long mk_par_blk(long n)
{
#if defined(MK_MUTATE) && MK_MUTATE == 211
    const long b = (n + MK_PAR_NBLK - 1) / MK_PAR_NBLK; /* planted fault: ceilings that overflow near LONG_MAX */
#else
    const long b = n / MK_PAR_NBLK + (n % MK_PAR_NBLK != 0); /* n <= 0 gives MK_PAR_BLK */
#endif
    return b > MK_PAR_BLK ? b : MK_PAR_BLK;
}
/* the number of such blocks (0 for n <= 0) */
static inline long mk_par_nblk(long n)
{
    const long b = mk_par_blk(n);
#if defined(MK_MUTATE) && MK_MUTATE == 211
    return n > 0 ? (n + b - 1) / b : 0;
#else
    return n > 0 ? n / b + (n % b != 0) : 0;
#endif
}
/* block k's items [*i0, *i1), for 0 <= k < mk_par_nblk(n) */
static inline void mk_par_range(long n, long k, long *i0, long *i1)
{
    const long b = mk_par_blk(n);
#if defined(MK_MUTATE) && MK_MUTATE == 211
    *i0 = k * b, *i1 = (k + 1) * b < n ? (k + 1) * b : n;
#else
    /* k b < n: no overflow, where (k + 1) b can pass LONG_MAX on the last block */
    *i0 = k * b, *i1 = b < n - *i0 ? *i0 + b : n;
#endif
}
/* stmt (a call of a block function on [i0, i1)) for every block of n items, the blocks spread over the threads */
#define MK_PAR_FOR_(n, i0, i1, stmt)                                                                                    \
    do {                                                                                                                 \
        const long mk_n_ = (n), mk_nb_ = mk_par_nblk(mk_n_);                                                             \
        MK_OMP_(omp parallel for schedule(static) if(mk_nb_ > 1))                                                        \
        for (long mk_k_ = 0; mk_k_ < mk_nb_; mk_k_++) {                                                                  \
            long i0, i1;                                                                                                 \
            mk_par_range(mk_n_, mk_k_, &i0, &i1);                                                                        \
            stmt;                                                                                                        \
        }                                                                                                                \
    } while (0)
/* out = the sum of part (a block function's partial on [i0, i1)) over the blocks of n items, added in block order; one
   block gives the block function's own sum */
#if defined(MK_MUTATE) && MK_MUTATE == 107
/* planted fault: the partials added by OpenMP's reduction(+:), whose order of adding follows the thread count */
#define MK_PAR_SUM_(n, i0, i1, part, out)                                                                                \
    do {                                                                                                                 \
        const long mk_n_ = (n), mk_nb_ = mk_par_nblk(mk_n_);                                                             \
        double mk_s_ = 0;                                                                                                \
        MK_OMP_(omp parallel for schedule(static) reduction(+:mk_s_))                                                    \
        for (long mk_k_ = 0; mk_k_ < mk_nb_; mk_k_++) {                                                                  \
            long i0, i1;                                                                                                 \
            mk_par_range(mk_n_, mk_k_, &i0, &i1);                                                                        \
            mk_s_ += (part);                                                                                             \
        }                                                                                                                \
        (out) = mk_s_;                                                                                                   \
    } while (0)
#else
#define MK_PAR_SUM_(n, i0, i1, part, out)                                                                                \
    do {                                                                                                                 \
        const long mk_n_ = (n), mk_nb_ = mk_par_nblk(mk_n_);                                                             \
        double mk_p_[MK_PAR_NBLK], mk_s_ = 0;                                                                            \
        MK_OMP_(omp parallel for schedule(static) if(mk_nb_ > 1))                                                        \
        for (long mk_k_ = 0; mk_k_ < mk_nb_; mk_k_++) {                                                                  \
            long i0, i1;                                                                                                 \
            mk_par_range(mk_n_, mk_k_, &i0, &i1);                                                                        \
            mk_p_[mk_k_] = (part);                                                                                       \
        }                                                                                                                \
        for (long mk_k_ = 0; mk_k_ < mk_nb_; mk_k_++) mk_s_ += mk_p_[mk_k_];                                             \
        (out) = mk_s_;                                                                                                   \
    } while (0)
#endif
#endif
