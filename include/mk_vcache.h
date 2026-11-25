/* mathkit: a triangle order for a GPU's post-transform vertex cache, and the cache misses of an order.

   mk_vcache_order reorders an indexed triangle list in place so that each triangle's vertices were, as far as it can
   manage, used by the triangles just before it: a GPU that keeps its latest vertex shader results then runs the shader
   fewer times. It is Forsyth's method (2006) with integer scores. A simulated cache holds the MK_VCACHE_SIZE (32)
   vertices used most recently, newest first. A vertex scores for its place p in that cache, 0.75 in the three newest
   places and (1 - (p - 3) / 29)^1.5 from place 3 to 31, 0 outside it, plus 2 / sqrt(n) for the n triangles still to come
   that use it (n over MK_VCACHE_VALENCE, 32, counts as 32); a triangle scores the sum over its three corners. The first
   triangle is the best of all; each next one is the best of the triangles still to come that use a vertex in the cache
   or one the last triangle pushed out of it; when there is none, the lowest-numbered triangle still to come. Equal scores
   go to the lowest-numbered triangle. The two terms are constant tables in 16.16 fixed point, rounded once from those
   formulas, so no floating point enters the order and the same input gives the same bytes on every machine and compiler.

   Each triangle keeps its three indices in their order, so its winding is unchanged, and the result is a permutation of
   the input's triangles. A triangle that repeats a vertex, and a triangle listed more than once, are kept like any other.
   Cost: each triangle placed rescores the triangles still to come around up to 35 vertices, so the time is linear in the
   triangle count where few triangles share a vertex (a body mesh of 1.77 M triangles: about 35 scores a triangle placed,
   0.55 s on one core of a desktop CPU), but a vertex with n corners in the list (a triangle naming it twice has two)
   costs O(n) on every step it stays in the cache, so a fan of n triangles round one vertex takes O(n^2). Scratch: about
   25 bytes a triangle and 12 a vertex.

   mk_vcache_acmr measures an order: the average cache misses per triangle (ACMR) under a FIFO cache of k vertices, the
   usual model of the hardware's cache. Each index is a lookup; a miss enters the vertex and, once k are held, drops the
   oldest. 3 is the worst; a large closed mesh, with about half as many vertices as triangles, cannot go below about 0.5.
   Scratch: 8 bytes a vertex. Header only, C23. */
#ifndef MK_VCACHE_H
#define MK_VCACHE_H
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MK_VCACHE_SIZE 32
#define MK_VCACHE_VALENCE 32
#define MK_VCACHE_MAXT 0x55555555u /* the most triangles: 3 nt must fit in 32 bits */

/* 65536 (0.75 for places 0 to 2, (1 - (p - 3) / 29)^1.5 for p from 3 to 31), rounded; place 32 is outside the cache */
static const int32_t mk_vcache_place_[MK_VCACHE_SIZE + 1] = {
    49152, 49152, 49152, 65536, 62176, 58875, 55634, 52456, 49340, 46289, 43303, 40384, 37534, 34755, 32047, 29414, 26857,
    24379, 21982, 19670, 17444, 15310, 13270, 11330, 9496, 7772, 6168, 4692, 3357, 2181, 1187, 420, 0};
/* 65536 x 2 / sqrt(n) for n triangles still to come, rounded (n 0 never scores) */
static const int32_t mk_vcache_left_[MK_VCACHE_VALENCE + 1] = {
    0, 131072, 92682, 75674, 65536, 58617, 53510, 49541, 46341, 43691, 41449, 39520, 37837, 36353, 35030, 33843, 32768,
    31790, 30894, 30070, 29309, 28602, 27945, 27330, 26755, 26214, 25705, 25225, 24770, 24339, 23930, 23541, 23170};

/* a vertex's score at place p (MK_VCACHE_SIZE: outside) with n triangles still to come */
static inline int32_t mk_vcache_score_(uint32_t p, uint32_t n)
{
#if defined(MK_MUTATE) && MK_MUTATE == 390
    (void)n;
    return mk_vcache_place_[p]; /* planted fault: the term for the triangles still to come dropped */
#elif defined(MK_MUTATE) && MK_MUTATE == 391
    return mk_vcache_place_[p < MK_VCACHE_SIZE ? p + 1 : p] + mk_vcache_left_[n < MK_VCACHE_VALENCE ? n : MK_VCACHE_VALENCE]; /* planted fault: the place one off */
#else
    return mk_vcache_place_[p] + mk_vcache_left_[n < MK_VCACHE_VALENCE ? n : MK_VCACHE_VALENCE];
#endif
}
/* a candidate as one key, its score above and its index's complement below, so the largest key is the best score and,
   of equal scores, the lowest index; 0 is no candidate (an index's complement is never 0 below MK_VCACHE_MAXT) */
static inline uint64_t mk_vcache_key_(int32_t s, uint32_t t)
{
#if defined(MK_MUTATE) && MK_MUTATE == 392
    return (uint64_t)(uint32_t)s << 32 | (t + 1); /* planted fault: equal scores go to the highest index */
#else
    return (uint64_t)(uint32_t)s << 32 | (uint32_t)~t;
#endif
}

/* Reorders the nt triangles of idx (3 nt vertex indices, each below nv) in place for the vertex cache. Returns 0, or -1
   for a bad argument (idx NULL with nt > 0, nt over MK_VCACHE_MAXT, an index not below nv) or no memory, with idx
   untouched. nt 0 returns 0. */
static inline int mk_vcache_order(uint32_t *idx, size_t nt, uint32_t nv)
{
    if (nt == 0) return 0;
    if (!idx || nt > MK_VCACHE_MAXT) return -1;
    const size_t ni = 3 * nt;
    for (size_t i = 0; i < ni; i++)
        if (idx[i] >= nv) return -1;
    uint32_t *start = malloc(sizeof *start * ((size_t)nv + 1)), *left = calloc(nv, sizeof *left), *adj = malloc(sizeof *adj * ni);
    uint32_t *out = malloc(sizeof *out * ni);
    int32_t *vs = malloc(sizeof *vs * nv);
    unsigned char *done = calloc(nt, 1);
    if (!start || !left || !adj || !out || !vs || !done) {
        free(start), free(left), free(adj), free(out), free(vs), free(done);
        return -1;
    }
    /* each vertex's triangles, in index order; left[v] how many are still to come (their list's live length) */
    for (size_t i = 0; i < ni; i++) left[idx[i]]++;
    start[0] = 0;
    for (uint32_t v = 0; v < nv; v++) start[v + 1] = start[v] + left[v], left[v] = 0;
    for (size_t t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) {
            const uint32_t v = idx[3 * t + k];
            adj[start[v] + left[v]++] = (uint32_t)t;
        }
    for (uint32_t v = 0; v < nv; v++) vs[v] = mk_vcache_score_(MK_VCACHE_SIZE, left[v]);
    uint64_t key = 0;
    for (size_t t = 0; t < nt; t++) {
        const uint64_t k = mk_vcache_key_(vs[idx[3 * t]] + vs[idx[3 * t + 1]] + vs[idx[3 * t + 2]], (uint32_t)t);
        key = k > key ? k : key;
    }
#if defined(MK_MUTATE) && MK_MUTATE == 392
    uint32_t best = (uint32_t)key - 1, cursor = 0;
#else
    uint32_t best = ~(uint32_t)key, cursor = 0;
#endif
    uint32_t cache[MK_VCACHE_SIZE + 3], nn[MK_VCACHE_SIZE + 3], nc = 0;
    for (size_t e = 0; e < nt; e++) {
        const uint32_t t = best, *tv = idx + 3 * (size_t)t;
        done[t] = 1;
        memcpy(out + 3 * e, tv, sizeof *out * 3);
        /* t off its vertices' lists; its vertices, each once, to the front of the cache */
        uint32_t m = 0;
        for (int k = 0; k < 3; k++) {
            const uint32_t v = tv[k], n = left[v];
            uint32_t *a = adj + start[v];
            for (uint32_t i = 0; i < n; i++)
                if (a[i] == t) {
                    a[i] = a[n - 1];
                    break;
                }
            left[v] = n - 1;
            if (!(m > 0 && nn[0] == v) && !(m > 1 && nn[1] == v)) nn[m++] = v;
        }
        const uint32_t m0 = m;
        for (uint32_t i = 0; i < nc; i++) {
            const uint32_t v = cache[i];
            if (v != nn[0] && (m0 < 2 || v != nn[1]) && (m0 < 3 || v != nn[2])) nn[m++] = v;
        }
        nc = m < MK_VCACHE_SIZE ? m : MK_VCACHE_SIZE;
        /* the new scores of the cache's vertices and of those pushed out of it */
        for (uint32_t i = 0; i < m; i++) {
            if (i < nc) cache[i] = nn[i];
            vs[nn[i]] = mk_vcache_score_(i < MK_VCACHE_SIZE ? i : MK_VCACHE_SIZE, left[nn[i]]);
        }
        /* the best of their triangles still to come */
        key = 0;
        for (uint32_t i = 0; i < m; i++) {
            const uint32_t *a = adj + start[nn[i]], n = left[nn[i]];
            for (uint32_t j = 0; j < n; j++) {
                const uint32_t *w = idx + 3 * (size_t)a[j];
                const uint64_t k = mk_vcache_key_(vs[w[0]] + vs[w[1]] + vs[w[2]], a[j]);
                key = k > key ? k : key;
            }
        }
#if defined(MK_MUTATE) && MK_MUTATE == 392
        if (key) best = (uint32_t)key - 1;
#else
        if (key) best = ~(uint32_t)key;
#endif
        else if (e + 1 < nt) { /* a dead end: the lowest-numbered triangle still to come */
            while (done[cursor]) cursor++;
            best = cursor;
        }
    }
    memcpy(idx, out, sizeof *out * ni);
    free(start), free(left), free(adj), free(out), free(vs), free(done);
    return 0;
}

/* The average cache misses per triangle of the order idx (nt triangles, indices below nv) under a FIFO cache of k
   vertices. Returns -1 for a bad argument (idx NULL with nt > 0, k below 1, an index not below nv) or no memory; 0 for
   nt 0. */
static inline double mk_vcache_acmr(const uint32_t *idx, size_t nt, uint32_t nv, int k)
{
    if (k < 1) return -1;
    if (nt == 0) return 0;
    if (!idx || nt > SIZE_MAX / 3) return -1;
    uint64_t *in = calloc(nv ? nv : 1, sizeof *in), miss = 0; /* the miss count at which each vertex entered; 0 never */
    if (!in) return -1;
    for (size_t i = 0; i < 3 * nt; i++) {
        const uint32_t v = idx[i];
        if (v >= nv) {
            free(in);
            return -1;
        }
#if defined(MK_MUTATE) && MK_MUTATE == 393
        if (in[v] && miss - in[v] <= (uint64_t)k) continue; /* planted fault: the cache one entry larger */
#else
        if (in[v] && miss - in[v] < (uint64_t)k) continue;
#endif
        in[v] = ++miss;
    }
    free(in);
    return (double)miss / (double)nt;
}
#endif
