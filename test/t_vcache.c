/* mk_vcache.h: mk_vcache_order on a shuffled regular grid and a shuffled sphere (each result's ACMR at caches 16, 32 and
   64 against the shuffled input's and under a stated bound; the result a permutation of the input's triangles, each with
   its corners in their order; the same bytes from two runs, and the bytes recorded for this version), equal scores going
   to the lowest index, triangles listed more than once and triangles that repeat a vertex, and the refusals with the list
   untouched; mk_vcache_acmr on orders counted by hand (a FIFO, not an LRU, of exactly k) and its refusals. */
#include "../include/mk_vcache.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "t.h"

static uint64_t next64(void)
{
    t_state ^= t_state << 13, t_state ^= t_state >> 7, t_state ^= t_state << 17;
    return t_state;
}
/* the triangles' order shuffled and each triangle's corners turned a random number of places (its winding kept) */
static void shuffle(uint32_t *T, size_t nt)
{
    for (size_t i = nt - 1; i > 0; i--) {
        const size_t j = (size_t)(next64() % (i + 1));
        uint32_t x[3];
        memcpy(x, T + 3 * i, sizeof x), memcpy(T + 3 * i, T + 3 * j, sizeof x), memcpy(T + 3 * j, x, sizeof x);
    }
    for (size_t i = 0; i < nt; i++) {
        const int r = (int)(next64() % 3);
        const uint32_t a = T[3 * i], b = T[3 * i + 1], c = T[3 * i + 2];
        if (r == 1) T[3 * i] = b, T[3 * i + 1] = c, T[3 * i + 2] = a;
        if (r == 2) T[3 * i] = c, T[3 * i + 1] = a, T[3 * i + 2] = b;
    }
}
/* a g x g grid of vertices, two triangles a cell */
static uint32_t *grid(uint32_t g, size_t *nt, uint32_t *nv)
{
    *nt = 2 * (size_t)(g - 1) * (g - 1), *nv = g * g;
    uint32_t *T = malloc(sizeof *T * 3 * *nt), n = 0;
    for (uint32_t i = 0; i + 1 < g; i++)
        for (uint32_t j = 0; j + 1 < g; j++) {
            const uint32_t a = i * g + j, b = a + 1, c = a + g, d = c + 1;
            T[n++] = a, T[n++] = c, T[n++] = b, T[n++] = b, T[n++] = c, T[n++] = d;
        }
    return T;
}
/* a closed sphere by latitude and longitude: two poles, rings rows of nl vertices each; each pole a fan of nl triangles */
static uint32_t *sphere(uint32_t rings, uint32_t nl, size_t *nt, uint32_t *nv)
{
    *nv = 2 + rings * nl, *nt = 2 * (size_t)nl * rings;
    uint32_t *T = malloc(sizeof *T * 3 * *nt), n = 0;
    const uint32_t S = *nv - 1;
    for (uint32_t j = 0; j < nl; j++) {
        const uint32_t j1 = (j + 1) % nl;
        T[n++] = 0, T[n++] = 1 + j, T[n++] = 1 + j1;
        T[n++] = S, T[n++] = 1 + (rings - 1) * nl + j1, T[n++] = 1 + (rings - 1) * nl + j;
    }
    for (uint32_t r = 0; r + 1 < rings; r++)
        for (uint32_t j = 0; j < nl; j++) {
            const uint32_t j1 = (j + 1) % nl, a = 1 + r * nl + j, b = 1 + r * nl + j1, c = a + nl, d = b + nl;
            T[n++] = a, T[n++] = c, T[n++] = b, T[n++] = b, T[n++] = c, T[n++] = d;
        }
    return T;
}
static int cmp_tri(const void *x, const void *y)
{
    const uint32_t *a = x, *b = y;
    for (int k = 0; k < 3; k++)
        if (a[k] != b[k]) return a[k] < b[k] ? -1 : 1;
    return 0;
}
/* the same triangles as ordered triples, each as many times */
static int same_triangles(const uint32_t *A, const uint32_t *B, size_t nt)
{
    uint32_t *a = malloc(sizeof *a * 3 * nt), *b = malloc(sizeof *b * 3 * nt);
    memcpy(a, A, sizeof *a * 3 * nt), memcpy(b, B, sizeof *b * 3 * nt);
    qsort(a, nt, sizeof *a * 3, cmp_tri), qsort(b, nt, sizeof *b * 3, cmp_tri);
    const int r = memcmp(a, b, sizeof *a * 3 * nt) == 0;
    free(a), free(b);
    return r;
}
/* FNV-1a over the indices as 4 little-endian bytes each */
static uint64_t hash(const uint32_t *T, size_t n)
{
    uint64_t h = 0xcbf29ce484222325u;
    for (size_t i = 0; i < n; i++)
        for (int k = 0; k < 4; k++) h = (h ^ ((T[i] >> (8 * k)) & 0xff)) * 0x100000001b3u;
    return h;
}

/* one mesh: reorder a copy twice, check the rows */
static void mesh_rows(const char *name, const uint32_t *T0, size_t nt, uint32_t nv, const double bound[3], uint64_t want)
{
    uint32_t *A = malloc(sizeof *A * 3 * nt), *B = malloc(sizeof *B * 3 * nt);
    memcpy(A, T0, sizeof *A * 3 * nt), memcpy(B, T0, sizeof *B * 3 * nt);
    const int ra = mk_vcache_order(A, nt, nv), rb = mk_vcache_order(B, nt, nv);
    static const int ks[3] = {16, 32, 64};
    double in[3], out[3];
    for (int i = 0; i < 3; i++) in[i] = mk_vcache_acmr(T0, nt, nv, ks[i]), out[i] = mk_vcache_acmr(A, nt, nv, ks[i]);
    CHECK(ra == 0 && rb == 0 && same_triangles(A, T0, nt), "%s: %zu triangles, %u vertices: the result a permutation of the input's triangles, each with its corners in order (rc %d)",
          name, nt, nv, ra);
    for (int i = 0; i < 3; i++)
        CHECK(out[i] <= bound[i] && out[i] < in[i], "%s: ACMR at cache %d %.4f, under %.2f (shuffled input %.4f)", name, ks[i], out[i], bound[i], in[i]);
    CHECK(memcmp(A, B, sizeof *A * 3 * nt) == 0, "%s: two runs give the same bytes", name);
    const uint64_t h = hash(A, 3 * nt);
    CHECK(h == want, "%s: the order is the one recorded for this version (hash %016llx, recorded %016llx)", name, (unsigned long long)h, (unsigned long long)want);
    free(A), free(B);
}

int main(void)
{
    { /* the scores as documented: each term rounded once from its formula, n over 32 counted as 32 */
        long off = 0;
        for (uint32_t p = 0; p <= MK_VCACHE_SIZE; p++)
            for (uint32_t n = 1; n <= 100; n++) {
                const double place = p < 3 ? 0.75 : p < MK_VCACHE_SIZE ? pow(1.0 - (double)(p - 3) / (MK_VCACHE_SIZE - 3), 1.5) : 0.0;
                const long want = lround(65536.0 * place) + lround(65536.0 * 2.0 / sqrt((double)(n < MK_VCACHE_VALENCE ? n : MK_VCACHE_VALENCE)));
                off += mk_vcache_score_(p, n) != want;
            }
        CHECK(off == 0, "a vertex's score is 65536 (0.75, or (1 - (p - 3) / 29)^1.5 at place p, 0 outside) + 65536 x 2 / sqrt(n), each rounded: %ld of %d differ", off,
              (MK_VCACHE_SIZE + 1) * 100);
    }
    /* The bounds are about 4 % above this version's results (grid 0.7254, 0.6977, 0.6696; sphere 0.7161, 0.6921, 0.6691);
       about 0.5 is the floor for a closed mesh, and an open grid's rim raises it */
    const double bound[3] = {0.75, 0.72, 0.70};
    { /* a shuffled 160 x 160 grid */
        size_t nt;
        uint32_t nv, *T = grid(160, &nt, &nv);
        shuffle(T, nt);
        mesh_rows("grid", T, nt, nv, bound, 0x161a696a9deb5a5bu);
        free(T);
    }
    { /* a shuffled closed sphere, 120 rings of 96: each pole is used by 96 triangles, past the score's table */
        size_t nt;
        uint32_t nv, *T = sphere(120, 96, &nt, &nv);
        shuffle(T, nt);
        mesh_rows("sphere", T, nt, nv, bound, 0x22698d363c71aa9eu);
        free(T);
    }
    { /* 1000 triangles with no vertex in common all score alike: the order stays the input's */
        enum { N = 1000 };
        uint32_t T[3 * N], U[3 * N];
        for (uint32_t i = 0; i < 3 * N; i++) T[i] = i;
        for (uint32_t i = 3 * N - 1; i > 0; i--) { /* the vertices numbered at random: each still in one triangle */
            const uint32_t j = (uint32_t)(next64() % (i + 1)), x = T[i];
            T[i] = T[j], T[j] = x;
        }
        memcpy(U, T, sizeof T);
        const int r = mk_vcache_order(U, N, 3 * N);
        CHECK(r == 0 && memcmp(U, T, sizeof T) == 0, "1000 triangles sharing no vertex, all scoring alike, keep their order (equal scores to the lowest index; rc %d)", r);
    }
    { /* an 8 x 8 grid with each triangle listed three times, and triangles that repeat a vertex, shuffled */
        size_t nt;
        uint32_t nv, *G = grid(9, &nt, &nv);
        const size_t n = 3 * nt + 4;
        uint32_t *T = malloc(sizeof *T * 3 * n), *U = malloc(sizeof *U * 3 * n);
        for (int c = 0; c < 3; c++) memcpy(T + 3 * nt * (size_t)c, G, sizeof *G * 3 * nt);
        const uint32_t extra[12] = {5, 5, 6, 7, 7, 7, 10, 11, 10, 40, 40, 40};
        memcpy(T + 9 * nt, extra, sizeof extra);
        shuffle(T, n);
        memcpy(U, T, sizeof *T * 3 * n);
        const int r = mk_vcache_order(U, n, nv);
        CHECK(r == 0 && same_triangles(U, T, n), "%zu triangles, each listed three times, and 4 that repeat a vertex: a permutation, corners in order (rc %d, ACMR at 32 %.4f from %.4f)", n,
              r, mk_vcache_acmr(U, n, nv, 32), mk_vcache_acmr(T, n, nv, 32));
        free(G), free(T), free(U);
    }
    { /* refusals, each with the list untouched; no triangles is no work */
        uint32_t T[6] = {0, 1, 2, 2, 1, 3}, U[6];
        memcpy(U, T, sizeof T);
        const int a = mk_vcache_order(U, 2, 3); /* index 3 with 3 vertices */
        const int a_same = memcmp(U, T, sizeof T) == 0;
        const int b = mk_vcache_order(U, 2, 0), c = mk_vcache_order(NULL, 1, 4), d = mk_vcache_order(U, (size_t)MK_VCACHE_MAXT + 1, 4);
        const int b_same = memcmp(U, T, sizeof T) == 0;
        const int e = mk_vcache_order(NULL, 0, 0), f = mk_vcache_order(U, 0, 4);
        CHECK(a == -1 && a_same && b == -1 && c == -1 && d == -1 && b_same && e == 0 && f == 0 && memcmp(U, T, sizeof T) == 0,
              "refused: an index not below nv (%d), nv 0 (%d), NULL (%d), nt over MK_VCACHE_MAXT (%d), the list untouched; nt 0 gives 0 (%d, %d)", a, b, c, d, e, f);
    }
    { /* mk_vcache_acmr counted by hand */
        const uint32_t A[6] = {0, 1, 2, 3, 0, 1}, B[9] = {0, 1, 2, 0, 3, 4, 0, 5, 6};
        const double a3 = mk_vcache_acmr(A, 2, 4, 3), a4 = mk_vcache_acmr(A, 2, 4, 4), b3 = mk_vcache_acmr(B, 3, 7, 3);
        CHECK(a3 == 3.0 && a4 == 2.0, "ACMR of 0 1 2, 3 0 1: 6 misses at cache 3 (%.4f a triangle), 4 at cache 4 (%.4f)", a3, a4);
        CHECK(fabs(b3 - 8.0 / 3) < 1e-15, "ACMR of 0 1 2, 0 3 4, 0 5 6 at cache 3: 8 misses, a FIFO (an LRU would keep 0: 7) (%.4f a triangle)", b3);
        const double r1 = mk_vcache_acmr(A, 2, 3, 3), r2 = mk_vcache_acmr(A, 2, 4, 0), r3 = mk_vcache_acmr(NULL, 1, 4, 3), r4 = mk_vcache_acmr(NULL, 0, 4, 3);
        CHECK(r1 == -1 && r2 == -1 && r3 == -1 && r4 == 0, "ACMR refuses an index not below nv (%g), cache 0 (%g), NULL (%g); nt 0 gives 0 (%g)", r1, r2, r3, r4);
    }
    return wrong;
}
