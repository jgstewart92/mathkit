/* mk_decimate: a sphere and a torus decimated within their distance bound (each input vertex near the result and the
   result near the input, judged here by brute force), closed with the Euler characteristic unchanged; a hollow ball
   whose two walls lie closer than the edges grow (no crossing pair, judged by brute force); a bumpy sphere's fold
   lengths; held vertices; a target; an open cap's rim; refusals. With a folder argument the results are written there
   for test/check.sh's meshcheck. */
/* every block the header allocates carries guard bytes on both sides, checked when it is freed (test-only): a write
   outside a block is a FAIL row here, where it would otherwise corrupt the heap and abort without a row (planted fault
   201, which leaves a needle's tip with no owner and writes before the owner table) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
enum { T_GUARD = 64 };
static long t_guard_bad;
static void *t_gd_malloc(size_t n)
{
    unsigned char *p = malloc(n + 2 * T_GUARD);
    if (!p) return NULL;
    memset(p, 0xFB, T_GUARD), memset(p + T_GUARD + n, 0xFB, T_GUARD);
    memcpy(p, &n, sizeof n);
    return p + T_GUARD;
}
static void t_gd_free(void *q)
{
    if (!q) return;
    unsigned char *p = (unsigned char *)q - T_GUARD;
    size_t n;
    memcpy(&n, p, sizeof n);
    int bad = 0;
    for (size_t i = sizeof n; i < T_GUARD; i++) bad |= p[i] != 0xFB;
    for (size_t i = 0; i < T_GUARD; i++) bad |= p[T_GUARD + n + i] != 0xFB;
    t_guard_bad += bad;
    free(p);
}
static void *t_gd_calloc(size_t a, size_t b)
{
    void *p = t_gd_malloc(a * b);
    if (p) memset(p, 0, a * b);
    return p;
}
static void *t_gd_realloc(void *q, size_t n)
{
    void *p = t_gd_malloc(n);
    if (p && q) {
        size_t o;
        memcpy(&o, (unsigned char *)q - T_GUARD, sizeof o);
        memcpy(p, q, o < n ? o : n);
    }
    if (p) t_gd_free(q);
    return p;
}
#define malloc(n) t_gd_malloc(n)
#define calloc(a, b) t_gd_calloc(a, b)
#define realloc(p, n) t_gd_realloc(p, n)
#define free(p) t_gd_free(p)
#include "../include/mk_dec.h"
#undef malloc
#undef calloc
#undef realloc
#undef free
#include "t.h"
#include "tmesh.h"

static int closed(int nv, const int *T, int nt, int *chi)
{
    long long *e = malloc(sizeof *e * 3 * (size_t)nt);
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) e[3 * t + k] = (long long)T[3 * t + k] * nv + T[3 * t + (k + 1) % 3];
    for (long g = 3L * nt / 2; g > 0; g /= 2)
        for (long i = g; i < 3L * nt; i++)
            for (long j = i; j >= g && e[j - g] > e[j]; j -= g) { const long long x = e[j]; e[j] = e[j - g], e[j - g] = x; }
    int ok = nt > 0;
    for (long i = 1; i < 3L * nt; i++) ok &= e[i] != e[i - 1];
    for (long i = 0; i < 3L * nt && ok; i++) {
        const long long r = (e[i] % nv) * nv + e[i] / nv;
        long lo = 0, hi = 3L * nt - 1, f = 0;
        while (lo <= hi) { const long m = (lo + hi) / 2; if (e[m] == r) { f = 1; break; } if (e[m] < r) lo = m + 1; else hi = m - 1; }
        ok &= f;
    }
    char *used = calloc((size_t)nv, 1);
    int nu = 0;
    for (int i = 0; i < 3 * nt; i++) nu += !used[T[i]], used[T[i]] = 1;
    *chi = nu - 3 * nt / 2 + nt;
    free(e), free(used);
    return ok;
}
/* the squared distance from p to triangle abc, by projection and the three edges (independent of mk_geom) */
static double seg2(const double *p, const double *a, const double *b)
{
    double ab[3], t = 0, l = 0, d = 0;
    for (int k = 0; k < 3; k++) ab[k] = b[k] - a[k], t += ab[k] * (p[k] - a[k]), l += ab[k] * ab[k];
    t = l > 0 ? fmin(1, fmax(0, t / l)) : 0;
    for (int k = 0; k < 3; k++) { const double e = a[k] + t * ab[k] - p[k]; d += e * e; }
    return d;
}
static double tri2(const double *p, const double *a, const double *b, const double *c)
{
    double e1[3], e2[3], n[3], best = fmin(seg2(p, a, b), fmin(seg2(p, b, c), seg2(p, c, a)));
    for (int k = 0; k < 3; k++) e1[k] = b[k] - a[k], e2[k] = c[k] - a[k];
    n[0] = e1[1] * e2[2] - e1[2] * e2[1], n[1] = e1[2] * e2[0] - e1[0] * e2[2], n[2] = e1[0] * e2[1] - e1[1] * e2[0];
    const double nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    if (nn > 0) {
        const double hh = ((p[0] - a[0]) * n[0] + (p[1] - a[1]) * n[1] + (p[2] - a[2]) * n[2]) / nn;
        double q[3], s[3];
        for (int k = 0; k < 3; k++) q[k] = p[k] - hh * n[k];
        const double *v[3] = {a, b, c};
        for (int i = 0; i < 3; i++) {
            const double *u = v[i], *w = v[(i + 1) % 3];
            const double eu[3] = {w[0] - u[0], w[1] - u[1], w[2] - u[2]}, qu[3] = {q[0] - u[0], q[1] - u[1], q[2] - u[2]};
            s[i] = (eu[1] * qu[2] - eu[2] * qu[1]) * n[0] + (eu[2] * qu[0] - eu[0] * qu[2]) * n[1] + (eu[0] * qu[1] - eu[1] * qu[0]) * n[2];
        }
        if ((s[0] >= 0 && s[1] >= 0 && s[2] >= 0) || (s[0] <= 0 && s[1] <= 0 && s[2] <= 0)) best = fmin(best, hh * hh * nn);
    }
    return best;
}
/* the largest distance of the input's vertices from the result, and of points on the result (corners, edge midpoints,
   centres) from the input, brute force */
static void drift(const double *V, int nv, const int *T, int nt, const int *To, int nto, double *in2out, double *out2in)
{
    double a = 0, b = 0;
    for (int i = 0; i < nv; i++) {
        double d = INFINITY;
        for (int t = 0; t < nto; t++) d = fmin(d, tri2(V + 3 * i, V + 3 * To[3 * t], V + 3 * To[3 * t + 1], V + 3 * To[3 * t + 2]));
        a = fmax(a, d);
    }
    for (int t = 0; t < nto; t++) {
        const double *p0 = V + 3 * To[3 * t], *p1 = V + 3 * To[3 * t + 1], *p2 = V + 3 * To[3 * t + 2];
        double s[4][3];
        for (int k = 0; k < 3; k++) s[0][k] = (p0[k] + p1[k]) / 2, s[1][k] = (p1[k] + p2[k]) / 2, s[2][k] = (p2[k] + p0[k]) / 2, s[3][k] = (p0[k] + p1[k] + p2[k]) / 3;
        for (int j = 0; j < 4; j++) {
            double d = INFINITY;
            for (int u = 0; u < nt; u++) d = fmin(d, tri2(s[j], V + 3 * T[3 * u], V + 3 * T[3 * u + 1], V + 3 * T[3 * u + 2]));
            b = fmax(b, d);
        }
    }
    *in2out = sqrt(a), *out2in = sqrt(b);
}
/* an edge of one triangle through the other's interior, triangles sharing no corner (Moller-Trumbore, strict) */
static int segtri(const double *p, const double *q, const double *a, const double *b, const double *c)
{
    double d[3], e1[3], e2[3], h[3], s[3], qq[3];
    for (int k = 0; k < 3; k++) d[k] = q[k] - p[k], e1[k] = b[k] - a[k], e2[k] = c[k] - a[k], s[k] = p[k] - a[k];
    h[0] = d[1] * e2[2] - d[2] * e2[1], h[1] = d[2] * e2[0] - d[0] * e2[2], h[2] = d[0] * e2[1] - d[1] * e2[0];
    const double det = e1[0] * h[0] + e1[1] * h[1] + e1[2] * h[2];
    if (fabs(det) < 1e-300) return 0;
    const double u = (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]) / det;
    if (u <= 0 || u >= 1) return 0;
    qq[0] = s[1] * e1[2] - s[2] * e1[1], qq[1] = s[2] * e1[0] - s[0] * e1[2], qq[2] = s[0] * e1[1] - s[1] * e1[0];
    const double v = (d[0] * qq[0] + d[1] * qq[1] + d[2] * qq[2]) / det, t = (e2[0] * qq[0] + e2[1] * qq[1] + e2[2] * qq[2]) / det;
    return v > 0 && u + v < 1 && t > 0 && t < 1;
}
static long crossings(const double *V, const int *T, int nt)
{
    long n = 0;
    for (int i = 0; i < nt; i++)
        for (int j = i + 1; j < nt; j++) {
            const int *a = T + 3 * i, *b = T + 3 * j;
            int sh = 0;
            for (int x = 0; x < 3; x++)
                for (int y = 0; y < 3; y++) sh |= a[x] == b[y];
            if (sh) continue;
            int c = 0;
            for (int k = 0; k < 3 && !c; k++)
                c = segtri(V + 3 * a[k], V + 3 * a[(k + 1) % 3], V + 3 * b[0], V + 3 * b[1], V + 3 * b[2]) ||
                    segtri(V + 3 * b[k], V + 3 * b[(k + 1) % 3], V + 3 * a[0], V + 3 * a[1], V + 3 * a[2]);
            n += c;
        }
    return n;
}
static int half(const double *p, void *ctx) { (void)ctx; return p[2] > 0; }

int main(int argc, char **argv)
{
    tm_dir = argc > 1 ? argv[1] : NULL;
    const double c0[3] = {0, 0, 0};
    mk_dec_opts o;
    mk_dec_report r;
    { /* a unit sphere, a bound of 0.005 */
        tm_mesh S = {0};
        tm_sphere(&S, 1, c0, 5);
        mk_dec_defaults(&o, 1), o.dist = 0.005, o.maxedge = 0, o.minalt = 1e-6;
        unsigned char *al = malloc((size_t)S.nv);
        int *To = malloc(sizeof(int) * 3 * (size_t)S.nt), nto = 0, chi = 0;
        const int e = mk_decimate(S.V, S.nv, S.T, S.nt, &o, al, To, &nto, &r);
        double a, b;
        drift(S.V, S.nv, S.T, S.nt, To, nto, &a, &b);
        int ok = e == 0 && closed(S.nv, To, nto, &chi);
        CHECK(ok && chi == 2 && r.nv < S.nv / 4 && r.refused[5] + r.refused[6] > 0, "a unit sphere, bound 0.005: %d -> %d vertices, closed, Euler %d; the bound refused %ld collapses", S.nv, r.nv, chi,
              r.refused[5] + r.refused[6]);
        CHECK(ok && a <= 0.005 * (1 + 1e-9) && b <= 0.005 * (1 + 1e-9), "every input vertex within %.5f of the result, the result within %.5f of the input (bound 0.005)", a, b);
        int bad = 0;
        for (int t = 0; t < 3 * nto; t++) bad += !al[To[t]];
        CHECK(ok && !bad && r.nt == nto && r.nv == nto / 2 + 2, "kept vertices are the input's own, the faces use only them (%d of %d)", r.nv, S.nv);
        tm_write("dec_sphere", S.V, S.nv, To, nto, NULL);
        /* the target, and held vertices */
        mk_dec_report r2;
        o.target = 2000;
        mk_decimate(S.V, S.nv, S.T, S.nt, &o, al, To, &nto, &r2);
        CHECK(r2.nv == 2000, "a target of 2000 vertices: %d", r2.nv);
        unsigned char *keep = calloc((size_t)S.nv, 1);
        for (int i = 0; i < S.nv; i += 7) keep[i] = i % 2 ? 1 : 3;
        o.target = 0, o.keep = keep;
        mk_decimate(S.V, S.nv, S.T, S.nt, &o, al, To, &nto, &r2);
        int lost = 0, nk = 0;
        for (int i = 0; i < S.nv; i++) nk += keep[i] != 0, lost += keep[i] && !al[i];
        CHECK(lost == 0 && r2.nv < r.nv + nk, "%d held vertices (1 and 3) all kept; %d vertices", nk, r2.nv);
        o.keep = NULL, free(keep);
        /* regions: the fold sums by a caller's function */
        o.nregion = 2, o.region = half;
        const int er = mk_decimate(S.V, S.nv, S.T, S.nt, &o, al, To, &nto, &r2);
        CHECK(er == 0 && r2.area0[0] > 0 && r2.area0[1] > 0 && fabs(r2.area0[0] - r2.area0[1]) < 0.03 * r2.area0[0], "two regions (z <= 0, z > 0): areas %.4f and %.4f", r2.area0[0], r2.area0[1]);
        o.nregion = 1, o.region = NULL;
        free(al), free(To), tm_free(&S);
    }
    { /* a torus, Euler 0 */
        tm_mesh S = {0};
        tm_torus(&S, 2, 0.7, 160, 64);
        mk_dec_defaults(&o, 1), o.dist = 0.01, o.maxedge = 0, o.minalt = 1e-6;
        unsigned char *al = malloc((size_t)S.nv);
        int *To = malloc(sizeof(int) * 3 * (size_t)S.nt), nto = 0, chi = 1;
        const int e = mk_decimate(S.V, S.nv, S.T, S.nt, &o, al, To, &nto, &r);
        double a, b;
        drift(S.V, S.nv, S.T, S.nt, To, nto, &a, &b);
        const int ok = e == 0 && closed(S.nv, To, nto, &chi);
        CHECK(ok && chi == 0 && r.nv < S.nv / 4 && a <= 0.01 * (1 + 1e-9) && b <= 0.01 * (1 + 1e-9), "a torus, bound 0.01: %d -> %d vertices, closed, Euler %d, drift %.5f and %.5f", S.nv, r.nv, chi, a, b);
        tm_write("dec_torus", S.V, S.nv, To, nto, NULL);
        free(al), free(To), tm_free(&S);
    }
    { /* a torus of 12 x 3: with no bound and no turn limit only the link condition keeps it a torus */
        tm_mesh S = {0};
        tm_torus(&S, 2, 0.7, 12, 3);
        mk_dec_defaults(&o, 1), o.dist = 0, o.maxedge = 0, o.minangle = 0, o.turn = -1, o.minalt = 0, o.clean = 0, o.fold = 0;
        unsigned char *al = malloc((size_t)S.nv);
        int *To = malloc(sizeof(int) * 3 * (size_t)S.nt), nto = 0, chi = 1;
        const int e = mk_decimate(S.V, S.nv, S.T, S.nt, &o, al, To, &nto, &r);
        CHECK(e == 0 && closed(S.nv, To, nto, &chi) && chi == 0 && r.refused[0] > 0, "a 12 x 3 torus, no bound, no turn limit: %d -> %d vertices, still a closed torus (Euler %d); %ld collapses refused by the link condition",
              S.nv, r.nv, chi, r.refused[0]);
        free(al), free(To), tm_free(&S);
    }
    { /* a hollow ball: an outer sphere and an inner one wound inward, 0.03 apart, so outer chords reach the inner */
        tm_mesh S = {0};
        tm_sphere(&S, 1, c0, 4);
        const int n1 = S.nv, t1 = S.nt;
        tm_sphere(&S, 0.97, c0, 4);
        for (int t = t1; t < S.nt; t++) { const int x = S.T[3 * t + 1]; S.T[3 * t + 1] = S.T[3 * t + 2], S.T[3 * t + 2] = x; }
        (void)n1;
        mk_dec_defaults(&o, 1), o.dist = 0.2, o.maxedge = 0, o.minangle = 0, o.minalt = 0, o.clean = 0;
        unsigned char *al = malloc((size_t)S.nv);
        int *To = malloc(sizeof(int) * 3 * (size_t)S.nt), nto = 0, chi = 0;
        const int e = mk_decimate(S.V, S.nv, S.T, S.nt, &o, al, To, &nto, &r);
        const long nx = e ? -1 : crossings(S.V, To, nto);
        CHECK(e == 0 && closed(S.nv, To, nto, &chi) && chi == 4 && nx == 0 && r.refused[7] > 0, "a hollow ball, walls 0.03 apart, bound 0.2: %d -> %d vertices, Euler %d, %ld crossing pairs; %ld collapses refused for a crossing",
              S.nv, r.nv, chi, nx, r.refused[7]);
        tm_write("dec_hollow", S.V, S.nv, To, nto, NULL);
        free(al), free(To), tm_free(&S);
    }
    { /* a bumpy sphere: fold lengths over 90 and over 80 degrees can only fall */
        tm_mesh S = {0};
        tm_sphere(&S, 1, c0, 5);
        unsigned long long st = 12345;
        for (int i = 0; i < S.nv; i++) {
            st ^= st << 13, st ^= st >> 7, st ^= st << 17;
            const double f = 1 + 0.02 * ((double)(st >> 11) * 0x1p-53 - 0.5);
            for (int k = 0; k < 3; k++) S.V[3 * i + k] *= f;
        }
        mk_dec_defaults(&o, 1), o.dist = 0.05, o.maxedge = 0, o.minangle = 0, o.minalt = 0, o.turn = -1, o.clean = 0;
        unsigned char *al = malloc((size_t)S.nv);
        int *To = malloc(sizeof(int) * 3 * (size_t)S.nt), nto = 0;
        const int e = mk_decimate(S.V, S.nv, S.T, S.nt, &o, al, To, &nto, &r);
        CHECK(e == 0 && r.fold1[0][0] <= r.fold0[0][0] && r.fold1[0][1] <= r.fold0[0][1] && r.refused[4] > 0,
              "a bumpy sphere, no turn limit: fold over 90 degrees %.3f -> %.3f, over 80 %.3f -> %.3f (none made); %ld collapses refused for a fold", r.fold0[0][0], r.fold1[0][0], r.fold0[0][1], r.fold1[0][1],
              r.refused[4]);
        free(al), free(To), tm_free(&S);
    }
    { /* an open cap (the sphere above z = 0.5): its rim moves only along the rim */
        tm_mesh S = {0}, C = {0};
        tm_sphere(&S, 1, c0, 5);
        int *map = malloc(sizeof(int) * (size_t)S.nv);
        for (int i = 0; i < S.nv; i++) map[i] = -1;
        for (int t = 0; t < S.nt; t++) {
            int in = 1;
            for (int k = 0; k < 3; k++) in &= S.V[3 * S.T[3 * t + k] + 2] > 0.5;
            if (!in) continue;
            int q[3];
            for (int k = 0; k < 3; k++) {
                const int v = S.T[3 * t + k];
                if (map[v] < 0) map[v] = tm_vert(&C, S.V[3 * v], S.V[3 * v + 1], S.V[3 * v + 2]);
                q[k] = map[v];
            }
            tm_tri(&C, q[0], q[1], q[2]);
        }
        mk_dec_defaults(&o, 1), o.dist = 0.005, o.maxedge = 0, o.minalt = 1e-6;
        unsigned char *al = malloc((size_t)C.nv);
        int *To = malloc(sizeof(int) * 3 * (size_t)C.nt), nto = 0;
        const int e = mk_decimate(C.V, C.nv, C.T, C.nt, &o, al, To, &nto, &r);
        /* rim vertices: ends of an edge used once; the result's must be the input's */
        unsigned char *rim0 = calloc((size_t)C.nv, 1), *rim1 = calloc((size_t)C.nv, 1);
        for (int pass = 0; pass < 2; pass++) {
            const int *F = pass ? To : C.T, nf = pass ? nto : C.nt;
            unsigned char *rm = pass ? rim1 : rim0;
            for (int t = 0; t < nf && !e; t++)
                for (int k = 0; k < 3; k++) {
                    const int a = F[3 * t + k], bb = F[3 * t + (k + 1) % 3];
                    int found = 0;
                    for (int u = 0; u < nf && !found; u++)
                        for (int j = 0; j < 3; j++) found |= F[3 * u + j] == bb && F[3 * u + (j + 1) % 3] == a;
                    if (!found) rm[a] = rm[bb] = 1;
                }
        }
        long offrim = 0;
        for (int i = 0; i < C.nv; i++) offrim += rim1[i] && !rim0[i];
        free(rim0), free(rim1);
        double a, b;
        drift(C.V, C.nv, C.T, C.nt, To, nto, &a, &b);
        CHECK(e == 0 && r.rim1 > 0 && r.rim1 < r.rim0 && a <= 0.005 * (1 + 1e-9) && b <= 0.005 * (1 + 1e-9) && offrim == 0, "an open cap: %d -> %d vertices, rim %ld -> %ld (every one of the input's rim, %ld not), drift %.5f and %.5f", C.nv, r.nv, r.rim0,
              r.rim1, offrim, a, b);
        free(map), free(al), free(To), tm_free(&S), tm_free(&C);
    }
    { /* a needle ear on an open grid's rim (its only triangle (w, v, u), u 1e-4 from v), every grid vertex held, a
         bound and the clean-up on: only the relaxed try of u onto v can take the needle, and no new triangle is left to
         own u, so the nearest remaining face of v's must own it (with none, the owner index is -1 and the write lands
         before the owner table) */
        tm_mesh G = {0};
        const int n = 6;
        for (int j = 0; j <= n; j++)
            for (int i = 0; i <= n; i++) tm_vert(&G, i, j, 0);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                const int a = j * (n + 1) + i, b = a + 1, c = a + n + 2, d = a + n + 1;
                tm_tri(&G, a, b, c), tm_tri(&G, a, c, d);
            }
        const int v = 2, w = 3, u = tm_vert(&G, 2 + 1e-4, -1e-5, 0);
        tm_tri(&G, w, v, u);
        unsigned char *keep = malloc((size_t)G.nv), *al = malloc((size_t)G.nv);
        for (int i = 0; i < G.nv; i++) keep[i] = i != u;
        mk_dec_defaults(&o, 1), o.dist = 0.01, o.maxedge = 0, o.minalt = 1e-3, o.keep = keep;
        int *To = malloc(sizeof(int) * 3 * (size_t)G.nt), nto = 0;
        const int e = mk_decimate(G.V, G.nv, G.T, G.nt, &o, al, To, &nto, &r);
        long lost = 0;
        for (int i = 0; i < G.nv; i++) lost += i != u && !al[i];
        CHECK(e == 0 && !al[u] && lost == 0 && nto == G.nt - 1 && r.clean_after == 1 && r.left_after == 0 && r.relaxed == 1 && r.maxown < 2e-5,
              "a needle ear on an open grid, every other vertex held: taken by the relaxed clean-up (%ld, %ld relaxed, %ld left), %d -> %d faces; its tip owned by v's nearest face, %.2g off",
              r.clean_after, r.relaxed, r.left_after, G.nt, nto, r.maxown);
        free(keep), free(al), free(To), tm_free(&G);
    }
    {
        double V[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
        int T[12] = {0, 2, 1, 0, 1, 3, 1, 2, 3, 0, 3, 2}, Tb[3] = {0, 1, 4}, Tr[3] = {0, 1, 1}, To[12], nto;
        unsigned char al[4];
        mk_dec_opts b;
        int k = 0, n = 0;
#define REF(call) (k++, n += (call) == -1)
        mk_dec_defaults(&b, 1);
        REF(mk_decimate(V, 0, T, 4, &b, al, To, &nto, &r));
        REF(mk_decimate(V, 4, T, 0, &b, al, To, &nto, &r));
        REF(mk_decimate(V, 4, Tb, 1, &b, al, To, &nto, &r));
        REF(mk_decimate(V, 4, Tr, 1, &b, al, To, &nto, &r));
        REF(mk_decimate(V, 4, T, 4, &b, NULL, To, &nto, &r));
        b.dist = -1, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), b.dist = 0.3;
        b.dist = NAN, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), b.dist = 0.3;
        b.fold = 180, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), b.fold = 80;
        b.turn = 2, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), b.turn = 0.2;
        b.nregion = 2, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), b.nregion = 9, b.region = half, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), b.nregion = 1, b.region = NULL;
        b.target = -1, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), b.target = 0;
        V[4] = NAN, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), V[4] = 0;
        V[3] = 1e300, REF(mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r)), V[3] = 1;
#undef REF
        const int e = mk_decimate(V, 4, T, 4, &b, al, To, &nto, &r);
        CHECK(n == k && e == 0 && r.nv == 4, "refusals: %d of %d (sizes, an index out of range, a repeated corner, no output, a bound negative or NaN, fold 180, turn 2, regions without a function or past 8, target -1, a NaN or 1e300 coordinate); a tetrahedron stays (%d vertices)",
              n, k, r.nv);
    }
    { /* the neighbour list taken right after a link test still holds the edge's two common vertices, which
         mk_dec_cheap_ stamps MARKS + 1: mk_dec_nbrs_ must not read them as seen (else they go unrated after every
         collapse, and about half of all neighbour lists come up short). On an octahedron, through the internals. */
        const double OV[18] = {1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};
        int OF[24] = {0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5}, nb[MK_DEC_MAXN];
        mk_dec_ c = {0};
        c.P = OV, c.nv = 6, c.nt = 8, c.nreg = 1, c.relax = 1, c.F = OF;
        c.VT = t_gd_calloc(6, sizeof *c.VT), c.NVT = t_gd_calloc(6, sizeof *c.NVT), c.CVT = t_gd_calloc(6, sizeof *c.CVT);
        c.MARK = t_gd_calloc(6, sizeof *c.MARK), c.deadv = t_gd_calloc(6, 1), c.bnd = t_gd_calloc(6, 1), c.N0 = t_gd_calloc(18, sizeof *c.N0);
        for (int t = 0; t < 8; t++)
            for (int k = 0; k < 3; k++) mk_dec_fan_add_(&c, OF[3 * t + k], t);
        const int ok = !c.oom && mk_dec_cheap_(&c, 0, 4), n = mk_dec_nbrs_(&c, 4, nb);
        CHECK(ok && n == 4, "an octahedron's edge passes the link test (%d), and its end's next neighbour list holds all %d of 4", ok, n);
        for (int i = 0; i < 6; i++) t_gd_free(c.VT[i]);
        t_gd_free(c.VT), t_gd_free(c.NVT), t_gd_free(c.CVT), t_gd_free(c.MARK), t_gd_free(c.deadv), t_gd_free(c.bnd), t_gd_free(c.N0);
    }
    CHECK(t_guard_bad == 0, "every block the decimation allocated kept its guard bytes: %ld written over", t_guard_bad);
    return wrong;
}
