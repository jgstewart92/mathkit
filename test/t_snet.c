/* mk_snet on known shapes (s-forms' rows, 2026-11-15, plus the kit's): a ball raw and relaxed, a ball split into two labels,
   a torus, a ball from a scalar field, two voxels touching at an edge, the refusals. With a folder argument it also writes the
   surfaces as OBJ files there, for test/check.sh to give to mesher's meshcheck. */
#include "../include/mk_snet.h"
#include "t.h"
#include <stdint.h>
#include <stdio.h>

/* one label's triangles: every directed edge once and its reverse once (closed and consistently wound), the Euler
   characteristic of the used vertices, the signed volume */
static void surf(const mk_snet *s, const int *T, int nt, int *closed, int *chi, double *vol)
{
    long ne = 3L * nt;
    long long *e = malloc(sizeof *e * (size_t)(ne ? ne : 1));
    char *used = calloc((size_t)s->nv + 1, 1);
    double v = 0;
    for (int t = 0; t < nt; t++) {
        for (int k = 0; k < 3; k++) e[3L * t + k] = (long long)T[3 * t + k] * s->nv + T[3 * t + (k + 1) % 3], used[T[3 * t + k]] = 1;
        const double *a = s->V + 3 * T[3 * t], *b = s->V + 3 * T[3 * t + 1], *c = s->V + 3 * T[3 * t + 2];
        v += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
    }
    for (long gap = ne / 2; gap > 0; gap /= 2) /* shell sort: the test stays independent of the header */
        for (long i = gap; i < ne; i++)
            for (long j = i; j >= gap && e[j - gap] > e[j]; j -= gap) { long long x = e[j]; e[j] = e[j - gap], e[j - gap] = x; }
    int ok = ne > 0;
    for (long i = 1; i < ne; i++) ok &= e[i] != e[i - 1];
    for (long i = 0; i < ne && ok; i++) {
        const long long a = e[i] / s->nv, b = e[i] % s->nv, r = b * s->nv + a;
        long lo = 0, hi = ne - 1, f = 0;
        while (lo <= hi) { long m = (lo + hi) / 2; if (e[m] == r) { f = 1; break; } if (e[m] < r) lo = m + 1; else hi = m - 1; }
        ok &= f;
    }
    int nvu = 0;
    for (int i = 0; i < s->nv; i++) nvu += used[i];
    *closed = ok, *chi = nvu - (int)(ne / 2) + nt, *vol = v;
    free(e), free(used);
}
static const char *odir;
static void obj(const char *name, const mk_snet *s, const int *T, int nt)
{
    if (!odir) return;
    char p[1024];
    snprintf(p, sizeof p, "%s/%s.obj", odir, name);
    FILE *f = fopen(p, "w");
    if (!f) return;
    for (int v = 0; v < s->nv; v++) fprintf(f, "v %.17g %.17g %.17g\n", s->V[3 * v], s->V[3 * v + 1], s->V[3 * v + 2]);
    for (int t = 0; t < nt; t++) fprintf(f, "f %d %d %d\n", T[3 * t] + 1, T[3 * t + 1] + 1, T[3 * t + 2] + 1);
    fclose(f);
}

int main(int argc, char **argv)
{
    odir = argc > 1 ? argv[1] : NULL;
    const int n = 28;
    const double R = 10, c0 = 13.5, org[3] = {-1, 2, 0.5}, h = 0.25, exact = 4.0 / 3 * T_PI * R * R * R * h * h * h;
    uint8_t *g = calloc((size_t)n * n * n, 1);
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                const double x = i + 0.5 - c0, y = j + 0.5 - c0, z = k + 0.5 - c0;
                g[((size_t)k * n + j) * n + i] = x * x + y * y + z * z < R * R ? 1 : 0;
            }
    double err[2], vol[2];
    for (int pass = 0; pass < 2; pass++) {
        mk_snet s;
        CHECK(mk_snet_build(g, n, n, n, org, h, pass ? 20 : 0, &s) == 0, "ball: built (%s)", pass ? "20 relaxations" : "raw");
        int *T = NULL, nt = mk_snet_faces(&s, g, 1, &T), closed, chi;
        surf(&s, T, nt, &closed, &chi, &vol[pass]);
        obj(pass ? "ball_relaxed" : "ball_raw", &s, T, nt);
        CHECK(closed, "ball (%s): %d triangles, every edge shared by two, wound one way", pass ? "relaxed" : "raw", nt);
        CHECK(chi == 2 && mk_snet_split(&s, g, 1) == 0, "ball (%s): Euler characteristic %d (a sphere's 2), no split cells", pass ? "relaxed" : "raw", chi);
        err[pass] = 0;
        int inside = 1;
        for (int v = 0; v < s.nv; v++) {
            const double *p = s.V + 3 * v;
            double r2 = 0;
            for (int a = 0; a < 3; a++) { const double q = (p[a] - org[a]) / h - c0; r2 += q * q; }
            err[pass] = fmax(err[pass], fabs(sqrt(r2) - R));
        }
        int found = 0;
        for (int k = 0; k <= n; k++) /* every vertex inside its own dual cell, found there by mk_snet_vertex */
            for (int j = 0; j <= n; j++)
                for (int i = 0; i <= n; i++) {
                    const int v = mk_snet_vertex(&s, i, j, k);
                    if (v < 0) continue;
                    found++;
                    const int ci[3] = {i, j, k};
                    for (int a = 0; a < 3; a++) { const double u = (s.V[3 * v + a] - org[a]) / h - (ci[a] - 0.5); inside &= u >= -1e-12 && u <= 1 + 1e-12; }
                }
        CHECK(inside && found == s.nv, "ball (%s): every vertex inside its dual cell (%d of %d found by cell)", pass ? "relaxed" : "raw", found, s.nv);
        CHECK(vol[pass] > 0 && fabs(vol[pass] / exact - 1) < 0.03, "ball (%s): volume %.4f of the sphere's %.4f (%+.2f %%), outward", pass ? "relaxed" : "raw", vol[pass], exact, 100 * (vol[pass] / exact - 1));
        free(T), mk_snet_free(&s);
    }
    CHECK(err[1] < err[0] && err[1] < 0.75, "relaxing brings the vertices to the sphere: largest radial error %.3f -> %.3f cells", err[0], err[1]);
    /* the same ball from its distance field: the crossings of iso put the raw vertices on the sphere */
    {
        float *phi = malloc(sizeof(float) * (size_t)n * n * n);
        for (int k = 0; k < n; k++)
            for (int j = 0; j < n; j++)
                for (int i = 0; i < n; i++) {
                    const double x = i + 0.5 - c0, y = j + 0.5 - c0, z = k + 0.5 - c0;
                    phi[((size_t)k * n + j) * n + i] = (float)(sqrt(x * x + y * y + z * z) - R);
                }
        mk_snet s;
        const int rb = mk_snet_build_iso(phi, n, n, n, org, h, 0, 0, &s);
        int *T = NULL, nt = mk_snet_faces_iso(&s, phi, 0, &T), closed, chi;
        double v, ei = 0;
        surf(&s, T, nt, &closed, &chi, &v);
        obj("ball_iso", &s, T, nt);
        for (int q = 0; q < s.nv; q++) {
            double r2 = 0;
            for (int a = 0; a < 3; a++) { const double u = (s.V[3 * q + a] - org[a]) / h - c0; r2 += u * u; }
            ei = fmax(ei, fabs(sqrt(r2) - R));
        }
        CHECK(rb == 0 && closed && chi == 2 && v > 0 && fabs(v / exact - 1) < 0.02 && ei < 0.1 && ei < err[0] / 4,
              "ball from a distance field (iso 0, raw): closed, Euler %d, volume %+.2f %%, largest radial error %.3f cells (labels raw: %.3f)", chi, 100 * (v / exact - 1), ei, err[0]);
        mk_snet t;
        phi[5] = NAN;
        CHECK(mk_snet_build_iso(phi, n, n, n, org, h, 0, 0, &t) == -1 && !t.V && t.nv == 0, "a NaN in phi is refused, the net left empty");
        free(T), mk_snet_free(&s), free(phi);
    }
    /* two labels touching: the ball split at x < c0 (1) and x >= c0 (2); their interface triangles are the same, wound oppositely */
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                uint8_t *p = g + ((size_t)k * n + j) * n + i;
                if (*p && i + 0.5 >= c0) *p = 2;
            }
    int faces12 = 0;
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i + 1 < n; i++) faces12 += g[((size_t)k * n + j) * n + i] == 1 && g[((size_t)k * n + j) * n + i + 1] == 2;
    mk_snet s;
    mk_snet_build(g, n, n, n, org, h, 20, &s);
    int *T1 = NULL, *T2 = NULL, n1 = mk_snet_faces(&s, g, 1, &T1), n2 = mk_snet_faces(&s, g, 2, &T2), shared = 0, c1, c2, x1, x2;
    double v1, v2;
    surf(&s, T1, n1, &c1, &x1, &v1), surf(&s, T2, n2, &c2, &x2, &v2);
    obj("half1", &s, T1, n1), obj("half2", &s, T2, n2);
    for (int a = 0; a < n1; a++) /* a triangle of 1 whose reverse (any rotation) is a triangle of 2 */
        for (int b = 0; b < n2; b++) {
            const int *p = T1 + 3 * a, *q = T2 + 3 * b;
            int m = 0;
            for (int r = 0; r < 3 && !m; r++) m = p[0] == q[r] && p[1] == q[(r + 2) % 3] && p[2] == q[(r + 1) % 3];
            if (m) { shared++; break; }
        }
    CHECK(c1 && c2 && x1 == 2 && x2 == 2, "two halves: each closed with Euler characteristic 2 (%d, %d)", x1, x2);
    CHECK(shared == 2 * faces12, "two halves: %d of label 1's triangles are label 2's reversed, the %d interface faces' two each", shared, faces12);
    CHECK(fabs((v1 + v2) / exact - 1) < 0.03, "two halves: their volumes sum to the sphere's within 3 %% (%.4f against %.4f)", v1 + v2, exact);
    free(T1), free(T2), mk_snet_free(&s);
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                const double x = i + 0.5 - c0, y = j + 0.5 - c0, z = k + 0.5 - c0, q = sqrt(x * x + y * y) - 8;
                g[((size_t)k * n + j) * n + i] = q * q + z * z < 3.5 * 3.5 ? 3 : 0;
            }
    mk_snet_build(g, n, n, n, org, h, 20, &s);
    int *T3 = NULL, n3 = mk_snet_faces(&s, g, 3, &T3), c3, x3;
    double v3;
    surf(&s, T3, n3, &c3, &x3, &v3);
    obj("torus", &s, T3, n3);
    CHECK(c3 && x3 == 0 && mk_snet_split(&s, g, 3) == 0, "torus: closed, Euler characteristic %d (a torus's 0), no split cells", x3);
    free(T3), mk_snet_free(&s);
    /* two voxels touching at an edge: their surface meets itself there (four triangles on one edge), and the split count says so */
    {
        uint8_t d[4 * 4 * 4] = {0};
        d[(1 * 4 + 1) * 4 + 1] = 1, d[(1 * 4 + 2) * 4 + 2] = 1;
        mk_snet e;
        mk_snet_build(d, 4, 4, 4, org, h, 0, &e);
        int *T = NULL, nt = mk_snet_faces(&e, d, 1, &T), cl, chi;
        double v;
        surf(&e, T, nt, &cl, &chi, &v);
        obj("diagonal", &e, T, nt);
        const int sp = mk_snet_split(&e, d, 1);
        CHECK(sp > 0, "two voxels sharing only an edge: %d split cells reported (their surface touches itself there)", sp);
        free(T), mk_snet_free(&e);
    }
    {
        mk_snet e;
        const double bad[3] = {0, NAN, 0};
        const int r1 = mk_snet_build(g, 0, n, n, org, h, 0, &e), r2 = mk_snet_build(g, n, n, n, org, 0, 0, &e), r3 = mk_snet_build(g, n, n, n, bad, h, 0, &e);
        const int r4 = mk_snet_build(g, n, n, n, org, h, -1, &e), r5 = mk_snet_build(g, INT_MAX, 2, 2, org, h, 0, &e);
        CHECK(r1 == -1 && r2 == -1 && r3 == -1 && r4 == -1 && r5 == -1 && !e.V && !e.cell && e.nv == 0,
              "refused, the net left empty: nx 0, h 0, a NaN origin, iters -1, nx INT_MAX");
    }
    free(g);
#ifndef T_OLD /* 1.2.2 reads zs = NULL */
    { /* the verification work N-A (Result 87; scratch/the verification work/mk-new/findings/demo.c) */
        const uint8_t lab1[1] = {1};
        const double o0[3] = {0, 0, 0};
        mk_snet s;
        const int r = mk_snet_build(lab1, 1, 1, 1, o0, -1, 0, &s);
        CHECK(r == -1 && mk_snet_vertex(&s, 0, 0, 0) == -1, "N-A: asking an empty net (a refusal's) for a vertex gives -1");
    }
#endif
    { /* review 2 #14 (1.6.1): faces asked of a grid other than the net's own found dual cells with no vertex and read V[-3] (1.6.0) */
        enum { G = 16 };
        static uint8_t ga[G * G * G], gb[G * G * G];
        for (int k = 0; k < G; k++)
            for (int j = 0; j < G; j++)
                for (int i = 0; i < G; i++) {
                    const double x = i - 7.5, y = j - 7.5, z = k - 7.5, r2 = x * x + y * y + z * z;
                    ga[(k * G + j) * G + i] = r2 < 9, gb[(k * G + j) * G + i] = r2 < 25;
                }
        const double o0[3] = {0, 0, 0};
        mk_snet s;
        int *Ta = NULL, *Tb = NULL, na = -9, nb = -9;
        if (!mk_snet_build(ga, G, G, G, o0, 1, 0, &s)) na = mk_snet_faces(&s, ga, 1, &Ta), nb = mk_snet_faces(&s, gb, 1, &Tb);
        CHECK(na > 0 && nb == -1 && !Tb, "review 2 #14: a ball's net asked for the faces of a larger ball's grid gives %d (-1 wanted; its own grid %d triangles)", nb, na);
        free(Ta), free(Tb), mk_snet_free(&s);
    }
    return wrong;
}
