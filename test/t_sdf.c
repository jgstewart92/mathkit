/* mk_sdf on known meshes: a sphere and a torus against their analytic distance (near the surface and far from it), the
   sign at every node of three nested shells, a random closed mesh against brute force, an open mesh (the winding number
   takes the odd rows), two overlapping spheres with and without MK_SDF_FILL, the refusals. */
#include "../include/mk_sdf.h"
#include "t.h"
#include <stdlib.h>
#include <string.h>

typedef struct { double *V; int nv, cap; int *T; int nt; } msh;
static int addv(msh *m, double x, double y, double z)
{
    if (m->nv == m->cap) m->V = realloc(m->V, sizeof(double) * 3 * (size_t)(m->cap = m->cap ? 2 * m->cap : 1024));
    m->V[3 * m->nv] = x, m->V[3 * m->nv + 1] = y, m->V[3 * m->nv + 2] = z;
    return m->nv++;
}
static int mid(msh *m, int a, int b, int *key, int *val, int nk)
{
    const int lo = a < b ? a : b, hi = a < b ? b : a;
    int h = (int)(((unsigned)lo * 73856093u ^ (unsigned)hi * 19349663u) % (unsigned)nk);
    while (key[2 * h] >= 0) {
        if (key[2 * h] == lo && key[2 * h + 1] == hi) return val[h];
        h = (h + 1) % nk;
    }
    double x = (m->V[3 * a] + m->V[3 * b]) / 2, y = (m->V[3 * a + 1] + m->V[3 * b + 1]) / 2, z = (m->V[3 * a + 2] + m->V[3 * b + 2]) / 2;
    const double l = sqrt(x * x + y * y + z * z);
    key[2 * h] = lo, key[2 * h + 1] = hi;
    return val[h] = addv(m, x / l, y / l, z / l);
}
/* a unit icosphere, outward */
static void icosphere(msh *m, int sub)
{
    const double p = (1 + sqrt(5.0)) / 2, s = sqrt(1 + p * p);
    const double v[12][3] = {{-1, p, 0}, {1, p, 0}, {-1, -p, 0}, {1, -p, 0}, {0, -1, p}, {0, 1, p}, {0, -1, -p}, {0, 1, -p}, {p, 0, -1}, {p, 0, 1}, {-p, 0, -1}, {-p, 0, 1}};
    static const int f[60] = {0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
                              3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1};
    memset(m, 0, sizeof *m);
    for (int i = 0; i < 12; i++) addv(m, v[i][0] / s, v[i][1] / s, v[i][2] / s);
    m->T = malloc(sizeof(int) * 60), m->nt = 20, memcpy(m->T, f, sizeof f);
    for (int r = 0; r < sub; r++) {
        const int nk = 8 * m->nt + 1;
        int *key = malloc(sizeof(int) * 2 * (size_t)nk), *val = malloc(sizeof(int) * (size_t)nk);
        for (int i = 0; i < 2 * nk; i++) key[i] = -1;
        int *T2 = malloc(sizeof(int) * 12 * (size_t)(unsigned)m->nt); /* unsigned: gcc -O3 warns on a negative nt */
        for (int t = 0; t < m->nt; t++) {
            const int a = m->T[3 * t], b = m->T[3 * t + 1], c = m->T[3 * t + 2], ab = mid(m, a, b, key, val, nk), bc = mid(m, b, c, key, val, nk), ca = mid(m, c, a, key, val, nk);
            const int q[12] = {a, ab, ca, b, bc, ab, c, ca, bc, ab, bc, ca};
            memcpy(T2 + 12 * t, q, sizeof q);
        }
        free(m->T), m->T = T2, m->nt *= 4, free(key), free(val);
    }
}
/* append a sphere of radius r at c (flip: wound inward) */
static void add_sphere(msh *out, const msh *u, double r, const double c[3], int flip)
{
    const int base = out->nv;
    for (int v = 0; v < u->nv; v++) addv(out, c[0] + r * u->V[3 * v], c[1] + r * u->V[3 * v + 1], c[2] + r * u->V[3 * v + 2]);
    out->T = realloc(out->T, sizeof(int) * 3 * (size_t)(out->nt + u->nt));
    for (int t = 0; t < u->nt; t++) {
        int *q = out->T + 3 * (size_t)(out->nt + t);
        q[0] = base + u->T[3 * t], q[1] = base + u->T[3 * t + (flip ? 2 : 1)], q[2] = base + u->T[3 * t + (flip ? 1 : 2)];
    }
    out->nt += u->nt;
}
static void torus(msh *m, double R, double r, int nu, int nw)
{
    memset(m, 0, sizeof *m);
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nw; j++) {
            const double u = 2 * T_PI * i / nu, w = 2 * T_PI * j / nw;
            addv(m, (R + r * cos(w)) * cos(u), (R + r * cos(w)) * sin(u), r * sin(w));
        }
    m->T = malloc(sizeof(int) * 6 * (size_t)nu * nw);
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nw; j++) {
            const int a = i * nw + j, b = ((i + 1) % nu) * nw + j, c = ((i + 1) % nu) * nw + (j + 1) % nw, d = i * nw + (j + 1) % nw;
            int *q = m->T + 3 * (size_t)m->nt;
            q[0] = a, q[1] = b, q[2] = c, q[3] = a, q[4] = c, q[5] = d, m->nt += 2;
        }
}
static void node(const mk_sdf *g, int i, int j, int k, double p[3]) { p[0] = g->lo[0] + i * g->h, p[1] = g->lo[1] + j * g->h, p[2] = g->lo[2] + k * g->h; }
static float at(const mk_sdf *g, int i, int j, int k) { return g->d[(size_t)i + (size_t)g->n[0] * ((size_t)j + (size_t)g->n[1] * (size_t)k)]; }

/* the error of the grid against an analytic distance, near the surface (within two cells) and beyond */
typedef double (*exact_fn)(const double *p);
static double sd_unit(const double *p) { return sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) - 1; }
static double sd_tor(const double *p) { const double q = sqrt(p[0] * p[0] + p[1] * p[1]) - 0.6; return sqrt(q * q + p[2] * p[2]) - 0.25; }
static void errs(const mk_sdf *g, exact_fn f, double *en, double *ef)
{
    *en = *ef = 0;
    for (int k = 0; k < g->n[2]; k++)
        for (int j = 0; j < g->n[1]; j++)
            for (int i = 0; i < g->n[0]; i++) {
                double p[3];
                node(g, i, j, k, p);
                const double ex = f(p), e = fabs(at(g, i, j, k) - ex);
                if (fabs(ex) < 2 * g->h) *en = fmax(*en, e);
                else *ef = fmax(*ef, e);
            }
}
/* how far the faces stray from the analytic surface: the largest |f| at the centroids and edge midpoints */
static double sag(const msh *m, exact_fn f)
{
    double s = 0;
    for (int t = 0; t < m->nt; t++) {
        const double *a = m->V + 3 * m->T[3 * t], *b = m->V + 3 * m->T[3 * t + 1], *c = m->V + 3 * m->T[3 * t + 2];
        const double pts[4][3] = {{(a[0] + b[0] + c[0]) / 3, (a[1] + b[1] + c[1]) / 3, (a[2] + b[2] + c[2]) / 3}, {(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2},
                                  {(b[0] + c[0]) / 2, (b[1] + c[1]) / 2, (b[2] + c[2]) / 2}, {(c[0] + a[0]) / 2, (c[1] + a[1]) / 2, (c[2] + a[2]) / 2}};
        for (int q = 0; q < 4; q++) s = fmax(s, fabs(f(pts[q])));
    }
    return s;
}

int main(void)
{
    msh u;
    icosphere(&u, 5);
    { /* the unit sphere (20480 faces) and a torus R 0.6 r 0.25 (65536 faces): the analytic distance, which the faces
         miss by their sag; the fast sweeps (Batty, SDFGen, 2007) within 0.1 h of it beyond the band, MK_SDF_EXACT
         within the sag everywhere */
        msh m;
        torus(&m, 0.6, 0.25, 256, 128);
        const double s1 = sag(&u, sd_unit), s2 = sag(&m, sd_tor);
        mk_sdf g, gx, t, tx;
        const int r = mk_sdf_mesh(u.V, u.nv, u.T, u.nt, 0.02, 0.3, 0, &g) | mk_sdf_mesh(u.V, u.nv, u.T, u.nt, 0.02, 0.3, MK_SDF_EXACT, &gx) |
                      mk_sdf_mesh(m.V, m.nv, m.T, m.nt, 0.015, 0.2, 0, &t) | mk_sdf_mesh(m.V, m.nv, m.T, m.nt, 0.015, 0.2, MK_SDF_EXACT, &tx);
        double a[2], b[2], c[2], d[2];
        errs(&g, sd_unit, &a[0], &a[1]), errs(&gx, sd_unit, &b[0], &b[1]), errs(&t, sd_tor, &c[0], &c[1]), errs(&tx, sd_tor, &d[0], &d[1]);
        CHECK(r == 0 && a[0] < s1 + 1e-6 && a[1] < 0.1 * g.h && b[0] < s1 + 1e-6 && b[1] < s1 + 1e-6 && g.odd_rows == 0,
              "unit sphere, h 0.02 (sag %.2e): largest error near the surface %.2e, beyond %.2e; exact %.2e and %.2e", s1, a[0], a[1], b[0], b[1]);
        CHECK(r == 0 && c[0] < s2 + 1e-6 && c[1] < 0.1 * t.h && d[0] < s2 + 1e-6 && d[1] < s2 + 1e-6 && t.odd_rows == 0,
              "torus R 0.6 r 0.25, h 0.015 (sag %.2e): largest error near the surface %.2e, beyond %.2e; exact %.2e and %.2e", s2, c[0], c[1], d[0], d[1]);
        mk_sdf_free(&g), mk_sdf_free(&gx), mk_sdf_free(&t), mk_sdf_free(&tx), free(m.V), free(m.T);
    }
    { /* three nested shells: a ball of radius 0.25 in the hole of a shell 0.5 .. 1 (its inner face wound inward). The
         sign at every node off the faces' sag against the analytic solid, and the winding number at 1 node in 31 */
        msh m = {0};
        const double c0[3] = {0, 0, 0};
        icosphere(&u, 4);
        add_sphere(&m, &u, 1, c0, 0), add_sphere(&m, &u, 0.5, c0, 1), add_sphere(&m, &u, 0.25, c0, 0);
        mk_sdf g;
        mk_sdf_mesh(m.V, m.nv, m.T, m.nt, 0.03, 0.15, 0, &g);
        long bad = 0, tot = 0, wbad = 0;
        for (int k = 0; k < g.n[2]; k++)
            for (int j = 0; j < g.n[1]; j++)
                for (int i = 0; i < g.n[0]; i++) {
                    double p[3];
                    node(&g, i, j, k, p);
                    const double rr = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                    if (fabs(rr - 1) < 6e-3 || fabs(rr - 0.5) < 3e-3 || fabs(rr - 0.25) < 1.5e-3) continue;
                    const int in = rr < 0.25 || (rr > 0.5 && rr < 1);
                    tot++, bad += in != (at(&g, i, j, k) < 0);
                    if (tot % 31 == 0) wbad += (fabs(mk_winding_number(m.V, m.T, m.nt, p)) > 0.5) != (at(&g, i, j, k) < 0);
                }
        CHECK(bad == 0 && wbad == 0 && g.odd_rows == 0, "three nested shells: the sign wrong at %ld of %ld nodes; the winding number disagrees at %ld of the %ld it was asked",
              bad, tot, wbad, tot / 31);
        mk_sdf_free(&g), free(m.V), free(m.T);
    }
    { /* a random closed mesh (a sphere of 1280 faces, each vertex pushed in or out by up to 30 %) against brute force:
         the distance and the winding number's sign at every node */
        msh m = {0};
        const double c0[3] = {0, 0, 0};
        icosphere(&u, 3);
        add_sphere(&m, &u, 1, c0, 0);
        for (int v = 0; v < m.nv; v++) {
            const double s = 1 + 0.3 * t_rand();
            for (int a = 0; a < 3; a++) m.V[3 * v + a] *= s;
        }
        mk_sdf g, gx;
        mk_sdf_mesh(m.V, m.nv, m.T, m.nt, 0.1, 0.4, 0, &g), mk_sdf_mesh(m.V, m.nv, m.T, m.nt, 0.1, 0.4, MK_SDF_EXACT, &gx);
        double e = 0, ex = 0;
        long sb = 0, exact = 0, tot = 0;
        for (int k = 0; k < g.n[2]; k++)
            for (int j = 0; j < g.n[1]; j++)
                for (int i = 0; i < g.n[0]; i++) {
                    double p[3], best = INFINITY;
                    node(&g, i, j, k, p);
                    for (int f = 0; f < m.nt; f++) {
                        double w[3], q[3];
                        best = fmin(best, mk_closest_pt_tri(p, m.V + 3 * m.T[3 * f], m.V + 3 * m.T[3 * f + 1], m.V + 3 * m.T[3 * f + 2], w, q));
                    }
                    best = sqrt(best);
                    const float d = at(&g, i, j, k), dx = at(&gx, i, j, k);
                    const int in = fabs(mk_winding_number(m.V, m.T, m.nt, p)) > 0.5;
                    e = fmax(e, fabs(fabsf(d) - best)), ex = fmax(ex, fabs(fabsf(dx) - best)), sb += (in != (d < 0)) + (in != (dx < 0));
                    exact += fabsf(d) == (float)best, tot++;
                }
        CHECK(e < 2e-3 * g.h && ex <= 1e-7 * (1 + e) && sb == 0,
              "a random closed mesh (1280 faces), %ld nodes against brute force: sweeps off by %.2e at most (%ld nodes not exact to float), exact mode %.2e; sign wrong at %ld (winding); %ld near nodes",
              tot, e, tot - exact, ex, sb, g.near);
        mk_sdf_free(&g), mk_sdf_free(&gx), free(m.V), free(m.T);
    }
    { /* an open mesh: a sphere with one face removed. Rows through the hole cross it an odd number of times; their sign
         comes from the winding number, so the inside still reads inside */
        msh m = {0};
        const double c0[3] = {0, 0, 0};
        icosphere(&u, 3);
        add_sphere(&m, &u, 1, c0, 0);
        int far = 0; /* the face that faces +x most, removed */
        for (int t = 1; t < m.nt; t++) if (m.V[3 * m.T[3 * t]] + m.V[3 * m.T[3 * t + 1]] + m.V[3 * m.T[3 * t + 2]] > m.V[3 * m.T[3 * far]] + m.V[3 * m.T[3 * far + 1]] + m.V[3 * m.T[3 * far + 2]]) far = t;
        memmove(m.T + 3 * far, m.T + 3 * far + 3, sizeof(int) * 3 * (size_t)(m.nt - far - 1)), m.nt--;
        mk_sdf g;
        mk_sdf_mesh(m.V, m.nv, m.T, m.nt, 0.04, 0.2, 0, &g);
        long bad = 0, tot = 0;
        for (int k = 0; k < g.n[2]; k++)
            for (int j = 0; j < g.n[1]; j++)
                for (int i = 0; i < g.n[0]; i++) {
                    double p[3];
                    node(&g, i, j, k, p);
                    const double rr = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
                    if (fabs(rr - 1) < 0.05) continue;
                    tot++, bad += (rr < 1) != (at(&g, i, j, k) < 0);
                }
        CHECK(g.odd_rows > 0 && bad == 0, "a sphere with one face missing: %ld rows odd on three rays, signed by the winding number (%ld more even on a second ray); the sign wrong at %ld of %ld nodes", g.odd_rows, g.reray, bad, tot);
        mk_sdf_free(&g), free(m.V), free(m.T);
    }
    { /* a region covered twice: a sphere of radius 0.4 inside one of radius 1, both wound outward (winding 2 inside the
         small one). Parity counts four crossings there and calls it outside; MK_SDF_FILL makes the enclosed pocket
         inside, and leaves the outside alone */
        msh m = {0};
        const double c0[3] = {0, 0, 0};
        icosphere(&u, 3);
        add_sphere(&m, &u, 1, c0, 0), add_sphere(&m, &u, 0.4, c0, 0);
        mk_sdf g0, g1;
        mk_sdf_mesh(m.V, m.nv, m.T, m.nt, 0.04, 0.2, 0, &g0), mk_sdf_mesh(m.V, m.nv, m.T, m.nt, 0.04, 0.2, MK_SDF_FILL, &g1);
        const int ic = (int)lround(-g0.lo[0] / g0.h), jc = (int)lround(-g0.lo[1] / g0.h), kc = (int)lround(-g0.lo[2] / g0.h);
        long badout = 0, nout = 0;
        for (int k = 0; k < g1.n[2]; k++)
            for (int j = 0; j < g1.n[1]; j++)
                for (int i = 0; i < g1.n[0]; i++) {
                    double p[3];
                    node(&g1, i, j, k, p);
                    if (sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) > 1.05) nout++, badout += at(&g1, i, j, k) < 0;
                }
        CHECK(at(&g0, ic, jc, kc) > 0 && at(&g1, ic, jc, kc) < 0 && g1.filled > 0 && g0.filled == 0 && badout == 0,
              "a region covered twice: its centre reads %+.3f by parity, %+.3f with MK_SDF_FILL (%ld nodes filled; %ld of %ld outside nodes made inside)",
              at(&g0, ic, jc, kc), at(&g1, ic, jc, kc), g1.filled, badout, nout);
        mk_sdf_free(&g0), mk_sdf_free(&g1), free(m.V), free(m.T);
    }
    { /* a square prism turned 45 degrees, |x| + |y| < r, |z| < 1, with r set so a diagonal of nodes lies 2.1e-6 inside
         its slanted faces: nearer than the rays' offset (1.2e-5), so their rays meet the face before the node. The
         winding number signs those nodes; parity alone would call them outside */
        const double r = 0.5 + 1e-6, h = 0.1;
        const double V[24] = {r, 0, -1, 0, r, -1, -r, 0, -1, 0, -r, -1, r, 0, 1, 0, r, 1, -r, 0, 1, 0, -r, 1};
        const int T[36] = {0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0, 4, 3, 4, 7, 0, 3, 2, 0, 2, 1, 4, 5, 6, 4, 6, 7};
        mk_sdf g;
        mk_sdf_mesh(V, 8, T, 12, h, 0.2, 0, &g);
        long bad = 0, tot = 0;
        for (int k = 0; k < g.n[2]; k++)
            for (int j = 0; j < g.n[1]; j++)
                for (int i = 0; i < g.n[0]; i++) {
                    double p[3];
                    node(&g, i, j, k, p);
                    const double s = r - fabs(p[0]) - fabs(p[1]), c = 1 - fabs(p[2]);
                    if (fabs(s) < 1e-9 || fabs(c) < 1e-9) continue; /* on a face */
                    tot++, bad += (s > 0 && c > 0) != (at(&g, i, j, k) < 0);
                }
        CHECK(bad == 0 && g.near > 0, "a prism with nodes 2.1e-6 inside its slanted faces: the sign wrong at %ld of %ld nodes (%ld near nodes, each signed by a ray of its own)", bad, tot, g.near);
        mk_sdf_free(&g);
    }
    { /* the refusals, before anything is allocated */
        mk_sdf g;
        const int tri[3] = {0, 1, 2}, bad[3] = {0, 1, 3};
        const double v3[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0}, vn[9] = {0, 0, 0, NAN, 0, 0, 0, 1, 0};
        /* r1: 1e10 nodes along x; r2: 3e19 nodes in all */
        const int r1 = mk_sdf_mesh(v3, 3, tri, 1, 1e-10, 0, 0, &g), r2 = mk_sdf_mesh(v3, 3, tri, 1, 5e-7, 0.5, 0, &g);
        const int r3 = mk_sdf_mesh(v3, 3, bad, 1, 0.1, 0, 0, &g), r4 = mk_sdf_mesh(vn, 3, tri, 1, 0.1, 0, 0, &g), r5 = mk_sdf_mesh(v3, 3, tri, 1, 0, 0, 0, &g);
        const int r6 = mk_sdf_mesh(v3, 3, tri, 0, 0.1, 0, 0, &g), r7 = mk_sdf_mesh(v3, 3, tri, 1, 0.1, -1, 0, &g), r8 = mk_sdf_mesh(v3, 3, tri, 1, 0.1, 0, 0, &g);
        CHECK(r1 == -1 && r2 == -1 && r3 == -1 && r4 == -1 && r5 == -1 && r6 == -1 && r7 == -1 && r8 == 0 && g.n[0] == 11,
              "refused: 1e10 nodes on an axis, 3e19 nodes in all, an index out of range, a NaN vertex, h 0, no faces, margin -1; a sound call works");
        mk_sdf_free(&g);
    }
    { /* trilinear sampling between nodes */
        mk_sdf g;
        icosphere(&u, 4);
        mk_sdf_mesh(u.V, u.nv, u.T, u.nt, 0.05, 0.3, 0, &g);
        double e = 0;
        for (int r = 0; r < 200; r++) {
            const double p[3] = {1.2 * t_rand(), 1.2 * t_rand(), 1.2 * t_rand()};
            e = fmax(e, fabs(mk_sdf_sample(&g, p) - sd_unit(p)));
        }
        CHECK(e < 0.02, "trilinear samples at 200 points against the sphere's distance: largest error %.2e (h 0.05)", e);
        mk_sdf_free(&g);
    }
    { /* MK_SDF_LATTICE: two grids of one mesh with margins 0.13 and 0.31 share their nodes, and with MK_SDF_EXACT their
         distances there agree bit for bit */
        mk_sdf a, b;
        icosphere(&u, 3);
        const int ra = mk_sdf_mesh(u.V, u.nv, u.T, u.nt, 0.05, 0.13, MK_SDF_LATTICE | MK_SDF_EXACT, &a), rb = mk_sdf_mesh(u.V, u.nv, u.T, u.nt, 0.05, 0.31, MK_SDF_LATTICE | MK_SDF_EXACT, &b);
        int on = 1;
        for (int c = 0; c < 3; c++) on &= fabs(a.lo[c] / 0.05 - round(a.lo[c] / 0.05)) < 1e-9 && fabs(b.lo[c] / 0.05 - round(b.lo[c] / 0.05)) < 1e-9;
        long diff = 0, shared = 0;
        const int s[3] = {(int)lround((a.lo[0] - b.lo[0]) / 0.05), (int)lround((a.lo[1] - b.lo[1]) / 0.05), (int)lround((a.lo[2] - b.lo[2]) / 0.05)};
        for (int k = 0; k < a.n[2]; k++)
            for (int j = 0; j < a.n[1]; j++)
                for (int i = 0; i < a.n[0]; i++) shared++, diff += at(&a, i, j, k) != at(&b, i + s[0], j + s[1], k + s[2]);
        CHECK(ra == 0 && rb == 0 && on && s[0] > 0 && diff == 0, "MK_SDF_LATTICE: node 0 on the lattice for margins 0.13 and 0.31; %ld of %ld shared nodes differ (exact mode)", diff, shared);
        mk_sdf_free(&a), mk_sdf_free(&b);
    }
    { /* mk_sdf_mesh_box: a box across the surface and past the whole grid's side, and a box just off the surface (every
         node outside, its nearest faces a cell or two beyond its side), on the lattice, against the whole grid at their
         shared nodes: with MK_SDF_EXACT every one bit for bit; plain, every sign and every value within the band. A box
         away from the mesh: every node within the sag of the analytic distance (MK_SDF_EXACT). The refusals. */
        const double hb = 0.05, box[2][6] = {{0.32, -0.41, -0.26, 1.33, 0.37, 0.52}, {1.08, -0.21, -0.19, 1.4, 0.22, 0.23}}, far[6] = {2, 2, 2, 2.5, 2.6, 2.7};
        mk_sdf w, wx, b[2], bx[2], fx;
        const int r = mk_sdf_mesh(u.V, u.nv, u.T, u.nt, hb, 0.2, MK_SDF_LATTICE, &w) | mk_sdf_mesh(u.V, u.nv, u.T, u.nt, hb, 0.2, MK_SDF_LATTICE | MK_SDF_EXACT, &wx) |
                      mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, box[0], MK_SDF_LATTICE, &b[0]) | mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, box[0], MK_SDF_LATTICE | MK_SDF_EXACT, &bx[0]) |
                      mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, box[1], MK_SDF_LATTICE, &b[1]) | mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, box[1], MK_SDF_LATTICE | MK_SDF_EXACT, &bx[1]) |
                      mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, far, MK_SDF_EXACT, &fx);
        long shared[2] = {0, 0}, dx = 0, dsign = 0, dband = 0, band[2] = {0, 0}, inside = 0;
        for (int q = 0; q < 2 && !r; q++) {
            int s[3];
            for (int c = 0; c < 3; c++) s[c] = (int)lround((b[q].lo[c] - w.lo[c]) / hb);
            for (int k = 0; k < b[q].n[2]; k++)
                for (int j = 0; j < b[q].n[1]; j++)
                    for (int i = 0; i < b[q].n[0]; i++) {
                        const int I = i + s[0], J = j + s[1], K = k + s[2];
                        if (I < 0 || J < 0 || K < 0 || I >= w.n[0] || J >= w.n[1] || K >= w.n[2]) continue;
                        const float vw = at(&w, I, J, K), vb = at(&b[q], i, j, k);
                        shared[q]++, dx += at(&wx, I, J, K) != at(&bx[q], i, j, k), dsign += (vw < 0) != (vb < 0), inside += vb < 0;
                        if (fabs(at(&wx, I, J, K)) < 1.9 * hb) band[q]++, dband += vw != vb;
                    }
        }
        double ef = 0;
        for (int k = 0; !r && k < fx.n[2]; k++)
            for (int j = 0; j < fx.n[1]; j++)
                for (int i = 0; i < fx.n[0]; i++) {
                    double p[3];
                    node(&fx, i, j, k, p);
                    ef = fmax(ef, fabs(at(&fx, i, j, k) - sd_unit(p)));
                }
        const double sg = sag(&u, sd_unit);
        CHECK(r == 0 && shared[0] > 1000 && shared[0] < (long)b[0].n[0] * b[0].n[1] * b[0].n[2] && shared[1] > 100 && inside > 0 && dx == 0 && dsign == 0 && band[0] > 100 && band[1] > 50 &&
                  dband == 0 && ef < sg + 1e-6,
              "two boxes (%d x %d x %d, %d x %d x %d nodes) against the whole grid at %ld and %ld shared nodes (%ld inside): %ld differ exact, %ld signs and %ld of %ld band values plain; a box away from the mesh: error %.2e (sag %.2e)",
              b[0].n[0], b[0].n[1], b[0].n[2], b[1].n[0], b[1].n[1], b[1].n[2], shared[0], shared[1], inside, dx, dsign, dband, band[0] + band[1], ef, sg);
        if (!r) mk_sdf_free(&w), mk_sdf_free(&wx), mk_sdf_free(&b[0]), mk_sdf_free(&bx[0]), mk_sdf_free(&b[1]), mk_sdf_free(&bx[1]), mk_sdf_free(&fx);
        mk_sdf z;
        int k = 0, nr = 0;
        const double bad1[6] = {1, 0, 0, 0, 1, 1}, bad2[6] = {0, 0, 0, NAN, 1, 1};
#define REF(call) (k++, nr += (call) == -1 && z.d == NULL && z.n[0] == 0)
        REF(mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, box[0], MK_SDF_FILL, &z));
        REF(mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, bad1, 0, &z));
        REF(mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, bad2, 0, &z));
        REF(mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, hb, NULL, 0, &z));
        REF(mk_sdf_mesh_box(u.V, u.nv, u.T, u.nt, 0, box[0], 0, &z));
#undef REF
        CHECK(nr == k, "mk_sdf_mesh_box refuses %d of %d: MK_SDF_FILL, hi below lo, a NaN corner, no box, h 0", nr, k);
    }
    free(u.V), free(u.T);
    { /* refusals that leave the grid empty, and the sampling of an empty grid; an open triangle whose re-cast rays
         cross nothing */
        double W[9] = {1e10, 1e10, 1e10, 1e10 + 1e-6, 1e10, 1e10, 1e10, 1e10 + 1e-6, 1e10};
        const int T1[3] = {0, 1, 2};
        mk_sdf g;
        memset(&g, 0x5a, sizeof g);
        const int r = mk_sdf_mesh(W, 3, T1, 1, 0x1p-1000, 0, MK_SDF_LATTICE, &g);
        CHECK(r == -1 && g.n[0] == 0 && g.d == NULL, "a lattice grid past the double range (h 2^-1000 at 1e10) is refused, the grid left empty (%d)", r);
        double V2[9] = {0, 0, 0, 1, 0, 0, 0, 1, 1};
        memset(&g, 0x5a, sizeof g);
        const int r2 = mk_sdf_mesh(V2, 3, T1, 1, 1e-9, 0, 0, &g);
        CHECK(r2 == -1 && g.n[0] == 0 && g.n[1] == 0 && g.n[2] == 0 && g.h == 0 && g.d == NULL, "a grid too large for memory is refused with the grid left empty, all zero (%d; n %d %d %d)", r2, g.n[0], g.n[1], g.n[2]);
        const int Tb[3] = {0, 1, 3};
        const int r3 = mk_sdf_mesh(V2, 3, Tb, 1, 0.1, 0, 0, &g);
        CHECK(r3 == -1 && isnan(mk_sdf_sample(&g, (double[3]){0, 0, 0})), "sampling the empty grid a refusal leaves reads NaN");
        double V3[9] = {0.5, 0, 0, 1, 0.5, 1, 0.25, 1, 1};
        const int r4 = mk_sdf_mesh(V3, 3, T1, 1, 0.25 / 3, 0.5, 0, &g);
        CHECK(r4 == 0 && g.reray + g.odd_rows > 0, "one open triangle (a re-cast ray that crosses nothing): 0, %ld rows re-cast, %ld odd", r4 ? -1L : (long)g.reray, r4 ? -1L : (long)g.odd_rows);
        if (!r4) mk_sdf_free(&g);
    }
    { /* a NaN coordinate samples NaN, not the clamped value at node 0; +inf clamps to the grid's edge */
        float d8[8] = {-1, 2, 3, 4, 5, 6, 7, 8};
        mk_sdf g = {.n = {2, 2, 2}, .h = 1, .d = d8};
        const double a = mk_sdf_sample(&g, (double[3]){NAN, 0, 0}), b = mk_sdf_sample(&g, (double[3]){0, 0, NAN}), c = mk_sdf_sample(&g, (double[3]){INFINITY, 0, 0});
        CHECK(isnan(a) && isnan(b) && c == 2, "sampling at a NaN coordinate gives NaN (%g, %g); +inf still clamps to the box (%g)", a, b, c);
    }
    return wrong;
}
