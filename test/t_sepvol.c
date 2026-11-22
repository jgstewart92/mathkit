/* mk_sepvol: two balls sunk into each other in one mesh (single mode), a ball moved out of a held one, two cubes 1 mm
   into each other (sharp edges, where a push along the averaged vertex normal leaves faces crossing), a third part
   beside the crease, crumpled and coarse meshes held against thin parts, held meshes that are refused, allocation
   failures, the inside test's corner cases, the volume change, the refusals. Units are metres. Whether a vertex is
   inside another part is judged here by the winding number of that part alone, and crossing pairs by this file's own
   segment test, not by the header's. The meshes in test/data are read by paths relative to the repository root. With a
   folder argument the results are written there as OBJ files, which test/check.sh checks with meshcheck. */
[[maybe_unused]] static int t_pairs = 1; /* MK_SEPVOL_PAIRS is read at each call: row 12 turns step 2 off with it */
#define MK_SEPVOL_PAIRS t_pairs
/* row 13's allocation failures (test-only): the header's allocations are counted, the one numbered t_fi_at fails, and
   while one is set to fail each fresh malloc block is filled with 0xA5, so a block read before it is written shows */
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
static long t_fi_n, t_fi_at = -1;
static void *t_fi_malloc(size_t n)
{
    if (t_fi_at >= 0 && t_fi_n++ == t_fi_at) return NULL;
    void *p = malloc(n);
    if (p && t_fi_at >= 0) memset(p, 0xA5, n);
    return p;
}
[[maybe_unused]] static void *t_fi_calloc(size_t a, size_t b) { return t_fi_at >= 0 && t_fi_n++ == t_fi_at ? NULL : calloc(a, b); }
[[maybe_unused]] static void *t_fi_realloc(void *p, size_t n) { return t_fi_at >= 0 && t_fi_n++ == t_fi_at ? NULL : realloc(p, n); }
#define malloc(n) t_fi_malloc(n)
#define calloc(a, b) t_fi_calloc(a, b)
#define realloc(p, n) t_fi_realloc(p, n)
#include "../include/mk_sepvol.h"
#undef malloc
#undef calloc
#undef realloc
#include "t.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { double *V; int *T, nv, nt; } mesh;
static mesh readobj(const char *path)
{
    mesh m = {0};
    FILE *f = fopen(path, "r");
    if (!f) return m;
    char l[512];
    int cv = 1024, ct = 2048;
    m.V = malloc(sizeof(double) * 3 * cv), m.T = malloc(sizeof(int) * 3 * ct);
    while (fgets(l, sizeof l, f)) {
        if (l[0] == 'v' && l[1] == ' ') {
            if (m.nv == cv) m.V = realloc(m.V, sizeof(double) * 3 * (cv *= 2));
            sscanf(l + 2, "%lf %lf %lf", m.V + 3 * m.nv, m.V + 3 * m.nv + 1, m.V + 3 * m.nv + 2), m.nv++;
        } else if (l[0] == 'f' && l[1] == ' ') {
            if (m.nt == ct) m.T = realloc(m.T, sizeof(int) * 3 * (ct *= 2));
            char *p = l + 2;
            for (int k = 0; k < 3; k++) {
                m.T[3 * m.nt + k] = (int)strtol(p, &p, 10) - 1;
                while (*p && *p != ' ') p++;
            }
            m.nt++;
        }
    }
    fclose(f);
    return m;
}
static const char *odir;
static void writeobj(const char *name, const double *V, int nv, const int *T, int nt)
{
    if (!odir) return;
    char p[1024];
    snprintf(p, sizeof p, "%s/%s.obj", odir, name);
    FILE *f = fopen(p, "w");
    if (!f) return;
    for (int i = 0; i < nv; i++) fprintf(f, "v %.17g %.17g %.17g\n", V[3 * i], V[3 * i + 1], V[3 * i + 2]);
    for (int t = 0; t < nt; t++) fprintf(f, "f %d %d %d\n", T[3 * t] + 1, T[3 * t + 1] + 1, T[3 * t + 2] + 1);
    fclose(f);
}
/* parts: connected components through faces; comp[v], returns the count */
static int parts(const int *T, int nt, int nv, int *comp)
{
    int *u = malloc(sizeof *u * nv);
    for (int i = 0; i < nv; i++) u[i] = i;
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) {
            int a = T[3 * t + k], b = T[3 * t + (k + 1) % 3];
            while (u[a] != a) a = u[a];
            while (u[b] != b) b = u[b];
            u[a] = b;
        }
    int n = 0;
    for (int i = 0; i < nv; i++) comp[i] = -1;
    int *id = malloc(sizeof *id * nv);
    for (int i = 0; i < nv; i++) id[i] = -1;
    for (int i = 0; i < nv; i++) {
        int r = i;
        while (u[r] != r) r = u[r];
        if (id[r] < 0) id[r] = n++;
        comp[i] = id[r];
    }
    free(u), free(id);
    return n;
}
/* vertices of one part inside another (winding number of that part's faces over 1/2) */
static long cross_inside(const double *V, int nv, const int *T, int nt)
{
    int *comp = malloc(sizeof *comp * nv), np = parts(T, nt, nv, comp);
    long n = 0;
    int *TP = malloc(sizeof *TP * 3 * nt);
    for (int p = 0; p < np; p++) {
        int m = 0;
        for (int t = 0; t < nt; t++)
            if (comp[T[3 * t]] == p) memcpy(TP + 3 * m, T + 3 * t, sizeof *TP * 3), m++;
        for (int i = 0; i < nv; i++)
            if (comp[i] != p && mk_winding_number(V, TP, m, V + 3 * i) > 0.5) n++;
    }
    free(comp), free(TP);
    return n;
}
static long held_inside(const double *V, int nv, const mesh *H)
{
    long n = 0;
    for (int i = 0; i < nv; i++) n += mk_winding_number(H->V, H->T, H->nt, V + 3 * i) > 0.5;
    return n;
}
/* the largest outward share of a move: move . n0 / |move| over the moved vertices, n0 the area-weighted normal */
static double outward(const double *V0, const double *V, int nv, const int *T, int nt)
{
    double *N = calloc((size_t)3 * nv, sizeof *N), w = -1;
    for (int t = 0; t < nt; t++) {
        const double *a = V0 + 3 * T[3 * t], *b = V0 + 3 * T[3 * t + 1], *c = V0 + 3 * T[3 * t + 2];
        double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]};
        for (int k = 0; k < 3; k++) for (int j = 0; j < 3; j++) N[3 * T[3 * t + k] + j] += n[j];
    }
    for (int i = 0; i < nv; i++) {
        double m[3], lm = 0, ln = 0, d = 0;
        for (int k = 0; k < 3; k++) m[k] = V[3 * i + k] - V0[3 * i + k], lm += m[k] * m[k], ln += N[3 * i + k] * N[3 * i + k], d += m[k] * N[3 * i + k];
        if (lm > 0 && ln > 0) w = fmax(w, d / sqrt(lm * ln));
    }
    free(N);
    return w;
}
/* a cube of side s at o split n x n a face, outward, welded within itself, appended to m */
static void cube(mesh *m, double ox, double oy, double oz, double s, int n)
{
    const int base = m->nv, side = n + 1;
    int *id = malloc(sizeof *id * side * side * side);
    for (int i = 0; i < side * side * side; i++) id[i] = -1;
    m->V = realloc(m->V, sizeof(double) * 3 * (m->nv + 6 * side * side));
    m->T = realloc(m->T, sizeof(int) * 3 * (m->nt + 12 * n * n));
    for (int f = 0; f < 6; f++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                int q[4];
                const int di[4] = {0, 1, 1, 0}, dj[4] = {0, 0, 1, 1};
                for (int c = 0; c < 4; c++) {
                    const int I = i + di[c], J = j + dj[c], e = (f % 2) * n, ax = f / 2;
                    const int X = ax == 0 ? e : ax == 1 ? J : I, Y = ax == 0 ? I : ax == 1 ? e : J, Z = ax == 0 ? J : ax == 1 ? I : e;
                    int *p = &id[(Z * side + Y) * side + X];
                    if (*p < 0) {
                        *p = m->nv - base;
                        double *v = m->V + 3 * m->nv++;
                        v[0] = ox + s * X / n, v[1] = oy + s * Y / n, v[2] = oz + s * Z / n;
                    }
                    q[c] = base + *p;
                }
                int *t = m->T + 3 * m->nt;
                if (f % 2 == 0) t[0] = q[0], t[1] = q[2], t[2] = q[1], t[3] = q[0], t[4] = q[3], t[5] = q[2];
                else t[0] = q[0], t[1] = q[1], t[2] = q[2], t[3] = q[0], t[4] = q[2], t[5] = q[3];
                m->nt += 2;
            }
    free(id);
}
static void append(mesh *m, const mesh *a, const double *d)
{
    m->V = realloc(m->V, sizeof(double) * 3 * (m->nv + a->nv)), m->T = realloc(m->T, sizeof(int) * 3 * (m->nt + a->nt));
    for (int i = 0; i < a->nv; i++) for (int k = 0; k < 3; k++) m->V[3 * (m->nv + i) + k] = a->V[3 * i + k] + d[k];
    for (int t = 0; t < 3 * a->nt; t++) m->T[3 * m->nt + t] = a->T[t] + m->nv;
    m->nv += a->nv, m->nt += a->nt;
}

/* an edge pq through triangle abc's interior: Moller-Trumbore, strict (meshcheck's test), not the header's */
static int t_segtri(const double *p, const double *q, const double *a, const double *b, const double *c)
{
    double d[3], e1[3], e2[3], s[3], h[3], g[3];
    for (int k = 0; k < 3; k++) d[k] = q[k] - p[k], e1[k] = b[k] - a[k], e2[k] = c[k] - a[k], s[k] = p[k] - a[k];
    h[0] = d[1] * e2[2] - d[2] * e2[1], h[1] = d[2] * e2[0] - d[0] * e2[2], h[2] = d[0] * e2[1] - d[1] * e2[0];
    const double det = e1[0] * h[0] + e1[1] * h[1] + e1[2] * h[2];
    const double ld = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]), l1 = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]),
                 l2 = sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
    if (!(fabs(det) > 1e-15 * ld * l1 * l2)) return 0;
    const double f = 1 / det, u = f * (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]);
    if (u <= 0 || u >= 1) return 0;
    g[0] = s[1] * e1[2] - s[2] * e1[1], g[1] = s[2] * e1[0] - s[0] * e1[2], g[2] = s[0] * e1[1] - s[1] * e1[0];
    const double v = f * (d[0] * g[0] + d[1] * g[1] + d[2] * g[2]);
    if (v <= 0 || u + v >= 1) return 0;
    const double t = f * (e2[0] * g[0] + e2[1] * g[1] + e2[2] * g[2]);
    return t > 0 && t < 1;
}
/* a mesh's own crossing pairs (two faces sharing no corner, an edge of one through the other); with comp, only two
   faces of one part. mk_geom's tree only finds the candidates. */
static long own_pairs(const double *V, const int *T, int nt, const int *comp)
{
    mk_bvh B;
    if (mk_bvh_build(&B, V, T, nt)) return -1;
    int cap = 1024, *buf = malloc(sizeof *buf * cap);
    long n = 0;
    for (int t = 0; t < nt; t++) {
        const int *X = T + 3 * t;
        mk_aabb bx;
        mk_tri_box_(V, T, t, &bx);
        int m = mk_bvh_box(&B, &bx, buf, cap);
        if (m > cap) free(buf), buf = malloc(sizeof *buf * (cap = m)), m = mk_bvh_box(&B, &bx, buf, cap);
        for (int i = 0; i < m; i++) {
            const int u = buf[i], *Y = T + 3 * u;
            if (u <= t || (comp && comp[X[0]] != comp[Y[0]])) continue;
            int sh = 0;
            for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) sh |= X[a] == Y[b];
            if (sh) continue;
            int c = 0;
            for (int k = 0; k < 3 && !c; k++)
                c = t_segtri(V + 3 * X[k], V + 3 * X[(k + 1) % 3], V + 3 * Y[0], V + 3 * Y[1], V + 3 * Y[2]) ||
                    t_segtri(V + 3 * Y[k], V + 3 * Y[(k + 1) % 3], V + 3 * X[0], V + 3 * X[1], V + 3 * X[2]);
            n += c;
        }
    }
    free(buf), mk_bvh_free(&B);
    return n;
}
/* crossing pairs between two meshes, by the same test */
static long pair_cross(const mesh *a, const mesh *b)
{
    long n = 0;
    for (int t = 0; t < a->nt; t++)
        for (int u = 0; u < b->nt; u++) {
            const double *P[3] = {a->V + 3 * a->T[3 * t], a->V + 3 * a->T[3 * t + 1], a->V + 3 * a->T[3 * t + 2]};
            const double *Q[3] = {b->V + 3 * b->T[3 * u], b->V + 3 * b->T[3 * u + 1], b->V + 3 * b->T[3 * u + 2]};
            int c = 0;
            for (int k = 0; k < 3 && !c; k++) c = t_segtri(P[k], P[(k + 1) % 3], Q[0], Q[1], Q[2]) || t_segtri(Q[k], Q[(k + 1) % 3], P[0], P[1], P[2]);
            n += c;
        }
    return n;
}
/* the crumpled cube: 10 mm, n x n a face, each face's inner vertices moved along its normal by up to amp at random (a
   height field over each face, so it crosses nothing) */
static mesh crumpled(int n, double amp)
{
    mesh m = {0};
    cube(&m, 0, 0, 0, 0.01, n);
    for (int i = 0; i < m.nv; i++) {
        double *v = m.V + 3 * i;
        int on = -1, cnt = 0;
        for (int k = 0; k < 3; k++) if (fabs(v[k]) < 1e-12 || fabs(v[k] - 0.01) < 1e-12) on = k, cnt++;
        if (cnt == 1) v[on] += amp * t_rand();
    }
    return m;
}
/* a box bar th x len x wd (x, y, z), turned tilt about y, centred on the cube's +x face at y = z = 5 mm, its lowest
   corner depth into it */
static mesh barn(double th, double len, double wd, double tilt, double depth, int n) /* n x n a face */
{
    mesh m = {0};
    cube(&m, 0, 0, 0, 1, n);
    double xmin = INFINITY;
    for (int i = 0; i < m.nv; i++) {
        double *v = m.V + 3 * i, x = v[0] * th, y = (v[1] - 0.5) * len + 0.005, z = (v[2] - 0.5) * wd;
        v[0] = x * cos(tilt) - z * sin(tilt), v[1] = y, v[2] = x * sin(tilt) + z * cos(tilt) + 0.005, xmin = fmin(xmin, v[0]);
    }
    for (int i = 0; i < m.nv; i++) m.V[3 * i] += 0.01 - depth - xmin;
    return m;
}
static mesh bar(double th, double len, double wd, double tilt, double depth) { return barn(th, len, wd, tilt, depth, 4); }
/* a box lo .. hi (metres), n x n a face */
static mesh box(const double *lo, const double *hi, int n)
{
    mesh m = {0};
    cube(&m, 0, 0, 0, 1, n);
    for (int i = 0; i < m.nv; i++) for (int k = 0; k < 3; k++) m.V[3 * i + k] = lo[k] + m.V[3 * i + k] * (hi[k] - lo[k]);
    return m;
}
static double maxmove(const double *V0, const double *V, int nv)
{
    double w = 0;
    for (int i = 0; i < nv; i++) {
        const double *a = V + 3 * i, *b = V0 + 3 * i;
        w = fmax(w, sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
    }
    return w;
}
/* a latitude-longitude sphere of radius 1 about c, nu x nw, outward, appended to m */
static void lat_sphere(mesh *m, int nu, int nw, const double *c)
{
    const int base = m->nv, nv = 2 + (nw - 1) * nu, nt = 2 * nu + 2 * (nw - 2) * nu;
    m->V = realloc(m->V, sizeof(double) * 3 * (m->nv + nv)), m->T = realloc(m->T, sizeof(int) * 3 * (m->nt + nt));
    double *v = m->V + 3 * m->nv;
    v[0] = c[0], v[1] = c[1], v[2] = c[2] + 1, v += 3;
    for (int i = 1; i < nw; i++)
        for (int j = 0; j < nu; j++, v += 3) {
            const double th = T_PI * i / nw, ph = 2 * T_PI * j / nu;
            v[0] = c[0] + sin(th) * cos(ph), v[1] = c[1] + sin(th) * sin(ph), v[2] = c[2] + cos(th);
        }
    v[0] = c[0], v[1] = c[1], v[2] = c[2] - 1;
    m->nv += nv;
    int *t = m->T + 3 * m->nt;
    const int south = base + nv - 1;
    for (int j = 0; j < nu; j++, t += 3) t[0] = base, t[1] = base + 1 + j, t[2] = base + 1 + (j + 1) % nu;
    for (int i = 0; i < nw - 2; i++)
        for (int j = 0; j < nu; j++, t += 6) {
            const int a = base + 1 + i * nu + j, b = base + 1 + (i + 1) * nu + j, cc = base + 1 + (i + 1) * nu + (j + 1) % nu, d = base + 1 + i * nu + (j + 1) % nu;
            t[0] = a, t[1] = b, t[2] = cc, t[3] = a, t[4] = cc, t[5] = d;
        }
    for (int j = 0; j < nu; j++, t += 3) t[0] = base + 1 + (nw - 2) * nu + j, t[1] = south, t[2] = base + 1 + (nw - 2) * nu + (j + 1) % nu;
    m->nt += nt;
}
/* vertices inside another part, each counted once, by the winding number of the other parts' faces alone (no ray of the
   header's); *twice: those inside two others */
static long others_inside(const double *V, int nv, const int *T, int nt, long *twice)
{
    int *comp = malloc(sizeof *comp * nv), np = parts(T, nt, nv, comp), *TP = malloc(sizeof *TP * 3 * nt);
    long n = 0, n2 = 0;
    for (int i = 0; i < nv; i++) {
        int cov = 0;
        for (int p = 0; p < np; p++) {
            if (p == comp[i]) continue;
            int m = 0;
            for (int t = 0; t < nt; t++)
                if (comp[T[3 * t]] == p) memcpy(TP + 3 * m, T + 3 * t, sizeof *TP * 3), m++;
            cov += mk_winding_number(V, TP, m, V + 3 * i) > 0.5;
        }
        n += cov >= 1, n2 += cov >= 2;
    }
    free(comp), free(TP);
    *twice = n2;
    return n;
}
/* a regular tetrahedron, corners at scale (+-1, +-1, +-1) about (cx, 0, 0), outward or wound inward, appended to m */
static void tetra(mesh *m, double scale, double cx, int inward)
{
    static const double c[12] = {-1, -1, -1, 1, 1, -1, -1, 1, 1, 1, -1, 1};
    static const int f[12] = {0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3};
    const int base = m->nv;
    m->V = realloc(m->V, sizeof(double) * 3 * (m->nv + 4)), m->T = realloc(m->T, sizeof(int) * 3 * (m->nt + 4));
    for (int i = 0; i < 4; i++) m->V[3 * m->nv] = scale * c[3 * i] + cx, m->V[3 * m->nv + 1] = scale * c[3 * i + 1], m->V[3 * m->nv + 2] = scale * c[3 * i + 2], m->nv++;
    for (int t = 0; t < 4; t++, m->nt++)
        for (int k = 0; k < 3; k++) m->T[3 * m->nt + k] = base + f[3 * t + (inward ? 2 - k : k)];
}
/* the cube [-1, 1]^3, one square a face, each split along a diagonal drawn from the xorshift state *g */
static void diag_cube(mesh *m, unsigned long long *g)
{
    const int base = m->nv, n = 2;
    int id[8];
    for (int i = 0; i < 8; i++) id[i] = -1;
    m->V = realloc(m->V, sizeof(double) * 3 * (m->nv + 8)), m->T = realloc(m->T, sizeof(int) * 3 * (m->nt + 12));
    for (int f = 0; f < 6; f++) {
        const int ax = f / 2, sg = f & 1 ? 1 : -1, u = (ax + 1) % 3, v = (ax + 2) % 3;
        int q[4];
        for (int c = 0; c < 4; c++) {
            int I[3];
            I[ax] = sg > 0, I[u] = c == 1 || c == 2, I[v] = c >= 2;
            const int key = (I[0] * n + I[1]) * n + I[2];
            if (id[key] < 0) {
                id[key] = m->nv - base;
                double *p = m->V + 3 * m->nv++;
                p[0] = -1 + 2.0 * I[0], p[1] = -1 + 2.0 * I[1], p[2] = -1 + 2.0 * I[2];
            }
            q[c] = base + id[key];
        }
        *g ^= *g >> 12, *g ^= *g << 25, *g ^= *g >> 27;
        const int flip = (int)((*g * 2685821657736338717ULL) & 1);
        int *t = m->T + 3 * m->nt;
        m->nt += 2;
        if (sg > 0) {
            if (flip) t[0] = q[0], t[1] = q[1], t[2] = q[2], t[3] = q[0], t[4] = q[2], t[5] = q[3];
            else t[0] = q[0], t[1] = q[1], t[2] = q[3], t[3] = q[1], t[4] = q[2], t[5] = q[3];
        } else {
            if (flip) t[0] = q[0], t[1] = q[2], t[2] = q[1], t[3] = q[0], t[4] = q[3], t[5] = q[2];
            else t[0] = q[0], t[1] = q[3], t[2] = q[1], t[3] = q[1], t[4] = q[3], t[5] = q[2];
        }
    }
}
/* a box lo .. lo + s, one square a face, outward or wound inward, appended to m */
static void box1(mesh *m, double x, double y, double z, double s, int inward)
{
    const int t0 = m->nt;
    cube(m, x, y, z, s, 1);
    if (inward) for (int t = t0; t < m->nt; t++) { const int a = m->T[3 * t + 1]; m->T[3 * t + 1] = m->T[3 * t + 2], m->T[3 * t + 2] = a; }
}
/* the cosine between faces a and b's normals, positions V */
static double face_cos(const double *V, const int *a, const int *b)
{
    double n[2][3];
    for (int s = 0; s < 2; s++) {
        const int *q = s ? b : a;
        const double *p0 = V + 3 * q[0], *p1 = V + 3 * q[1], *p2 = V + 3 * q[2];
        const double u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]}, w[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
        n[s][0] = u[1] * w[2] - u[2] * w[1], n[s][1] = u[2] * w[0] - u[0] * w[2], n[s][2] = u[0] * w[1] - u[1] * w[0];
    }
    const double l = sqrt((n[0][0] * n[0][0] + n[0][1] * n[0][1] + n[0][2] * n[0][2]) * (n[1][0] * n[1][0] + n[1][1] * n[1][1] + n[1][2] * n[1][2]));
    return l > 0 ? (n[0][0] * n[1][0] + n[0][1] * n[1][1] + n[0][2] * n[1][2]) / l : 1;
}
/* over each pair of faces sharing an edge (brute force): the edges folded over 90 degrees at V, their count (what: 0)
   and length (what: 1); with V0, only those that met at under 90 degrees at V0 (at 90 or more there: left out) */
static double folded(const double *V0, const double *V, const int *T, int nt, int what)
{
    double s = 0;
    for (int t = 0; t < nt; t++)
        for (int u = t + 1; u < nt; u++) {
            int sh = 0, e[2] = {0, 0};
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    if (T[3 * t + a] == T[3 * u + b] && sh < 2) e[sh++] = T[3 * t + a];
            if (sh != 2 || !(face_cos(V, T + 3 * t, T + 3 * u) < 0) || (V0 && !(face_cos(V0, T + 3 * t, T + 3 * u) > 1e-9))) continue;
            const double *p = V + 3 * e[0], *q = V + 3 * e[1];
            s += what ? sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2])) : 1;
        }
    return s;
}
static long new_folds(const double *V0, const double *V, const int *T, int nt) { return (long)folded(V0, V, T, nt, 0); }
[[maybe_unused]] static double fold_length(const double *V, const int *T, int nt) { return folded(NULL, V, T, nt, 1); }
/* two meshes as one file, for meshcheck */
static void writepair(const char *name, const mesh *a, const mesh *b)
{
    mesh pr = {0};
    const double z[3] = {0, 0, 0};
    append(&pr, a, z), append(&pr, b, z);
    writeobj(name, pr.V, pr.nv, pr.T, pr.nt);
    free(pr.V), free(pr.T);
}

int main(int argc, char **argv)
{
    odir = argc > 1 ? argv[1] : NULL;
    const double gap = 5e-5, step = 2.5e-4; /* metres: 0.05 mm (the gap meshcheck's probe needs), 0.25 mm */
    const int rings = 3, rounds = 30;
    mk_sepvol_report r;

    /* 1. twoball, single mode */
    mesh tb = readobj("test/data/sepvol_twoball.obj");
    CHECK(tb.nv == 1284 && tb.nt == 2560, "twoball read: %d vertices, %d faces", tb.nv, tb.nt);
    double *tb0 = malloc(sizeof(double) * 3 * tb.nv);
    memcpy(tb0, tb.V, sizeof(double) * 3 * tb.nv);
    int *tbT0 = malloc(sizeof(int) * 3 * tb.nt);
    memcpy(tbT0, tb.T, sizeof(int) * 3 * tb.nt);
    const long tbin0 = cross_inside(tb.V, tb.nv, tb.T, tb.nt);
    int rc = mk_sepvol(tb.V, tb.nv, tb.T, tb.nt, NULL, 0, NULL, 0, gap, step, rings, rounds, &r);
    const long tbin1 = cross_inside(tb.V, tb.nv, tb.T, tb.nt);
    CHECK(rc == 0 && r.inside0 > 50 && r.inside == 0 && r.cross0 > 50 && r.cross == 0, "twoball: %ld vertices inside the other ball, %ld after %d rounds (winding: %ld, %ld); crossing face pairs %ld, %ld after", r.inside0, r.inside, r.rounds, tbin0, tbin1, r.cross0, r.cross);
    CHECK(tbin0 > 50 && tbin1 == 0, "twoball by the winding number of each ball: %ld inside before, %ld after", tbin0, tbin1);
    CHECK(!memcmp(tbT0, tb.T, sizeof(int) * 3 * tb.nt), "twoball: faces unchanged");
    const double vv = mk_sepvol_volume_(tb0, tb.T, tb.nt);
    CHECK(r.vol0 == vv && r.vol1 < r.vol0 && r.vol1 > 0.97 * r.vol0, "twoball volume %.4f ml -> %.4f ml (%.2f %%)", r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1));
    CHECK(r.worst <= rounds * step && r.worst > 0 && r.moved > 0, "twoball: %ld vertices moved, worst %.4f mm", r.moved, r.worst * 1e3);
    const double tbo = outward(tb0, tb.V, tb.nv, tb.T, tb.nt);
    CHECK(tbo < 0.01, "twoball: no move outward along a vertex's normal (largest share %.3f)", tbo);
    CHECK(r.status == 0 && r.self0 == 0 && r.self1 == 0, "twoball: status %d; own crossing pairs within one ball %ld, %ld after (the balls crossing each other is what is undone, not a self-crossing)", r.status, r.self0, r.self1);
    writeobj("twoball", tb.V, tb.nv, tb.T, tb.nt);

    /* 2. ballB moved out of a held ballA */
    mesh A = readobj("test/data/sepvol_ballA.obj"), Bm = readobj("test/data/sepvol_ballB.obj");
    CHECK(A.nv == 642 && Bm.nv == 642, "ballA, ballB read: %d, %d vertices", A.nv, Bm.nv);
    double *A0 = malloc(sizeof(double) * 3 * A.nv), *B0 = malloc(sizeof(double) * 3 * Bm.nv);
    int *AT0 = malloc(sizeof(int) * 3 * A.nt);
    memcpy(A0, A.V, sizeof(double) * 3 * A.nv), memcpy(AT0, A.T, sizeof(int) * 3 * A.nt), memcpy(B0, Bm.V, sizeof(double) * 3 * Bm.nv);
    const long hin0 = held_inside(Bm.V, Bm.nv, &A);
    rc = mk_sepvol(Bm.V, Bm.nv, Bm.T, Bm.nt, A.V, A.nv, A.T, A.nt, gap, step, rings, rounds, &r);
    const long hin1 = held_inside(Bm.V, Bm.nv, &A);
    CHECK(rc == 0 && r.inside0 > 20 && r.inside == 0 && r.rounds >= 2 && r.cross0 > 20 && r.cross == 0, "ballB out of held ballA: %ld inside, %ld after %d rounds; crossing face pairs %ld, %ld after", r.inside0, r.inside, r.rounds, r.cross0, r.cross);
    CHECK(hin0 > 20 && hin1 == 0, "ballB by ballA's winding number: %ld inside before, %ld after", hin0, hin1);
    CHECK(!memcmp(A0, A.V, sizeof(double) * 3 * A.nv) && !memcmp(AT0, A.T, sizeof(int) * 3 * A.nt), "ballA unchanged byte for byte");
    CHECK(r.status == 0 && !r.held_open && !r.held_flipped && !r.held_cross && r.held_vol > 0 && !r.self0 && !r.self1, "ballA taken as held: no open, flipped or crossing edges (%ld, %ld, %ld), volume %.4f ml; ballB's own crossing pairs %ld, %ld after", r.held_open, r.held_flipped, r.held_cross, r.held_vol * 1e6, r.self0, r.self1);
    CHECK(r.vol1 < r.vol0 && r.vol1 > 0.9 * r.vol0, "ballB volume %.4f ml -> %.4f ml (%.2f %%)", r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1));
    { /* faces, not only vertices, clear of ballA: each of ballB's face centres (where meshcheck probes) at least gap/4
         from it */
        double dmin = INFINITY;
        for (int t = 0; t < Bm.nt; t++) {
            double g[3] = {0, 0, 0}, w[3], q[3];
            for (int k = 0; k < 3; k++) for (int j = 0; j < 3; j++) g[j] += Bm.V[3 * Bm.T[3 * t + k] + j] / 3;
            if (mk_winding_number(A.V, A.T, A.nt, g) > 0.5) { dmin = -1; break; }
            for (int u = 0; u < A.nt; u++)
                dmin = fmin(dmin, mk_closest_pt_tri(g, A.V + 3 * A.T[3 * u], A.V + 3 * A.T[3 * u + 1], A.V + 3 * A.T[3 * u + 2], w, q));
        }
        dmin = dmin < 0 ? dmin : sqrt(dmin);
        CHECK(dmin >= 0.25 * gap, "ballB's face centres: the nearest %.4f mm from ballA (gap %.3f mm)", dmin * 1e3, gap * 1e3);
    }
    const double ho = outward(B0, Bm.V, Bm.nv, Bm.T, Bm.nt);
    CHECK(ho < 0.01, "ballB: no move outward along a vertex's normal (largest share %.3f)", ho);
    {
        mesh pr = {0};
        const double z[3] = {0, 0, 0};
        append(&pr, &A, z), append(&pr, &Bm, z);
        writeobj("pair", pr.V, pr.nv, pr.T, pr.nt);
        free(pr.V), free(pr.T);
    }

    /* 3. two cubes 10 mm, 8 x 8 a face, the second 1 mm into the first (sharp edges, where a push along the averaged
       vertex normal would leave faces crossing); the volume lost against the 1 mm x 8 mm x 8 mm overlap */
    mesh C = {0};
    cube(&C, 0, 0, 0, 0.01, 8), cube(&C, 0.009, 0.002, 0.002, 0.01, 8);
    writeobj("cubes_before", C.V, C.nv, C.T, C.nt);
    double *C0 = malloc(sizeof(double) * 3 * C.nv);
    memcpy(C0, C.V, sizeof(double) * 3 * C.nv);
    const long cin0 = cross_inside(C.V, C.nv, C.T, C.nt);
    const double ov = 1e-3 * 8e-3 * 8e-3; /* the overlap, m^3 */
    rc = mk_sepvol(C.V, C.nv, C.T, C.nt, NULL, 0, NULL, 0, gap, step, rings, rounds, &r);
    const long cin1 = cross_inside(C.V, C.nv, C.T, C.nt);
    CHECK(C.nv == 772 && C.nt == 1536, "cubes: %d vertices, %d faces", C.nv, C.nt);
    CHECK(rc == 0 && r.inside0 > 0 && r.inside == 0 && r.cross0 > 0 && r.cross == 0, "cubes: %ld vertices inside the other cube, %ld after %d rounds; crossing face pairs %ld, %ld after", r.inside0, r.inside, r.rounds, r.cross0, r.cross);
    CHECK(cin0 > 0 && cin1 == 0, "cubes by the winding number of each cube: %ld inside before, %ld after", cin0, cin1);
    CHECK(r.vol0 - r.vol1 >= ov && r.vol0 - r.vol1 < 1.5 * ov, "cubes volume %.4f ml -> %.4f ml: %.4f ml lost, %.2f times the 1 mm overlap (%.4f ml)", r.vol0 * 1e6, r.vol1 * 1e6, (r.vol0 - r.vol1) * 1e6, (r.vol0 - r.vol1) / ov, ov * 1e6);
    const double co = outward(C0, C.V, C.nv, C.T, C.nt);
    CHECK(co < 0.01, "cubes: no move outward along a vertex's normal (largest share %.3f)", co);
    writeobj("cubes", C.V, C.nv, C.T, C.nt);
    { /* one round from the start: each vertex inside the other cube moves the full step (no stall), and the vertices on
         the first cube's x = 10 mm edges go straight back along -x (edge-aware), not along their 45 degree normal. Each
         push is smoothed toward its pushed neighbours' (MK_SEPVOL_SMOOTH passes): within that many rings (1.25 mm each)
         of the second cube's side faces, where the pushes turn to -y or -z, the edge vertices' pushes turn part way too
         (up to 34 degrees); farther away they stay on -x. */
        memcpy(C.V, C0, sizeof(double) * 3 * C.nv);
        mk_sepvol(C.V, C.nv, C.T, C.nt, NULL, 0, NULL, 0, gap, step, rings, 1, &r);
        int ne = 0, nok = 0, nin = 0, nfull = 0, nf = 0, nfok = 0;
        double worst = 0, fworst = 0;
        const double band = 0.002 + (MK_SEPVOL_SMOOTH + 1.5) * 0.00125; /* beyond the smoothing's reach */
        for (int i = 0; i < C.nv; i++) {
            const double *p = C0 + 3 * i, m[3] = {C.V[3 * i] - p[0], C.V[3 * i + 1] - p[1], C.V[3 * i + 2] - p[2]};
            const double l = sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
            const int inB = p[0] > 0.009 - 1e-9 && p[1] > 0.002 && p[1] < 0.012 && p[2] > 0.002 && p[2] < 0.012 && p[0] < 0.0100001;
            if (i < 386 && inB) nin++, nfull += l > 0.99 * step;
            if (i < 386 && fabs(p[0] - 0.01) < 1e-12 && (fabs(p[1] - 0.01) < 1e-12 || fabs(p[2] - 0.01) < 1e-12) && p[1] > 0.0025 && p[2] > 0.0025) {
                ne++;
                const double ang = acos(fmin(1, -m[0] / l)) * 180 / T_PI;
                worst = fmax(worst, ang), nok += ang < 40;
                if (p[1] > band && p[2] > band) nf++, fworst = fmax(fworst, ang), nfok += ang < 10;
            }
        }
        CHECK(nin > 40 && nfull == nin, "cubes, one round: %d of the first cube's %d vertices inside the second moved the full step", nfull, nin);
        CHECK(ne >= 10 && nok == ne, "cubes, one round: %d of %d vertices on the x = 10 mm edges moved nearer -x than their 45 deg normal (worst %.1f deg off -x)", nok, ne, worst);
        CHECK(nf >= 5 && nfok == nf, "cubes, one round: %d of %d of them past %.2f mm from the second cube's side faces moved along -x (worst %.1f deg off)", nfok, nf, band * 1e3, fworst);
    }

    /* 4. a third part beside the crease of two overlapping balls: nothing moves into it, and it does not move */
    {
        mesh T3 = {0};
        const double z[3] = {0, 0, 0};
        mesh a = readobj("test/data/sepvol_ballA.obj");
        double lo[3] = {1e9, 1e9, 1e9}, hi[3] = {-1e9, -1e9, -1e9}, ctr[3];
        for (int i = 0; i < a.nv; i++) for (int k = 0; k < 3; k++) lo[k] = fmin(lo[k], a.V[3 * i + k]), hi[k] = fmax(hi[k], a.V[3 * i + k]);
        for (int k = 0; k < 3; k++) ctr[k] = 0.5 * (lo[k] + hi[k]);
        const double R = 0.5 * (hi[0] - lo[0]), Ry = 0.5 * (hi[1] - lo[1]), sh[3] = {1.4 * R, 0, 0};
        /* ballA is an ellipsoid (semi-axes R, Ry, R: 5, 7, 5 mm). The crease of A and A + sh lies at x = ctr + 0.7 R,
           out to y = 0.714 Ry; A scaled by 0.3 and centred at (0.7 R, 1.14 Ry) sits beside it, outside both (its
           nearest point reads 1.03 to 1.05 in either's ellipsoid norm), 0.13 Ry above the crease */
        const double rs = 0.3, sh3[3] = {0.7 * R, 1.14 * Ry, 0};
        mesh s = {0};
        append(&s, &a, z);
        for (int i = 0; i < s.nv; i++) for (int k = 0; k < 3; k++) s.V[3 * i + k] = ctr[k] + rs * (s.V[3 * i + k] - ctr[k]);
        append(&T3, &a, z), append(&T3, &a, sh), append(&T3, &s, sh3);
        const int third = 2 * a.nv;
        double *T30 = malloc(sizeof(double) * 3 * T3.nv);
        memcpy(T30, T3.V, sizeof(double) * 3 * T3.nv);
        const long in3a = cross_inside(T3.V, T3.nv, T3.T, T3.nt);
        rc = mk_sepvol(T3.V, T3.nv, T3.T, T3.nt, NULL, 0, NULL, 0, gap, step, rings, rounds, &r);
        const long in3b = cross_inside(T3.V, T3.nv, T3.T, T3.nt);
        CHECK(rc == 0 && in3a > 0 && in3b == 0 && r.inside == 0, "three parts: %ld vertices inside another part before, %ld after (%d rounds)", in3a, in3b, r.rounds);
        CHECK(!memcmp(T30 + 3 * third, T3.V + 3 * third, sizeof(double) * 3 * (T3.nv - third)), "three parts: the third (beside the crease) unchanged");
        writeobj("three", T3.V, T3.nv, T3.T, T3.nt);
        free(T30), free(T3.V), free(T3.T), free(s.V), free(s.T), free(a.V), free(a.T);
    }

    /* 6. a dense surface crumpled at its edges' scale, held against a sharp part: moves of half an edge can fold it
       into own crossing pairs between faces one edge apart, which the smoothing and the self-crossing guard (steps 3
       and 4) prevent. Here a crumpled cube (10 mm, 30 x 30 a face, its faces' inner vertices up to 0.1 mm off along the
       face normal) is held against a 1 mm square bar turned 45 degrees, its edge 0.4 mm in (gap 0.1 mm, step 0.25 mm,
       4 rings, 30 rounds). */
    {
        mesh K = crumpled(30, 1e-4), R = bar(0.001, 0.006, 0.001, T_PI / 4, 4e-4);
        double *R0 = malloc(sizeof(double) * 3 * R.nv);
        memcpy(R0, R.V, sizeof(double) * 3 * R.nv);
        const long k0 = own_pairs(K.V, K.T, K.nt, NULL), kr0 = pair_cross(&K, &R);
        { /* for meshcheck's negative control: the cube and the bar as they come in */
            mesh pr = {0};
            const double z[3] = {0, 0, 0};
            append(&pr, &K, z), append(&pr, &R, z);
            writeobj("crumpled_before", pr.V, pr.nv, pr.T, pr.nt);
            free(pr.V), free(pr.T);
        }
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, 1e-4, 2.5e-4, 4, 30, &r);
        const long k1 = own_pairs(K.V, K.T, K.nt, NULL), kr1 = pair_cross(&K, &R), kin1 = held_inside(K.V, K.nv, &R);
        CHECK(K.nv == 5402 && K.nt == 10800 && k0 == 0 && kr0 > 100, "crumpled cube: %d vertices, %d faces, %ld own crossing pairs; %ld crossing pairs with the bar", K.nv, K.nt, k0, kr0);
        CHECK(rc == 0 && k1 == 0, "crumpled cube against the bar: own crossing pairs %ld -> %ld (the test's count)", k0, k1);
        CHECK(rc == 0 && r.inside == 0 && r.cross == 0 && kin1 == 0 && kr1 == 0, "crumpled cube out of the bar: %ld vertices inside, %ld after %d rounds (winding: %ld); crossing pairs %ld, %ld after (the test's count: %ld)", r.inside0, r.inside, r.rounds, kin1, r.cross0, r.cross, kr1);
        CHECK(!memcmp(R0, R.V, sizeof(double) * 3 * R.nv), "the bar unchanged byte for byte");
        CHECK(r.vol1 < r.vol0 && r.vol1 > 0.99 * r.vol0 && r.worst < 0.002, "crumpled cube volume %.4f ml -> %.4f ml (%.2f %%), worst move %.3f mm", r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1), r.worst * 1e3);
        CHECK(r.status == 0 && r.self0 == k0 && r.self1 == k1, "crumpled cube: the report's own crossing pairs %ld, %ld after (the test's: %ld, %ld); %ld moves halved and kept, %ld refused", r.self0, r.self1, k0, k1, r.shortened, r.refused);
        CHECK(r.rounds <= 10 && r.refused <= 20, "crumpled cube: separated in %d rounds with %ld moves refused (at most 10 and 20)", r.rounds, r.refused);
        { /* for meshcheck: the cube and the bar as one file */
            mesh pr = {0};
            const double z[3] = {0, 0, 0};
            append(&pr, &K, z), append(&pr, &R, z);
            writeobj("crumpled", pr.V, pr.nv, pr.T, pr.nt);
            free(pr.V), free(pr.T);
        }
        free(K.V), free(K.T), free(R.V), free(R.T), free(R0);
    }

    /* 6b. deeper: a 0.6 mm bar, its edge 0.8 mm into a cube crumpled up to 0.15 mm. Denting that far folds the crumples
       even with the pushes smoothed; the guard halves or refuses those moves, and the run never hands back a state
       worse than the input. The cube's vertices inside the bar are past its middle
       (the bar is thinner than the dent), so their pushes point outward and step 3 would take them away; the held exit
       sends them back the way they came, and the cube comes clear. */
    {
        mesh K = crumpled(30, 1.5e-4), R = bar(0.0006, 0.006, 0.0006, T_PI / 4, 8e-4);
        double *R0 = malloc(sizeof(double) * 3 * R.nv);
        memcpy(R0, R.V, sizeof(double) * 3 * R.nv);
        const long k0 = own_pairs(K.V, K.T, K.nt, NULL), kr0 = pair_cross(&K, &R), kin0 = held_inside(K.V, K.nv, &R);
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, 1e-4, 2.5e-4, 4, 30, &r);
        const long k1 = own_pairs(K.V, K.T, K.nt, NULL), kr1 = pair_cross(&K, &R), kin1 = held_inside(K.V, K.nv, &R);
        CHECK(rc == 0 && k0 == 0 && k1 == 0, "crumpled cube, bar 0.8 mm in: own crossing pairs %ld -> %ld (the test's count)", k0, k1);
        CHECK(rc == 0 && kr0 > 100 && kr1 <= kr0 && kin1 <= kin0, "crumpled cube, bar 0.8 mm in: no worse than it came in: crossing pairs with the bar %ld -> %ld, vertices inside it %ld -> %ld (winding), %d rounds", kr0, kr1, kin0, kin1, r.rounds);
        CHECK(!memcmp(R0, R.V, sizeof(double) * 3 * R.nv) && r.vol1 > 0.98 * r.vol0, "crumpled cube, bar 0.8 mm in: the bar unchanged byte for byte; volume %.4f ml -> %.4f ml (%.2f %%)", r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1));
        CHECK(r.status == 0 && r.self0 == 0 && r.self1 == 0 && r.inside <= r.inside0 && r.cross <= r.cross0, "crumpled cube, bar 0.8 mm in: the report's own crossing pairs %ld, %ld after; %ld moves halved and kept, %ld refused; inside %ld -> %ld, crossing pairs %ld -> %ld", r.self0, r.self1, r.shortened, r.refused, r.inside0, r.inside, r.cross0, r.cross);
        CHECK(rc == 0 && kr1 == 0 && kin1 == 0 && r.inside == 0 && r.cross == 0, "crumpled cube, bar 0.8 mm in, clear of it by the held exit: crossing pairs %ld -> %ld, vertices inside %ld -> %ld (the report: %ld, %ld)", kr0, kr1, kin0, kin1, r.inside, r.cross);
        free(K.V), free(K.T), free(R.V), free(R.T), free(R0);
    }

    /* 7. a held mesh that cannot say inside from outside (open, crossing itself, inside out, one face turned) is
       refused: -2, V as it came in, the report's status and held_ fields saying why. */
    {
        mesh a = readobj("test/data/sepvol_ballA.obj"), b = readobj("test/data/sepvol_ballB.obj");
        double *b0 = malloc(sizeof(double) * 3 * b.nv);
        memcpy(b0, b.V, sizeof(double) * 3 * b.nv);
        mesh u = {0}; /* two overlapping copies of ballA, 2 mm apart */
        const double z[3] = {0, 0, 0}, sh[3] = {0.002, 0, 0};
        append(&u, &a, z), append(&u, &a, sh);
        int *flip = malloc(sizeof(int) * 3 * a.nt), *one = malloc(sizeof(int) * 3 * a.nt);
        memcpy(one, a.T, sizeof(int) * 3 * a.nt);
        for (int t = 0; t < a.nt; t++) flip[3 * t] = a.T[3 * t], flip[3 * t + 1] = a.T[3 * t + 2], flip[3 * t + 2] = a.T[3 * t + 1];
        one[1] = a.T[2], one[2] = a.T[1];
        struct { const char *what; const double *V; int nv; const int *T; int nt; int bit; } cs[] = {
            {"open (ballA without its last 3 faces)", a.V, a.nv, a.T, a.nt - 3, 1},
            {"crossing itself (two overlapping copies of ballA)", u.V, u.nv, u.T, u.nt, 4},
            {"inside out (ballA, every face turned)", a.V, a.nv, flip, a.nt, 2},
            {"one face turned", a.V, a.nv, one, a.nt, 2},
        };
        for (int c = 0; c < 4; c++) {
            rc = mk_sepvol(b.V, b.nv, b.T, b.nt, cs[c].V, cs[c].nv, cs[c].T, cs[c].nt, gap, step, rings, rounds, &r);
            CHECK(rc == -2 && !memcmp(b.V, b0, sizeof(double) * 3 * b.nv), "held mesh %s: refused (%d), ballB as it came in", cs[c].what, rc);
            /* each case's open, flipped and crossing counts: 0 or not, and 3 exactly */
            const long want[4][3] = {{1, 0, 0}, {0, 0, 1}, {0, 0, 0}, {0, 3, 0}};
            const int ok = r.status == cs[c].bit && (!want[c][0] == !r.held_open) && (want[c][1] == 3 ? r.held_flipped == 3 : r.held_flipped == 0) &&
                           (!want[c][2] == !r.held_cross) && r.inside0 == 0 && r.moved == 0 && (c != 2 || r.held_vol < 0);
            CHECK(ok, "held mesh %s: status %d; open edges %ld, edges run one way %ld, crossing pairs %ld, volume %.4f ml", cs[c].what, r.status, r.held_open, r.held_flipped, r.held_cross, r.held_vol * 1e6);
        }
        free(a.V), free(a.T), free(b.V), free(b.T), free(b0), free(u.V), free(u.T), free(flip), free(one);
    }

    /* 8. single mode on a part that already crosses itself: twoball, one vertex of the first ball (on its side away
       from the second) pushed through that ball's far side. The part's own crossings are reported and never made more;
       the two balls still come apart (a tenth of the vertices inside, or fewer, are left). */
    {
        mesh s = readobj("test/data/sepvol_twoball.obj");
        int *comp = malloc(sizeof *comp * s.nv);
        parts(s.T, s.nt, s.nv, comp);
        double c0[3] = {0, 0, 0}, c1[3] = {0, 0, 0};
        int n0 = 0, n1 = 0;
        for (int i = 0; i < s.nv; i++) {
            double *c = comp[i] == comp[0] ? c0 : c1;
            for (int k = 0; k < 3; k++) c[k] += s.V[3 * i + k];
            if (comp[i] == comp[0]) n0++; else n1++;
        }
        for (int k = 0; k < 3; k++) c0[k] /= n0, c1[k] /= n1;
        int far = -1;
        double best = -INFINITY;
        for (int i = 0; i < s.nv; i++) { /* the first ball's vertex farthest from the second's centre */
            if (comp[i] != comp[0]) continue;
            const double d = (s.V[3 * i] - c1[0]) * (s.V[3 * i] - c1[0]) + (s.V[3 * i + 1] - c1[1]) * (s.V[3 * i + 1] - c1[1]) + (s.V[3 * i + 2] - c1[2]) * (s.V[3 * i + 2] - c1[2]);
            if (d > best) best = d, far = i;
        }
        /* through the centre and out, 0.3 radius past the far side */
        for (int k = 0; k < 3; k++) s.V[3 * far + k] = c0[k] + 1.3 * (c0[k] - s.V[3 * far + k]);
        const long s0 = own_pairs(s.V, s.T, s.nt, comp), in0 = cross_inside(s.V, s.nv, s.T, s.nt);
        rc = mk_sepvol(s.V, s.nv, s.T, s.nt, NULL, 0, NULL, 0, gap, step, rings, rounds, &r);
        const long s1 = own_pairs(s.V, s.T, s.nt, comp), in1 = cross_inside(s.V, s.nv, s.T, s.nt);
        CHECK(rc == 0 && s0 > 0 && s1 <= s0, "a part crossing itself, single mode: its own crossing pairs %ld -> %ld (the test's count, one part's faces)", s0, s1);
        CHECK((r.status & MK_SEPVOL_SELF_CROSSES) && r.self0 == s0 && r.self1 == s1, "a part crossing itself: status %d (self-crossing reported), the report's own pairs %ld, %ld after (the test's: %ld, %ld)", r.status, r.self0, r.self1, s0, s1);
        CHECK(in0 > 50 && in1 <= in0 / 10, "a part crossing itself beside the second ball: vertices of one ball inside the other %ld -> %ld (the winding number)", in0, in1);
        free(s.V), free(s.T), free(comp);
    }

    /* 9. a runaway to avoid: a coarse mesh held against a part thinner than its faces, with a large step. Step 2 can
       push corners that lie behind one face of the thin part through it, round after round, until the mesh turns inside
       out; the self-crossing guard and the best state bound that, and the held exit (step 1) clears it. Here a cube
       10 mm, 4 x 4 a face, against a plate 0.4 mm thick, 12 mm square, 24 x 24 a face, turned 30 degrees, its edge
       1.5 mm into the cube; gap 0.2 mm, step 3 mm, 3 rings, 30 rounds. */
    {
        mesh K = {0}, R = barn(4e-4, 0.012, 0.012, T_PI / 6, 1.5e-3, 24);
        cube(&K, 0, 0, 0, 0.01, 4);
        double *K0 = malloc(sizeof(double) * 3 * K.nv);
        memcpy(K0, K.V, sizeof(double) * 3 * K.nv);
        const long k0 = own_pairs(K.V, K.T, K.nt, NULL), kr0 = pair_cross(&K, &R), kin0 = held_inside(K.V, K.nv, &R);
        writepair("coarse_before", &K, &R);
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, 2e-4, 3e-3, 3, 30, &r);
        const long k1 = own_pairs(K.V, K.T, K.nt, NULL), kr1 = pair_cross(&K, &R), kin1 = held_inside(K.V, K.nv, &R);
        const double w = maxmove(K0, K.V, K.nv);
        CHECK(K.nv == 98 && R.nv == 3458 && k0 == 0 && kr0 > 100, "coarse cube (%d vertices) against a thin plate (%d): %ld own crossing pairs, %ld with the plate, %ld vertices inside it", K.nv, R.nv, k0, kr0, kin0);
        CHECK(rc == 0 && k1 == 0 && r.vol1 > 0.8 * r.vol0 && w <= 3e-3 && kr1 <= kr0, "coarse cube against a thin plate, no runaway: own crossing pairs %ld -> %ld, with the plate %ld -> %ld, volume %.4f -> %.4f ml (%.2f %%), worst move %.3f mm (step 3 mm)", k0, k1, kr0, kr1, r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1), w * 1e3);
        CHECK(rc == 0 && kr1 == 0 && kin1 == 0 && r.inside == 0 && r.cross == 0, "coarse cube clear of the thin plate by the held exit: crossing pairs %ld -> %ld, vertices inside %ld -> %ld (the report: %ld, %ld), %d rounds", kr0, kr1, kin0, kin1, r.inside, r.cross, r.rounds);
        writepair("coarse", &K, &R);
        free(K.V), free(K.T), free(R.V), free(R.T), free(K0);
    }

    /* 10. past the middle of a thin held part: a vertex whose nearest exit lies on the part's far side has a push that
       points outward, which step 3 takes away, so it would slide sideways and stall. A cube 10 mm, 20 x 20 a face, its
       +x face 0.7 mm into a plate 1 mm thick and 5.8 mm square (16 x 16 a face, off the cube's grid lines), the test's
       settings. The vertices leave by the side they came in: the volume lost is one to two times the overlap, and no
       move is as long as twice the depth plus the gap. */
    {
        mesh K = {0};
        cube(&K, 0, 0, 0, 0.01, 20);
        const double lo[3] = {0.0093, 0.0021, 0.0021}, hi[3] = {0.0103, 0.0079, 0.0079}, ov = 0.7e-3 * 5.8e-3 * 5.8e-3;
        mesh R = box(lo, hi, 16);
        double *K0 = malloc(sizeof(double) * 3 * K.nv), *R0 = malloc(sizeof(double) * 3 * R.nv);
        memcpy(K0, K.V, sizeof(double) * 3 * K.nv), memcpy(R0, R.V, sizeof(double) * 3 * R.nv);
        const long kin0 = held_inside(K.V, K.nv, &R), kr0 = pair_cross(&K, &R);
        writepair("plate_before", &K, &R);
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, gap, step, rings, rounds, &r);
        const long k1 = own_pairs(K.V, K.T, K.nt, NULL), kr1 = pair_cross(&K, &R), kin1 = held_inside(K.V, K.nv, &R);
        const double w = maxmove(K0, K.V, K.nv);
        CHECK(rc == 0 && kin0 > 100 && kin1 == 0 && kr1 == 0 && k1 == 0 && r.inside == 0 && r.cross == 0, "cube 0.7 mm into a 1 mm plate (past its middle): vertices inside %ld -> %ld, crossing pairs with it %ld -> %ld, own %ld; %d rounds", kin0, kin1, kr0, kr1, k1, r.rounds);
        CHECK(rc == 0 && r.vol0 - r.vol1 >= ov && r.vol0 - r.vol1 < 2 * ov && w < 2 * (0.7e-3 + gap) && !memcmp(R0, R.V, sizeof(double) * 3 * R.nv), "cube out of the plate the way it came: volume %.4f -> %.4f ml, %.2f times the overlap; worst move %.3f mm (the depth 0.7 mm); the plate unchanged", r.vol0 * 1e6, r.vol1 * 1e6, (r.vol0 - r.vol1) / ov, w * 1e3);
        writepair("plate", &K, &R);
        free(K.V), free(K.T), free(R.V), free(R.T), free(K0), free(R0);
    }

    /* 11. where going back is wrong: a bar 4.4 mm wide through a slab 3 mm thick (20 mm square, 20 x 20 a face), its
       top 0.5 mm above the slab's. The slab's top vertices inside the bar have their nearest exit ahead (the bar's
       top), but the bar runs on through the slab: along their inward normal the slab's own bottom comes before the
       bar's far side, so their nearest exit stays. Without that check they would be pushed down through the slab; the
       row allows a volume loss under 15 %. The slab's vertices are moved up to 0.1 um at random, so no two of its faces
       lie exactly in one plane: a side folded flat onto itself in its own plane reads as crossing to this file's
       segment test (rounding on a face's edge), not to the header's. */
    {
        const double lo[3] = {0, 0, 0}, hi[3] = {0.02, 0.02, 0.003}, blo[3] = {0.0079, -0.005, -0.008}, bhi[3] = {0.0123, 0.025, 0.0035};
        mesh S = box(lo, hi, 20), R = box(blo, bhi, 8);
        for (int i = 0; i < 3 * S.nv; i++) S.V[i] += 1e-7 * t_rand();
        const long sin0 = held_inside(S.V, S.nv, &R);
        rc = mk_sepvol(S.V, S.nv, S.T, S.nt, R.V, R.nv, R.T, R.nt, gap, step, rings, rounds, &r);
        const long s1 = own_pairs(S.V, S.T, S.nt, NULL), sin1 = held_inside(S.V, S.nv, &R);
        CHECK(rc == 0 && sin0 > 100 && s1 == 0 && r.vol1 > 0.85 * r.vol0 && sin1 <= sin0, "bar through a slab: the slab's volume %.4f -> %.4f ml (%.2f %%), own crossing pairs %ld, vertices inside the bar %ld -> %ld", r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1), s1, sin0, sin1);
        free(S.V), free(S.T), free(R.V), free(R.T);
    }

    /* 12. MK_SEPVOL_PAIRS: the plate of row 9, square to x and 0.4 mm thick, runs right through the coarse cube 1.5 mm
       in from its face, with no vertex of the cube inside it. No vertex move can clear that; step 2 squeezes the cube
       round it (about a third of its volume). With the switch off nothing moves; the crossing pairs are reported. */
    {
        mesh K = {0}, R = barn(4e-4, 0.012, 0.012, 0, 1.5e-3, 24);
        cube(&K, 0, 0, 0, 0.01, 4);
        double *K0 = malloc(sizeof(double) * 3 * K.nv);
        memcpy(K0, K.V, sizeof(double) * 3 * K.nv);
        t_pairs = 0;
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, 2e-4, 3e-3, 3, 30, &r);
        t_pairs = 1;
        CHECK(rc == 0 && r.inside0 == 0 && r.cross0 > 100 && r.cross == r.cross0 && !memcmp(K0, K.V, sizeof(double) * 3 * K.nv), "a plate through the coarse cube, step 2 off: crossing pairs %ld -> %ld, %ld vertices inside, the cube unchanged byte for byte", r.cross0, r.cross, r.inside0);
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, 2e-4, 3e-3, 3, 30, &r);
        CHECK(rc == 0 && r.cross < r.cross0 && r.vol1 < 0.9 * r.vol0, "the same, step 2 on: crossing pairs %ld -> %ld, the cube squeezed round the plate, %.4f -> %.4f ml (%.2f %%)", r.cross0, r.cross, r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1));
        free(K.V), free(K.T), free(R.V), free(R.T), free(K0);
    }

    /* 13. no memory: -1 with V as it came in and the report zeroed. Two cubes 2 x 2 a face, 1 mm into each other,
       single mode, 3 rounds; each of the call's allocations made to fail in turn, from the first (before the saved copy
       of V is filled, so copying that block back over V would show as 0xA5 bytes) to those of the last round (after V
       has moved). */
    {
        mesh m = {0};
        cube(&m, 0, 0, 0, 0.01, 2), cube(&m, 0.009, 0.002, 0.002, 0.01, 2);
        double *m0 = malloc(sizeof(double) * 3 * m.nv);
        memcpy(m0, m.V, sizeof(double) * 3 * m.nv);
        const mk_sepvol_report zero = {0};
        t_fi_at = LONG_MAX, t_fi_n = 0; /* counted, none failing */
        rc = mk_sepvol(m.V, m.nv, m.T, m.nt, NULL, 0, NULL, 0, gap, step, rings, 3, &r);
        const long nalloc = t_fi_n;
        const int rounds0 = r.rounds, moved = r.moved > 0 && r.rounds >= 2;
        const long moved0 = r.moved;
        long nok = 0, nm1 = 0;
        for (long k = 0; k < nalloc; k++) {
            memcpy(m.V, m0, sizeof(double) * 3 * m.nv);
            memset(&r, 0x5A, sizeof r);
            t_fi_at = k, t_fi_n = 0;
            const int rk = mk_sepvol(m.V, m.nv, m.T, m.nt, NULL, 0, NULL, 0, gap, step, rings, 3, &r);
            nm1 += rk == -1;
            nok += rk == -1 && !memcmp(m.V, m0, sizeof(double) * 3 * m.nv) && !memcmp(&r, &zero, sizeof r);
        }
        t_fi_at = -1;
        CHECK(rc == 0 && moved && nalloc > 20 && nok == nalloc, "allocation failures: each of the call's %ld allocations failing in turn: %ld return -1 (no memory), %ld of them with V exactly as it came in and the report zeroed (without failures: %d rounds, %ld vertices moved)",
              nalloc, nm1, nok, rounds0, moved0);
        free(m.V), free(m.T), free(m0);
    }
    /* 14. a part thinner than eps (a thousandth of the mean edge) is inside nothing. A cube 0.2 mm across, 4 m from a
       cube 2 m across (eps 0.96 mm): the point eps inside each of its vertices lies outside it, outside every part,
       where even parity would read "inside another part". The winding number (0, not 2) decides: nothing is inside,
       nothing moves. */
    {
        mesh m = {0};
        cube(&m, -1, -1, -1, 2, 2), cube(&m, 5 - 1e-4, -1e-4, -1e-4, 2e-4, 1);
        double *m0 = malloc(sizeof(double) * 3 * m.nv);
        memcpy(m0, m.V, sizeof(double) * 3 * m.nv);
        rc = mk_sepvol(m.V, m.nv, m.T, m.nt, NULL, 0, NULL, 0, 0, 0.1, 0, 3, &r);
        CHECK(rc == 0 && r.inside0 == 0 && r.moved == 0 && r.cross0 == 0 && !memcmp(m.V, m0, sizeof(double) * 3 * m.nv), "a cube 0.2 mm across beside one 2 m across: %ld vertices inside another part, %ld moved, %d rounds",
              r.inside0, r.moved, r.rounds);
        free(m.V), free(m.T), free(m0);
    }
    /* 15. row 14's winding number made cheap: on a closed mesh mk_sepvol_two_ reads the clean rays' facings summed
       (exactly the winding number there) where mk_winding_number costs every triangle at every point. At the point eps
       inside every vertex of twoball, the two cubes and row 14's pair, the decision on the facings (closed 1) is the
       one on the winding number (closed 0); the even-parity points include ones covered twice and by nothing */
    {
        /* mk_sepvol_two_'s three ray directions */
        static const double dirs[3][3] = {{0.5773502691896258, 0.5773502691896258, 0.5773502691896258},
                                          {-0.3015113445777636, 0.9045340337332909, 0.3015113445777636},
                                          {0.2672612419124244, -0.5345224838248488, 0.8017837257372732}};
        mesh ms[3] = {readobj("test/data/sepvol_twoball.obj"), {0}, {0}};
        cube(&ms[1], 0, 0, 0, 0.01, 8), cube(&ms[1], 0.009, 0.002, 0.002, 0.01, 8);
        cube(&ms[2], -1, -1, -1, 2, 2), cube(&ms[2], 5 - 1e-4, -1e-4, -1e-4, 2e-4, 1);
        long npt = 0, dis = 0, nin = 0, e0 = 0, e2 = 0, closed = 0;
        for (int c = 0; c < 3; c++) {
            mesh *m = &ms[c];
            mk_bvh B;
            double *N = malloc(sizeof(double) * 3 * m->nv), el = 0;
            if (mk_bvh_build(&B, m->V, m->T, m->nt)) { dis++; continue; }
            mk_sepvol_normals_(m->V, m->nv, m->T, m->nt, N);
            for (int t = 0; t < m->nt; t++) el += sqrt(pow(m->V[3 * m->T[3 * t]] - m->V[3 * m->T[3 * t + 1]], 2) + pow(m->V[3 * m->T[3 * t] + 1] - m->V[3 * m->T[3 * t + 1] + 1], 2) + pow(m->V[3 * m->T[3 * t] + 2] - m->V[3 * m->T[3 * t + 1] + 2], 2));
            const double eps = 1e-3 * el / m->nt;
            closed += mk_sepvol_closed_(m->T, m->nt) == 1;
            for (int i = 0; i < m->nv; i++, npt++) {
                const double p[3] = {m->V[3 * i] - eps * N[3 * i], m->V[3 * i + 1] - eps * N[3 * i + 1], m->V[3 * i + 2] - eps * N[3 * i + 2]};
                const int a = mk_sepvol_two_(&B, p, 1), b = mk_sepvol_two_(&B, p, 0);
                dis += a != b, nin += a;
                int odd = 0, even = 0, w0 = 0, nw = 0;
                for (int k = 0; k < 3; k++) {
                    int edge, w, h = mk_bvh_count_wind_(&B, p, dirs[k], &edge, &w);
                    if (!edge) odd += h & 1, even += !(h & 1), w0 = nw++ ? w0 : w;
                }
                if (even && !odd) e0 += w0 == 0, e2 += w0 == 2;
            }
            mk_bvh_free(&B), free(N), free(m->V), free(m->T);
        }
        CHECK(closed == 3 && dis == 0 && nin > 0 && e0 > 0 && e2 > 0, "the facings for the winding number (closed meshes): at %ld points (twoball, two cubes, the 0.2 mm and 2 m cubes) %ld decisions unlike mk_winding_number's; %ld inside another part; even parity with facings 2 at %ld points, 0 at %ld",
              npt, dis, nin, e2, e0);
    }
    /* 16. the facings only on a closed mesh. A cube 1 across, 4 x 4 a face, and an open sheet (two triangles, 11
       across, facing +z) through it at z = 0.7: below the sheet each of mk_sepvol_two_'s rays leaves the cube and
       passes the sheet, facings 2, while the winding number is 1 and the sheet's solid angle, under 1.5. The mesh has
       open edges, so the winding number decides */
    {
        mesh m = {0};
        cube(&m, 0, 0, 0, 1, 4);
        const int nv0 = m.nv, nt0 = m.nt;
        m.V = realloc(m.V, sizeof(double) * 3 * (m.nv + 4)), m.T = realloc(m.T, sizeof(int) * 3 * (m.nt + 2));
        const double sv[4][3] = {{-5, -5, 0.7}, {6, -5, 0.7}, {6, 6, 0.7}, {-5, 6, 0.7}};
        memcpy(m.V + 3 * m.nv, sv, sizeof sv), m.nv += 4;
        const int st[6] = {nv0, nv0 + 1, nv0 + 2, nv0, nv0 + 2, nv0 + 3};
        memcpy(m.T + 3 * m.nt, st, sizeof st), m.nt += 2;
        const int cc = mk_sepvol_closed_(m.T, nt0), cs = mk_sepvol_closed_(m.T, m.nt);
        mk_bvh B;
        double *N = malloc(sizeof(double) * 3 * m.nv), el = 0;
        const int rb = mk_bvh_build(&B, m.V, m.T, m.nt);
        mk_sepvol_normals_(m.V, m.nv, m.T, m.nt, N);
        for (int t = 0; t < m.nt; t++) el += sqrt(pow(m.V[3 * m.T[3 * t]] - m.V[3 * m.T[3 * t + 1]], 2) + pow(m.V[3 * m.T[3 * t] + 1] - m.V[3 * m.T[3 * t + 1] + 1], 2) + pow(m.V[3 * m.T[3 * t] + 2] - m.V[3 * m.T[3 * t + 1] + 2], 2));
        const double eps = 1e-3 * el / m.nt;
        int nf = 0, nw = 0, nwind = 0;
        double wmax = 0;
        for (int i = 0; i < nv0 && !rb; i++) {
            if (m.V[3 * i + 2] > 0.7) continue;
            const double p[3] = {m.V[3 * i] - eps * N[3 * i], m.V[3 * i + 1] - eps * N[3 * i + 1], m.V[3 * i + 2] - eps * N[3 * i + 2]};
            const double w = mk_winding_number(m.V, m.T, m.nt, p);
            nf += mk_sepvol_two_(&B, p, 1), nw += mk_sepvol_two_(&B, p, cs), nwind += w > 1.5, wmax = fmax(wmax, w);
        }
        CHECK(!rb && cc == 1 && cs == 0 && nf > 50 && nw == nwind && nwind == 0, "an open sheet through a cube: the cube closed (%d), with the sheet not (%d); below the sheet the facings read %d vertices inside, the winding number (at most %.3f) %d, mk_sepvol_two_ %d",
              cc, cs, nf, wmax, nwind, nw);
        if (!rb) mk_bvh_free(&B);
        free(N), free(m.V), free(m.T);
    }

    /* 17. a point covered three times has odd parity, yet lies inside another part. Three unit spheres with a common
       region, 16 x 10: the vertices inside another part, counted by the winding number of the other parts alone (159,
       47 of them inside two others), against the header's count */
    {
        mesh m = {0};
        for (int p = 0; p < 3; p++) {
            const double c[3] = {0.5 * cos(2 * T_PI * p / 3), 0.5 * sin(2 * T_PI * p / 3), 0.01 * p};
            lat_sphere(&m, 16, 10, c);
        }
        long twice;
        const long ref = others_inside(m.V, m.nv, m.T, m.nt, &twice);
        double *V = malloc(sizeof(double) * 3 * m.nv);
        memcpy(V, m.V, sizeof(double) * 3 * m.nv);
        rc = mk_sepvol(V, m.nv, m.T, m.nt, NULL, 0, NULL, 0, 0.01, 0.05, 2, 1, &r);
        CHECK(rc == 0 && ref == 159 && twice == 47 && r.inside0 == ref, "three spheres: %ld vertices inside another part by the others' winding number, %ld of them inside two others; the header reads %ld", ref, twice, r.inside0);
        /* 19. the same with a mesh that is not closed: one open triangle far off, so the facings cannot stand for the
           winding number */
        const double farv[9] = {50, 0, 0, 51, 0, 0, 50, 1, 0};
        const int nv0 = m.nv, ot[3] = {nv0, nv0 + 1, nv0 + 2};
        m.V = realloc(m.V, sizeof(double) * 3 * (m.nv + 3)), m.T = realloc(m.T, sizeof(int) * 3 * (m.nt + 1));
        memcpy(m.V + 3 * m.nv, farv, sizeof farv), memcpy(m.T + 3 * m.nt, ot, sizeof ot), m.nv += 3, m.nt++;
        V = realloc(V, sizeof(double) * 3 * m.nv);
        memcpy(V, m.V, sizeof(double) * 3 * m.nv);
        rc = mk_sepvol(V, m.nv, m.T, m.nt, NULL, 0, NULL, 0, 0.01, 0.05, 2, 1, &r);
        CHECK(rc == 0 && mk_sepvol_closed_(m.T, m.nt) == 0 && r.inside0 == ref, "three spheres and an open triangle far off (the mesh not closed): the header reads %ld inside (%ld)", r.inside0, ref);
        free(V), free(m.V), free(m.T);
    }
    /* 18. regular nested tetrahedra, one to four shells, none crossing: every vertex of every inner shell is inside
       another part */
    {
        int ok = 0;
        long got[5] = {0};
        for (int n = 1; n <= 4; n++) {
            mesh m = {0};
            for (int j = 0; j < n; j++) tetra(&m, 4 - j, 0, 0);
            rc = mk_sepvol(m.V, m.nv, m.T, m.nt, NULL, 0, NULL, 0, 0, 1e-4, 0, 1, &r);
            got[n] = r.inside0, ok += rc == 0 && r.inside0 == 4 * (n - 1) && r.cross0 == 0;
            free(m.V), free(m.T);
        }
        CHECK(ok == 4, "nested tetrahedra, 1 to 4 shells: %ld, %ld, %ld, %ld vertices inside another part (0, 4, 8, 12 wanted)", got[1], got[2], got[3], got[4]);
    }
    /* 20. each part of the held mesh wound outward, unless it is a cavity inside another. A held mesh with a separate
       tetrahedron wound inward (refused, though its total volume is positive); a hollow cube (an inner shell wound
       inward) with M in its wall and in its cavity; an outward island in the cavity; an inward shell in the cavity */
    {
        mesh H = {0}, M = {0};
        tetra(&H, 2, 0, 0), tetra(&H, 0.5, 10, 1), tetra(&M, 0.1, 10, 0);
        double h0v[24];
        memcpy(h0v, H.V, sizeof h0v);
        rc = mk_sepvol(M.V, M.nv, M.T, M.nt, H.V, H.nv, H.T, H.nt, 0, 1e-4, 0, 1, &r);
        CHECK(rc == -2 && r.status == MK_SEPVOL_HELD_FLIPPED && r.held_inward == 1 && r.held_vol > 0 && !memcmp(h0v, H.V, sizeof h0v), "a held mesh with a separate tetrahedron wound inward: returns %d, status %d, %ld parts wound inward, total volume %.3f (positive)", rc, r.status, r.held_inward, r.held_vol);
        free(H.V), free(H.T), free(M.V), free(M.T);
        mesh Ho = {0}, W = {0}, C = {0};
        box1(&Ho, 0, 0, 0, 4, 0), box1(&Ho, 1, 1, 1, 2, 1);
        box1(&W, 0.25, 1.75, 1.75, 0.5, 0), box1(&C, 1.75, 1.75, 1.75, 0.5, 0);
        rc = mk_sepvol(W.V, W.nv, W.T, W.nt, Ho.V, Ho.nv, Ho.T, Ho.nt, 0.01, 0.05, 0, 3, &r);
        const int okw = rc == 0 && r.held_inward == 0 && !(r.status & MK_SEPVOL_HELD_FLIPPED) && r.inside0 == 8;
        const long inw = r.inside0;
        rc = mk_sepvol(C.V, C.nv, C.T, C.nt, Ho.V, Ho.nv, Ho.T, Ho.nt, 0.01, 0.05, 0, 3, &r);
        CHECK(okw && rc == 0 && r.status == 0 && r.inside0 == 0, "a hollow cube held (inner shell wound inward): M in its wall %ld of 8 vertices inside, M in its cavity %ld, status %d", inw, r.inside0, r.status);
        mesh Is = {0}, In = {0};
        box1(&Is, 0, 0, 0, 4, 0), box1(&Is, 1, 1, 1, 2, 1), box1(&Is, 1.5, 1.5, 1.5, 1, 0);
        box1(&In, 0, 0, 0, 4, 0), box1(&In, 1, 1, 1, 2, 1), box1(&In, 1.5, 1.5, 1.5, 1, 1);
        rc = mk_sepvol(W.V, W.nv, W.T, W.nt, Is.V, Is.nv, Is.T, Is.nt, 0.01, 0.05, 0, 1, &r);
        const int oki = rc == 0 && !(r.status & MK_SEPVOL_HELD_FLIPPED);
        rc = mk_sepvol(W.V, W.nv, W.T, W.nt, In.V, In.nv, In.T, In.nt, 0.01, 0.05, 0, 1, &r);
        CHECK(oki && rc == -2 && (r.status & MK_SEPVOL_HELD_FLIPPED) && r.held_inward == 1, "an outward island in the cavity: taken; an inward shell in the cavity (no solid around it): returns %d, %ld parts wound inward", rc, r.held_inward);
        free(Ho.V), free(Ho.T), free(W.V), free(W.T), free(C.V), free(C.T), free(Is.V), free(Is.T), free(In.V), free(In.T);
    }
    /* 21. two bars 6 x 0.2 x 0.2 crossed like an X, no vertex of either inside the other: moving vertices cannot clear
       them, and squeezing them through costs a large share of their volume. The run stops, hands back the best state
       before a part lost more than MK_SEPVOL_CRUSH, and says the parts still cross */
    {
        mesh m = {0};
        unsigned long long g = 7;
        diag_cube(&m, &g); /* [-1, 1]^3, each square split along a pseudo-random diagonal (xorshift, seed 7) */
        for (int i = 0; i < m.nv; i++) m.V[3 * i] *= 3, m.V[3 * i + 1] *= 0.1, m.V[3 * i + 2] *= 0.1;
        const int v0 = m.nv;
        diag_cube(&m, &g);
        for (int i = v0; i < m.nv; i++) m.V[3 * i] *= 0.1, m.V[3 * i + 1] *= 3, m.V[3 * i + 2] = 0.1 * m.V[3 * i + 2] + 0.05;
        double *V = malloc(sizeof(double) * 3 * m.nv);
        memcpy(V, m.V, sizeof(double) * 3 * m.nv);
        rc = mk_sepvol(V, m.nv, m.T, m.nt, NULL, 0, NULL, 0, 0.01, 0.05, 1, 60, &r);
        CHECK(rc == 0 && r.part_loss <= MK_SEPVOL_CRUSH && r.vol1 >= (1 - MK_SEPVOL_CRUSH) * r.vol0 && r.inside == 0, "two crossed bars: volume %.4f -> %.4f, the most any part lost %.1f %% (limit %.0f %%), %ld of %ld crossing pairs left after %d rounds", r.vol0, r.vol1, 100 * r.part_loss, 100 * (double)MK_SEPVOL_CRUSH, r.cross, r.cross0, r.rounds);
        CHECK((r.status & MK_SEPVOL_CRUSHED) && (r.status & MK_SEPVOL_STILL_CROSSES) && r.cross > 0, "two crossed bars: the status says the run stopped (%d) and the parts still cross (%d): status %d", (r.status & MK_SEPVOL_CRUSHED) != 0, (r.status & MK_SEPVOL_STILL_CROSSES) != 0, r.status);
        memcpy(V, m.V, sizeof(double) * 3 * m.nv);
        rc = mk_sepvol(V, m.nv, m.T, m.nt, NULL, 0, NULL, 0, 0.01, 0.05, 1, 2, &r);
        CHECK(rc == 0 && r.status == MK_SEPVOL_STILL_CROSSES && r.cross > 0, "two crossed bars, two rounds only: pairs left (%ld) and nothing crushed: status %d", r.cross, r.status);
        memcpy(tb.V, tb0, sizeof(double) * 3 * tb.nv);
        rc = mk_sepvol(tb.V, tb.nv, tb.T, tb.nt, NULL, 0, NULL, 0, gap, step, rings, rounds, &r);
        CHECK(rc == 0 && r.status == 0 && r.cross == 0 && r.part_loss > 0 && r.part_loss < 0.1, "twoball, a separation that succeeds: status %d, the most any part lost %.2f %%", r.status, 100 * r.part_loss);
        free(V), free(m.V), free(m.T);
    }

    /* 22. (held mode) a box 6.05 mm thick and 7.32 mm square (3 x 3 a face), 3.61 mm into a cube 10 mm across
       (18 x 18 a face), off its middle; gap 0.1 mm, step 0.649 mm, 3 rings, 22 rounds. Pushing out of a part this thick
       folds the cube's face. The state handed back has no more inside vertices plus crossing pairs than the input, no
       edge folded over 90 degrees that met at under 90 at the start (the cube's own edges, at 90, may fold further),
       and the report says what is left */
    {
        mesh K = {0};
        cube(&K, 0, 0, 0, 0.01, 18);
        const double lo[3] = {0.01 - 3.607e-3, 6.469e-3 - 3.659e-3, 0.005 - 3.659e-3}, hi[3] = {0.01 - 3.607e-3 + 6.048e-3, 6.469e-3 + 3.659e-3, 0.005 + 3.659e-3};
        mesh R = box(lo, hi, 3);
        double *K0 = malloc(sizeof(double) * 3 * K.nv);
        memcpy(K0, K.V, sizeof(double) * 3 * K.nv);
        const long kin0 = held_inside(K.V, K.nv, &R), kr0 = pair_cross(&K, &R);
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, 1e-4, 6.49e-4, 3, 22, &r);
        const long kin1 = held_inside(K.V, K.nv, &R), kr1 = pair_cross(&K, &R), k1 = own_pairs(K.V, K.T, K.nt, NULL), nf = new_folds(K0, K.V, K.T, K.nt);
        CHECK(rc == 0 && kin0 > 100 && kr0 > 100 && kin1 + kr1 <= kin0 + kr0 && r.inside + r.cross <= r.inside0 + r.cross0 && r.inside < r.inside0 && r.inside == kin1 && r.cross == kr1, "a thick box 3.6 mm into a cube: never more inside plus crossing than it came in: vertices inside %ld -> %ld, crossing pairs %ld -> %ld (the report: %ld -> %ld, %ld -> %ld)", kin0, kin1, kr0, kr1, r.inside0, r.inside, r.cross0, r.cross);
        CHECK(rc == 0 && nf == 0 && k1 == 0, "a thick box 3.6 mm into a cube: %ld edges folded over 90 degrees that met at under 90 (the test's count), own crossing pairs %ld", nf, k1);
        const double f0 = fold_length(K0, K.T, K.nt), f1 = fold_length(K.V, K.T, K.nt);
        CHECK(fabs(r.fold0 - f0) <= 1e-9 * fmax(f0, f1) && fabs(r.fold1 - f1) <= 1e-9 * fmax(f0, f1) && f1 > 0, "a thick box 3.6 mm into a cube: the report's fold lengths %.4f -> %.4f mm (the test's: %.4f -> %.4f mm)", r.fold0 * 1e3, r.fold1 * 1e3, f0 * 1e3, f1 * 1e3);
        CHECK(r.inside > 0 && r.cross > 0 && (r.status & MK_SEPVOL_STILL_INSIDE) && (r.status & MK_SEPVOL_STILL_CROSSES) && !(r.status & MK_SEPVOL_PAIRS_ROSE) == !(r.cross > r.cross0), "a thick box 3.6 mm into a cube: vertices (%ld) and pairs (%ld, %ld at the start) left, and the status says so: %d", r.inside, r.cross, r.cross0, r.status);
        free(K.V), free(K.T), free(R.V), free(R.T), free(K0);
    }
    { /* where the push raises the crossing pairs: a box 5.70 mm thick and 7.82 mm square (12 x 12 a face), 4.39 mm into a
         cube 15 x 15 a face; step 0.846 mm, 1 ring, 25 rounds. The fewest inside vertices plus crossing pairs come with
         more crossing pairs than the input, and the status says so (the fewest inside vertices alone: 132 -> 30 inside,
         190 -> 362 crossing pairs, more in sum than the input) */
        mesh K = {0};
        cube(&K, 0, 0, 0, 0.01, 15);
        const double lo[3] = {0.01 - 4.386e-3, 3.977e-3 - 3.912e-3, 0.005 - 3.912e-3}, hi[3] = {0.01 - 4.386e-3 + 5.698e-3, 3.977e-3 + 3.912e-3, 0.005 + 3.912e-3};
        mesh R = box(lo, hi, 12);
        const long kin0 = held_inside(K.V, K.nv, &R), kr0 = pair_cross(&K, &R);
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, 1e-4, 8.46e-4, 1, 25, &r);
        const long kin1 = held_inside(K.V, K.nv, &R), kr1 = pair_cross(&K, &R);
        CHECK(rc == 0 && kin0 > 100 && kin1 + kr1 <= kin0 + kr0 && r.inside + r.cross <= r.inside0 + r.cross0 && kin1 < kin0, "a box 4.4 mm into a cube, the push adding crossing pairs: vertices inside %ld -> %ld, crossing pairs %ld -> %ld (the report: %ld -> %ld, %ld -> %ld)", kin0, kin1, kr0, kr1, r.inside0, r.inside, r.cross0, r.cross);
        CHECK(r.cross > r.cross0 && (r.status & MK_SEPVOL_PAIRS_ROSE), "a box 4.4 mm into a cube: the crossing pairs rose (%ld -> %ld), and the status says so: %d", r.cross0, r.cross, r.status);
        free(K.V), free(K.T), free(R.V), free(R.T);
    }
    /* 23. (held mode) the moves smoothed in space: a box 5.6 mm thick and 4.7 mm square (2 x 2 a face), 1.5 mm into
       a cube 10 mm across (19 x 19 a face); gap 0.1 mm, step 0.816 mm, 2 rings, 15 rounds. Each vertex's own push folds
       the face past 90 degrees, so the guard holds the moves back; smoothed, the patch moves together and clears */
    {
        mesh K = {0};
        cube(&K, 0, 0, 0, 0.01, 19);
        const double lo[3] = {0.01 - 1.526e-3, 5.018e-3 - 2.329e-3, 0.005 - 2.329e-3}, hi[3] = {0.01 - 1.526e-3 + 5.614e-3, 5.018e-3 + 2.329e-3, 0.005 + 2.329e-3};
        mesh R = box(lo, hi, 2);
        double *K0 = malloc(sizeof(double) * 3 * K.nv);
        memcpy(K0, K.V, sizeof(double) * 3 * K.nv);
        const long kin0 = held_inside(K.V, K.nv, &R), kr0 = pair_cross(&K, &R);
        rc = mk_sepvol(K.V, K.nv, K.T, K.nt, R.V, R.nv, R.T, R.nt, 1e-4, 8.16e-4, 2, 15, &r);
        const long kin1 = held_inside(K.V, K.nv, &R), kr1 = pair_cross(&K, &R), k1 = own_pairs(K.V, K.T, K.nt, NULL), nf = new_folds(K0, K.V, K.T, K.nt);
        CHECK(rc == 0 && kin0 > 30 && kr0 > 30 && kin1 == 0 && kr1 == 0 && r.inside == 0 && r.cross == 0 && nf == 0 && k1 == 0, "a box 1.5 mm into a cube, the moves smoothed in space: vertices inside %ld -> %ld, crossing pairs %ld -> %ld, %ld new folds over 90 degrees, own crossing pairs %ld; %d rounds", kin0, kin1, kr0, kr1, nf, k1, r.rounds);
        const double ov = 1.526e-3 * 4.658e-3 * 4.658e-3;
        CHECK(rc == 0 && r.vol0 - r.vol1 >= ov && r.vol0 - r.vol1 < 3 * ov, "a box 1.5 mm into a cube: volume %.4f -> %.4f ml, %.2f times the overlap (the smoothing spreads the dent; at most 3)", r.vol0 * 1e6, r.vol1 * 1e6, (r.vol0 - r.vol1) / ov);
        free(K.V), free(K.T), free(R.V), free(R.T), free(K0);
    }
    /* 23b. row 11's check on M's own surface, with the moves smoothed: a bar 1.9 mm wide through row 11's slab,
       its top 0.35 mm above the slab's (the slab jittered by its own xorshift, seed 11). Without the check the slab's
       top is pushed down through the slab and it loses some 10 % of its volume; with it, under 5 % */
    {
        const double lo[3] = {0, 0, 0}, hi[3] = {0.02, 0.02, 0.003}, blo[3] = {0.0101 - 0.00095, -0.005, -0.008}, bhi[3] = {0.0101 + 0.00095, 0.025, 0.00335};
        mesh S = box(lo, hi, 20), R = box(blo, bhi, 8);
        unsigned long long g = 11;
        for (int i = 0; i < 3 * S.nv; i++) {
            g ^= g << 13, g ^= g >> 7, g ^= g << 17;
            S.V[i] += 1e-7 * ((double)(g >> 11) * 0x1.0p-52 - 1.0);
        }
        const long sin0 = held_inside(S.V, S.nv, &R), sr0 = pair_cross(&S, &R);
        rc = mk_sepvol(S.V, S.nv, S.T, S.nt, R.V, R.nv, R.T, R.nt, gap, step, rings, rounds, &r);
        const long s1 = own_pairs(S.V, S.T, S.nt, NULL), sin1 = held_inside(S.V, S.nv, &R), sr1 = pair_cross(&S, &R);
        CHECK(rc == 0 && sin0 > 100 && s1 == 0 && sin1 + sr1 <= sin0 + sr0 && r.vol1 > 0.95 * r.vol0, "a bar 1.9 mm wide through a slab: its volume %.4f -> %.4f ml (%.2f %%), own crossing pairs %ld, vertices inside the bar %ld -> %ld, crossing pairs %ld -> %ld", r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1), s1, sin0, sin1, sr0, sr1);
        free(S.V), free(S.T), free(R.V), free(R.T);
    }
    /* 24. several held meshes in one call, each kept separate. A ball 5 mm across (60 x 30) against two boxes that
       overlap each other, one through its +x side and one through its +y side, and an open box near its -x side. Their
       union crosses itself, so mk_sepvol refuses it (-2); mk_sepvol_multi takes the two boxes, each its own inside, and
       leaves the open one out with its report saying why */
    {
        mesh M = {0};
        const double c0[3] = {0, 0, 0};
        lat_sphere(&M, 60, 30, c0);
        for (int i = 0; i < 3 * M.nv; i++) M.V[i] *= 5e-3;
        const double alo[3] = {3.5e-3, -2e-3, -3e-3}, ahi[3] = {9e-3, 4.5e-3, 3e-3}, blo[3] = {-2e-3, 3.8e-3, -3e-3}, bhi[3] = {4.5e-3, 9e-3, 3e-3};
        const double clo[3] = {-9e-3, -2e-3, -2e-3}, chi[3] = {-4.2e-3, 2e-3, 2e-3};
        mesh A = box(alo, ahi, 8), Bx = box(blo, bhi, 8), C = box(clo, chi, 4), U = {0};
        C.nt--; /* open: a face short */
        const double z[3] = {0, 0, 0};
        append(&U, &A, z), append(&U, &Bx, z);
        double *M0 = malloc(sizeof(double) * 3 * M.nv), *W = malloc(sizeof(double) * 3 * M.nv), *A0 = malloc(sizeof(double) * 3 * A.nv), *B0 = malloc(sizeof(double) * 3 * Bx.nv);
        int *MT0 = malloc(sizeof(int) * 3 * M.nt);
        memcpy(M0, M.V, sizeof(double) * 3 * M.nv), memcpy(W, M.V, sizeof(double) * 3 * M.nv), memcpy(MT0, M.T, sizeof(int) * 3 * M.nt);
        memcpy(A0, A.V, sizeof(double) * 3 * A.nv), memcpy(B0, Bx.V, sizeof(double) * 3 * Bx.nv);
        rc = mk_sepvol(W, M.nv, M.T, M.nt, U.V, U.nv, U.T, U.nt, 1e-4, 2.5e-4, 3, 30, &r);
        CHECK(rc == -2 && (r.status & MK_SEPVOL_HELD_CROSSES) && !memcmp(W, M0, sizeof(double) * 3 * M.nv), "the two boxes as one held mesh: refused (%d, status %d, %ld crossing pairs of its own)", rc, r.status, r.held_cross);
        const long ia0 = held_inside(M.V, M.nv, &A), ib0 = held_inside(M.V, M.nv, &Bx), ca0 = pair_cross(&M, &A), cb0 = pair_cross(&M, &Bx);
        const double *HV[3] = {A.V, Bx.V, C.V};
        const int *HT[3] = {A.T, Bx.T, C.T}, hnv[3] = {A.nv, Bx.nv, C.nv}, hnt[3] = {A.nt, Bx.nt, C.nt};
        mk_sepvol_report hr[3];
        rc = mk_sepvol_multi(M.V, M.nv, M.T, M.nt, 3, HV, hnv, HT, hnt, 1e-4, 2.5e-4, 3, 30, &r, hr);
        const long ia1 = held_inside(M.V, M.nv, &A), ib1 = held_inside(M.V, M.nv, &Bx), ca1 = pair_cross(&M, &A), cb1 = pair_cross(&M, &Bx), k1 = own_pairs(M.V, M.T, M.nt, NULL);
        CHECK(rc == 0 && ia0 > 50 && ib0 > 50 && ia1 == 0 && ib1 == 0 && ca1 == 0 && cb1 == 0 && k1 == 0, "the ball out of both boxes in one call: inside %ld, %ld -> %ld, %ld; crossing pairs %ld, %ld -> %ld, %ld (the test's counts); own crossing pairs %ld", ia0, ib0, ia1, ib1, ca0, cb0, ca1, cb1, k1);
        CHECK(hr[0].inside0 == ia0 && hr[1].inside0 == ib0 && hr[0].cross0 == ca0 && hr[1].cross0 == cb0 && hr[0].inside == 0 && hr[1].inside == 0 && hr[0].cross == 0 && hr[1].cross == 0 && r.inside0 == ia0 + ib0 && r.cross0 == ca0 + cb0 && r.inside == 0 && r.cross == 0 && hr[0].rounds > 0 && hr[1].rounds > 0, "each box's report: inside %ld -> %ld and %ld -> %ld, crossing pairs %ld -> %ld and %ld -> %ld, rounds %d and %d; summed %ld -> %ld, %ld -> %ld", hr[0].inside0, hr[0].inside, hr[1].inside0, hr[1].inside, hr[0].cross0, hr[0].cross, hr[1].cross0, hr[1].cross, hr[0].rounds, hr[1].rounds, r.inside0, r.inside, r.cross0, r.cross);
        CHECK(hr[2].status == MK_SEPVOL_HELD_OPEN && hr[2].held_open > 0 && hr[2].rounds == 0 && (r.status & MK_SEPVOL_HELD_OPEN) && !(r.status & ~MK_SEPVOL_HELD_OPEN), "the open box left out: its status %d (%ld open edges), the call's %d", hr[2].status, hr[2].held_open, r.status);
        CHECK(!memcmp(A0, A.V, sizeof(double) * 3 * A.nv) && !memcmp(B0, Bx.V, sizeof(double) * 3 * Bx.nv) && !memcmp(MT0, M.T, sizeof(int) * 3 * M.nt) && r.vol1 < r.vol0 && r.vol1 > 0.8 * r.vol0 && r.self0 == 0 && r.self1 == 0, "the boxes unchanged byte for byte, the ball's faces too; its volume %.4f -> %.4f ml (%.2f %%), own pairs %ld -> %ld", r.vol0 * 1e6, r.vol1 * 1e6, 100 * (r.vol1 / r.vol0 - 1), r.self0, r.self1);
        /* refusals: no held mesh, a NULL one, and every one left out */
        memcpy(W, M.V, sizeof(double) * 3 * M.nv);
        const double *HN[1] = {NULL};
        const int nn = mk_sepvol_multi(W, M.nv, M.T, M.nt, 0, HV, hnv, HT, hnt, 1e-4, 2.5e-4, 3, 30, &r, hr) == -1 && mk_sepvol_multi(W, M.nv, M.T, M.nt, 1, HN, hnv, HT, hnt, 1e-4, 2.5e-4, 3, 30, &r, hr) == -1 &&
                       mk_sepvol_multi(W, M.nv, M.T, M.nt, 2, HV, hnv, HT, hnt, 1e-4, 0, 3, 30, &r, hr) == -1;
        rc = mk_sepvol_multi(W, M.nv, M.T, M.nt, 1, HV + 2, hnv + 2, HT + 2, hnt + 2, 1e-4, 2.5e-4, 3, 30, &r, hr);
        CHECK(nn && rc == -2 && hr[0].status == MK_SEPVOL_HELD_OPEN && !memcmp(W, M.V, sizeof(double) * 3 * M.nv), "mk_sepvol_multi refuses no held mesh, a NULL one and step 0 (-1), and returns -2 with M unchanged when every held mesh is left out (status %d)", hr[0].status);
        free(M.V), free(M.T), free(A.V), free(A.T), free(Bx.V), free(Bx.T), free(C.V), free(C.T), free(U.V), free(U.T), free(M0), free(W), free(A0), free(B0), free(MT0);
    }
    /* 25. a push out of one held mesh can take M into another: a cube 10 mm across (17 x 17 a face), a box 3.44
       mm thick 3.16 mm into its +x face (5.23 mm square, 7 x 7 a face) and one 2.98 mm thick 3.76 mm in (6.86 mm square,
       8 x 8 a face) beside it, overlapping it; step 0.648 mm, 2 rings, 4 rounds. Each push sinks the face further into
       the other box; the state handed back has no more inside vertices plus crossing pairs, summed over both, than the
       start (the last state: 394 -> 513), and the reports agree with the test's counts */
    {
        mesh K = {0};
        cube(&K, 0, 0, 0, 0.01, 17);
        const double alo[3] = {0.01 - 3.157e-3, 6.447e-3 - 2.6155e-3, 4.905e-3 - 2.6155e-3}, ahi[3] = {0.01 - 3.157e-3 + 3.438e-3, 6.447e-3 + 2.6155e-3, 4.905e-3 + 2.6155e-3};
        const double blo[3] = {0.01 - 3.758e-3, 5.830e-3 - 3.4295e-3, 7.246e-3 - 3.4295e-3}, bhi[3] = {0.01 - 3.758e-3 + 2.975e-3, 5.830e-3 + 3.4295e-3, 7.246e-3 + 3.4295e-3};
        mesh A = box(alo, ahi, 7), Bx = box(blo, bhi, 8);
        const double *HV[2] = {A.V, Bx.V};
        const int *HT[2] = {A.T, Bx.T}, hnv[2] = {A.nv, Bx.nv}, hnt[2] = {A.nt, Bx.nt};
        const long ia0 = held_inside(K.V, K.nv, &A), ib0 = held_inside(K.V, K.nv, &Bx), ca0 = pair_cross(&K, &A), cb0 = pair_cross(&K, &Bx);
        mk_sepvol_report hr[2];
        rc = mk_sepvol_multi(K.V, K.nv, K.T, K.nt, 2, HV, hnv, HT, hnt, 1e-4, 6.48e-4, 2, 4, &r, hr);
        const long ia1 = held_inside(K.V, K.nv, &A), ib1 = held_inside(K.V, K.nv, &Bx), ca1 = pair_cross(&K, &A), cb1 = pair_cross(&K, &Bx);
        CHECK(rc == 0 && ia0 > 20 && ib0 > 20 && ia1 + ib1 + ca1 + cb1 <= ia0 + ib0 + ca0 + cb0 && hr[0].inside == ia1 && hr[1].inside == ib1 && hr[0].cross == ca1 && hr[1].cross == cb1 && r.inside + r.cross <= r.inside0 + r.cross0, "two held boxes beside each other: inside %ld, %ld -> %ld, %ld; crossing pairs %ld, %ld -> %ld, %ld (the test's counts, the reports' the same); summed %ld -> %ld", ia0, ib0, ia1, ib1, ca0, cb0, ca1, cb1, ia0 + ib0 + ca0 + cb0, ia1 + ib1 + ca1 + cb1);
        CHECK(!(r.status & MK_SEPVOL_PAIRS_ROSE) == !(r.cross > r.cross0) && !(hr[0].status & MK_SEPVOL_PAIRS_ROSE) == !(ca1 > ca0) && !(hr[1].status & MK_SEPVOL_PAIRS_ROSE) == !(cb1 > cb0), "two held boxes: MK_SEPVOL_PAIRS_ROSE where the crossing pairs rose: the call's %d, each box's %d and %d", (r.status & MK_SEPVOL_PAIRS_ROSE) != 0, (hr[0].status & MK_SEPVOL_PAIRS_ROSE) != 0, (hr[1].status & MK_SEPVOL_PAIRS_ROSE) != 0);
        free(K.V), free(K.T), free(A.V), free(A.T), free(Bx.V), free(Bx.T);
    }
    /* 26. no memory, held mode and several held meshes: as row 13, each allocation of a call failing in turn. A
       cube 4 x 4 a face against two boxes through its +x face, 3 rounds; mk_sepvol with the first (its smoothing in
       space and its folds allocate too), then mk_sepvol_multi with both: -1, V as it came in, the reports zeroed */
    {
        mesh m = {0};
        cube(&m, 0, 0, 0, 0.01, 4);
        const double alo[3] = {0.0092, 0.001, 0.001}, ahi[3] = {0.012, 0.004, 0.009}, blo[3] = {0.0094, 0.006, 0.001}, bhi[3] = {0.012, 0.009, 0.009};
        mesh A = box(alo, ahi, 3), Bx = box(blo, bhi, 3);
        double *m0 = malloc(sizeof(double) * 3 * m.nv);
        memcpy(m0, m.V, sizeof(double) * 3 * m.nv);
        const double *HV[2] = {A.V, Bx.V};
        const int *HT[2] = {A.T, Bx.T}, hnv[2] = {A.nv, Bx.nv}, hnt[2] = {A.nt, Bx.nt};
        const mk_sepvol_report zero = {0};
        mk_sepvol_report hr[2];
        long nalloc[2], nok[2] = {0, 0}, moved[2];
        for (int f = 0; f < 2; f++) {
            memcpy(m.V, m0, sizeof(double) * 3 * m.nv);
            t_fi_at = LONG_MAX, t_fi_n = 0;
            rc = f ? mk_sepvol_multi(m.V, m.nv, m.T, m.nt, 2, HV, hnv, HT, hnt, gap, step, rings, 3, &r, hr) : mk_sepvol(m.V, m.nv, m.T, m.nt, A.V, A.nv, A.T, A.nt, gap, step, rings, 3, &r);
            nalloc[f] = rc == 0 ? t_fi_n : -1, moved[f] = r.moved;
            for (long k = 0; k < nalloc[f]; k++) {
                memcpy(m.V, m0, sizeof(double) * 3 * m.nv);
                memset(&r, 0x5A, sizeof r), memset(hr, 0x5A, sizeof hr);
                t_fi_at = k, t_fi_n = 0;
                const int rk = f ? mk_sepvol_multi(m.V, m.nv, m.T, m.nt, 2, HV, hnv, HT, hnt, gap, step, rings, 3, &r, hr) : mk_sepvol(m.V, m.nv, m.T, m.nt, A.V, A.nv, A.T, A.nt, gap, step, rings, 3, &r);
                nok[f] += rk == -1 && !memcmp(m.V, m0, sizeof(double) * 3 * m.nv) && !memcmp(&r, &zero, sizeof r) && (!f || (!memcmp(hr, &zero, sizeof zero) && !memcmp(hr + 1, &zero, sizeof zero)));
            }
            t_fi_at = -1;
        }
        CHECK(nalloc[0] > 20 && nalloc[1] > nalloc[0] && moved[0] > 0 && moved[1] > 0 && nok[0] == nalloc[0] && nok[1] == nalloc[1], "allocation failures, held mode: %ld of mk_sepvol's %ld and %ld of mk_sepvol_multi's %ld return -1 with V exactly as it came in and the reports zeroed (without failures %ld and %ld vertices moved)", nok[0], nalloc[0], nok[1], nalloc[1], moved[0], moved[1]);
        free(m.V), free(m.T), free(m0), free(A.V), free(A.T), free(Bx.V), free(Bx.T);
    }
    /* 27. the volume given back where there is room: ballB pushed out of held ballA lost volume; mk_sepvol_restore
       moves its vertices 1 mm or more from ballA out along their normals toward the volume it had. Nothing comes nearer
       ballA than 0.5 mm, and nothing crosses */
    {
        mesh a = readobj("test/data/sepvol_ballA.obj"), b = readobj("test/data/sepvol_ballB.obj");
        rc = mk_sepvol(b.V, b.nv, b.T, b.nt, a.V, a.nv, a.T, a.nt, gap, step, rings, rounds, &r);
        const double v0 = r.vol0, lost = r.vol0 - r.vol1, depth = 1e-3;
        double *b0 = malloc(sizeof(double) * 3 * b.nv);
        memcpy(b0, b.V, sizeof(double) * 3 * b.nv);
        const double *HV[1] = {a.V};
        const int *HT[1] = {a.T}, hnv[1] = {a.nv}, hnt[1] = {a.nt};
        const int rr = mk_sepvol_restore(b.V, b.nv, b.T, b.nt, 1, HV, hnv, HT, hnt, v0, depth, &r);
        double dmin = INFINITY; /* the moved vertices' distance from ballA, before their move */
        for (int i = 0; i < b.nv; i++) {
            if (!memcmp(b.V + 3 * i, b0 + 3 * i, sizeof(double) * 3)) continue;
            for (int u = 0; u < a.nt; u++) {
                double w[3], q[3];
                dmin = fmin(dmin, mk_closest_pt_tri(b0 + 3 * i, a.V + 3 * a.T[3 * u], a.V + 3 * a.T[3 * u + 1], a.V + 3 * a.T[3 * u + 2], w, q));
            }
        }
        dmin = sqrt(dmin);
        const long in1 = held_inside(b.V, b.nv, &a), c1 = pair_cross(&b, &a), k1 = own_pairs(b.V, b.T, b.nt, NULL), nf = new_folds(b0, b.V, b.T, b.nt);
        CHECK(rc == 0 && rr == 0 && lost > 0 && fabs(r.vol1 - v0) < 0.05 * lost && r.moved > 0 && dmin >= depth, "ballB's volume given back: %.4f ml lost in the push, %.4f ml short after (%d steps, %ld vertices moved, the nearest %.3f mm from ballA)", lost * 1e6, (v0 - r.vol1) * 1e6, r.rounds, r.moved, dmin * 1e3);
        CHECK(rr == 0 && in1 == 0 && c1 == 0 && k1 == 0 && nf == 0 && r.inside == 0 && r.cross == 0 && r.self1 == 0, "ballB after: vertices inside ballA %ld, crossing pairs %ld, own %ld, new folds over 90 degrees %ld (the report: %ld, %ld, %ld)", in1, c1, k1, nf, r.inside, r.cross, r.self1);
        free(a.V), free(a.T), free(b.V), free(b.T), free(b0);
    }
    { /* where the corners are far from a held part but a face is not: a cube one square a face, a plate 2 mm square and
         0.1 mm thick 0.3 mm off the middle of its +x face. Its corners are 4 mm and more from the plate;
         given 20 % more volume, that face would sweep over the plate: its corners stay, the other four move */
        mesh k = {0}, p = {0};
        box1(&k, 0, 0, 0, 0.01, 0);
        const double lo[3] = {0.0103, 0.004, 0.004}, hi[3] = {0.0104, 0.006, 0.006};
        p = box(lo, hi, 2);
        double *k0 = malloc(sizeof(double) * 3 * k.nv);
        memcpy(k0, k.V, sizeof(double) * 3 * k.nv);
        const double *HV[1] = {p.V};
        const int *HT[1] = {p.T}, hnv[1] = {p.nv}, hnt[1] = {p.nt};
        rc = mk_sepvol_restore(k.V, k.nv, k.T, k.nt, 1, HV, hnv, HT, hnt, 1.2e-6, 1e-3, &r);
        int still = 0; /* the +x face's corners where they were */
        for (int i = 0; i < k.nv; i++) still += k0[3 * i] > 0.005 && !memcmp(k0 + 3 * i, k.V + 3 * i, sizeof(double) * 3);
        CHECK(rc == 0 && pair_cross(&k, &p) == 0 && r.cross == 0 && r.cross0 == 0 && still == 4 && fabs(r.vol1 - 1.2e-6) < 1e-3 * 1.2e-6 && r.self1 == 0, "a face that would go through a held plate: %ld crossing pairs, its %d corners where they were, the rest moved: volume %.4f -> %.4f ml (1.2 wanted; %d steps)", pair_cross(&k, &p), still, r.vol0 * 1e6, r.vol1 * 1e6, r.rounds);
        memcpy(k.V, k0, sizeof(double) * 3 * k.nv);
        /* refusals: depth not positive, the volume not finite, a held mesh open */
        const int n1 = mk_sepvol_restore(k.V, k.nv, k.T, k.nt, 1, HV, hnv, HT, hnt, 1.2e-6, 0, &r) == -1 && mk_sepvol_restore(k.V, k.nv, k.T, k.nt, 1, HV, hnv, HT, hnt, NAN, 1e-3, &r) == -1;
        p.nt--;
        const int hnt2[1] = {p.nt};
        rc = mk_sepvol_restore(k.V, k.nv, k.T, k.nt, 1, HV, hnv, HT, hnt2, 1.2e-6, 1e-3, &r);
        CHECK(n1 && rc == -2 && (r.status & MK_SEPVOL_HELD_OPEN) && !memcmp(k0, k.V, sizeof(double) * 3 * k.nv), "mk_sepvol_restore refuses depth 0 and a NaN volume (-1), and an open held mesh (%d, status %d), the cube unchanged", rc, r.status);
        rc = mk_sepvol_restore(k.V, k.nv, k.T, k.nt, 0, NULL, NULL, NULL, NULL, 1.2e-6, 1e-3, &r);
        CHECK(rc == 0 && fabs(r.vol1 - 1.2e-6) < 1e-3 * 1.2e-6 && r.self1 == 0, "the cube, no held mesh: %.4f -> %.4f ml (1.2 wanted) in %d steps", r.vol0 * 1e6, r.vol1 * 1e6, r.rounds);
        free(k.V), free(k.T), free(p.V), free(p.T), free(k0);
    }

    /* 5. refusals: V as it came in */
    {
        double V[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0}, V0[9];
        int T[3] = {0, 1, 2}, Tb[3] = {0, 1, 3};
        memcpy(V0, V, sizeof V);
        int n = 0, k = 0;
#define REF(call) (k++, n += (call) == -1 && !memcmp(V, V0, sizeof V))
        REF(mk_sepvol(V, 0, T, 1, NULL, 0, NULL, 0, gap, step, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 0, NULL, 0, NULL, 0, gap, step, 3, 5, &r));
        REF(mk_sepvol(V, 3, Tb, 1, NULL, 0, NULL, 0, gap, step, 3, 5, &r));
        REF(mk_sepvol(V, -1, T, 1, NULL, 0, NULL, 0, gap, step, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, NAN, step, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, -1e-6, step, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, INFINITY, step, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, gap, 0, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, gap, INFINITY, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, gap, NAN, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, gap, step, -1, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, gap, step, 3, 0, &r));
        REF(mk_sepvol(V, 3, T, 1, V0, 3, NULL, 1, gap, step, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, V0, 3, Tb, 1, gap, step, 3, 5, &r));
        REF(mk_sepvol(V, 3, T, 1, V0, 0, T, 1, gap, step, 3, 5, &r));
        V[4] = NAN, V0[4] = NAN; /* compared bytewise: the NaN stays as it came */
        REF(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, gap, step, 3, 5, &r));
        V[4] = V0[4] = 0;
        double H[9] = {0, 0, 0, 1, 0, 0, 0, INFINITY, 0};
        REF(mk_sepvol(V, 3, T, 1, H, 3, T, 1, gap, step, 3, 5, &r));
#undef REF
        CHECK(n == k, "refusals: %d of %d return -1 with V unchanged (sizes, an index out of range, NaN and infinite input, gap, step, rings, rounds, a bad held part)", n, k);
        CHECK(mk_sepvol(V, 3, T, 1, NULL, 0, NULL, 0, gap, step, 3, 5, NULL) == 0, "one lone triangle, no report: 0, nothing to move");
    }
    free(tb.V), free(tb.T), free(tb0), free(tbT0), free(A.V), free(A.T), free(Bm.V), free(Bm.T), free(A0), free(B0), free(AT0), free(C.V), free(C.T), free(C0);
    return wrong;
}
