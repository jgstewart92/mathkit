/* Test meshes, independent of the headers: an icosphere, a subdivided box and a torus, each appended to a tm_mesh
   (V: 3 doubles a vertex, T: 3 indices a face, grown with realloc; tm_free frees both). tm_volume is the signed volume
   (positive for an outward closed mesh). tm_write writes an OBJ file into tm_dir; with tm_dir NULL it does nothing. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { double *V; int *T; int nv, nt; } tm_mesh;
[[maybe_unused]] static void tm_free(tm_mesh *m) { free(m->V), free(m->T), memset(m, 0, sizeof *m); }
[[maybe_unused]] static int tm_vert(tm_mesh *m, double x, double y, double z)
{
    m->V = realloc(m->V, sizeof(double) * 3 * (size_t)(m->nv + 1));
    m->V[3 * m->nv] = x, m->V[3 * m->nv + 1] = y, m->V[3 * m->nv + 2] = z;
    return m->nv++;
}
[[maybe_unused]] static void tm_tri(tm_mesh *m, int a, int b, int c)
{
    m->T = realloc(m->T, sizeof(int) * 3 * (size_t)(m->nt + 1));
    m->T[3 * m->nt] = a, m->T[3 * m->nt + 1] = b, m->T[3 * m->nt + 2] = c, m->nt++;
}
/* an icosphere of radius R about c, S subdivisions (20 * 4^S faces), outward, appended to m */
[[maybe_unused]] static void tm_sphere(tm_mesh *m, double R, const double *c, int S)
{
    const double t = (1 + sqrt(5.0)) / 2;
    const double v0[12][3] = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
    const int f0[20][3] = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8}, {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
    int nv = 12, nt = 20;
    double *V = malloc(sizeof(double) * 3 * 12);
    int *T = malloc(sizeof f0);
    for (int i = 0; i < 12; i++) {
        const double l = sqrt(v0[i][0] * v0[i][0] + v0[i][1] * v0[i][1] + v0[i][2] * v0[i][2]);
        for (int k = 0; k < 3; k++) V[3 * i + k] = v0[i][k] / l;
    }
    memcpy(T, f0, sizeof f0);
    for (int s = 0; s < S; s++) {
        int *key = malloc(sizeof(int) * 6 * (size_t)nt), *val = malloc(sizeof(int) * 3 * (size_t)nt), *head = malloc(sizeof(int) * (size_t)nv), *nxt = malloc(sizeof(int) * 3 * (size_t)nt), nk = 0;
        int *T2 = malloc(sizeof(int) * 12 * (size_t)nt);
        for (int i = 0; i < nv; i++) head[i] = -1;
        V = realloc(V, sizeof(double) * 3 * (size_t)(nv + 3 * nt));
        for (int q = 0; q < nt; q++) {
            int mm[3];
            for (int e = 0; e < 3; e++) {
                const int a = T[3 * q + e], b = T[3 * q + (e + 1) % 3], lo = a < b ? a : b, hi = a < b ? b : a;
                int f = -1;
                for (int z = head[lo]; z >= 0; z = nxt[z])
                    if (key[2 * z + 1] == hi) { f = val[z]; break; }
                if (f < 0) {
                    double p[3], l = 0;
                    for (int k = 0; k < 3; k++) p[k] = V[3 * a + k] + V[3 * b + k], l += p[k] * p[k];
                    l = sqrt(l);
                    for (int k = 0; k < 3; k++) V[3 * nv + k] = p[k] / l;
                    key[2 * nk] = lo, key[2 * nk + 1] = hi, val[nk] = f = nv++, nxt[nk] = head[lo], head[lo] = nk, nk++;
                }
                mm[e] = f;
            }
            const int a = T[3 * q], b = T[3 * q + 1], cc = T[3 * q + 2], nn[12] = {a, mm[0], mm[2], b, mm[1], mm[0], cc, mm[2], mm[1], mm[0], mm[1], mm[2]};
            memcpy(T2 + 12 * q, nn, sizeof nn);
        }
        free(T), free(key), free(val), free(head), free(nxt);
        T = T2, nt *= 4;
    }
    const int base = m->nv;
    for (int i = 0; i < nv; i++) tm_vert(m, c[0] + R * V[3 * i], c[1] + R * V[3 * i + 1], c[2] + R * V[3 * i + 2]);
    for (int q = 0; q < nt; q++) tm_tri(m, base + T[3 * q], base + T[3 * q + 1], base + T[3 * q + 2]);
    free(V), free(T);
}
/* a box lo .. lo + s, n x n a face, outward, welded within itself, appended to m */
[[maybe_unused]] static void tm_box(tm_mesh *m, const double *lo, const double *s, int n)
{
    const int side = n + 1, base = m->nv;
    int *id = malloc(sizeof(int) * (size_t)(side * side * side));
    for (int i = 0; i < side * side * side; i++) id[i] = -1;
    for (int f = 0; f < 6; f++)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                int q[4];
                const int di[4] = {0, 1, 1, 0}, dj[4] = {0, 0, 1, 1};
                for (int c = 0; c < 4; c++) {
                    const int I = i + di[c], J = j + dj[c], e = (f % 2) * n, ax = f / 2;
                    const int X = ax == 0 ? e : ax == 1 ? J : I, Y = ax == 0 ? I : ax == 1 ? e : J, Z = ax == 0 ? J : ax == 1 ? I : e;
                    int *p = &id[(Z * side + Y) * side + X];
                    if (*p < 0) *p = tm_vert(m, lo[0] + s[0] * X / n, lo[1] + s[1] * Y / n, lo[2] + s[2] * Z / n) - base;
                    q[c] = base + *p;
                }
                if (f % 2 == 0) tm_tri(m, q[0], q[2], q[1]), tm_tri(m, q[0], q[3], q[2]);
                else tm_tri(m, q[0], q[1], q[2]), tm_tri(m, q[0], q[2], q[3]);
            }
    free(id);
}
/* a torus about z: tube radius r at radius R, nu x nw quads of two faces, outward, appended to m */
[[maybe_unused]] static void tm_torus(tm_mesh *m, double R, double r, int nu, int nw)
{
    const int base = m->nv;
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nw; j++) {
            const double u = 2 * 3.14159265358979323846 * i / nu, w = 2 * 3.14159265358979323846 * j / nw;
            tm_vert(m, (R + r * cos(w)) * cos(u), (R + r * cos(w)) * sin(u), r * sin(w));
        }
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nw; j++) {
            const int a = base + i * nw + j, b = base + ((i + 1) % nu) * nw + j, c = base + ((i + 1) % nu) * nw + (j + 1) % nw, d = base + i * nw + (j + 1) % nw;
            tm_tri(m, a, b, c), tm_tri(m, a, c, d);
        }
}
[[maybe_unused]] static double tm_volume(const double *V, const int *T, int nt)
{
    double v = 0;
    for (int t = 0; t < nt; t++) {
        const double *a = V + 3 * T[3 * t], *b = V + 3 * T[3 * t + 1], *c = V + 3 * T[3 * t + 2];
        v += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
    }
    return v;
}
static const char *tm_dir;
[[maybe_unused]] static void tm_write(const char *name, const double *V, int nv, const int *T, int nt, const unsigned char *alive)
{
    if (!tm_dir) return;
    char p[1024];
    snprintf(p, sizeof p, "%s/%s.obj", tm_dir, name);
    FILE *f = fopen(p, "w");
    if (!f) return;
    int *map = malloc(sizeof(int) * (size_t)(nv ? nv : 1)), k = 0;
    for (int i = 0; i < nv; i++) map[i] = !alive || alive[i] ? k++ : -1;
    for (int i = 0; i < nv; i++)
        if (map[i] >= 0) fprintf(f, "v %.17g %.17g %.17g\n", V[3 * i], V[3 * i + 1], V[3 * i + 2]);
    for (int t = 0; t < nt; t++) fprintf(f, "f %d %d %d\n", map[T[3 * t]] + 1, map[T[3 * t + 1]] + 1, map[T[3 * t + 2]] + 1);
    fclose(f), free(map);
}
