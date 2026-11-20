#define _POSIX_C_SOURCE 199309L
/* speed: CG on the 64^3 Laplacian (none, Jacobi, IC(0)); BVH queries per second against brute force. `make bench` */
#include "../include/mk_sparse.h"
#include "../include/mk_geom.h"
#include <stdio.h>
#include <time.h>

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static unsigned long long st = 88172645463325252ULL;
static double rnd(void) { st ^= st << 13, st ^= st >> 7, st ^= st << 17; return (double)(st >> 11) * 0x1.0p-52 - 1.0; }

int main(void)
{
    enum { M = 64 };
    int n = M * M * M, k = 0;
    int *ri = malloc(sizeof(int) * 14 * (size_t)n), *ci = ri + 7 * (size_t)n;
    double *v = malloc(sizeof(double) * 7 * (size_t)n), *b = malloc(sizeof(double) * (size_t)n), *x = malloc(sizeof(double) * (size_t)n);
    for (int z = 0; z < M; z++)
        for (int y = 0; y < M; y++)
            for (int xx = 0; xx < M; xx++) {
                int r = (z * M + y) * M + xx;
                ri[k] = r, ci[k] = r, v[k++] = 6;
                if (xx > 0) ri[k] = r, ci[k] = r - 1, v[k++] = -1;
                if (xx < M - 1) ri[k] = r, ci[k] = r + 1, v[k++] = -1;
                if (y > 0) ri[k] = r, ci[k] = r - M, v[k++] = -1;
                if (y < M - 1) ri[k] = r, ci[k] = r + M, v[k++] = -1;
                if (z > 0) ri[k] = r, ci[k] = r - M * M, v[k++] = -1;
                if (z < M - 1) ri[k] = r, ci[k] = r + M * M, v[k++] = -1;
            }
    mk_csr A;
    mk_csr_from_triplets(&A, n, k, ri, ci, v);
    for (int i = 0; i < n; i++) b[i] = rnd();
    const char *pn[3] = {"none", "Jacobi", "IC(0)"};
    int pk[3] = {MK_PC_NONE, MK_PC_JACOBI, MK_PC_IC0};
    printf("CG, 64^3 Laplacian (%d unknowns, %d entries), tol 1e-8, one thread:\n", n, A.rowptr[n]);
    for (int p = 0; p < 3; p++) {
        for (int i = 0; i < n; i++) x[i] = 0;
        int it;
        double rr, t0 = now();
        int r = mk_cg(&A, b, x, pk[p], 1e-8, 5000, &it, &rr);
        double dt = now() - t0;
        printf("  %-7s %s  %4d iterations  %7.3f s  %6.2f ms an iteration  residual %.2g\n", pn[p], r ? "FAIL" : "ok  ", it, dt, 1e3 * dt / it, rr);
    }
    mk_csr_free(&A), free(ri), free(v), free(b), free(x);

    enum { NU = 512, NW = 256 };
    int nv = NU * NW, nt = 2 * NU * NW;
    double *V = malloc(sizeof(double) * 3 * (size_t)nv);
    int *T = malloc(sizeof(int) * 3 * (size_t)nt);
    for (int i = 0; i < NU; i++)
        for (int j = 0; j < NW; j++) {
            double u = 6.283185307179586 * i / NU, w = 6.283185307179586 * j / NW;
            double *q = V + 3 * (size_t)(i * NW + j);
            q[0] = (2 + 0.7 * cos(w)) * cos(u), q[1] = (2 + 0.7 * cos(w)) * sin(u), q[2] = 0.7 * sin(w);
        }
    for (int i = 0, t = 0; i < NU; i++)
        for (int j = 0; j < NW; j++) {
            int a = i * NW + j, bb = ((i + 1) % NU) * NW + j, c = ((i + 1) % NU) * NW + (j + 1) % NW, d = i * NW + (j + 1) % NW;
            T[3 * t] = a, T[3 * t + 1] = bb, T[3 * t + 2] = c, t++;
            T[3 * t] = a, T[3 * t + 1] = c, T[3 * t + 2] = d, t++;
        }
    mk_bvh B;
    double t0 = now();
    mk_bvh_build(&B, V, T, nt);
    double tb = now() - t0;
    printf("BVH on a %d-triangle torus: built in %.3f s, %d nodes\n", nt, tb, B.nnode);
    enum { Q = 200000, QB = 200 };
    double sink = 0;
    t0 = now();
    for (int q = 0; q < Q; q++) {
        double p[3] = {3 * rnd(), 3 * rnd(), 3 * rnd()}, c[3];
        int tr;
        sink += mk_bvh_nearest(&B, p, INFINITY, &tr, c);
    }
    double tn = now() - t0;
    t0 = now();
    for (int q = 0; q < QB; q++) {
        double p[3] = {3 * rnd(), 3 * rnd(), 3 * rnd()}, best = INFINITY;
        for (int s = 0; s < nt; s++) {
            double w[3], c[3], d2 = mk_closest_pt_tri(p, V + 3 * (size_t)T[3 * s], V + 3 * (size_t)T[3 * s + 1], V + 3 * (size_t)T[3 * s + 2], w, c);
            if (d2 < best) best = d2;
        }
        sink += best;
    }
    double tnb = now() - t0;
    t0 = now();
    int hits = 0;
    for (int q = 0; q < Q; q++) {
        double o[3] = {3 * rnd(), 3 * rnd(), 3 * rnd()}, d[3] = {rnd(), rnd(), rnd()}, th;
        hits += mk_bvh_first_hit(&B, o, d, 0, INFINITY, &th, NULL) >= 0;
    }
    double tr = now() - t0;
    t0 = now();
    for (int q = 0; q < QB; q++) {
        double o[3] = {3 * rnd(), 3 * rnd(), 3 * rnd()}, d[3] = {rnd(), rnd(), rnd()}, bt = INFINITY, th;
        mk_ray r = mk_ray_make(o, d);
        for (int s = 0; s < nt; s++)
            if (mk_ray_tri(&r, V + 3 * (size_t)T[3 * s], V + 3 * (size_t)T[3 * s + 1], V + 3 * (size_t)T[3 * s + 2], 0, bt, &th, NULL)) bt = th;
        sink += bt < INFINITY;
    }
    double trb = now() - t0;
    printf("  nearest point: %.0f a second with the tree, %.0f brute force: %.0f times\n", Q / tn, QB / tnb, (Q / tn) / (QB / tnb));
    printf("  first hit:     %.0f a second with the tree, %.0f brute force: %.0f times (%d of %d rays hit)\n", Q / tr, QB / trb, (Q / tr) / (QB / trb), hits, Q);
    mk_bvh_free(&B), free(V), free(T);
    return sink < 0;
}
