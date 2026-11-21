/* mathkit's threads (mk_par.h). Built with -fopenmp, every threaded routine runs at 1, 2 and 12 threads and must give
   the same bytes. Each result's hash is printed on a BITS line, and check.sh requires a build without OpenMP to print
   the same lines. Also: the batch queries against their single calls, the sdf cases reaching the stages they are meant
   to (odd rows, near nodes, several tiles), CG's block sums against one serial sum (written out here), the block
   arithmetic near LONG_MAX, and mk_welch from two threads on one plan (a row, no BITS line). */
#include "../include/mk_sparse.h"
#include "../include/mk_sdf.h"
#include "../include/mk_hhd.h"
#include "../include/mk_spec.h"
#include "t.h"
#include <limits.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#define NRUN 3
#else
#define NRUN 1
#endif
static const int nthr[3] = {1, 2, 12};
static void threads(int t)
{
#ifdef _OPENMP
    omp_set_num_threads(t);
#else
    (void)t;
#endif
}
static uint64_t fnv(const void *p, size_t n, uint64_t h)
{
    const unsigned char *b = p;
    for (size_t i = 0; i < n; i++) h = (h ^ b[i]) * 1099511628211ULL;
    return h;
}
#define FNV0 1469598103934665603ULL
/* the row for one routine: the same hash at every thread count; the BITS line for check.sh's comparison of builds */
static void same(const char *what, const uint64_t *h)
{
#ifdef _OPENMP
    CHECK(h[0] == h[1] && h[1] == h[2], "%s: the same bytes at 1, 2 and 12 threads (%016llx %016llx %016llx)", what, (unsigned long long)h[0],
          (unsigned long long)h[1], (unsigned long long)h[2]);
#endif
    printf("BITS  %-44s %016llx\n", what, (unsigned long long)h[0]);
}

/* the 7-point Laplacian on an m^3 grid: Dirichlet, periodic, or with upwind convection c along x (nonsymmetric) */
static void lap3(mk_csr *A, int m, int periodic, double c)
{
    const int n = m * m * m;
    int *ri = malloc(sizeof(int) * 14 * (size_t)n), *ci = ri + 7 * (size_t)n, k = 0;
    double *v = malloc(sizeof(double) * 7 * (size_t)n);
    for (int z = 0; z < m; z++)
        for (int y = 0; y < m; y++)
            for (int x = 0; x < m; x++) {
                const int r = (z * m + y) * m + x, nb[6][3] = {{x - 1, y, z}, {x + 1, y, z}, {x, y - 1, z}, {x, y + 1, z}, {x, y, z - 1}, {x, y, z + 1}};
                ri[k] = r, ci[k] = r, v[k++] = 6 + c;
                for (int e = 0; e < 6; e++) {
                    int a = nb[e][0], b = nb[e][1], d = nb[e][2];
                    if (periodic) a = (a + m) % m, b = (b + m) % m, d = (d + m) % m;
                    else if (a < 0 || b < 0 || d < 0 || a >= m || b >= m || d >= m) continue;
                    ri[k] = r, ci[k] = (d * m + b) * m + a, v[k++] = e == 0 ? -1 - c : -1;
                }
            }
    mk_csr_from_triplets(A, n, k, ri, ci, v);
    free(ri), free(v);
}
/* CG without a preconditioner, mk_cg's steps with every dot product one serial sum: the reference for the block sums */
static double dot1(const double *a, const double *b, int n)
{
    double s = 0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}
static int cg140(const mk_csr *A, const double *b, double *x, double tol, int maxit, int *iters)
{
    const int n = A->n;
    double *r = malloc(sizeof(double) * 4 * (size_t)n), *p = r + n, *q = p + n, *u = q + n, big = 0;
    int e = 0, it = 0, ret = 1;
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int k = A->rowptr[i]; k < A->rowptr[i + 1]; k++) s += A->val[k] * x[A->col[k]];
        r[i] = b[i] - s, big = fmax(big, fabs(b[i]));
    }
    if (big > 0 && isfinite(big)) frexp(big, &e);
    for (int i = 0; i < n; i++) r[i] = ldexp(r[i], -e), u[i] = ldexp(b[i], -e);
    const double sc = ldexp(1.0, e);
    double bn = sqrt(dot1(u, u, n)), rn = sqrt(dot1(r, r, n));
    memcpy(p, r, sizeof *p * (size_t)n);
    double rz = dot1(r, r, n), rz0 = rz;
    if (rn / bn <= tol) ret = 0;
    while (ret == 1 && it < maxit) {
        for (int i = 0; i < n; i++) {
            double s = 0;
            for (int k = A->rowptr[i]; k < A->rowptr[i + 1]; k++) s += A->val[k] * p[A->col[k]];
            q[i] = s;
        }
        const double pq = dot1(p, q, n), a = rz / pq;
        for (int i = 0; i < n; i++) x[i] += (a * sc) * p[i], r[i] -= a * q[i];
        rn = sqrt(dot1(r, r, n)), it++;
        if (rn / bn <= tol) { ret = 0; break; }
        const double rz2 = dot1(r, r, n);
        if (!(rz2 > 0) || (rz2 < DBL_MIN && rz2 < rz0 * 0x1p-100)) break;
        const double beta = rz2 / rz;
        for (int i = 0; i < n; i++) p[i] = r[i] + beta * p[i];
        rz = rz2;
    }
    *iters = it;
    free(r);
    return ret;
}

typedef struct { double *V; int nv, cap; int *T; int nt, tcap; } msh;
static int addv(msh *m, double x, double y, double z)
{
    if (m->nv == m->cap) m->V = realloc(m->V, sizeof(double) * 3 * (size_t)(m->cap = m->cap ? 2 * m->cap : 1024));
    m->V[3 * m->nv] = x, m->V[3 * m->nv + 1] = y, m->V[3 * m->nv + 2] = z;
    return m->nv++;
}
static void addt(msh *m, int a, int b, int c)
{
    if (m->nt == m->tcap) m->T = realloc(m->T, sizeof(int) * 3 * (size_t)(m->tcap = m->tcap ? 2 * m->tcap : 1024));
    m->T[3 * m->nt] = a, m->T[3 * m->nt + 1] = b, m->T[3 * m->nt + 2] = c, m->nt++;
}
static void torus(msh *m, int nu, int nw, double R, double r)
{
    memset(m, 0, sizeof *m);
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nw; j++) {
            const double u = 2 * T_PI * i / nu, w = 2 * T_PI * j / nw;
            addv(m, (R + r * cos(w)) * cos(u), (R + r * cos(w)) * sin(u), r * sin(w));
        }
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nw; j++) {
            const int a = i * nw + j, b = ((i + 1) % nu) * nw + j, c = ((i + 1) % nu) * nw + (j + 1) % nw, d = i * nw + (j + 1) % nw;
            addt(m, a, b, c), addt(m, a, c, d);
        }
}
/* a UV sphere of radius 1, outward, with the faces whose centre has x > cut left out (a hole facing +x; cut >= 1: closed) */
static void sphere(msh *m, int nu, int nw, double cut)
{
    memset(m, 0, sizeof *m);
    addv(m, 0, 0, 1);
    for (int i = 1; i < nw; i++)
        for (int j = 0; j < nu; j++) {
            const double th = T_PI * i / nw, ph = 2 * T_PI * j / nu;
            addv(m, sin(th) * cos(ph), sin(th) * sin(ph), cos(th));
        }
    const int south = addv(m, 0, 0, -1);
    int tri[3 * 8192], k = 0;
    for (int j = 0; j < nu; j++) tri[k++] = 0, tri[k++] = 1 + j, tri[k++] = 1 + (j + 1) % nu;
    for (int i = 0; i < nw - 2; i++)
        for (int j = 0; j < nu; j++) {
            const int a = 1 + i * nu + j, b = 1 + (i + 1) * nu + j, c = 1 + (i + 1) * nu + (j + 1) % nu, d = 1 + i * nu + (j + 1) % nu;
            tri[k++] = a, tri[k++] = b, tri[k++] = c, tri[k++] = a, tri[k++] = c, tri[k++] = d;
        }
    for (int j = 0; j < nu; j++) tri[k++] = 1 + (nw - 2) * nu + j, tri[k++] = south, tri[k++] = 1 + (nw - 2) * nu + (j + 1) % nu;
    for (int t = 0; t < k; t += 3)
        if ((m->V[3 * tri[t]] + m->V[3 * tri[t + 1]] + m->V[3 * tri[t + 2]]) / 3 <= cut) addt(m, tri[t], tri[t + 1], tri[t + 2]);
}
/* an axis-aligned box [0, s]^3, outward, two triangles a face */
static void box(msh *m, double s)
{
    memset(m, 0, sizeof *m);
    for (int i = 0; i < 8; i++) addv(m, i & 1 ? s : 0, i & 2 ? s : 0, i & 4 ? s : 0);
    static const int f[36] = {0, 2, 3, 0, 3, 1, 4, 5, 7, 4, 7, 6, 0, 1, 5, 0, 5, 4, 2, 6, 7, 2, 7, 3, 0, 4, 6, 0, 6, 2, 1, 3, 7, 1, 7, 5};
    for (int t = 0; t < 36; t += 3) addt(m, f[t], f[t + 1], f[t + 2]);
}
static uint64_t sdfhash(const msh *m, double h, double margin, int flags, mk_sdf *keep)
{
    mk_sdf g;
    if (mk_sdf_mesh(m->V, m->nv, m->T, m->nt, h, margin, flags, &g)) return 0;
    uint64_t s = fnv(g.d, sizeof(float) * (size_t)g.n[0] * g.n[1] * g.n[2], FNV0);
    const long c[4] = {g.odd_rows, g.reray, g.near, g.filled};
    s = fnv(c, sizeof c, fnv(g.n, sizeof g.n, fnv(g.lo, sizeof g.lo, s)));
    if (keep) *keep = g;
    else mk_sdf_free(&g);
    return s;
}

int main(void)
{
    uint64_t hs[3];
    /* sparse: 24^3 = 13824 unknowns, four blocks of rows */
    enum { M = 24 };
    const int n = M * M * M;
    mk_csr A, Pd, C;
    lap3(&A, M, 0, 0), lap3(&Pd, M, 1, 0), lap3(&C, M, 0, 2.5);
    double *b = malloc(sizeof(double) * 4 * (size_t)n), *x = b + n, *y = x + n, *bp = y + n, mean = 0;
    for (int i = 0; i < n; i++) b[i] = t_rand(), bp[i] = t_rand(), mean += bp[i];
    for (int i = 0; i < n; i++) bp[i] -= mean / n;
    for (int r = 0; r < NRUN; r++) { /* every output cleared before each run, so a left-out entry cannot keep the last run's value */
        for (int i = 0; i < n; i++) y[i] = NAN;
        threads(nthr[r]), mk_csr_matvec(&A, b, y), hs[r] = fnv(y, sizeof(double) * (size_t)n, FNV0);
    }
    same("mk_csr_matvec (24^3 Laplacian)", hs);
    static const struct { const char *name; int bicg, pc, periodic; } sv[] = {
        {"mk_cg none", 0, MK_PC_NONE, 0}, {"mk_cg Jacobi", 0, MK_PC_JACOBI, 0}, {"mk_cg block-Jacobi", 0, MK_PC_BJACOBI3, 0},
        {"mk_cg IC(0)", 0, MK_PC_IC0, 0}, {"mk_cg mean-free (periodic)", 0, MK_PC_NONE | MK_PC_MEANFREE, 1},
        {"mk_bicgstab none (convection-diffusion)", 1, MK_PC_NONE, 0}, {"mk_bicgstab Jacobi (convection-diffusion)", 1, MK_PC_JACOBI, 0}};
    for (size_t s = 0; s < sizeof sv / sizeof *sv; s++) {
        for (int r = 0; r < NRUN; r++) {
            threads(nthr[r]);
            int it = -1;
            double rr = -1;
            memset(x, 0, sizeof(double) * (size_t)n);
            const mk_csr *K = sv[s].bicg ? &C : sv[s].periodic ? &Pd : &A;
            const double *rhs = sv[s].periodic ? bp : b;
            const int ret = sv[s].bicg ? mk_bicgstab(K, rhs, x, sv[s].pc, 1e-10, 3000, &it, &rr) : mk_cg(K, rhs, x, sv[s].pc, 1e-10, 3000, &it, &rr);
            const double tail[3] = {ret, it, rr};
            hs[r] = fnv(tail, sizeof tail, fnv(x, sizeof(double) * (size_t)n, FNV0));
            if (r == 0 && ret) printf("note: %s returned %d\n", sv[s].name, ret);
        }
        same(sv[s].name, hs);
    }
    /* CG with block sums against CG with one serial sum: close and in about the same iterations, but not the same bits,
       which shows the block sums are in use */
    threads(nthr[NRUN - 1]);
    int it5, it4;
    memset(x, 0, sizeof(double) * (size_t)n), memset(y, 0, sizeof(double) * (size_t)n);
    const int r5 = mk_cg(&A, b, x, MK_PC_NONE, 1e-10, 3000, &it5, NULL), r4 = cg140(&A, b, y, 1e-10, 3000, &it4);
    double dmax = 0, xmax = 0;
    int same_bits = 1;
    for (int i = 0; i < n; i++) dmax = fmax(dmax, fabs(x[i] - y[i])), xmax = fmax(xmax, fabs(y[i])), same_bits &= x[i] == y[i];
    CHECK(!r5 && !r4 && dmax <= 1e-9 * xmax && abs(it5 - it4) <= 2 && !same_bits,
          "CG, 24^3 Laplacian, tol 1e-10: block sums against one serial sum: x differs by %.2g of its largest (%s), %d and %d iterations", dmax / xmax,
          same_bits ? "bit for bit the same" : "bits differ", it5, it4);
    mk_csr_free(&A), mk_csr_free(&Pd), mk_csr_free(&C), free(b);

    /* sdf: a torus (several tiles), an open sphere (odd rows: the winding number), a box with nodes 1e-5 h off its faces (near nodes) */
    msh tor, open, bx;
    torus(&tor, 96, 48, 0.6, 0.25), sphere(&open, 48, 24, 0.97), box(&bx, (10 + 1e-5) * 0.05);
    const struct { const char *name; const msh *m; double h, margin; int flags; } sc[] = {
        {"mk_sdf_mesh torus", &tor, 0.04, 0.3, 0}, {"mk_sdf_mesh torus MK_SDF_EXACT", &tor, 0.04, 0.3, MK_SDF_EXACT},
        {"mk_sdf_mesh torus MK_SDF_FILL | MK_SDF_LATTICE", &tor, 0.04, 0.3, MK_SDF_FILL | MK_SDF_LATTICE}, {"mk_sdf_mesh open sphere (odd rows)", &open, 0.05, 0.2, 0},
        {"mk_sdf_mesh box (near nodes)", &bx, 0.05, 0.25, 0}};
    mk_sdf keep[5];
    for (size_t s = 0; s < sizeof sc / sizeof *sc; s++) {
        for (int r = 0; r < NRUN; r++) threads(nthr[r]), hs[r] = sdfhash(sc[s].m, sc[s].h, sc[s].margin, sc[s].flags, r ? NULL : &keep[s]);
        same(sc[s].name, hs);
    }
    CHECK(keep[0].n[1] > 3 * MK_SDF_TILE && keep[0].n[2] > 2 * MK_SDF_TILE && keep[3].odd_rows > 0 && keep[4].near > 0,
          "the sdf cases reach what they test: the torus grid %d x %d x %d (%d x %d tiles), the open sphere %ld odd rows, the box %ld near nodes", keep[0].n[0],
          keep[0].n[1], keep[0].n[2], (keep[0].n[1] + MK_SDF_TILE - 1) / MK_SDF_TILE, (keep[0].n[2] + MK_SDF_TILE - 1) / MK_SDF_TILE, keep[3].odd_rows, keep[4].near);
    for (int s = 0; s < 5; s++) mk_sdf_free(&keep[s]);

    /* batch queries on a torus: the same bytes at every thread count, and the single calls' answers */
    msh tq;
    torus(&tq, 128, 64, 0.6, 0.25);
    mk_bvh B;
    mk_bvh_build(&B, tq.V, tq.T, tq.nt);
    enum { NP = 20000, NW = 400 };
    double *P = malloc(sizeof(double) * 3 * NP), *d2 = malloc(sizeof(double) * NP), *q = malloc(sizeof(double) * 3 * NP), *w = malloc(sizeof(double) * NP);
    int *tri = malloc(sizeof(int) * NP), *in = malloc(sizeof(int) * NP);
    for (int i = 0; i < NP; i++) P[3 * i] = 1.1 * t_rand(), P[3 * i + 1] = 1.1 * t_rand(), P[3 * i + 2] = 0.4 * t_rand();
    for (int r = 0; r < NRUN; r++) {
        threads(nthr[r]);
        for (int i = 0; i < NP; i++) d2[i] = NAN, tri[i] = -2, q[3 * i] = q[3 * i + 1] = q[3 * i + 2] = NAN;
        mk_bvh_nearest_many(&B, P, NP, 0.04, d2, tri, q);
        hs[r] = fnv(q, sizeof(double) * 3 * NP, fnv(tri, sizeof(int) * NP, fnv(d2, sizeof(double) * NP, FNV0)));
    }
    same("mk_bvh_nearest_many (20000 points, maxd2 0.04)", hs);
    int bad = 0, none = 0;
    for (int i = 0; i < NP; i++) {
        double c[3] = {NAN, NAN, NAN};
        int t;
        const double e = mk_bvh_nearest(&B, P + 3 * i, 0.04, &t, c);
        none += t < 0;
        bad += memcmp(&e, &d2[i], sizeof e) || t != tri[i] || memcmp(c, q + 3 * i, sizeof c);
    }
    CHECK(bad == 0 && none > 0 && none < NP, "mk_bvh_nearest_many: every point's distance, triangle and point the single call's (%d differ; %d beyond the bound)", bad, none);
    for (int r = 0; r < NRUN; r++) {
        for (int i = 0; i < NP; i++) in[i] = -1;
        threads(nthr[r]), mk_bvh_inside_many(&B, P, NP, in), hs[r] = fnv(in, sizeof(int) * NP, FNV0);
    }
    same("mk_bvh_inside_many (20000 points)", hs);
    int nin = 0;
    bad = 0;
    for (int i = 0; i < NP; i++) bad += in[i] != mk_bvh_inside(&B, P + 3 * i), nin += in[i];
    CHECK(bad == 0 && nin > 0 && nin < NP, "mk_bvh_inside_many: the single calls' answers (%d differ; %d inside)", bad, nin);
    for (int r = 0; r < NRUN; r++) {
        for (int i = 0; i < NW; i++) w[i] = NAN;
        threads(nthr[r]), mk_winding_number_many(tq.V, tq.T, tq.nt, P, NW, w), hs[r] = fnv(w, sizeof(double) * NW, FNV0);
    }
    same("mk_winding_number_many (400 points)", hs);
    bad = 0;
    for (int i = 0; i < NW; i++) {
        const double e = mk_winding_number(tq.V, tq.T, tq.nt, P + 3 * i);
        bad += memcmp(&e, &w[i], sizeof e) != 0;
    }
    CHECK(bad == 0, "mk_winding_number_many: the single calls' values bit for bit (%d differ)", bad);
    mk_bvh_free(&B), free(P), free(d2), free(q), free(w), free(tri), free(in);

    { /* mk_hhd_split: its transforms by column blocks and its solves by rows, spread over the threads */
        enum { HX = 256, HY = 144 };
        mk_hhd H;
        double *hu = malloc(sizeof(double) * HX * HY), *hv = malloc(sizeof(double) * HX * HY);
        for (int k = 0; k < HX * HY; k++) hu[k] = t_rand(), hv[k] = t_rand();
        int hok = !mk_hhd_init(&H, HX, HY);
        for (int r = 0; r < NRUN && hok; r++) {
            for (int k = 0; k < HX * HY; k++) H.gx[k] = H.rx[k] = H.hx[k] = NAN;
            threads(nthr[r]), hok &= !mk_hhd_split(&H, hu, hv);
            hs[r] = fnv(H.gx, sizeof(double) * (6 * HX * HY + 2 * (HX + 1) * (HY + 1)), FNV0); /* the six parts and both potentials: one block */
        }
        CHECK(hok, "mk_hhd_split 256 x 144 ran at every thread count");
        same("mk_hhd_split (256 x 144)", hs);
        mk_hhd_free(&H), free(hu), free(hv);
    }

    { /* the block arithmetic for n near LONG_MAX, where the plain ceiling (n + k - 1) / k would overflow a long
         (n = 9223371998217895936 among them). The test's own ceilings are in unsigned long, which cannot wrap there */
        volatile long nv[] = {LONG_MAX, LONG_MAX - 1, LONG_MAX - 255, LONG_MAX - 256, 9223371998217895936L, LONG_MAX / 2 + 1, (1L << 62) + 1};
        int ok = 1, nok = 0;
        for (size_t c = 0; c < sizeof nv / sizeof *nv; c++) {
            const long n = nv[c], b = mk_par_blk(n), nb = mk_par_nblk(n);
            unsigned long ub = ((unsigned long)n + MK_PAR_NBLK - 1) / MK_PAR_NBLK;
            if (ub < MK_PAR_BLK) ub = MK_PAR_BLK;
            const unsigned long unb = ((unsigned long)n + ub - 1) / ub;
            int k_ok = b > 0 && (unsigned long)b == ub && nb >= 1 && nb <= MK_PAR_NBLK && (unsigned long)nb == unb;
            long at = 0;
            for (long k = 0; k < nb && k < MK_PAR_NBLK && k_ok; k++) {
                long i0, i1;
                mk_par_range(n, k, &i0, &i1);
                k_ok &= i0 == at && i1 > i0 && i1 - i0 <= b, at = i1;
            }
            k_ok &= at == n;
            ok &= k_ok, nok += k_ok;
        }
        CHECK(ok, "mk_par near LONG_MAX: %d of %d counts from 2^62 + 1 to LONG_MAX cut into at most %d blocks of the exact ceiling's length, tiling [0, n)", nok,
              (int)(sizeof nv / sizeof *nv), MK_PAR_NBLK);
    }
#ifdef _OPENMP
    { /* two threads share one const mk_spec plan, each with its own signal, 100 rounds: every answer equals the
         serial one, so mk_stft_frame keeps its windowed frame out of the plan */
        enum { NS = 1 << 15, NF = 1024, NB = NF / 2 + 1, RR = 100 };
        static double sx[2][NS], ref[2][NB];
        double *got = malloc(sizeof(double) * 2 * NB); /* not a fixed-size array: gcc -O3 cannot bound mk_welch's nb and warns */
        for (int k = 0; k < 2; k++)
            for (int i = 0; i < NS; i++) sx[k][i] = (k ? 3 : 1) * t_rand() + (k ? sin(0.3 * i) : 0);
        mk_spec sp;
        int bad = 0, nt2 = 0;
        if (!got || mk_spec_init(&sp, NF, MK_WIN_HANN, 0, 1000)) bad = -1;
        else {
            for (int k = 0; k < 2; k++) mk_welch(&sp, sx[k], NS, NF / 2, 1, MK_SPEC_DENSITY, ref[k]);
            for (int r = 0; r < RR; r++) {
#pragma omp parallel num_threads(2)
                {
                    const int k = omp_get_thread_num();
                    if (k == 1) nt2 = 1;
                    mk_welch(&sp, sx[k], NS, NF / 2, 1, MK_SPEC_DENSITY, got + k * NB);
                }
                for (int k = 0; k < 2; k++) bad += memcmp(got + k * NB, ref[k], sizeof ref[k]) != 0;
            }
            mk_spec_free(&sp);
        }
        CHECK(bad == 0 && nt2, "mk_welch from 2 threads on one plan, %d rounds: %d answers differ from the serial ones", RR, bad);
        free(got);
    }
#endif
    free(tor.V), free(tor.T), free(open.V), free(open.T), free(bx.V), free(bx.T), free(tq.V), free(tq.T);
    return wrong;
}
