/* Regression rows for the verification work's fuzz findings against 1.5.3, patched in 1.5.4 (its demos: scratch/the verification work/mk-153f/findings/demo_*.c). Each
   row is the finding's exact call and fails on 1.5.3's headers. The calls that hung (O-1, O-1b) or ran for most of a minute (D-2) run in
   a child process with a limit of T_LIM seconds of CPU, so that a hang is a FAIL row here and not the driver's timeout. Those run first:
   a child must not inherit a parent that has already used OpenMP. */
#define _GNU_SOURCE /* fork, setitimer */
#include <float.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>
/* G-5's allocations: the one numbered t_al_fail fails; while t_al_on is set a block freed twice is counted and not freed again */
static long t_al_n, t_al_fail = -1, t_al_twice;
static int t_al_on, t_al_ng;
static void *t_al_gone[256];
static void *t_al_took(void *p)
{
    if (t_al_on) for (int i = 0; i < t_al_ng; i++) if (t_al_gone[i] == p) t_al_gone[i] = t_al_gone[--t_al_ng];
    return p;
}
static void *t_al_malloc(size_t n) { return t_al_took(t_al_on && t_al_n++ == t_al_fail ? NULL : malloc(n)); }
static void *t_al_calloc(size_t a, size_t b) { return t_al_took(t_al_on && t_al_n++ == t_al_fail ? NULL : calloc(a, b)); }
static void *t_al_realloc(void *q, size_t n) { return t_al_took(t_al_on && t_al_n++ == t_al_fail ? NULL : realloc(q, n)); }
static void t_al_free(void *p)
{
    if (!p) return;
    if (t_al_on) {
        for (int i = 0; i < t_al_ng; i++) if (t_al_gone[i] == p) { t_al_twice++; return; }
        if (t_al_ng < 256) t_al_gone[t_al_ng++] = p;
    }
    free(p);
}
#define malloc(n) t_al_malloc(n)
#define calloc(a, b) t_al_calloc(a, b)
#define realloc(p, n) t_al_realloc(p, n)
#define free(p) t_al_free(p)
#include "../include/mk_graft.h"
#undef malloc
#undef calloc
#undef realloc
#undef free
#include "../include/mk_dec.h"
#include "../include/mk_ovl.h"
#include "../include/mk_sdf.h"
#include "t.h"

enum { T_LIM = 5 };
/* fn(io) in a child limited to sec seconds of CPU (and a wall limit, should it sleep); io comes back through a pipe. 0: it returned,
   1: over the limit, 2: it died otherwise */
static int t_limited(void (*fn)(void *), void *io, size_t nio, int sec)
{
    int pfd[2];
    fflush(stdout);
    if (pipe(pfd)) return 2;
    const pid_t pid = fork();
    if (pid < 0) return 2;
    if (pid == 0) {
        close(pfd[0]);
        const struct itimerval it = {{0, 0}, {sec, 0}};
        setitimer(ITIMER_VIRTUAL, &it, NULL), alarm(4 * (unsigned)sec + 10);
        fn(io);
        _exit(write(pfd[1], io, nio) == (ssize_t)nio ? 0 : 2);
    }
    close(pfd[1]);
    size_t got = 0;
    for (ssize_t n; got < nio && (n = read(pfd[0], (char *)io + got, nio - got)) > 0;) got += (size_t)n;
    close(pfd[0]);
    int st = 0;
    waitpid(pid, &st, 0);
    if (WIFEXITED(st) && WEXITSTATUS(st) == 0 && got == nio) return 0;
    return WIFSIGNALED(st) && (WTERMSIG(st) == SIGVTALRM || WTERMSIG(st) == SIGALRM) ? 1 : 2;
}

/* an axis-aligned box [lo, hi]: 8 corners, 12 outward faces */
static void box(double *V, int *T, const double *lo, const double *hi)
{
    for (int c = 0; c < 8; c++) V[3 * c] = c & 1 ? hi[0] : lo[0], V[3 * c + 1] = c & 2 ? hi[1] : lo[1], V[3 * c + 2] = c & 4 ? hi[2] : lo[2];
    static const int F[36] = {0, 2, 1, 1, 2, 3, 4, 5, 6, 5, 7, 6, 0, 1, 4, 1, 5, 4, 2, 6, 3, 3, 6, 7, 0, 4, 2, 2, 4, 6, 1, 3, 5, 3, 7, 5};
    memcpy(T, F, sizeof F);
}
typedef struct { double lo0[3], hi0[3], lo1[3], hi1[3], h; int axis, flags, e; mk_ovl_report r; } ovl_case;
static void ovl_run(void *p)
{
    ovl_case *c = p;
    double VA[24], VB[24];
    int TA[36], TB[36];
    box(VA, TA, c->lo0, c->hi0), box(VB, TB, c->lo1, c->hi1);
    c->e = mk_overlap_volume(VA, 8, TA, 12, VB, 8, TB, 12, c->h, c->axis, c->flags, &c->r);
}
/* two boxes; the same box twice with lo1 NULL */
static ovl_case ovl(const double *lo0, const double *hi0, const double *lo1, const double *hi1, double h, int flags)
{
    ovl_case c = {.h = h, .axis = 2, .flags = flags, .e = -9};
    memcpy(c.lo0, lo0, sizeof c.lo0), memcpy(c.hi0, hi0, sizeof c.hi0);
    memcpy(c.lo1, lo1 ? lo1 : lo0, sizeof c.lo1), memcpy(c.hi1, hi1 ? hi1 : hi0, sizeof c.hi1);
    ovl_run(&c);
    return c;
}
static int near(double a, double b, double tol) { return fabs(a / b - 1) < tol; }

typedef struct { double *V; int *T, nv, nt; } mesh;
static int dv(mesh *m, double x, double y, double z)
{
    m->V = realloc(m->V, sizeof(double) * 3 * (size_t)(m->nv + 1));
    m->V[3 * (size_t)m->nv] = x, m->V[3 * (size_t)m->nv + 1] = y, m->V[3 * (size_t)m->nv + 2] = z;
    return m->nv++;
}
static void dt(mesh *m, int a, int b, int c)
{
    m->T = realloc(m->T, sizeof(int) * 3 * (size_t)(m->nt + 1));
    m->T[3 * (size_t)m->nt] = a, m->T[3 * (size_t)m->nt + 1] = b, m->T[3 * (size_t)m->nt + 2] = c, m->nt++;
}
typedef struct { int rc, nv; } dec_out;
/* D-2: a fine grid (edge 1 mm, 1800 triangles) beside one large tilted triangle 128 cells a side */
static void d2_run(void *p)
{
    dec_out *o = p;
    mesh m = {0};
    const int n = 30;
    for (int j = 0; j <= n; j++) for (int i = 0; i <= n; i++) dv(&m, 1e-3 * i, 1e-3 * j, 0);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) { const int a = j * (n + 1) + i; dt(&m, a, a + 1, a + n + 2), dt(&m, a, a + n + 2, a + n + 1); }
    const double L = 128 * 1e-3;
    const int a = dv(&m, 0.05, 0, 0), b = dv(&m, 0.05 + L, L, 0), c = dv(&m, 0.05, L, L);
    dt(&m, a, b, c);
    mk_dec_opts op;
    mk_dec_defaults(&op, 1e-3), op.maxedge = 2e-3;
    unsigned char *al = malloc((size_t)m.nv);
    int *To = malloc(sizeof(int) * 3 * (size_t)m.nt), nto = 0;
    mk_dec_report r;
    o->rc = mk_decimate(m.V, m.nv, m.T, m.nt, &op, al, To, &nto, &r), o->nv = o->rc ? -1 : r.nv;
    free(al), free(To), free(m.V), free(m.T);
}
static int dec_call(mesh *m, const mk_dec_opts *op, mk_dec_report *r)
{
    unsigned char *al = malloc((size_t)m->nv);
    int *To = malloc(sizeof(int) * 3 * (size_t)m->nt), nto = 0;
    const int e = mk_decimate(m->V, m->nv, m->T, m->nt, op, al, To, &nto, r);
    free(al), free(To);
    return e;
}
static double psi_inf(const double *p, void *ctx) /* the plane z = 0, +inf past z = lim */
{
    const double lim = *(const double *)ctx;
    return p[2] > lim ? INFINITY : p[2];
}
static double plane_z(const double *p, void *ctx) { (void)ctx; return p[2]; }
static void sdf_cube(double *V, int *T, double s)
{
    for (int q = 0; q < 8; q++) V[3 * q] = q & 1 ? s : -s, V[3 * q + 1] = q & 2 ? s : -s, V[3 * q + 2] = q & 4 ? s : -s;
    static const int F[36] = {0, 2, 3, 0, 3, 1, 4, 5, 7, 4, 7, 6, 0, 1, 5, 0, 5, 4, 2, 6, 7, 2, 7, 3, 0, 4, 6, 0, 6, 2, 1, 3, 7, 1, 7, 5};
    memcpy(T, F, sizeof F);
}

int main(void)
{
    { /* O-1, O-1b: the lattice's double counters stopped moving past 2^53. The coordinates there round to 256 (O-1) and 2^-52 (O-1b), wider
         than the spacing h, so the lines bunch: the volume is the definition's within 7 % (272 lines for 256, 65 x 65 for 64 x 64) */
        const double lo[3] = {0x1p60, 0, 0}, hi[3] = {0x1p60 + 1024, 1024, 1024}, lo2[3] = {1, 1, 1}, hi2[3] = {1 + 0x1p-48, 1 + 0x1p-48, 1 + 0x1p-48};
        ovl_case c = {.h = 64, .axis = 2, .e = -9};
        memcpy(c.lo0, lo, sizeof lo), memcpy(c.hi0, hi, sizeof hi), memcpy(c.lo1, lo, sizeof lo), memcpy(c.hi1, hi, sizeof hi);
        int t = t_limited(ovl_run, &c, sizeof c, T_LIM);
        CHECK(t == 0 && c.e == 0 && near(c.r.both, 0x1p30, 0.07), "O-1: two boxes 1024 wide at x = 2^60, h = 64: %s, both %.6g (2^30 = %.6g)", t ? "no return in the limit" : "returned", c.r.both, 0x1p30);
        ovl_case d = {.h = 0x1p-54, .axis = 2, .e = -9};
        memcpy(d.lo0, lo2, sizeof lo2), memcpy(d.hi0, hi2, sizeof hi2), memcpy(d.lo1, lo2, sizeof lo2), memcpy(d.hi1, hi2, sizeof hi2);
        t = t_limited(ovl_run, &d, sizeof d, T_LIM);
        CHECK(t == 0 && d.e == 0 && near(d.r.both, 0x1p-144, 0.07), "O-1b: two boxes 2^-48 wide at (1, 1, 1), h = 2^-54: %s, both %.6g (2^-144 = %.6g)", t ? "no return in the limit" : "returned", d.r.both, 0x1p-144);
    }
    { /* D-2: the large triangle was entered in every grid cell of its box */
        dec_out o = {-9, -1};
        const int t = t_limited(d2_run, &o, sizeof o, T_LIM);
        CHECK(t == 0 && o.rc == 0 && o.nv > 0, "D-2: 1800 triangles of 1 mm beside one triangle 128 cells a side: %s in %d s of CPU, rc %d, %d vertices left", t ? "no return" : "returned", T_LIM, o.rc, o.nv);
    }
    { /* O-4: the pad vanished far along the axis, and the face at the lines' origin was lost */
        const double lo[3] = {0, 0, 0x1p60}, hi[3] = {1, 1, 0x1p60 + 0x1p24};
        const ovl_case c = ovl(lo, hi, NULL, NULL, 0.125, 0);
        CHECK(c.e == 0 && c.r.both == 0x1p24 && c.r.odd == 0, "O-4: two boxes [0,1]^2 x [2^60, 2^60 + 2^24], h = 1/8: both %.17g (2^24), %ld odd lines", c.r.both, c.r.odd);
    }
    { /* O-5, O-5b: h * h overflowed: 0 * inf = NaN where nothing was found */
        const double lo[3] = {0, 0, 0}, hi[3] = {1, 1, 1};
        const ovl_case c = ovl(lo, hi, NULL, NULL, 1e200, MK_OVL_ALL);
        CHECK(c.e == 0 && c.r.both == 0 && c.r.vol_a == 0, "O-5: two unit boxes, h = 1e200 (no line meets them): both %g, vol_a %g (0)", c.r.both, c.r.vol_a);
        const double p = 0.5371 * 1e160, q = 0.5213 * 1e160, s = 0x1p-52 * 1e160 * 8;
        const double lo2[3] = {p - s, q - s, 0}, hi2[3] = {p + s, q + s, 1e-100};
        const ovl_case d = ovl(lo2, hi2, NULL, NULL, 1e160, MK_OVL_ALL);
        CHECK(d.e == 0 && near(d.r.both, 1e220, 1e-12), "O-5b: a box around one line, h = 1e160: both %.17g (1e-100 h^2 = 1e220)", d.r.both);
    }
    { /* O-6: h * h subnormal, its bits lost */
        const double h = 0x1.123456789abcdp-525, W = 6 * h, E = 0x1p950, def = (double)(36.0L * h * h * 0x1p951L);
        const double lo[3] = {0, 0, -E}, hi[3] = {W, W, E};
        const ovl_case c = ovl(lo, hi, NULL, NULL, h, 0);
        CHECK(c.e == 0 && near(c.r.both, def, 1e-12), "O-6: 36 lines, h = 0x1.123456789abcdp-525: both %.17g against %.17g", c.r.both, def);
    }
    { /* O-7, O-7b: the line origin at the low end of both meshes lost the far mesh's crossings */
        const double lo0[3] = {-2, -2, -0x1p60}, hi0[3] = {2, 2, 1}, lo1[3] = {-0.5, -0.5, -0.5}, hi1[3] = {0.5, 0.5, 0.5}, lo2[3] = {-2, -2, -0x1p33};
        const ovl_case c = ovl(lo0, hi0, lo1, hi1, 0.125, 0), d = ovl(lo2, hi0, lo1, hi1, 0.125, 0);
        CHECK(c.e == 0 && near(c.r.both, 1, 1e-12) && c.r.odd == 0, "O-7: a unit box inside one reaching z = -2^60: both %.17g (1), %ld odd lines", c.r.both, c.r.odd);
        CHECK(d.e == 0 && near(d.r.both, 1, 1e-12), "O-7b: the same reaching -2^33: both %.17g (1)", d.r.both);
    }
    { /* H-end: the hit search left out a crossing at either end of its range */
        double V[24], t[8];
        int T[36];
        const double lo[3] = {0, 0, 0}, hi[3] = {1, 1, 1}, o[3] = {0.3, 0.4, 0};
        box(V, T, lo, hi);
        mk_bvh B;
        const int e = mk_bvh_build(&B, V, T, 12), n = e ? -1 : mk_ovl_hits_(&B, o, 2, 1, t, 8);
        if (!e) mk_bvh_free(&B);
        CHECK(n == 2, "H-end: a unit box's crossings at t = 0 and t = 1 of the line (0.3, 0.4) from z = 0, length 1: %d found (2)", n);
    }
    { /* D-5: a positive bound under about 1.5e-162 squared to 0, which means no bound */
        mesh m = {0};
        const double p = 1.618033988749895, l = sqrt(1 + p * p);
        const double P[12][3] = {{-1, p, 0}, {1, p, 0}, {-1, -p, 0}, {1, -p, 0}, {0, -1, p}, {0, 1, p}, {0, -1, -p}, {0, 1, -p}, {p, 0, -1}, {p, 0, 1}, {-p, 0, -1}, {-p, 0, 1}};
        const int F[20][3] = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                              {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
        for (int i = 0; i < 12; i++) dv(&m, P[i][0] / l, P[i][1] / l, P[i][2] / l);
        for (int i = 0; i < 20; i++) dt(&m, F[i][0], F[i][1], F[i][2]);
        mk_dec_opts o;
        mk_dec_report r;
        memset(&o, 0, sizeof o), o.dist = 1e-170, o.turn = -1, o.nregion = 1;
        const int e = dec_call(&m, &o, &r);
        CHECK(e == 0 && r.nv == 12, "D-5: a unit icosahedron under the bound 1e-170: %d of its 12 vertices left (12)", e ? -1 : r.nv);
        free(m.V), free(m.T);
    }
    { /* D-6: maxown counted a collapse the crossing test then refused */
        mesh m = {0};
        const int c0 = dv(&m, 0, 0, 0), c1 = dv(&m, 1, 0, 0), c2 = dv(&m, 1, 1, 0), c3 = dv(&m, 0, 1, 0), u = dv(&m, 0.5, 0.5, 0.1);
        dt(&m, u, c0, c1), dt(&m, u, c1, c2), dt(&m, u, c2, c3), dt(&m, u, c3, c0);
        const int a = dv(&m, 0.45, 0.48, 0.05), b = dv(&m, 0.55, 0.5, -0.3), c = dv(&m, 0.5, 0.56, 0.04);
        dt(&m, a, b, c);
        mk_dec_opts o;
        mk_dec_report r;
        memset(&o, 0, sizeof o), o.dist = 0.2, o.turn = -1, o.nregion = 1;
        const int e = dec_call(&m, &o, &r);
        CHECK(e == 0 && r.nv == 8 && r.maxown == 0, "D-6: a tent over a small triangle, nothing removable: %d of 8 vertices left, maxown %g (0)", e ? -1 : r.nv, e ? -1.0 : r.maxown);
        free(m.V), free(m.T);
    }
    { /* G-1: a finite lo whose components' sum overflows was refused */
        const int n[3] = {2, 2, 2};
        const double lo[3] = {DBL_MAX, DBL_MAX, 0};
        const float F[8] = {-1, 1, 1, 1, 1, 1, 1, 1};
        mk_gmesh o;
        const int r = mk_mtet(F, n, lo, 1.0, &o);
        CHECK(r == 0 && o.nv > 0, "G-1: mk_mtet with lo (DBL_MAX, DBL_MAX, 0), finite: returns %d, %d vertices", r, o.nv);
        mk_gmesh_free(&o);
    }
    { /* G-2: a far node gave infinite coordinates */
        const int n[3] = {4, 2, 2};
        const double lo[3] = {0, 0, 0};
        const float F[16] = {1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        mk_gmesh o;
        const int r = mk_mtet(F, n, lo, DBL_MAX / 2, &o);
        long nf = 0;
        for (long i = 0; !r && i < 3 * (long)o.nv; i++) nf += !isfinite(o.V[i]);
        CHECK(r == -1 || (r == 0 && nf == 0), "G-2: mk_mtet with h = DBL_MAX / 2, a far corner past DBL_MAX: returns %d, %ld coordinates not finite", r, nf);
        mk_gmesh_free(&o);
    }
    { /* G-3: a seam psi of +inf on part of the kept side gave the patch NaN cut points */
        static const double V[8 * 3] = {-1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1};
        static const int T[12 * 3] = {0, 3, 2, 0, 2, 1, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6, 1, 2, 6, 1, 6, 5, 0, 4, 7, 0, 7, 3};
        double lim = 0.05;
        mk_graft g;
        const int ei = mk_graft_init(&g, V, 8, T, 12, psi_inf, &lim, 0.0, 1e-9);
        const int n[3] = {13, 13, 13};
        const double lo[3] = {-1.5, -1.5, -1.5}, h = 0.25;
        float *F = malloc(sizeof *F * 13 * 13 * 13);
        for (int k = 0; k < 13; k++)
            for (int j = 0; j < 13; j++)
                for (int i = 0; i < 13; i++) {
                    const double x = lo[0] + i * h, y = lo[1] + j * h, z = lo[2] + k * h;
                    F[(k * 13 + j) * 13 + i] = (float)(sqrt(x * x + y * y + z * z) - 1.1);
                }
        mk_gmesh P = {0};
        const int r = ei ? -9 : mk_graft_patch(&g, F, n, lo, h, &P);
        long nf = 0;
        for (long i = 0; !r && i < 3 * (long)P.nv; i++) nf += !isfinite(P.V[i]);
        CHECK(r == 0 && P.nv > 0 && nf == 0, "G-3: a cube grafted across a seam function of +inf past z = 0.05: patch returns %d, %d vertices, %ld coordinates not finite", r, P.nv, nf);
        if (!r) mk_gmesh_free(&P);
        if (!ei) mk_graft_free(&g);
        free(F);
    }
    { /* G-5: mk_mtet's map freed its keys and values when one malloc failed, and mk_mtet freed them again */
        const int n[3] = {2, 2, 2};
        const double lo[3] = {0, 0, 0};
        const float F[8] = {-1, 1, 1, 1, 1, 1, 1, 1};
        mk_gmesh o;
        t_al_on = 1, t_al_n = 0, t_al_fail = 1, t_al_twice = 0, t_al_ng = 0;
        const int r = mk_mtet(F, n, lo, 1.0, &o);
        t_al_on = 0, t_al_fail = -1;
        CHECK(r == -1 && t_al_twice == 0, "G-5: mk_mtet with its second allocation failing: returns %d, %ld blocks freed twice", r, t_al_twice);
    }
    { /* R-1: mk_graft_refine took coordinates that are not finite */
        static const double bad[2] = {NAN, INFINITY};
        int nref = 0;
        for (int b = 0; b < 2; b++) {
            double V[12] = {0, 0, -1, 1, 0, -1, 1, 0, 1, 0, 0, 1};
            const int T[6] = {0, 1, 2, 0, 2, 3};
            V[6] = bad[b];
            mk_gmesh o;
            const int r = mk_graft_refine(V, 4, T, 2, plane_z, NULL, 0.3, 0.5, &o);
            nref += r == -1;
            if (r >= 0) mk_gmesh_free(&o);
        }
        CHECK(nref == 2, "R-1: mk_graft_refine with a NaN, then an infinite coordinate: %d of 2 refused", nref);
    }
    { /* G1: mk_aabb_ray's pad ignored the ray origin's size, so a ray from about 1.7e12 away lost a hit on the unit cube */
        const double V[24] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1};
        const int T[36] = {0, 2, 1, 1, 2, 3, 4, 5, 6, 5, 7, 6, 0, 1, 4, 1, 5, 4, 2, 6, 3, 3, 6, 7, 0, 4, 2, 2, 4, 6, 1, 3, 5, 3, 7, 5};
        const double o[3] = {0x1.7f6f8a0efdp+40, 0x1.1736b51d31p+40, 0x1.004bd17526p+40};
        const double d[3] = {-0x1.7f6f8a0efcfdep+40, -0x1.1736b51d30fffp+40, -0x1.004bd17526p+40};
        mk_bvh B;
        int edge = -1, wind = -9;
        const int e = mk_bvh_build(&B, V, T, 12), n = e ? -1 : mk_bvh_count_wind_(&B, o, d, &edge, &wind);
        if (!e) mk_bvh_free(&B);
        CHECK(n == 2 && edge == 0 && wind == 0, "G1: a ray from 1.7e12 away through the closed unit cube: %d hits (2), facings summed %d (0, the point is outside)", n, wind);
    }
    { /* S-1: a box refused on axis 1 or 2 left axis 0's fields in *g */
        double V[24];
        int T[36], ne = 0;
        sdf_cube(V, T, 1);
        static const mk_sdf z;
        for (int c = 0; c < 3; c++) {
            double b[6] = {-0.5, -0.5, -0.5, 0.5, 0.5, 0.5};
            if (c < 2) b[1 + c] = NAN;
            else b[4] = -0.6; /* hi below lo on axis 1 */
            mk_sdf g;
            memset(&g, 0x5a, sizeof g);
            const int r = mk_sdf_mesh_box(V, 8, T, 12, 0.1, b, c == 2 ? MK_SDF_LATTICE : 0, &g);
            ne += r == -1 && !memcmp(&g, &z, sizeof g);
        }
        CHECK(ne == 3, "S-1: a box refused on axis 1 or 2 (NaN, NaN, hi below lo): -1 with *g empty in %d of 3", ne);
    }
    { /* S-2: a triangle more than INT_MAX cells from node 0: (int) of its cell was undefined, its crossings lost, every node read outside */
        double V[24];
        int T[36];
        const double s = 0x1p33;
        sdf_cube(V, T, s);
        const double b[6] = {-2, -2, -2, 2, 2, 2};
        mk_sdf g;
        const int r = mk_sdf_mesh_box(V, 8, T, 12, 1, b, 0, &g);
        long out = 0, N = 0;
        if (!r) {
            N = (long)g.n[0] * g.n[1] * g.n[2];
            for (long i = 0; i < N; i++) out += g.d[i] > 0;
            mk_sdf_free(&g);
        }
        CHECK(r == 0 && N == 125 && out == 0, "S-2: a closed cube 2^34 cells across around a grid of %ld nodes: %ld read outside (0)", N, out);
    }
    return wrong;
}
