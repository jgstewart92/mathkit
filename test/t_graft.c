/* mk_mtet and the graft: marching tetrahedra on a sphere's and a torus's distance grids (closed, outward, the volume
   against 4/3 pi and 2 pi^2 R r^2), then a dent carved into a unit icosphere S inside a 45-degree cap: the field from
   mk_sdf_mesh, the patch thinned with mk_decimate between patch and stitch, stitched back; the result closed with S's
   own triangles kept outside the cap, its volume against S's less the analytic lens of the carved ball, nothing inside
   that ball. Then the patch stitched unthinned (on that S and a coarse one), mk_graft_refine, a hole in the patch by
   the seam refused, and the kept piece of most volume either sign (S inward with a speck outward; a large inward piece
   beside an outward body: 2). With a folder argument the results are written there for test/check.sh's meshcheck. */
#include "../include/mk_dec.h"
#include "../include/mk_graft.h"
#include "../include/mk_sdf.h"
#include "t.h"
#include "tmesh.h"

/* closed and consistently wound: every directed edge once, its reverse once; and the Euler characteristic */
static int closed(const double *V, int nv, const int *T, int nt, int *chi)
{
    (void)V;
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
/* the cap: 45 degrees about +z, psi < 0 inside it; the dent: a ball of radius 0.3 about the pole */
static const double CAP = 0.7853981633974483, BALL = 0.3;
static double psi(const double *p, void *ctx)
{
    (void)ctx;
    const double l = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    return l > 0 ? acos(fmax(-1, fmin(1, p[2] / l))) - CAP : -CAP;
}
static double lensvol(double R, double r, double d) { return T_PI * (R + r - d) * (R + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * R + 6 * r * R - 3 * R * R) / (12 * d); }

int main(int argc, char **argv)
{
    tm_dir = argc > 1 ? argv[1] : NULL;
    const double c0[3] = {0, 0, 0};
    { /* mk_mtet on analytic fields */
        const int n[3] = {61, 61, 61};
        const double lo[3] = {-1.5, -1.5, -1.5}, h = 0.05;
        float *F = malloc(sizeof(float) * 61 * 61 * 61), *G = malloc(sizeof(float) * 61 * 61 * 61);
        for (int k = 0; k < 61; k++)
            for (int j = 0; j < 61; j++)
                for (int i = 0; i < 61; i++) {
                    const double x = lo[0] + i * h, y = lo[1] + j * h, z = lo[2] + k * h;
                    F[(k * 61 + j) * 61 + i] = (float)(sqrt(x * x + y * y + z * z) - 1);
                    G[(k * 61 + j) * 61 + i] = (float)(hypot(hypot(x, y) - 0.9, z) - 0.35);
                }
        mk_gmesh s, t;
        int cs, ct;
        const int r1 = mk_mtet(F, n, lo, h, &s), r2 = mk_mtet(G, n, lo, h, &t);
        const double vs = tm_volume(s.V, s.T, s.nt), vt = tm_volume(t.V, t.T, t.nt), es = 4.0 / 3 * T_PI, et = 2 * T_PI * T_PI * 0.9 * 0.35 * 0.35;
        const int ks = r1 == 0 && closed(s.V, s.nv, s.T, s.nt, &cs), kt = r2 == 0 && closed(t.V, t.nv, t.T, t.nt, &ct);
        CHECK(ks && cs == 2 && fabs(vs / es - 1) < 0.01, "marching tetrahedra on a unit sphere's field, h 0.05: closed, Euler %d, volume %.4f against %.4f (%.2f %%)", cs, vs, es, 100 * (vs / es - 1));
        CHECK(kt && ct == 0 && fabs(vt / et - 1) < 0.02, "on a torus's: closed, Euler %d, volume %.4f against %.4f (%.2f %%)", ct, vt, et, 100 * (vt / et - 1));
        tm_write("mtet_sphere", s.V, s.nv, s.T, s.nt, NULL), tm_write("mtet_torus", t.V, t.nv, t.T, t.nt, NULL);
        int k = 0, nr = 0;
        const int n1[3] = {1, 61, 61}, nb[3] = {1 << 16, 1 << 16, 1 << 16};
        mk_gmesh z;
#define REF(call) (k++, nr += (call) == -1 && z.V == NULL)
        REF(mk_mtet(F, n1, lo, h, &z));
        REF(mk_mtet(F, n, lo, 0, &z));
        REF(mk_mtet(F, n, lo, NAN, &z));
        REF(mk_mtet(NULL, n, lo, h, &z));
        REF(mk_mtet(F, nb, lo, h, &z));
#undef REF
        CHECK(nr == k, "mk_mtet refuses %d of %d: a size under 2, h 0 or NaN, no field, nodes past INT_MAX", nr, k);
        mk_gmesh_free(&s), mk_gmesh_free(&t), free(F), free(G);
    }
    /* the graft: S a unit icosphere; the field S's own distance (mk_sdf_mesh) carved by a ball at the pole */
    tm_mesh S = {0};
    tm_sphere(&S, 1, c0, 5);
    mk_sdf g;
    const double h = 0.02;
    int r = mk_sdf_mesh(S.V, S.nv, S.T, S.nt, h, 0.06, 0, &g);
    const size_t N = (size_t)g.n[0] * g.n[1] * g.n[2];
    for (size_t q = 0; q < N && !r; q++) {
        const int i = (int)(q % (size_t)g.n[0]), j = (int)(q / (size_t)g.n[0] % (size_t)g.n[1]), k = (int)(q / ((size_t)g.n[0] * g.n[1]));
        const double x = g.lo[0] + i * h, y = g.lo[1] + j * h, z = g.lo[2] + k * h - 1;
        g.d[q] = (float)fmax(g.d[q], BALL - sqrt(x * x + y * y + z * z));
    }
    mk_graft G;
    const int ri = r ? -1 : mk_graft_init(&G, S.V, S.nv, S.T, S.nt, psi, NULL, 2 * h, 1e-6);
    mk_gmesh P;
    const int rp = ri ? -1 : mk_graft_patch(&G, g.d, g.n, g.lo, g.h, &P);
    CHECK(r == 0 && ri == 0 && rp == 0 && P.nt > 1000, "the patch: mk_sdf_mesh's grid at h 0.02 carved, its surface by marching tetrahedra kept in the cap: %d vertices, %d triangles", rp ? -1 : P.nv,
          rp ? -1 : P.nt);
    /* between patch and stitch: the patch thinned as an open surface, its rim moving only along itself */
    mk_dec_opts o;
    mk_dec_defaults(&o, 1);
    o.dist = 0.2 * h, o.maxedge = 0.1, o.minalt = 1e-4 * h;
    unsigned char *alive = malloc((size_t)P.nv);
    int *To = malloc(sizeof(int) * 3 * (size_t)P.nt), nto = 0;
    mk_dec_report dr;
    const int rd = mk_decimate(P.V, P.nv, P.T, P.nt, &o, alive, To, &nto, &dr);
    CHECK(rd == 0 && dr.nv < P.nv / 2 && dr.rim1 > 0, "the patch thinned by mk_decimate: %d -> %d vertices, rim %ld -> %ld", P.nv, dr.nv, dr.rim0, dr.rim1);
    int *map = malloc(sizeof(int) * (size_t)P.nv), kv = 0;
    double *PV = malloc(sizeof(double) * 3 * (size_t)P.nv);
    for (int i = 0; i < P.nv; i++)
        if (alive[i]) map[i] = kv, memcpy(PV + 3 * kv, P.V + 3 * i, 3 * sizeof *PV), kv++;
    for (int t = 0; t < 3 * nto; t++) To[t] = map[To[t]];
    mk_gmesh O;
    mk_graft_report gr;
    const int rs = mk_graft_stitch(&G, PV, kv, To, nto, 2.5 * h, &O, &gr);
    int chi = 0;
    const int ok = rs == 0 && closed(O.V, O.nv, O.T, O.nt, &chi);
    CHECK(ok && chi == 2 && gr.dup == 0 && gr.lone == 0 && gr.pieces == 1, "stitched: closed and consistently wound, Euler %d, %d piece; %ld rim vertices put into S's rim, %ld welded, largest move %.5f (h %.2f)", chi,
          gr.pieces, gr.inserted, gr.welded, gr.snap, h);
    const double vexp = tm_volume(S.V, S.T, S.nt) - lensvol(1, BALL, 1);
    CHECK(ok && fabs(gr.vol1 / vexp - 1) < 3e-3, "the volume %.5f against S's less the ball's lens %.5f (%.3f %%)", gr.vol1, vexp, 100 * (gr.vol1 / vexp - 1));
    long kept = 0, own = 0, inball = 0, bad1 = 0;
    for (int i = 0; ok && i < O.nv; i++) {
        const double *p = O.V + 3 * i, z = p[2] - 1;
        inball += sqrt(p[0] * p[0] + p[1] * p[1] + z * z) < BALL - 0.5 * h;
        if (O.keep[i] == 1) { /* S's own vertex: exactly one of S's */
            int f = 0;
            for (int j = 0; j < S.nv && !f; j++) f = !memcmp(p, S.V + 3 * j, 3 * sizeof *p);
            bad1 += !f, own++;
        }
    }
    for (int t = 0; ok && t < S.nt; t++) { /* each S triangle beyond the seam's reach is in the result as it was */
        int out = 1;
        for (int k = 0; k < 3; k++) out &= psi(S.V + 3 * S.T[3 * t + k], NULL) > 0.05;
        if (!out) continue;
        int f = 0;
        for (int u = 0; u < O.nt && !f; u++) {
            int same = 1;
            for (int k = 0; k < 3; k++) same &= !memcmp(O.V + 3 * O.T[3 * u + k], S.V + 3 * S.T[3 * t + k], 3 * sizeof(double));
            f = same;
        }
        kept += f, kept -= !f * 1000000;
    }
    CHECK(ok && bad1 == 0 && own > 0 && kept > 10000, "S's own vertices (%ld) are S's, digit for digit; %ld of its triangles away from the seam kept as they were", own, kept);
    CHECK(ok && inball == 0, "nothing of the result inside the carved ball (%ld vertices)", inball);
    tm_write("graft", O.V, O.nv, O.T, O.nt, NULL);
    { /* the patch stitched as it comes, not thinned: its dense rim wiggles back on S's edges, so some rim vertices come
         out of order and must be merged; closed, one piece, the expected volume, on this S and on a coarse one
         (3 subdivisions, edges about 6 h) */
        mk_gmesh O2 = {0};
        mk_graft_report r2 = {0}, r3 = {0};
        const int s2 = mk_graft_stitch(&G, P.V, P.nv, P.T, P.nt, 2.5 * h, &O2, &r2);
        int c2 = 0, c3 = 0;
        const int k2 = s2 == 0 && closed(O2.V, O2.nv, O2.T, O2.nt, &c2);
        tm_mesh S3 = {0};
        tm_sphere(&S3, 1, c0, 3);
        mk_sdf g3;
        int s3 = mk_sdf_mesh(S3.V, S3.nv, S3.T, S3.nt, h, 0.06, 0, &g3), k3 = 0;
        for (size_t q = 0; !s3 && q < (size_t)g3.n[0] * g3.n[1] * g3.n[2]; q++) {
            const int i = (int)(q % (size_t)g3.n[0]), j = (int)(q / (size_t)g3.n[0] % (size_t)g3.n[1]), k = (int)(q / ((size_t)g3.n[0] * g3.n[1]));
            const double x = g3.lo[0] + i * h, y = g3.lo[1] + j * h, z = g3.lo[2] + k * h - 1;
            g3.d[q] = (float)fmax(g3.d[q], BALL - sqrt(x * x + y * y + z * z));
        }
        mk_graft G3;
        mk_gmesh P3 = {0}, O3 = {0};
        const int i3 = s3 ? -1 : mk_graft_init(&G3, S3.V, S3.nv, S3.T, S3.nt, psi, NULL, 2 * h, 1e-6);
        if (!s3) s3 = i3 ? -9 : mk_graft_patch(&G3, g3.d, g3.n, g3.lo, g3.h, &P3) ? -8 : mk_graft_stitch(&G3, P3.V, P3.nv, P3.T, P3.nt, 2.5 * h, &O3, &r3);
        k3 = s3 == 0 && closed(O3.V, O3.nv, O3.T, O3.nt, &c3);
        const double v3 = tm_volume(S3.V, S3.T, S3.nt) - lensvol(1, BALL, 1);
        long x2 = -1, x3 = -1;
        mk_bvh B2;
        if (k2 && !mk_bvh_build(&B2, O2.V, O2.T, O2.nt)) x2 = mk_bvh_crossings(&B2, NULL), mk_bvh_free(&B2);
        if (k3 && !mk_bvh_build(&B2, O3.V, O3.T, O3.nt)) x3 = mk_bvh_crossings(&B2, NULL), mk_bvh_free(&B2);
        CHECK(k2 && c2 == 2 && r2.dup == 0 && r2.lone == 0 && r2.pieces == 1 && r2.merged > 0 && fabs(r2.vol1 / vexp - 1) < 3e-3 && k3 && c3 == 2 && r3.dup == 0 &&
                  r3.lone == 0 && r3.pieces == 1 && r3.merged > 0 && fabs(r3.vol1 / v3 - 1) < 1e-2,
              "the patch not thinned (%d and %d triangles): closed, Euler %d and %d, rc %d and %d, dup %ld and %ld, lone %ld and %ld; %ld and %ld rim vertices out of order merged; volume %.3f %% and %.3f %% off; crossing pairs by the seam %ld and %ld (not judged: thinning lays them flat)",
              P.nt, P3.nt, c2, c3, s2, s3, r2.dup, r3.dup, r2.lone, r3.lone, r2.merged, r3.merged, 100 * (r2.vol1 / vexp - 1), 100 * (r3.vol1 / v3 - 1), x2, x3);
        { /* mk_graft_refine: the coarse S refined where the seam cuts it (edges to 2 h within 2 h of it), then grafted
             as before */
            mk_gmesh R4 = {0}, P4 = {0}, O4 = {0};
            mk_graft G4;
            mk_graft_report r4 = {0};
            const int rr = s3 ? -9 : mk_graft_refine(S3.V, S3.nv, S3.T, S3.nt, psi, NULL, 2 * h, 2 * h, &R4);
            int cr = 0, c4 = 0;
            const int kr = rr > 0 && closed(R4.V, R4.nv, R4.T, R4.nt, &cr);
            long left = 0, kept = 0, away = 0;
            for (int t = 0; kr && t < R4.nt; t++)
                for (int k = 0; k < 3; k++) {
                    const double *a = R4.V + 3 * R4.T[3 * t + k], *b = R4.V + 3 * R4.T[3 * t + (k + 1) % 3], pa = psi(a, NULL), pb = psi(b, NULL);
                    left += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])) > 2 * h && fmin(pa, pb) < 2 * h && fmax(pa, pb) > -2 * h;
                }
            for (int t = 0; kr && t < S3.nt; t++) { /* S's triangles beyond the band: kept, the same corners in order */
                double lo = INFINITY, hi = -INFINITY;
                for (int k = 0; k < 3; k++) lo = fmin(lo, psi(S3.V + 3 * S3.T[3 * t + k], NULL)), hi = fmax(hi, psi(S3.V + 3 * S3.T[3 * t + k], NULL));
                if (lo < 2 * h && hi > -2 * h) continue;
                away++;
                for (int u = 0; u < R4.nt; u++)
                    if (!memcmp(R4.T + 3 * u, S3.T + 3 * t, 3 * sizeof(int))) { kept++; break; }
            }
            const double vr = kr ? tm_volume(R4.V, R4.T, R4.nt) : 0, vs3 = tm_volume(S3.V, S3.T, S3.nt);
            const int i4 = kr ? mk_graft_init(&G4, R4.V, R4.nv, R4.T, R4.nt, psi, NULL, 2 * h, 1e-6) : -1;
            const int s4 = i4 ? -9 : mk_graft_patch(&G4, g3.d, g3.n, g3.lo, g3.h, &P4) ? -8 : mk_graft_stitch(&G4, P4.V, P4.nv, P4.T, P4.nt, 2.5 * h, &O4, &r4);
            const int k4 = s4 == 0 && closed(O4.V, O4.nv, O4.T, O4.nt, &c4);
            long x4 = -1;
            mk_bvh B4;
            if (k4 && !mk_bvh_build(&B4, O4.V, O4.T, O4.nt)) x4 = mk_bvh_crossings(&B4, NULL), mk_bvh_free(&B4);
            CHECK(kr && cr == 2 && left == 0 && kept == away && away > 500 && fabs(vr / vs3 - 1) < 1e-12 && k4 && c4 == 2 && fabs(r4.vol1 / v3 - 1) < 1e-2,
                  "mk_graft_refine on the coarse S: %d -> %d triangles in %d passes, closed, Euler %d, volume %.1e off, %ld long edges left by the seam, %ld of %ld triangles away from it as they were; grafted: rc %d, Euler %d, rim %d -> %d edges, %ld crossing pairs (unrefined %ld)",
                  S3.nt, R4.nt, rr, cr, vr / vs3 - 1, left, kept, away, s4, c4, r3.rim_s, r4.rim_s, x4, x3);
            if (!i4) mk_graft_free(&G4);
            mk_gmesh_free(&R4), mk_gmesh_free(&P4), mk_gmesh_free(&O4);
        }
        /* a hole in the patch by the seam: one triangle off its rim, its corners within 1.2 h of the seam, taken out */
        int hole = -1;
        for (int t = 0; t < P.nt && hole < 0; t++) {
            int near = 1;
            for (int k = 0; k < 3; k++) near &= !P.keep[P.T[3 * t + k]] && fabs(psi(P.V + 3 * P.T[3 * t + k], NULL)) < 1.2 * h;
            if (near) hole = t;
        }
        int *Th = malloc(sizeof(int) * 3 * (size_t)P.nt), nh = 0;
        for (int t = 0; t < P.nt; t++)
            if (t != hole) memcpy(Th + 3 * nh++, P.T + 3 * t, 3 * sizeof(int));
        mk_gmesh Oh = {0};
        mk_graft_report rh = {0};
        const int sh = hole < 0 ? -9 : mk_graft_stitch(&G, P.V, P.nv, Th, nh, 2.5 * h, &Oh, &rh);
        CHECK(hole >= 0 && sh != 0 && rh.loops_s == 1 && rh.loops_p == 2 && rh.unpaired == 1,
              "a patch with a one-triangle hole by the seam: not passed (rc %d), rim loops S %d, patch %d, one without a partner (%d)", sh, rh.loops_s, rh.loops_p, rh.unpaired);
        /* the piece kept: S wound inward (a cavity's shell) with its patch, and a small closed speck wound outward
           beside them; and a large closed piece wound inward beside an outward body */
        int *Tr = malloc(sizeof(int) * 3 * (size_t)S.nt);
        for (int t = 0; t < S.nt; t++) Tr[3 * t] = S.T[3 * t], Tr[3 * t + 1] = S.T[3 * t + 2], Tr[3 * t + 2] = S.T[3 * t + 1];
        float *Fn = malloc(sizeof(float) * N);
        for (size_t q = 0; q < N; q++) Fn[q] = -g.d[q];
        mk_graft Gi;
        mk_gmesh Pi = {0}, Oi = {0}, Ob = {0};
        mk_graft_report ri = {0};
        tm_mesh sp = {0};
        tm_sphere(&sp, 0.05, c0, 1);
        const int ii = mk_graft_init(&Gi, S.V, S.nv, Tr, S.nt, psi, NULL, 2 * h, 1e-6);
        int si = ii ? -9 : mk_graft_patch(&Gi, Fn, g.n, g.lo, g.h, &Pi) ? -8 : 0, ci = 0;
        if (!si) {
            tm_mesh A = {0};
            for (int i = 0; i < Pi.nv; i++) tm_vert(&A, Pi.V[3 * i], Pi.V[3 * i + 1], Pi.V[3 * i + 2]);
            for (int t = 0; t < Pi.nt; t++) tm_tri(&A, Pi.T[3 * t], Pi.T[3 * t + 1], Pi.T[3 * t + 2]);
            for (int i = 0; i < sp.nv; i++) tm_vert(&A, sp.V[3 * i], sp.V[3 * i + 1], sp.V[3 * i + 2]);
            for (int t = 0; t < sp.nt; t++) tm_tri(&A, Pi.nv + sp.T[3 * t], Pi.nv + sp.T[3 * t + 1], Pi.nv + sp.T[3 * t + 2]);
            si = mk_graft_stitch(&Gi, A.V, A.nv, A.T, A.nt, 2.5 * h, &Oi, &ri);
            tm_free(&A);
        }
        const int ki = si == 0 && closed(Oi.V, Oi.nv, Oi.T, Oi.nt, &ci);
        CHECK(ki && ci == 2 && ri.pieces == 2 && ri.inverted == 1 && fabs(ri.vol1 / -vexp - 1) < 3e-3,
              "S wound inward with its patch, a speck wound outward beside: the body kept (volume %.5f, Euler %d, rc %d), not the speck; %d pieces, %d wound against S", ri.vol1, ci,
              si, ri.pieces, ri.inverted);
        /* beside the outward body, a closed piece wound inward: radius 1.5 (kept: of most volume), 0.5 (dropped, 12 %
           of the body's volume), 0.05 (dropped, 0.01 %: the specks of rounding a clean stitch leaves are of that
           order). The stitch returns 2 when the kept piece is wound against S or such pieces hold over 1 % of its
           volume, so 2, 2 and 0 */
        const double rad[3] = {1.5, 0.5, 0.05};
        int sb[3];
        mk_graft_report rq[3];
        for (int c = 0; c < 3; c++) {
            tm_mesh B = {0}, bs = {0};
            tm_sphere(&bs, rad[c], c0, 2);
            for (int t = 0; t < bs.nt; t++) { const int x = bs.T[3 * t + 1]; bs.T[3 * t + 1] = bs.T[3 * t + 2], bs.T[3 * t + 2] = x; }
            for (int i = 0; i < kv; i++) tm_vert(&B, PV[3 * i], PV[3 * i + 1], PV[3 * i + 2]);
            for (int t = 0; t < nto; t++) tm_tri(&B, To[3 * t], To[3 * t + 1], To[3 * t + 2]);
            for (int i = 0; i < bs.nv; i++) tm_vert(&B, bs.V[3 * i], bs.V[3 * i + 1], bs.V[3 * i + 2]);
            for (int t = 0; t < bs.nt; t++) tm_tri(&B, kv + bs.T[3 * t], kv + bs.T[3 * t + 1], kv + bs.T[3 * t + 2]);
            sb[c] = mk_graft_stitch(&G, B.V, B.nv, B.T, B.nt, 2.5 * h, &Ob, &rq[c]);
            mk_gmesh_free(&Ob), tm_free(&B), tm_free(&bs);
        }
        CHECK(sb[0] == 2 && rq[0].vol1 < 0 && fabs(rq[0].dropped_vol / gr.vol1 - 1) < 1e-9 && sb[1] == 2 && fabs(rq[1].vol1 / gr.vol1 - 1) < 1e-9 && rq[1].inverted_vol < -0.01 * rq[1].vol1 &&
                  sb[2] == 0 && fabs(rq[2].vol1 / gr.vol1 - 1) < 1e-9 && rq[2].inverted == 1 && rq[2].inverted_vol < 0 && rq[2].inverted_vol > -1e-3 * rq[2].vol1,
              "a closed piece wound inward beside the body, radius 1.5, 0.5, 0.05: rc %d, %d, %d (2: kept %.3f, or dropped %.4f against S over 1 %% of the body; the speck, %.2g, reported only); %d, %d, %d inverted",
              sb[0], sb[1], sb[2], rq[0].vol1, rq[1].inverted_vol, rq[2].inverted_vol, rq[0].inverted, rq[1].inverted, rq[2].inverted);
        tm_free(&sp), free(Tr), free(Fn), free(Th), mk_gmesh_free(&Pi), mk_gmesh_free(&Oi), mk_gmesh_free(&Ob), mk_gmesh_free(&Oh);
        if (!ii) mk_graft_free(&Gi);
        if (!i3) mk_graft_free(&G3);
        mk_gmesh_free(&O2), mk_gmesh_free(&O3), mk_gmesh_free(&P3), mk_sdf_free(&g3), tm_free(&S3);
    }
    {
        mk_gmesh z;
        mk_graft_report zr;
        int k = 0, nr = 0;
        const int Tb[3] = {0, 1, 99999999};
#define REF(call) (k++, nr += (call) == -1)
        REF(mk_graft_stitch(&G, PV, kv, To, nto, -1, &z, &zr));
        REF(mk_graft_stitch(&G, PV, kv, To, nto, NAN, &z, &zr));
        REF(mk_graft_stitch(&G, PV, kv, Tb, 1, 2.5 * h, &z, &zr));
        REF(mk_graft_stitch(&G, PV, kv, To, nto, 1e-9, &z, &zr)); /* the rims too far apart for that snap */
        mk_graft G2;
        REF(mk_graft_init(&G2, S.V, S.nv, S.T, S.nt, NULL, NULL, h, 0));
        REF(mk_graft_init(&G2, S.V, S.nv, S.T, S.nt, psi, NULL, -1, 0));
        REF(mk_graft_init(&G2, S.V, 2, S.T, S.nt, psi, NULL, h, 0));
#undef REF
        CHECK(nr == k, "graft refusals: %d of %d (snap negative, NaN or too small for the rims, an index out of range; no psi, near negative, an index past nv)", nr, k);
    }
    mk_gmesh_free(&O), mk_gmesh_free(&P), free(alive), free(To), free(map), free(PV), mk_graft_free(&G), mk_sdf_free(&g), tm_free(&S);
    return wrong;
}
