/* mk_mtet and the graft. mk_mtet on a sphere's and a torus's distance grids: closed, no vertex pinched, no zero-area
   triangle, the Euler characteristic, and the volume against its law (low by (h/R)^2 / 2 for a ball, (h/r)^2 / 6 for a
   tube of radius r); a sphere cut by the box, open only on the box's faces; a field 0 at a plane of nodes, no zero-area
   triangle. Then a dent carved into a unit icosphere S inside a 45-degree cap: the field from mk_sdf_mesh, the patch
   thinned with mk_decimate between patch and stitch, stitched back; the result closed, rc 0 with no fault, S's own
   triangles kept outside the cap, its volume against S's less the analytic lens of the carved ball, nothing inside that
   ball. Then the patch stitched unthinned (closed on that S; on a coarse one refused, rc 5, its rim triangles crossing
   S's), mk_graft_refine, a hole in the patch by the seam refused with a reason, the kept piece of most volume either
   sign, rc 2 before rc 1, a split (rc 3), a zero-area triangle (rc 4), a seam just off one of S's vertices (no sliver),
   a pinched vertex (rc 1), and the refusals with their reasons. With a folder argument the results are written there
   for test/check.sh's meshcheck. */
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
static int ufind(int *p, int x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; }
/* vertices whose triangles do not all join round them, edge to edge (each vertex's triangles joined where two share a
   second corner) */
static long pinched(int nv, const int *T, int nt)
{
    int *off = calloc((size_t)nv + 1, sizeof *off), *L = malloc(sizeof *L * 3 * (size_t)nt), *p = malloc(sizeof *p * 3 * (size_t)nt);
    for (int i = 0; i < 3 * nt; i++) off[T[i] + 1]++;
    for (int v = 0; v < nv; v++) off[v + 1] += off[v];
    int *fill = malloc(sizeof *fill * (size_t)(nv ? nv : 1));
    memcpy(fill, off, sizeof *fill * (size_t)nv);
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) L[fill[T[3 * t + k]]++] = t;
    long np = 0;
    for (int v = 0; v < nv; v++) {
        const int a = off[v], n = off[v + 1] - a;
        for (int i = 0; i < n; i++) p[a + i] = a + i;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                int sh = 0;
                for (int x = 0; x < 3; x++)
                    for (int y = 0; y < 3; y++) sh += T[3 * L[a + i] + x] != v && T[3 * L[a + i] + x] == T[3 * L[a + j] + y];
                if (sh) p[ufind(p, a + i)] = ufind(p, a + j);
            }
        int g = 0;
        for (int i = 0; i < n; i++) g += ufind(p, a + i) == a + i;
        np += g > 1;
    }
    free(off), free(L), free(p), free(fill);
    return np;
}
/* the smallest |2A| of the triangles; *nz: how many have |2A| under z */
static double min2a(const double *V, const int *T, int nt, double z, long *nz)
{
    double m = INFINITY;
    *nz = 0;
    for (int t = 0; t < nt; t++) {
        const double *a = V + 3 * T[3 * t], *b = V + 3 * T[3 * t + 1], *c = V + 3 * T[3 * t + 2];
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double x = u[1] * w[2] - u[2] * w[1], y = u[2] * w[0] - u[0] * w[2], q = u[0] * w[1] - u[1] * w[0], A2 = sqrt(x * x + y * y + q * q);
        m = fmin(m, A2), *nz += !(A2 >= z);
    }
    return m;
}
static long xpairs(const double *V, const int *T, int nt)
{
    mk_bvh B;
    long x = -1;
    if (!mk_bvh_build(&B, V, T, nt)) x = mk_bvh_crossings(&B, NULL), mk_bvh_free(&B);
    return x;
}
/* the cap: 45 degrees about +z, psi < 0 inside it; the dent: a ball of radius 0.3 about the pole */
static const double CAP = 0.7853981633974483, BALL = 0.3;
static double psi(const double *p, void *ctx)
{
    (void)ctx;
    const double l = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    return l > 0 ? acos(fmax(-1, fmin(1, p[2] / l))) - CAP : -CAP;
}
static double psi_off(const double *p, void *ctx) { return psi(p, NULL) - *(const double *)ctx; } /* the seam moved by *ctx */
static double lensvol(double R, double r, double d) { return T_PI * (R + r - d) * (R + r - d) * (d * d + 2 * d * r - 3 * r * r + 2 * d * R + 6 * r * R - 3 * R * R) / (12 * d); }
/* a patch thinned as the callers do (mk_decimate on the open surface, its rim moving only along itself): PV, PT
   malloc'd */
static int thin(const mk_gmesh *P, double h, double **PV, int **PT, int *npv, int *npt)
{
    mk_dec_opts o;
    mk_dec_defaults(&o, 1);
    o.dist = 0.2 * h, o.maxedge = 0.1, o.minalt = 1e-4 * h;
    unsigned char *alive = malloc((size_t)P->nv);
    int *To = malloc(sizeof(int) * 3 * (size_t)P->nt), nto = 0, *map = malloc(sizeof(int) * (size_t)P->nv), kv = 0;
    double *V = malloc(sizeof(double) * 3 * (size_t)P->nv);
    mk_dec_report dr;
    const int rd = mk_decimate(P->V, P->nv, P->T, P->nt, &o, alive, To, &nto, &dr);
    for (int i = 0; !rd && i < P->nv; i++)
        if (alive[i]) map[i] = kv, memcpy(V + 3 * kv, P->V + 3 * i, 3 * sizeof *V), kv++;
    for (int t = 0; !rd && t < 3 * nto; t++) To[t] = map[To[t]];
    free(alive), free(map);
    *PV = V, *PT = To, *npv = kv, *npt = nto;
    return rd;
}
/* the field carved by the ball at the pole, on S's own distance grid */
static int carved(const tm_mesh *S, double h, mk_sdf *g)
{
    const int r = mk_sdf_mesh(S->V, S->nv, S->T, S->nt, h, 0.06, 0, g);
    for (size_t q = 0; !r && q < (size_t)g->n[0] * g->n[1] * g->n[2]; q++) {
        const int i = (int)(q % (size_t)g->n[0]), j = (int)(q / (size_t)g->n[0] % (size_t)g->n[1]), k = (int)(q / ((size_t)g->n[0] * g->n[1]));
        const double x = g->lo[0] + i * h, y = g->lo[1] + j * h, z = g->lo[2] + k * h - 1;
        g->d[q] = (float)fmax(g->d[q], BALL - sqrt(x * x + y * y + z * z));
    }
    return r;
}

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
        long zs, zt;
        const long ps = ks ? pinched(s.nv, s.T, s.nt) : -1, pt = kt ? pinched(t.nv, t.T, t.nt) : -1;
        const double ms = ks ? min2a(s.V, s.T, s.nt, 1e-6 * h * h, &zs) : 0, mt = kt ? min2a(t.V, t.T, t.nt, 1e-6 * h * h, &zt) : 0;
        CHECK(ks && cs == 2 && ps == 0 && zs == 0 && kt && ct == 0 && pt == 0 && zt == 0,
              "marching tetrahedra on a unit sphere's and a torus's fields, h 0.05: closed, Euler %d and %d, pinched vertices %ld and %ld, smallest |2A| %.2g h^2 and %.2g h^2", cs, ct, ps, pt, ms / (h * h),
              mt / (h * h));
        /* the surface sinks by about h^2 H / 6 on average (H the mean curvature; measured from R/h 5 to 40 on offset grids, the
           law within 6 %): a ball's volume low by (h/R)^2 / 2, a tube's by (h/r)^2 / 6. Up to 1.6 times that on this grid,
           aligned with the part: low, by under twice the law */
        const double ls = -0.5 * (h / 1) * (h / 1), lt = -(h / 0.35) * (h / 0.35) / 6;
        CHECK(ks && vs < es && vs / es - 1 > 2 * ls && kt && vt < et && vt / et - 1 > 2 * lt,
              "the volume: sphere %.4f against %.4f (%.3f %%, the law %.3f %%), torus %.4f against %.4f (%.3f %%, the law %.3f %%)", vs, es, 100 * (vs / es - 1), 100 * ls, vt, et, 100 * (vt / et - 1),
              100 * lt);
        tm_write("mtet_sphere", s.V, s.nv, s.T, s.nt, NULL), tm_write("mtet_torus", t.V, t.nv, t.T, t.nt, NULL);
        { /* the unit sphere in a box from -0.5 to 1.5: cut by three of its faces, open there and only there */
            const int nb[3] = {41, 41, 41};
            const double lb[3] = {-0.5, -0.5, -0.5}, hb[3] = {lb[0] + 40 * h, lb[1] + 40 * h, lb[2] + 40 * h};
            mk_gmesh o;
            float *B = malloc(sizeof(float) * 41 * 41 * 41);
            for (int k = 0; k < 41; k++)
                for (int j = 0; j < 41; j++)
                    for (int i = 0; i < 41; i++) B[(k * 41 + j) * 41 + i] = F[((k + 20) * 61 + (j + 20)) * 61 + (i + 20)];
            const int rc = mk_mtet(B, nb, lb, h, &o);
            long open = 0, off = 0, zb = 0;
            if (!rc) {
                long long *e = malloc(sizeof *e * 3 * (size_t)o.nt);
                for (int u = 0; u < o.nt; u++)
                    for (int k = 0; k < 3; k++) e[3 * u + k] = (long long)o.T[3 * u + k] * o.nv + o.T[3 * u + (k + 1) % 3];
                for (int u = 0; u < 3 * o.nt; u++) {
                    const long long rv = (e[u] % o.nv) * o.nv + e[u] / o.nv;
                    int f = 0;
                    for (int w = 0; w < 3 * o.nt && !f; w++) f = e[w] == rv;
                    if (f) continue;
                    open++;
                    const double *a = o.V + 3 * (e[u] / o.nv), *b = o.V + 3 * (e[u] % o.nv);
                    int face = 0;
                    for (int x = 0; x < 3; x++) face |= (a[x] == lb[x] && b[x] == lb[x]) || (a[x] == hb[x] && b[x] == hb[x]);
                    off += !face;
                }
                free(e);
                min2a(o.V, o.T, o.nt, 1e-6 * h * h, &zb);
            }
            CHECK(rc == 0 && open > 100 && off == 0 && zb == 0, "the sphere cut by a box from -0.5 to 1.5: %ld open edges, %ld of them off the box's faces; %ld zero-area triangles", open, off, zb);
            mk_gmesh_free(&o), free(B);
        }
        { /* a field 0 on a plane of nodes (a node at 0 counts outside): each crossing at that node, held off it */
            const int np[3] = {12, 12, 12};
            const double lp[3] = {0, 0, 0};
            float *P = malloc(sizeof(float) * 12 * 12 * 12);
            for (int k = 0; k < 12; k++)
                for (int j = 0; j < 12; j++)
                    for (int i = 0; i < 12; i++) P[(k * 12 + j) * 12 + i] = (float)(k - 6);
            mk_gmesh o;
            const int rp = mk_mtet(P, np, lp, h, &o);
            long zp = -1;
            const double mp = rp ? 0 : min2a(o.V, o.T, o.nt, 1e-6 * h * h, &zp);
            CHECK(rp == 0 && o.nt > 100 && zp == 0, "a field 0 on a plane of nodes: %d triangles, %ld of zero area, the smallest |2A| %.3g h^2", rp ? -1 : o.nt, zp, mp / (h * h));
            mk_gmesh_free(&o), free(P);
        }
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
    int r = carved(&S, h, &g);
    const size_t N = (size_t)g.n[0] * g.n[1] * g.n[2];
    mk_graft G;
    const int ri = r ? -1 : mk_graft_init(&G, S.V, S.nv, S.T, S.nt, psi, NULL, 2 * h, 1e-6);
    mk_gmesh P;
    const int rp = ri ? -1 : mk_graft_patch(&G, g.d, g.n, g.lo, g.h, &P);
    CHECK(r == 0 && ri == 0 && rp == 0 && P.nt > 1000, "the patch: mk_sdf_mesh's grid at h 0.02 carved, its surface by marching tetrahedra kept in the cap: %d vertices, %d triangles", rp ? -1 : P.nv,
          rp ? -1 : P.nt);
    /* between patch and stitch: the patch thinned as an open surface, its rim moving only along itself */
    double *PV = NULL;
    int *To = NULL, kv = 0, nto = 0;
    const int rd = thin(&P, h, &PV, &To, &kv, &nto);
    CHECK(rd == 0 && kv < P.nv / 2, "the patch thinned by mk_decimate: %d -> %d vertices", P.nv, kv);
    mk_gmesh O;
    mk_graft_report gr;
    const int rs = mk_graft_stitch(&G, PV, kv, To, nto, 2.5 * h, &O, &gr);
    int chi = 0;
    const int ok = rs == 0 && closed(O.V, O.nv, O.T, O.nt, &chi);
    long zo = 0;
    const long po = ok ? pinched(O.nv, O.T, O.nt) : -1, xo = ok ? xpairs(O.V, O.T, O.nt) : -1;
    const double mo = ok ? min2a(O.V, O.T, O.nt, 1e-6 * h * h, &zo) : 0;
    CHECK(ok && chi == 2 && gr.dup == 0 && gr.lone == 0 && gr.pieces == 1 && gr.reason == 0 && gr.faults == 0 && po == 0 && zo == 0 && xo == 0 && gr.zero_area == 0 && gr.cross1 == 0,
          "stitched: rc %d, no fault (%#x); closed and consistently wound, Euler %d, %d piece, %ld pinched; %ld rim vertices put into S's rim, %ld welded, largest move %.5f (h %.2f); crossing pairs %ld "
          "(the stitch's count %ld, S's %ld), smallest |2A| %.2g h^2",
          rs, gr.faults, chi, gr.pieces, po, gr.inserted, gr.welded, gr.snap, h, xo, gr.cross1, gr.cross0, mo / (h * h));
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
         out of order and are merged, and its inside near the rim follows the rim onto S's. Closed, one piece, rc 0, with no
         crossing pair, on this S and on a coarse one (3 subdivisions, edges about 6 h) */
        mk_gmesh O2 = {0};
        mk_graft_report r2 = {0}, r3 = {0};
        const int s2 = mk_graft_stitch(&G, P.V, P.nv, P.T, P.nt, 2.5 * h, &O2, &r2);
        int c2 = 0, c3 = 0;
        const int k2 = s2 == 0 && closed(O2.V, O2.nv, O2.T, O2.nt, &c2);
        tm_mesh S3 = {0};
        tm_sphere(&S3, 1, c0, 3);
        mk_sdf g3;
        int s3 = carved(&S3, h, &g3), k3 = 0;
        mk_graft G3;
        mk_gmesh P3 = {0}, O3 = {0};
        const int i3 = s3 ? -1 : mk_graft_init(&G3, S3.V, S3.nv, S3.T, S3.nt, psi, NULL, 2 * h, 1e-6);
        const int p3 = s3 || i3 ? -9 : mk_graft_patch(&G3, g3.d, g3.n, g3.lo, g3.h, &P3);
        if (!s3) s3 = i3 ? -9 : p3 ? -8 : mk_graft_stitch(&G3, P3.V, P3.nv, P3.T, P3.nt, 2.5 * h, &O3, &r3);
        k3 = s3 == 0 && closed(O3.V, O3.nv, O3.T, O3.nt, &c3);
        const double v3 = tm_volume(S3.V, S3.T, S3.nt) - lensvol(1, BALL, 1);
        const long x2 = k2 ? xpairs(O2.V, O2.T, O2.nt) : -1, x3 = k3 ? xpairs(O3.V, O3.T, O3.nt) : -1;
        CHECK(k2 && c2 == 2 && r2.dup == 0 && r2.lone == 0 && r2.pieces == 1 && r2.merged > 0 && r2.faults == 0 && x2 == 0 && fabs(r2.vol1 / vexp - 1) < 3e-3 && k3 && c3 == 2 && r3.dup == 0 &&
                  r3.lone == 0 && r3.pieces == 1 && r3.faults == 0 && x3 == 0 && fabs(r3.vol1 / v3 - 1) < 1e-2,
              "the patch not thinned (%d and %d triangles): rc %d and %d, closed, Euler %d and %d, dup %ld and %ld, lone %ld and %ld; %ld and %ld rim vertices out of order merged; volume %.3f %% and %.3f %% off; crossing pairs %ld and %ld (the stitch's count %ld and %ld, S's %ld and %ld)",
              P.nt, P3.nt, s2, s3, c2, c3, r2.dup, r3.dup, r2.lone, r3.lone, r2.merged, r3.merged, 100 * (r2.vol1 / vexp - 1), 100 * (r3.vol1 / v3 - 1), x2, x3, r2.cross1, r3.cross1, r2.cross0, r3.cross0);
        { /* mk_graft_refine: the coarse S refined where the seam cuts it (edges to 2 h within 2 h of it), then grafted
             as before, the patch thinned */
            mk_gmesh R4 = {0}, P4 = {0}, O4 = {0};
            mk_graft G4;
            mk_graft_report r4 = {0};
            const int rr = p3 ? -9 : mk_graft_refine(S3.V, S3.nv, S3.T, S3.nt, psi, NULL, 2 * h, 2 * h, &R4);
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
            double *V4 = NULL;
            int *T4 = NULL, n4v = 0, n4t = 0;
            const int s4 = i4 ? -9 : mk_graft_patch(&G4, g3.d, g3.n, g3.lo, g3.h, &P4) ? -8 : thin(&P4, h, &V4, &T4, &n4v, &n4t) ? -7 : mk_graft_stitch(&G4, V4, n4v, T4, n4t, 2.5 * h, &O4, &r4);
            const int k4 = s4 == 0 && closed(O4.V, O4.nv, O4.T, O4.nt, &c4);
            const long x4 = k4 ? xpairs(O4.V, O4.T, O4.nt) : -1;
            CHECK(kr && cr == 2 && left == 0 && kept == away && away > 500 && fabs(vr / vs3 - 1) < 1e-12 && k4 && c4 == 2 && x4 == 0 && fabs(r4.vol1 / v3 - 1) < 1e-2,
                  "mk_graft_refine on the coarse S: %d -> %d triangles in %d passes, closed, Euler %d, volume %.1e off, %ld long edges left by the seam, %ld of %ld triangles away from it as they were; grafted with the patch thinned: rc %d, Euler %d, rim %d -> %d edges, %ld crossing pairs",
                  S3.nt, R4.nt, rr, cr, vr / vs3 - 1, left, kept, away, s4, c4, r3.rim_s, r4.rim_s, x4);
            if (!i4) mk_graft_free(&G4);
            mk_gmesh_free(&R4), mk_gmesh_free(&P4), mk_gmesh_free(&O4), free(V4), free(T4);
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
        CHECK(hole >= 0 && sh != 0 && rh.reason != 0 && rh.loops_s == 1 && rh.loops_p == 2 && rh.unpaired == 1,
              "a patch with a one-triangle hole by the seam: not passed (rc %d, reason %d), rim loops S %d, patch %d, one without a partner (%d)", sh, rh.reason, rh.loops_s, rh.loops_p, rh.unpaired);
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
           order); then radius 1.5 with one triangle in twice, inverted and open. The stitch returns 2 when the kept
           piece is wound against S or such pieces hold over 1 % of its volume, before 1 (open): so 2, 2, 0 and 2 */
        const double rad[4] = {1.5, 0.5, 0.05, 1.5};
        int sb[4];
        mk_graft_report rq[4];
        for (int c = 0; c < 4; c++) {
            tm_mesh B = {0}, bs = {0};
            tm_sphere(&bs, rad[c], c0, 2);
            for (int t = 0; t < bs.nt; t++) { const int x = bs.T[3 * t + 1]; bs.T[3 * t + 1] = bs.T[3 * t + 2], bs.T[3 * t + 2] = x; }
            for (int i = 0; i < kv; i++) tm_vert(&B, PV[3 * i], PV[3 * i + 1], PV[3 * i + 2]);
            for (int t = 0; t < nto; t++) tm_tri(&B, To[3 * t], To[3 * t + 1], To[3 * t + 2]);
            for (int i = 0; i < bs.nv; i++) tm_vert(&B, bs.V[3 * i], bs.V[3 * i + 1], bs.V[3 * i + 2]);
            for (int t = 0; t < bs.nt; t++) tm_tri(&B, kv + bs.T[3 * t], kv + bs.T[3 * t + 1], kv + bs.T[3 * t + 2]);
            if (c == 3) tm_tri(&B, kv + bs.T[0], kv + bs.T[1], kv + bs.T[2]); /* a triangle twice: its edges used twice */
            sb[c] = mk_graft_stitch(&G, B.V, B.nv, B.T, B.nt, 2.5 * h, &Ob, &rq[c]);
            mk_gmesh_free(&Ob), tm_free(&B), tm_free(&bs);
        }
        CHECK(sb[0] == 2 && rq[0].vol1 < 0 && fabs(rq[0].dropped_vol / gr.vol1 - 1) < 1e-9 && sb[1] == 2 && fabs(rq[1].vol1 / gr.vol1 - 1) < 1e-9 && rq[1].inverted_vol < -0.01 * rq[1].vol1 &&
                  sb[2] == 0 && fabs(rq[2].vol1 / gr.vol1 - 1) < 1e-9 && rq[2].inverted == 1 && rq[2].inverted_vol < 0 && rq[2].inverted_vol > -1e-3 * rq[2].vol1,
              "a closed piece wound inward beside the body, radius 1.5, 0.5, 0.05: rc %d, %d, %d (2: kept %.3f, or dropped %.4f against S over 1 %% of the body; the speck, %.2g, reported only); %d, %d, %d inverted",
              sb[0], sb[1], sb[2], rq[0].vol1, rq[1].inverted_vol, rq[2].inverted_vol, rq[0].inverted, rq[1].inverted, rq[2].inverted);
        CHECK(sb[3] == MK_GRAFT_INVERTED && rq[3].dup > 0 && rq[3].faults >> MK_GRAFT_OPEN & 1 && rq[3].faults >> MK_GRAFT_INVERTED & 1,
              "the kept piece wound inward and open (a triangle in twice): rc %d, 2 before 1 (faults %#x, dup %ld)", sb[3], rq[3].faults, rq[3].dup);
        tm_free(&sp), free(Tr), free(Fn), free(Th), mk_gmesh_free(&Pi), mk_gmesh_free(&Oi), mk_gmesh_free(&Ob), mk_gmesh_free(&Oh);
        if (!ii) mk_graft_free(&Gi);
        if (!i3) mk_graft_free(&G3);
        mk_gmesh_free(&O2), mk_gmesh_free(&O3), mk_gmesh_free(&P3), mk_sdf_free(&g3), tm_free(&S3);
    }
    { /* a split: beside the body a closed outward piece the stitch drops, of radius 0.3 (2.7 % of S's volume: rc 3)
         and 0.05 (0.01 %: rc 0) */
        const double rad[2] = {0.3, 0.05}, cs[3] = {0, 0, -3};
        int sx[2];
        mk_graft_report rx[2];
        for (int c = 0; c < 2; c++) {
            tm_mesh B = {0}, bs = {0};
            tm_sphere(&bs, rad[c], cs, 2);
            for (int i = 0; i < kv; i++) tm_vert(&B, PV[3 * i], PV[3 * i + 1], PV[3 * i + 2]);
            for (int t = 0; t < nto; t++) tm_tri(&B, To[3 * t], To[3 * t + 1], To[3 * t + 2]);
            for (int i = 0; i < bs.nv; i++) tm_vert(&B, bs.V[3 * i], bs.V[3 * i + 1], bs.V[3 * i + 2]);
            for (int t = 0; t < bs.nt; t++) tm_tri(&B, kv + bs.T[3 * t], kv + bs.T[3 * t + 1], kv + bs.T[3 * t + 2]);
            mk_gmesh Ox;
            sx[c] = mk_graft_stitch(&G, B.V, B.nv, B.T, B.nt, 2.5 * h, &Ox, &rx[c]);
            mk_gmesh_free(&Ox), tm_free(&B), tm_free(&bs);
        }
        CHECK(sx[0] == MK_GRAFT_SPLIT && rx[0].ndropped == 1 && rx[0].dropped_abs > 0.01 * rx[0].vol0 && fabs(rx[0].vol1 / gr.vol1 - 1) < 1e-9 && sx[1] == 0 && rx[1].ndropped == 1,
              "a closed piece dropped beside the body, radius 0.3 and 0.05: rc %d and %d (3: dropped %.4f, %.2f %% of S's volume; %.2g, %.3f %%)", sx[0], sx[1], rx[0].dropped_abs,
              100 * rx[0].dropped_abs / rx[0].vol0, rx[1].dropped_abs, 100 * rx[1].dropped_abs / rx[1].vol0);
    }
    { /* a zero-area triangle in the patch, away from its rim: one triangle split at the middle of an edge, the third
         piece flat on that edge. The stitch returns 4 */
        int tz = -1;
        for (int t = 0; t < nto && tz < 0; t++) {
            int deep = 1;
            for (int k = 0; k < 3; k++) deep &= psi(PV + 3 * To[3 * t + k], NULL) < -0.2;
            if (deep) tz = t;
        }
        tm_mesh B = {0};
        for (int i = 0; i < kv; i++) tm_vert(&B, PV[3 * i], PV[3 * i + 1], PV[3 * i + 2]);
        int sz = -9;
        mk_graft_report rz = {0};
        if (tz >= 0) {
            const int a = To[3 * tz], b = To[3 * tz + 1], c = To[3 * tz + 2];
            const int m = tm_vert(&B, 0.5 * (PV[3 * a] + PV[3 * b]), 0.5 * (PV[3 * a + 1] + PV[3 * b + 1]), 0.5 * (PV[3 * a + 2] + PV[3 * b + 2]));
            for (int t = 0; t < nto; t++)
                if (t != tz) tm_tri(&B, To[3 * t], To[3 * t + 1], To[3 * t + 2]);
            tm_tri(&B, a, m, c), tm_tri(&B, m, b, c), tm_tri(&B, a, b, m);
            mk_gmesh Oz;
            sz = mk_graft_stitch(&G, B.V, B.nv, B.T, B.nt, 2.5 * h, &Oz, &rz);
            mk_gmesh_free(&Oz);
        }
        CHECK(sz == MK_GRAFT_ZERO && rz.zero_area == 1 && rz.dup == 0 && rz.lone == 0 && rz.pinched == 0, "a patch triangle of zero area: rc %d (4), %ld new triangles of zero area", sz, rz.zero_area);
        tm_free(&B);
    }
    { /* the seam just off one of S's vertices (the one nearest it, moved onto it, then nudged 1e-6 off): a cut of its
         triangles within 1e-6 of it would leave slivers of zero area. Grafted as above: rc 0, none under 1e-6 h^2 */
        int vs = 0;
        for (int i = 1; i < S.nv; i++)
            if (fabs(psi(S.V + 3 * i, NULL)) < fabs(psi(S.V + 3 * vs, NULL))) vs = i;
        double off = psi(S.V + 3 * vs, NULL);
        mk_graft Gs;
        mk_gmesh Ps = {0}, Os = {0};
        mk_graft_report rv = {0};
        double *Vs = NULL;
        int *Ts = NULL, nsv = 0, nst = 0;
        const int is = mk_graft_init(&Gs, S.V, S.nv, S.T, S.nt, psi_off, &off, 2 * h, 1e-6);
        const int ss = is ? -9 : mk_graft_patch(&Gs, g.d, g.n, g.lo, g.h, &Ps) ? -8 : thin(&Ps, h, &Vs, &Ts, &nsv, &nst) ? -7 : mk_graft_stitch(&Gs, Vs, nsv, Ts, nst, 2.5 * h, &Os, &rv);
        long zs = -1;
        int cs = 0;
        const int ks = ss == 0 && closed(Os.V, Os.nv, Os.T, Os.nt, &cs);
        const double ms = ks ? min2a(Os.V, Os.T, Os.nt, 1e-6 * h * h, &zs) : 0;
        CHECK(ks && cs == 2 && zs == 0 && rv.zero_area == 0 && rv.faults == 0, "the seam through S's vertex %d (psi there %.1e before the nudge): rc %d, closed, Euler %d, %ld triangles of zero area, the smallest |2A| %.2g h^2",
              vs, psi_off(S.V + 3 * vs, &off), ss, cs, zs, ms / (h * h));
        if (!is) mk_graft_free(&Gs);
        mk_gmesh_free(&Ps), mk_gmesh_free(&Os), free(Vs), free(Ts);
    }
    { /* a pinched vertex: a small closed ball inside the body touching the patch at one of its vertices, sharing it.
         Closed edge by edge and one piece, but not a manifold: rc 1 */
        int vp = -1;
        for (int i = 0; i < kv && vp < 0; i++)
            if (psi(PV + 3 * i, NULL) < -0.2) vp = i;
        tm_mesh B = {0}, bs = {0};
        int sp = -9;
        mk_graft_report rpn = {0};
        if (vp >= 0) {
            const double *x = PV + 3 * vp, C[3] = {0, 0, 1}, rs_ = 0.004;
            double nin[3], l = 0;
            for (int k = 0; k < 3; k++) nin[k] = x[k] - C[k], l += nin[k] * nin[k];
            l = sqrt(l);
            const double cb[3] = {x[0] + rs_ * nin[0] / l, x[1] + rs_ * nin[1] / l, x[2] + rs_ * nin[2] / l};
            tm_sphere(&bs, rs_, cb, 1);
            int near = 0;
            for (int i = 1; i < bs.nv; i++) {
                double di = 0, dn = 0;
                for (int k = 0; k < 3; k++) di += (bs.V[3 * i + k] - x[k]) * (bs.V[3 * i + k] - x[k]), dn += (bs.V[3 * near + k] - x[k]) * (bs.V[3 * near + k] - x[k]);
                if (di < dn) near = i;
            }
            for (int i = 0; i < kv; i++) tm_vert(&B, PV[3 * i], PV[3 * i + 1], PV[3 * i + 2]);
            for (int t = 0; t < nto; t++) tm_tri(&B, To[3 * t], To[3 * t + 1], To[3 * t + 2]);
            int *id = malloc(sizeof(int) * (size_t)bs.nv);
            for (int i = 0; i < bs.nv; i++) id[i] = i == near ? vp : tm_vert(&B, bs.V[3 * i], bs.V[3 * i + 1], bs.V[3 * i + 2]);
            for (int t = 0; t < bs.nt; t++) tm_tri(&B, id[bs.T[3 * t]], id[bs.T[3 * t + 1]], id[bs.T[3 * t + 2]]);
            free(id);
            mk_gmesh Op;
            sp = mk_graft_stitch(&G, B.V, B.nv, B.T, B.nt, 2.5 * h, &Op, &rpn);
            const long pp = sp >= 0 ? pinched(Op.nv, Op.T, Op.nt) : -1;
            CHECK(sp == MK_GRAFT_OPEN && rpn.pinched == 1 && pp == 1 && rpn.dup == 0 && rpn.lone == 0 && rpn.pieces == 1,
                  "a ball touching the patch at one vertex inside the body: rc %d (1), %d piece, dup %ld, lone %ld, pinched %ld (the test's count %ld)", sp, rpn.pieces, rpn.dup, rpn.lone, rpn.pinched, pp);
            mk_gmesh_free(&Op);
        } else CHECK(0, "a pinched vertex: no patch vertex deep in the cap");
        tm_free(&B), tm_free(&bs);
    }
    { /* a crossing: a patch edge by the rim moved past the seam, one end out of S and one in, so it pierces S's kept
         part. Closed, one piece, but crossing S: rc 5, the stitch's count of crossing pairs the whole result's (S has
         none) */
        int vx = -1, vy = -1;
        for (int i = 0; i < kv && vx < 0; i++) {
            const double s = psi(PV + 3 * i, NULL);
            if (s < -1.5 * h && s > -3 * h) vx = i;
        }
        for (int t = 0; vx >= 0 && vy < 0 && t < nto; t++) /* a neighbour of vx's, as deep in the cap */
            for (int k = 0; k < 3; k++)
                if (To[3 * t + k] == vx) {
                    const int w = To[3 * t + (k + 1) % 3];
                    if (psi(PV + 3 * w, NULL) < -1.5 * h) vy = w;
                }
        int sxc = -9;
        mk_graft_report rxc = {0};
        long xx = -1;
        if (vy >= 0) {
            tm_mesh B = {0};
            for (int i = 0; i < kv; i++) tm_vert(&B, PV[3 * i], PV[3 * i + 1], PV[3 * i + 2]);
            for (int t = 0; t < nto; t++) tm_tri(&B, To[3 * t], To[3 * t + 1], To[3 * t + 2]);
            for (int q = 0; q < 2; q++) {
                double *x = B.V + 3 * (q ? vy : vx);
                const double th = acos(x[2] / sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2])) + 0.2, ph = atan2(x[1], x[0]), rr = q ? 0.85 : 1.1;
                x[0] = rr * sin(th) * cos(ph), x[1] = rr * sin(th) * sin(ph), x[2] = rr * cos(th);
            }
            mk_gmesh Ox;
            sxc = mk_graft_stitch(&G, B.V, B.nv, B.T, B.nt, 2.5 * h, &Ox, &rxc);
            if (sxc >= 0) xx = xpairs(Ox.V, Ox.T, Ox.nt);
            mk_gmesh_free(&Ox), tm_free(&B);
        }
        CHECK(sxc == MK_GRAFT_CROSS && rxc.cross0 == 0 && rxc.cross1 > 0 && rxc.cross1 == xx && rxc.dup == 0 && rxc.lone == 0 && rxc.pieces == 1,
              "a patch edge moved through S beyond the seam: rc %d (5), crossing pairs %ld -> %ld (the whole result's %ld), closed, %d piece", sxc, rxc.cross0, rxc.cross1, xx, rxc.pieces);
    }
    {
        mk_gmesh z;
        mk_graft_report zr;
        int k = 0, nr = 0, na = 0;
        const int Tb[3] = {0, 1, 99999999};
#define REF(call) (k++, nr += (call) == -1, na += zr.reason == MK_GRAFT_E_ARGS)
        REF(mk_graft_stitch(&G, PV, kv, To, nto, -1, &z, &zr));
        REF(mk_graft_stitch(&G, PV, kv, To, nto, NAN, &z, &zr));
        REF(mk_graft_stitch(&G, PV, kv, Tb, 1, 2.5 * h, &z, &zr));
#undef REF
        const int rsn = mk_graft_stitch(&G, PV, kv, To, nto, 1e-9, &z, &zr); /* the rims too far apart for that snap */
        k++, nr += rsn == -1;
        const int snapr = zr.reason == MK_GRAFT_E_SNAP && zr.snap > 1e-9 && zr.snap < 2.5 * h;
        mk_graft G2;
#define REF(call) (k++, nr += (call) == -1)
        REF(mk_graft_init(&G2, S.V, S.nv, S.T, S.nt, NULL, NULL, h, 0));
        REF(mk_graft_init(&G2, S.V, S.nv, S.T, S.nt, psi, NULL, -1, 0));
        REF(mk_graft_init(&G2, S.V, 2, S.T, S.nt, psi, NULL, h, 0));
#undef REF
        CHECK(nr == k && na == 3 && snapr,
              "graft refusals: %d of %d (snap negative, NaN or too small for the rims, an index out of range; no psi, near negative, an index past nv); reasons: %d of 3 bad arguments, the snap %d (%.2g, its distance)", nr,
              k, na, zr.reason, zr.snap);
    }
    mk_gmesh_free(&O), mk_gmesh_free(&P), free(To), free(PV), mk_graft_free(&G), mk_sdf_free(&g), tm_free(&S);
    return wrong;
}
