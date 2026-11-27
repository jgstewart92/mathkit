/* t_cut: mk_cut.h. mk_cut_loops on a closed U tube (genus 0, rounded ends): three cuts, the left leg's plane also crossing
   the right leg (that loop stays whole), the pieces closed and outward, their volumes summing to the whole, each cap shared
   wound opposite, each cap vertex on its plane, the source map reproducing every new vertex bit for bit; a cut 2 um off a
   ring of vertices, which snaps; a cone whose plane must move one step to pass; a flat tube's cap shape. mk_cut_clip on a
   ball (one plane, then three), on a torus through its axis (two loops, two halves) and through its equator (nested,
   refused). The probe on a tube inside a tube. Each refusal. With a folder argument it writes OBJ files for meshcheck. */
#include "../include/mk_cut.h"
#include "t.h"
#include "tmesh.h"
#include <stdint.h>
#define K 48
#define M 400
static double V[3 * (K * M + 2)];
static int T[3 * (2 * K * M + 2 * K)];
static int nv, nt;
static void centre(double s, double *c, double *tg)
{
    /* s in [0,3]: up the left leg (x=-0.05, z 0 -> 0.1), round the bend (radius 0.05, centre (0,0,0.1)), down the right leg */
    if (s < 1) c[0] = -0.05, c[1] = 0, c[2] = 0.1 * s, tg[0] = 0, tg[1] = 0, tg[2] = 1;
    else if (s < 2) { const double a = T_PI * (s - 1); c[0] = -0.05 * cos(a), c[1] = 0, c[2] = 0.1 + 0.05 * sin(a); tg[0] = sin(a), tg[1] = 0, tg[2] = cos(a); }
    else c[0] = 0.05, c[1] = 0, c[2] = 0.1 - 0.1 * (s - 2), tg[0] = 0, tg[1] = 0, tg[2] = -1;
}
static double vol(const double *VV, const int *TT, int n) { return tm_volume(VV, TT, n); }
static int cmp64(const void *a, const void *b) { const int64_t x = *(const int64_t *)a, y = *(const int64_t *)b; return x < y ? -1 : x > y; }
/* closed and consistently wound: every directed edge once, and its reverse once */
static int oriented(const int *TT, int n, int nvv)
{
    int64_t *k = malloc(sizeof *k * 3 * (size_t)n);
    for (int t = 0; t < n; t++) for (int e = 0; e < 3; e++) k[3 * t + e] = (int64_t)TT[3 * t + e] * nvv + TT[3 * t + (e + 1) % 3];
    qsort(k, 3 * (size_t)n, sizeof *k, cmp64);
    int ok = n > 0;
    for (size_t i = 0; ok && i < 3 * (size_t)n; i++) {
        if (i + 1 < 3 * (size_t)n && k[i + 1] == k[i]) ok = 0;
        const int64_t r = (k[i] % nvv) * nvv + k[i] / nvv;
        if (!bsearch(&r, k, 3 * (size_t)n, sizeof *k, cmp64)) ok = 0;
    }
    free(k);
    return ok;
}
/* the smallest twice-area over squared edge lengths of the triangles (0 for one of zero area) */
static double minshape(const double *VV, const int *TT, int n)
{
    double q = INFINITY;
    for (int t = 0; t < n; t++) {
        const double *a = VV + 3 * TT[3 * t], *b = VV + 3 * TT[3 * t + 1], *c = VV + 3 * TT[3 * t + 2];
        const double e1[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, e2[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]}, e3[3] = {c[0] - b[0], c[1] - b[1], c[2] - b[2]};
        const double x = e1[1] * e2[2] - e1[2] * e2[1], y = e1[2] * e2[0] - e1[0] * e2[2], z = e1[0] * e2[1] - e1[1] * e2[0];
        const double l2 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2] + e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2] + e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2];
        q = fmin(q, sqrt(x * x + y * y + z * z) / l2);
    }
    return q;
}
/* every piece of o: oriented, outward; their volumes' sum against the split surface's (relative: the caps cancel); *dv
   the split surface's against the input's vtot (relative: the snapped vertices' moves); written as NAME0, NAME1, ... */
static int pieces(const mk_cut_out *o, double vtot, double *rel, double *dv, const char *name)
{
    int ok = 1;
    double vs = 0;
    unsigned char *used = calloc((size_t)o->nv, 1);
    for (int g = 0; g < o->nreg; g++) {
        int *TT, m;
        if (mk_cut_piece(o, g, &TT, &m)) { ok = 0; continue; }
        const double v = vol(o->V, TT, m);
        ok &= oriented(TT, m, o->nv) && v > 0;
        vs += v;
        if (name) {
            char nm[64];
            memset(used, 0, (size_t)o->nv);
            for (int i = 0; i < 3 * m; i++) used[TT[i]] = 1;
            snprintf(nm, sizeof nm, "%s%d", name, g);
            tm_write(nm, o->V, o->nv, TT, m, used);
        }
        free(TT);
    }
    free(used);
    const double vsplit = vol(o->V, o->T, o->nt);
    *rel = fabs(vs - vsplit) / vsplit, *dv = fabs(vsplit - vtot) / vtot;
    return ok;
}
/* each cap triangle is in its loop's two pieces, as it is in the - piece and reversed in the + piece */
static int caps_shared(const mk_cut_out *o)
{
    int ok = o->ncapt > 0;
    for (int s = 0; ok && s < 2; s++)
        for (int l = 0; ok && l < o->nloop; l++) {
            int *TT, m;
            if (mk_cut_piece(o, o->side[2 * l + s], &TT, &m)) return 0;
            int64_t *k = malloc(sizeof *k * (size_t)(m ? m : 1));
            for (int t = 0; t < m; t++) { /* each triangle by its corners rotated smallest first */
                const int *c = TT + 3 * t, r = c[0] < c[1] ? (c[0] < c[2] ? 0 : 2) : (c[1] < c[2] ? 1 : 2);
                k[t] = ((int64_t)c[r] * o->nv + c[(r + 1) % 3]) * o->nv + c[(r + 2) % 3];
            }
            qsort(k, (size_t)m, sizeof *k, cmp64);
            for (int t = 0; ok && t < o->ncapt; t++) {
                if (o->capk[t] != l) continue;
                const int *c = o->capT + 3 * t;
                const int w[3] = {c[0], s ? c[2] : c[1], s ? c[1] : c[2]}, r = w[0] < w[1] ? (w[0] < w[2] ? 0 : 2) : (w[1] < w[2] ? 1 : 2);
                const int64_t q = ((int64_t)w[r] * o->nv + w[(r + 1) % 3]) * o->nv + w[(r + 2) % 3];
                if (!bsearch(&q, k, (size_t)m, sizeof *k, cmp64)) ok = 0;
            }
            free(k), free(TT);
        }
    return ok;
}
/* the largest distance of a cap vertex from its plane (through Pc, normal N) */
static double capoff(const mk_cut_out *o, const double *N)
{
    double w = 0;
    for (int t = 0; t < o->ncapt; t++) {
        const int k = o->loopk[o->capk[t]];
        const double *n = N + 3 * k, nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        for (int j = 0; j < 3; j++) {
            const double *x = o->V + 3 * o->capT[3 * t + j];
            w = fmax(w, fabs((x[0] - o->Pc[3 * k]) * n[0] + (x[1] - o->Pc[3 * k + 1]) * n[1] + (x[2] - o->Pc[3 * k + 2]) * n[2]) / nn);
        }
    }
    return w;
}
/* every new vertex from its source: an edge of the input crossing its plane, t in [0.001, 0.999], and V[a] + t (V[b] -
   V[a]) its position bit for bit; *nnew the new vertices */
static int srcmap(const mk_cut_out *o, const double *VV, const int *TT, int ntt, int *nnew)
{
    int64_t *e = malloc(sizeof *e * 3 * (size_t)ntt);
    for (int t = 0; t < ntt; t++) for (int j = 0; j < 3; j++) { const int64_t a = TT[3 * t + j], b = TT[3 * t + (j + 1) % 3]; e[3 * t + j] = (a < b ? a : b) * o->nv0 + (a < b ? b : a); }
    qsort(e, 3 * (size_t)ntt, sizeof *e, cmp64);
    int ok = 1;
    *nnew = o->nv - o->nv0;
    for (int i = 0; i < *nnew; i++) {
        const int a = o->src[2 * i], b = o->src[2 * i + 1];
        const double t = o->frac[i];
        const int64_t key = (int64_t)(a < b ? a : b) * o->nv0 + (a < b ? b : a);
        if (a < 0 || b < 0 || a >= o->nv0 || b >= o->nv0 || !bsearch(&key, e, 3 * (size_t)ntt, sizeof *e, cmp64) || !(t >= 0.001 && t <= 0.999)) { ok = 0; continue; }
        for (int j = 0; j < 3; j++) if (VV[3 * a + j] + t * (VV[3 * b + j] - VV[3 * a + j]) != o->V[3 * (o->nv0 + i) + j]) ok = 0;
    }
    free(e);
    return ok;
}
/* a closed tube along z from 0 to 1, rings of k points of radius ra(z) (x) and rb(z) (y), capped by a point at each end */
static void tube(tm_mesh *m, int rings, int k, double r0, double slope, double flat)
{
    for (int i = 0; i < rings; i++) {
        const double z = (double)i / (rings - 1), r = r0 + slope * z;
        for (int j = 0; j < k; j++) { const double a = 2 * T_PI * j / k; tm_vert(m, r * cos(a), flat * r * sin(a), z); }
    }
    const int p0 = tm_vert(m, 0, 0, -r0), p1 = tm_vert(m, 0, 0, 1 + r0);
    for (int i = 0; i + 1 < rings; i++)
        for (int j = 0; j < k; j++) {
            const int a = i * k + j, b = i * k + (j + 1) % k, c = (i + 1) * k + (j + 1) % k, d = (i + 1) * k + j;
            tm_tri(m, a, b, c), tm_tri(m, a, c, d);
        }
    for (int j = 0; j < k; j++) tm_tri(m, p0, (j + 1) % k, j), tm_tri(m, p1, (rings - 1) * k + j, (rings - 1) * k + (j + 1) % k);
}
static void outward(tm_mesh *m)
{
    if (tm_volume(m->V, m->T, m->nt) < 0)
        for (int t = 0; t < m->nt; t++) { const int x = m->T[3 * t + 1]; m->T[3 * t + 1] = m->T[3 * t + 2], m->T[3 * t + 2] = x; }
}
int main(int argc, char **argv)
{
    tm_dir = argc > 1 ? argv[1] : NULL;
    mk_cut_opts om, o1;
    mk_cut_defaults(&om, 1e-3), mk_cut_defaults(&o1, 1);
    CHECK(om.snap == 5e-6 && om.shift == 1e-4 && om.step == 1e-5 && o1.snap == 0.005 && o1.shift == 0.1 && o1.step == 0.01,
          "defaults: snap %g, shift %g in steps of %g in metres (mm 1e-3); %g, %g, %g in mm", om.snap, om.shift, om.step, o1.snap, o1.shift, o1.step);
    const double r = 0.01;
    for (int i = 0; i < M; i++) {
        double c[3], tg[3];
        centre(3.0 * i / (M - 1), c, tg);
        const double u[3] = {0, 1, 0}, w[3] = {tg[1] * u[2] - tg[2] * u[1], tg[2] * u[0] - tg[0] * u[2], tg[0] * u[1] - tg[1] * u[0]};
        for (int j = 0; j < K; j++) {
            const double a = 2 * T_PI * j / K;
            for (int q = 0; q < 3; q++) V[3 * nv + q] = c[q] + r * (cos(a) * u[q] + sin(a) * w[q]);
            nv++;
        }
    }
    const int p0 = nv, p1 = nv + 1;
    { double c[3], tg[3]; centre(0, c, tg); V[3 * nv] = c[0], V[3 * nv + 1] = c[1], V[3 * nv + 2] = c[2] - r; nv++; centre(3, c, tg); V[3 * nv] = c[0], V[3 * nv + 1] = c[1], V[3 * nv + 2] = c[2] - r; nv++; }
    for (int i = 0; i + 1 < M; i++)
        for (int j = 0; j < K; j++) {
            const int a = i * K + j, b = i * K + (j + 1) % K, c = (i + 1) * K + (j + 1) % K, d = (i + 1) * K + j;
            T[3 * nt] = a, T[3 * nt + 1] = b, T[3 * nt + 2] = c, nt++;
            T[3 * nt] = a, T[3 * nt + 1] = c, T[3 * nt + 2] = d, nt++;
        }
    for (int j = 0; j < K; j++) {
        T[3 * nt] = p0, T[3 * nt + 1] = (j + 1) % K, T[3 * nt + 2] = j, nt++;
        T[3 * nt] = p1, T[3 * nt + 1] = (M - 1) * K + j, T[3 * nt + 2] = (M - 1) * K + (j + 1) % K, nt++;
    }
    const double V0 = vol(V, T, nt);
    if (V0 < 0) for (int t = 0; t < nt; t++) { const int x = T[3 * t + 1]; T[3 * t + 1] = T[3 * t + 2], T[3 * t + 2] = x; }
    const double vtot = fabs(V0);
    CHECK(oriented(T, nt, nv), "U tube: %d vertices, %d triangles, closed and wound one way, volume %.6e", nv, nt, vtot);
    mk_cut_out o;
    const double P[9] = {-0.05, 0, 0.05, 0, 0, 0.15, 0.05, 0, 0.03}, N[9] = {0, 0, 1, 1, 0, 0, 0, 0, -1}, R[3] = {0.02, 0.02, 0.02};
    {
        const int rc = mk_cut_loops(V, nv, T, nt, 3, P, N, R, &om, &o);
        CHECK(rc == 0 && o.nreg == 4 && o.nloop == 3 && o.code == 0 && o.bad == -1, "three cuts: rc %d (%s), regions %d, loops %d, cap triangles %d", rc, rc ? o.reason : "ok", o.nreg, o.nloop, o.ncapt);
        if (!rc) {
            double rel, dv;
            const int ok = pieces(&o, vtot, &rel, &dv, "cut_u");
            CHECK(ok && rel < 1e-12 && dv < 1e-12, "the 4 pieces closed, wound one way and outward, their volumes summing to the whole: %.2e relative (the input's %.2e)", rel, dv);
            CHECK(caps_shared(&o), "each cap in its two pieces, as it is on the - side and reversed on the + side");
            int near = 1;
            for (int l = 0; l < 3; l++) { /* the loop cut is the one round its cut's point: its crossings' centre within r of P */
                double c[3] = {0, 0, 0};
                int n = 0;
                for (int t = 0; t < o.ncapt; t++) if (o.capk[t] == l) for (int j = 0; j < 3; j++) { for (int q = 0; q < 3; q++) c[q] += o.V[3 * o.capT[3 * t + j] + q]; n++; }
                near &= n > 0 && hypot(hypot(c[0] / n - P[3 * l], c[1] / n - P[3 * l + 1]), c[2] / n - P[3 * l + 2]) < r;
            }
            CHECK(near && o.loopn[0] <= 2 * K + 2, "each loop the one round its point; the left cut's plane leaves the right leg whole (loop %d steps, a ring of %d)", o.loopn[0], K);
            int nnew;
            const int sm = srcmap(&o, V, T, nt, &nnew);
            CHECK(sm && nnew == o.nv - nv && nnew > 0, "the source map: each of %d new vertices V[a] + t (V[b] - V[a]) bit for bit, (a, b) an input edge, t in [0.001, 0.999]", nnew);
            const double off = capoff(&o, N);
            CHECK(off < 1e-15, "every cap vertex on its plane: %.1e at most", off);
            const double q = minshape(o.V, o.capT, o.ncapt);
            CHECK(q > 0.002, "the caps' worst triangle, best-shaped ears first: twice its area over its squared edges %.4f, over 0.002 (an equilateral 0.29; the first valid ear each time 0.0006)", q);
        }
        mk_cut_free(&o);
    }
    { /* snapping: a cut 2 um off ring 50 of the left leg: the ring's vertices go onto the plane, the loop through them */
        const double zr = 0.1 * 3.0 * 50 / (M - 1);
        const double Ps[3] = {-0.05, 0, zr + 2e-6}, Ns[3] = {0, 0, 1}, Rr[1] = {0.02};
        const int rs = mk_cut_loops(V, nv, T, nt, 1, Ps, Ns, Rr, &om, &o);
        double rel = 1, dv = 1, moved = 0;
        int ok = 0, za = 0;
        if (!rs) {
            ok = pieces(&o, vtot, &rel, &dv, NULL);
            for (int g = 0; g < o.nreg; g++) { int *TT, m; mk_cut_piece(&o, g, &TT, &m); za += minshape(o.V, TT, m) < 1e-6; free(TT); }
            for (int i = 0; i < o.nsnapped; i++) { const int v = o.snapped[i]; moved = fmax(moved, hypot(hypot(o.V[3 * v] - V[3 * v], o.V[3 * v + 1] - V[3 * v + 1]), o.V[3 * v + 2] - V[3 * v + 2])); }
        }
        CHECK(!rs && o.nreg == 2 && ok && za == 0 && rel < 1e-12 && dv < 1e-6 && o.loopn[0] == 2 * K && o.nv == nv && o.nsnapped == K && moved < om.snap,
              "a cut 2 um off a ring snaps onto it: rc %d, regions %d, closed %d, no sliver %d, volume %.1e relative (the input's %.1e), loop %d steps, %d new vertices, %d snapped, moved %.2g",
              rs, o.nreg, ok, za == 0, rel, dv, rs ? 0 : o.loopn[0], rs ? -1 : o.nv - nv, o.nsnapped, moved);
        const double off = rs ? 1 : capoff(&o, Ns);
        CHECK(!rs && off < 1e-15, "the snapped ring lies on the plane: %.1e at most", off);
        mk_cut_free(&o);
    }
    { /* the refusals on the U tube */
        const double Rs[3] = {0.005, 0.02, 0.02};
        int n = mk_cut_loops(V, nv, T, nt, 3, P, N, Rs, &om, &o);
        CHECK(n == MK_CUT_FAR && o.bad == 0 && !o.V && !o.T, "a search radius under the tube's: MK_CUT_FAR, out empty (%s)", o.reason);
        mk_cut_free(&o);
        const double P2[6] = {-0.05, 0, 0.05, -0.05, 0, 0.051}, N2[6] = {0, 0, 1, 0.6, 0, 0.8}, R2[2] = {0.02, 0.02};
        n = mk_cut_loops(V, nv, T, nt, 2, P2, N2, R2, &om, &o);
        CHECK(n == MK_CUT_CROSS && o.bad == 1, "two crossing cuts: MK_CUT_CROSS (%s)", o.reason);
        mk_cut_free(&o);
        const double P3[3] = {0.5, 0.5, 0.5};
        n = mk_cut_loops(V, nv, T, nt, 1, P3, N, R, &om, &o);
        CHECK(n == MK_CUT_NOLOOP, "a point off the mesh: MK_CUT_NOLOOP (%s)", o.reason);
        mk_cut_free(&o);
        n = mk_cut_clip(V, nv, T, nt, 1, P3, N, &om, &o);
        CHECK(n == MK_CUT_NOLOOP, "a clip plane that misses the mesh: MK_CUT_NOLOOP (%s)", o.reason);
        mk_cut_free(&o);
        const int t0 = T[0];
        T[0] = nv;
        n = mk_cut_loops(V, nv, T, nt, 1, P, N, R, &om, &o);
        T[0] = t0;
        const double v0 = V[1];
        V[1] = NAN;
        const int n2 = mk_cut_loops(V, nv, T, nt, 1, P, N, R, &om, &o);
        V[1] = v0;
        const double Z[3] = {0, 0, 0};
        const int n3 = mk_cut_loops(V, nv, T, nt, 1, P, Z, R, &om, &o);
        const int n4 = mk_cut_clip(V, nv, T, nt, 0, P, N, &om, &o);
        CHECK(n == -1 && n2 == -1 && n3 == -1 && n4 == -1, "bad arguments: an index out of range %d, a NaN coordinate %d, a zero normal %d, no cut %d (all -1)", n, n2, n3, n4);
        n = mk_cut_loops(V, nv, T, nt - 1, 1, P, N, R, &om, &o);
        CHECK(n == MK_CUT_OPEN, "a triangle dropped: MK_CUT_OPEN (%s)", o.reason);
        mk_cut_free(&o);
    }
    { /* a cone (radius 0.2 + 0.5 z): at P its loop reaches 0.002 past rmax, one step down it passes */
        tm_mesh c = {0};
        tube(&c, 40, K, 0.2, 0.5, 1);
        outward(&c);
        const double Pc[3] = {0, 0, 0.5035}, Nc[3] = {0, 0, 1}, Rc[1] = {0.2 + 0.5 * 0.5035 - 0.002};
        const int rc = mk_cut_loops(c.V, c.nv, c.T, c.nt, 1, Pc, Nc, Rc, &o1, &o);
        const double want = Pc[2] + (-1 * 1) * o1.step * 1.0;
        CHECK(rc == 0 && o.Pc[2] == want && o.Pc[0] == 0 && o.Pc[1] == 0 && o.loopr[0] <= Rc[0], "a cone: the plane moved one step down to pass rmax: rc %d, Pc z %.6f (want %.6f), loop %.5f <= %.5f", rc, rc ? 0 : o.Pc[2], want, rc ? 0 : o.loopr[0], Rc[0]);
        mk_cut_free(&o);
        tm_free(&c);
    }
    { /* the same cone, shift 0.6 of a step: no whole step fits in shift, so the plane may not move and the cut is refused */
        tm_mesh c = {0};
        tube(&c, 40, K, 0.2, 0.5, 1);
        outward(&c);
        const double Pc[3] = {0, 0, 0.5035}, Nc[3] = {0, 0, 1}, Rc[1] = {0.2 + 0.5 * 0.5035 - 0.002};
        mk_cut_opts os = o1;
        os.shift = 0.6 * os.step;
        const int rc = mk_cut_loops(c.V, c.nv, c.T, c.nt, 1, Pc, Nc, Rc, &os, &o);
        CHECK(rc == MK_CUT_FAR, "a cone, shift 0.6 of a step: refused, not moved a step past shift: rc %d (want %d), Pc z %.6f", rc, MK_CUT_FAR, rc ? Pc[2] : o.Pc[2]);
        mk_cut_free(&o);
        tm_free(&c);
    }
    { /* ear clipping of a triangle with three points along each side (collinear to rounding): it used to find no ear */
        const double C[3][2] = {{0.62771525095162328, 0.13010466858942163}, {-0.09252613225740336, 0.097665699433785863}, {-0.032775963062944724, -1.2495572287837899}};
        double q[12][2];
        int tri[10][3], n = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 4; j++) q[n][0] = C[i][0] + j / 4.0 * (C[(i + 1) % 3][0] - C[i][0]), q[n][1] = C[i][1] + j / 4.0 * (C[(i + 1) % 3][1] - C[i][1]), n++;
        const int rc = mk_cut_ears_(12, (const double (*)[2])q, tri);
        double a = 0, worst = INFINITY;
        int used[12] = {0}, all = 1;
        for (int t = 0; !rc && t < 10; t++) {
            const double *x = q[tri[t][0]], *y = q[tri[t][1]], *z = q[tri[t][2]];
            const double tw = (y[0] - x[0]) * (z[1] - x[1]) - (y[1] - x[1]) * (z[0] - x[0]);
            a += tw / 2, worst = fmin(worst, tw);
            for (int j = 0; j < 3; j++) used[tri[t][j]] = 1;
        }
        for (int i = 0; i < 12; i++) all &= used[i];
        const double want = ((C[1][0] - C[0][0]) * (C[2][1] - C[0][1]) - (C[1][1] - C[0][1]) * (C[2][0] - C[0][0])) / 2;
        CHECK(!rc && fabs(a - want) < 1e-14 && worst > 0 && all, "ear clipping, a triangle with three points a side: rc %d (0: clipped), area %.15g of %.15g, every corner used %d, no triangle turned %d", rc, a, want, all, worst > 0);
    }
    { /* a ball clipped: one plane, two halves; then three planes, four pieces */
        tm_mesh b = {0};
        const double c0[3] = {0, 0, 0};
        tm_sphere(&b, 1, c0, 4);
        outward(&b);
        const double vb = tm_volume(b.V, b.T, b.nt);
        const double Pb[9] = {0, 0, 0.1234567, 0, 0, -0.5, 0, 0, 0.6}, Nb[9] = {0, 0, 1, 0, 0, 1, 0, 0, 1};
        int rc = mk_cut_clip(b.V, b.nv, b.T, b.nt, 1, Pb, Nb, &o1, &o);
        double rel = 1, dv = 1;
        int ok = !rc && pieces(&o, vb, &rel, &dv, "cut_ball");
        int plus = 1;
        if (!rc) for (int t = 0; t < o.nt; t++) for (int j = 0; j < 3; j++) { const int v = o.T[3 * t + j]; if (v < b.nv && o.V[3 * v + 2] > Pb[2] + o1.snap) plus &= o.reg[t] == o.side[1]; }
        CHECK(ok && o.nloop == 1 && o.nreg == 2 && rel < 1e-12 && plus && caps_shared(&o), "a ball clipped by one plane: rc %d, loops %d, pieces %d, closed and outward, volumes %.1e relative, the + piece above it", rc, o.nloop, o.nreg, rel);
        CHECK(!rc && dv < o1.snap * 4 * T_PI / vb && o.nsnapped > 0, "its %d snapped vertices (snap %g) moved the volume by %.1e of it, under snap times the area over the volume", o.nsnapped, o1.snap, dv);
        mk_cut_free(&o);
        rc = mk_cut_clip(b.V, b.nv, b.T, b.nt, 3, Pb, Nb, &o1, &o);
        ok = !rc && pieces(&o, vb, &rel, &dv, NULL);
        CHECK(ok && o.nloop == 3 && o.nreg == 4 && rel < 1e-12 && caps_shared(&o), "a ball clipped by three planes in one call: rc %d, loops %d, pieces %d, volumes %.1e relative", rc, o.nloop, o.nreg, rel);
        mk_cut_free(&o);
        /* the ball wound inward: the pieces follow it */
        for (int t = 0; t < b.nt; t++) { const int x = b.T[3 * t + 1]; b.T[3 * t + 1] = b.T[3 * t + 2], b.T[3 * t + 2] = x; }
        rc = mk_cut_clip(b.V, b.nv, b.T, b.nt, 1, Pb, Nb, &o1, &o);
        int inw = !rc && o.nreg == 2;
        double vs = 0;
        for (int g = 0; inw && g < o.nreg; g++) {
            int *TT, m;
            mk_cut_piece(&o, g, &TT, &m);
            const double v = vol(o.V, TT, m);
            inw &= oriented(TT, m, o.nv) && v < 0;
            vs += v;
            free(TT);
        }
        const double vsplit = rc ? 1 : vol(o.V, o.T, o.nt);
        CHECK(inw && fabs(vs - vsplit) < 1e-12 * fabs(vsplit) && caps_shared(&o), "the ball wound inward, clipped: rc %d, pieces %d, each wound one way and inward, volumes %.1e relative", rc, o.nreg, fabs(vs - vsplit) / fabs(vsplit));
        mk_cut_free(&o);
        tm_free(&b);
    }
    { /* a torus (R 1, r 0.3) */
        tm_mesh t = {0};
        tm_torus(&t, 1, 0.3, 48, 24);
        outward(&t);
        const double vt = tm_volume(t.V, t.T, t.nt);
        const double Pt[3] = {0.0123, 0, 0}, Nt[3] = {1, 0, 0};
        int rc = mk_cut_clip(t.V, t.nv, t.T, t.nt, 1, Pt, Nt, &o1, &o);
        double rel = 1, dv = 1;
        int ok = !rc && pieces(&o, vt, &rel, &dv, "cut_torus");
        CHECK(ok && o.nloop == 2 && o.nreg == 2 && rel < 1e-12 && dv < 1e-12 && caps_shared(&o), "a torus clipped through its axis: rc %d, loops %d, pieces %d (two halves, each with two caps), volumes %.1e relative", rc, o.nloop, o.nreg, rel);
        int nnew;
        CHECK(!rc && srcmap(&o, t.V, t.T, t.nt, &nnew) && nnew > 0, "its source map: %d new vertices bit for bit", rc ? 0 : nnew);
        if (!rc && tm_dir) { /* the negative control for meshcheck: a half with its + side cap not reversed */
            int *TT, m;
            const int g = o.side[1];
            mk_cut_piece(&o, g, &TT, &m);
            unsigned char *used = calloc((size_t)o.nv, 1);
            for (int i = 0; i < 3 * m; i++) used[TT[i]] = 1;
            int n0 = 0;
            for (int q = 0; q < o.nt; q++) n0 += o.reg[q] == g;
            for (int q = n0; q < m; q++) { const int x = TT[3 * q + 1]; TT[3 * q + 1] = TT[3 * q + 2], TT[3 * q + 2] = x; }
            tm_write("cut_flip", o.V, o.nv, TT, m, used);
            free(used), free(TT);
        }
        mk_cut_free(&o);
        const double Pl[3] = {0.0123, 1, 0}, Rl[1] = {0.5};
        rc = mk_cut_loops(t.V, t.nv, t.T, t.nt, 1, Pl, Nt, Rl, &o1, &o);
        CHECK(rc == MK_CUT_REGIONS, "one loop of the two through its axis does not separate: MK_CUT_REGIONS (%s)", o.reason);
        mk_cut_free(&o);
        const double Pe[3] = {0, 0, 0.0123}, Ne[3] = {0, 0, 1};
        rc = mk_cut_clip(t.V, t.nv, t.T, t.nt, 1, Pe, Ne, &o1, &o);
        CHECK(rc == MK_CUT_NESTED && o.bad == 0, "a torus clipped through its equator (its inner loop inside its outer): MK_CUT_NESTED (%s)", o.reason);
        mk_cut_free(&o);
        tm_free(&t);
    }
    { /* the snap where moving onto the plane would change the cut: snap 0.001 */
        mk_cut_opts os = o1;
        os.snap = 0.001;
        tm_mesh t = {0};
        tm_torus(&t, 1, 0.3, 48, 24);
        outward(&t);
        /* a plane 0.0005 off the torus's inner equator vertex (0.7, 0, 0), a saddle: on the plane, the loop's two passes
           beside it would meet at it */
        int sv = -1;
        for (int i = 0; i < t.nv; i++) if (fabs(t.V[3 * i] - 0.7) < 1e-12 && fabs(t.V[3 * i + 1]) < 1e-12 && fabs(t.V[3 * i + 2]) < 1e-12) sv = i;
        const double Ps[3] = {0.7005, 0, 0}, Ns[3] = {1, 0, 0};
        int rc = mk_cut_clip(t.V, t.nv, t.T, t.nt, 1, Ps, Ns, &os, &o);
        double rel = 1, dv = 1;
        int ok = !rc && pieces(&o, tm_volume(t.V, t.T, t.nt), &rel, &dv, NULL), listed = 0;
        for (int i = 0; !rc && i < o.nsnapped; i++) listed |= o.snapped[i] == sv;
        CHECK(sv >= 0 && ok && o.nloop == 1 && o.nreg == 2 && rel < 1e-12 && !listed && caps_shared(&o),
              "a plane %.4f off a saddle vertex (snap %g): it keeps its place, the loop passes it on both sides: rc %d, loops %d, pieces %d, closed, volumes %.1e (%s)", 0.0005, os.snap, rc, o.nloop, o.nreg, rel, rc ? o.reason : "ok");
        mk_cut_free(&o);
        tm_free(&t);
        /* a plane 0.0005 inside a ball's vertex, every neighbour on the far side: a dimple thinner than snap, not cut */
        tm_mesh b = {0};
        const double c0[3] = {0, 0, 0};
        tm_sphere(&b, 1, c0, 3);
        outward(&b);
        const double *v0 = b.V, Pd[3] = {0.9995 * v0[0], 0.9995 * v0[1], 0.9995 * v0[2]}, Nd[3] = {-v0[0], -v0[1], -v0[2]};
        rc = mk_cut_clip(b.V, b.nv, b.T, b.nt, 1, Pd, Nd, &os, &o);
        CHECK(rc == MK_CUT_NOLOOP, "a plane 0.0005 short of a ball's vertex, all its neighbours beyond it: no loop, MK_CUT_NOLOOP (%s)", o.reason);
        mk_cut_free(&o);
        tm_free(&b);
    }
    { /* a small closed tube inside a big one (two pieces): the probe */
        static double V2[3 * 2 * (K * 20 + 2)];
        static int T2[3 * 2 * (2 * K * 20 + 2 * K)];
        int n2 = 0, t2 = 0;
        for (int c = 0; c < 2; c++) {
            const double rr = c ? 0.005 : 0.03, x0 = c ? 0.015 : 0, z0 = c ? 0.02 : 0, z1 = c ? 0.08 : 0.1;
            const int b = n2;
            for (int i = 0; i < 20; i++) for (int j = 0; j < K; j++) { const double a = 2 * T_PI * j / K; V2[3 * n2] = x0 + rr * cos(a), V2[3 * n2 + 1] = rr * sin(a), V2[3 * n2 + 2] = z0 + (z1 - z0) * i / 19; n2++; }
            const int pa = n2, pb = n2 + 1;
            V2[3 * n2] = x0, V2[3 * n2 + 1] = 0, V2[3 * n2 + 2] = z0 - rr; n2++;
            V2[3 * n2] = x0, V2[3 * n2 + 1] = 0, V2[3 * n2 + 2] = z1 + rr; n2++;
            for (int i = 0; i + 1 < 20; i++) for (int j = 0; j < K; j++) {
                const int a = b + i * K + j, bb = b + i * K + (j + 1) % K, cc = b + (i + 1) * K + (j + 1) % K, d = b + (i + 1) * K + j;
                T2[3 * t2] = a, T2[3 * t2 + 1] = bb, T2[3 * t2 + 2] = cc, t2++; T2[3 * t2] = a, T2[3 * t2 + 1] = cc, T2[3 * t2 + 2] = d, t2++;
            }
            for (int j = 0; j < K; j++) { T2[3 * t2] = pa, T2[3 * t2 + 1] = b + (j + 1) % K, T2[3 * t2 + 2] = b + j, t2++; T2[3 * t2] = pb, T2[3 * t2 + 1] = b + 19 * K + j, T2[3 * t2 + 2] = b + 19 * K + (j + 1) % K, t2++; }
        }
        mk_cut_probe pr;
        const int ri = mk_cut_probe_init(&pr, V2, n2, T2, t2, &om);
        double rl, ar;
        int m, nest;
        const double pz[3] = {-0.015, 0, 0.05}, nz[3] = {0, 0, 1}, pin[3] = {0.015, 0, 0.05}, far[3] = {-0.02, 0.01, 0.05};
        const int r1 = mk_cut_probe_loop(&pr, pz, nz, &rl, &ar, &m, &nest);
        const int in1 = mk_cut_probe_inside(&pr, pin), in2 = mk_cut_probe_inside(&pr, far);
        CHECK(!ri && !r1 && nest == 1 && in1 && in2, "the probe: a loop round a point between the tubes holds the small tube's loop: rc %d, nested %d (loop %.1f mm), both points inside it", r1, nest, rl * 1e3);
        const int r2 = mk_cut_probe_loop(&pr, pin, nz, &rl, &ar, &m, &nest);
        const int in3 = mk_cut_probe_inside(&pr, pin), in4 = mk_cut_probe_inside(&pr, far);
        CHECK(!r2 && nest == 0 && rl < 0.006 && in3 && !in4, "a loop round the small tube's own axis: rc %d, nested %d (loop %.1f mm), its centre inside, the far point not", r2, nest, rl * 1e3);
        const int r3 = mk_cut_probe_loop(&pr, P, nz, &rl, &ar, &m, &nest);
        CHECK(r3 == 1, "a probe plane round a point outside every loop: 1");
        const double Pn[3] = {-0.015, 0, 0.05}, Rn[1] = {0.05};
        const int r4 = mk_cut_loops(V2, n2, T2, t2, 1, Pn, nz, Rn, &om, &o);
        CHECK(r4 == MK_CUT_NESTED && o.bad == 0, "mk_cut_loops there: MK_CUT_NESTED (%s)", o.reason);
        mk_cut_free(&o);
        mk_cut_probe_free(&pr);
        int *TT, mm;
        CHECK(mk_cut_piece(&o, 0, &TT, &mm) == -1 && !TT, "mk_cut_piece of an empty result: -1");
    }
    { /* a closed tube whose cross-section is a figure eight (the surface crosses itself): its loop is not simple */
        tm_mesh f = {0};
        for (int i = 0; i < 20; i++)
            for (int j = 0; j < K; j++) { const double a = 2 * T_PI * (j + 0.5) / K; tm_vert(&f, cos(a), sin(a) * cos(a), i / 19.0); }
        const int q0 = tm_vert(&f, 0.5, 0.2, -0.3), q1 = tm_vert(&f, 0.5, 0.2, 1.3);
        for (int i = 0; i + 1 < 20; i++)
            for (int j = 0; j < K; j++) {
                const int a = i * K + j, b = i * K + (j + 1) % K, c = (i + 1) * K + (j + 1) % K, d = (i + 1) * K + j;
                tm_tri(&f, a, b, c), tm_tri(&f, a, c, d);
            }
        for (int j = 0; j < K; j++) tm_tri(&f, q0, (j + 1) % K, j), tm_tri(&f, q1, 19 * K + j, 19 * K + (j + 1) % K);
        const double Pf[3] = {0.5, 0.2, 0.5123}, Nf[3] = {0, 0, 1};
        const int rc = mk_cut_clip(f.V, f.nv, f.T, f.nt, 1, Pf, Nf, &o1, &o);
        CHECK(rc == MK_CUT_CAP && o.bad == 0, "a tube of figure-eight section (crossing itself) clipped: MK_CUT_CAP (%s)", o.reason);
        mk_cut_free(&o);
        tm_free(&f);
    }
    { /* never a piece left open: a snap large against the edges near the cut is refused (MK_CUT_SPLIT), not returned open */
        tm_mesh c = {0};
        tube(&c, 20, 24, 0.2, 0.3, 1), outward(&c);
        mk_cut_opts ob;
        mk_cut_defaults(&ob, 1);
        ob.snap = 0.055, ob.shift = 0.01, ob.step = 0.001;
        const double P[3] = {0, 0, 0.262}, N[3] = {0.3, 0.1, 1};
        mk_cut_out o;
        const int rc = mk_cut_clip(c.V, c.nv, c.T, c.nt, 1, P, N, &ob, &o);
        CHECK(rc == MK_CUT_SPLIT && o.nreg == 0, "a snap of 0.055 against edges of about 0.05: refused MK_CUT_SPLIT (%d), not returned open: %s", rc, o.reason);
        if (rc == 0) mk_cut_free(&o);
        /* one triangle wound the other way round: not consistently wound, MK_CUT_OPEN */
        const int x = c.T[3 * 7 + 1];
        c.T[3 * 7 + 1] = c.T[3 * 7 + 2], c.T[3 * 7 + 2] = x;
        const double P2[3] = {0, 0, 0.5}, N2[3] = {0, 0, 1};
        const int r2 = mk_cut_clip(c.V, c.nv, c.T, c.nt, 1, P2, N2, &o1, &o);
        CHECK(r2 == MK_CUT_OPEN, "one triangle flipped: refused MK_CUT_OPEN (%d): %s", r2, o.reason);
        if (r2 == 0) mk_cut_free(&o);
        tm_free(&c);
    }
    return wrong;
}
