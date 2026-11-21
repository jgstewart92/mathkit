/* mk_geom against independent references written here: closest point on a triangle (projection and the three edges),
   ray and triangle (Moller and Trumbore 1997), BVH queries against brute force over every triangle, inside or out
   against the analytic sphere and torus, crossing pairs against a signed-volume segment test; then edge cases, one row
   each, against exact values or brute force. */
#include "../include/mk_geom.h"
#include "t.h"
#include "tmesh.h"

/* segment pq strictly through triangle abc by signed volumes (p and q strictly either side of its plane, pq strictly
   inside the three edge wedges), not Moller-Trumbore */
static double orient_ref(const double *a, const double *b, const double *c, const double *d)
{
    const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]}, w[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
    return u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
}
[[maybe_unused]] static int segtri_ref(const double *p, const double *q, const double *a, const double *b, const double *c)
{
    const double s1 = orient_ref(a, b, c, p), s2 = orient_ref(a, b, c, q);
    if (!((s1 > 0 && s2 < 0) || (s1 < 0 && s2 > 0))) return 0;
    const double e0 = orient_ref(p, q, a, b), e1 = orient_ref(p, q, b, c), e2 = orient_ref(p, q, c, a);
    return (e0 > 0 && e1 > 0 && e2 > 0) || (e0 < 0 && e1 < 0 && e2 < 0);
}
static double seg_d2(const double *p, const double *a, const double *b)
{
    double ab[3], ap[3], t = 0, l = 0, d = 0;
    for (int k = 0; k < 3; k++) ab[k] = b[k] - a[k], ap[k] = p[k] - a[k], t += ab[k] * ap[k], l += ab[k] * ab[k];
    t = l > 0 ? fmin(1, fmax(0, t / l)) : 0;
    for (int k = 0; k < 3; k++) { double e = a[k] + t * ab[k] - p[k]; d += e * e; }
    return d;
}
static double brute_d2(const double *p, const double *a, const double *b, const double *c)
{
    double e1[3], e2[3], n[3], ap[3];
    for (int k = 0; k < 3; k++) e1[k] = b[k] - a[k], e2[k] = c[k] - a[k], ap[k] = p[k] - a[k];
    n[0] = e1[1] * e2[2] - e1[2] * e2[1], n[1] = e1[2] * e2[0] - e1[0] * e2[2], n[2] = e1[0] * e2[1] - e1[1] * e2[0];
    double nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2], best = fmin(seg_d2(p, a, b), fmin(seg_d2(p, b, c), seg_d2(p, c, a)));
    if (nn > 0) {
        double h = (ap[0] * n[0] + ap[1] * n[1] + ap[2] * n[2]) / nn, q[3];
        for (int k = 0; k < 3; k++) q[k] = p[k] - h * n[k];
        double s = 1, prev = 0;
        const double *v[3] = {a, b, c};
        for (int i = 0; i < 3; i++) { /* is q inside: the same side of all three edges */
            const double *u = v[i], *w = v[(i + 1) % 3];
            double eu[3], qu[3], cr[3];
            for (int k = 0; k < 3; k++) eu[k] = w[k] - u[k], qu[k] = q[k] - u[k];
            cr[0] = eu[1] * qu[2] - eu[2] * qu[1], cr[1] = eu[2] * qu[0] - eu[0] * qu[2], cr[2] = eu[0] * qu[1] - eu[1] * qu[0];
            double sg = cr[0] * n[0] + cr[1] * n[1] + cr[2] * n[2];
            if (i && sg * prev < 0) s = 0;
            prev = sg;
        }
        if (s) best = fmin(best, h * h * nn);
    }
    return best;
}
/* ray o + t d against triangle abc, t > 0 (Moller and Trumbore 1997): the reference for mk_ray_tri */
static int mt(const double *o, const double *d, const double *a, const double *b, const double *c, double *t)
{
    double e1[3], e2[3], pv[3], tv[3], qv[3];
    for (int k = 0; k < 3; k++) e1[k] = b[k] - a[k], e2[k] = c[k] - a[k], tv[k] = o[k] - a[k];
    pv[0] = d[1] * e2[2] - d[2] * e2[1], pv[1] = d[2] * e2[0] - d[0] * e2[2], pv[2] = d[0] * e2[1] - d[1] * e2[0];
    double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
    if (fabs(det) < 1e-300) return 0;
    double u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) / det;
    qv[0] = tv[1] * e1[2] - tv[2] * e1[1], qv[1] = tv[2] * e1[0] - tv[0] * e1[2], qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
    double v = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) / det;
    *t = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) / det;
    return u >= 0 && v >= 0 && u + v <= 1 && *t > 0;
}

static double V[60000];
static int T[120000], nv, nt;
static void torus(int nu, int nw, double R, double r)
{
    nv = nt = 0;
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nw; j++) {
            double u = 2 * T_PI * i / nu, w = 2 * T_PI * j / nw;
            V[3 * nv] = (R + r * cos(w)) * cos(u), V[3 * nv + 1] = (R + r * cos(w)) * sin(u), V[3 * nv + 2] = r * sin(w), nv++;
        }
    for (int i = 0; i < nu; i++)
        for (int j = 0; j < nw; j++) {
            int a = i * nw + j, b = ((i + 1) % nu) * nw + j, c = ((i + 1) % nu) * nw + (j + 1) % nw, d = i * nw + (j + 1) % nw;
            T[3 * nt] = a, T[3 * nt + 1] = b, T[3 * nt + 2] = c, nt++;
            T[3 * nt] = a, T[3 * nt + 1] = c, T[3 * nt + 2] = d, nt++;
        }
}
static void sphere(int nu, int nw) /* a UV sphere of radius 1 with poles, outward */
{
    nv = nt = 0;
    V[0] = 0, V[1] = 0, V[2] = 1, nv = 1;
    for (int i = 1; i < nw; i++)
        for (int j = 0; j < nu; j++) {
            double th = T_PI * i / nw, ph = 2 * T_PI * j / nu;
            V[3 * nv] = sin(th) * cos(ph), V[3 * nv + 1] = sin(th) * sin(ph), V[3 * nv + 2] = cos(th), nv++;
        }
    V[3 * nv] = 0, V[3 * nv + 1] = 0, V[3 * nv + 2] = -1;
    int south = nv++;
    for (int j = 0; j < nu; j++) T[3 * nt] = 0, T[3 * nt + 1] = 1 + j, T[3 * nt + 2] = 1 + (j + 1) % nu, nt++;
    for (int i = 0; i < nw - 2; i++)
        for (int j = 0; j < nu; j++) {
            int a = 1 + i * nu + j, b = 1 + (i + 1) * nu + j, c = 1 + (i + 1) * nu + (j + 1) % nu, d = 1 + i * nu + (j + 1) % nu;
            T[3 * nt] = a, T[3 * nt + 1] = b, T[3 * nt + 2] = c, nt++;
            T[3 * nt] = a, T[3 * nt + 1] = c, T[3 * nt + 2] = d, nt++;
        }
    for (int j = 0; j < nu; j++) {
        int a = 1 + (nw - 2) * nu + j, b = 1 + (nw - 2) * nu + (j + 1) % nu;
        T[3 * nt] = a, T[3 * nt + 1] = south, T[3 * nt + 2] = b, nt++;
    }
}

int main(void)
{
    int cbad = 0;
    double cw = 0;
    for (int k = 0; k < 200000; k++) {
        double p[3], a[3], b[3], c[3], w[3], q[3];
        for (int i = 0; i < 3; i++) p[i] = 2 * t_rand(), a[i] = t_rand(), b[i] = t_rand(), c[i] = k % 7 ? t_rand() : a[i] + 1e-3 * t_rand();
        double d1 = mk_closest_pt_tri(p, a, b, c, w, q), d2 = brute_d2(p, a, b, c);
        cw = fmax(cw, fabs(d1 - d2));
        if (fabs(d1 - d2) > 1e-12 * (1 + d2)) cbad++;
    }
    CHECK(cbad == 0, "closest point on a triangle against an independent projection-and-edges method, 200000 cases with slivers (%d differ, worst %.2g)", cbad, cw);
    int rbad = 0, rcount = 0;
    for (int k = 0; k < 200000; k++) {
        double o[3], d[3], a[3], b[3], c[3], t1 = 0, t2 = 0, br[3];
        for (int i = 0; i < 3; i++) o[i] = 2 * t_rand(), d[i] = t_rand(), a[i] = t_rand(), b[i] = t_rand(), c[i] = t_rand();
        mk_ray r = mk_ray_make(o, d);
        int h1 = mk_ray_tri(&r, a, b, c, 0, INFINITY, &t1, br), h2 = mt(o, d, a, b, c, &t2);
        if (h1 && fmin(br[0], fmin(br[1], br[2])) < 1e-9) continue; /* too near an edge for two methods to agree */
        rcount += h1;
        if (h1 != h2 || (h1 && fabs(t1 - t2) > 1e-9 * (1 + t2))) rbad++;
    }
    CHECK(rbad == 0 && rcount > 1000, "watertight ray-triangle against Moller-Trumbore, 200000 random rays (%d hits, %d differ)", rcount, rbad);

    sphere(64, 32);
    int miss = 0, mtmiss = 0;
    for (int t = 0; t < nt; t++)
        for (int e = 0; e < 4; e++) { /* aim from the centre at each vertex and each edge midpoint */
            double d[3], o[3] = {0, 0, 0}, tt;
            const double *a = V + 3 * T[3 * t + e % 3], *b = V + 3 * T[3 * t + (e + 1) % 3];
            for (int k = 0; k < 3; k++) d[k] = e == 3 ? a[k] : 0.5 * (a[k] + b[k]);
            int hit = 0, hm = 0;
            mk_ray r = mk_ray_make(o, d);
            for (int s = 0; s < nt; s++) {
                hit += mk_ray_tri(&r, V + 3 * T[3 * s], V + 3 * T[3 * s + 1], V + 3 * T[3 * s + 2], 0, INFINITY, &tt, NULL);
                hm += mt(o, d, V + 3 * T[3 * s], V + 3 * T[3 * s + 1], V + 3 * T[3 * s + 2], &tt);
            }
            miss += hit == 0, mtmiss += hm == 0;
        }
    CHECK(miss == 0, "rays from inside aimed exactly at every edge and vertex of a 3968-triangle sphere: %d slip through (Moller-Trumbore: %d)", miss, mtmiss);

    int fails = 0, nhits = 0, nearbad = 0, anybad = 0, inbad = 0, wbad = 0;
    for (int mesh = 0; mesh < 2; mesh++) {
        if (mesh == 0) sphere(48, 24);
        else torus(64, 32, 2, 0.7);
        mk_bvh B;
        mk_bvh_build(&B, V, T, nt);
        for (int k = 0; k < 4000; k++) {
            double o[3], d[3], bt = INFINITY, tb;
            int bi = -1;
            for (int i = 0; i < 3; i++) o[i] = 3 * t_rand(), d[i] = t_rand();
            mk_ray r = mk_ray_make(o, d);
            for (int s = 0; s < nt; s++) {
                double th;
                if (mk_ray_tri(&r, V + 3 * T[3 * s], V + 3 * T[3 * s + 1], V + 3 * T[3 * s + 2], 0, INFINITY, &th, NULL) && th < bt) bt = th, bi = s;
            }
            int hi = mk_bvh_first_hit(&B, o, d, 0, INFINITY, &tb, NULL);
            if (hi != bi || (bi >= 0 && tb != bt)) fails++;
            nhits += bi >= 0;
            if (mk_bvh_any_hit(&B, o, d, 0, 1.5) != (bi >= 0 && bt < 1.5)) anybad++;
            double p[3] = {3 * t_rand(), 3 * t_rand(), 3 * t_rand()}, q[3], bd = INFINITY;
            int bn = -1, tn;
            for (int s = 0; s < nt; s++) {
                double w[3], c[3], d2 = mk_closest_pt_tri(p, V + 3 * T[3 * s], V + 3 * T[3 * s + 1], V + 3 * T[3 * s + 2], w, c);
                if (d2 < bd) bd = d2, bn = s;
            }
            double nd = mk_bvh_nearest(&B, p, INFINITY, &tn, q);
            if (tn != bn || nd != bd) nearbad++;
            /* inside: parity against the analytic solid, away from the facets, and against the winding number */
            double rr = mesh == 0 ? sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) : 0;
            double tq = mesh == 1 ? hypot(hypot(p[0], p[1]) - 2, p[2]) : 0;
            int truth = mesh == 0 ? (rr < 0.99 ? 1 : rr > 1.0 ? 0 : -1) : (tq < 0.69 ? 1 : tq > 0.7 ? 0 : -1);
            if (truth >= 0 && k < 1500) {
                int in = mk_bvh_inside(&B, p);
                double wn = mk_winding_number(V, T, nt, p);
                if (in != truth) inbad++;
                if ((wn > 0.5) != truth) wbad++;
            }
        }
        mk_bvh_free(&B);
    }
    CHECK(fails == 0 && nhits > 1000, "BVH first hit equals brute force on a sphere and a torus: 8000 rays, %d hits, %d differ", nhits, fails);
    CHECK(anybad == 0, "BVH any-hit within t 1.5 agrees with brute force (%d differ)", anybad);
    CHECK(nearbad == 0, "BVH nearest point equals brute force: 8000 points, same triangle and distance (%d differ)", nearbad);
    CHECK(inbad == 0 && wbad == 0, "inside or out: ray parity and the winding number both match the analytic sphere and torus (%d, %d wrong of 3000)", inbad, wbad);
    mk_bvh E;
    double tt, qq[3] = {0, 0, 0};
    int tn;
    mk_bvh_build(&E, V, T, 0);
    CHECK(mk_bvh_first_hit(&E, (double[3]){0, 0, 0}, (double[3]){1, 0, 0}, 0, INFINITY, &tt, NULL) == -1 && mk_bvh_nearest(&E, qq, 5, &tn, qq) == 5 && tn == -1, "an empty tree answers no hit and nothing near");
    mk_bvh_free(&E);

    /* edge cases, one row each: interval ends, non-finite and huge coordinates, deep trees, ties, slivers */
    {
        static const double V1[] = {0, 0, 1, 1, 0, 1, 0, 1, 1};
        static const int T1[] = {0, 1, 2};
        mk_bvh B;
        mk_bvh_build(&B, V1, T1, 1);
        double o[3] = {0.2, 0.2, 0}, d[3] = {0, 0, 1}, t = -1;
        int h = mk_bvh_first_hit(&B, o, d, 0, 1, &t, NULL), a = mk_bvh_any_hit(&B, o, d, 0, 1);
        CHECK(h == -1 && a == 0, "a hit at exactly t = tmax is outside the open interval (first hit %d, any hit %d)", h, a);
        mk_bvh_free(&B);
    }
    {
        int ok = 0;
        for (int mode = 0; mode < 4; mode++) {
            double W[3 * 24];
            int TT[3 * 8];
            for (int t = 0; t < 8; t++)
                for (int j = 0; j < 3; j++) {
                    int v = 3 * t + j;
                    W[3 * v] = t + (j == 1), W[3 * v + 1] = (j == 2), W[3 * v + 2] = 0, TT[3 * t + j] = v;
                }
            if (mode == 0) W[0] = NAN;
            if (mode == 1) W[0] = INFINITY;
            if (mode == 2) W[0] = W[3] = W[6] = 1e308;
            if (mode == 3) for (int v = 21; v < 24; v++) W[3 * v] = 1e308;
            mk_bvh B;
            const int r = mk_bvh_build(&B, W, TT, 8);
            double o[3] = {3.2, 0.2, -1}, d[3] = {0, 0, 1}, th;
            ok += (r == 0 || r == -1) && (r || mk_bvh_first_hit(&B, o, d, 0, INFINITY, &th, NULL) == 3);
            if (!r) mk_bvh_free(&B);
        }
        CHECK(ok == 4, "a tree over a NaN, an infinite or a 1e308 vertex builds (or refuses) and still finds the good triangles (%d of 4)", ok);
    }
    {
        int nd = 0;
        static double W[9 * 6006];
        static int TT[3 * 6006];
        for (int ax = 0; ax < 6; ax++)
            for (int e = -500; e <= 500; e++) {
                double x = ldexp(ax & 1 ? -1 : 1, e), s = fabs(x) * 1e-3, P[9] = {0};
                int a = ax / 2, b = (a + 1) % 3;
                P[a] = x, P[3 + a] = x + s, P[6 + a] = x, P[6 + b] = s;
                for (int k = 0; k < 9; k++) W[9 * nd + k] = P[k];
                TT[3 * nd] = 3 * nd, TT[3 * nd + 1] = 3 * nd + 1, TT[3 * nd + 2] = 3 * nd + 2, nd++;
            }
        mk_bvh B;
        mk_bvh_build(&B, W, TT, nd);
        int fh = 0, ch = 0, nn = 0;
        for (int t = 0; t < nd; t += 3) {
            double c[3], e1[3], e2[3], n[3], o[3], d[3];
            for (int k = 0; k < 3; k++) c[k] = (W[9 * t + k] + W[9 * t + 3 + k] + W[9 * t + 6 + k]) / 3, e1[k] = W[9 * t + 3 + k] - W[9 * t + k], e2[k] = W[9 * t + 6 + k] - W[9 * t + k];
            n[0] = e1[1] * e2[2] - e1[2] * e2[1], n[1] = e1[2] * e2[0] - e1[0] * e2[2], n[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]), sc = fabs(c[0]) + fabs(c[1]) + fabs(c[2]);
            for (int k = 0; k < 3; k++) o[k] = c[k] + sc * n[k] / nl, d[k] = -n[k] / nl;
            mk_ray r = mk_ray_make(o, d);
            double bt = INFINITY, th, q[3], bd = INFINITY;
            int bi = -1, bc = 0, bn = -1, tn;
            for (int s = 0; s < nd; s++) {
                double tt, w[3], qq[3];
                if (mk_ray_tri(&r, W + 3 * TT[3 * s], W + 3 * TT[3 * s + 1], W + 3 * TT[3 * s + 2], 0, INFINITY, &tt, NULL)) { bc++; if (tt < bt) bt = tt, bi = s; }
                double d2 = mk_closest_pt_tri(c, W + 3 * TT[3 * s], W + 3 * TT[3 * s + 1], W + 3 * TT[3 * s + 2], w, qq);
                if (d2 < bd) bd = d2, bn = s;
            }
            fh += mk_bvh_first_hit(&B, o, d, 0, INFINITY, &th, NULL) != bi;
            ch += mk_bvh_count_hits(&B, o, d) != bc;
            double ndd = mk_bvh_nearest(&B, c, INFINITY, &tn, q);
            nn += tn != bn || ndd != bd;
        }
        CHECK(fh == 0 && ch == 0 && nn == 0, "6006 triangles from 2^-500 to 2^500 (a deep tree): 2002 queries, first hit %d, hit count %d, nearest %d differ from brute force", fh, ch, nn);
        mk_bvh_free(&B);
    }
    {
        static const double C8[] = {-1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1};
        static const int Ta[] = {0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 6, 3, 6, 2, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5};
        static const int Tb[] = {0, 2, 1, 0, 3, 2, 4, 5, 7, 5, 6, 7, 0, 1, 5, 0, 5, 4, 3, 7, 2, 2, 7, 6, 0, 4, 7, 0, 7, 3, 1, 2, 6, 1, 6, 5};
        int in = 0;
        for (int m = 0; m < 2; m++) {
            mk_bvh B;
            mk_bvh_build(&B, C8, m ? Tb : Ta, 12);
            in += mk_bvh_inside(&B, (double[3]){0, 0, 0});
            mk_bvh_free(&B);
        }
        CHECK(in == 2, "the centre of a cube is inside, with either diagonal on its faces (rays through edges don't vote): %d of 2", in);
    }
    {
        mk_bvh B;
        const int r = mk_bvh_build(&B, V, T, 1 << 30), r2 = mk_bvh_build(&B, V, T, -2147483643);
        CHECK(r == -1 && r2 == -1, "a tree of 2^30 or a negative count of triangles is refused before 2 nt + 1 overflows (%d, %d)", r, r2);
    }
    {
        double o[3] = {-0x1.854316adb7bdp+0, 0x1.7c7cae7cb20dfp+0, 0x1.058a117099cfp-3}, d[3] = {0, -2, 1};
        double a[3] = {0x1p-42, 0, 0}, b[3] = {0x1.004189374bc6ap-42, 0, 0}, c[3] = {0x1p-42, 0x1.0624dd2f1a9fcp-52, 0}, t, br[3];
        mk_ray r = mk_ray_make(o, d);
        const int h = mk_ray_tri(&r, a, b, c, 0, INFINITY, &t, br);
        CHECK(h == 0, "a ray 1.5 from a triangle of 2^-52 is no hit (%d)", h);
    }
    {
        double p[3] = {-0x1.a60414e4274d8p-3, -0x1.a15ed8c36a34p-5, -0x1.b6ef39b22ebb8p-2}, a[3] = {-0x1.391f6886a282p-1, 0x1.a28cd13d467dcp-2, -0x1.164dfaf3d1da8p-1};
        double b[3] = {0x1.133ee6bbe6934p-2, -0x1.32a5d083beb92p-1, -0x1.883c48ab7b0bp-3}, c[3] = {-0x1.a87a926d20e5cp-2, 0x1.7737744192f6cp-3, -0x1.dbed1386be45ep-2};
        double w[3], q[3], d = mk_closest_pt_tri(p, a, b, c, w, q);
        double p2[3] = {-4, -7, 8}, a2[3] = {7, -6, 3}, b2[3] = {7, 6, 7}, c2[3] = {7, 0x1.a1a286ea9c3b8p+0, 0x1.62cd8b3e37afap+2}, w2[3], q2[3];
        double d2 = mk_closest_pt_tri(p2, a2, b2, c2, w2, q2);
        CHECK(fabs(d - 0.0020183485327915739) <= 1e-15 && fabs(d2 - 146.6) <= 1e-12, "closest point on two slivers: %.17g (exact 0.0020183485327915739), %.17g (exact 146.6)", d, d2);
    }
    {
        double p[3] = {-8, 2, 7}, a[3] = {2, -5, -4}, c[3] = {0x1.fffffffffffffp+0, -5, -4}, w[3], q[3];
        double d = mk_closest_pt_tri(p, a, a, c, w, q), e = c[0] - a[0], t = fmin(1, fmax(0, (p[0] - a[0]) / e)), s = 0;
        for (int k = 0; k < 3; k++) { double x = a[k] + (k == 0 ? t * e : 0) - p[k]; s += x * x; }
        CHECK(fabs(d - s) <= 1e-12 * s && w[0] == w[0] && w[1] == w[1] && w[2] == w[2], "closest point on a triangle with a repeated corner: %.17g (the segment: %.17g), weights not NaN", d, s);
    }
    {
        double o[3] = {0x1.14699deba4d1cp+0, 0x1.f2738091721f4p-1, 0x1.1d966d08fd674p+0}, d[3] = {-0x1.1473a42358bc4p+0, -0x1.f273836cf2dcp-1, -0x1.c2db25017bf6ep-1};
        double a[3] = {0x1.cfbe9f82ab97p+29, 0x1.21041d856fecep+2, -0x1.d8b7e558bb0bap-31}, b[3] = {0x1.7a11b7654584p+35, -0x1.968e0aa823c26p+6, 0x1.a8cb856ba0c14p-27};
        double c[3] = {-0x1.9764ea4020f28p+34, 0x1.78ac2f5dc9f2p-2, -0x1.f4f9f1e772d54p-6}, t = -1, br[3];
        mk_ray r = mk_ray_make(o, d);
        const int h = mk_ray_tri(&r, a, b, c, 0, INFINITY, &t, br);
        CHECK(h && fabs(t - 1.2705559762972542) < 1e-6, "a 5e10-wide triangle 1.5 from the origin: hit %d at t %.17g (exact 1.2705559762972542)", h, t);
    }
    {
        double o[3] = {0x1.eceec5b961cb7p+0, -0x1.20a08cfe2cf0ep-1, -0x1.f8eae8d2b34adp+0}, d[3] = {-0x1.b27ac90e4303cp-2, 0, -0x1.d3f80b09f012p-1};
        double W[9] = {-0x1.64606b710a6b2p-21, -0x1.22bfc6c7ab27p-5, 0x1.65268a1e06dbp+37, 0x1.0c555bffdd6p-41, 0x1.b81c95c70b5dp-25, -0x1.53211865b5c32p+28,
                       -0x1.6cb7aeb8eb4bp-26, -0x1.ad42ec26c46d2p+35, 0x1.b39486183aed8p+8};
        const int T1[3] = {0, 1, 2};
        const double tmin = 0x1.2267035f6e37cp+2;
        mk_ray r = mk_ray_make(o, d);
        double t = -1, tb = -1;
        const int h = mk_ray_tri(&r, W, W + 3, W + 6, tmin, INFINITY, &t, NULL);
        mk_bvh B;
        mk_bvh_build(&B, W, T1, 1);
        const int hb = mk_bvh_first_hit(&B, o, d, tmin, INFINITY, &tb, NULL);
        mk_bvh_free(&B);
        mk_aabb bx;
        double te;
        mk_aabb_empty(&bx), mk_aabb_add(&bx, W), mk_aabb_add(&bx, W + 3), mk_aabb_add(&bx, W + 6);
        const int in = h ? mk_aabb_ray(&bx, &r, tmin, t, &te) : -1;
        CHECK(h && in == 1 && hb == 0 && tb == t, "a thin triangle 1e11 wide: its own box holds the hit (%d), and the tree finds it (%d, %d)", in, h, hb);
    }
    { /* two faces of a scaled cube tie exactly at p (a case found by a random search): the lower index wins */
        static const double Q[] = {0x0p+0, -0x1.1p+11, 0x1.8p+8, 0x0p+0, -0x1.ep+10, 0x1.8p+8, 0x0p+0, -0x1.ep+10, 0x1.4p+9, 0x0p+0, -0x1.1p+11, 0x1.4p+9,
                                   0x1p+8, -0x1.1p+11, 0x1.8p+8, 0x1p+8, -0x1.ep+10, 0x1.8p+8, 0x1p+8, -0x1.ep+10, 0x1.4p+9, 0x1p+8, -0x1.1p+11, 0x1.4p+9};
        static const int QT[] = {0, 3, 1, 1, 3, 2, 4, 5, 7, 5, 6, 7, 0, 7, 3, 0, 4, 7, 1, 2, 6, 1, 6, 5, 0, 1, 4, 4, 1, 5, 3, 7, 2, 7, 6, 2};
        const double p[3] = {0x1.a343726fd65b3p-1, 0x1.82eadc5d5d982p+0, 0x1.5db8bcff5315fp+0};
        double bd = INFINITY, q[3];
        int bn = -1, tn;
        for (int s = 0; s < 12; s++) {
            double w[3], c[3], d2 = mk_closest_pt_tri(p, Q + 3 * QT[3 * s], Q + 3 * QT[3 * s + 1], Q + 3 * QT[3 * s + 2], w, c);
            if (d2 < bd) bd = d2, bn = s;
        }
        mk_bvh B;
        mk_bvh_build(&B, Q, QT, 12);
        const double nd = mk_bvh_nearest(&B, p, INFINITY, &tn, q);
        mk_bvh_free(&B);
        CHECK(tn == bn && nd == bd && bn == 7, "an exact tie in the nearest distance goes to the lower index, as in brute force: the tree %d, brute force %d", tn, bn);
    }
    { /* mk_bvh_box against every triangle's box, on the torus, with boxes that only touch counted as meeting */
        torus(64, 32, 2, 0.7);
        mk_bvh B;
        mk_bvh_build(&B, V, T, nt);
        static int out[8192];
        int bad = 0, tot = 0;
        for (int k = 0; k < 2000; k++) {
            mk_aabb q;
            double c[3] = {3 * t_rand(), 3 * t_rand(), t_rand()}, h = k % 4 ? 0.3 * fabs(t_rand()) : 0;
            for (int i = 0; i < 3; i++) q.lo[i] = c[i] - h, q.hi[i] = c[i] + h;
            if (k % 5 == 0) { /* a box whose low x face is exactly a triangle's high x */
                mk_aabb tb;
                mk_tri_box_(V, T, k % nt, &tb);
                q.lo[0] = tb.hi[0], q.hi[0] = tb.hi[0] + 0.1;
                for (int i = 1; i < 3; i++) q.lo[i] = tb.lo[i], q.hi[i] = tb.hi[i];
            }
            int n = mk_bvh_box(&B, &q, out, 8192), m = 0;
            for (int s = 0; s < nt; s++) {
                mk_aabb tb;
                mk_tri_box_(V, T, s, &tb);
                int meet = 1;
                for (int i = 0; i < 3; i++) meet &= tb.lo[i] <= q.hi[i] && q.lo[i] <= tb.hi[i];
                if (!meet) continue;
                m++;
                int found = 0;
                for (int j = 0; j < n && j < 8192; j++) found |= out[j] == s;
                bad += !found;
            }
            bad += n != m, tot += m;
        }
        CHECK(bad == 0 && tot > 1000 && mk_bvh_box(&B, &(mk_aabb){{9, 9, 9}, {10, 10, 10}}, NULL, 0) == 0, "box query against every triangle's box: 2000 boxes, %d triangles found, %d wrong", tot, bad);
        mk_bvh_free(&B);
    }
    { /* mk_bvh_crossings against brute force over every pair, on two spheres sunk into each other as one mesh: all
         pairs, and only those with a face marked (centroid x > 0.6, about where they cross) */
        tm_mesh M = {0};
        const double c0[3] = {0, 0, 0}, c1[3] = {1.2, 0.1, 0.05};
        tm_sphere(&M, 1, c0, 3), tm_sphere(&M, 1, c1, 3);
        unsigned char *mark = malloc((size_t)M.nt);
        for (int t = 0; t < M.nt; t++) mark[t] = (M.V[3 * M.T[3 * t]] + M.V[3 * M.T[3 * t + 1]] + M.V[3 * M.T[3 * t + 2]]) / 3 > 0.6;
        long all = 0, marked = 0;
        for (int a = 0; a < M.nt; a++)
            for (int b = a + 1; b < M.nt; b++) {
                const int *X = M.T + 3 * a, *Y = M.T + 3 * b;
                int sh = 0, x = 0;
                for (int i = 0; i < 9; i++) sh |= X[i / 3] == Y[i % 3];
                for (int k = 0; k < 3 && !sh && !x; k++)
                    x = segtri_ref(M.V + 3 * X[k], M.V + 3 * X[(k + 1) % 3], M.V + 3 * Y[0], M.V + 3 * Y[1], M.V + 3 * Y[2]) ||
                        segtri_ref(M.V + 3 * Y[k], M.V + 3 * Y[(k + 1) % 3], M.V + 3 * X[0], M.V + 3 * X[1], M.V + 3 * X[2]);
                all += x, marked += x && (mark[a] || mark[b]);
            }
        mk_bvh B;
        const int rb = mk_bvh_build(&B, M.V, M.T, M.nt);
        const long ca = rb ? -2 : mk_bvh_crossings(&B, NULL), cm = rb ? -2 : mk_bvh_crossings(&B, mark);
        if (!rb) mk_bvh_free(&B);
        CHECK(ca == all && cm == marked && marked > 0 && marked < all, "crossing pairs of two spheres sunk into each other: %ld (brute force %ld); with a face marked %ld (brute force %ld)", ca, all, cm, marked);
        free(mark), tm_free(&M);
    }
    { /* the hits' facings summed along a ray (mk_bvh_count_wind_, which mk_sepvol uses as the winding number on a
         closed mesh) against mk_winding_number, on the same two spheres (covered twice where they overlap): at 3000
         points in a box round them, every ray with no edge hit sums to the winding number rounded, and counts what
         mk_bvh_count_hits counts */
        tm_mesh M = {0};
        const double c0[3] = {0, 0, 0}, c1[3] = {1.2, 0.1, 0.05};
        tm_sphere(&M, 1, c0, 3), tm_sphere(&M, 1, c1, 3);
        mk_bvh B;
        const int rb = mk_bvh_build(&B, M.V, M.T, M.nt);
        int bad = rb != 0, nray = 0, cls[3] = {0, 0, 0};
        double worst = 0;
        for (int k = 0; k < 3000 && !rb; k++) {
            const double p[3] = {0.6 + 1.8 * t_rand(), 1.2 * t_rand(), 1.2 * t_rand()}, wn = mk_winding_number(M.V, M.T, M.nt, p);
            const double d[3] = {t_rand(), t_rand(), t_rand()};
            int edge, w, c = mk_bvh_count_wind_(&B, p, d, &edge, &w);
            if (edge) continue;
            nray++, worst = fmax(worst, fabs(wn - w));
            bad += fabs(wn - w) > 1e-6 || c != mk_bvh_count_hits(&B, p, d) || (c - w) % 2;
            if (w >= 0 && w <= 2) cls[w]++;
        }
        if (!rb) mk_bvh_free(&B);
        CHECK(bad == 0 && nray > 2900 && cls[0] > 100 && cls[1] > 100 && cls[2] > 100, "the facings summed along a ray: the winding number at %d of %d points (0: %d, 1: %d, 2: %d; largest difference %.2g), the hits mk_bvh_count_hits counts",
              nray - bad, nray, cls[0], cls[1], cls[2], worst);
        tm_free(&M);
    }
    { /* where no ray votes (each of mk_bvh_inside's three fixed rays hits an edge or a vertex), the winding number
         decides by its magnitude, not its sign, so a closed mesh wound inward still reads inside, as the parity says
         elsewhere. A tetrahedron with a vertex on each fixed ray from the origin, both windings. */
        /* mk_bvh_inside's three ray directions */
        static const double dirs[3][3] = {{0.5773502691896258, 0.5773502691896258, 0.5773502691896258},
                                          {-0.3015113445777636, 0.9045340337332909, 0.3015113445777636},
                                          {0.2672612419124244, -0.5345224838248488, 0.8017837257372732}};
        double V[12];
        for (int k = 0; k < 3; k++)
            for (int a = 0; a < 3; a++) V[3 * k + a] = 2 * dirs[k][a];
        for (int a = 0; a < 3; a++) V[9 + a] = -(dirs[0][a] + dirs[1][a] + dirs[2][a]);
        const int To[12] = {0, 1, 2, 0, 3, 1, 1, 3, 2, 2, 3, 0}, Ti[12] = {0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3};
        const double o[3] = {0, 0, 0};
        mk_bvh Bo, Bi;
        mk_bvh_build(&Bo, V, To, 4), mk_bvh_build(&Bi, V, Ti, 4);
        int ed = 0;
        for (int k = 0; k < 3; k++) {
            int e;
            mk_bvh_count_wind_(&Bo, o, dirs[k], &e, NULL);
            ed += e != 0;
        }
        const double wo = mk_winding_number(V, To, 4, o), wi = mk_winding_number(V, Ti, 4, o);
        const int io = mk_bvh_inside(&Bo, o), ii = mk_bvh_inside(&Bi, o);
        CHECK(ed == 3 && wo > 0.5 && wi < -0.5 && io == 1 && ii == 1, "no ray votes: the origin inside a tetrahedron whose vertices sit on the 3 rays (%d through a vertex), wound either way: winding %g and %g, inside %d and %d",
              ed, wo, wi, io, ii);
        mk_bvh_free(&Bo), mk_bvh_free(&Bi);
    }
    { /* mk_bvh_free zeroes nt, so a query after a free (or after a build that failed for memory) finds no tree instead
         of reading a NULL node */
        static const double V[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
        static const int T[3] = {0, 1, 2};
        const double o[3] = {0.2, 0.2, 1}, d[3] = {0, 0, -1};
        double t;
        mk_bvh B;
        mk_bvh_build(&B, V, T, 1);
        const int h0 = mk_bvh_first_hit(&B, o, d, 0, INFINITY, &t, NULL);
        mk_bvh_free(&B);
        const int nt = B.nt, h1 = nt ? -9 : mk_bvh_first_hit(&B, o, d, 0, INFINITY, &t, NULL);
        CHECK(h0 == 0 && nt == 0 && h1 == -1, "a freed tree has nt %d and a query on it finds nothing (%d; before the free %d)", nt, h1, h0);
    }
    return wrong;
}
