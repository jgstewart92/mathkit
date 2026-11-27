/* mathkit: a closed triangle mesh cut along planes into closed pieces that join exactly. Header only, C23.
   Each cut is a plane, a point P on it and a normal N. Two calls choose which of the plane's loops (the closed curves in
   which it meets the surface) are cut:
     mk_cut_loops  only the loop that winds round P in the plane, the one of least area if several do; other loops of the
                   same plane (other coils of a tube) stay whole. It must reach no farther than rmax from P, and no other
                   loop of the plane may lie inside it. If the plane through P fails, it is moved along N by step, 2 step,
                   ... up to shift, each way (+ first), and the first place that passes is used (out->Pc).
     mk_cut_clip   every loop of each plane (a clip: the caller keeps whichever pieces it wants); no other loop of the
                   plane may lie inside one.
   Then, in both: a vertex nearer the plane than snap is moved onto it (along N), and a loop through it passes through
   it, when that leaves one run of + neighbours round it. Where the plane passes beside a vertex twice (a saddle: two
   runs), on the plane both passes would meet at it and pinch the loop, so it is left where it is; a vertex with only +
   neighbours joins their side (a dimple thinner than snap is not cut). Every other crossed edge gets one new vertex at
   the plane, held between 0.001 and 0.999 of its edge (MADE UP), so no triangle of near-zero area is made beside the loop. The crossed triangles are split on the loops, the surface falls
   into regions bounded by the loops (a flood fill across every edge but the loops'), and each loop gets one cap, ear
   clipped in its plane on the loop's own vertices, the best-shaped valid ear first (area over the squared edge lengths,
   so a run of near-collinear points is fanned from a corner, not clipped into slivers). A cap is wound so that the region
   on the - side of its plane takes it as it is and the one on the + side reversed (+N for an outward input), so two
   pieces share the cap's vertices and triangles exactly, wound as the input, and their volumes sum to the whole. mk_cut_piece gives a region's triangles with its caps.
   The source of each new vertex: it is V[a] + t (V[b] - V[a]) of the input's V, computed in that form, a and b the ends
   of its crossed edge (out->src[2 i], out->src[2 i + 1], t = out->frac[i] for output vertex nv + i), so a caller carries
   weights or any per-vertex value along by the same sum. A snapped vertex keeps its index (out->snapped lists them); it
   moved by under snap.
   Refusals (return codes over 0, with out->reason; out->bad the cut it names, or -1): MK_CUT_OPEN the input is not a
   closed, consistently wound manifold (an edge not in two triangles, or in two used the same way round); MK_CUT_NOLOOP
   no loop winds round P (mk_cut_clip: the plane misses the mesh); MK_CUT_FAR the loop reaches past rmax; MK_CUT_NESTED
   another loop of the plane lies inside the loop (its cap would run through that tube); MK_CUT_CROSS two cuts' loops
   through one triangle or one snapped vertex; MK_CUT_REGIONS (mk_cut_loops) regions other than ncut + 1 (a loop that
   does not separate, or a mesh not in one piece); MK_CUT_CAP a loop not a simple polygon in its plane; MK_CUT_SPLIT the
   split surface or a piece not closed (a snap large against the edges near a cut), or a piece inside out (a cut's point
   within shift of the surface): every piece is checked before 0 is returned. -1: a bad argument (a coordinate
   not finite, an index out of range, a zero normal, rmax not over 0) or no memory. On any return but 0, out holds only
   reason, code and bad.
   A cap is flat to rounding, but where a new vertex is held at 0.001 or 0.999 of its edge: then it lies off the plane by
   under 0.001 of the edge's height across it (beside a saddle left in place, or on an edge over 1000 snap long).
   Not checked here: a cap crossing the surface or another cap (a plane grazing the surface beside the loop). Run a mesh
   check on the pieces. (A checker that reads exactly coplanar triangles near each other as crossing will find such
   pairs among a cap's own triangles written at full precision.)
   The input: a closed manifold, consistently wound (each edge in two triangles, used once each way). Single-threaded;
   memory about 90 bytes a triangle and 8 a vertex for each cut (a 13.9 M-triangle mesh in 8 cuts: 1.8 GB beside its
   input). Units are the caller's; every MADE UP default is mk_cut_defaults'. */
#ifndef MK_CUT_H
#define MK_CUT_H
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { MK_CUT_OK = 0, MK_CUT_OPEN = 1, MK_CUT_NOLOOP = 2, MK_CUT_FAR = 3, MK_CUT_NESTED = 4, MK_CUT_CROSS = 5, MK_CUT_REGIONS = 6, MK_CUT_CAP = 7, MK_CUT_SPLIT = 8 };

typedef struct {
    double snap;  /* a vertex nearer a plane than this is moved onto it (the top comment says when) */
    double shift; /* mk_cut_loops: the farthest a plane moves along N to find a loop it may cut (0: never moved) */
    double step;  /* in steps of this */
} mk_cut_opts;
/* Defaults for a mesh in which one millimetre is mm units (metres: 1e-3). All MADE UP: snap 5 um (a vertex that near a
   plane, left where it is, gives a corner triangle of order (snap / edge)^2 of its triangle's area), shift 0.1 mm in
   steps of 10 um. */
static inline void mk_cut_defaults(mk_cut_opts *o, double mm)
{
    memset(o, 0, sizeof *o);
    o->snap = 0.005 * mm, o->shift = 0.1 * mm, o->step = 0.01 * mm;
}

/* The result. Loops: in mk_cut_loops loop k is cut k's; in mk_cut_clip a plane's loops in the order found. */
typedef struct {
    double *V; int nv, nv0;   /* the input's vertices (snapped ones moved), then the new ones; nv0 the input's count */
    int *T; int nt;           /* the split surface */
    int *reg; int nreg;       /* per triangle of T: its region */
    int *capT; int ncapt;     /* cap triangles, as the - side's region takes them (+N of their plane for an outward input) */
    int *capk;                /* per cap triangle: its loop */
    int nloop;
    int *loopk;               /* per loop: its cut (plane) */
    int *side;                /* per loop: [2 l] the region on the - side, [2 l + 1] on the + side */
    int *loopn;               /* per loop: its steps (crossed triangles) */
    double *loopr, *looparea; /* per loop: its farthest point from its cut's Pc, and its area in the plane */
    double *margin;           /* per loop: the least fraction along a crossed edge (0 where it passes a snapped vertex) */
    double *Pc;               /* per cut: the plane point used (P, or P moved along N) */
    int *src; double *frac;   /* per new vertex nv0 + i: its edge's ends src[2 i], src[2 i + 1] and t along it */
    int *snapped; int nsnapped; /* the input vertices moved onto a plane */
    int code, bad;            /* the return code; the cut a refusal names (-1: none) */
    char reason[200];
} mk_cut_out;

static inline void mk_cut_free(mk_cut_out *o)
{
    free(o->V), free(o->T), free(o->reg), free(o->capT), free(o->capk), free(o->loopk), free(o->side), free(o->loopn);
    free(o->loopr), free(o->looparea), free(o->margin), free(o->Pc), free(o->src), free(o->frac), free(o->snapped);
    memset(o, 0, sizeof *o);
}

typedef struct { int64_t key; int t; } mk_cut_he_;
static inline int mk_cut_hecmp_(const void *a, const void *b)
{
    const int64_t x = ((const mk_cut_he_ *)a)->key, y = ((const mk_cut_he_ *)b)->key;
    return x < y ? -1 : x > y;
}
/* neighbours across each triangle's edges (nb[3 t + e], e the edge from corner e to e + 1), -1 where an edge is not in
   exactly two triangles: 0, 1 no memory, 2 not a closed manifold */
static inline int mk_cut_nbrs_(const int *T, int nt, int nv, int *nb)
{
    mk_cut_he_ *h = malloc(sizeof *h * 3 * (size_t)nt);
    if (!h) return 1;
    for (int t = 0; t < nt; t++)
        for (int e = 0; e < 3; e++) {
            const int64_t a = T[3 * t + e], b = T[3 * t + (e + 1) % 3];
            h[3 * t + e] = (mk_cut_he_){(a < b ? a : b) * (int64_t)nv + (a < b ? b : a), 3 * t + e};
        }
    qsort(h, 3 * (size_t)nt, sizeof *h, mk_cut_hecmp_);
    int bad = 0;
    for (size_t i = 0; i < 3 * (size_t)nt; i++) nb[i] = -1;
    for (size_t i = 0; i < 3 * (size_t)nt;) {
        size_t j = i;
        while (j < 3 * (size_t)nt && h[j].key == h[i].key) j++;
#if defined(MK_MUTATE) && MK_MUTATE == 444
        if (j - i == 2) /* 444: two triangles using an edge the same way round accepted */
#else
        if (j - i == 2 && T[h[i].t] != T[h[i + 1].t]) /* the two uses run opposite ways: consistently wound */
#endif
            nb[h[i].t] = h[i + 1].t / 3, nb[h[i + 1].t] = h[i].t / 3;
        else bad = 1;
        i = j;
    }
    free(h);
#if defined(MK_MUTATE) && MK_MUTATE == 441
    bad = 0; /* 441: an open mesh not refused */
#endif
    return bad ? 2 : 0;
}
static inline int mk_cut_intcmp_(const void *a, const void *b) { const int x = *(const int *)a, y = *(const int *)b; return (x > y) - (x < y); }
static inline double mk_cut_dot_(const double *a, const double *b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
static inline double mk_cut_t_(double da, double db) { double t = da / (da - db); return t < 0.001 ? 0.001 : t > 0.999 ? 0.999 : t; }
static inline int mk_cut_root_(int *p, int a) { while (p[a] != a) a = p[a] = p[p[a]]; return a; }
/* the unit normal nz of n and a basis (u, w) of the plane, u x w = nz */
static inline void mk_cut_basis_(const double *n, double *nz, double *u, double *w)
{
    const double nn = sqrt(mk_cut_dot_(n, n));
    nz[0] = n[0] / nn, nz[1] = n[1] / nn, nz[2] = n[2] / nn;
    const double a[3] = {fabs(nz[0]) < 0.9 ? 1 : 0, fabs(nz[0]) < 0.9 ? 0 : 1, 0};
    u[0] = a[1] * nz[2] - a[2] * nz[1], u[1] = a[2] * nz[0] - a[0] * nz[2], u[2] = a[0] * nz[1] - a[1] * nz[0];
    const double ul = sqrt(mk_cut_dot_(u, u));
    u[0] /= ul, u[1] /= ul, u[2] /= ul;
    w[0] = nz[1] * u[2] - nz[2] * u[1], w[1] = nz[2] * u[0] - nz[0] * u[2], w[2] = nz[0] * u[1] - nz[1] * u[0];
}
/* ear clipping of a simple polygon q[n] (2-D, counter-clockwise) into tri[n - 2][3], indices into q: 0, or 1 if it is not
   simple, -1 for no memory. The best-shaped valid ear each time (twice its area over its squared edge lengths); only reflex corners can lie
   inside an ear. */
static inline int mk_cut_ears0_(int n, const double (*q)[2], int (*tri)[3])
{
    int *id = malloc(sizeof *id * (size_t)n), m = n, nt = 0;
    if (!id) return -1;
    for (int i = 0; i < n; i++) id[i] = i;
    while (m > 3) {
        int bi = -1;
        double bq = 0;
        for (int i = 0; i < m; i++) {
            const int a = id[(i + m - 1) % m], b = id[i], c = id[(i + 1) % m];
            const double cr = (q[b][0] - q[a][0]) * (q[c][1] - q[a][1]) - (q[b][1] - q[a][1]) * (q[c][0] - q[a][0]);
            if (cr <= 0) continue;
            const double l2 = (q[b][0] - q[a][0]) * (q[b][0] - q[a][0]) + (q[b][1] - q[a][1]) * (q[b][1] - q[a][1]) + (q[c][0] - q[b][0]) * (q[c][0] - q[b][0]) +
                              (q[c][1] - q[b][1]) * (q[c][1] - q[b][1]) + (q[a][0] - q[c][0]) * (q[a][0] - q[c][0]) + (q[a][1] - q[c][1]) * (q[a][1] - q[c][1]);
            const double qual = cr / l2;
            if (qual <= bq) continue;
            int ear = 1;
            for (int j = 0; ear && j < m; j++) {
                const int pa = id[(j + m - 1) % m], p = id[j], pc = id[(j + 1) % m];
                if (p == a || p == b || p == c) continue;
                const double rc = (q[p][0] - q[pa][0]) * (q[pc][1] - q[pa][1]) - (q[p][1] - q[pa][1]) * (q[pc][0] - q[pa][0]);
                if (rc > 0) continue; /* a convex corner cannot be inside an ear */
                const double d1 = (q[b][0] - q[a][0]) * (q[p][1] - q[a][1]) - (q[b][1] - q[a][1]) * (q[p][0] - q[a][0]);
                const double d2 = (q[c][0] - q[b][0]) * (q[p][1] - q[b][1]) - (q[c][1] - q[b][1]) * (q[p][0] - q[b][0]);
                const double d3 = (q[a][0] - q[c][0]) * (q[p][1] - q[c][1]) - (q[a][1] - q[c][1]) * (q[p][0] - q[c][0]);
                if (d1 >= 0 && d2 >= 0 && d3 >= 0) ear = 0;
            }
            if (ear) bi = i, bq = qual;
#if defined(MK_MUTATE) && MK_MUTATE == 438
            if (bi >= 0) break; /* 438: the first valid ear, not the best-shaped */
#endif
        }
        if (bi < 0) break;
        tri[nt][0] = id[(bi + m - 1) % m], tri[nt][1] = id[bi], tri[nt][2] = id[(bi + 1) % m], nt++;
        memmove(id + bi, id + bi + 1, sizeof *id * (size_t)(m - 1 - bi));
        m--;
    }
    if (m == 3) tri[nt][0] = id[0], tri[nt][1] = id[1], tri[nt][2] = id[2], nt++;
    free(id);
    return nt == n - 2 ? 0 : 1;
}
/* tri[nt] covers q[n] exactly: every triangle counter-clockwise, their areas summing to the polygon's */
static inline int mk_cut_ears_ok_(int n, const double (*q)[2], const int (*tri)[3], int nt)
{
    double a = 0, s = 0;
    for (int i = 0; i < n; i++) a += q[i][0] * q[(i + 1) % n][1] - q[(i + 1) % n][0] * q[i][1];
    for (int t = 0; t < nt; t++) {
        const double *x = q[tri[t][0]], *y = q[tri[t][1]], *z = q[tri[t][2]];
        const double c = (y[0] - x[0]) * (z[1] - x[1]) - (y[1] - x[1]) * (z[0] - x[0]);
        if (!(c > 0)) return 0;
        s += c;
    }
    return fabs(s - a) <= 1e-9 * fabs(a);
}
/* The clipping above, kept where it covers the polygon; where it fails or folds (a polygon with corners collinear to
   rounding, as where a loop crosses a flat fan or grid: a point on an ear's side blocks that ear, and a triangle with
   three points a side had none), those flat corners (turning by under 1e-9) are taken out, the rest clipped, and each
   run put back fanned from the opposite corner of the triangle on its side, so every corner is used and each side of
   the polygon is a side of one triangle. */
static inline int mk_cut_ears_(int n, const double (*q)[2], int (*tri)[3])
{
    const int e = mk_cut_ears0_(n, q, tri);
    if (e < 0 || (!e && mk_cut_ears_ok_(n, q, (const int (*)[3])tri, n - 2))) return e;
    int *nx = malloc(sizeof *nx * (size_t)n), *pv = malloc(sizeof *pv * (size_t)n), *id = malloc(sizeof *id * (size_t)n), left = n, m = 0, rc = 1;
    char *gone = calloc((size_t)n, 1);
    double (*r)[2] = malloc(sizeof *r * (size_t)n);
    if (!nx || !pv || !id || !gone || !r) { rc = -1; goto done; }   /* no memory (the core's code too: 1 not simple, -1 no memory) */
    for (int i = 0; i < n; i++) nx[i] = (i + 1) % n, pv[i] = (i + n - 1) % n;
    for (int changed = 1; changed && left > 3;) {
        changed = 0;
        for (int i = 0; i < n && left > 3; i++) {
            if (gone[i]) continue;
            const double *a = q[pv[i]], *b = q[i], *c = q[nx[i]];
            const double cr = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            const double l1 = (b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]), l2 = (c[0] - b[0]) * (c[0] - b[0]) + (c[1] - b[1]) * (c[1] - b[1]);
            if (fabs(cr) <= 1e-9 * sqrt(l1 * l2) && (c[0] - b[0]) * (b[0] - a[0]) + (c[1] - b[1]) * (b[1] - a[1]) > 0)
                gone[i] = 1, nx[pv[i]] = nx[i], pv[nx[i]] = pv[i], left--, changed = 1;
        }
    }
    double area = 0;
    for (int i = 0; i < n; i++) if (!gone[i]) id[m] = i, r[m][0] = q[i][0], r[m][1] = q[i][1], m++;
    for (int i = 0; i < m; i++) area += r[i][0] * r[(i + 1) % m][1] - r[(i + 1) % m][0] * r[i][1];
    if (m < 3 || !(area > 0)) goto done;
    const int e0 = mk_cut_ears0_(m, (const double (*)[2])r, tri);
    if (e0) { rc = e0; goto done; }
    int nt = m - 2;
    for (int t = 0; t < nt; t++) for (int j = 0; j < 3; j++) tri[t][j] = id[tri[t][j]];
    for (int s = 0; s < m; s++) {
        const int u = id[s], v = id[(s + 1) % m];
        if ((u + 1) % n == v) continue;
        int t = 0, j = 3;
        for (; t < nt; t++) {
            for (j = 0; j < 3 && !(tri[t][j] == u && tri[t][(j + 1) % 3] == v); j++) {}
            if (j < 3) break;
        }
        if (t == nt) goto done;
        const int w = tri[t][(j + 2) % 3];
        tri[t][0] = u, tri[t][1] = (u + 1) % n, tri[t][2] = w;
        for (int x = (u + 1) % n; x != v; x = (x + 1) % n) tri[nt][0] = x, tri[nt][1] = (x + 1) % n, tri[nt][2] = w, nt++;
    }
    rc = nt == n - 2 && mk_cut_ears_ok_(n, q, (const int (*)[3])tri, nt) ? 0 : 1;
done:
    free(nx), free(pv), free(id), free(gone), free(r);
    return rc;
}
/* the vertices' distances to the plane (p, unit nz). A vertex nearer than snap gets 0, which reads as on the plane and,
   in the sides, as the - side, when that leaves one run of + neighbours round it (or none). Two runs, a saddle the plane
   passes beside: on the plane both passes would meet at it and pinch the loop, so it keeps its own distance (it is left
   where it is). Only + neighbours: it joins their side (a dimple thinner than snap is not cut). Decided again until no
   vertex changes, since a vertex's side counts in its neighbours' runs. 0, or -1 for no memory. */
typedef struct { int v, t, j; } mk_cut_vt_;
static inline int mk_cut_vtcmp_(const void *a, const void *b)
{
    const mk_cut_vt_ *x = a, *y = b;
    return x->v != y->v ? (x->v > y->v) - (x->v < y->v) : (x->t > y->t) - (x->t < y->t);
}
static inline int mk_cut_dist_(const double *V, int nv, const int *T, int nt, const double *p, const double *nz, double snap, double *d)
{
#if defined(MK_MUTATE) && MK_MUTATE == 429
    snap = 0; /* 429: no snapping */
#endif
    long ns = 0;
    for (int i = 0; i < nv; i++) {
        const double r[3] = {V[3 * i] - p[0], V[3 * i + 1] - p[1], V[3 * i + 2] - p[2]};
        d[i] = mk_cut_dot_(r, nz);
        if (fabs(d[i]) < snap) d[i] = 0, ns++;
    }
#if defined(MK_MUTATE) && MK_MUTATE == 442
    ns = 0; /* 442: every vertex within snap onto the plane, saddles too */
#endif
    if (!ns) return 0;
    long np = 0, cp = 6 * ns + 16;
    mk_cut_vt_ *a = malloc(sizeof *a * (size_t)cp);
    if (!a) return -1;
    for (int t = 0; t < nt; t++)
        for (int j = 0; j < 3; j++) {
            if (d[T[3 * t + j]] != 0) continue;
            if (np == cp) {
                mk_cut_vt_ *a2 = realloc(a, sizeof *a * (size_t)(cp *= 2));
                if (!a2) { free(a); return -1; }
                a = a2;
            }
            a[np++] = (mk_cut_vt_){T[3 * t + j], t, j};
        }
    qsort(a, (size_t)np, sizeof *a, mk_cut_vtcmp_);
    for (int pass = 0, changed = 1; changed && pass < 64; pass++) {
        changed = 0;
        for (long i = 0; i < np;) {
            long e = i;
            while (e < np && a[e].v == a[i].v) e++;
            const int v = a[i].v;
            if (a[i].j >= 0 && d[v] == 0) {
                int runs = 0, plus = 0, minus = 0;
                for (long q = i; q < e; q++) {
                    const int t = a[q].t, j = a[q].j, su = d[T[3 * t + (j + 1) % 3]] > 0, sw = d[T[3 * t + (j + 2) % 3]] > 0;
                    runs += su != sw, plus += su + sw, minus += !su + !sw;
                }
                const double r[3] = {V[3 * v] - p[0], V[3 * v + 1] - p[1], V[3 * v + 2] - p[2]}, dv = mk_cut_dot_(r, nz);
                if (runs > 2) d[v] = dv;
                else if (plus && !minus) d[v] = dv > 0 ? dv : snap;
                if (d[v] != 0) {
                    changed = 1;
                    for (long q = i; q < e; q++) a[q].j = -1; /* decided */
                }
            }
            i = e;
        }
    }
    free(a);
    return 0;
}
static inline int mk_cut_on_(double d) { return d == 0; }
#if defined(MK_MUTATE) && MK_MUTATE == 434
enum { MK_CUT_PER_ = 3 }; /* (434: room for a triangle split twice) */
#else
enum { MK_CUT_PER_ = 2 }; /* triangles a loop step adds: three for its one */
#endif

/* A plane's loops, traced through the crossed triangles: loop l's triangles in order are tri[off[l] .. off[l + 1]), its
   points (each step's exit crossing, in the plane's basis about p) pt[off[l] ..], its signed area and its winding angle
   round p, its farthest point from p. */
typedef struct {
    int nl, cl, ntri, ctri;
    int *off, *tri;
    double (*pt)[2];
    double *area, *wind, *r;
    double u[3], w[3];
} mk_cut_lps_;
static inline void mk_cut_lps_free_(mk_cut_lps_ *L) { free(L->off), free(L->tri), free(L->pt), free(L->area), free(L->wind), free(L->r), memset(L, 0, sizeof *L); }
/* every loop of the plane (p, n); d[nv] gets the vertices' distances to it; seen[nt] (zeroed once) and *stamp are the
   caller's scratch, kept between calls. 0, or -1 for no memory. */
static inline int mk_cut_trace_(const double *V, int nv, const int *T, int nt, const int *nb, const double *p, const double *n, double snap, double *d, int *seen, int *stamp, mk_cut_lps_ *L)
{
    mk_cut_lps_free_(L);
    (*stamp)++;
    double nz[3];
    mk_cut_basis_(n, nz, L->u, L->w);
    if (mk_cut_dist_(V, nv, T, nt, p, nz, snap, d)) return -1;
    L->off = malloc(sizeof *L->off), L->cl = 0;
    if (!L->off) return -1;
    L->off[0] = 0;
    for (int t0 = 0; t0 < nt; t0++) {
        const int a = T[3 * t0], b = T[3 * t0 + 1], c = T[3 * t0 + 2];
        const int sa = d[a] > 0, sb = d[b] > 0, sc = d[c] > 0;
        if (sa == sb && sb == sc) continue;
        if (seen[t0] == *stamp) continue;
        /* the loop through t0: in by one crossed edge, out by the other */
        int m = 0, t = t0, from = -1, ok = 1;
        const int base = L->ntri;
        double area = 0, wind = 0, rr = 0, prev[2] = {0, 0}, first[2] = {0, 0};
        for (;;) {
            if (seen[t] == *stamp) { ok = t == t0; break; }
            seen[t] = *stamp;
            if (base + m == L->ctri) {
                const int c2 = 2 * L->ctri + 64;
                int *nt2 = realloc(L->tri, sizeof *nt2 * (size_t)c2);
                double (*np2)[2] = realloc(L->pt, sizeof *np2 * (size_t)c2);
                if (nt2) L->tri = nt2;
                if (np2) L->pt = np2;
                if (!nt2 || !np2) return -1;
                L->ctri = c2;
            }
            L->tri[base + m] = t, m++;
            int e_out = -1;
            for (int e = 0; e < 3; e++) {
                const int x = T[3 * t + e], y = T[3 * t + (e + 1) % 3];
                if ((d[x] > 0) != (d[y] > 0) && nb[3 * t + e] != from) { e_out = e; break; }
            }
            if (e_out < 0) { ok = 0; break; }
            const int x = T[3 * t + e_out], y = T[3 * t + (e_out + 1) % 3];
            const double tt = mk_cut_t_(d[x], d[y]);
            const double q[3] = {V[3 * x] + tt * (V[3 * y] - V[3 * x]) - p[0], V[3 * x + 1] + tt * (V[3 * y + 1] - V[3 * x + 1]) - p[1], V[3 * x + 2] + tt * (V[3 * y + 2] - V[3 * x + 2]) - p[2]};
            const double q2[2] = {mk_cut_dot_(q, L->u), mk_cut_dot_(q, L->w)};
            const double r = sqrt(q2[0] * q2[0] + q2[1] * q2[1]);
            if (r > rr) rr = r;
            L->pt[base + m - 1][0] = q2[0], L->pt[base + m - 1][1] = q2[1];
            if (m > 1) { area += 0.5 * (prev[0] * q2[1] - prev[1] * q2[0]); wind += atan2(prev[0] * q2[1] - prev[1] * q2[0], prev[0] * q2[0] + prev[1] * q2[1]); }
            else first[0] = q2[0], first[1] = q2[1];
            prev[0] = q2[0], prev[1] = q2[1];
            from = t, t = nb[3 * t + e_out];
            if (t < 0) { ok = 0; break; }
            if (t == t0) break;
        }
        if (!ok || m < 3) continue; /* (its triangles stay seen, its slots are reused) */
        area += 0.5 * (prev[0] * first[1] - prev[1] * first[0]);
        wind += atan2(prev[0] * first[1] - prev[1] * first[0], prev[0] * first[0] + prev[1] * first[1]);
        if (L->nl == L->cl) {
            const int c2 = 2 * L->cl + 8;
            int *o2 = realloc(L->off, sizeof *o2 * (size_t)(c2 + 1));
            if (o2) L->off = o2;
            double *a2 = realloc(L->area, sizeof *a2 * (size_t)c2);
            if (a2) L->area = a2;
            double *w2 = realloc(L->wind, sizeof *w2 * (size_t)c2);
            if (w2) L->wind = w2;
            double *r2 = realloc(L->r, sizeof *r2 * (size_t)c2);
            if (r2) L->r = r2;
            if (!o2 || !a2 || !w2 || !r2) return -1;
            L->cl = c2;
        }
        L->area[L->nl] = area, L->wind[L->nl] = wind, L->r[L->nl] = rr;
        L->ntri = base + m, L->off[++L->nl] = L->ntri;
    }
    return 0;
}
/* is the point x (2-D) inside loop l's polygon? (even-odd) */
static inline int mk_cut_inpoly_(const mk_cut_lps_ *L, int l, const double *x)
{
    const double (*bp)[2] = (const double (*)[2])(L->pt + L->off[l]);
    const int nbp = L->off[l + 1] - L->off[l];
    int in = 0;
    for (int i = 0, j = nbp - 1; i < nbp; j = i++)
        if ((bp[i][1] > x[1]) != (bp[j][1] > x[1]) && x[0] < (bp[j][0] - bp[i][0]) * (x[1] - bp[i][1]) / (bp[j][1] - bp[i][1]) + bp[i][0]) in = !in;
    return in;
}
/* of L's loops, the one of least area that winds round p (-1: none), and how many others start inside it */
static inline int mk_cut_pick_(const mk_cut_lps_ *L, int *nested)
{
    int best = -1;
    double best_area = INFINITY;
    for (int l = 0; l < L->nl; l++) {
#if !(defined(MK_MUTATE) && MK_MUTATE == 428)
        if (fabs(L->wind[l]) < 3.14159265358979323846) continue; /* p is not inside it (428: any loop) */
#endif
#if defined(MK_MUTATE) && MK_MUTATE == 427
        if (best < 0 || fabs(L->area[l]) > best_area) best_area = fabs(L->area[l]), best = l; /* 427: the largest */
#else
        if (fabs(L->area[l]) < best_area) best_area = fabs(L->area[l]), best = l;
#endif
    }
    *nested = 0;
    for (int l = 0; best >= 0 && l < L->nl; l++)
        if (l != best) *nested += mk_cut_inpoly_(L, best, L->pt[L->off[l]]);
#if defined(MK_MUTATE) && MK_MUTATE == 432
    *nested = 0; /* 432: a nested loop not refused */
#endif
    return best;
}

/* A probe for callers that try many planes on one mesh: the neighbours once, then mk_cut_probe_loop per plane: 0 if a
   loop winds round p (its farthest point *r, its area *area, its steps *m, the other loops inside it *nested), 1 if none,
   -1 for no memory; mk_cut_probe_inside after it: is the point q (on that plane) inside that loop? */
typedef struct {
    const double *V; const int *T; int nv, nt, stamp; double snap;
    int *nb, *seen; double *d;
    mk_cut_lps_ L; int pick;
    double p[3];
} mk_cut_probe;
static inline void mk_cut_probe_free(mk_cut_probe *pr) { free(pr->nb), free(pr->seen), free(pr->d), mk_cut_lps_free_(&pr->L); memset(pr, 0, sizeof *pr); }
/* 0, -1 (bad argument or no memory), MK_CUT_OPEN */
static inline int mk_cut_probe_init(mk_cut_probe *pr, const double *V, int nv, const int *T, int nt, const mk_cut_opts *o)
{
    memset(pr, 0, sizeof *pr);
    if (!V || !T || !o || nv < 3 || nt < 4 || nt > INT_MAX / 3) return -1;
    pr->V = V, pr->T = T, pr->nv = nv, pr->nt = nt, pr->snap = o->snap, pr->pick = -1;
    pr->nb = malloc(sizeof(int) * 3 * (size_t)nt), pr->seen = calloc((size_t)nt, sizeof(int)), pr->d = malloc(sizeof(double) * (size_t)nv);
    if (!pr->nb || !pr->seen || !pr->d) { mk_cut_probe_free(pr); return -1; }
    const int rc = mk_cut_nbrs_(T, nt, nv, pr->nb);
    if (rc) { mk_cut_probe_free(pr); return rc == 2 ? MK_CUT_OPEN : -1; }
    return 0;
}
static inline int mk_cut_probe_loop(mk_cut_probe *pr, const double *p, const double *n, double *r, double *area, int *m, int *nested)
{
    *r = 0, *area = 0, *m = 0, *nested = 0, pr->pick = -1;
    memcpy(pr->p, p, sizeof pr->p);
    if (mk_cut_trace_(pr->V, pr->nv, pr->T, pr->nt, pr->nb, p, n, pr->snap, pr->d, pr->seen, &pr->stamp, &pr->L)) return -1;
    const int b = mk_cut_pick_(&pr->L, nested);
    if (b < 0) return 1;
    pr->pick = b, *r = pr->L.r[b], *area = fabs(pr->L.area[b]), *m = pr->L.off[b + 1] - pr->L.off[b];
    return 0;
}
static inline int mk_cut_probe_inside(const mk_cut_probe *pr, const double *q)
{
    if (pr->pick < 0) return 0;
    const double r[3] = {q[0] - pr->p[0], q[1] - pr->p[1], q[2] - pr->p[2]}, x[2] = {mk_cut_dot_(r, pr->L.u), mk_cut_dot_(r, pr->L.w)};
    return mk_cut_inpoly_(&pr->L, pr->pick, x);
}

/* Both calls: the cuts' loops chosen (every = 0: mk_cut_loops, 1: mk_cut_clip), then the split, the regions and the caps. */
/* six times a triangle's signed volume against the origin */
static inline double mk_cut_tvol_(const double *V, const int *t)
{
    const double *a = V + 3 * (size_t)t[0], *b = V + 3 * (size_t)t[1], *c = V + 3 * (size_t)t[2];
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}
static inline int mk_cut_piece(const mk_cut_out *o, int g, int **Tg, int *ntg);
static inline int mk_cut_run_(const double *V, int nv, const int *T, int nt, int ncut, const double *P, const double *N, const double *rmax, const mk_cut_opts *opt, int every, mk_cut_out *o)
{
    memset(o, 0, sizeof *o);
    o->bad = -1;
    int rc = -1, nloop = 0, cloop = 0, nb_blk = 0, nloopv = 0;
    int *nb = NULL, *claim = NULL, *seen = NULL, *lv = NULL, *blkt = NULL, *loopbase = NULL, *nlt = NULL, *mark = NULL, **lt = NULL, stamp = 0;
    double *d = NULL, **dk = NULL;
    mk_cut_he_ *blk = NULL;
    mk_cut_lps_ L = {0};
#define MK_CUT_FAIL_(code, cut, ...) do { rc = (code), o->bad = (cut); snprintf(o->reason, sizeof o->reason, __VA_ARGS__); goto out; } while (0)
    if (!V || !T || !P || !N || !opt || (!every && !rmax) || ncut < 1 || nv < 3 || nt < 4 || nt > INT_MAX / 3 || !(opt->snap >= 0) || !(opt->shift >= 0) || !(opt->step >= 0))
        MK_CUT_FAIL_(-1, -1, "a bad argument");
    for (size_t i = 0; i < 3 * (size_t)nt; i++)
        if (T[i] < 0 || T[i] >= nv) MK_CUT_FAIL_(-1, -1, "a bad argument: triangle %zu names vertex %d of %d", i / 3, T[i], nv);
#if !(defined(MK_MUTATE) && MK_MUTATE == 440)
    for (size_t i = 0; i < 3 * (size_t)nv; i++)
        if (!isfinite(V[i])) MK_CUT_FAIL_(-1, -1, "a bad argument: vertex %zu not finite", i / 3); /* (440: accepted) */
#endif
    for (int k = 0; k < ncut; k++) {
        const double *n = N + 3 * k;
        if (!isfinite(P[3 * k]) || !isfinite(P[3 * k + 1]) || !isfinite(P[3 * k + 2]) || !isfinite(n[0]) || !isfinite(n[1]) || !isfinite(n[2]) || !(mk_cut_dot_(n, n) > 0) || (!every && !(rmax[k] > 0)))
            MK_CUT_FAIL_(-1, k, "a bad argument: cut %d's point, normal or rmax", k);
    }
    nb = malloc(sizeof *nb * 3 * (size_t)nt), claim = malloc(sizeof *claim * (size_t)nt), seen = calloc((size_t)nt, sizeof *seen);
    d = malloc(sizeof *d * (size_t)nv), dk = calloc((size_t)ncut, sizeof *dk);
    o->Pc = malloc(sizeof *o->Pc * 3 * (size_t)ncut);
    if (!nb || !claim || !seen || !d || !dk || !o->Pc) MK_CUT_FAIL_(-1, -1, "out of memory");
    rc = mk_cut_nbrs_(T, nt, nv, nb);
    if (rc == 2) MK_CUT_FAIL_(MK_CUT_OPEN, -1, "not a closed manifold: an edge not in two triangles");
    if (rc) MK_CUT_FAIL_(-1, -1, "out of memory");
    for (int t = 0; t < nt; t++) claim[t] = -1;
    for (int k = 0; k < ncut; k++) {
        double nz0[3], u0[3], w0[3];
        mk_cut_basis_(N + 3 * k, nz0, u0, w0);
        /* the loops this cut takes, as indices into L, and L's plane point */
        int pick[1] = {-1}, *take = pick, ntake = 0;
        if (every) {
            memcpy(o->Pc + 3 * k, P + 3 * k, 3 * sizeof *P);
            if (mk_cut_trace_(V, nv, T, nt, nb, P + 3 * k, N + 3 * k, opt->snap, d, seen, &stamp, &L)) MK_CUT_FAIL_(-1, -1, "out of memory");
            if (!L.nl) MK_CUT_FAIL_(MK_CUT_NOLOOP, k, "cut %d: its plane meets the mesh in no loop", k);
            for (int l = 0; l < L.nl; l++) {
                int nest = 0;
                for (int j = 0; j < L.nl; j++)
                    if (j != l) nest += mk_cut_inpoly_(&L, l, L.pt[L.off[j]]);
#if defined(MK_MUTATE) && MK_MUTATE == 432
                nest = 0;
#endif
                if (nest) MK_CUT_FAIL_(MK_CUT_NESTED, k, "cut %d: %d other loop(s) of its plane lie inside its loop %d: its cap would run through another tube", k, nest, l);
            }
            ntake = L.nl;
#if defined(MK_MUTATE) && MK_MUTATE == 436
            ntake = 1; /* 436: only the plane's first loop */
#endif
            take = NULL;
        } else {
            /* the plane moved along N (0, +step, -step, +2 step, ...) to the first place where a loop winds round its point,
               reaches no farther than rmax and holds no other loop */
            int f_found = 0, f_nested = 0;
            double f_r = 0;
            /* whole steps up to shift, never past it (a hair over for rounding: 0.1 / 0.01 gives 10) */
            const int jmax = opt->step > 0 ? 2 * (int)fmin(floor(opt->shift / opt->step * (1 + 1e-9)), 1e6) : 0;
#if defined(MK_MUTATE) && MK_MUTATE == 439
            const int jlim = 0 * jmax; /* 439: the plane never moved */
#else
            const int jlim = jmax;
#endif
            for (int j = 0; j <= jlim && pick[0] < 0; j++) {
                const double del = (j & 1 ? 1 : -1) * ((j + 1) / 2) * opt->step;
                const double pq[3] = {P[3 * k] + del * nz0[0], P[3 * k + 1] + del * nz0[1], P[3 * k + 2] + del * nz0[2]};
                if (mk_cut_trace_(V, nv, T, nt, nb, pq, N + 3 * k, opt->snap, d, seen, &stamp, &L)) MK_CUT_FAIL_(-1, -1, "out of memory");
                int nst = 0;
                const int b = mk_cut_pick_(&L, &nst);
                if (j == 0) f_found = b >= 0, f_nested = nst, f_r = b >= 0 ? L.r[b] : 0;
#if defined(MK_MUTATE) && MK_MUTATE == 433
                if (b < 0 || nst) continue; /* 433: rmax not enforced */
#else
                if (b < 0 || nst || L.r[b] > rmax[k]) continue;
#endif
                pick[0] = b, ntake = 1;
                memcpy(o->Pc + 3 * k, pq, sizeof pq);
            }
            if (pick[0] < 0) {
                if (!f_found) MK_CUT_FAIL_(MK_CUT_NOLOOP, k, "cut %d: no loop of its plane winds round its point", k);
                if (f_nested) MK_CUT_FAIL_(MK_CUT_NESTED, k, "cut %d: %d other loop(s) of its plane lie inside its loop: its cap would run through another tube", k, f_nested);
                MK_CUT_FAIL_(MK_CUT_FAR, k, "cut %d: its loop reaches %.6g from its point, over its rmax %.6g", k, f_r, rmax[k]);
            }
        }
        dk[k] = malloc(sizeof *dk[k] * (size_t)nv);
        if (!dk[k]) MK_CUT_FAIL_(-1, -1, "out of memory");
        if (mk_cut_dist_(V, nv, T, nt, o->Pc + 3 * k, nz0, opt->snap, dk[k])) MK_CUT_FAIL_(-1, -1, "out of memory"); /* the chosen plane's distances, for the split */
        for (int i = 0; i < ntake; i++) {
            const int l = take ? take[i] : i, m = L.off[l + 1] - L.off[l];
            if (nloop == cloop) {
                cloop = 2 * cloop + 8;
                int **lt2 = realloc(lt, sizeof *lt2 * (size_t)cloop);
                if (lt2) lt = lt2;
                int *n2 = realloc(nlt, sizeof *n2 * (size_t)cloop);
                if (n2) nlt = n2;
                int *k2 = realloc(o->loopk, sizeof *k2 * (size_t)cloop);
                if (k2) o->loopk = k2;
                double *r2 = realloc(o->loopr, sizeof *r2 * (size_t)cloop);
                if (r2) o->loopr = r2;
                double *a2 = realloc(o->looparea, sizeof *a2 * (size_t)cloop);
                if (a2) o->looparea = a2;
                double *g2 = realloc(o->margin, sizeof *g2 * (size_t)cloop);
                if (g2) o->margin = g2;
                if (!lt2 || !n2 || !k2 || !r2 || !a2 || !g2) MK_CUT_FAIL_(-1, -1, "out of memory");
            }
            lt[nloop] = malloc(sizeof **lt * (size_t)m);
            if (!lt[nloop]) MK_CUT_FAIL_(-1, -1, "out of memory");
            memcpy(lt[nloop], L.tri + L.off[l], sizeof **lt * (size_t)m);
            nlt[nloop] = m, o->loopk[nloop] = k, o->loopr[nloop] = L.r[l], o->looparea[nloop] = fabs(L.area[l]);
            /* the least fraction along a crossed edge */
            double mg = 1;
            for (int s = 0; s < m && mg > 0; s++) {
                const int t = lt[nloop][s], tn = lt[nloop][(s + 1) % m];
                int e = 0;
                while (e < 3 && nb[3 * t + e] != tn) e++;
                if (e == 3) { mg = 0; break; }
                const int x = T[3 * t + e], y = T[3 * t + (e + 1) % 3];
                const double tt = dk[k][x] / (dk[k][x] - dk[k][y]), g = tt < 1 - tt ? tt : 1 - tt;
                if (g < mg) mg = g > 0 ? g : 0; /* (a snapped end: -0 or 0) */
            }
            o->margin[nloop] = mg;
            nloop++;
            for (int s = 0; s < m; s++) {
                const int t = lt[nloop - 1][s];
#if !(defined(MK_MUTATE) && MK_MUTATE == 434)
                if (claim[t] >= 0 && o->loopk[claim[t]] != k) MK_CUT_FAIL_(MK_CUT_CROSS, k, "cut %d: its loop crosses cut %d's (a triangle on both)", k, o->loopk[claim[t]]);
#endif
                claim[t] = nloop - 1;
            }
        }
    }
    o->nloop = nloop;
    mk_cut_lps_free_(&L);
    free(seen), seen = NULL;
    /* the split: one vertex per crossed edge, each crossed triangle into three */
    for (int l = 0; l < nloop; l++) nloopv += nlt[l];
    if ((size_t)nv + (size_t)nloopv > INT_MAX / 3 || (size_t)nt + MK_CUT_PER_ * (size_t)nloopv > INT_MAX / 3) MK_CUT_FAIL_(-1, -1, "a bad argument: the result too large");
    o->V = malloc(sizeof *o->V * 3 * (size_t)(nv + nloopv)), o->T = malloc(sizeof *o->T * 3 * (size_t)(nt + MK_CUT_PER_ * nloopv));
    o->src = malloc(sizeof *o->src * 2 * (size_t)(nloopv ? nloopv : 1)), o->frac = malloc(sizeof *o->frac * (size_t)(nloopv ? nloopv : 1));
    o->snapped = malloc(sizeof *o->snapped * (size_t)(nloopv ? nloopv : 1)), mark = calloc((size_t)nv, sizeof *mark);
    lv = malloc(sizeof *lv * (size_t)(nloopv ? nloopv : 1));                 /* per loop step: the vertex on its exit edge */
    blk = malloc(sizeof *blk * (size_t)(nloopv ? nloopv : 1));               /* the loops' edges, each with its loop */
    blkt = malloc(sizeof *blkt * 2 * (size_t)(nloopv ? nloopv : 1));
    loopbase = malloc(sizeof *loopbase * (size_t)(nloop + 1));
    o->loopn = malloc(sizeof *o->loopn * (size_t)(nloop ? nloop : 1)), o->side = malloc(sizeof *o->side * 2 * (size_t)(nloop ? nloop : 1));
    if (!o->V || !o->T || !o->src || !o->frac || !o->snapped || !mark || !lv || !blk || !blkt || !loopbase || !o->loopn || !o->side) MK_CUT_FAIL_(-1, -1, "out of memory");
    memcpy(o->V, V, sizeof *V * 3 * (size_t)nv);
    o->nv = o->nv0 = nv;
    loopbase[0] = 0;
    for (int l = 0; l < nloop; l++) loopbase[l + 1] = loopbase[l] + nlt[l], o->loopn[l] = nlt[l];
    /* each loop's vertices: step i's exit edge is the edge it shares with step i + 1 */
    for (int l = 0; l < nloop; l++) {
        const int k = o->loopk[l];
        const double *dd = dk[k];
        for (int i = 0; i < nlt[l]; i++) {
            const int t = lt[l][i], tn = lt[l][(i + 1) % nlt[l]];
            int e = 0;
            while (e < 3 && nb[3 * t + e] != tn) e++;
            if (e == 3) MK_CUT_FAIL_(MK_CUT_SPLIT, k, "cut %d: loop step %d is not joined to the next", k, i);
            const int x = T[3 * t + e], y = T[3 * t + (e + 1) % 3];
            const int sv = mk_cut_on_(dd[x]) ? x : mk_cut_on_(dd[y]) ? y : -1;
            if (sv >= 0) {
                const double *n = N + 3 * k, nn = sqrt(mk_cut_dot_(n, n));
                const double r[3] = {V[3 * sv] - o->Pc[3 * k], V[3 * sv + 1] - o->Pc[3 * k + 1], V[3 * sv + 2] - o->Pc[3 * k + 2]}, h = mk_cut_dot_(r, n) / nn;
#if !(defined(MK_MUTATE) && MK_MUTATE == 430)
                for (int j = 0; j < 3; j++) o->V[3 * sv + j] = V[3 * sv + j] - h * n[j] / nn; /* (430: left where it is) */
#else
                (void)h;
#endif
                if (mark[sv] && mark[sv] != k + 1)   /* two cuts' loops through one snapped vertex: they meet there */
                    MK_CUT_FAIL_(MK_CUT_CROSS, k, "cut %d: its loop meets cut %d's at vertex %d", k, mark[sv] - 1, sv);
                if (!mark[sv]) mark[sv] = k + 1, o->snapped[o->nsnapped++] = sv;
                lv[loopbase[l] + i] = sv;
                continue;
            }
            const double tt = mk_cut_t_(dd[x], dd[y]);
            for (int j = 0; j < 3; j++) o->V[3 * o->nv + j] = V[3 * x + j] + tt * (V[3 * y + j] - V[3 * x + j]);
            const int q = o->nv - nv;
            o->src[2 * q] = x, o->src[2 * q + 1] = y;
#if defined(MK_MUTATE) && MK_MUTATE == 437
            o->frac[q] = 1 - tt; /* 437: the fraction from the other end */
#else
            o->frac[q] = tt;
#endif
            lv[loopbase[l] + i] = o->nv++;
        }
    }
    free(mark), mark = NULL;
    /* the triangles: the uncrossed as they are; a crossed one, its lone corner a (alone on its side), its entry point on
       edge (c, a) and its exit point on edge (a, b), into (a, pab, pca) and the quad (pab, b, c, pca) as two */
    int ot = 0;
    for (int t = 0; t < nt; t++)
        if (claim[t] < 0) memcpy(o->T + 3 * ot, T + 3 * t, 3 * sizeof *T), ot++;
    for (int l = 0; l < nloop; l++) {
        const int k = o->loopk[l];
        for (int i = 0; i < nlt[l]; i++) {
            const int t = lt[l][i], tp = lt[l][(i + nlt[l] - 1) % nlt[l]], tn = lt[l][(i + 1) % nlt[l]];
            const int vin = lv[loopbase[l] + (i + nlt[l] - 1) % nlt[l]], vout = lv[loopbase[l] + i];
            int ein = 0, eout = 0;
            while (ein < 3 && nb[3 * t + ein] != tp) ein++;
            while (eout < 3 && nb[3 * t + eout] != tn) eout++;
            /* the lone corner is the one the two crossed edges share: edges ein (ein -> ein + 1) and eout */
            int a = -1;
            for (int j = 0; j < 3; j++)
                if ((j == ein || j == (ein + 1) % 3) && (j == eout || j == (eout + 1) % 3)) a = j;
            if (a < 0 || ein == eout) MK_CUT_FAIL_(MK_CUT_SPLIT, k, "cut %d: a loop triangle with one crossed edge", k);
            const int A = T[3 * t + a], B = T[3 * t + (a + 1) % 3], C = T[3 * t + (a + 2) % 3];
            const int pab = eout == a ? vout : vin, pca = eout == a ? vin : vout; /* the cut points on edges (A, B) and (C, A) */
            const int tri[3][3] = {{A, pab, pca}, {pab, B, C}, {pab, C, pca}};
            for (int j = 0; j < 3; j++) {
                if (tri[j][0] == tri[j][1] || tri[j][1] == tri[j][2] || tri[j][2] == tri[j][0]) continue; /* (a snapped corner) */
                memcpy(o->T + 3 * ot, tri[j], sizeof tri[j]), ot++;
            }
            if (pab != pca) blk[nb_blk++] = (mk_cut_he_){pab < pca ? (int64_t)pab * (nv + nloopv) + pca : (int64_t)pca * (nv + nloopv) + pab, l};
        }
    }
    o->nt = ot;
    /* the regions: a union over every shared edge but the loops' */
    qsort(blk, (size_t)nb_blk, sizeof *blk, mk_cut_hecmp_);
    for (int i = 0; i < 2 * nb_blk; i++) blkt[i] = -1;
    {
        mk_cut_he_ *h = malloc(sizeof *h * 3 * (size_t)(ot ? ot : 1));
        int *par = malloc(sizeof *par * (size_t)(ot ? ot : 1));
        if (!h || !par) { free(h), free(par); MK_CUT_FAIL_(-1, -1, "out of memory"); }
        const int64_t NV = nv + nloopv;
        for (int t = 0; t < ot; t++) {
            par[t] = t;
            for (int e = 0; e < 3; e++) { const int64_t a = o->T[3 * t + e], b = o->T[3 * t + (e + 1) % 3]; h[3 * t + e] = (mk_cut_he_){(a < b ? a : b) * NV + (a < b ? b : a), t}; }
        }
        qsort(h, 3 * (size_t)ot, sizeof *h, mk_cut_hecmp_);
        int bad = 0;
        for (size_t i = 0; i < 3 * (size_t)ot;) {
            size_t j = i;
            while (j < 3 * (size_t)ot && h[j].key == h[i].key) j++;
            if (j - i != 2) bad = 1;
            else {
                int lo = 0, hi = nb_blk - 1, isblk = 0;
                while (lo <= hi) {
                    const int mid = (lo + hi) / 2;
                    if (blk[mid].key == h[i].key) { isblk = 1, blkt[2 * mid] = h[i].t, blkt[2 * mid + 1] = h[i + 1].t; break; }
                    if (blk[mid].key < h[i].key) lo = mid + 1; else hi = mid - 1;
                }
                if (!isblk) { const int ra = mk_cut_root_(par, h[i].t), rb = mk_cut_root_(par, h[i + 1].t); if (ra != rb) par[ra] = rb; }
            }
            i = j;
        }
        free(h);
        if (bad) { free(par); MK_CUT_FAIL_(MK_CUT_SPLIT, -1, "the split surface has an edge not in two triangles"); }
        o->reg = malloc(sizeof *o->reg * (size_t)(ot ? ot : 1));
        int *rid = malloc(sizeof *rid * (size_t)(ot ? ot : 1));
        if (!o->reg || !rid) { free(par), free(rid); MK_CUT_FAIL_(-1, -1, "out of memory"); }
        for (int t = 0; t < ot; t++) rid[t] = -1;
        o->nreg = 0;
        for (int t = 0; t < ot; t++) { const int r = mk_cut_root_(par, t); if (rid[r] < 0) rid[r] = o->nreg++; o->reg[t] = rid[r]; }
        free(rid), free(par);
        /* each loop's two sides from the triangles on its edges: a triangle's third corner off the plane tells its side (the
           split's own pieces can all be gone where the loop runs through snapped vertices) */
        for (int s = 0; s < 2 * nloop; s++) o->side[s] = -1;
        int clash = 0;
        for (int m2 = 0; m2 < nb_blk; m2++) {
            const int l = blk[m2].t, k = o->loopk[l];
            const int64_t ea = blk[m2].key / NV, eb = blk[m2].key % NV;
            for (int q = 0; q < 2; q++) {
                const int t = blkt[2 * m2 + q];
                if (t < 0) continue;
                int x = -1;
                for (int j = 0; j < 3; j++) if (o->T[3 * t + j] != ea && o->T[3 * t + j] != eb) x = o->T[3 * t + j];
                if (x < 0 || x >= nv || mk_cut_on_(dk[k][x])) continue;
                const int sd = 2 * l + (dk[k][x] > 0);
                if (o->side[sd] < 0) o->side[sd] = o->reg[t]; else if (o->side[sd] != o->reg[t]) clash = 1;
            }
        }
        for (int s = 0; s < 2 * nloop; s++)
            if (o->side[s] < 0) MK_CUT_FAIL_(MK_CUT_NOLOOP, o->loopk[s / 2], "cut %d: its plane only touches the mesh (a face in the plane on its %s side)", o->loopk[s / 2], s % 2 ? "+" : "-");
        if (clash) MK_CUT_FAIL_(MK_CUT_SPLIT, -1, "a loop's side holds two regions");
    }
    /* the caps: each loop's vertices, ear clipped in its plane, wound so the - side's region takes it as it is (+N for an
       outward input): its first edge runs the other way in that region's triangle on it */
    o->capT = malloc(sizeof *o->capT * 3 * (size_t)(nloopv ? nloopv : 1)), o->capk = malloc(sizeof *o->capk * (size_t)(nloopv ? nloopv : 1));
    if (!o->capT || !o->capk) MK_CUT_FAIL_(-1, -1, "out of memory");
    for (int l = 0; l < nloop; l++) {
        const int k = o->loopk[l];
        int m = nlt[l];
        double nz[3], u[3], w[3];
        mk_cut_basis_(N + 3 * k, nz, u, w);
        double (*q)[2] = malloc(sizeof *q * (size_t)m);
        int (*tri)[3] = malloc(sizeof *tri * (size_t)m), *ids = malloc(sizeof *ids * (size_t)m);
        if (!q || !tri || !ids) { free(q), free(tri), free(ids); MK_CUT_FAIL_(-1, -1, "out of memory"); }
        double area = 0;
        int mm = 0;
        for (int i = 0; i < m; i++) {
            const int v = lv[loopbase[l] + i];
            if (mm && ids[mm - 1] == v) continue;
            const double r[3] = {o->V[3 * v] - o->Pc[3 * k], o->V[3 * v + 1] - o->Pc[3 * k + 1], o->V[3 * v + 2] - o->Pc[3 * k + 2]};
            q[mm][0] = mk_cut_dot_(r, u), q[mm][1] = mk_cut_dot_(r, w), ids[mm++] = v;
        }
        while (mm > 1 && ids[mm - 1] == ids[0]) mm--;
        m = mm;
        for (int i = 0; i < m; i++) area += q[i][0] * q[(i + 1) % m][1] - q[(i + 1) % m][0] * q[i][1];
        if (area < 0)
            for (int i = 0; i < m / 2; i++) {
                const double t0 = q[i][0], t1 = q[i][1];
                const int ti = ids[i];
                q[i][0] = q[m - 1 - i][0], q[i][1] = q[m - 1 - i][1], ids[i] = ids[m - 1 - i];
                q[m - 1 - i][0] = t0, q[m - 1 - i][1] = t1, ids[m - 1 - i] = ti;
            }
        int pinch = 0;   /* a loop through one vertex twice touches itself there: not a simple polygon */
        int *sorted = m > 0 ? malloc(sizeof *sorted * (size_t)m) : NULL;
        if (m > 0 && !sorted) { free(q), free(tri), free(ids); MK_CUT_FAIL_(-1, -1, "out of memory"); }
        if (m > 0) {
            memcpy(sorted, ids, sizeof *sorted * (size_t)m);
            qsort(sorted, (size_t)m, sizeof *sorted, mk_cut_intcmp_);
            for (int i = 1; i < m; i++) pinch |= sorted[i] == sorted[i - 1];
            free(sorted);
        }
        const int ec = m < 3 || pinch ? 1 : mk_cut_ears_(m, (const double (*)[2])q, tri);
        if (ec) {
            free(q), free(tri), free(ids);
            if (ec < 0) MK_CUT_FAIL_(-1, -1, "out of memory");
            MK_CUT_FAIL_(MK_CUT_CAP, k, "cut %d: its loop %d is not a simple polygon in its plane (ear clipping failed)", k, l);
        }
        int flip = 0;
        {
            const int64_t a = ids[0], b = ids[1], key = (a < b ? a : b) * (int64_t)(nv + nloopv) + (a < b ? b : a);
            int lo = 0, hi = nb_blk - 1;
            while (lo <= hi) {
                const int mid = (lo + hi) / 2;
                if (blk[mid].key == key) {
                    for (int q = 0; q < 2; q++) {
                        const int t = blkt[2 * mid + q];
                        if (t < 0 || o->reg[t] != o->side[2 * l]) continue;
                        for (int j = 0; j < 3; j++) if (o->T[3 * t + j] == a && o->T[3 * t + (j + 1) % 3] == b) flip = 1;
                    }
                    break;
                }
                if (blk[mid].key < key) lo = mid + 1; else hi = mid - 1;
            }
        }
#if defined(MK_MUTATE) && MK_MUTATE == 443
        flip = 0; /* 443: always +N, whatever the input's winding */
#endif
        for (int i = 0; i < m - 2; i++) o->capT[3 * o->ncapt] = ids[tri[i][0]], o->capT[3 * o->ncapt + 1] = ids[tri[i][flip ? 2 : 1]], o->capT[3 * o->ncapt + 2] = ids[tri[i][flip ? 1 : 2]], o->capk[o->ncapt++] = l;
        free(q), free(tri), free(ids);
    }
#if !(defined(MK_MUTATE) && MK_MUTATE == 435)
    if (!every && o->nreg != ncut + 1) MK_CUT_FAIL_(MK_CUT_REGIONS, -1, "%d regions for %d cuts (a loop does not separate, or the mesh is not one piece)", o->nreg, ncut);
#endif
#if !(defined(MK_MUTATE) && MK_MUTATE == 445)
    { /* each piece closed, consistently wound and wound as the input (positive volume for an outward input): a snap large
         against the edges near a cut can leave a piece open, and a point within shift of the surface a piece inside out */
        double whole = 0;
        for (int t = 0; t < nt; t++) whole += mk_cut_tvol_(V, T + 3 * (size_t)t);
        for (int g = 0; g < o->nreg; g++) {
            int *Tg = NULL, ntg = 0;
            if (mk_cut_piece(o, g, &Tg, &ntg)) MK_CUT_FAIL_(-1, -1, "no memory");
            int *pnb = malloc(sizeof *pnb * 3 * (size_t)(ntg ? ntg : 1));
            const int pc = pnb ? mk_cut_nbrs_(Tg, ntg, o->nv, pnb) : 1;
            double vg = 0;
            for (int t = 0; t < ntg; t++) vg += mk_cut_tvol_(o->V, Tg + 3 * (size_t)t);
            free(pnb), free(Tg);
            if (pc == 1) MK_CUT_FAIL_(-1, -1, "no memory");
            if (pc) MK_CUT_FAIL_(MK_CUT_SPLIT, -1, "piece %d is not closed (a snap of %g against short edges near a cut: lower it)", g, opt->snap);
            if (vg * whole <= 0) MK_CUT_FAIL_(MK_CUT_SPLIT, -1, "piece %d is inside out (volume %g against the whole's %g: a cut's point at the surface)", g, vg, whole);
        }
    }
#endif
    rc = 0;
out:
#undef MK_CUT_FAIL_
    for (int l = 0; lt && l < nloop; l++) free(lt[l]);
    for (int k = 0; dk && k < ncut; k++) free(dk[k]);
    free(lt), free(nlt), free(dk), free(nb), free(claim), free(seen), free(d), free(lv), free(blk), free(blkt), free(loopbase), free(mark);
    mk_cut_lps_free_(&L);
    o->code = rc;
    if (rc) {
        char why[sizeof o->reason];
        const int bad = o->bad;
        memcpy(why, o->reason, sizeof why);
        mk_cut_free(o);
        memcpy(o->reason, why, sizeof why), o->code = rc, o->bad = bad;
    }
    return rc;
}
/* mk_cut_loops: ncut cuts at P[3 k] with normals N[3 k] (any length), loop k within rmax[k] of its point; options o (from
   mk_cut_defaults). Returns 0 with out filled (free it with mk_cut_free), or a refusal or -1 (the top comment). */
static inline int mk_cut_loops(const double *V, int nv, const int *T, int nt, int ncut, const double *P, const double *N, const double *rmax, const mk_cut_opts *o, mk_cut_out *out)
{
    return mk_cut_run_(V, nv, T, nt, ncut, P, N, rmax, o, 0, out);
}
/* mk_cut_clip: every loop of each of the nplane planes (P[3 k], N[3 k]); out->nloop loops, out->loopk[l] each one's plane. */
static inline int mk_cut_clip(const double *V, int nv, const int *T, int nt, int nplane, const double *P, const double *N, const mk_cut_opts *o, mk_cut_out *out)
{
    return mk_cut_run_(V, nv, T, nt, nplane, P, N, NULL, o, 1, out);
}
/* Region g as a closed mesh: its triangles of out->T in their order, then its caps (on a loop's - side as they are, on its
   + side reversed), indices into out->V. *Tg is malloc'd (free it), *ntg its triangles. 0, or -1 (g out of range, no
   memory). */
static inline int mk_cut_piece(const mk_cut_out *o, int g, int **Tg, int *ntg)
{
    *Tg = NULL, *ntg = 0;
    if (!o || g < 0 || g >= o->nreg) return -1;
    int n = 0;
    for (int t = 0; t < o->nt; t++) n += o->reg[t] == g;
    for (int t = 0; t < o->ncapt; t++) n += o->side[2 * o->capk[t]] == g || o->side[2 * o->capk[t] + 1] == g;
    int *R = malloc(sizeof *R * 3 * (size_t)(n ? n : 1)), m = 0;
    if (!R) return -1;
    for (int t = 0; t < o->nt; t++)
        if (o->reg[t] == g) memcpy(R + 3 * m++, o->T + 3 * t, 3 * sizeof *R);
    for (int t = 0; t < o->ncapt; t++) {
        const int l = o->capk[t], *c = o->capT + 3 * t;
        if (o->side[2 * l] == g) R[3 * m] = c[0], R[3 * m + 1] = c[1], R[3 * m + 2] = c[2], m++;
#if defined(MK_MUTATE) && MK_MUTATE == 431
        else if (o->side[2 * l + 1] == g) R[3 * m] = c[0], R[3 * m + 1] = c[1], R[3 * m + 2] = c[2], m++; /* 431: not reversed */
#else
        else if (o->side[2 * l + 1] == g) R[3 * m] = c[0], R[3 * m + 1] = c[2], R[3 * m + 2] = c[1], m++;
#endif
    }
    *Tg = R, *ntg = m;
    return 0;
}
#endif
