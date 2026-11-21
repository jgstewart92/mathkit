/* mathkit: parts of closed outward triangle meshes that sink into each other, moved apart by volume, the mesh kept as
   it is (the same vertices in the same order, the same faces). One call, mk_sepvol, with a report. Header only, C23.
   Two modes. Single (HV NULL): one mesh M whose parts overlap (toes pressed together, two balls sunk into each other).
   A vertex is inside another part when a point just inside it is covered at least twice: the winding number there over
   1 (on a closed mesh the clean rays' facings summed, on any other mk_winding_number). Parity does not decide it: a
   part thinner than that point's depth gives even parity though covered by nothing, and a point covered three times
   gives odd parity. Both sides give way, each by half. Held (HV given): only M moves, by the whole depth; the held mesh
   H is read, never written (a bone held, a muscle or film moved out of it).
   Each round has four steps.
   1. Every vertex inside another part moves toward the nearest point of that part's surface (mk_bvh_nearest_lab):
      edge-aware. At a 90 degree edge the averaged vertex normal points along the diagonal, and a push along it leaves
      crossing pairs; the nearest point is on the face the vertex is actually behind. In single mode a vertex's own
      surface is left out of its search: its patch of inside vertices and MK_SEPVOL_OWN rings around it. In held mode a
      vertex whose nearest exit lies ahead along its own outward normal has gone past the middle of a part thinner than
      M there (a blade, a membrane): step 3 would take the push's outward part away and leave it sliding sideways, round
      after round. It leaves by the side it came in instead, along its inward normal to where that ray leaves H, unless
      M's own surface comes first on that ray (H runs through M there, and going back would push it through M): then its
      nearest exit stays.
   2. Where no vertex is inside, faces can still cross (an edge through a face between two outside corners) or sit
      nearer than gap. For each such pair (mk_bvh_box finds them), a corner behind the other face's plane, or within
      gap/2 of it, and over it, goes to gap in front along that face's normal; in a crossing with no corner over the
      other face (an X), every corner behind the other's plane does. Near pairs are faces facing each other (normals
      over 120 degrees apart), so a part's own curved surface is never one. MK_SEPVOL_PAIRS 0 leaves this step out: the
      pairs are counted and reported, never pushed. A held part that runs through M with no vertex of M inside it (a
      ligament through a coarse muscle, a rib through a sheet) cannot be cleared by moving vertices; this step then
      squeezes M round it (a 0.4 mm plate through a cube of 2.5 mm faces takes a third of the cube's volume).
   3. The moves (depth plus gap, halved in single mode; at most step) are smoothed: each pushed vertex's push turns
      toward the mean of its pushed neighbours', MK_SEPVOL_SMOOTH passes, keeping its own length, so a crumpled patch
      moves together instead of tearing. They spread `rings` rings to the vertices that did not move (half their
      neighbours' mean a ring, so a patch moves, not a spike), and no vertex moves outward along its own normal: that
      part is removed, so nothing grows into a third part. Near another part's edge, smoothing turns a push part way
      toward the other face's (up to 34 degrees off its own face's normal on two cubes).
   4. No move may make M cross itself. The round's moves are tried together; a face that would cross one of M's own
      faces it did not cross at the round's start (any of M's faces in held mode; one of its own part's in single mode,
      where faces of two parts crossing is what is being undone) has its moving corners' moves halved, up to
      MK_SEPVOL_HALVE times, then rolled back (refused), until no such pair is left. Pairs that crossed before are left
      as they are: an input that crosses itself is reported (MK_SEPVOL_SELF_CROSSES) and never made worse. The guard
      holds for the doubles handed back; faces it let come within a rounding of each other can cross once the
      coordinates are written with fewer digits, so write them with %.17g. Steps 3 and 4 work together: on a surface
      crumpled at the edge's scale, unsmoothed moves of half an edge fold it into crossing pairs, and the guard alone
      then stalls.
   mk_sepvol hands back the state with the fewest inside vertices plus crossing pairs it measured, the later on a tie: a
   push the guard stalls can end worse than a round before it, and then the input or an earlier round is returned.
   Rounds repeat until no vertex is inside, no faces cross and none are near, or `rounds` is spent. A run also ends when
   a round's moves are all refused, after MK_SEPVOL_STALL rounds without a better state, or when a part of a closed M
   has lost more than MK_SEPVOL_CRUSH of its volume while vertices are inside or pairs cross (MK_SEPVOL_CRUSHED; the
   best state before that is handed back). Some crossings no vertex move can clear: two long bars crossed like an X have
   no vertex inside each other, and step 2 only squeezes them. Whenever pairs are left at the end,
   MK_SEPVOL_STILL_CROSSES says so.
   The held mesh must define inside and outside: closed, every edge of two faces running it opposite ways, no crossing
   pairs, positive volume, and each connected part wound outward unless it lies inside another part, where it is a
   cavity (a hollow solid). Otherwise mk_sepvol returns -2 with the report's status and held_ fields saying why.
   No vertex moves more than step in a round, so none more than rounds times step in all. Lengths are in the caller's
   unit (gap and step in it). mk_closest_pt_tri's face weights can be thousands of ulps off on a triangle of aspect 1e5
   (about 1e-12 of the triangle's size), far under any gap, so this cannot change a push.
   Experimental: its results may change between releases. Single-threaded. */
#ifndef MK_SEPVOL_H
#define MK_SEPVOL_H
#include "mk_geom.h"
#ifndef MK_SEPVOL_OWN
#define MK_SEPVOL_OWN 3 /* rings of a vertex's own surface, around its patch, that its search leaves out */
#endif
#ifndef MK_SEPVOL_HALVE
#define MK_SEPVOL_HALVE 4 /* halvings of a move that would make M cross itself, before it is refused */
#endif
#ifndef MK_SEPVOL_STALL
#define MK_SEPVOL_STALL 8 /* rounds with no fewer inside vertices plus crossing pairs before the run ends */
#endif
#ifndef MK_SEPVOL_SMOOTH
#define MK_SEPVOL_SMOOTH 2 /* passes turning each push toward its pushed neighbours' mean direction (step 3) */
#endif
#ifndef MK_SEPVOL_PAIRS
#define MK_SEPVOL_PAIRS 1 /* step 2 on (1) or off (0: crossing and near pairs counted, never pushed). Any int
                             expression, read at each call, so a caller may name its own variable and choose per call */
#endif
#ifndef MK_SEPVOL_CRUSH
#define MK_SEPVOL_CRUSH 0.35 /* the share of a part's volume a run may take while vertices are inside or pairs cross.
                                Any double expression, read at each call; 1 or more turns the stop off */
#endif
/* the report's status bits. The held ones make mk_sepvol return -2 without moving anything. */
#define MK_SEPVOL_HELD_OPEN 1     /* the held mesh has edges not shared by exactly two faces (holes, or 3+ faces) */
#define MK_SEPVOL_HELD_FLIPPED 2  /* its faces disagree on their facing (an edge run one way by both), a part is wound
                                     inward with no part around it, or it encloses no volume */
#define MK_SEPVOL_HELD_CROSSES 4  /* its faces cross each other (overlapping parts, or a part through itself) */
#define MK_SEPVOL_SELF_CROSSES 8  /* M crossed itself at the start (see self0): it ran, and never added to them */
#define MK_SEPVOL_STILL_CROSSES 16 /* faces still cross in the state handed back (rep->cross): parts not separated */
#define MK_SEPVOL_CRUSHED 32      /* stopped: going on would have taken more than MK_SEPVOL_CRUSH of a part */
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int rounds;        /* rounds run */
    long inside0;      /* M's vertices inside another part (held: inside H) before */
    long inside;       /* and after */
    long cross0, cross; /* crossing face pairs (two faces that share no corner, an edge of one through the other),
                           before and after */
    long moved;        /* vertices moved (whose returned position differs from the input's) */
    double worst;      /* the largest total move of a vertex */
    double vol0, vol1; /* M's enclosed volume before and after (signed: positive for an outward mesh) */
    int status;        /* MK_SEPVOL_* bits, 0 when nothing was found */
    long self0, self1; /* M's own crossing pairs before and after (held: any two faces; single: two of one part) */
    long refused;      /* moves rolled back, summed over rounds: MK_SEPVOL_HALVE halvings still made M cross itself */
    long shortened;    /* vertex moves halved at least once and then kept, summed over rounds */
    /* the held mesh: open or shared-by-3+ edges, edges run the same way by both faces, crossing pairs */
    long held_open, held_flipped, held_cross;
    double held_vol;   /* its enclosed volume (signed) */
    long held_inward;  /* the held mesh's parts wound inward with no other part around them (MK_SEPVOL_HELD_FLIPPED) */
    double part_loss;  /* the largest share of a part's volume lost in the state handed back (0 when M is not closed) */
} mk_sepvol_report;

static inline double mk_sepvol_volume_(const double *V, const int *T, int nt)
{
    double s = 0;
    for (int t = 0; t < nt; t++) {
        const double *a = V + 3 * (size_t)T[3 * (size_t)t], *b = V + 3 * (size_t)T[3 * (size_t)t + 1], *c = V + 3 * (size_t)T[3 * (size_t)t + 2];
        s += a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    return s / 6;
}
static inline int mk_sepvol_find_(int *uf, int a)
{
    while (uf[a] != a) a = uf[a] = uf[uf[a]];
    return a;
}
static inline int mk_sepvol_ok_(const double *V, int nv, const int *T, int nt)
{
    if (nv < 1 || nt < 1 || (size_t)nv > SIZE_MAX / 64 || (size_t)nt > SIZE_MAX / 64) return 0;
    for (size_t q = 0; q < 3 * (size_t)nv; q++) if (!isfinite(V[q])) return 0;
    for (size_t q = 0; q < 3 * (size_t)nt; q++) if (T[q] < 0 || T[q] >= nv) return 0;
    return 1;
}
/* p covered by two parts of one mesh, p being just inside a part: the winding number at p is 2 or more. On a closed
   mesh (closed 1: every edge run as often each way) it is the clean rays' facings summed (mk_bvh_count_wind_), when
   they agree, whatever their parity: three parts over p read odd, and even can mean covered by nothing. Otherwise, or
   when no ray is clean, mk_winding_number, O(triangles), over 1.5; except that on a mesh not closed, odd parity on
   every clean ray and no ray's facings over 1 reads covered once, which keeps that cost from every vertex of an open
   mesh. mk_bvh_inside is not used here: its winding fallback reads a point inside two parts as inside. */
static inline int mk_sepvol_two_(const mk_bvh *B, const double *p, int closed)
{
    static const double dirs[3][3] = {{0.5773502691896258, 0.5773502691896258, 0.5773502691896258},
                                      {-0.3015113445777636, 0.9045340337332909, 0.3015113445777636},
                                      {0.2672612419124244, -0.5345224838248488, 0.8017837257372732}};
    int odd = 0, even = 0, nw = 0, w0 = 0, wmax = INT_MIN, same = 1;
    for (int k = 0; k < 3; k++) {
        int edge, w, c = mk_bvh_count_wind_(B, p, dirs[k], &edge, &w);
        if (edge) continue;
        odd += c & 1, even += !(c & 1);
        if (nw++) same &= w == w0;
        else w0 = w;
        if (w > wmax) wmax = w;
    }
#if defined(MK_MUTATE) && MK_MUTATE == 213
    if ((void)closed, (void)same, even && !odd) return 1; /* 213: even read as covered twice */
#endif
#if defined(MK_MUTATE) && MK_MUTATE == 216
    if (odd && !even) return 0; /* 216: odd read as covered once before the facings are looked at */
#endif
    if (closed && nw && same) return w0 > 1;
#if defined(MK_MUTATE) && MK_MUTATE == 217
    if (odd && !even) return 0; /* 217: on an open mesh odd parity is covered once, whatever the facings */
#else
    if (odd && !even && wmax <= 1) return 0;
#endif
    return mk_winding_number(B->V, B->T, B->nt, p) > 1.5;
}
static inline int mk_sepvol_u64cmp_(const void *a, const void *b)
{
    const uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}
/* every edge of the mesh run as often one way as the other (a repeated corner's edge aside): its winding number is then
   the facings summed along any ray with no edge hit. 1, 0, or -1 for no memory */
static inline int mk_sepvol_closed_(const int *T, int nt)
{
    uint64_t *e = malloc(sizeof *e * 3 * (size_t)nt);
    if (!e) return -1;
    size_t ne = 0;
    for (size_t t = 0; t < (size_t)nt; t++)
        for (int k = 0; k < 3; k++) {
            const uint64_t a = (uint64_t)T[3 * t + k], b = (uint64_t)T[3 * t + (k + 1) % 3];
            /* the edge, and which way this face runs it */
            if (a != b) e[ne++] = a < b ? (a << 32 | b) << 1 : (b << 32 | a) << 1 | 1;
        }
    qsort(e, ne, sizeof *e, mk_sepvol_u64cmp_);
    int ok = 1;
    for (size_t i = 0, j; i < ne && ok; i = j) {
        long run = 0;
        for (j = i; j < ne && e[j] >> 1 == e[i] >> 1; j++) run += e[j] & 1 ? -1 : 1;
        ok = run == 0;
    }
    free(e);
#if defined(MK_MUTATE) && MK_MUTATE == 215
    ok = 1; /* 215: every mesh taken as closed */
#endif
    return ok;
}
/* which of M's vertices are inside another part: the count, flags in in[] */
static inline long mk_sepvol_mark_(const double *V, int nv, const double *N, const mk_bvh *B, int held, int closed, double eps, char *in)
{
    long n = 0;
    for (int i = 0; i < nv; i++) {
        const double *v = V + 3 * (size_t)i, *u = N + 3 * (size_t)i, p[3] = {v[0] - eps * u[0], v[1] - eps * u[1], v[2] - eps * u[2]};
        in[i] = held ? mk_bvh_inside(B, v) : mk_sepvol_two_(B, p, closed);
        n += in[i];
    }
    return n;
}
static inline void mk_sepvol_normals_(const double *V, int nv, const int *T, int nt, double *N)
{
    memset(N, 0, sizeof *N * 3 * (size_t)nv);
    for (int t = 0; t < nt; t++) {
        const double *a = V + 3 * (size_t)T[3 * (size_t)t], *b = V + 3 * (size_t)T[3 * (size_t)t + 1], *c = V + 3 * (size_t)T[3 * (size_t)t + 2];
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double n[3] = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]};
        for (int k = 0; k < 3; k++) for (int j = 0; j < 3; j++) N[3 * (size_t)T[3 * (size_t)t + k] + j] += n[j];
    }
    for (int i = 0; i < nv; i++) {
        double *u = N + 3 * (size_t)i, l = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        if (l > 0) u[0] /= l, u[1] /= l, u[2] /= l;
    }
}

static inline void mk_sepvol_fnormal_(const double *V, const int *q, double *n)
{
    const double *a = V + 3 * (size_t)q[0], *b = V + 3 * (size_t)q[1], *c = V + 3 * (size_t)q[2];
    const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    n[0] = u[1] * w[2] - u[2] * w[1], n[1] = u[2] * w[0] - u[0] * w[2], n[2] = u[0] * w[1] - u[1] * w[0];
    const double l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (l > 0) n[0] /= l, n[1] /= l, n[2] /= l;
}
/* an edge of triangle P (in VP) through triangle Q (in VQ) */
static inline int mk_sepvol_edge_through_(const double *VP, const int *P, const double *VQ, const int *Q)
{
    for (int e = 0; e < 3; e++) {
        const double *u = VP + 3 * (size_t)P[e], *v = VP + 3 * (size_t)P[(e + 1) % 3], d[3] = {v[0] - u[0], v[1] - u[1], v[2] - u[2]};
        const mk_ray ry = mk_ray_make(u, d);
        double t, b[3];
        if (mk_ray_tri(&ry, VQ + 3 * (size_t)Q[0], VQ + 3 * (size_t)Q[1], VQ + 3 * (size_t)Q[2], 0, 1, &t, b)) return 1;
    }
    return 0;
}
/* step 2's rule for one pair: each corner of X (in VX) behind Y's plane, or nearer it than gap/2, and over Y (any:
   every corner behind it), is pushed along Y's normal nY to gap in front, by `share` of the way (1/2 when both sides
   give way), with `sign`. acc and cnt gather pushes per vertex. Returns the number of such corners. */
static inline int mk_sepvol_corners_(const double *VX, const int *X, const double *VY, const int *Y, const double *nY, double gap, int any,
                                      double share, double sign, const int *to, double *acc, int *cnt)
{
    int n = 0;
    for (int k = 0; k < 3; k++) {
        const double *c = VX + 3 * (size_t)X[k], *a = VY + 3 * (size_t)Y[0];
        const double s = (c[0] - a[0]) * nY[0] + (c[1] - a[1]) * nY[1] + (c[2] - a[2]) * nY[2];
        if (!(s < 0.5 * gap)) continue;
        double w[3], q[3];
        const double d2 = mk_closest_pt_tri(c, a, VY + 3 * (size_t)Y[1], VY + 3 * (size_t)Y[2], w, q);
        if (any ? !(s < 0) : d2 > s * s * (1 + 1e-9) + 1e-30) continue; /* not over Y: nearest point on an edge */
        n++;
        if (!acc) continue;
        const double m = share * (gap - s);
        for (int j = 0; j < 3; j++) {
            const int v = to ? to[j] : X[k]; /* to: the three corners of the moving face, when X's corner can't move */
            if (to == NULL && j) break;
            for (int i = 0; i < 3; i++) acc[3 * (size_t)v + i] += sign * m * nY[i];
            cnt[v]++;
        }
    }
    return n;
}
/* crossing pairs of M's faces with each other (O M's own tree, VH NULL) or with H's (O H's tree): their count; faces
   that share a corner are not a pair. Faces of opposite facing (normals at over 120 degrees) nearer than gap/2 are a
   near pair: counted in *nnear. With acc, the pushes of each pair are gathered (mk_sepvol_corners_). -1: no memory. */
static inline long mk_sepvol_pairs_(const double *V, const int *T, int nt, const mk_bvh *O, const double *VH, const int *TH, double gap,
                                    double *acc, int *cnt, long *nnear)
{
    int cap = 64, *buf = malloc(sizeof *buf * (size_t)cap);
    long np = 0;
    *nnear = 0;
    if (!buf) return -1;
    const double *VO = VH ? VH : V;
    const int *TO = VH ? TH : T;
    for (int t = 0; t < nt; t++) {
        const int *X = T + 3 * (size_t)t;
        mk_aabb bx;
        mk_tri_box_(V, T, t, &bx);
        for (int k = 0; k < 3; k++) bx.lo[k] -= 0.5 * gap, bx.hi[k] += 0.5 * gap;
        int n = mk_bvh_box(O, &bx, buf, cap);
        if (n > cap) {
            free(buf), cap = n;
            if (!(buf = malloc(sizeof *buf * (size_t)cap))) return -1;
            n = mk_bvh_box(O, &bx, buf, cap);
        }
        double nX[3];
        mk_sepvol_fnormal_(V, X, nX);
        for (int i = 0; i < n; i++) {
            const int u = buf[i], *Y = TO + 3 * (size_t)u;
            if (!VH && u <= t) continue;
            if (!VH && (X[0] == Y[0] || X[0] == Y[1] || X[0] == Y[2] || X[1] == Y[0] || X[1] == Y[1] || X[1] == Y[2] || X[2] == Y[0] ||
                        X[2] == Y[1] || X[2] == Y[2]))
                continue;
            double nY[3];
            mk_sepvol_fnormal_(VO, Y, nY);
            const int cross = mk_sepvol_edge_through_(V, X, VO, Y) || mk_sepvol_edge_through_(VO, Y, V, X);
#if defined(MK_MUTATE) && MK_MUTATE == 155
            if (!cross) continue; /* 155: near pairs ignored */
#endif
            if (!cross && !(nX[0] * nY[0] + nX[1] * nY[1] + nX[2] * nY[2] < -0.5)) continue;
            int pushed = 0;
#if defined(MK_MUTATE) && MK_MUTATE == 154
            const int fall = 0; /* 154: no fallback for an X */
#else
            const int fall = cross;
#endif
            for (int any = 0; any <= fall && !pushed; any++) { /* a crossing pair with no corner over the other face (an
                                                                  X): every corner behind the other's plane */
                if (VH) { /* M's corners out of H's face, and M's face off H's corners (H stays) */
                    pushed = mk_sepvol_corners_(V, X, VO, Y, nY, gap, any, 1, 1, NULL, acc, cnt);
                    pushed += mk_sepvol_corners_(VO, Y, V, X, nX, gap, any, 1, -1, X, acc, cnt);
                } else {
                    pushed = mk_sepvol_corners_(V, X, VO, Y, nY, gap, any, 0.5, 1, NULL, acc, cnt);
                    pushed += mk_sepvol_corners_(VO, Y, V, X, nX, gap, any, 0.5, 1, NULL, acc, cnt);
                }
            }
            np += cross, *nnear += !cross && pushed;
        }
    }
    free(buf);
    return np;
}

/* the box query into a growing buffer; -1: no memory */
static inline int mk_sepvol_box_(const mk_bvh *B, const mk_aabb *q, int **buf, int *cap)
{
    int n = mk_bvh_box(B, q, *buf, *cap);
    if (n > *cap) {
        int *b = realloc(*buf, sizeof *b * (size_t)n);
        if (!b) return -1;
        *buf = b, *cap = n, n = mk_bvh_box(B, q, *buf, *cap);
    }
    return n;
}
/* a mesh's own crossing pairs (mk_tri_cross, as mk_bvh_crossings counts them), each counted once; with part, only pairs
   of faces of one part. B: its tree, over positions each within grow of V's. -1: no memory */
static inline long mk_sepvol_self_(const double *V, const int *T, int nt, const mk_bvh *B, double grow, const int *part)
{
    int cap = 64, *buf = malloc(sizeof *buf * (size_t)cap);
    long np = 0;
    if (!buf) return -1;
    for (int t = 0; t < nt; t++) {
        const int *X = T + 3 * (size_t)t;
        mk_aabb bx;
        mk_tri_box_(V, T, t, &bx);
        for (int k = 0; k < 3; k++) bx.lo[k] -= grow, bx.hi[k] += grow;
        const int n = mk_sepvol_box_(B, &bx, &buf, &cap);
        if (n < 0) { free(buf); return -1; }
        for (int i = 0; i < n; i++) {
            const int u = buf[i], *Y = T + 3 * (size_t)u;
            if (u <= t || (part && part[X[0]] != part[Y[0]])) continue;
            np += mk_tri_cross(V, X, Y);
        }
    }
    free(buf);
    return np;
}
/* the held mesh's connected parts wound inward (negative volume, about the part's own corner, past rounding) that lie
   inside no other part: the winding number of the rest of the mesh at one of the part's vertices is not over 1/2. An
   inward part inside another is a cavity, a hollow solid's inner wall. The count, -1: no memory */
static inline long mk_sepvol_inward_(const double *HV, const int *HT, int hnt)
{
    int nv = 0;
    for (size_t q = 0; q < 3 * (size_t)hnt; q++) if (HT[q] >= nv) nv = HT[q] + 1;
    int *uf = malloc(sizeof *uf * (size_t)nv), *rest = malloc(sizeof *rest * 3 * (size_t)hnt);
    double *vol = calloc((size_t)nv, sizeof *vol), *mag = calloc((size_t)nv, sizeof *mag);
    long bad = -1;
    if (!uf || !rest || !vol || !mag) goto out;
    bad = 0;
    for (int i = 0; i < nv; i++) uf[i] = i;
    for (size_t t = 0; t < (size_t)hnt; t++)
        for (int k = 0; k < 2; k++) uf[mk_sepvol_find_(uf, HT[3 * t + k])] = mk_sepvol_find_(uf, HT[3 * t + k + 1]);
    for (size_t t = 0; t < (size_t)hnt; t++) {
        const int r = mk_sepvol_find_(uf, HT[3 * t]);
        const double *o = HV + 3 * (size_t)r, *a = HV + 3 * (size_t)HT[3 * t], *b = HV + 3 * (size_t)HT[3 * t + 1], *c = HV + 3 * (size_t)HT[3 * t + 2];
        const double u[3] = {a[0] - o[0], a[1] - o[1], a[2] - o[2]}, v[3] = {b[0] - o[0], b[1] - o[1], b[2] - o[2]}, w[3] = {c[0] - o[0], c[1] - o[1], c[2] - o[2]};
        const double d = u[0] * (v[1] * w[2] - v[2] * w[1]) + u[1] * (v[2] * w[0] - v[0] * w[2]) + u[2] * (v[0] * w[1] - v[1] * w[0]);
        vol[r] += d, mag[r] += fabs(d);
    }
    for (int r = 0; r < nv; r++) {
        if (uf[r] != r || !(vol[r] < -1e-12 * mag[r])) continue;
#if defined(MK_MUTATE) && MK_MUTATE == 218
        continue; /* 218: only the total's sign looked at */
#endif
#if defined(MK_MUTATE) && MK_MUTATE == 219
        bad++; /* 219: a cavity refused like a part wound inward on its own */
        continue;
#endif
        int n = 0;
        for (size_t t = 0; t < (size_t)hnt; t++)
            if (mk_sepvol_find_(uf, HT[3 * t]) != r) memcpy(rest + 3 * (size_t)n++, HT + 3 * t, sizeof *rest * 3);
        bad += !(n && mk_winding_number(HV, rest, n, HV + 3 * (size_t)r) > 0.5);
    }
out:
    free(uf), free(rest), free(vol), free(mag);
    return bad;
}
/* can the held mesh say inside from outside: its edges, facing, crossing pairs and volume into r; the status bits found
   (0: usable), -1: no memory */
static inline int mk_sepvol_heldcheck_(const double *HV, const int *HT, int hnt, const mk_bvh *BH, mk_sepvol_report *r)
{
    uint64_t *e = malloc(sizeof *e * 3 * (size_t)hnt);
    if (!e) return -1;
    for (size_t t = 0; t < (size_t)hnt; t++)
        for (int k = 0; k < 3; k++) {
            const uint64_t a = (uint64_t)HT[3 * t + k], b = (uint64_t)HT[3 * t + (k + 1) % 3];
            /* the edge, and which way this face runs it */
            e[3 * t + k] = (a < b ? (a << 32 | b) << 1 : (b << 32 | a) << 1 | 1);
        }
    qsort(e, 3 * (size_t)hnt, sizeof *e, mk_sepvol_u64cmp_);
    for (size_t i = 0, j; i < 3 * (size_t)hnt; i = j) {
        for (j = i + 1; j < 3 * (size_t)hnt && e[j] >> 1 == e[i] >> 1;) j++;
        if (j - i != 2) r->held_open++;
        else if ((e[i] & 1) == (e[i + 1] & 1)) r->held_flipped++;
    }
    free(e);
    r->held_vol = mk_sepvol_volume_(HV, HT, hnt);
    if ((r->held_cross = mk_sepvol_self_(HV, HT, hnt, BH, 0, NULL)) < 0) return -1;
    if ((r->held_inward = mk_sepvol_inward_(HV, HT, hnt)) < 0) return -1;
    int st = 0;
    if (r->held_open) st |= MK_SEPVOL_HELD_OPEN;
    if (r->held_flipped || r->held_inward || !(r->held_vol > 0)) st |= MK_SEPVOL_HELD_FLIPPED;
    if (r->held_cross) st |= MK_SEPVOL_HELD_CROSSES;
#if defined(MK_MUTATE) && MK_MUTATE == 157
    st = 0; /* 157: every held mesh taken */
#endif
    return st;
}
/* step 4: the moves D tried together (s[i] = 1 for each moving vertex, 0 otherwise; V the positions at the round's
   start, B M's tree over positions each within grow of V's). A face that would cross one of M's faces it did not cross
   at the start (with part: only one of its own part's) halves its moving corners' moves, MK_SEPVOL_HALVE times, then
   refuses them (lvl[i] > MK_SEPVOL_HALVE), until no such pair is left. Vn: the positions after. 0, or -1 for no
   memory. */
static inline int mk_sepvol_guard_(const double *V, double *Vn, const double *D, int nv, const int *T, int nt, const mk_bvh *B, double grow,
                                   const int *part, double *s, unsigned char *lvl, char *ch, char *bad)
{
    double maxm = 0;
    for (int i = 0; i < nv; i++) {
        const double *e = D + 3 * (size_t)i;
        lvl[i] = 0, ch[i] = s[i] > 0;
        if (ch[i]) maxm = fmax(maxm, sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]));
    }
    if (nv < 1) return 0;
    memcpy(Vn, V, sizeof *Vn * 3 * (size_t)nv);
    int cap = 64, *buf = malloc(sizeof *buf * (size_t)cap);
    if (!buf) return -1;
    for (;;) {
        for (int i = 0; i < nv; i++)
            if (ch[i]) for (int k = 0; k < 3; k++) Vn[3 * (size_t)i + k] = V[3 * (size_t)i + k] + s[i] * D[3 * (size_t)i + k];
        memset(bad, 0, (size_t)nv);
        long nbad = 0;
        for (int t = 0; t < nt; t++) {
            const int *X = T + 3 * (size_t)t;
            if (!ch[X[0]] && !ch[X[1]] && !ch[X[2]]) continue;
            mk_aabb bx; /* a candidate's box after its move lies within maxm of its box before, which the tree holds */
            mk_tri_box_(Vn, T, t, &bx);
            for (int k = 0; k < 3; k++) bx.lo[k] -= maxm + grow, bx.hi[k] += maxm + grow;
            const int n = mk_sepvol_box_(B, &bx, &buf, &cap);
            if (n < 0) { free(buf); return -1; }
            for (int i = 0; i < n; i++) {
                const int u = buf[i], *Y = T + 3 * (size_t)u;
#if defined(MK_MUTATE) && MK_MUTATE == 158
                if ((void)part, u == t) continue; /* 158: faces of two parts guarded too, in single mode */
#else
                if (u == t || (part && part[X[0]] != part[Y[0]])) continue;
#endif
                if (!mk_tri_cross(Vn, X, Y) || mk_tri_cross(V, X, Y)) continue;
                for (int k = 0; k < 3; k++) {
                    if (s[X[k]] > 0) bad[X[k]] = 1;
                    if (s[Y[k]] > 0) bad[Y[k]] = 1;
                }
                nbad++;
            }
        }
#if defined(MK_MUTATE) && MK_MUTATE == 156
        nbad = 0; /* 156: no guard, every move whole */
#endif
        if (!nbad) break;
        for (int i = 0; i < nv; i++) {
            ch[i] = bad[i];
            if (bad[i]) s[i] = ++lvl[i] > MK_SEPVOL_HALVE ? 0 : 0.5 * s[i];
        }
    }
    free(buf);
    return 0;
}

/* step 3's smoothing: each pushed vertex's push (fix: step 1's and step 2's) turns toward the mean of its own and its
   pushed neighbours' pushes, MK_SEPVOL_SMOOTH times, keeping its own length. acc, na: scratch. */
static inline void mk_sepvol_smooth_(double *D, const char *fix, int nv, const int *T, int nt, double *acc, int *na)
{
    for (int it = 0; it < MK_SEPVOL_SMOOTH && nv > 0; it++) {
        memcpy(acc, D, sizeof *acc * 3 * (size_t)nv);
        for (int i = 0; i < nv; i++) na[i] = 1;
        for (int t = 0; t < nt; t++)
            for (int k = 0; k < 3; k++) {
                const int v = T[3 * (size_t)t + k];
                if (!fix[v]) continue;
                for (int j = 1; j < 3; j++) {
                    const int u = T[3 * (size_t)t + (k + j) % 3];
                    if (!fix[u]) continue;
                    for (int c = 0; c < 3; c++) acc[3 * (size_t)v + c] += D[3 * (size_t)u + c];
                    na[v]++;
                }
            }
        for (int i = 0; i < nv; i++) {
            if (!fix[i] || na[i] < 2) continue;
            double *e = D + 3 * (size_t)i, x[3], l0 = 0, l1 = 0;
            for (int c = 0; c < 3; c++) x[c] = acc[3 * (size_t)i + c] / na[i], l0 += e[c] * e[c], l1 += x[c] * x[c];
            if (l1 > 0) for (int c = 0; c < 3; c++) e[c] = x[c] * sqrt(l0 / l1);
        }
    }
}

/* each part's volume into pv[its root vertex] (part[v]: the root), about the root's own place, so a part far from the
   origin keeps its digits; 0 for a vertex that is no root */
static inline void mk_sepvol_partvol_(const double *V, const int *T, int nt, const int *part, int nv, double *pv)
{
    memset(pv, 0, sizeof *pv * (size_t)nv);
    for (size_t t = 0; t < (size_t)nt; t++) {
        const int r = part[T[3 * t]];
        const double *o = V + 3 * (size_t)r, *a = V + 3 * (size_t)T[3 * t], *b = V + 3 * (size_t)T[3 * t + 1], *c = V + 3 * (size_t)T[3 * t + 2];
        const double u[3] = {a[0] - o[0], a[1] - o[1], a[2] - o[2]}, v[3] = {b[0] - o[0], b[1] - o[1], b[2] - o[2]}, w[3] = {c[0] - o[0], c[1] - o[1], c[2] - o[2]};
        pv[r] += (u[0] * (v[1] * w[2] - v[2] * w[1]) + u[1] * (v[2] * w[0] - v[0] * w[2]) + u[2] * (v[0] * w[1] - v[1] * w[0])) / 6;
    }
}
/* the largest share of a part's volume lost against pv0 */
static inline double mk_sepvol_loss_(const double *pv0, const double *pv, int nv)
{
    double l = 0;
    for (int i = 0; i < nv; i++)
        if (pv0[i] != 0 && 1 - pv[i] / pv0[i] > l) l = 1 - pv[i] / pv0[i];
    return l;
}

/* step 1 in held mode: v, inside H, with its nearest exit dir at distance *d. When that exit lies ahead along v's own
   outward normal u, v has gone past the middle of a part thinner than M there, and the step 3 rule takes the outward
   part of the push away: it would slide sideways round after round. It leaves by the side it came in instead, along -u
   to where that ray leaves H, if it leaves H before it leaves M (MB: M's tree, its positions within tm / 2 of V; tm
   also skips v's own faces). If M's surface comes first, H runs through M there and going back would push v through M
   itself: the nearest exit stays. th: the ray's start. */
static inline void mk_sepvol_entry_(const mk_bvh *H, const mk_bvh *MB, const double *v, const double *u, double th, double tm, double *d, double *dir)
{
#if defined(MK_MUTATE) && MK_MUTATE == 193
    return; /* 193: the nearest exit always */
#endif
    if (!(dir[0] * u[0] + dir[1] * u[1] + dir[2] * u[2] > 0)) return;
    const double w[3] = {-u[0], -u[1], -u[2]};
    double t, b[3];
    if (mk_bvh_first_hit(H, v, w, th, INFINITY, &t, b) < 0) return;
#if !(defined(MK_MUTATE) && MK_MUTATE == 194) /* 194: back through H whatever lies behind */
    double s;
    if (tm < t && mk_bvh_first_hit(MB, v, w, tm, t, &s, b) >= 0) return;
#else
    (void)MB, (void)tm;
#endif
    *d = t, dir[0] = w[0], dir[1] = w[1], dir[2] = w[2];
}

/* Moves M's vertices (V[3 nv], faces T[3 nt]) in place, out of the other parts of M (HV NULL) or out of the held mesh H
   (HV[3 hnv], HT[3 hnt]). gap >= 0, step > 0, rings >= 0 as above, rounds >= 1. Returns 0 when it ran (rep->inside and
   rep->cross say whether anything is left: rounds may run out first, every move be refused, or MK_SEPVOL_CRUSH stop it;
   rep->status says so); -1 for a bad mesh (a size under 1, an index out of range, a coordinate not finite), gap or step
   not finite or out of range, rings < 0, rounds < 1, or no memory (rep zeroed); -2 for a held mesh that cannot say
   inside from outside (rep->status and the held_ fields say why). On -1 and -2, V is as it came in. rep may be NULL. */
static inline int mk_sepvol(double *V, int nv, const int *T, int nt, const double *HV, int hnv, const int *HT, int hnt, double gap,
                            double step, int rings, int rounds, mk_sepvol_report *rep)
{
    mk_sepvol_report r = {0};
    if (rep) *rep = r;
    const int held = HV != NULL;
    if (!mk_sepvol_ok_(V, nv, T, nt) || (held && (!HT || !mk_sepvol_ok_(HV, hnv, HT, hnt))) || !(gap >= 0) || !isfinite(gap) ||
        !(step > 0) || !isfinite(step) || rings < 0 || rounds < 1)
        return -1;
    double *N = malloc(sizeof *N * 3 * (size_t)nv), *D = malloc(sizeof *D * 3 * (size_t)nv), *acc = malloc(sizeof *acc * 3 * (size_t)nv);
    double *V0 = malloc(sizeof *V0 * 3 * (size_t)nv), *Vr = malloc(sizeof *Vr * 3 * (size_t)nv), *Vn = malloc(sizeof *Vn * 3 * (size_t)nv);
    double *Vb = malloc(sizeof *Vb * 3 * (size_t)nv), *Vk = malloc(sizeof *Vk * 3 * (size_t)nv);
    double *sc = malloc(sizeof *sc * (size_t)nv), *pv0 = malloc(sizeof *pv0 * (size_t)nv), *pv = malloc(sizeof *pv * (size_t)nv);
#if !(defined(MK_MUTATE) && MK_MUTATE == 212) /* 212: V0 filled after every allocation is checked */
    /* before any goto out: the restore there copies V0 back */
    if (V0) memcpy(V0, V, sizeof *V * 3 * (size_t)nv);
#endif
    int *na = malloc(sizeof *na * (size_t)nv), *uf = malloc(sizeof *uf * (size_t)nv), *vl = malloc(sizeof *vl * (size_t)nv);
    int *lab = held ? NULL : malloc(sizeof *lab * (size_t)nt), *pc = malloc(sizeof *pc * (size_t)nv), *part = malloc(sizeof *part * (size_t)nv);
    double *P = malloc(sizeof *P * 3 * (size_t)nv);
    char *in = malloc((size_t)nv), *fix = malloc((size_t)nv), *ch = malloc((size_t)nv), *bad = malloc((size_t)nv);
    unsigned char *lvl = malloc((size_t)nv);
    mk_bvh B = {0}, BH = {0};
    int ret = -1, closed = 0;
    if (!N || !D || !acc || !V0 || !Vr || !Vn || !Vb || !Vk || !sc || !pv0 || !pv || !na || !uf || !vl || (!held && !lab) || !part || !in || !fix || !pc ||
        !P || !ch || !bad || !lvl)
        goto out;
#if defined(MK_MUTATE) && MK_MUTATE == 212
    memcpy(V0, V, sizeof *V * 3 * (size_t)nv);
#endif
    if (held) {
        if (mk_bvh_build(&BH, HV, HT, hnt)) goto out;
        const int st = mk_sepvol_heldcheck_(HV, HT, hnt, &BH, &r);
        if (st < 0) goto out;
        if (st) { /* nothing moved: the report says why */
            r.status = st, ret = -2;
            goto out;
        }
    }
    /* M's parts: joined through faces; whether M is closed, for mk_sepvol_two_ and for what a part's volume means */
    if ((closed = mk_sepvol_closed_(T, nt)) < 0) goto out;
    for (int i = 0; i < nv; i++) uf[i] = i;
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) uf[mk_sepvol_find_(uf, T[3 * (size_t)t + k])] = mk_sepvol_find_(uf, T[3 * (size_t)t + (k + 1) % 3]);
    for (int i = 0; i < nv; i++) part[i] = mk_sepvol_find_(uf, i);
    /* whose crossings with itself count as M's own: any two faces held, one part's in single mode */
    const int *pself = held ? NULL : part;
    mk_sepvol_partvol_(V, T, nt, part, nv, pv0);
    const double crush = MK_SEPVOL_CRUSH;
    /* eps: a thousandth of the mean edge, for the point just inside a vertex */
    double el = 0;
    for (int t = 0; t < nt; t++) {
        const double *a = V + 3 * (size_t)T[3 * (size_t)t], *b = V + 3 * (size_t)T[3 * (size_t)t + 1];
        el += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
    }
    const double eps = 1e-3 * el / nt;
    r.vol0 = mk_sepvol_volume_(V, T, nt);
    long best = LONG_MAX;
    int better = 0; /* the last round that measured fewer than any before */
    for (int round = 0;; round++) {
        mk_sepvol_normals_(V, nv, T, nt, N);
        memcpy(Vr, V, sizeof *V * 3 * (size_t)nv); /* the round's start */
        /* M's tree, over Vb. Single mode asks it where V is inside: rebuilt every round. Held mode only looks for M's
           own crossings in it: kept while no vertex is more than a mean edge from where the tree has it, and its boxes
           grown by that */
        double grow = 0;
        if (held && round) for (int i = 0; i < nv; i++) {
            const double *a = V + 3 * (size_t)i, *b = Vb + 3 * (size_t)i;
            grow = fmax(grow, sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
        }
        if (!held || !round || grow > el / nt) {
            memcpy(Vb, V, sizeof *V * 3 * (size_t)nv), grow = 0;
            mk_bvh_free(&B);
            if (mk_bvh_build(&B, Vb, T, nt)) goto out;
        }
        if (round == 0) {
#if defined(MK_MUTATE) && MK_MUTATE == 191
            /* 191: two parts crossing counted as a self-crossing */
            if ((r.self0 = mk_sepvol_self_(Vr, T, nt, &B, 0, NULL)) < 0) goto out;
#else
            if ((r.self0 = mk_sepvol_self_(Vr, T, nt, &B, 0, pself)) < 0) goto out;
#endif
            if (r.self0) r.status |= MK_SEPVOL_SELF_CROSSES;
        }
        const mk_bvh *O = held ? &BH : &B;
        const long nin = mk_sepvol_mark_(V, nv, N, O, held, !held && closed, eps, in);
        if (round == 0) r.inside0 = nin;
        memset(P, 0, sizeof *P * 3 * (size_t)nv), memset(pc, 0, sizeof *pc * (size_t)nv);
        long nnear;
#if defined(MK_MUTATE) && MK_MUTATE == 195
        const int pairs = 1; /* 195: the switch read as on whatever it says */
#else
        const int pairs = (MK_SEPVOL_PAIRS) != 0;
#endif
        const long np = mk_sepvol_pairs_(V, T, nt, O, held ? HV : NULL, held ? HT : NULL, gap, pairs ? P : NULL, pc, &nnear);
        if (np < 0) goto out;
        if (round == 0) r.cross0 = np;
        /* a state that took too much of a part while what is left is not cleared: the best state before it stands */
        if (closed && (nin || np)) {
            mk_sepvol_partvol_(V, T, nt, part, nv, pv);
#if !(defined(MK_MUTATE) && MK_MUTATE == 220) /* 220: squeezed as far as the rounds go */
            if (mk_sepvol_loss_(pv0, pv, nv) > crush) {
                r.status |= MK_SEPVOL_CRUSHED;
                break;
            }
#else
            (void)crush;
#endif
        }
        if (nin + np <= best) { /* the state kept: the fewest inside vertices and crossing pairs so far, the later on a
                                   tie. A stalled push can end worse than a round before it */
            if (nin + np < best) better = round;
            best = nin + np, r.inside = nin, r.cross = np;
            memcpy(Vk, V, sizeof *V * 3 * (size_t)nv);
        }
        if (round - better >= MK_SEPVOL_STALL) break; /* no better state in that many rounds: stalled */
#if defined(MK_MUTATE) && MK_MUTATE == 153
        if ((!nin && !np && !nnear) || round == rounds || round == 1) break; /* 153: one round only */
#else
        if ((!nin && !np && !nnear) || round == rounds) break;
#endif
        r.rounds = round + 1;
        if (!held) { /* the inside vertices' patches: joined through faces with at least two inside corners */
            for (int i = 0; i < nv; i++) uf[i] = i;
            for (int t = 0; t < nt; t++) {
                const int *q = T + 3 * (size_t)t;
                for (int k = 0; k < 3; k++) {
                    const int a = q[k], b = q[(k + 1) % 3];
                    if (in[a] && in[b]) uf[mk_sepvol_find_(uf, a)] = mk_sepvol_find_(uf, b);
                }
            }
            /* each patch and MK_SEPVOL_OWN rings around it is the vertex's own surface, left out of its search; every
               other face is a candidate (label nv, never a patch's) */
            for (int i = 0; i < nv; i++) vl[i] = in[i] ? mk_sepvol_find_(uf, i) : -1;
            for (int g = 0; g < MK_SEPVOL_OWN; g++) {
                memcpy(na, vl, sizeof *na * (size_t)nv);
                for (int t = 0; t < nt; t++) {
                    const int *q = T + 3 * (size_t)t;
                    for (int k = 0; k < 3; k++)
                        if (na[q[k]] >= 0) for (int j = 0; j < 3; j++) if (vl[q[j]] < 0) vl[q[j]] = na[q[k]];
                }
            }
            for (int t = 0; t < nt; t++) {
                const int *q = T + 3 * (size_t)t;
                lab[t] = in[q[0]] ? vl[q[0]] : in[q[1]] ? vl[q[1]] : in[q[2]] ? vl[q[2]] : vl[q[0]] >= 0 ? vl[q[0]] : vl[q[1]] >= 0 ? vl[q[1]] : vl[q[2]] >= 0 ? vl[q[2]] : nv;
            }
        }
        memset(D, 0, sizeof *D * 3 * (size_t)nv);
        for (int i = 0; i < nv; i++) {
            fix[i] = 0;
            if (!in[i]) continue;
            const double *v = V + 3 * (size_t)i;
            double q[3], dir[3], d;
            int tri = -1;
#if defined(MK_MUTATE) && MK_MUTATE == 151
            /* 151: a push along the vertex's own inward normal, as far as the first hit */
            double tt, bc[3];
            const double nd[3] = {-N[3 * (size_t)i], -N[3 * (size_t)i + 1], -N[3 * (size_t)i + 2]};
            tri = mk_bvh_first_hit(O, v, nd, 1e-9 * el / nt, INFINITY, &tt, bc);
            if (tri < 0) continue;
            d = tt, dir[0] = nd[0], dir[1] = nd[1], dir[2] = nd[2];
            (void)q;
#else
            int skip = held ? 0 : mk_sepvol_find_(uf, i);
#if defined(MK_MUTATE) && MK_MUTATE == 152
            skip = INT_MAX; /* 152: its own patch not left out */
#endif
            const double d2 = held ? mk_bvh_nearest(O, v, INFINITY, &tri, q) : mk_bvh_nearest_lab(O, v, INFINITY, lab, skip, &tri, q);
            if (tri < 0 || !(d2 > 0)) continue;
            d = sqrt(d2);
            for (int k = 0; k < 3; k++) dir[k] = (q[k] - v[k]) / d;
            if (held) mk_sepvol_entry_(O, &B, v, N + 3 * (size_t)i, 1e-9 * el / nt, 2 * grow + 1e-6 * el / nt, &d, dir);
#endif
            double a = (held ? d : 0.5 * d) + gap;
            if (a > step) a = step;
            for (int k = 0; k < 3; k++) D[3 * (size_t)i + k] = a * dir[k];
            fix[i] = 1;
        }
        for (int i = 0; i < nv; i++) /* step 2's pushes, on vertices not inside: crossing and near faces */
            if (!fix[i] && pc[i]) {
                double *e = D + 3 * (size_t)i, l = 0;
                for (int k = 0; k < 3; k++) e[k] = P[3 * (size_t)i + k] / pc[i], l += e[k] * e[k];
                l = sqrt(l);
                if (l > step) for (int k = 0; k < 3; k++) e[k] *= step / l;
                fix[i] = 1;
            }
#if !(defined(MK_MUTATE) && MK_MUTATE == 159) /* 159: no smoothing, the guard alone */
        mk_sepvol_smooth_(D, fix, nv, T, nt, acc, na);
#endif
        for (int ring = 0; ring < rings; ring++) {
            memset(acc, 0, sizeof *acc * 3 * (size_t)nv), memset(na, 0, sizeof *na * (size_t)nv);
            for (int t = 0; t < nt; t++)
                for (int k = 0; k < 3; k++) {
                    const int v = T[3 * (size_t)t + k];
                    if (fix[v]) continue;
                    for (int j = 1; j < 3; j++) {
                        const int u = T[3 * (size_t)t + (k + j) % 3];
                        for (int c = 0; c < 3; c++) acc[3 * (size_t)v + c] += D[3 * (size_t)u + c];
                        na[v]++;
                    }
                }
            for (int i = 0; i < nv; i++)
                if (!fix[i] && na[i]) for (int c = 0; c < 3; c++) D[3 * (size_t)i + c] = 0.5 * acc[3 * (size_t)i + c] / na[i];
        }
        for (int i = 0; i < nv; i++) {
            double *e = D + 3 * (size_t)i;
            const double *u = N + 3 * (size_t)i, out = e[0] * u[0] + e[1] * u[1] + e[2] * u[2];
            if (out > 0) e[0] -= out * u[0], e[1] -= out * u[1], e[2] -= out * u[2]; /* never out of its own part */
            sc[i] = e[0] != 0 || e[1] != 0 || e[2] != 0;
        }
        if (mk_sepvol_guard_(Vr, Vn, D, nv, T, nt, &B, grow, pself, sc, lvl, ch, bad)) goto out;
        long moved = 0;
        for (int i = 0; i < nv; i++) {
            if (!(D[3 * (size_t)i] != 0 || D[3 * (size_t)i + 1] != 0 || D[3 * (size_t)i + 2] != 0)) continue;
            if (lvl[i] > MK_SEPVOL_HALVE) { r.refused++; continue; }
            r.shortened += lvl[i] > 0, moved++;
            memcpy(V + 3 * (size_t)i, Vn + 3 * (size_t)i, sizeof *V * 3);
        }
        if (!moved) break; /* every move refused: another round would try the same ones */
    }
#if !(defined(MK_MUTATE) && MK_MUTATE == 192) /* 192: the last state returned, not the best */
    memcpy(V, Vk, sizeof *V * 3 * (size_t)nv);
#else
    (void)Vk;
#endif
    double grow = 0; /* the tree has V within grow */
    for (int i = 0; i < nv; i++) {
        const double *a = V + 3 * (size_t)i, *b = Vb + 3 * (size_t)i;
        grow = fmax(grow, sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
    }
    if ((r.self1 = mk_sepvol_self_(V, T, nt, &B, grow, pself)) < 0) goto out;
    for (int i = 0; i < nv; i++) {
        const double *a = V + 3 * (size_t)i, *b = V0 + 3 * (size_t)i;
        r.moved += a[0] != b[0] || a[1] != b[1] || a[2] != b[2];
        r.worst = fmax(r.worst, sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
    }
    r.vol1 = mk_sepvol_volume_(V, T, nt);
    if (closed) mk_sepvol_partvol_(V, T, nt, part, nv, pv), r.part_loss = mk_sepvol_loss_(pv0, pv, nv);
#if !(defined(MK_MUTATE) && MK_MUTATE == 221) /* 221: crossing pairs left, and the status silent about them */
    if (r.cross > 0) r.status |= MK_SEPVOL_STILL_CROSSES;
#endif
    ret = 0;
out:
    if (ret) memcpy(V, V0 ? V0 : V, V0 ? sizeof *V * 3 * (size_t)nv : 0);
    mk_bvh_free(&B), mk_bvh_free(&BH);
    free(N), free(D), free(acc), free(V0), free(Vr), free(Vn), free(Vb), free(Vk), free(sc), free(pv0), free(pv), free(na), free(uf), free(vl), free(lab), free(part), free(in), free(fix),
        free(pc), free(P), free(ch), free(bad), free(lvl);
    if (rep && ret != -1) *rep = r;
    return ret;
}
#endif
