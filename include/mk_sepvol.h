/* mathkit: parts of closed outward triangle meshes that sink into each other, moved apart by volume, the mesh kept as
   it is (the same vertices in the same order, the same faces). mk_sepvol, with a report; mk_sepvol_multi for several
   held meshes in one call, mk_sepvol_restore to give volume back after a push. Header only, C23.
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
      toward the other face's (up to 34 degrees off its own face's normal on two cubes). In held mode the moves are then
      smoothed in space, before that last rule: each vertex within 3 sigma of a moving one takes
      the mean of the moves of every vertex of M within 3 sigma, moving or not and whatever part it is on, weighted by a
      Gaussian of the distance; sigma is MK_SEPVOL_SIGMA times the larger of step and M's mean edge. The move is
      smoothed, not the shape: a patch dents as a whole instead of creasing, and pieces close in space but far apart on
      the mesh move alike. A mean is no longer than its longest move, so step still bounds a move.
   4. No move may make M cross itself. The round's moves are tried together; a face that would cross one of M's own
      faces it did not cross at the round's start (any of M's faces in held mode; one of its own part's in single mode,
      where faces of two parts crossing is what is being undone) has its moving corners' moves halved, up to
      MK_SEPVOL_HALVE times, then rolled back (refused), until no such pair is left. Pairs that crossed before are left
      as they are: an input that crosses itself is reported (MK_SEPVOL_SELF_CROSSES) and never made worse. In held mode
      a face also may not fold over 90 degrees against a neighbour across an edge that it met at under 90 at the round's
      start (an edge at 90 or more, as on a cube, is left as it is): faces that share a corner are never a crossing pair,
      and a face turned over onto its neighbour reads, just inside it, as outside the mesh. So every edge folded over 90
      degrees at the end was at 90 or more at the start (rep->fold0, fold1: the length of those over 90). The guard
      holds for the doubles handed back; faces it let come within a rounding of each other can cross once the
      coordinates are written with fewer digits, so write them with %.17g. Steps 3 and 4 work together: on a surface
      crumpled at the edge's scale, unsmoothed moves of half an edge fold it into crossing pairs, and the guard alone
      then stalls.
   mk_sepvol hands back the state it measured with the fewest inside vertices plus crossing pairs (the input is one, so
   the sum never ends above the input's): a push the guard stalls can end worse than a round before it, and then the
   input or an earlier round comes back. Of states with the same sum the earliest, a later one being only more squeezed;
   once nothing is left, the latest (near pairs further apart). A push can trade one count for the other: when the
   crossing pairs end above the input's, MK_SEPVOL_PAIRS_ROSE says so. MK_SEPVOL_STILL_INSIDE and
   MK_SEPVOL_STILL_CROSSES say what is left.
   Rounds repeat until no vertex is inside, no faces cross and none are near, or `rounds` is spent. A run also ends when
   a round's moves are all refused, after MK_SEPVOL_STALL rounds without fewer in sum, or when a part of a closed M
   has lost more than MK_SEPVOL_CRUSH of its volume while vertices are inside or pairs cross (MK_SEPVOL_CRUSHED; the
   best state before that is handed back). Some crossings no vertex move can clear: two long bars crossed like an X have
   no vertex inside each other, and step 2 only squeezes them. mk_sepvol_multi takes several held meshes in one call.
   The held mesh must define inside and outside: closed, every edge of two faces running it opposite ways, no crossing
   pairs, positive volume, and each connected part wound outward unless it lies inside another part, where it is a
   cavity (a hollow solid). Otherwise mk_sepvol returns -2 with the report's status and held_ fields saying why.
   No vertex moves more than step in a round, so none more than rounds times step in all. Lengths are in the caller's
   unit (gap and step in it). mk_closest_pt_tri's face weights can be thousands of ulps off on a triangle of aspect 1e5
   (about 1e-12 of the triangle's size), far under any gap, so this cannot change a push.
   Single-threaded. */
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
#define MK_SEPVOL_STALL 8 /* rounds with no fewer inside vertices plus crossing pairs, summed, before the run ends */
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
#ifndef MK_SEPVOL_SIGMA
#define MK_SEPVOL_SIGMA 1.0 /* held mode, step 3: the moves smoothed in space by a Gaussian of sigma this times the larger of
                               step and M's mean edge; 0 turns it off. Any double expression, read at each call */
#endif
/* the report's status bits. The held ones make mk_sepvol return -2 without moving anything. */
#define MK_SEPVOL_HELD_OPEN 1     /* the held mesh has edges not shared by exactly two faces (holes, or 3+ faces) */
#define MK_SEPVOL_HELD_FLIPPED 2  /* its faces disagree on their facing (an edge run one way by both), a part is wound
                                     inward with no part around it, or it encloses no volume */
#define MK_SEPVOL_HELD_CROSSES 4  /* its faces cross each other (overlapping parts, or a part through itself) */
#define MK_SEPVOL_SELF_CROSSES 8  /* M crossed itself at the start (see self0): it ran, and never added to them */
#define MK_SEPVOL_STILL_CROSSES 16 /* faces still cross in the state handed back (rep->cross): parts not separated */
#define MK_SEPVOL_CRUSHED 32      /* stopped: going on would have taken more than MK_SEPVOL_CRUSH of a part */
#define MK_SEPVOL_STILL_INSIDE 64 /* vertices are still inside in the state handed back (rep->inside) */
#define MK_SEPVOL_PAIRS_ROSE 128  /* the state handed back has more crossing pairs than the input (rep->cross > cross0),
                                     though inside vertices plus crossing pairs are no more */
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
    double fold0, fold1; /* the length of M's edges folded over 90 degrees (its two faces' normals more than 90 degrees apart),
                            before and after, in the caller's unit */
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
   refuses them (lvl[i] > MK_SEPVOL_HALVE), until no such pair is left. So does a face that would fold over 90 degrees
   against a neighbour across an edge (nb: the faces across each face's three edges, -1 for none) it met at under 90
   degrees at the start (an edge at 90 or more is left as it is): two faces that share a corner are never a crossing
   pair, and a face turned over onto its neighbour reads, just inside it, as outside the mesh. Vn: the positions after. 0, or -1 for no memory. */
static inline int mk_sepvol_guard_(const double *V, double *Vn, const double *D, int nv, const int *T, int nt, const mk_bvh *B, double grow,
                                   const int *part, const int *nb, double *s, unsigned char *lvl, char *ch, char *bad)
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
#if !(defined(MK_MUTATE) && MK_MUTATE == 343) /* 343: folds against a neighbour let through */
            for (int k = 0; nb && k < 3; k++) {
                const int u = nb[3 * (size_t)t + k];
                if (u < 0) continue;
                const int *Y = T + 3 * (size_t)u;
                double a[3], b[3];
                mk_sepvol_fnormal_(Vn, X, a), mk_sepvol_fnormal_(Vn, Y, b);
                if (!(a[0] * b[0] + a[1] * b[1] + a[2] * b[2] < 0)) continue;
                mk_sepvol_fnormal_(V, X, a), mk_sepvol_fnormal_(V, Y, b);
                if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] < 1e-9) continue; /* 90 degrees or more at the start: left as it is */
                for (int c = 0; c < 3; c++) {
                    if (s[X[c]] > 0) bad[X[c]] = 1;
                    if (s[Y[c]] > 0) bad[Y[c]] = 1;
                }
                nbad++;
            }
#else
            (void)nb;
#endif
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

/* a grid of cells R across over V: cell q's vertices are a list in its bucket (head[hn], hn a power of 2; next[nv]) */
static inline void mk_sepvol_cell_(const double *p, const double *lo, double R, long long *c)
{
    for (int k = 0; k < 3; k++) {
        const double q = floor((p[k] - lo[k]) / R);
        c[k] = q < 0x1p52 ? (long long)q : (1LL << 52); /* far points share the last cell: only the distance test decides */
    }
}
static inline int mk_sepvol_bucket_(const long long *c, int hn)
{
    const uint64_t h = (uint64_t)c[0] * 0x9E3779B97F4A7C15u ^ (uint64_t)c[1] * 0xC2B2AE3D27D4EB4Fu ^ (uint64_t)c[2] * 0x165667B19E3779F9u;
    return (int)((h ^ h >> 29) & (uint64_t)(hn - 1));
}
/* the distinct buckets of the 27 cells round p's: their count, into b[27] */
static inline int mk_sepvol_near_(const double *p, const double *lo, double R, int hn, int *b)
{
    long long c[3], q[3];
    mk_sepvol_cell_(p, lo, R, c);
    int n = 0;
    for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++)
            for (int dz = -1; dz <= 1; dz++) {
                q[0] = c[0] + dx, q[1] = c[1] + dy, q[2] = c[2] + dz;
                const int h = mk_sepvol_bucket_(q, hn);
                int seen = 0;
                for (int i = 0; i < n && !seen; i++) seen = b[i] == h;
                if (!seen) b[n++] = h;
            }
    return n;
}
/* step 3 in held mode: the round's moves D smoothed in space: each vertex within 3 sigma of a
   moving one takes the mean of the moves of every vertex of M within 3 sigma, moving or not, weighted by a Gaussian of
   their distance (sigma), whatever part they are on, so pieces close in space move together however far apart they are
   on the mesh. A mean is no longer than the longest move in it, so no move grows past step. V the positions; head[hn]
   (hn a power of 2 >= nv) and next[nv] the grid; acc[3 nv] scratch, tch[nv] marks. */
static inline void mk_sepvol_gauss_(const double *V, int nv, double *D, double sigma, int *head, int *next, int hn, double *acc, char *tch)
{
    if (nv < 1) return;
    const double R = 3 * sigma, R2 = R * R, k = -0.5 / (sigma * sigma);
    double lo[3] = {INFINITY, INFINITY, INFINITY};
    for (int i = 0; i < nv; i++) for (int c = 0; c < 3; c++) lo[c] = fmin(lo[c], V[3 * (size_t)i + c]);
    for (int h = 0; h < hn; h++) head[h] = -1;
    for (int i = nv - 1; i >= 0; i--) { /* each list in index order */
        long long c[3];
        mk_sepvol_cell_(V + 3 * (size_t)i, lo, R, c);
        const int h = mk_sepvol_bucket_(c, hn);
        next[i] = head[h], head[h] = i;
    }
    memset(acc, 0, sizeof *acc * 3 * (size_t)nv), memset(tch, 0, (size_t)nv);
    int b[27];
    for (int j = 0; j < nv; j++) { /* each move spread over its neighbours */
        const double *e = D + 3 * (size_t)j, *p = V + 3 * (size_t)j;
        if (e[0] == 0 && e[1] == 0 && e[2] == 0) continue;
        const int nb = mk_sepvol_near_(p, lo, R, hn, b);
        for (int q = 0; q < nb; q++)
            for (int i = head[b[q]]; i >= 0; i = next[i]) {
                const double *x = V + 3 * (size_t)i, r2 = (x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]) + (x[2] - p[2]) * (x[2] - p[2]);
                if (!(r2 <= R2)) continue;
                const double w = exp(k * r2);
                for (int c = 0; c < 3; c++) acc[3 * (size_t)i + c] += w * e[c];
                tch[i] = 1;
            }
    }
    for (int i = 0; i < nv; i++) { /* over the weights of every vertex near it */
        if (!tch[i]) continue;
        const double *p = V + 3 * (size_t)i;
        double ws = 0;
        const int nb = mk_sepvol_near_(p, lo, R, hn, b);
        for (int q = 0; q < nb; q++)
            for (int j = head[b[q]]; j >= 0; j = next[j]) {
                const double *x = V + 3 * (size_t)j, r2 = (x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]) + (x[2] - p[2]) * (x[2] - p[2]);
                if (r2 <= R2) ws += exp(k * r2);
            }
        for (int c = 0; c < 3; c++) D[3 * (size_t)i + c] = acc[3 * (size_t)i + c] / ws;
    }
}
/* the faces across each face's edges: nb[3 t + k] the one face other than t on edge T[3 t + k], T[3 t + (k + 1) % 3],
   -1 when the edge has no other face or more than one. 0, or -1 for no memory */
static inline int mk_sepvol_adj_(const int *T, int nt, int nv, int *nb)
{
    int *off = calloc((size_t)nv + 1, sizeof *off), *inc = malloc(sizeof *inc * 3 * (size_t)nt);
    if (!off || !inc) { free(off), free(inc); return -1; }
    for (size_t q = 0; q < 3 * (size_t)nt; q++) off[T[q] + 1]++;
    for (int i = 0; i < nv; i++) off[i + 1] += off[i];
    for (size_t q = 0; q < 3 * (size_t)nt; q++) inc[off[T[q]]++] = (int)(q / 3);
    for (int i = nv; i > 0; i--) off[i] = off[i - 1];
    off[0] = 0;
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) {
            const int a = T[3 * (size_t)t + k], b = T[3 * (size_t)t + (k + 1) % 3];
            int other = -1, n = 0;
            for (int j = off[a]; j < off[a + 1] && a != b; j++) {
                const int u = inc[j], *Y = T + 3 * (size_t)u;
                if (u != t && (Y[0] == b || Y[1] == b || Y[2] == b)) other = u, n++;
            }
            nb[3 * (size_t)t + k] = n == 1 ? other : -1;
        }
    free(off), free(inc);
    return 0;
}
/* the length of M's edges folded over 90 degrees: edges of two faces (nb) whose normals are more than 90 degrees apart
   (meshq's f90 is it over the area) */
static inline double mk_sepvol_fold_(const double *V, const int *T, int nt, const int *nb)
{
    double L = 0;
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) {
            const int u = nb[3 * (size_t)t + k];
            if (u < t) continue; /* each edge once, from its lower face */
            double a[3], b[3];
            mk_sepvol_fnormal_(V, T + 3 * (size_t)t, a), mk_sepvol_fnormal_(V, T + 3 * (size_t)u, b);
#if defined(MK_MUTATE) && MK_MUTATE == 345
            if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] < 0.5) { /* 345: folds over 60 degrees counted */
#else
            if (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] < 0) {
#endif
                const double *p = V + 3 * (size_t)T[3 * (size_t)t + k], *q = V + 3 * (size_t)T[3 * (size_t)t + (k + 1) % 3];
                L += sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]));
            }
        }
    return L;
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
    if (rep) memset(rep, 0, sizeof *rep); /* padding too: a caller may compare the bytes */
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
    int *lab = held ? NULL : malloc(sizeof *lab * (size_t)nt), *pc = malloc(sizeof *pc * (size_t)nv), *part = calloc((size_t)nv, sizeof *part);
    double *P = malloc(sizeof *P * 3 * (size_t)nv);
    char *in = malloc((size_t)nv), *fix = malloc((size_t)nv), *ch = malloc((size_t)nv), *bad = malloc((size_t)nv);
    unsigned char *lvl = malloc((size_t)nv);
    int *head = NULL, *next = NULL, *nb = NULL, hn = 1;
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
    /* each face's neighbours across its edges: for the guard's folds and the folds reported */
    if (!(nb = malloc(sizeof *nb * 3 * (size_t)nt)) || mk_sepvol_adj_(T, nt, nv, nb)) goto out;
    r.fold0 = mk_sepvol_fold_(V, T, nt, nb);
    /* held mode's smoothing in space (step 3): its sigma, and the grid of cells 3 sigma across it looks for neighbours in */
    const double sig = held ? (MK_SEPVOL_SIGMA) * fmax(step, el / nt) : 0;
    if (sig > 0) {
        while (hn < nv && hn < 1 << 30) hn *= 2;
        if (!(head = malloc(sizeof *head * (size_t)hn)) || !(next = malloc(sizeof *next * (size_t)nv))) goto out;
    }
    long best = LONG_MAX;
    long kept = LONG_MAX;
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
        /* the state kept: the fewest inside vertices plus crossing pairs so far (the input is one, so the sum never
           ends above the input's); the earliest on a tie, the latest once nothing is left. A stalled push can end worse
           than a round before it */
        if (nin + np < best) better = round, best = nin + np;
#if defined(MK_MUTATE) && MK_MUTATE == 340
        const long key = nin; /* 340: the fewest inside vertices, whatever the crossing pairs did */
#else
        const long key = nin + np;
#endif
        if (key < kept || (key == kept && !key)) {
            kept = key, r.inside = nin, r.cross = np;
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
#if !(defined(MK_MUTATE) && MK_MUTATE == 341) /* 341: no smoothing in space */
        if (sig > 0) mk_sepvol_gauss_(V, nv, D, sig, head, next, hn, acc, ch);
#endif
        for (int i = 0; i < nv; i++) {
            double *e = D + 3 * (size_t)i;
            const double *u = N + 3 * (size_t)i, out = e[0] * u[0] + e[1] * u[1] + e[2] * u[2];
            if (out > 0) e[0] -= out * u[0], e[1] -= out * u[1], e[2] -= out * u[2]; /* never out of its own part */
            sc[i] = e[0] != 0 || e[1] != 0 || e[2] != 0;
        }
        if (mk_sepvol_guard_(Vr, Vn, D, nv, T, nt, &B, grow, pself, held ? nb : NULL, sc, lvl, ch, bad)) goto out;
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
#if !(defined(MK_MUTATE) && MK_MUTATE == 342) /* 342: vertices left inside, and the status silent about them */
    if (r.inside > 0) r.status |= MK_SEPVOL_STILL_INSIDE;
#endif
#if !(defined(MK_MUTATE) && MK_MUTATE == 348) /* 348: more crossing pairs than the input's, and the status silent */
    if (r.cross > r.cross0) r.status |= MK_SEPVOL_PAIRS_ROSE;
#endif
    r.fold1 = mk_sepvol_fold_(V, T, nt, nb);
    ret = 0;
out:
    if (ret) memcpy(V, V0 ? V0 : V, V0 ? sizeof *V * 3 * (size_t)nv : 0);
    mk_bvh_free(&B), mk_bvh_free(&BH);
    free(N), free(D), free(acc), free(V0), free(Vr), free(Vn), free(Vb), free(Vk), free(sc), free(pv0), free(pv), free(na), free(uf), free(vl), free(lab), free(part), free(in), free(fix),
        free(pc), free(P), free(ch), free(bad), free(lvl), free(head), free(next), free(nb);
    if (rep && ret != -1) *rep = r;
    return ret;
}

/* crossing pairs of M's faces F[0, nf) (all nf = nt faces with F NULL) with H's (BH its tree over HV, HT), counted as
   mk_sepvol_pairs_ counts them; M's positions V. -1: no memory */
static inline long mk_sepvol_xcount_(const double *V, const int *T, const int *F, int nf, const mk_bvh *BH, const double *HV, const int *HT)
{
    int cap = 64, *buf = malloc(sizeof *buf * (size_t)cap);
    long np = 0;
    if (!buf) return -1;
    for (int f = 0; f < nf; f++) {
        const int t = F ? F[f] : f, *X = T + 3 * (size_t)t;
        mk_aabb bx;
        mk_tri_box_(V, T, t, &bx);
        const int n = mk_sepvol_box_(BH, &bx, &buf, &cap);
        if (n < 0) { free(buf); return -1; }
        for (int i = 0; i < n; i++) {
            const int *Y = HT + 3 * (size_t)buf[i];
            np += mk_sepvol_edge_through_(V, X, HV, Y) || mk_sepvol_edge_through_(HV, Y, V, X);
        }
    }
    free(buf);
    return np;
}
static inline int mk_sepvol_boxes_meet_(const mk_aabb *a, const mk_aabb *b)
{
    for (int k = 0; k < 3; k++) if (a->hi[k] < b->lo[k] || b->hi[k] < a->lo[k]) return 0;
    return 1;
}
/* M's vertices S[0, ns) inside H (its tree BH, its box hb), at positions V */
static inline long mk_sepvol_xin_(const double *V, const int *S, int ns, const mk_bvh *BH, const mk_aabb *hb)
{
    long n = 0;
    for (int s = 0; s < ns; s++) {
        const double *v = V + 3 * (size_t)(S ? S[s] : s);
        if (mk_aabb_dist2(hb, v) == 0) n += mk_bvh_inside(BH, v);
    }
    return n;
}

#ifndef MK_SEPVOL_PASSES
#define MK_SEPVOL_PASSES 2 /* mk_sepvol_multi: passes over the held meshes */
#endif
/* Several held meshes in one call, each kept separate: M moved out of held mesh k (HV[k] with hnv[k] vertices, HT[k] with
   hnt[k] faces, k < nh) by mk_sepvol, one mesh at a time, the one with the most inside vertices plus crossing pairs first
   MK_SEPVOL_PASSES passes over the meshes that still have something left. Each is its own
   inside and outside: no union is formed, so held meshes that overlap or touch are taken, where their union would cross
   itself. A push out of one mesh can take M into another; the state handed back is the one after a push (or the input)
   with the fewest inside vertices plus crossing pairs summed over the meshes, so the sum never ends above the input's;
   the earliest on a tie, the latest once nothing is left. A mesh whose last push moved nothing is not pushed again
   until something else has moved. A held mesh that cannot say inside from outside is left out: its report's status and
   held_ fields say why, and rep->status carries its bits.
   rep: inside0, inside, cross0, cross summed over the held meshes (a vertex inside two counts twice); M's own crossings,
   volume, part_loss, folds, moved and worst against the input; rounds, refused and shortened summed over the pushes;
   status as mk_sepvol's for the sums (MK_SEPVOL_PAIRS_ROSE: the summed crossing pairs above the input's), with the held
   bits of any mesh left out. hrep[k] (may be NULL): inside0 and cross0 with mesh k at the start, inside and cross in the
   state handed back, status as mk_sepvol's for those (PAIRS_ROSE: its own crossing pairs rose), rounds, refused and
   shortened summed over its pushes, its held_ fields. Returns 0 when it ran, -1 as mk_sepvol (and for nh < 1 or a NULL mesh), -2
   when every held mesh was left out. On -1 and -2 V is as it came in. */
static inline int mk_sepvol_multi(double *V, int nv, const int *T, int nt, int nh, const double *const *HV, const int *hnv, const int *const *HT,
                                  const int *hnt, double gap, double step, int rings, int rounds, mk_sepvol_report *rep, mk_sepvol_report *hrep)
{
    mk_sepvol_report r = {0};
    if (rep) memset(rep, 0, sizeof *rep); /* padding too: a caller may compare the bytes */
    if (hrep && nh > 0) memset(hrep, 0, sizeof *hrep * (size_t)nh);
    if (nh < 1 || !HV || !hnv || !HT || !hnt || !mk_sepvol_ok_(V, nv, T, nt) || !(gap >= 0) || !isfinite(gap) || !(step > 0) || !isfinite(step) ||
        rings < 0 || rounds < 1)
        return -1;
    for (int k = 0; k < nh; k++)
        if (!HV[k] || !HT[k] || !mk_sepvol_ok_(HV[k], hnv[k], HT[k], hnt[k])) return -1;
    mk_bvh *BH = calloc((size_t)nh, sizeof *BH), B = {0};
    mk_aabb *hb = malloc(sizeof *hb * (size_t)nh);
    mk_sepvol_report *hr = calloc((size_t)nh, sizeof *hr);
    long *in0 = malloc(sizeof *in0 * (size_t)nh), *cr0 = malloc(sizeof *cr0 * (size_t)nh), *in = malloc(sizeof *in * (size_t)nh), *cr = malloc(sizeof *cr * (size_t)nh);
    long *ink = malloc(sizeof *ink * (size_t)nh), *crk = malloc(sizeof *crk * (size_t)nh);
    int *ord = malloc(sizeof *ord * (size_t)nh), *S = malloc(sizeof *S * (size_t)nv), *F = malloc(sizeof *F * (size_t)nt), *uf = malloc(sizeof *uf * (size_t)nv),
        *part = calloc((size_t)nv, sizeof *part), *nb = malloc(sizeof *nb * 3 * (size_t)nt);
    char *ok = calloc((size_t)nh, 1), *mv = malloc((size_t)nv), *still = calloc((size_t)nh, 1);
    double *V0 = malloc(sizeof *V0 * 3 * (size_t)nv), *Vs = malloc(sizeof *Vs * 3 * (size_t)nv), *Vk = malloc(sizeof *Vk * 3 * (size_t)nv),
           *pv0 = malloc(sizeof *pv0 * (size_t)nv), *pv = malloc(sizeof *pv * (size_t)nv);
    int ret = -1, nok = 0, closed = 0;
#if !(defined(MK_MUTATE) && MK_MUTATE == 346) /* 346: V0 filled after every allocation is checked */
    if (V0) memcpy(V0, V, sizeof *V * 3 * (size_t)nv); /* before any goto out: the restore there copies V0 back */
#endif
    if (!BH || !hb || !hr || !in0 || !cr0 || !in || !cr || !ink || !crk || !ord || !S || !F || !uf || !part || !nb || !ok || !mv || !still || !V0 || !Vs || !Vk || !pv0 ||
        !pv)
        goto out;
#if defined(MK_MUTATE) && MK_MUTATE == 346
    memcpy(V0, V, sizeof *V * 3 * (size_t)nv);
#endif
    for (int k = 0; k < nh; k++) { /* the held meshes: each checked, its tree, its counts with M at the start */
        if (mk_bvh_build(BH + k, HV[k], HT[k], hnt[k])) goto out;
        const int st = mk_sepvol_heldcheck_(HV[k], HT[k], hnt[k], BH + k, hr + k);
        if (st < 0) goto out;
        hr[k].status = st, r.status |= st, ok[k] = !st, nok += !st;
        mk_aabb_empty(hb + k);
        for (int i = 0; i < hnv[k]; i++) mk_aabb_add(hb + k, HV[k] + 3 * (size_t)i);
        in0[k] = cr0[k] = 0;
        if (ok[k]) {
            in0[k] = mk_sepvol_xin_(V, NULL, nv, BH + k, hb + k);
            if ((cr0[k] = mk_sepvol_xcount_(V, T, NULL, nt, BH + k, HV[k], HT[k])) < 0) goto out;
        }
        in[k] = ink[k] = in0[k], cr[k] = crk[k] = cr0[k], hr[k].inside0 = in0[k], hr[k].cross0 = cr0[k];
        r.inside0 += in0[k], r.cross0 += cr0[k];
        ord[k] = k;
    }
    if (!nok) { ret = -2; goto out; }
    for (int a = 1; a < nh; a++) /* the most inside vertices plus crossing pairs first; ties in the caller's order */
        for (int b = a; b > 0 && in0[ord[b]] + cr0[ord[b]] > in0[ord[b - 1]] + cr0[ord[b - 1]]; b--) {
            const int s = ord[b];
            ord[b] = ord[b - 1], ord[b - 1] = s;
        }
    if ((closed = mk_sepvol_closed_(T, nt)) < 0) goto out;
    for (int i = 0; i < nv; i++) uf[i] = i;
    for (int t = 0; t < nt; t++)
        for (int k = 0; k < 3; k++) uf[mk_sepvol_find_(uf, T[3 * (size_t)t + k])] = mk_sepvol_find_(uf, T[3 * (size_t)t + (k + 1) % 3]);
    for (int i = 0; i < nv; i++) part[i] = mk_sepvol_find_(uf, i);
    mk_sepvol_partvol_(V, T, nt, part, nv, pv0);
    r.vol0 = mk_sepvol_volume_(V, T, nt);
    if (mk_sepvol_adj_(T, nt, nv, nb)) goto out;
    r.fold0 = mk_sepvol_fold_(V, T, nt, nb);
    if (mk_bvh_build(&B, V, T, nt)) goto out;
    if ((r.self0 = mk_sepvol_self_(V, T, nt, &B, 0, NULL)) < 0) goto out;
    mk_bvh_free(&B);
    memcpy(Vk, V, sizeof *V * 3 * (size_t)nv);
    long kept = r.inside0 + r.cross0;
    for (int pass = 0; pass < MK_SEPVOL_PASSES; pass++) {
        int pushed = 0;
        for (int o = 0; o < nh; o++) {
            const int k = ord[o];
            if (!ok[k] || !(in[k] + cr[k]) || still[k]) continue; /* still: its last push moved nothing, and nothing has since */
            memcpy(Vs, V, sizeof *V * 3 * (size_t)nv);
            mk_sepvol_report q;
            const int rc = mk_sepvol(V, nv, T, nt, HV[k], hnv[k], HT[k], hnt[k], gap, step, rings, rounds, &q);
            if (rc == -1) goto out;
            hr[k].rounds += q.rounds, hr[k].refused += q.refused, hr[k].shortened += q.shortened;
            r.rounds += q.rounds, r.refused += q.refused, r.shortened += q.shortened, r.status |= q.status & MK_SEPVOL_CRUSHED;
            int ns = 0, nf = 0; /* what moved, and the faces it moved */
            for (int i = 0; i < nv; i++) {
                const double *a = V + 3 * (size_t)i, *b = Vs + 3 * (size_t)i;
                mv[i] = a[0] != b[0] || a[1] != b[1] || a[2] != b[2];
                if (mv[i]) S[ns++] = i;
            }
            if (!ns) { still[k] = 1; continue; }
            pushed = 1, memset(still, 0, (size_t)nh);
            mk_aabb mb;
            mk_aabb_empty(&mb);
            for (int t = 0; t < nt; t++) /* the faces it moved, and the box of all they span before and after */
                if (mv[T[3 * (size_t)t]] || mv[T[3 * (size_t)t + 1]] || mv[T[3 * (size_t)t + 2]]) {
                    F[nf++] = t;
                    for (int c = 0; c < 3; c++) mk_aabb_add(&mb, V + 3 * (size_t)T[3 * (size_t)t + c]), mk_aabb_add(&mb, Vs + 3 * (size_t)T[3 * (size_t)t + c]);
                }
            for (int j = 0; j < nh; j++) { /* each other held mesh: its counts after the push, from what moved */
                if (j == k || !ok[j] || !mk_sepvol_boxes_meet_(&mb, hb + j)) continue;
                const long a1 = mk_sepvol_xcount_(V, T, F, nf, BH + j, HV[j], HT[j]), a0 = mk_sepvol_xcount_(Vs, T, F, nf, BH + j, HV[j], HT[j]);
                if (a1 < 0 || a0 < 0) goto out;
                in[j] += mk_sepvol_xin_(V, S, ns, BH + j, hb + j) - mk_sepvol_xin_(Vs, S, ns, BH + j, hb + j), cr[j] += a1 - a0;
            }
            in[k] = q.inside, cr[k] = q.cross;
            long sum = 0;
            for (int j = 0; j < nh; j++) sum += in[j] + cr[j];
#if defined(MK_MUTATE) && MK_MUTATE == 344
            sum = 0; /* 344: the last state handed back, whatever it did to the other held meshes */
#endif
            if (sum < kept || (sum == kept && !sum)) {
                kept = sum;
                memcpy(Vk, V, sizeof *V * 3 * (size_t)nv), memcpy(ink, in, sizeof *in * (size_t)nh), memcpy(crk, cr, sizeof *cr * (size_t)nh);
            }
        }
        if (!pushed) break;
    }
    memcpy(V, Vk, sizeof *V * 3 * (size_t)nv);
    for (int k = 0; k < nh; k++) {
        hr[k].inside = ink[k], hr[k].cross = crk[k], r.inside += ink[k], r.cross += crk[k];
        if (ok[k] && crk[k]) hr[k].status |= MK_SEPVOL_STILL_CROSSES;
        if (ok[k] && ink[k]) hr[k].status |= MK_SEPVOL_STILL_INSIDE;
        if (ok[k] && crk[k] > cr0[k]) hr[k].status |= MK_SEPVOL_PAIRS_ROSE;
    }
    if (mk_bvh_build(&B, V, T, nt)) goto out;
    if ((r.self1 = mk_sepvol_self_(V, T, nt, &B, 0, NULL)) < 0) goto out;
    if (r.self0) r.status |= MK_SEPVOL_SELF_CROSSES;
    if (r.cross) r.status |= MK_SEPVOL_STILL_CROSSES;
    if (r.inside) r.status |= MK_SEPVOL_STILL_INSIDE;
    if (r.cross > r.cross0) r.status |= MK_SEPVOL_PAIRS_ROSE;
    for (int i = 0; i < nv; i++) {
        const double *a = V + 3 * (size_t)i, *b = V0 + 3 * (size_t)i;
        r.moved += a[0] != b[0] || a[1] != b[1] || a[2] != b[2];
        r.worst = fmax(r.worst, sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
    }
    r.vol1 = mk_sepvol_volume_(V, T, nt);
    if (closed) mk_sepvol_partvol_(V, T, nt, part, nv, pv), r.part_loss = mk_sepvol_loss_(pv0, pv, nv);
    r.fold1 = mk_sepvol_fold_(V, T, nt, nb);
    ret = 0;
out:
    if (ret == -1 && V0) memcpy(V, V0, sizeof *V * 3 * (size_t)nv);
    mk_bvh_free(&B);
    for (int k = 0; BH && k < nh; k++) mk_bvh_free(BH + k);
    if (hrep && ret != -1) memcpy(hrep, hr, sizeof *hr * (size_t)nh);
    if (rep && ret != -1) *rep = r;
    free(BH), free(hb), free(hr), free(in0), free(cr0), free(in), free(cr), free(ink), free(crk), free(ord), free(S), free(F), free(uf), free(part), free(nb), free(ok),
        free(mv), free(still), free(V0), free(Vs), free(Vk), free(pv0), free(pv);
    return ret;
}

/* M's volume given back toward vol where there is room. The movers: M's vertices outside every
   held mesh (k < nh; nh 0: none) none of whose faces has a held face within depth of its box. They move out along their
   normals, all by one distance, so that M's enclosed volume meets vol: MK_SEPVOL_RESTORE steps of Newton's method (the
   rate: the faces they touch, each against their normals), each at most depth / 2, the movers chosen again at each
   step. A face moving less than depth stays clear of every held mesh, so none is crossed, entered or swept over: no
   vertex ends inside a held mesh, and no face crosses one, that did not at the start. Each step goes through held mode's
   guard (no new own crossing pair, no new fold over 90 degrees); where it holds moves back the volume falls short of
   vol. rep: vol0, vol1, moved, worst, self0, self1, fold0, fold1, inside0, inside, cross0, cross (summed over the held
   meshes), rounds (the steps run), refused, shortened; status MK_SEPVOL_SELF_CROSSES when M crossed itself at the start.
   Returns 0 when it ran, -1 as mk_sepvol_multi (and for depth not positive or vol not finite), -2 when a held mesh
   cannot say inside from outside (rep->status has its bits). On -1 and -2 V is as it came in. */
#ifndef MK_SEPVOL_RESTORE
#define MK_SEPVOL_RESTORE 4 /* mk_sepvol_restore's Newton steps */
#endif
static inline int mk_sepvol_restore(double *V, int nv, const int *T, int nt, int nh, const double *const *HV, const int *hnv, const int *const *HT,
                                    const int *hnt, double vol, double depth, mk_sepvol_report *rep)
{
    mk_sepvol_report r = {0};
    if (rep) memset(rep, 0, sizeof *rep); /* padding too: a caller may compare the bytes */
    if (nh < 0 || (nh && (!HV || !hnv || !HT || !hnt)) || !mk_sepvol_ok_(V, nv, T, nt) || !isfinite(vol) || !(depth > 0) || !isfinite(depth)) return -1;
    for (int k = 0; k < nh; k++)
        if (!HV[k] || !HT[k] || !mk_sepvol_ok_(HV[k], hnv[k], HT[k], hnt[k])) return -1;
    mk_bvh *BH = calloc((size_t)(nh ? nh : 1), sizeof *BH), B = {0};
    mk_aabb *hb = malloc(sizeof *hb * (size_t)(nh ? nh : 1));
    double *V0 = malloc(sizeof *V0 * 3 * (size_t)nv), *N = malloc(sizeof *N * 3 * (size_t)nv), *D = malloc(sizeof *D * 3 * (size_t)nv),
           *Vn = malloc(sizeof *Vn * 3 * (size_t)nv), *s = malloc(sizeof *s * (size_t)nv);
    int *nb = malloc(sizeof *nb * 3 * (size_t)nt);
    char *ch = malloc((size_t)nv), *bad = malloc((size_t)nv);
    int cap = 64, *buf = malloc(sizeof *buf * (size_t)cap);
    unsigned char *lvl = malloc((size_t)nv);
    int ret = -1;
    if (V0) memcpy(V0, V, sizeof *V * 3 * (size_t)nv); /* before any goto out: the restore there copies V0 back */
    if (!BH || !hb || !V0 || !N || !D || !Vn || !s || !nb || !ch || !bad || !buf || !lvl) goto out;
    for (int k = 0; k < nh; k++) {
        if (mk_bvh_build(BH + k, HV[k], HT[k], hnt[k])) goto out;
        mk_sepvol_report hr = {0};
        const int st = mk_sepvol_heldcheck_(HV[k], HT[k], hnt[k], BH + k, &hr);
        if (st < 0) goto out;
        if (st) { r.status |= st, ret = -2; goto out; }
        mk_aabb_empty(hb + k);
        for (int i = 0; i < hnv[k]; i++) mk_aabb_add(hb + k, HV[k] + 3 * (size_t)i);
        r.inside0 += mk_sepvol_xin_(V, NULL, nv, BH + k, hb + k);
        const long c = mk_sepvol_xcount_(V, T, NULL, nt, BH + k, HV[k], HT[k]);
        if (c < 0) goto out;
        r.cross0 += c;
    }
    if (mk_sepvol_adj_(T, nt, nv, nb)) goto out;
    r.fold0 = mk_sepvol_fold_(V, T, nt, nb), r.vol0 = mk_sepvol_volume_(V, T, nt);
    if (mk_bvh_build(&B, V, T, nt)) goto out;
    if ((r.self0 = mk_sepvol_self_(V, T, nt, &B, 0, NULL)) < 0) goto out;
    if (r.self0) r.status |= MK_SEPVOL_SELF_CROSSES;
    for (int it = 0; it < MK_SEPVOL_RESTORE; it++) {
        const double v = mk_sepvol_volume_(V, T, nt);
        if (!(v != vol)) break;
        mk_sepvol_normals_(V, nv, T, nt, N);
        long nm = 0;
        for (int i = 0; i < nv; i++) { /* the movers: outside every held mesh */
            const double *p = V + 3 * (size_t)i;
            int ok = 1;
            for (int k = 0; k < nh && ok; k++) ok = !(mk_aabb_dist2(hb + k, p) == 0 && mk_bvh_inside(BH + k, p));
            s[i] = ok;
        }
        for (int t = 0; t < nt; t++) { /* and on no face with a held face within depth of its box: a face moving less than
                                          depth stays clear of every held mesh, so none is crossed, entered or swept over */
            const int *X = T + 3 * (size_t)t;
            if (!s[X[0]] && !s[X[1]] && !s[X[2]]) continue;
            mk_aabb fb;
            mk_tri_box_(V, T, t, &fb);
            for (int c = 0; c < 3; c++) fb.lo[c] -= depth, fb.hi[c] += depth;
            int near = 0;
            for (int k = 0; k < nh && !near; k++) {
                if (!mk_sepvol_boxes_meet_(&fb, hb + k)) continue;
                const int n = mk_sepvol_box_(BH + k, &fb, &buf, &cap);
                if (n < 0) goto out;
                near = n > 0;
            }
#if defined(MK_MUTATE) && MK_MUTATE == 347
            near = 0; /* 347: faces near a held mesh move too */
#endif
            if (near) s[X[0]] = s[X[1]] = s[X[2]] = 0;
        }
        for (int i = 0; i < nv; i++) nm += s[i] > 0;
        double area = 0;
        for (int t = 0; t < nt; t++) {
            const int *X = T + 3 * (size_t)t;
            if (!s[X[0]] && !s[X[1]] && !s[X[2]]) continue;
            const double *a = V + 3 * (size_t)X[0], *b = V + 3 * (size_t)X[1], *c = V + 3 * (size_t)X[2];
            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            const double n[3] = {u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]};
            double f = 0; /* the face's rate of volume: half its area vector against its moving corners' mean normal */
            for (int c = 0; c < 3; c++)
                if (s[X[c]]) f += (n[0] * N[3 * (size_t)X[c]] + n[1] * N[3 * (size_t)X[c] + 1] + n[2] * N[3 * (size_t)X[c] + 2]) / 6;
            area += f;
        }
        if (!nm || !(area > 0)) break;
        double d = (vol - v) / area;
        if (fabs(d) > 0.5 * depth) d = copysign(0.5 * depth, d);
        for (int i = 0; i < nv; i++) for (int c = 0; c < 3; c++) D[3 * (size_t)i + c] = s[i] ? d * N[3 * (size_t)i + c] : 0;
        mk_bvh_free(&B);
        if (mk_bvh_build(&B, V, T, nt)) goto out;
        if (mk_sepvol_guard_(V, Vn, D, nv, T, nt, &B, 0, NULL, nb, s, lvl, ch, bad)) goto out;
        long moved = 0;
        for (int i = 0; i < nv; i++) {
            if (!(D[3 * (size_t)i] != 0 || D[3 * (size_t)i + 1] != 0 || D[3 * (size_t)i + 2] != 0)) continue;
            if (!(s[i] > 0)) { r.refused++; continue; }
            r.shortened += lvl[i] > 0, moved++;
            memcpy(V + 3 * (size_t)i, Vn + 3 * (size_t)i, sizeof *V * 3);
        }
        r.rounds = it + 1;
        if (!moved) break;
    }
    for (int k = 0; k < nh; k++) { /* the counts after */
        r.inside += mk_sepvol_xin_(V, NULL, nv, BH + k, hb + k);
        const long c = mk_sepvol_xcount_(V, T, NULL, nt, BH + k, HV[k], HT[k]);
        if (c < 0) goto out;
        r.cross += c;
    }
    mk_bvh_free(&B);
    if (mk_bvh_build(&B, V, T, nt)) goto out;
    if ((r.self1 = mk_sepvol_self_(V, T, nt, &B, 0, NULL)) < 0) goto out;
    for (int i = 0; i < nv; i++) {
        const double *a = V + 3 * (size_t)i, *b = V0 + 3 * (size_t)i;
        r.moved += a[0] != b[0] || a[1] != b[1] || a[2] != b[2];
        r.worst = fmax(r.worst, sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2])));
    }
    r.vol1 = mk_sepvol_volume_(V, T, nt), r.fold1 = mk_sepvol_fold_(V, T, nt, nb);
    ret = 0;
out:
    if (ret && V0) memcpy(V, V0, sizeof *V * 3 * (size_t)nv);
    mk_bvh_free(&B);
    for (int k = 0; BH && k < nh; k++) mk_bvh_free(BH + k);
    if (rep && ret != -1) *rep = r;
    free(BH), free(hb), free(V0), free(N), free(D), free(Vn), free(s), free(nb), free(ch), free(bad), free(buf), free(lvl);
    return ret;
}
#endif
