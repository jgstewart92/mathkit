/* mathkit: a part rebuilt closed and clear of the meshes around it (its room). Header only, C23.
   mk_refit takes a part (the source: any triangle mesh, closed or not, in one or many pieces, crossing itself or not) and
   an array of obstacles, and returns a closed, uncrossed surface as near the source's shape as the room allows: at least
   a clearance off every obstacle, every piece of the source kept (a cavity only on request). Units are the caller's.
   Field. On the lattice h Z^3, in tiles of MK_REFIT_TILE^3 cells, only where the source is (tiles within a band of its
   triangles, and tiles whose centre is inside it): F = max(Fs, max over obstacles of (gap - e)). Fs is the source's
   signed distance d (mk_bvh_nearest within the band; negative inside), and for a piece with open edges |d| - s, a
   shell of half-thickness s round it (an open piece has no inside); a closed piece thinner than `thin` (twice its
   volume over its area) gets min(d, |d| - h), so it stays at least two cells thick. e is the distance to an obstacle,
   negative inside a closed one; an open obstacle has no inside and is kept clear of by |e| alone. Fs's sign comes from
   the source's closed pieces: with MK_REFIT_PARITY, mk_bvh_inside on all of them as one mesh (where two pieces overlap
   the parity is even, so the overlap reads outside and the solid is hollowed there); with MK_REFIT_UNION, inside any
   one piece (mk_bvh_inside on each), so a region two pieces cover stays solid, and a cavity a piece encloses is filled.
   A piece covering a region twice by folding through itself is hollowed there either way. Obstacles' inside is
   mk_bvh_inside's. The same field, cheaper: a brick of 8^3 nodes with no triangle within the band of any of its nodes
   takes the band's value and one inside test for the whole brick (no surface crosses it), for the source and for each
   obstacle; along a row a node's inside state is the node before's when that one lay over h (1 + 1e-6) from every
   triangle (the step between them cannot reach the surface, and the node is not within rounding of it). A node on the surface (Fs 0) is inside where its sign test says so.
   Surface: mk_mtet on each tile (in lattice units, so neighbouring tiles' vertices on their shared faces are the same
   bits), the tiles joined in their order (the same mesh at any thread count) and welded by position; then pieces of
   volume under vmin dropped and inward pieces (voids) dropped unless keep_voids; then mk_decimate (dist dec_h h).
   Contact: mk_refit_check, the same test the retry loop and the caller's delivery use: per obstacle whose box lies
   within `reach` of the part's, face pairs crossing between the two (an edge of one strictly through a face of the
   other, as mk_tri_cross), vertices of each inside the other (only against a closed mesh: an open one has no inside),
   and the least distance from each one's vertices to the other's surface (within reach). In contact when any count is
   over 0 or the gap is under gap_ok. On a contact the build is made again with the gap gap_step wider, up to `tries`
   builds. With a target volume the part grows (Fs - g) into free room toward it, g from 0 to grow_max; of the clear
   builds the largest at most the target is kept (none grows past it, but one without growth may already be).
   Closed or open is the caller's flag on each obstacle, and must be decided on the mesh as given: an edge used once (or
   three times) makes a mesh open. Welding equal positions first pinches a closed mesh with two vertices at one point
   (two branches of a nerve that touch) into one with edges used four times, and an obstacle wrongly read as open loses
   its inside: a part can then be built, and pass the check, lying inside it (in one anatomical model 34 of its 2,359
   parts, most of them nerves, read open that way). mk_refit_open_edges counts edge uses on the mesh as given.
   Grid size: h is hmin, raised by 7 % steps until the band's nodes, about A (2 V / A + 8 h) / (2 h^3) for area A and
   closed volume V, fit in nact; or the caller's h.
   Threads (with -fopenmp): the tiles, the inside test of tile centres and the checker's loops; the result is the same
   bytes at any thread count. Single tile memory: (MK_REFIT_TILE + 1)^3 floats a thread.
   Every MADE UP default is mk_refit_defaults'; each is a choice, not a measured or derived value. */
#ifndef MK_REFIT_H
#define MK_REFIT_H
#include "mk_dec.h"
#include "mk_geom.h"
#include "mk_graft.h"
#include "mk_par.h"
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum { MK_REFIT_PARITY = 0, MK_REFIT_UNION = 1 };
/* status: clear; in contact after every try; not closed or crossing itself (or no volume); nothing left */
enum { MK_REFIT_OK = 0, MK_REFIT_CONTACT = 1, MK_REFIT_BAD = 2, MK_REFIT_EMPTY = 3 };
enum { MK_REFIT_TILE = 64 };

typedef struct {
    double h;        /* > 0: the lattice spacing, no budget */
    double hmin;     /* the finest spacing the budget starts from */
    double nact;     /* the budget: nodes in the band */
    double shell;    /* an open piece's shell: half-thickness max(shell, shell_h h) */
    double shell_h;
    double thin;     /* a closed piece thinner than this (2 V / A) gets a shell of h */
    double gap;      /* the clearance built: gap + gap_h h */
    double gap_h;
    double gap_ok;   /* the checker's least gap */
    double gap_step; /* added to the gap after a try in contact */
    double reach;    /* obstacles within this of the part's box are checked; distances beyond it are not measured */
    int tries;       /* builds at most, without a target */
    int tries_target;/* builds at most, with a target */
    double vmin;     /* pieces of less volume are dropped */
    int keep_voids;  /* 1: inward pieces (cavities) kept */
    double dec_h;    /* mk_decimate's dist, times h; 0: not decimated */
    double maxedge;  /* mk_decimate's longest new edge */
    double mm;       /* one millimetre in the caller's units (mk_dec_defaults' other bounds) */
    double target;   /* > 0: the volume to grow toward, never past */
    double grow_max; /* the most the part may grow */
    double grow_tol; /* the target is met within this share of it */
    int sign;        /* MK_REFIT_PARITY or MK_REFIT_UNION */
} mk_refit_opts;

/* Defaults for a mesh in which one millimetre is mm units (metres: 1e-3). All MADE UP: hmin 0.08 mm, nact 14e6, shell
   0.15 mm and 1.2, thin 1 mm, gap 0.07 mm + 0.6 h, gap_ok 0.05 mm, gap_step 0.05 mm, reach 2 mm, tries 3 (7 with a
   target), vmin 0.5 mm^3, dec_h 0.2, maxedge 1.5 mm, grow_max 1 mm, grow_tol 0.03, voids dropped, MK_REFIT_UNION. */
static inline void mk_refit_defaults(mk_refit_opts *o, double mm)
{
    memset(o, 0, sizeof *o);
    o->hmin = 0.08 * mm, o->nact = 14e6, o->shell = 0.15 * mm, o->shell_h = 1.2, o->thin = 1 * mm;
    o->gap = 0.07 * mm, o->gap_h = 0.6, o->gap_ok = 0.05 * mm, o->gap_step = 0.05 * mm, o->reach = 2 * mm;
    o->tries = 3, o->tries_target = 7, o->vmin = 0.5 * mm * mm * mm, o->dec_h = 0.2, o->maxedge = 1.5 * mm, o->mm = mm;
    o->grow_max = 1 * mm, o->grow_tol = 0.03, o->sign = MK_REFIT_UNION;
}

/* An obstacle (or any mesh to check): its arrays (borrowed, they must outlive it), whether it is closed, its tree and
   box. Made by mk_refit_obj_init, freed by mk_refit_obj_free. */
typedef struct {
    const double *V;
    const int *T;
    int nv, nt, closed;
    mk_bvh B;
    double lo[3], hi[3];
} mk_refit_obj;
static inline void mk_refit_obj_free(mk_refit_obj *m) { mk_bvh_free(&m->B), memset(m, 0, sizeof *m); }
/* Returns 0, or -1 (no face, an index out of range, a coordinate not finite, no memory) with *m empty. */
static inline int mk_refit_obj_init(mk_refit_obj *m, const double *V, int nv, const int *T, int nt, int closed)
{
    memset(m, 0, sizeof *m);
    if (!V || !T || nv < 1 || nt < 1) return -1;
    for (size_t q = 0; q < 3 * (size_t)nt; q++)
        if (T[q] < 0 || T[q] >= nv) return -1;
    for (int k = 0; k < 3; k++) m->lo[k] = INFINITY, m->hi[k] = -INFINITY;
    for (size_t v = 0; v < (size_t)nv; v++)
        for (int k = 0; k < 3; k++) {
            const double x = V[3 * v + k];
            if (!isfinite(x)) return -1;
            m->lo[k] = fmin(m->lo[k], x), m->hi[k] = fmax(m->hi[k], x);
        }
    if (mk_bvh_build(&m->B, V, T, nt)) { memset(m, 0, sizeof *m); return -1; }
    m->V = V, m->T = T, m->nv = nv, m->nt = nt, m->closed = closed != 0;
    return 0;
}

static inline int mk_refit_cmpu64_(const void *x, const void *y)
{
    const uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y;
    return (a > b) - (a < b);
}
/* Edge uses of (T, nt) on edges not used by exactly two faces: 0 for a closed mesh (each edge counted by how often it
   is used, so one used three times adds 3). The mesh as given: two vertices at one point are two vertices. -1 for no
   memory or nt < 0. */
static inline long mk_refit_open_edges(const int *T, int nt)
{
    if (nt < 0 || (nt && !T)) return -1;
    const size_t ne = 3 * (size_t)nt;
    uint64_t *e = malloc(sizeof *e * (ne ? ne : 1));
    if (!e) return -1;
    for (size_t t = 0; t < (size_t)nt; t++)
        for (int k = 0; k < 3; k++) {
            const uint64_t a = (uint32_t)T[3 * t + k], b = (uint32_t)T[3 * t + (k + 1) % 3];
            e[3 * t + k] = a < b ? a << 32 | b : b << 32 | a;
        }
    qsort(e, ne, sizeof *e, mk_refit_cmpu64_);
    long bad = 0;
    for (size_t i = 0; i < ne;) {
        size_t j = i;
        while (j < ne && e[j] == e[i]) j++;
#if defined(MK_MUTATE) && MK_MUTATE == 409
        if (j - i == 1) bad += 1; /* 409: an edge used three times is not open */
#else
        if (j - i != 2) bad += (long)(j - i);
#endif
        i = j;
    }
    free(e);
    return bad;
}

/* One pair's contact: crossing face pairs, a's vertices inside b, b's inside a, the least gap (reach when nothing is
   nearer), contact (1 or 0), tested (0: the boxes lie farther than reach apart, so nothing was measured). */
typedef struct {
    long cross, a_in_b, b_in_a;
    double gap;
    int contact, tested;
} mk_refit_pair;
/* a's vertices within reach of b's box: those inside b (b closed), and the least squared distance of the others */
static inline void mk_refit_side_(const mk_refit_obj *a, const mk_refit_obj *b, double reach, long *in, double *g2)
{
    const long n = a->nv, nb = mk_par_nblk(n);
    long cnt[MK_PAR_NBLK];
    double mn[MK_PAR_NBLK];
    MK_OMP_(omp parallel for schedule(dynamic, 1) if(nb > 1))
    for (long k = 0; k < nb; k++) {
        long i0, i1, c = 0;
        double m = reach * reach;
        mk_par_range(n, k, &i0, &i1);
        for (long i = i0; i < i1; i++) {
            const double *p = a->V + 3 * i;
            int out = 0;
            for (int x = 0; x < 3; x++) out |= p[x] > b->hi[x] + reach || p[x] < b->lo[x] - reach;
            if (out) continue;
#if defined(MK_MUTATE) && MK_MUTATE == 410
            if (mk_bvh_inside(&b->B, p)) { c++; continue; } /* 410: inside an open mesh counted */
#else
            if (b->closed && mk_bvh_inside(&b->B, p)) { c++; continue; }
#endif
            int t;
            double q[3];
            const double d2 = mk_bvh_nearest(&b->B, p, reach * reach, &t, q);
            if (d2 < m) m = d2;
        }
        cnt[k] = c, mn[k] = m;
    }
    *in = 0, *g2 = reach * reach;
    for (long k = 0; k < nb; k++) *in += cnt[k], *g2 = fmin(*g2, mn[k]);
}
/* face X of mesh VA crosses face Y of mesh VB: an edge of one strictly through the other (two meshes share no corner) */
static inline int mk_refit_tri_cross_(const double *VA, const int *X, const double *VB, const int *Y)
{
    const double *a0 = VA + 3 * (size_t)X[0], *a1 = VA + 3 * (size_t)X[1], *a2 = VA + 3 * (size_t)X[2];
    const double *b0 = VB + 3 * (size_t)Y[0], *b1 = VB + 3 * (size_t)Y[1], *b2 = VB + 3 * (size_t)Y[2];
    const double *A[3] = {a0, a1, a2}, *Bq[3] = {b0, b1, b2};
    for (int k = 0; k < 3; k++)
#if defined(MK_MUTATE) && MK_MUTATE == 420
        if (mk_seg_tri_strict(A[k], A[(k + 1) % 3], b0, b1, b2) || (0 && Bq[k])) return 1; /* 420: only a's edges through b's faces */
#else
        if (mk_seg_tri_strict(A[k], A[(k + 1) % 3], b0, b1, b2) || mk_seg_tri_strict(Bq[k], Bq[(k + 1) % 3], a0, a1, a2)) return 1;
#endif
    return 0;
}
/* squared distance between segments p1 q1 and p2 q2 (Ericson 2005, 5.1.9) */
static inline double mk_refit_segseg2_(const double *p1, const double *q1, const double *p2, const double *q2)
{
    double d1[3], d2[3], r[3], s, t;
    for (int k = 0; k < 3; k++) d1[k] = q1[k] - p1[k], d2[k] = q2[k] - p2[k], r[k] = p1[k] - p2[k];
    const double a = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2], e = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2], f = d2[0] * r[0] + d2[1] * r[1] + d2[2] * r[2];
    if (!(a > 0) && !(e > 0)) s = t = 0;
    else if (!(a > 0)) s = 0, t = fmin(1, fmax(0, f / e));
    else {
        const double c = d1[0] * r[0] + d1[1] * r[1] + d1[2] * r[2];
        if (!(e > 0)) t = 0, s = fmin(1, fmax(0, -c / a));
        else {
            const double b = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2], dn = a * e - b * b;
            s = dn > 0 ? fmin(1, fmax(0, (b * f - c * e) / dn)) : 0;
            t = (b * s + f) / e;
            if (t < 0) t = 0, s = fmin(1, fmax(0, -c / a));
            else if (t > 1) t = 1, s = fmin(1, fmax(0, (b - c) / a));
        }
    }
    double x = 0;
    for (int k = 0; k < 3; k++) {
        const double y = p1[k] + d1[k] * s - p2[k] - d2[k] * t;
        x += y * y;
    }
    return x;
}
/* squared distance between faces X of VA and Y of VB that do not cross: the least of corner to face and edge to edge */
static inline double mk_refit_facegap2_(const double *VA, const int *X, const double *VB, const int *Y)
{
    const double *A[3] = {VA + 3 * (size_t)X[0], VA + 3 * (size_t)X[1], VA + 3 * (size_t)X[2]}, *B[3] = {VB + 3 * (size_t)Y[0], VB + 3 * (size_t)Y[1], VB + 3 * (size_t)Y[2]};
    double m = INFINITY, w[3], q[3];
    for (int k = 0; k < 3; k++) {
        m = fmin(m, mk_closest_pt_tri(A[k], B[0], B[1], B[2], w, q)), m = fmin(m, mk_closest_pt_tri(B[k], A[0], A[1], A[2], w, q));
        for (int j = 0; j < 3; j++) m = fmin(m, mk_refit_segseg2_(A[k], A[(k + 1) % 3], B[j], B[(j + 1) % 3]));
    }
    return m;
}
/* the crossing pairs between a and b, and (gap_ok > 0) the least squared face gap of pairs within gap_ok (*fg2, else
   untouched); -1 for no memory */
static inline long mk_refit_cross_(const mk_refit_obj *a, const mk_refit_obj *b, double gap_ok, double *fg2)
{
    const long n = a->nt, nb = mk_par_nblk(n);
    long cnt[MK_PAR_NBLK];
    double fmn[MK_PAR_NBLK];
    int fail = 0;
    MK_OMP_(omp parallel for schedule(dynamic, 1) if(nb > 1))
    for (long k = 0; k < nb; k++) {
        long i0, i1, c = 0;
        double fm = INFINITY;
        int cap = 64, *buf = malloc(sizeof *buf * (size_t)cap);
        mk_par_range(n, k, &i0, &i1);
        cnt[k] = 0, fmn[k] = INFINITY;
        if (!buf) {
            MK_OMP_(omp atomic write)
            fail = 1;
            continue;
        }
        for (long t = i0; t < i1; t++) {
            mk_aabb bx;
            mk_aabb_empty(&bx);
            for (int j = 0; j < 3; j++) mk_aabb_add(&bx, a->V + 3 * (size_t)a->T[3 * t + j]);
            for (int x = 0; x < 3; x++) bx.lo[x] -= gap_ok, bx.hi[x] += gap_ok;
            int meet = 1;
            for (int x = 0; x < 3; x++) meet &= bx.lo[x] <= b->hi[x] && b->lo[x] <= bx.hi[x];
            if (!meet) continue;
            int m = mk_bvh_box(&b->B, &bx, buf, cap);
            if (m > cap) {
                int *b2 = realloc(buf, sizeof *b2 * (size_t)m);
                if (!b2) {
                    MK_OMP_(omp atomic write)
                    fail = 1;
                    break;
                }
                buf = b2, cap = m, m = mk_bvh_box(&b->B, &bx, buf, cap);
            }
            for (int i = 0; i < m; i++) {
                const int x = mk_refit_tri_cross_(a->V, a->T + 3 * t, b->V, b->T + 3 * (size_t)buf[i]);
                c += x;
                if (!x && gap_ok > 0) fm = fmin(fm, mk_refit_facegap2_(a->V, a->T + 3 * t, b->V, b->T + 3 * (size_t)buf[i]));
            }
        }
        free(buf);
        cnt[k] = c, fmn[k] = fm;
    }
    if (fail) return -1;
    long s = 0;
    for (long k = 0; k < nb; k++) s += cnt[k];
    if (fg2)
        for (long k = 0; k < nb; k++) *fg2 = fmin(*fg2, fmn[k]);
    return s;
}
/* The contact test of mesh a against mesh b (each from mk_refit_obj_init), under o's reach and gap_ok (o NULL:
   mk_refit_defaults in metres). Returns 1 in contact, 0 clear (or too far apart to test), -1 for bad arguments or no
   memory, with *r filled. "Inside" counts only against a closed mesh. */
static inline int mk_refit_check(const mk_refit_obj *a, const mk_refit_obj *b, const mk_refit_opts *o, mk_refit_pair *r)
{
    mk_refit_opts od;
    if (!o) mk_refit_defaults(&od, 1e-3), o = &od;
    if (!r) return -1;
    memset(r, 0, sizeof *r);
    r->gap = o->reach;
    if (!a || !b || !a->V || !b->V || !(o->reach >= 0) || !(o->gap_ok >= 0)) return -1;
#if defined(MK_MUTATE) && MK_MUTATE == 421
    for (int k = 0; k < 3; k++)
        if (b->lo[k] > a->hi[k] || a->lo[k] > b->hi[k]) return 0; /* 421: boxes apart by less than reach not tested */
#else
    for (int k = 0; k < 3; k++)
        if (b->lo[k] > a->hi[k] + o->reach || a->lo[k] > b->hi[k] + o->reach) return 0;
#endif
    r->tested = 1;
    double ga, gb;
    mk_refit_side_(a, b, o->reach, &r->a_in_b, &ga);
    mk_refit_side_(b, a, o->reach, &r->b_in_a, &gb);
    double fg = o->reach * o->reach;
    r->cross = mk_refit_cross_(a, b, o->gap_ok, &fg);
    if (r->cross < 0) return -1;
#if defined(MK_MUTATE) && MK_MUTATE == 424
    r->gap = sqrt(fmin(ga, gb)); /* 424: vertex to surface only, so near faces with no vertex near read clear */
#else
    r->gap = sqrt(fmin(fmin(ga, gb), fg));
#endif
#if defined(MK_MUTATE) && MK_MUTATE == 411
    r->contact = r->cross > 0 || r->a_in_b > 0 || r->b_in_a > 0; /* 411: the gap not looked at */
#else
    r->contact = r->cross > 0 || r->a_in_b > 0 || r->b_in_a > 0 || r->gap < o->gap_ok;
#endif
    return r->contact;
}

/* What a refit did. pieces_in: the source's pieces (faces joined through shared vertex indices), open_in of them with
   open edges; vol_in: the closed pieces' volume; pieces_out: the result's pieces, dropped: the volume of those dropped
   under vmin; lost: source pieces with no part of the result within 3 (h + gap) of any of their vertices; obstacles:
   those within reach of the source's box; contacts and least_gap (reach if nothing nearer): the result's against the
   room, worst the obstacle of the last contact, or of the least gap (-1: none); h, gap, grow: as built; tries: builds
   made; tiles and nodes: of the kept build; open_out, cross_out: the result's open edge uses and crossing pairs;
   status: MK_REFIT_*. */
typedef struct {
    int pieces_in, open_in, pieces_out, lost, obstacles, contacts, worst, tries, status;
    double vol_in, vol_out, dropped, least_gap, h, gap, grow;
    long tiles, nodes, open_out, cross_out;
} mk_refit_report;

/* the source taken apart: piece of each vertex and face, each piece's open edges, signed volume, area; closed pieces'
   faces first in Tc, a tree over them (PARITY) or one each (UNION) */
typedef struct {
    int nc, ncl, *vlab, *flab, *Tc, *cstart, *ccount;
    long *co;
    double *cv, *ca, *s, (*cbox)[6];
    mk_bvh B, CB, *PB;
    int sign;
} mk_refit_src_;
static inline void mk_refit_src_free_(mk_refit_src_ *S)
{
    if (S->PB)
        for (int c = 0; c < S->ncl; c++) mk_bvh_free(&S->PB[c]);
    mk_bvh_free(&S->B), mk_bvh_free(&S->CB);
    free(S->vlab), free(S->flab), free(S->Tc), free(S->cstart), free(S->ccount), free(S->co), free(S->cv), free(S->ca), free(S->s), free(S->cbox), free(S->PB);
    memset(S, 0, sizeof *S);
}
static inline int mk_refit_uf_(int *p, int x)
{
    while (p[x] != x) p[x] = p[p[x]], x = p[x];
    return x;
}
static inline int mk_refit_insrc_(const mk_refit_src_ *S, const double *p)
{
    if (!S->ncl) return 0;
#if defined(MK_MUTATE) && MK_MUTATE == 412
    if (0) { /* 412: UNION read as PARITY: two pieces' overlap hollow */
#else
    if (S->sign == MK_REFIT_UNION) {
#endif
        for (int c = 0; c < S->ncl; c++) {
            const double *b = S->cbox[c];
            if (p[0] < b[0] || p[1] < b[1] || p[2] < b[2] || p[0] > b[3] || p[1] > b[4] || p[2] > b[5]) continue;
            if (mk_bvh_inside(&S->PB[c], p)) return 1;
        }
        return 0;
    }
    return mk_bvh_inside(&S->CB, p);
}
/* the source's pieces as given (the shells are set once h is known); -1 for no memory */
static inline int mk_refit_src_init_(mk_refit_src_ *S, const double *V, int nv, const int *T, int nt, const double *lo, int sign)
{
    memset(S, 0, sizeof *S);
    S->sign = sign;
    int *par = malloc(sizeof *par * (size_t)nv), *id = malloc(sizeof *id * (size_t)nv);
    S->vlab = malloc(sizeof *S->vlab * (size_t)nv), S->flab = malloc(sizeof *S->flab * (size_t)nt);
    uint64_t *e = malloc(sizeof *e * 3 * (size_t)nt);
    if (!par || !id || !S->vlab || !S->flab || !e) goto fail;
    for (int i = 0; i < nv; i++) par[i] = i, id[i] = -1;
    for (int t = 0; t < nt; t++) {
        const int a = mk_refit_uf_(par, T[3 * (size_t)t]);
        for (int k = 1; k < 3; k++) {
            const int b = mk_refit_uf_(par, T[3 * (size_t)t + k]);
            if (a != b) par[b] = a;
        }
    }
    for (int t = 0; t < nt; t++) {
        const int r = mk_refit_uf_(par, T[3 * (size_t)t]);
        if (id[r] < 0) id[r] = S->nc++;
        S->flab[t] = id[r];
    }
    for (int v = 0; v < nv; v++) S->vlab[v] = id[mk_refit_uf_(par, v)];
    const int nc = S->nc;
    S->co = calloc((size_t)nc, sizeof *S->co), S->cv = calloc((size_t)nc, sizeof *S->cv), S->ca = calloc((size_t)nc, sizeof *S->ca);
    S->s = calloc((size_t)nc, sizeof *S->s), S->cstart = calloc((size_t)nc + 1, sizeof *S->cstart), S->ccount = calloc((size_t)nc, sizeof *S->ccount);
    S->Tc = malloc(sizeof *S->Tc * 3 * (size_t)nt), S->cbox = malloc(sizeof *S->cbox * (size_t)nc);
    if (!S->co || !S->cv || !S->ca || !S->s || !S->cstart || !S->ccount || !S->Tc || !S->cbox) goto fail;
    for (size_t t = 0; t < (size_t)nt; t++)
        for (int k = 0; k < 3; k++) {
            const uint64_t a = (uint32_t)T[3 * t + k], b = (uint32_t)T[3 * t + (k + 1) % 3];
            e[3 * t + k] = a < b ? a << 32 | b : b << 32 | a;
        }
    qsort(e, 3 * (size_t)nt, sizeof *e, mk_refit_cmpu64_);
    for (size_t i = 0, ne = 3 * (size_t)nt; i < ne;) {
        size_t j = i;
        while (j < ne && e[j] == e[i]) j++;
        if (j - i != 2) S->co[S->vlab[e[i] >> 32]] += (long)(j - i);
        i = j;
    }
    for (int t = 0; t < nt; t++) {
        double a[3], b[3], c[3];
        for (int k = 0; k < 3; k++)
            a[k] = V[3 * (size_t)T[3 * (size_t)t] + k] - lo[k], b[k] = V[3 * (size_t)T[3 * (size_t)t + 1] + k] - lo[k], c[k] = V[3 * (size_t)T[3 * (size_t)t + 2] + k] - lo[k];
        const double v6 = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double x = u[1] * w[2] - u[2] * w[1], y = u[2] * w[0] - u[0] * w[2], z = u[0] * w[1] - u[1] * w[0];
        S->cv[S->flab[t]] += v6 / 6, S->ca[S->flab[t]] += 0.5 * sqrt(x * x + y * y + z * z);
    }
    /* closed pieces first (their old order), then the open ones */
    int *rank = id; /* reused: the new number of each piece */
    for (int c = 0, r = 0; c < nc; c++)
        if (!S->co[c]) rank[c] = r++;
    S->ncl = 0;
    for (int c = 0; c < nc; c++) S->ncl += !S->co[c];
    for (int c = 0, r = S->ncl; c < nc; c++)
        if (S->co[c]) rank[c] = r++;
    for (int t = 0; t < nt; t++) S->ccount[rank[S->flab[t]]]++;
    for (int c = 0; c < nc; c++) S->cstart[c + 1] = S->cstart[c] + S->ccount[c];
    for (int c = 0; c < nc; c++) par[c] = S->cstart[c]; /* par reused: each piece's next slot */
    for (int t = 0; t < nt; t++) memcpy(S->Tc + 3 * (size_t)par[rank[S->flab[t]]]++, T + 3 * (size_t)t, 3 * sizeof *T);
    for (int c = 0; c < nc; c++)
        for (int k = 0; k < 3; k++) S->cbox[c][k] = INFINITY, S->cbox[c][3 + k] = -INFINITY;
    for (int v = 0; v < nv; v++) {
        if (S->vlab[v] < 0) continue; /* a vertex of no face */
        double *b = S->cbox[rank[S->vlab[v]]];
        for (int k = 0; k < 3; k++) b[k] = fmin(b[k], V[3 * (size_t)v + k]), b[3 + k] = fmax(b[3 + k], V[3 * (size_t)v + k]);
    }
    if (mk_bvh_build(&S->B, V, T, nt)) goto fail;
    if (S->ncl) {
        if (S->sign == MK_REFIT_UNION) {
            S->PB = calloc((size_t)S->ncl, sizeof *S->PB);
            if (!S->PB) goto fail;
            for (int c = 0; c < S->ncl; c++)
                if (mk_bvh_build(&S->PB[c], V, S->Tc + 3 * (size_t)S->cstart[c], S->ccount[c])) goto fail;
        } else if (mk_bvh_build(&S->CB, V, S->Tc, S->cstart[S->ncl]))
            goto fail;
    }
    free(par), free(id), free(e);
    return 0;
fail:
    free(par), free(id), free(e), mk_refit_src_free_(S);
    return -1;
}

/* the welded surface of the tiles: a hash on the three coordinates' bits */
static inline uint64_t mk_refit_hash_(const double *p)
{
    uint64_t k[3];
    memcpy(k, p, sizeof k);
    uint64_t h = k[0] * 0x9E3779B97F4A7C15u ^ k[1] * 0xC2B2AE3D27D4EB4Fu ^ k[2] * 0x165667B19E3779F9u;
    h ^= h >> 29, h *= 0xBF58476D1CE4E5B9u, h ^= h >> 32;
    return h;
}
/* (V, nv, T, nt) welded in place by exact position (the vertices renumbered in their order), faces with a repeated
   corner dropped */
static inline int mk_refit_weld_(double *V, int *nv, int *T, int *nt)
{
    size_t cap = 16;
    while (cap < 2 * (size_t)*nv + 16) cap *= 2;
    int *tab = malloc(sizeof *tab * cap), *map = malloc(sizeof *map * ((size_t)*nv + 1)), n = 0;
    if (!tab || !map) { free(tab), free(map); return -1; }
    memset(tab, 255, sizeof *tab * cap);
    for (int v = 0; v < *nv; v++) {
        const double *p = V + 3 * (size_t)v;
        size_t s = mk_refit_hash_(p) & (cap - 1);
        while (tab[s] >= 0 && memcmp(V + 3 * (size_t)tab[s], p, 3 * sizeof *p)) s = (s + 1) & (cap - 1);
        if (tab[s] < 0) {
            memmove(V + 3 * (size_t)n, p, 3 * sizeof *p);
            tab[s] = n++;
        }
        map[v] = tab[s];
    }
    int t2 = 0;
    for (int t = 0; t < *nt; t++) {
        const int a = map[T[3 * (size_t)t]], b = map[T[3 * (size_t)t + 1]], c = map[T[3 * (size_t)t + 2]];
        if (a == b || b == c || a == c) continue;
        T[3 * (size_t)t2] = a, T[3 * (size_t)t2 + 1] = b, T[3 * (size_t)t2 + 2] = c, t2++;
    }
    free(tab), free(map);
    *nv = n, *nt = t2;
    return 0;
}
/* keep only the faces with keep[t]; vertices renumbered in their order */
static inline int mk_refit_compact_(double *V, int *nv, int *T, int *nt, const unsigned char *keep)
{
    int *map = malloc(sizeof *map * ((size_t)*nv + 1)), t2 = 0, n = 0;
    if (!map) return -1;
    for (int t = 0; t < *nt; t++)
        if (!keep || keep[t]) memmove(T + 3 * (size_t)t2++, T + 3 * (size_t)t, 3 * sizeof *T);
    memset(map, 255, sizeof *map * ((size_t)*nv + 1));
    for (size_t x = 0; x < 3 * (size_t)t2; x++) map[T[x]] = 0;
    for (int v = 0; v < *nv; v++)
        if (!map[v]) map[v] = n, memmove(V + 3 * (size_t)n, V + 3 * (size_t)v, 3 * sizeof *V), n++;
    for (size_t x = 0; x < 3 * (size_t)t2; x++) T[x] = map[T[x]];
    free(map);
    *nv = n, *nt = t2;
    return 0;
}
static inline double mk_refit_vol6_(const double *V, const int *T, int t, const double *o)
{
    double a[3], b[3], c[3];
    for (int k = 0; k < 3; k++)
        a[k] = V[3 * (size_t)T[3 * (size_t)t] + k] - o[k], b[k] = V[3 * (size_t)T[3 * (size_t)t + 1] + k] - o[k], c[k] = V[3 * (size_t)T[3 * (size_t)t + 2] + k] - o[k];
    return a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
}

/* one build: the field on the tiles, the surface, pieces, decimation; out empty when nothing is solid. -1: no memory */
static inline int mk_refit_build_(const double *V, const int *T, int nt, const mk_refit_src_ *S, const mk_refit_obj *room, int nroom, const mk_refit_opts *o,
                                  const double *lo, const double *hi, double h, double gap, double grow, mk_gmesh *out, long *ntiles, long *nnodes, double *dropped)
{
    enum { TC = MK_REFIT_TILE, N1 = MK_REFIT_TILE + 1, NB = (MK_REFIT_TILE + 8) / 8 };
    memset(out, 0, sizeof *out);
    *ntiles = *nnodes = 0, *dropped = 0;
    /* hc: a node's inside state is carried to the next only from over hc off every triangle, so the next lies at least
       1e-6 h off it, where the sign test cannot be turned by rounding: a node within rounding of a face is tested on its
       own, the same in each tile that holds it */
#if defined(MK_MUTATE) && MK_MUTATE == 423
    const double band = 3 * h + gap + grow, pad = band + 2 * h, hc = h; /* 423: carried right up to h */
#else
    const double band = 3 * h + gap + grow, pad = band + 2 * h, hc = h * (1 + 1e-6);
#endif
    long L0[3], Lt[3];
    for (int k = 0; k < 3; k++) {
        L0[k] = (long)floor((lo[k] - pad) / h), Lt[k] = ((long)ceil((hi[k] + pad) / h) - L0[k]) / TC + 1;
        if (Lt[k] < 1 || Lt[k] > (1L << 20)) return -1;
    }
    const long ntile = Lt[0] * Lt[1] * Lt[2];
    unsigned char *use = calloc((size_t)ntile, 1);
    if (!use) return -1;
    for (int t = 0; t < nt; t++) {
        long a[3], b[3];
        for (int k = 0; k < 3; k++) {
            double x = INFINITY, y = -INFINITY;
            for (int j = 0; j < 3; j++) {
                const double z = V[3 * (size_t)T[3 * (size_t)t + j] + k];
                x = fmin(x, z), y = fmax(y, z);
            }
            a[k] = ((long)floor((x - pad) / h) - L0[k]) / TC, b[k] = ((long)ceil((y + pad) / h) - L0[k]) / TC;
            if (a[k] < 0) a[k] = 0;
            if (b[k] >= Lt[k]) b[k] = Lt[k] - 1;
        }
        for (long z = a[2]; z <= b[2]; z++)
            for (long y = a[1]; y <= b[1]; y++)
                for (long x = a[0]; x <= b[0]; x++) use[x + Lt[0] * (y + Lt[1] * z)] = 1;
    }
    if (S->ncl) {
        MK_OMP_(omp parallel for schedule(dynamic, 16))
        for (long i = 0; i < ntile; i++) {
            if (use[i]) continue;
            const long t3[3] = {i % Lt[0], i / Lt[0] % Lt[1], i / Lt[0] / Lt[1]};
            double p[3];
            for (int k = 0; k < 3; k++) p[k] = ((double)(L0[k] + t3[k] * TC) + 0.5 * TC) * h;
            if (mk_refit_insrc_(S, p)) use[i] = 2;
        }
    }
    long nu = 0;
    for (long i = 0; i < ntile; i++) nu += use[i] != 0;
    long *ul = malloc(sizeof *ul * (size_t)(nu ? nu : 1));
    mk_gmesh *sl = calloc((size_t)(nu ? nu : 1), sizeof *sl);
    long *nn = calloc((size_t)(nu ? nu : 1), sizeof *nn);
    if (!ul || !sl || !nn) { free(use), free(ul), free(sl), free(nn); return -1; }
    for (long i = 0, j = 0; i < ntile; i++)
        if (use[i]) ul[j++] = i;
    free(use);
    int fail = 0;
    MK_OMP_(omp parallel)
    {
        float *G = malloc(sizeof *G * (size_t)N1 * N1 * N1);
        int *tn = malloc(sizeof *tn * (size_t)(nroom + 1)), *wl = malloc(sizeof *wl * (size_t)(nroom + 1));
        signed char *bw = malloc((size_t)NB * NB * NB * (size_t)(nroom + 1)), *wc = malloc((size_t)(nroom + 1));
        MK_OMP_(omp for schedule(dynamic, 1))
        for (long u = 0; u < nu; u++) {
            if (!G || !tn || !wl || !bw || !wc) {
                MK_OMP_(omp atomic write)
                fail = 1;
                continue;
            }
            const long i = ul[u], t3[3] = {i % Lt[0], i / Lt[0] % Lt[1], i / Lt[0] / Lt[1]};
            long og[3];
            double tl[3], th[3];
            int ntn = 0, any = 0;
            for (int k = 0; k < 3; k++) og[k] = L0[k] + t3[k] * TC, tl[k] = (double)og[k] * h, th[k] = (double)(og[k] + TC) * h;
            for (int j = 0; j < nroom; j++) {
                int m = 1;
                for (int k = 0; k < 3; k++) m &= room[j].lo[k] - band <= th[k] && tl[k] <= room[j].hi[k] + band;
                if (m) tn[ntn++] = j;
            }
            /* bricks of 8^3 nodes: one with no triangle within the band of any of its nodes is far as a whole, with one
               inside state (no surface crosses it); its nodes then take the value a node-by-node query would give */
            signed char sb[NB * NB * NB];
            double bc[NB * NB * NB][4];
            for (int z = 0; z < NB; z++)
                for (int y = 0; y < NB; y++)
                    for (int x = 0; x < NB; x++) {
                        const int bi = x + NB * (y + NB * z), xyz[3] = {x, y, z};
                        double r2 = 0, *c = bc[bi];
                        for (int k = 0; k < 3; k++) {
                            const int n0 = 8 * xyz[k], nk = N1 - n0 < 8 ? N1 - n0 : 8;
                            const double ex = 0.5 * (nk - 1);
                            c[k] = ((double)(og[k] + n0) + ex) * h, r2 += ex * ex;
                        }
                        c[3] = sqrt(r2) * h * (1 + 1e-9) + 1e-300;
                        int t;
                        double q[3];
                        const double rr = band + c[3];
                        (void)mk_bvh_nearest(&S->B, c, rr * rr, &t, q);
#if defined(MK_MUTATE) && MK_MUTATE == 414
                        sb[bi] = t >= 0 ? -1 : 0; /* 414: a far brick read as outside */
#else
                        sb[bi] = t >= 0 ? -1 : (signed char)mk_refit_insrc_(S, c);
#endif
                    }
            memset(bw, -2, (size_t)NB * NB * NB * (size_t)(ntn ? ntn : 1));
            for (int c = 0; c < N1; c++)
                for (int b = 0; b < N1; b++) {
                    int sc = -1; /* the source's inside state carried from the node before, when that one lay over h from it */
                    for (int j = 0; j < ntn; j++) wl[j] = -2;
                    for (int a = 0; a < N1; a++) {
                        const double p[3] = {(double)(og[0] + a) * h, (double)(og[1] + b) * h, (double)(og[2] + c) * h};
                        const int bi = (a >> 3) + NB * ((b >> 3) + NB * (c >> 3));
                        double q[3], d;
                        int t;
                        if (sb[bi] >= 0) d = sb[bi] ? -band : band, sc = sb[bi];
                        else {
                            const double d2 = mk_bvh_nearest(&S->B, p, band * band, &t, q);
#if defined(MK_MUTATE) && MK_MUTATE == 414
                            const int in = t < 0 ? 0 : sc >= 0 ? sc : mk_refit_insrc_(S, p); /* 414: a far node read as outside */
#else
                            const int in = sc >= 0 ? sc : mk_refit_insrc_(S, p);
#endif
                            if (t < 0) d = in ? -band : band, sc = in;
                            else {
                                const double dd = sqrt(d2), s = S->s[S->flab[t]], sd = in ? -dd : dd;
                                d = (s > 0 ? fmin(sd, dd - s) : sd) - grow, sc = dd > hc ? in : -1;
                            }
                        }
                        if (d < 2 * h)
                            for (int j = 0; j < ntn; j++) {
                                const mk_refit_obj *w = room + tn[j];
                                int outb = 0;
                                for (int k = 0; k < 3; k++) outb |= p[k] < w->lo[k] - band || p[k] > w->hi[k] + band;
                                if (outb) continue;
                                signed char *st = bw + (size_t)bi * (size_t)ntn + (size_t)j;
                                if (*st == -2) {
                                    int tt;
                                    const double rr = band + bc[bi][3];
                                    (void)mk_bvh_nearest(&w->B, bc[bi], rr * rr, &tt, q);
#if defined(MK_MUTATE) && MK_MUTATE == 415
                                    *st = tt >= 0 ? -1 : (signed char)mk_bvh_inside(&w->B, bc[bi]); /* 415: an open obstacle given an inside */
#else
                                    *st = tt >= 0 ? -1 : w->closed ? (signed char)mk_bvh_inside(&w->B, bc[bi]) : 0;
#endif
                                }
                                double e;
                                if (*st >= 0) e = *st ? -band : band, wc[j] = *st;
                                else {
                                    int tt;
                                    const double e2 = mk_bvh_nearest(&w->B, p, band * band, &tt, q);
                                    e = tt >= 0 ? sqrt(e2) : band;
#if defined(MK_MUTATE) && MK_MUTATE == 415
                                    if (1) {
#else
                                    if (w->closed) {
#endif
                                        const int in = wl[j] == a - 1 && wc[j] >= 0 ? wc[j] : mk_bvh_inside(&w->B, p);
                                        wc[j] = e > hc ? (signed char)in : -1;
                                        if (in) e = -e;
                                    } else
                                        wc[j] = 0;
                                }
                                wl[j] = a;
                                d = fmax(d, gap - e);
                            }
                        float g = (float)(fmax(-band, fmin(band, d)) / h);
                        if (g == 0) g = signbit(d) ? -1e-6f : 1e-6f; /* a node on the surface: inside where its sign says so */
                        G[(size_t)a + (size_t)N1 * ((size_t)b + (size_t)N1 * (size_t)c)] = g;
                        any |= g < 0;
                    }
                }
            nn[u] = (long)N1 * N1 * N1;
            if (!any) continue;
            const int nd[3] = {N1, N1, N1};
            const double z3[3] = {0, 0, 0};
            if (mk_mtet(G, nd, z3, 1.0, &sl[u])) {
                MK_OMP_(omp atomic write)
                fail = 1;
                continue;
            }
            for (size_t x = 0; x < (size_t)sl[u].nv; x++)
                for (int k = 0; k < 3; k++) sl[u].V[3 * x + k] += (double)(t3[k] * TC); /* lattice units from L0, the same bits in each tile */
        }
        free(G), free(tn), free(wl), free(bw), free(wc);
    }
    long tv = 0, tt = 0;
    for (long u = 0; u < nu; u++) tv += sl[u].nv, tt += sl[u].nt, *nnodes += nn[u];
    *ntiles = nu;
    if (!fail && (tv > INT_MAX / 2 || tt > INT_MAX / 4)) fail = 1;
    double *RV = NULL;
    int *RT = NULL;
    if (!fail && tt) {
        RV = malloc(sizeof *RV * 3 * (size_t)tv), RT = malloc(sizeof *RT * 3 * (size_t)tt);
        if (!RV || !RT) fail = 1;
    }
    int rnv = 0, rnt = 0;
    for (long u0 = 0; u0 < nu; u0++) {
#if defined(MK_MUTATE) && MK_MUTATE == 419
        const long u = MK_PAR_NTHREADS_() > 1 ? nu - 1 - u0 : u0; /* 419: the tiles joined in an order that follows the thread count */
#else
        const long u = u0;
#endif
        if (!fail && sl[u].nt) {
            memcpy(RV + 3 * (size_t)rnv, sl[u].V, sizeof *RV * 3 * (size_t)sl[u].nv);
            for (size_t x = 0; x < 3 * (size_t)sl[u].nt; x++) RT[3 * (size_t)rnt + x] = sl[u].T[x] + rnv;
            rnv += sl[u].nv, rnt += sl[u].nt;
        }
        mk_gmesh_free(&sl[u]);
    }
    free(sl), free(ul), free(nn);
    if (fail) { free(RV), free(RT); return -1; }
    if (!rnt) { free(RV), free(RT); return 0; }
#if !(defined(MK_MUTATE) && MK_MUTATE == 416) /* 416: the tiles not welded: the surface open on every tile face */
    if (mk_refit_weld_(RV, &rnv, RT, &rnt)) { free(RV), free(RT); return -1; }
#endif
    for (size_t x = 0; x < (size_t)rnv; x++)
        for (int k = 0; k < 3; k++) RV[3 * x + k] = ((double)L0[k] + RV[3 * x + k]) * h;
    /* pieces: kept over vmin, voids dropped unless asked */
    int *par = malloc(sizeof *par * (size_t)rnv), *pid = malloc(sizeof *pid * (size_t)rnv), np = 0;
    unsigned char *keep = malloc((size_t)rnt + 1);
    if (!par || !pid || !keep) { free(par), free(pid), free(keep), free(RV), free(RT); return -1; }
    for (int v = 0; v < rnv; v++) par[v] = v, pid[v] = -1;
    for (int t = 0; t < rnt; t++) {
        const int a = mk_refit_uf_(par, RT[3 * (size_t)t]);
        for (int k = 1; k < 3; k++) {
            const int b = mk_refit_uf_(par, RT[3 * (size_t)t + k]);
            if (a != b) par[b] = a;
        }
    }
    for (int t = 0; t < rnt; t++) {
        const int r = mk_refit_uf_(par, RT[3 * (size_t)t]);
        if (pid[r] < 0) pid[r] = np++;
    }
    double *pv = calloc((size_t)np + 1, sizeof *pv);
    if (!pv) { free(par), free(pid), free(keep), free(RV), free(RT); return -1; }
    for (int t = 0; t < rnt; t++) pv[pid[mk_refit_uf_(par, RT[3 * (size_t)t])]] += mk_refit_vol6_(RV, RT, t, lo) / 6;
    for (int c = 0; c < np; c++)
        if (pv[c] > 0 && pv[c] < o->vmin) *dropped += pv[c];
    for (int t = 0; t < rnt; t++) {
        const double v = pv[pid[mk_refit_uf_(par, RT[3 * (size_t)t])]];
#if defined(MK_MUTATE) && MK_MUTATE == 417
        keep[t] = v >= o->vmin || v < 0; /* 417: voids kept */
#else
        keep[t] = v >= o->vmin || (v < 0 && o->keep_voids);
#endif
    }
    free(par), free(pid), free(pv);
    if (mk_refit_compact_(RV, &rnv, RT, &rnt, keep)) { free(keep), free(RV), free(RT); return -1; }
    free(keep);
    if (!rnt) { free(RV), free(RT); return 0; }
    if (o->dec_h > 0) {
        mk_dec_opts op;
        mk_dec_defaults(&op, o->mm);
        op.dist = o->dec_h * h, op.maxedge = o->maxedge;
        unsigned char *al = malloc((size_t)rnv);
        int *To = malloc(sizeof *To * 3 * (size_t)rnt), nto = 0;
        if (!al || !To) { free(al), free(To), free(RV), free(RT); return -1; }
        if (mk_decimate(RV, rnv, RT, rnt, &op, al, To, &nto, NULL) == 0) {
            memcpy(RT, To, sizeof *To * 3 * (size_t)nto), rnt = nto;
            if (mk_refit_compact_(RV, &rnv, RT, &rnt, NULL)) { free(al), free(To), free(RV), free(RT); return -1; }
        }
        free(al), free(To);
    }
    out->V = RV, out->T = RT, out->nv = rnv, out->nt = rnt;
    return 0;
}

/* Refits the source (V[3 nv], T[3 nt]) clear of room[0 .. nroom) under o (NULL: mk_refit_defaults in metres) into out
   (malloc'd, free with mk_gmesh_free; empty when nothing was built). Returns rep->status (MK_REFIT_OK, MK_REFIT_CONTACT,
   MK_REFIT_BAD, MK_REFIT_EMPTY), or -1 for bad arguments (no face, an index out of range, a coordinate not finite, an
   option out of range) or no memory, with out empty. rep may be NULL. */
static inline int mk_refit(const double *V, int nv, const int *T, int nt, const mk_refit_obj *room, int nroom, const mk_refit_opts *o, mk_gmesh *out,
                           mk_refit_report *rep)
{
    mk_refit_opts od;
    mk_refit_report R;
    memset(&R, 0, sizeof R), R.worst = -1, R.status = MK_REFIT_EMPTY;
    if (rep) *rep = R;
    if (!out) return -1;
    memset(out, 0, sizeof *out);
    if (!o) mk_refit_defaults(&od, 1e-3), o = &od;
    if (!V || !T || nv < 1 || nt < 1 || nv > INT_MAX / 4 || nt > INT_MAX / 4 || nroom < 0 || (nroom && !room)) return -1;
    if (!(o->h >= 0) || !(o->hmin > 0) || !(o->nact > 0) || !(o->shell >= 0) || !(o->shell_h >= 0) || !(o->thin >= 0) || !(o->gap >= 0) || !(o->gap_h >= 0) ||
        !(o->gap_ok >= 0) || !(o->gap_step >= 0) || !(o->reach >= 0) || o->tries < 1 || o->tries_target < 1 || !(o->vmin >= 0) || !(o->dec_h >= 0) ||
        !(o->maxedge >= 0) || !(o->mm > 0) || !(o->target >= 0) || !(o->grow_max >= 0) || !(o->grow_tol >= 0) ||
        (o->sign != MK_REFIT_PARITY && o->sign != MK_REFIT_UNION) || !isfinite(o->h) || !isfinite(o->hmin) || !isfinite(o->gap) || !isfinite(o->target) ||
        !isfinite(o->grow_max) || !isfinite(o->reach))
        return -1;
    for (size_t q = 0; q < 3 * (size_t)nt; q++)
        if (T[q] < 0 || T[q] >= nv) return -1;
    double lo[3], hi[3];
    for (int k = 0; k < 3; k++) lo[k] = INFINITY, hi[k] = -INFINITY;
    for (size_t v = 0; v < (size_t)nv; v++)
        for (int k = 0; k < 3; k++)
            if (!isfinite(V[3 * v + k])) return -1;
#if defined(MK_MUTATE) && MK_MUTATE == 426
    for (size_t q = 0; q < 3 * (size_t)nv; q++) /* 426: the box of every vertex, used or not */
#else
    for (size_t q = 0; q < 3 * (size_t)nt; q++) /* the box of the vertices faces use: an unused vertex far off costs nothing */
#endif
        for (int k = 0; k < 3; k++) {
#if defined(MK_MUTATE) && MK_MUTATE == 426
            const double x = V[q / 3 * 3 + k];
#else
            const double x = V[3 * (size_t)T[q] + k];
#endif
            lo[k] = fmin(lo[k], x), hi[k] = fmax(hi[k], x);
        }
    /* h from the budget, before the shells (which need it) */
    double A = 0, vol = 0;
    mk_refit_src_ S;
    if (mk_refit_src_init_(&S, V, nv, T, nt, lo, o->sign)) return -1;
    for (int c = 0; c < S.nc; c++) {
        A += S.ca[c];
        if (!S.co[c]) vol += fabs(S.cv[c]);
    }
    double h = o->h;
    if (!(h > 0) && !(isfinite(A) && isfinite(vol))) { mk_refit_src_free_(&S); return -1; } /* sizes past double's range */
    if (!(h > 0)) {
        const double tau = A > 0 ? 2 * vol / A : 0;
        for (h = o->hmin;; h *= 1.07)
#if defined(MK_MUTATE) && MK_MUTATE == 422
            if (1 || tau) break; /* 422: hmin whatever the budget */
#else
            if (0.5 * A * (tau + 8 * h) / (h * h * h) <= o->nact) break;
#endif
    }
    for (int c = 0; c < S.nc; c++) {
        const double thick = S.co[c] || S.ca[c] <= 0 ? 0 : 2 * fabs(S.cv[c]) / S.ca[c];
#if defined(MK_MUTATE) && MK_MUTATE == 413
        S.s[c] = S.co[c] ? fmax(o->shell, o->shell_h * h) : 0;
        (void)thick;
#else
        S.s[c] = S.co[c] ? fmax(o->shell, o->shell_h * h) : thick < o->thin ? h : 0;
#endif
    }
    R.pieces_in = S.nc, R.open_in = S.nc - S.ncl, R.vol_in = vol, R.h = h;
    for (int j = 0; j < nroom; j++) {
        int m = 1;
        for (int k = 0; k < 3; k++) m &= room[j].lo[k] - o->reach <= hi[k] && lo[k] <= room[j].hi[k] + o->reach;
        R.obstacles += m;
    }
    double gap = o->gap + o->gap_h * h, grow = 0;
    const double target = o->target;
    const int ntry = target > 0 ? o->tries_target : o->tries;
    mk_gmesh best = {0};
    int bok = 0, ret = 0;
    double bv = 0, g_fit = 0, g_over = INFINITY; /* the largest grow at or under the target, the least over it */
    mk_refit_report B = R;
    for (int tryn = 0; tryn < ntry; tryn++) {
        mk_gmesh m;
        long ntl, nnd;
        double drop;
        if (mk_refit_build_(V, T, nt, &S, room, nroom, o, lo, hi, h, gap, grow, &m, &ntl, &nnd, &drop)) { ret = -1; break; }
        R.tries = tryn + 1;
        if (!m.nt) {
            mk_gmesh_free(&m);
            if (!best.nt) B = R, B.gap = gap, B.grow = grow, B.tiles = ntl, B.nodes = nnd, B.dropped = drop, B.status = MK_REFIT_EMPTY;
            break;
        }
        mk_refit_report C = R;
        C.gap = gap, C.grow = grow, C.tiles = ntl, C.nodes = nnd, C.dropped = drop, C.least_gap = o->reach;
        mk_refit_obj ob;
        if (mk_refit_obj_init(&ob, m.V, m.nv, m.T, m.nt, 1)) { mk_gmesh_free(&m); ret = -1; break; }
        C.open_out = mk_refit_open_edges(m.T, m.nt), C.cross_out = mk_bvh_crossings(&ob.B, NULL);
        if (C.open_out < 0 || C.cross_out < 0) { mk_refit_obj_free(&ob), mk_gmesh_free(&m); ret = -1; break; }
        double v6 = 0, area = 0;
        for (int t = 0; t < m.nt; t++) {
            v6 += mk_refit_vol6_(m.V, m.T, t, lo);
            const double *a = m.V + 3 * (size_t)m.T[3 * (size_t)t], *b = m.V + 3 * (size_t)m.T[3 * (size_t)t + 1], *c = m.V + 3 * (size_t)m.T[3 * (size_t)t + 2];
            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            const double x = u[1] * w[2] - u[2] * w[1], y = u[2] * w[0] - u[0] * w[2], z = u[0] * w[1] - u[1] * w[0];
            area += 0.5 * sqrt(x * x + y * y + z * z);
        }
        C.vol_out = v6 / 6;
        { /* pieces of the result */
            int *par = malloc(sizeof *par * (size_t)m.nv), np = 0;
            if (!par) { mk_refit_obj_free(&ob), mk_gmesh_free(&m); ret = -1; break; }
            for (int v = 0; v < m.nv; v++) par[v] = v;
            for (int t = 0; t < m.nt; t++) {
                const int a = mk_refit_uf_(par, m.T[3 * (size_t)t]);
                for (int k = 1; k < 3; k++) {
                    const int b = mk_refit_uf_(par, m.T[3 * (size_t)t + k]);
                    if (a != b) par[b] = a;
                }
            }
            for (int v = 0; v < m.nv; v++) np += par[v] == v;
            free(par);
            C.pieces_out = np;
        }
        for (int j = 0; j < nroom; j++) {
            mk_refit_pair pr;
            const int r = mk_refit_check(&ob, room + j, o, &pr);
            if (r < 0) { ret = -1; break; }
            if (!pr.tested) continue;
            if (pr.gap < C.least_gap) {
                C.least_gap = pr.gap;
                if (!C.contacts) C.worst = j;
            }
            if (r) C.contacts++, C.worst = j;
        }
        if (ret < 0) { mk_refit_obj_free(&ob), mk_gmesh_free(&m); break; }
        { /* source pieces left with nothing near */
            unsigned char *near = calloc((size_t)S.nc + 1, 1);
            if (!near) { mk_refit_obj_free(&ob), mk_gmesh_free(&m); ret = -1; break; }
            const double r2 = (3 * (h + gap) + grow) * (3 * (h + gap) + grow); /* a grown part sits farther out */
            for (int v = 0; v < nv; v++) {
                const int c = S.vlab[v];
                if (c < 0 || near[c]) continue;
                int t;
                double q[3];
                (void)mk_bvh_nearest(&ob.B, V + 3 * (size_t)v, r2, &t, q);
                if (t >= 0) near[c] = 1;
            }
            for (int c = 0; c < S.nc; c++) C.lost += !near[c];
            free(near);
        }
        mk_refit_obj_free(&ob);
        C.status = C.open_out || C.cross_out || C.vol_out <= 0 ? MK_REFIT_BAD : C.contacts ? MK_REFIT_CONTACT : MK_REFIT_OK;
        const int okk = C.status == MK_REFIT_OK && (target <= 0 || C.vol_out <= target * 1.0001);
        if (!best.nt || (okk && (!bok || C.vol_out > bv)) || (!bok && !okk)) {
            mk_gmesh_free(&best);
            best = m, bok = okk, bv = C.vol_out, B = C;
        } else
            mk_gmesh_free(&m);
        if (C.contacts) {
#if !(defined(MK_MUTATE) && MK_MUTATE == 418) /* 418: no wider gap after a contact */
            gap += o->gap_step;
#endif
            continue;
        }
        if (target <= 0) break;
        const double dv = target - C.vol_out;
        if (okk) g_fit = fmax(g_fit, grow);
        else if (C.status == MK_REFIT_OK) g_over = fmin(g_over, grow);
#if defined(MK_MUTATE) && MK_MUTATE == 425
        if (fabs(dv) < o->grow_tol * target || (dv > 0 && grow >= o->grow_max)) break; /* 425: stops on a build over the target */
#else
        if ((okk && dv < o->grow_tol * target) || (dv > 0 && grow >= o->grow_max)) break; /* met from under, or no more room */
#endif
        const double g0 = grow;
        grow += area > 0 ? dv / area : 0;
        if (grow >= g_over || grow <= g_fit) grow = isfinite(g_over) ? 0.5 * (g_fit + g_over) : grow; /* over once: halve back */
        if (grow < 0) grow = 0;
        if (grow > o->grow_max) grow = o->grow_max;
        if (grow == g0) break;
    }
    B.tries = R.tries;
    mk_refit_src_free_(&S);
    if (ret < 0) { mk_gmesh_free(&best); return -1; }
    *out = best;
    if (rep) *rep = B;
    return B.status;
}
#endif
