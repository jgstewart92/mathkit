/* mk_refit and its checker on shapes with known answers, in millimetres: a ball refit against a ball it overlaps; a
   ball beside one 0.03 apart, one inside, crossing, and an open sheet through it (the checker's counts, inside only
   against a closed mesh); edge uses; a thin open sheet and a thin closed slab given their shells; two balls as one
   source (both pieces kept) and two overlapping balls (UNION solid, PARITY hollow at the overlap); a hollow ball's
   cavity (dropped, or kept on request); an open obstacle through a ball (clear on both sides, nothing taken as inside
   it); a first build in contact and the wider gap after it; a target volume, never passed; the budget's h against its
   formula; refusals; with OpenMP, the same bytes at 1, 2 and 12 threads. Writes OBJ files for meshcheck. */
#include "../include/mk_refit.h"
#include "t.h"
#include "tmesh.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* every result's own checks: closed, no crossing pair, outward (positive volume) */
static int clean(const mk_gmesh *m)
{
    mk_bvh B;
    if (!m->nt || mk_bvh_build(&B, m->V, m->T, m->nt)) return 0;
    const long x = mk_bvh_crossings(&B, NULL);
    mk_bvh_free(&B);
    return mk_refit_open_edges(m->T, m->nt) == 0 && x == 0 && tm_volume(m->V, m->T, m->nt) > 0;
}
/* the least distance from m's vertices to o's surface, and o's to m's, by brute force over faces */
static double gapbf(const double *V, int nv, const double *W, const int *T, int nt)
{
    double g = INFINITY;
    for (int i = 0; i < nv; i++)
        for (int t = 0; t < nt; t++) {
            double w[3], q[3];
            const double d = mk_closest_pt_tri(V + 3 * i, W + 3 * T[3 * t], W + 3 * T[3 * t + 1], W + 3 * T[3 * t + 2], w, q);
            if (d < g) g = d;
        }
    return sqrt(g);
}
[[maybe_unused]] static uint64_t hashm(const mk_gmesh *m)
{
    uint64_t h = 1469598103934665603u;
    const unsigned char *p = (const unsigned char *)m->V;
    for (size_t i = 0; i < 24 * (size_t)m->nv; i++) h = (h ^ p[i]) * 1099511628211u;
    p = (const unsigned char *)m->T;
    for (size_t i = 0; i < 12 * (size_t)m->nt; i++) h = (h ^ p[i]) * 1099511628211u;
    return h;
}
/* an open square sheet in the plane z = z0, side L about (cx, cy), n x n quads */
static void sheet(tm_mesh *m, double cx, double cy, double z0, double L, int n)
{
    const int base = m->nv;
    for (int j = 0; j <= n; j++)
        for (int i = 0; i <= n; i++) tm_vert(m, cx - L / 2 + L * i / n, cy - L / 2 + L * j / n, z0);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            const int a = base + j * (n + 1) + i, b = a + 1, c = a + n + 2, d = a + n + 1;
            tm_tri(m, a, b, c), tm_tri(m, a, c, d);
        }
}

int main(int argc, char **argv)
{
    tm_dir = argc > 1 ? argv[1] : NULL;
    mk_refit_opts o;
    mk_refit_defaults(&o, 1.0);
    o.h = 0.1, o.thin = 0; /* no shell for thin closed pieces but where it is tested: volumes then compare with the shapes' own */
    const double c0[3] = {0, 0, 0}, c1[3] = {1.5, 0, 0}, R = 1;
    { /* the defaults scale with mm */
        mk_refit_opts m3;
        mk_refit_defaults(&m3, 1e-3);
        CHECK(fabs(m3.hmin / 0.08e-3 - 1) < 1e-15 && fabs(m3.gap / 0.07e-3 - 1) < 1e-15 && fabs(m3.vmin / 5e-10 - 1) < 1e-15 && fabs(m3.reach / 2e-3 - 1) < 1e-15 && fabs(m3.thin / 1e-3 - 1) < 1e-15 && o.vmin == 0.5 && o.sign == MK_REFIT_UNION && o.tries == 3 && o.tries_target == 7,
              "defaults in metres: hmin %g, gap %g, vmin %g m^3 (0.5 mm^3), reach %g; in mm vmin %g; UNION, 3 tries (7 with a target)", m3.hmin, m3.gap, m3.vmin, m3.reach, o.vmin);
    }
    tm_mesh A = {0}, Bm = {0};
    tm_sphere(&A, R, c0, 4), tm_sphere(&Bm, R, c1, 4);
    mk_refit_obj ob, oa;
    CHECK(mk_refit_obj_init(&ob, Bm.V, Bm.nv, Bm.T, Bm.nt, 1) == 0 && mk_refit_obj_init(&oa, A.V, A.nv, A.T, A.nt, 1) == 0, "obstacle objects made");
    { /* edge uses */
        tm_mesh O = {0}, P = {0};
        tm_sphere(&O, R, c0, 2);
        const long e0 = mk_refit_open_edges(O.T, O.nt);
        O.nt--;
        const long e1 = mk_refit_open_edges(O.T, O.nt);
        tm_tri(&O, O.T[0], O.T[1], O.T[2]), tm_tri(&O, O.T[3 * (O.nt - 2)], O.T[3 * (O.nt - 2) + 1], O.T[3 * (O.nt - 2) + 2]);
        const long e2 = mk_refit_open_edges(O.T, O.nt);
        const double ca[3] = {0, 0, 0}, cb[3] = {2, 0, 0};
        tm_sphere(&P, R, ca, 2), tm_sphere(&P, R, cb, 2); /* two balls touching at (1, 0, 0), each with its own vertex there */
        const long e3 = mk_refit_open_edges(P.T, P.nt);
        CHECK(e0 == 0 && e1 == 3 && e2 >= 9 && e3 == 0, "edge uses on no edge used other than twice: closed ball %ld (0), a face less %ld (3), faces doubled %ld (an edge used three times counts 3 times), two balls touching at a point as given %ld (0)", e0, e1, e2, e3);
        tm_free(&O), tm_free(&P);
    }
    { /* the checker on pairs with known answers */
        tm_mesh G = {0}, In = {0}, X = {0}, Sh = {0};
        const double cg[3] = {2.03, 0, 0}, ci[3] = {0.2, 0, 0}, cx[3] = {1.0, 0, 0};
        tm_sphere(&G, R, cg, 4), tm_sphere(&In, 0.3, ci, 3), tm_sphere(&X, R, cx, 4);
        sheet(&Sh, 0, 0, 0.01, 3, 12);
        mk_refit_obj og, oi, ox, os, osc;
        mk_refit_obj_init(&og, G.V, G.nv, G.T, G.nt, 1), mk_refit_obj_init(&oi, In.V, In.nv, In.T, In.nt, 1), mk_refit_obj_init(&ox, X.V, X.nv, X.T, X.nt, 1);
        mk_refit_obj_init(&os, Sh.V, Sh.nv, Sh.T, Sh.nt, 0), mk_refit_obj_init(&osc, Sh.V, Sh.nv, Sh.T, Sh.nt, 1);
        mk_refit_pair p0, p1, p2, p3, p4, p5;
        const int r0 = mk_refit_check(&oa, &og, &o, &p0), r1 = mk_refit_check(&oa, &oi, &o, &p1), r2 = mk_refit_check(&oa, &ox, &o, &p2);
        const int r3 = mk_refit_check(&oa, &os, &o, &p3), r4 = mk_refit_check(&os, &oa, &o, &p4), r5 = mk_refit_check(&oa, &osc, &o, &p5);
        const double gb = fmin(gapbf(A.V, A.nv, G.V, G.T, G.nt), gapbf(G.V, G.nv, A.V, A.T, A.nt));
        CHECK(r0 == 1 && p0.tested && p0.cross == 0 && p0.a_in_b == 0 && p0.b_in_a == 0 && fabs(p0.gap - gb) < 1e-12 && p0.gap < 0.05,
              "two balls 0.03 apart: contact by the gap alone, %.6f (brute force %.6f, under gap_ok 0.05)", p0.gap, gb);
        CHECK(r1 == 1 && p1.b_in_a == In.nv && p1.a_in_b == 0 && p1.cross == 0, "a small ball inside: all %ld of its %d vertices inside, none of the big one's in it, no crossing", p1.b_in_a, In.nv);
        CHECK(r2 == 1 && p2.cross > 0 && p2.a_in_b > 0 && p2.b_in_a > 0, "two balls 1.0 apart: %ld crossing pairs, %ld and %ld vertices inside", p2.cross, p2.a_in_b, p2.b_in_a);
        CHECK(r3 == 1 && r4 == 1 && p3.cross > 0 && p4.cross == p3.cross && p3.a_in_b == 0 && p4.b_in_a == 0 && p4.a_in_b > 0,
              "an open sheet through the ball: %ld crossing pairs either way round; nothing inside the open sheet (%ld, %ld); %ld of its vertices inside the ball", p3.cross, p3.a_in_b, p4.b_in_a, p4.a_in_b);
        CHECK(r5 == 1 && p5.a_in_b > 0, "the same sheet flagged closed: %ld of the ball's vertices read inside it (why the flag matters)", p5.a_in_b);
        tm_mesh Far = {0};
        const double cf[3] = {5, 0, 0};
        tm_sphere(&Far, R, cf, 2);
        mk_refit_obj of;
        mk_refit_obj_init(&of, Far.V, Far.nv, Far.T, Far.nt, 1);
        mk_refit_pair pf;
        const int rf = mk_refit_check(&oa, &of, &o, &pf);
        CHECK(rf == 0 && !pf.tested && pf.gap == o.reach, "a ball 3 mm off: past reach, not tested, gap = reach %.1f", pf.gap);
        { /* a needle through a plate: only the needle's edges pierce the plate's two faces */
            tm_mesh Pl = {0}, Nd = {0};
            sheet(&Pl, 0, 0, 0, 2, 1);
            const double nlo[3] = {0.25, -0.25, -1}, ns[3] = {0.1, 0.1, 2}; /* off the plate's diagonal */
            tm_box(&Nd, nlo, ns, 1);
            mk_refit_obj opl, ond;
            mk_refit_obj_init(&opl, Pl.V, Pl.nv, Pl.T, Pl.nt, 0), mk_refit_obj_init(&ond, Nd.V, Nd.nv, Nd.T, Nd.nt, 1);
            mk_refit_pair q1, q2;
            mk_refit_check(&opl, &ond, &o, &q1), mk_refit_check(&ond, &opl, &o, &q2);
            CHECK(q1.cross > 0 && q1.cross == q2.cross, "a needle through a plate of two faces: %ld and %ld crossing pairs, either way round", q1.cross, q2.cross);
            mk_refit_obj_free(&opl), mk_refit_obj_free(&ond), tm_free(&Pl), tm_free(&Nd);
        }
        mk_refit_obj_free(&og), mk_refit_obj_free(&oi), mk_refit_obj_free(&ox), mk_refit_obj_free(&os), mk_refit_obj_free(&osc), mk_refit_obj_free(&of);
        tm_free(&G), tm_free(&In), tm_free(&X), tm_free(&Sh), tm_free(&Far);
    }
    { /* a ball refit against a ball it overlaps (centres 1.5 apart): closed, clear of it, the cut off */
        mk_gmesh m;
        mk_refit_report r;
        const int s = mk_refit(A.V, A.nv, A.T, A.nt, &ob, 1, &o, &m, &r);
        mk_refit_obj om;
        mk_refit_pair p;
        mk_refit_obj_init(&om, m.V, m.nv, m.T, m.nt, 1);
        const int c = mk_refit_check(&om, &ob, &o, &p);
        const double va = tm_volume(A.V, A.T, A.nt), v = tm_volume(m.V, m.T, m.nt), g = gapbf(m.V, m.nv, Bm.V, Bm.T, Bm.nt);
        /* the ball less its lens with the other grown by the gap 0.13 (radius 1.13, 1.5 apart): 0.595 */
        CHECK(s == MK_REFIT_OK && clean(&m) && c == 0 && r.contacts == 0 && r.pieces_out == 1 && r.lost == 0 && g >= o.gap_ok && v < va - 0.45 && v > va - 0.75 && fabs(r.vol_out - v) < 1e-9 * va,
              "a ball against a ball it overlaps: OK, closed and uncrossed, 1 piece, clear (least gap %.4f, report %.4f), volume %.4f of %.4f", g, r.least_gap, v, va);
        CHECK(r.h == 0.1 && fabs(r.gap - (0.07 + 0.06)) < 1e-12 && r.tries == 1 && r.obstacles == 1 && r.pieces_in == 1 && r.open_in == 0 && fabs(r.vol_in - va) < 1e-9 * va,
              "its report: h %.3f, gap %.3f (0.07 + 0.6 h), %d try, %d obstacle, %d piece in", r.h, r.gap, r.tries, r.obstacles, r.pieces_in);
        tm_write("refit_ball", m.V, m.nv, m.T, m.nt, NULL);
        mk_refit_obj_free(&om), mk_gmesh_free(&m);
    }
    { /* a thin open sheet: a closed slab of its shell's thickness; a thin closed slab kept two cells thick */
        tm_mesh Sh = {0}, Sl = {0};
        mk_refit_opts o1 = o;
        o1.thin = 1;
        sheet(&Sh, 0, 0, 0, 2, 8);
        const double slo[3] = {-1, -1, -0.04}, ss[3] = {2, 2, 0.08};
        tm_box(&Sl, slo, ss, 4);
        mk_gmesh m, n;
        mk_refit_report r, q;
        const int s = mk_refit(Sh.V, Sh.nv, Sh.T, Sh.nt, NULL, 0, &o1, &m, &r), t = mk_refit(Sl.V, Sl.nv, Sl.T, Sl.nt, NULL, 0, &o1, &n, &q);
        const double v = tm_volume(m.V, m.T, m.nt), w = tm_volume(n.V, n.T, n.nt);
        /* the sheet's shell, the sheet swept by a ball of radius max(0.15, 0.12): 2 x 2 x 0.3 + its rim 8 x pi 0.15^2 / 2 = 1.48;
           the slab swept by a ball of radius h: 0.32 + 0.2 (4 + 0.32) + pi 0.01 (4.08) + 4/3 pi 0.001 = 1.32 */
        CHECK(s == MK_REFIT_OK && clean(&m) && r.open_in == 1 && r.pieces_out == 1 && v > 1.3 && v < 1.6, "an open sheet 2 x 2: OK, closed, 1 piece, volume %.4f (its shell 1.48)", v);
        CHECK(t == MK_REFIT_OK && clean(&n) && q.pieces_out == 1 && q.lost == 0 && w > 1.15 && w < 1.45, "a closed slab 0.08 thick, under h: OK, closed, 1 piece, volume %.4f (with its shell of h 1.32; without, it falls apart)", w);
        tm_write("refit_sheet", m.V, m.nv, m.T, m.nt, NULL), tm_write("refit_slab", n.V, n.nv, n.T, n.nt, NULL);
        mk_gmesh_free(&m), mk_gmesh_free(&n), tm_free(&Sh), tm_free(&Sl);
    }
    { /* two balls apart as one source: both kept; two overlapping: UNION solid, PARITY hollow where they overlap */
        tm_mesh Two = {0}, Ov = {0};
        const double ca[3] = {0, 0, 0}, cb[3] = {3, 0, 0}, cc[3] = {1.0, 0, 0};
        tm_sphere(&Two, R, ca, 4), tm_sphere(&Two, 0.6, cb, 4);
        tm_sphere(&Ov, R, ca, 4), tm_sphere(&Ov, R, cc, 4);
        mk_gmesh m, u, p;
        mk_refit_report r, ru, rp;
        mk_refit_opts ou = o, op = o; /* voids kept, so a hollow shows */
        ou.keep_voids = 1, op.keep_voids = 1, op.sign = MK_REFIT_PARITY;
        const int s = mk_refit(Two.V, Two.nv, Two.T, Two.nt, NULL, 0, &o, &m, &r);
        const int su = mk_refit(Ov.V, Ov.nv, Ov.T, Ov.nt, NULL, 0, &ou, &u, &ru), sp = mk_refit(Ov.V, Ov.nv, Ov.T, Ov.nt, NULL, 0, &op, &p, &rp);
        const double lens = T_PI * (4 * R + 1.0) * (2 * R - 1.0) * (2 * R - 1.0) / 12, uni = 2 * 4.0 / 3 * T_PI - lens;
        const double vu = tm_volume(u.V, u.T, u.nt), vp = tm_volume(p.V, p.T, p.nt);
        CHECK(s == MK_REFIT_OK && clean(&m) && r.pieces_in == 2 && r.pieces_out == 2 && r.lost == 0, "two balls 3 apart (radii 1 and 0.6) as one source: %d pieces in, %d out, %d lost", r.pieces_in, r.pieces_out, r.lost);
        CHECK(su == MK_REFIT_OK && clean(&u) && ru.pieces_out == 1 && fabs(vu / uni - 1) < 0.03, "two balls overlapping by a lens of %.4f, UNION: one solid piece (%d, status %d, clean %d), volume %.4f against their union %.4f", lens, ru.pieces_out, su, clean(&u), vu, uni);
        CHECK(sp == MK_REFIT_OK && clean(&p) && vp < vu - 0.5 * lens, "the same, PARITY: the overlap read outside, volume %.4f (UNION %.4f, lens %.4f)", vp, vu, lens);
        tm_write("refit_union", u.V, u.nv, u.T, u.nt, NULL);
        mk_gmesh_free(&m), mk_gmesh_free(&u), mk_gmesh_free(&p), tm_free(&Two), tm_free(&Ov);
    }
    { /* a hollow ball (outer shell outward, inner inward): PARITY keeps the cavity in the field; the void is dropped, or
         kept with keep_voids. A ball much larger than the band: its inside far from the surface stays solid. */
        tm_mesh H = {0};
        tm_sphere(&H, 1.5, c0, 4);
        const int b = H.nv;
        tm_sphere(&H, 0.8, c0, 4);
        for (int t = 0; t < H.nt; t++)
            if (H.T[3 * t] >= b) { const int x = H.T[3 * t + 1]; H.T[3 * t + 1] = H.T[3 * t + 2], H.T[3 * t + 2] = x; }
        mk_refit_opts op = o;
        op.sign = MK_REFIT_PARITY;
        mk_gmesh m, k;
        mk_refit_report r, rk;
        const int s = mk_refit(H.V, H.nv, H.T, H.nt, NULL, 0, &op, &m, &r);
        op.keep_voids = 1;
        const int sk = mk_refit(H.V, H.nv, H.T, H.nt, NULL, 0, &op, &k, &rk);
        const double vo = 4.0 / 3 * T_PI * 1.5 * 1.5 * 1.5, vi = 4.0 / 3 * T_PI * 0.8 * 0.8 * 0.8, v = tm_volume(m.V, m.T, m.nt), w = tm_volume(k.V, k.T, k.nt);
        CHECK(s == MK_REFIT_OK && clean(&m) && r.pieces_out == 1 && fabs(v / vo - 1) < 0.03, "a hollow ball, voids dropped: 1 piece, volume %.4f (the outer ball %.4f)", v, vo);
        CHECK(sk == MK_REFIT_OK && rk.pieces_out == 2 && fabs(w / (vo - vi) - 1) < 0.04, "kept: 2 pieces, volume %.4f (the shell %.4f)", w, vo - vi);
        mk_gmesh_free(&m), mk_gmesh_free(&k), tm_free(&H);
        tm_mesh Big = {0};
        tm_sphere(&Big, 4, c0, 4);
        mk_gmesh g;
        mk_refit_report rg;
        mk_refit_opts o2 = o;
        o2.h = 0.1, o2.keep_voids = 1; /* a hollow would show; 2 x 2 x 2 tiles, welded */
        const int sg = mk_refit(Big.V, Big.nv, Big.T, Big.nt, NULL, 0, &o2, &g, &rg);
        const double vb = tm_volume(Big.V, Big.T, Big.nt), vg = tm_volume(g.V, g.T, g.nt);
        CHECK(sg == MK_REFIT_OK && clean(&g) && rg.pieces_out == 1 && rg.tiles == 8 && fabs(vg / vb - 1) < 0.02,
              "a ball of radius 4 at h 0.1, far nodes inside it: 1 solid piece over %ld tiles, closed where they meet, volume %.3f of %.3f", rg.tiles, vg, vb);
        mk_gmesh_free(&g), tm_free(&Big);
    }
    { /* an open obstacle through the ball: the ball cut in two, clear of the sheet on both sides; nothing taken as inside it */
        tm_mesh Sh = {0};
        sheet(&Sh, 0, 0, 0.003, 3, 24);
        mk_refit_obj os;
        mk_refit_obj_init(&os, Sh.V, Sh.nv, Sh.T, Sh.nt, 0);
        mk_gmesh m;
        mk_refit_report r;
        const int s = mk_refit(A.V, A.nv, A.T, A.nt, &os, 1, &o, &m, &r);
        const double va = tm_volume(A.V, A.T, A.nt), v = tm_volume(m.V, m.T, m.nt), g = gapbf(m.V, m.nv, Sh.V, Sh.T, Sh.nt);
        CHECK(s == MK_REFIT_OK && clean(&m) && r.pieces_out == 2 && g >= o.gap_ok && v > va - 1.1 && v < va - 0.5,
              "a ball through an open sheet: OK, 2 pieces, clear of it (least gap %.4f), volume %.4f of %.4f (a slot 0.26 wide through the middle: about 0.82)", g, v, va);
        tm_write("refit_slot", m.V, m.nv, m.T, m.nt, NULL);
        mk_refit_obj_free(&os), mk_gmesh_free(&m), tm_free(&Sh);
    }
    { /* a first build in contact: no clearance built, so the checker's 0.05 is missed; the next try is wider */
        mk_refit_opts oc = o;
        oc.gap = 0, oc.gap_h = 0, oc.gap_step = 0.06, oc.dec_h = 0;
        mk_gmesh m;
        mk_refit_report r;
        const int s = mk_refit(A.V, A.nv, A.T, A.nt, &ob, 1, &oc, &m, &r);
        CHECK(s == MK_REFIT_OK && r.tries >= 2 && r.contacts == 0 && r.least_gap >= oc.gap_ok && fabs(r.gap - 0.06 * (r.tries - 1)) < 1e-12,
              "no clearance built: %d tries, the kept one at gap %.2f clear (least gap %.4f)", r.tries, r.gap, r.least_gap);
        mk_gmesh_free(&m);
        oc.tries = 1;
        const int s1 = mk_refit(A.V, A.nv, A.T, A.nt, &ob, 1, &oc, &m, &r);
        CHECK(s1 == MK_REFIT_CONTACT && r.contacts == 1 && r.worst == 0, "the same with one try: CONTACT, %d contact with obstacle %d (least gap %.4f)", r.contacts, r.worst, r.least_gap);
        mk_gmesh_free(&m);
    }
    { /* a target volume: the cut ball grows into free room toward the whole ball's volume, never past it */
        mk_refit_opts ot = o;
        const double va = tm_volume(A.V, A.T, A.nt);
        mk_gmesh m0, m;
        mk_refit_report r0, r;
        mk_refit(A.V, A.nv, A.T, A.nt, &ob, 1, &o, &m0, &r0);
        ot.target = va;
        const int s = mk_refit(A.V, A.nv, A.T, A.nt, &ob, 1, &ot, &m, &r);
        const double v0 = tm_volume(m0.V, m0.T, m0.nt), v = tm_volume(m.V, m.T, m.nt);
        CHECK(s == MK_REFIT_OK && clean(&m) && r.contacts == 0 && v > v0 + 0.1 && v <= va * 1.0001 && r.grow > 0 && r.grow <= ot.grow_max && r.tries >= 2,
              "a target of the whole ball %.4f: grown %.3f in %d tries, volume %.4f (without: %.4f), never past it, clear", va, r.grow, r.tries, v, v0);
        ot.target = 0.5 * va;
        mk_gmesh_free(&m);
        const int s2 = mk_refit(A.V, A.nv, A.T, A.nt, &ob, 1, &ot, &m, &r);
        CHECK(s2 == MK_REFIT_OK && r.grow == 0 && r.vol_out == r0.vol_out && r.tries == 1, "a target under what is built: no growth, the result as without (%.4f), %d try", r.vol_out, r.tries);
        mk_gmesh_free(&m0), mk_gmesh_free(&m);
    }
    { /* the checker's gap is face to face: an edge 0.02 over a plate with every vertex 1 or more away */
        tm_mesh P = {0}, E = {0};
        sheet(&P, 0, 0, 0, 2, 1);
        tm_vert(&E, -3, 0, 0.02), tm_vert(&E, 3, 0, 0.02), tm_vert(&E, 0, 0, 1), tm_tri(&E, 0, 1, 2);
        mk_refit_obj op, oe;
        mk_refit_pair p;
        mk_refit_obj_init(&op, P.V, P.nv, P.T, P.nt, 0), mk_refit_obj_init(&oe, E.V, E.nv, E.T, E.nt, 0);
        const int c = mk_refit_check(&op, &oe, &o, &p);
        CHECK(c == 1 && p.cross == 0 && fabs(p.gap - 0.02) < 1e-12, "a face's edge 0.02 over a plate, every vertex 1 or more off: contact, gap %.6f (0.02), %ld crossing pairs", p.gap, p.cross);
        mk_refit_obj_free(&op), mk_refit_obj_free(&oe), tm_free(&P), tm_free(&E);
    }
    { /* a plate 0.4 thick cutting into a ball of radius 3, gap = gap_ok = 0.05 at h 0.3: clear of it face to face, checked by
         sampling the result's faces against the plate's surface */
        tm_mesh Bl = {0}, Pl = {0};
        tm_sphere(&Bl, 3, c0, 5);
        const double x0 = 1.5, x1 = 8, y0 = -6, y1 = 6, z0 = -0.2, z1 = 0.2;
        const double cx[8][3] = {{x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0}, {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1}};
        static const int f[12][3] = {{0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5}, {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
        for (int i = 0; i < 8; i++) tm_vert(&Pl, cx[i][0], cx[i][1], cx[i][2]);
        for (int i = 0; i < 12; i++) tm_tri(&Pl, f[i][0], f[i][1], f[i][2]);
        mk_refit_obj opl;
        mk_refit_obj_init(&opl, Pl.V, Pl.nv, Pl.T, Pl.nt, 1);
        mk_refit_opts og = o;
        og.h = 0.3, og.gap = og.gap_ok = 0.05, og.gap_h = 0, og.tries = 8; /* thinning moves the surface up to 0.2 h: more tries widen the gap */
        mk_gmesh m;
        mk_refit_report r;
        const int s = mk_refit(Bl.V, Bl.nv, Bl.T, Bl.nt, &opl, 1, &og, &m, &r);
        double g = INFINITY;
        for (int t = 0; t < m.nt; t++)
            for (int i = 0; i <= 12; i++)
                for (int j = 0; i + j <= 12; j++) {
                    double p[3], q[3];
                    int tt;
                    for (int k = 0; k < 3; k++) p[k] = (m.V[3 * m.T[3 * t] + k] * (12 - i - j) + m.V[3 * m.T[3 * t + 1] + k] * i + m.V[3 * m.T[3 * t + 2] + k] * j) / 12;
                    g = fmin(g, sqrt(mk_bvh_nearest(&opl.B, p, 1, &tt, q)));
                }
        CHECK(s == MK_REFIT_OK && clean(&m) && g >= og.gap_ok && r.least_gap >= og.gap_ok, "a plate cutting into a ball, gap 0.05 at h 0.3, 8 tries: OK in %d tries, least gap %.4f reported, %.4f sampled on its faces", r.tries, r.least_gap, g);
        mk_gmesh_free(&m);
        { /* a growth target 1.3 times the ball's volume, the room far: met from under, never past it */
            mk_refit_obj ofar;
            const double cf[3] = {40, 0, 0};
            tm_mesh Fb = {0};
            tm_sphere(&Fb, 1, cf, 3);
            mk_refit_obj_init(&ofar, Fb.V, Fb.nv, Fb.T, Fb.nt, 1);
            mk_refit_opts ot = og;
            const double vb = tm_volume(Bl.V, Bl.T, Bl.nt);
            ot.target = 1.3 * vb;
            const int s2 = mk_refit(Bl.V, Bl.nv, Bl.T, Bl.nt, &ofar, 1, &ot, &m, &r);
            const double v = tm_volume(m.V, m.T, m.nt);
            CHECK(s2 == MK_REFIT_OK && clean(&m) && v <= ot.target * 1.0001 && v >= ot.target * (1 - ot.grow_tol), "a target 1.3 times the ball (%.3f): %.3f in %d tries, grown %.4f", ot.target, v, r.tries, r.grow);
            mk_gmesh_free(&m), mk_refit_obj_free(&ofar), tm_free(&Fb);
        }
        { /* a vertex no face uses, far off, changes nothing */
            mk_gmesh m0, m1;
            mk_refit_report r0, r1;
            const int s0 = mk_refit(Bl.V, Bl.nv, Bl.T, Bl.nt, &opl, 1, &og, &m0, &r0);
            tm_vert(&Bl, 1e7, 1e7, 1e7);
            const int s1 = mk_refit(Bl.V, Bl.nv, Bl.T, Bl.nt, &opl, 1, &og, &m1, &r1);
            CHECK(s0 == MK_REFIT_OK && s1 == MK_REFIT_OK && hashm(&m0) == hashm(&m1) && r0.tiles == r1.tiles, "a vertex of no face at 1e7: the same bytes (%016llx), %ld tiles either way", (unsigned long long)hashm(&m1), r1.tiles);
            mk_gmesh_free(&m0), mk_gmesh_free(&m1);
        }
        mk_refit_obj_free(&opl), tm_free(&Bl), tm_free(&Pl);
    }
    { /* the budget: h from hmin by 7 % steps until 0.5 A (2 V / A + 8 h) / h^3 fits nact */
        mk_refit_opts ob2 = o;
        ob2.h = 0, ob2.hmin = 0.02, ob2.nact = 2e5;
        mk_gmesh m;
        mk_refit_report r;
        mk_refit(A.V, A.nv, A.T, A.nt, NULL, 0, &ob2, &m, &r);
        double area = 0;
        for (int t = 0; t < A.nt; t++) {
            const double *a = A.V + 3 * A.T[3 * t], *b = A.V + 3 * A.T[3 * t + 1], *c = A.V + 3 * A.T[3 * t + 2];
            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            area += 0.5 * sqrt((u[1] * w[2] - u[2] * w[1]) * (u[1] * w[2] - u[2] * w[1]) + (u[2] * w[0] - u[0] * w[2]) * (u[2] * w[0] - u[0] * w[2]) + (u[0] * w[1] - u[1] * w[0]) * (u[0] * w[1] - u[1] * w[0]));
        }
        const double tau = 2 * tm_volume(A.V, A.T, A.nt) / area;
        double h = 0.02;
        while (0.5 * area * (tau + 8 * h) / (h * h * h) > 2e5) h *= 1.07;
        CHECK(fabs(r.h - h) < 1e-12 * h && r.h > 0.02, "the budget's h: %.5f (by its formula %.5f, from hmin 0.02 at 2e5 nodes)", r.h, h);
        mk_gmesh_free(&m);
    }
    { /* refusals */
        mk_gmesh m;
        mk_refit_report r;
        mk_refit_opts bad = o;
        int n = 0, k = 0;
        double nanv[9] = {0, 0, 0, 1, 0, 0, 0, NAN, 0};
        const int T1[3] = {0, 1, 2}, Tb[3] = {0, 1, 3};
#define REF(call) (k++, n += (call) == -1 && m.nt == 0 && m.V == NULL)
        REF(mk_refit(A.V, A.nv, A.T, 0, NULL, 0, &o, &m, &r));
        REF(mk_refit(A.V, 3, Tb, 1, NULL, 0, &o, &m, &r));
        REF(mk_refit(nanv, 3, T1, 1, NULL, 0, &o, &m, &r));
        REF(mk_refit(A.V, A.nv, A.T, A.nt, NULL, -1, &o, &m, &r));
        bad.sign = 7;
        REF(mk_refit(A.V, A.nv, A.T, A.nt, NULL, 0, &bad, &m, &r));
        bad = o, bad.h = NAN;
        REF(mk_refit(A.V, A.nv, A.T, A.nt, NULL, 0, &bad, &m, &r));
        bad = o, bad.tries = 0;
        REF(mk_refit(A.V, A.nv, A.T, A.nt, NULL, 0, &bad, &m, &r));
        bad = o, bad.hmin = 0, bad.h = 0;
        REF(mk_refit(A.V, A.nv, A.T, A.nt, NULL, 0, &bad, &m, &r));
        mk_refit_obj x;
        k++, n += mk_refit_obj_init(&x, nanv, 3, T1, 1, 1) == -1;
        k++, n += mk_refit_obj_init(&x, A.V, 3, Tb, 1, 1) == -1;
        mk_refit_pair p;
        k++, n += mk_refit_check(&oa, NULL, &o, &p) == -1;
        k++, n += mk_refit_open_edges(NULL, 3) == -1;
        CHECK(n == k, "refusals: %d of %d (no face, an index out of range, NaN, a negative room, a bad sign, h NaN, no tries, no h at all, bad objects, no obstacle, no faces)", n, k);
    }
#ifdef _OPENMP
    { /* the same bytes at 1, 2 and 12 threads */
        tm_mesh Ov = {0};
        const double ca[3] = {0, 0, 0}, cc[3] = {1.2, 0.3, 0};
        tm_sphere(&Ov, R, ca, 4), tm_sphere(&Ov, 0.7, cc, 4);
        mk_refit_opts ot = o;
        ot.h = 0.03; /* several tiles */
        uint64_t hs[3];
        long tiles = 0;
        const int nth[3] = {1, 2, 12};
        for (int i = 0; i < 3; i++) {
            omp_set_num_threads(nth[i]);
            mk_gmesh m;
            mk_refit_report r;
            mk_refit(Ov.V, Ov.nv, Ov.T, Ov.nt, &ob, 1, &ot, &m, &r);
            hs[i] = hashm(&m), tiles = r.tiles;
            mk_gmesh_free(&m);
        }
        CHECK(hs[0] == hs[1] && hs[0] == hs[2] && tiles > 4, "the same bytes at 1, 2 and 12 threads (%016llx, %ld tiles)", (unsigned long long)hs[0], tiles);
        tm_free(&Ov);
    }
#endif
    mk_refit_obj_free(&ob), mk_refit_obj_free(&oa), tm_free(&A), tm_free(&Bm);
    return wrong;
}
