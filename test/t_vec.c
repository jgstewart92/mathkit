/* mk_vec.h: quaternion and transform identities on 1000 random rotations, integration at a constant rate against the
   exact turn, normalization across the exponent range, NaN propagation, mk_lerp at its ends, and the version string. */
#include "../include/mk_vec.h"
#include "../include/mk_version.h"
#include <string.h>
#include "t.h"

static double vd(mk_v3 a, mk_v3 b) { return mk_v3_dist(a, b); }
static mk_v3 rv(void) { return mk_v3_make(t_rand(), t_rand(), t_rand()); }

int main(void)
{
    mk_v3 r = mk_q_rotate(mk_q_axis_angle(mk_v3_make(0, 0, 1), T_PI / 2), mk_v3_make(1, 0, 0));
    CHECK(vd(r, mk_v3_make(0, 1, 0)) < 1e-15, "x turned 90 deg about z is y (error %.2g)", vd(r, mk_v3_make(0, 1, 0)));
    double worst_c = 0, worst_m = 0, worst_u = 0, worst_x = 0, worst_l = 0;
    for (int k = 0; k < 1000; k++) {
        mk_quat a = mk_q_from_rotvec(mk_v3_scale(rv(), 3)), b = mk_q_from_rotvec(mk_v3_scale(rv(), 3));
        mk_v3 v = rv();
        worst_c = fmax(worst_c, vd(mk_q_rotate(mk_q_mul(a, b), v), mk_q_rotate(a, mk_q_rotate(b, v))));
        double R[3][3];
        mk_q_to_mat(a, R);
        mk_v3 mv = {R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z, R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z,
                    R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z};
        mk_quat a2 = mk_q_from_mat(R);
        worst_m = fmax(worst_m, fmax(vd(mv, mk_q_rotate(a, v)), vd(mk_q_rotate(a2, v), mk_q_rotate(a, v))));
        worst_u = fmax(worst_u, vd(mk_q_unrotate(a, mk_q_rotate(a, v)), v));
        mk_xf x = {a, rv()};
        worst_x = fmax(worst_x, vd(mk_xf_apply(mk_xf_compose(mk_xf_inv(x), x), v), v));
        mk_quat h = mk_q_slerp(mk_q_identity(), a, 0.5);
        double full = 2 * acos(fmin(1, fabs(a.w))), half = 2 * acos(fmin(1, fabs(h.w)));
        worst_l = fmax(worst_l, fabs(half - full / 2));
    }
    CHECK(worst_c < 1e-14, "rotate(a b, v) = rotate(a, rotate(b, v)), 1000 random (worst %.2g)", worst_c);
    CHECK(worst_m < 1e-14, "matrix form and back agree with the rotation (worst %.2g)", worst_m);
    CHECK(worst_u < 1e-14, "unrotate undoes rotate (worst %.2g)", worst_u);
    CHECK(worst_x < 1e-14, "inverse transform composed with it is the identity (worst %.2g)", worst_x);
    CHECK(worst_l < 1e-12, "slerp halfway turns half the angle (worst %.2g)", worst_l);
    mk_v3 w = mk_v3_make(0.3, -0.7, 0.2);
    mk_quat q = mk_q_identity();
    for (int i = 0; i < 1000; i++) q = mk_q_integrate(q, w, 0.001);
    mk_quat e = mk_q_from_rotvec(w);
    double ie = vd(mk_q_rotate(q, mk_v3_make(1, 2, 3)), mk_q_rotate(e, mk_v3_make(1, 2, 3)));
    CHECK(ie < 1e-12, "1000 steps at constant omega equal one turn of omega t (error %.2g)", ie);
    CHECK(mk_v3_len(mk_v3_norm(mk_v3_make(0, 0, 0))) == 0, "the unit vector of zero is zero, not NaN");
    CHECK(strcmp(MK_VERSION, "1.6.4") == 0 && MK_VERSION_MAJOR == 1 && MK_VERSION_MINOR == 6 && MK_VERSION_PATCH == 4, "mathkit version %s", MK_VERSION);
    { /* a zero axis gives the identity */
        mk_quat z = mk_q_axis_angle(mk_v3_make(0, 0, 0), 1);
        CHECK(z.w == 1 && z.x == 0 && z.y == 0 && z.z == 0, "a zero axis gives the identity (got %g, %g, %g, %g)", z.w, z.x, z.y, z.z);
    }
    { /* normalization without overflow or underflow: (3, 4, 0) x 2^e, e from -1070 to 1020; a quaternion near 1e-170 */
        int ok = 1;
        for (int e = -1070; e <= 1020; e += 10) {
            mk_v3 u = mk_v3_norm((mk_v3){3 * ldexp(1, e), 4 * ldexp(1, e), 0});
            ok &= fabs(u.x - 0.6) < 1e-15 && fabs(u.y - 0.8) < 1e-15 && u.z == 0;
        }
        CHECK(ok, "v3_norm of (3, 4, 0) x 2^e is (0.6, 0.8, 0) for e from -1070 to 1020");
        mk_quat q = mk_q_norm((mk_quat){1e-170, 1e-170, 0, 0});
        CHECK(fabs(q.w - sqrt(0.5)) < 1e-15 && fabs(q.x - sqrt(0.5)) < 1e-15, "q_norm of a quaternion near 1e-170 keeps its direction (w %g, x %g)", q.w, q.x);
    }
    { /* a NaN quaternion, rate or matrix entry stays NaN, so a diverged attitude never reads as level; zero normalizes
         to the identity */
        const mk_quat q = mk_q_axis_angle((mk_v3){0, 0, 1}, 0.7), a = mk_q_norm((mk_quat){NAN, 0, 0, 0}), b = mk_q_integrate(q, (mk_v3){NAN, 0, 0}, 0.01);
        const double R[3][3] = {{NAN, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        const mk_quat c = mk_q_from_mat(R), z = mk_q_norm((mk_quat){0, 0, 0, 0});
        CHECK(isnan(a.w) && isnan(b.w) && isnan(c.w) && z.w == 1 && z.x == 0, "q_norm(NaN, 0, 0, 0) w %g, q_integrate with a NaN rate w %g, q_from_mat with a NaN entry w %g; zero still the identity (%g)",
              a.w, b.w, c.w, z.w);
    }
    { /* a zero axis with a NaN or infinite angle gives NaN, not the identity */
        const mk_quat a = mk_q_axis_angle((mk_v3){0, 0, 0}, NAN), b = mk_q_axis_angle((mk_v3){0, 0, 0}, INFINITY), c = mk_q_axis_angle((mk_v3){0, 0, 0}, 0.7);
        CHECK(isnan(a.w) && isnan(b.w) && c.w == 1 && c.x == 0, "q_axis_angle about a zero axis: angle NaN w %g, inf w %g; 0.7 still the identity (w %g)", a.w, b.w, c.w);
    }
    { /* mk_lerp gives exactly a at t 0 and b at t 1, even beside an infinite end ((1 - t) a + t b gives 0 inf = NaN) */
        const double a = mk_lerp(1, INFINITY, 0), b = mk_lerp(1, -INFINITY, 0), c = mk_lerp(-0.0, 1, 0), d = mk_lerp(1, INFINITY, 1);
        CHECK(a == 1 && b == 1 && c == 0 && signbit(c) && d == INFINITY, "lerp(1, inf, 0) %g, lerp(1, -inf, 0) %g, lerp(-0, 1, 0) %g (exactly a), lerp(1, inf, 1) %g", a, b,
              c, d);
    }
    return wrong;
}
