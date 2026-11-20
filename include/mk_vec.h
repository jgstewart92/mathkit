/* mathkit: 3-vectors, unit quaternions (w, x, y, z) and rigid transforms, double precision. Header only, C23. */
#ifndef MK_VEC_H
#define MK_VEC_H
#include <math.h>

typedef struct { double x, y, z; } mk_v3;
typedef struct { double w, x, y, z; } mk_quat;
typedef struct { mk_quat r; mk_v3 t; } mk_xf; /* p' = R p + t */

static inline mk_v3 mk_v3_make(double x, double y, double z) { return (mk_v3){x, y, z}; }
static inline mk_v3 mk_v3_add(mk_v3 a, mk_v3 b) { return (mk_v3){a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline mk_v3 mk_v3_sub(mk_v3 a, mk_v3 b) { return (mk_v3){a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline mk_v3 mk_v3_scale(mk_v3 a, double s) { return (mk_v3){a.x * s, a.y * s, a.z * s}; }
static inline mk_v3 mk_v3_mad(mk_v3 a, mk_v3 b, double s) { return (mk_v3){a.x + b.x * s, a.y + b.y * s, a.z + b.z * s}; }
static inline double mk_v3_dot(mk_v3 a, mk_v3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline mk_v3 mk_v3_cross(mk_v3 a, mk_v3 b)
{
    return (mk_v3){a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
static inline double mk_v3_len2(mk_v3 a) { return mk_v3_dot(a, a); }
static inline double mk_v3_len(mk_v3 a) { return sqrt(mk_v3_dot(a, a)); }
static inline double mk_v3_dist(mk_v3 a, mk_v3 b) { return mk_v3_len(mk_v3_sub(a, b)); }
/* unit vector; a zero vector comes back unchanged */
static inline mk_v3 mk_v3_norm(mk_v3 a)
{
    double m = fmax(fabs(a.x), fmax(fabs(a.y), fabs(a.z)));
    if (m > 0 && m < INFINITY) { /* by a power of two, exact: the same result wherever |a|^2 did not under- or overflow */
        int e = ilogb(m);
        a = (mk_v3){ldexp(a.x, -e), ldexp(a.y, -e), ldexp(a.z, -e)};
    }
    double l = mk_v3_len(a);
    return l > 0 ? mk_v3_scale(a, 1.0 / l) : a;
}
/* a + t (b - a), exactly a at t = 0 and b at t = 1 (1.6.2: a also beside an infinite b, which gave NaN, and -0 kept); where b - a
   overflows, (1 - t) a + t b (1.6.0, the review's S1 and I1). The kit's one linear interpolation: mk_interp's and mk_stats' use it. */
static inline double mk_lerp(double a, double b, double t)
{
    const double d = b - a;
    if (t == 0) return a;
    if (t == 1) return b;
    return isfinite(d) ? a + t * d : (1 - t) * a + t * b;
}
static inline mk_v3 mk_v3_lerp(mk_v3 a, mk_v3 b, double t) { return (mk_v3){mk_lerp(a.x, b.x, t), mk_lerp(a.y, b.y, t), mk_lerp(a.z, b.z, t)}; }

static inline mk_quat mk_q_identity(void) { return (mk_quat){1, 0, 0, 0}; }
static inline mk_quat mk_q_mul(mk_quat a, mk_quat b)
{
    return (mk_quat){a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
#if defined(MK_MUTATE) && MK_MUTATE == 2
                     a.w * b.x + a.x * b.w - a.y * b.z + a.z * b.y,
#else
                     a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
#endif
                     a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
                     a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}
static inline mk_quat mk_q_conj(mk_quat q) { return (mk_quat){q.w, -q.x, -q.y, -q.z}; }
static inline double mk_q_dot(mk_quat a, mk_quat b) { return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z; }
static inline mk_quat mk_q_norm(mk_quat q)
{
    double m = fmax(fmax(fabs(q.w), fabs(q.x)), fmax(fabs(q.y), fabs(q.z)));
    if (m > 0 && m < INFINITY) {
        int e = ilogb(m);
        q = (mk_quat){ldexp(q.w, -e), ldexp(q.x, -e), ldexp(q.y, -e), ldexp(q.z, -e)};
    }
    double l = sqrt(mk_q_dot(q, q));
    return l == 0 ? mk_q_identity() : (mk_quat){q.w / l, q.x / l, q.y / l, q.z / l}; /* 1.6.1: a NaN comes back NaN (1.6.0: identity) */
}
/* rotation by ang radians about axis (any length; a zero axis gives the identity, 1.6.0, for a finite angle: 1.6.2, NaN for a NaN or
   infinite one, as about any other axis) */
static inline mk_quat mk_q_axis_angle(mk_v3 axis, double ang)
{
    axis = mk_v3_norm(axis);
    if (axis.x == 0 && axis.y == 0 && axis.z == 0 && isfinite(ang)) return mk_q_identity();
    double s = sin(ang / 2);
    return (mk_quat){cos(ang / 2), axis.x * s, axis.y * s, axis.z * s};
}
/* rotation vector r = axis * angle (the exponential map); exact at every angle, normalized at small ones */
static inline mk_quat mk_q_from_rotvec(mk_v3 r)
{
    double a = mk_v3_len(r);
    if (a < 1e-12) return mk_q_norm((mk_quat){1, r.x / 2, r.y / 2, r.z / 2});
    return mk_q_axis_angle(mk_v3_scale(r, 1 / a), a);
}
/* v' = v + 2w (u x v) + 2 u x (u x v), u the vector part: no matrix */
static inline mk_v3 mk_q_rotate(mk_quat q, mk_v3 v)
{
    mk_v3 u = {q.x, q.y, q.z};
    mk_v3 t = mk_v3_scale(mk_v3_cross(u, v), 2);
#if defined(MK_MUTATE) && MK_MUTATE == 1
    return mk_v3_add(mk_v3_sub(v, mk_v3_scale(t, q.w)), mk_v3_cross(u, t));
#else
    return mk_v3_add(mk_v3_add(v, mk_v3_scale(t, q.w)), mk_v3_cross(u, t));
#endif
}
static inline mk_v3 mk_q_unrotate(mk_quat q, mk_v3 v) { return mk_q_rotate(mk_q_conj(q), v); }
/* rotation matrix, row-major: R[i][j] */
static inline void mk_q_to_mat(mk_quat q, double R[3][3])
{
    double w = q.w, x = q.x, y = q.y, z = q.z;
    R[0][0] = 1 - 2 * (y * y + z * z), R[0][1] = 2 * (x * y - w * z), R[0][2] = 2 * (x * z + w * y);
    R[1][0] = 2 * (x * y + w * z), R[1][1] = 1 - 2 * (x * x + z * z), R[1][2] = 2 * (y * z - w * x);
    R[2][0] = 2 * (x * z - w * y), R[2][1] = 2 * (y * z + w * x), R[2][2] = 1 - 2 * (x * x + y * y);
}
/* a rotation matrix to a unit quaternion (Shepperd: the largest of trace and diagonal picks the branch) */
static inline mk_quat mk_q_from_mat(const double R[3][3])
{
    double q[4], tr = R[0][0] + R[1][1] + R[2][2];
    if (tr > 0) {
        double s = 2 * sqrt(tr + 1);
        q[0] = s / 4, q[1] = (R[2][1] - R[1][2]) / s, q[2] = (R[0][2] - R[2][0]) / s, q[3] = (R[1][0] - R[0][1]) / s;
    } else {
        int i = R[1][1] > R[0][0] ? 1 : 0;
        if (R[2][2] > R[i][i]) i = 2;
        int j = (i + 1) % 3, k = (i + 2) % 3;
        double s = 2 * sqrt(1 + R[i][i] - R[j][j] - R[k][k]);
        q[0] = (R[k][j] - R[j][k]) / s;
        q[1 + i] = s / 4, q[1 + j] = (R[j][i] + R[i][j]) / s, q[1 + k] = (R[k][i] + R[i][k]) / s;
    }
    return mk_q_norm((mk_quat){q[0], q[1], q[2], q[3]});
}
/* the shorter arc from a (t = 0) to b (t = 1) */
static inline mk_quat mk_q_slerp(mk_quat a, mk_quat b, double t)
{
    double d = mk_q_dot(a, b);
    if (d < 0) b = (mk_quat){-b.w, -b.x, -b.y, -b.z}, d = -d;
    if (d > 0.9995)
        return mk_q_norm((mk_quat){a.w + t * (b.w - a.w), a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), a.z + t * (b.z - a.z)});
    double th = acos(d), s = sin(th), wa = sin((1 - t) * th) / s, wb = sin(t * th) / s;
    return (mk_quat){wa * a.w + wb * b.w, wa * a.x + wb * b.x, wa * a.y + wb * b.y, wa * a.z + wb * b.z};
}
/* q after turning at angular velocity omega (world axes) for dt: exact for constant omega */
static inline mk_quat mk_q_integrate(mk_quat q, mk_v3 omega, double dt)
{
    return mk_q_norm(mk_q_mul(mk_q_from_rotvec(mk_v3_scale(omega, dt)), q));
}

static inline mk_v3 mk_xf_apply(mk_xf a, mk_v3 p) { return mk_v3_add(mk_q_rotate(a.r, p), a.t); }
static inline mk_xf mk_xf_compose(mk_xf a, mk_xf b) { return (mk_xf){mk_q_mul(a.r, b.r), mk_v3_add(mk_q_rotate(a.r, b.t), a.t)}; }
static inline mk_xf mk_xf_inv(mk_xf a)
{
    mk_quat c = mk_q_conj(a.r);
    return (mk_xf){c, mk_v3_scale(mk_q_rotate(c, a.t), -1)};
}
#endif
