/* mk_hhd on known cases: a pinhole camera's motion field, from Longuet-Higgins and Prazdny 1980 (The interpretation of
   a moving retinal image, Proc. R. Soc. Lond. B 208), over a plane facing the camera, a tilted plane and a box room;
   the focus fit's residual; then the split's own properties: three parts orthogonal and adding up, each part's divergence
   and curl, a field built from known potentials split back into them (along either axis, and through Bluestein's
   transform), the refusals, the sum within 4 ulps and |cos| within 1e-10 at sizes from 2 x 2 to 3840 x 2160 (with
   OpenMP, the same bytes at 1 to 16 threads), and 1920 x 1080. The camera's numbers are arbitrary but
   plausible: a 256 x 144 frame, focal length 1000 px (a 14.6 x 8.2 degree view) or 200 px for the room (65 x 39
   degrees), walking forward at 2 m/s toward a wall 4 m off, heading 0.1 m/s right and 0.05 m/s up, rolling at 0.8 rad/s
   (a quick head tilt), panning (yaw) at 0.02 rad/s. The camera frame: x right, y down (image rows), z forward; the flow
   in pixels per second at pixel centres. */
#include "../include/mk_hhd.h"
#include "t.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    double f, U, V, W, A, B, C; /* focal length (px); translation (m/s); angular velocity about x, y, z: pitch, yaw, roll (rad/s) */
    int scene;                  /* 0 a plane facing the camera at Z0, 1 a tilted plane, 2 a box room */
    double Z0, p, q;            /* the tilted plane: 1 / Z = (1 + p x / f + q y / f) / Z0 */
} cam;
static double depth(const cam *c, double x, double y)
{
    if (c->scene == 0) return c->Z0;
    if (c->scene == 1) return c->Z0 / (1 + c->p * x / c->f + c->q * y / c->f);
    /* the room: back wall 6 m ahead, side walls 1.8 m left and 2.2 m right, floor 1.2 m down, ceiling 1.5 m up */
    double Z = 6.0, t;
    if (x > 0 && (t = 2.2 * c->f / x) < Z) Z = t;
    if (x < 0 && (t = -1.8 * c->f / x) < Z) Z = t;
    if (y > 0 && (t = 1.2 * c->f / y) < Z) Z = t;
    if (y < 0 && (t = -1.5 * c->f / y) < Z) Z = t;
    return Z;
}
/* Longuet-Higgins and Prazdny's equations; (x, y) from the principal point, the frame's centre */
static void flow(const cam *c, int nx, int ny, double *u, double *v)
{
    const double cx = (nx - 1) / 2.0, cy = (ny - 1) / 2.0, f = c->f;
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            const double x = i - cx, y = j - cy, Z = depth(c, x, y);
            u[j * nx + i] = (-f * c->U + x * c->W) / Z + c->A * x * y / f - c->B * (f + x * x / f) + c->C * y;
            v[j * nx + i] = (-f * c->V + y * c->W) / Z + c->A * (f + y * y / f) - c->B * x * y / f - c->C * x;
        }
}
static double dot(const double *ax, const double *ay, const double *bx, const double *by, size_t n)
{
    double s = 0;
    for (size_t k = 0; k < n; k++) s += ax[k] * bx[k] + ay[k] * by[k];
    return s;
}
/* the largest |<a, b>| / sum |v|^2 over the three pairs */
static double orth(const mk_hhd *H, const double *u, const double *v)
{
    const size_t n = (size_t)H->nx * H->ny;
    const double vv = dot(u, v, u, v, n);
    return fmax(fabs(dot(H->gx, H->gy, H->rx, H->ry, n)), fmax(fabs(dot(H->gx, H->gy, H->hx, H->hy, n)), fabs(dot(H->rx, H->ry, H->hx, H->hy, n)))) / vv;
}
/* the largest |box div| (what = 0) or |box curl| (1) of a field over the inner corners */
static double boxop(int nx, int ny, const double *u, const double *v, int what)
{
    double m = 0;
    for (int b = 1; b < ny; b++)
        for (int a = 1; a < nx; a++) {
            const int ll = (b - 1) * nx + a - 1, lr = ll + 1, ul = ll + nx, ur = ul + 1;
            const double d = what ? 0.5 * (v[lr] - v[ll] + v[ur] - v[ul]) - 0.5 * (u[ul] - u[ll] + u[ur] - u[lr])
                                  : 0.5 * (u[lr] - u[ll] + u[ur] - u[ul]) + 0.5 * (v[ul] - v[ll] + v[ur] - v[lr]);
            m = fmax(m, fabs(d));
        }
    return m;
}
static double absmax(const double *a, size_t n)
{
    double m = 0;
    for (size_t k = 0; k < n; k++) m = fmax(m, fabs(a[k]));
    return m;
}
static double now(void)
{
    struct timespec t;
    timespec_get(&t, TIME_UTC);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}
/* the focus fit's residual recomputed from the parts: sqrt(sum ((x - x0) v'y - (y - y0) v'x)^2 / sum |v'|^2), v' the
   field less a swirl at rate rot about the frame's centre, (fx, fy) the focus from the centre */
static double foe_rms_ref(const mk_hhd *H, double rot, double fx, double fy)
{
    const double cx = (H->nx - 1) / 2.0, cy = (H->ny - 1) / 2.0;
    double num = 0, den = 0;
    for (int j = 0; j < H->ny; j++)
        for (int i = 0; i < H->nx; i++) {
            const size_t k = (size_t)j * H->nx + i;
            const double x = i - cx, y = j - cy;
            const double vx = H->gx[k] + H->rx[k] + H->hx[k] + rot * y, vy = H->gy[k] + H->ry[k] + H->hy[k] - rot * x;
            const double m = (x - fx) * vy - (y - fy) * vx;
            num += m * m, den += vx * vx + vy * vy;
        }
    return sqrt(num / den);
}
enum { NX = 256, NY = 144 };
static double U[NX * NY], V[NX * NY];
/* split a camera's flow; the readout about the frame's centre */
static mk_hhd_read cam_read(mk_hhd *H, cam c)
{
    mk_hhd_read R;
    flow(&c, NX, NY, U, V);
    if (mk_hhd_split(H, U, V)) memset(&R, 0, sizeof R), R.expand = NAN;
    else mk_hhd_readout(H, NAN, NAN, &R);
    R.foe[0] -= (NX - 1) / 2.0, R.foe[1] -= (NY - 1) / 2.0; /* from the principal point, as the heading f (U, V) / W is */
    return R;
}

int main(void)
{
    mk_hhd H;
    CHECK(mk_hhd_init(&H, NX, NY) == 0, "mk_hhd_init 256 x 144");
    const cam base = {1000, 0.1, -0.05, 2.0, 0, 0, 0.8, 0, 4.0, 0, 0};
    const double hx = base.f * base.U / base.W, hy = base.f * base.V / base.W, a = base.W / base.Z0, C = base.C, B = 0.02;

    /* the plane facing the camera: forward, heading off centre, roll */
    mk_hhd_read R = cam_read(&H, base);
    CHECK(hypot(R.foe[0] - hx, R.foe[1] - hy) < 1e-6, "plane, forward + roll: the focus of expansion at the heading (%.6f, %.6f) px, f (U, V) / W = (%.1f, %.1f)", R.foe[0], R.foe[1], hx, hy);
    CHECK(fabs(R.expand - a) < 1e-12 * a, "plane: the expansion rate W / Z = %.3f 1/s (read %.15f)", a, R.expand);
    CHECK(fabs(R.rot + C) < 1e-12, "plane: the rotational part's rate is the roll, -C = %.2f rad/s in (x right, y down) (read %.15f)", -C, R.rot);
    CHECK(fabs(R.trans[0] + base.f * base.U / base.Z0) < 1e-9 && fabs(R.trans[1] + base.f * base.V / base.Z0) < 1e-9,
          "plane: the harmonic part's mean is the heading's sideways flow -f (U, V) / Z = (%.1f, %.1f) px/s (read %.6f, %.6f)", -base.f * base.U / base.Z0,
          -base.f * base.V / base.Z0, R.trans[0], R.trans[1]);
    CHECK(orth(&H, U, V) < 1e-12, "plane: the three parts orthogonal (worst |<a, b>| / |v|^2 %.2g)", orth(&H, U, V));
    const double q0 = foe_rms_ref(&H, R.rot, R.foe[0], R.foe[1]);
    CHECK(R.foe_rms < 1e-6 && fabs(R.foe_rms - q0) < 1e-9, "plane, forward + roll: the focus fit's residual foe_rms %.2g px, its lines all through the focus (the test's own sum %.2g)", R.foe_rms, q0);
    {
        const mk_hhd_read R0 = R;
        cam c = base;
        c.B = B, R = cam_read(&H, c);
        /* the pan's flow, -f B - B x^2 / f across, has no divergence or curl on average over a frame centred on the axis;
           its near-uniform part moves the flow's zero by f B Z / W to first order, as one plane's flow cannot tell a pan
           from a turned heading (Longuet-Higgins and Prazdny) */
        const double sh = c.f * B * c.Z0 / c.W;
        CHECK(fabs(R.expand - R0.expand) < 1e-12 && fabs(R.rot - R0.rot) < 1e-12, "plane + yaw %.2f rad/s: the expansion and roll rates unchanged (%.2g, %.2g)", B,
              R.expand - R0.expand, R.rot - R0.rot);
        CHECK(fabs(R.foe[0] - (hx + sh)) < 0.1 && fabs(R.foe[1] - hy) < 0.1,
              "plane + yaw: the focus read moves to heading + f B Z / W = (%.1f, %.1f) px (read %.3f, %.3f): physics, not separable in one plane", hx + sh, hy, R.foe[0], R.foe[1]);
        const double mx2 = (NX * NX - 1) / 12.0; /* the mean x^2 over the frame's pixel centres */
        CHECK(fabs(R.trans[0] - R0.trans[0] + c.f * B * (1 + mx2 / (c.f * c.f))) < 1e-9 && fabs(R.trans[1] - R0.trans[1]) < 1e-9,
              "plane + yaw: the pan lands in the harmonic part's mean, -f B (1 + <x^2> / f^2) = %.4f px/s (read %.4f)", -c.f * B * (1 + mx2 / (c.f * c.f)), R.trans[0] - R0.trans[0]);
        const double q = foe_rms_ref(&H, R.rot, R.foe[0], R.foe[1]);
        CHECK(R.foe_rms > 0.05 && fabs(R.foe_rms - q) < 1e-9 * R.foe_rms,
              "plane + yaw: a residual foe_rms %.3f px against %.2g px without the pan, while the focus moved %.0f px: the pan's near-uniform part leaves none, only its second-order part (the test's own sum %.3f)",
              R.foe_rms, R0.foe_rms, sh, q);
    }
    /* where each motion lands, alone: the parts' shares of its energy */
    double froll[3], ffwd[3];
    {
        cam c = base;
        c.U = c.V = c.W = 0, R = cam_read(&H, c), memcpy(froll, R.e, sizeof froll);
        CHECK(fabs(R.rot + C) < 1e-12 && R.e[0] < 1e-20, "roll alone: rate %.15f; shares gradient %.2g, rotational %.4f, harmonic %.4f (the swirl crossing the frame is harmonic)", R.rot,
              R.e[0], R.e[1], R.e[2]);
        CHECK(isnan(R.foe[0]), "roll alone: no focus of expansion (NaN)");
        c = base, c.U = c.V = 0, c.C = 0, R = cam_read(&H, c), memcpy(ffwd, R.e, sizeof ffwd);
        CHECK(hypot(R.foe[0], R.foe[1]) < 1e-9 && R.e[1] < 1e-20, "forward alone: focus (%.2g, %.2g), shares gradient %.4f, rotational %.2g, harmonic %.4f", R.foe[0], R.foe[1], R.e[0], R.e[1],
              R.e[2]);
        CHECK(fabs(ffwd[0] - froll[1]) < 1e-9 && fabs(ffwd[2] - froll[2]) < 1e-9,
              "the frame treats an expansion and a swirl alike: forward's gradient share %.6f = roll's rotational share %.6f", ffwd[0], froll[1]);
        static const double fs[3] = {1000, 300, 150}, lo[3] = {0.999, 0.99, 0.95};
        double prev = 1;
        for (int k = 0; k < 3; k++) {
            c = base, c.U = c.V = c.W = 0, c.C = 0, c.B = B, c.f = fs[k], R = cam_read(&H, c);
            CHECK(R.e[2] > lo[k] && R.e[2] <= prev, "yaw alone, f %.0f px (%.0f degrees across): harmonic share %.5f, gradient %.5f, rotational %.5f (over %.3f, and less as the view widens)",
                  fs[k], 2 * atan(NX / 2.0 / fs[k]) * 180 / T_PI, R.e[2], R.e[0], R.e[1], lo[k]);
            prev = R.e[2];
        }
        c = base, c.U = 0, c.W = 0, c.C = 0, c.V = 0.2, R = cam_read(&H, c);
        CHECK(R.e[2] > 1 - 1e-12 && isnan(R.foe[0]) && fabs(R.trans[1] + c.f * c.V / c.Z0) < 1e-9,
              "sliding sideways alone: all harmonic (%.15f), no focus, mean %.1f px/s", R.e[2], R.trans[1]);
    }
    /* the tilted plane: the translation's flow over a slanted depth has curl of its own, W (q x0 - p y0) / (2 f Z0) on average,
       and a divergence other than 2 W / Z0 */
    {
        cam c = base;
        c.scene = 1, c.p = 0.8, c.q = -0.5, R = cam_read(&H, c);
        const double e = a * (1 - (c.p * hx + c.q * hy) / (2 * c.f)), r = -C + c.W * (c.q * hx - c.p * hy) / (2 * c.f * c.Z0);
        CHECK(fabs(R.expand - e) < 1e-12, "tilted plane: the expansion rate is the frame's mean, (W / Z0) (1 - (p x0 + q y0) / 2 f) = %.6f (read %.15f)", e, R.expand);
        CHECK(fabs(R.rot - r) < 1e-12, "tilted plane: the rotational rate is the roll plus the slant's shear, %.6f (read %.15f)", r, R.rot);
        CHECK(hypot(R.foe[0] - hx, R.foe[1] - hy) < 0.5, "tilted plane: the focus within %.3f px of the heading (that shear left in the swirl taken out)",
              hypot(R.foe[0] - hx, R.foe[1] - hy));
        CHECK(orth(&H, U, V) < 1e-12, "tilted plane: the parts orthogonal (%.2g)", orth(&H, U, V));
    }
    /* the box room, 65 degrees across: walls, floor and ceiling at their own depths, creases between them */
    {
        cam c = base;
        c.scene = 2, c.f = 200, c.U = 0.02, c.V = -0.01, R = cam_read(&H, c);
        const double gx = c.f * c.U / c.W, gy = c.f * c.V / c.W;
        CHECK(hypot(R.foe[0] - gx, R.foe[1] - gy) < 0.5, "box room, forward + roll: the focus within %.3f px of the heading (%.1f, %.1f)", hypot(R.foe[0] - gx, R.foe[1] - gy), gx, gy);
        CHECK(fabs(R.rot + C) < 0.02 * C, "box room: the rotational rate %.6f within 2 %% of -C (the room's own shear in it, %.2f %%)", R.rot, 100 * (R.rot + C) / C);
        CHECK(orth(&H, U, V) < 1e-12, "box room: the parts orthogonal (%.2g)", orth(&H, U, V));
        const size_t n = (size_t)NX * NY;
        double sum = 0;
        for (size_t k = 0; k < n; k++) sum = fmax(sum, fmax(fabs(H.gx[k] + H.rx[k] + H.hx[k] - U[k]), fabs(H.gy[k] + H.ry[k] + H.hy[k] - V[k])));
        CHECK(sum < 1e-13 * fmax(absmax(U, n), absmax(V, n)), "box room: the parts add up to the field (worst %.2g px/s)", sum);
        const double sc = fmax(absmax(U, n), absmax(V, n));
        const double dr = boxop(NX, NY, H.rx, H.ry, 0), cg = boxop(NX, NY, H.gx, H.gy, 1), dh = boxop(NX, NY, H.hx, H.hy, 0), ch = boxop(NX, NY, H.hx, H.hy, 1);
        CHECK(dr < 1e-13 * sc && cg < 1e-13 * sc, "box room: the rotational part has no divergence (%.2g), the gradient part no curl (%.2g)", dr, cg);
        CHECK(dh < 1e-10 * sc && ch < 1e-10 * sc, "box room: the harmonic part has neither (%.2g, %.2g)", dh, ch);
    }
    mk_hhd_free(&H);

    /* known potentials split back: phi0, psi0 random at the inner corners, h0 = a uniform shift + a pure strain (x, -y) / 100;
       37 x 23 takes a pass of radix 23 by its pairs (2 ny = 46), 64 x 37 its transforms along x (2 nx = 128) and its
       solves along y, 67 x 71 Bluestein's chirp z (2 ny = 142 = 2 71) */
    static const int sz[4][2] = {{37, 23}, {64, 37}, {67, 71}, {NX, NY}};
    for (int t = 0; t < 4; t++) {
        const int nx = sz[t][0], ny = sz[t][1], nv1 = nx + 1;
        mk_hhd K;
        if (mk_hhd_init(&K, nx, ny)) { CHECK(0, "mk_hhd_init %d x %d", nx, ny); continue; }
        double *p0 = calloc((size_t)nv1 * (ny + 1), sizeof(double)), *s0 = calloc((size_t)nv1 * (ny + 1), sizeof(double));
        double *u = malloc(sizeof(double) * (size_t)nx * ny), *v = malloc(sizeof(double) * (size_t)nx * ny);
        for (int b = 1; b < ny; b++)
            for (int a = 1; a < nx; a++) p0[b * nv1 + a] = t_rand(), s0[b * nv1 + a] = t_rand();
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                const int c = j * nx + i, o = j * nv1 + i, p = o + nv1;
                const double fx = 0.5 * (p0[o + 1] - p0[o] + p0[p + 1] - p0[p]), fy = 0.5 * (p0[p] - p0[o] + p0[p + 1] - p0[o + 1]);
                const double sx = 0.5 * (s0[o + 1] - s0[o] + s0[p + 1] - s0[p]), sy = 0.5 * (s0[p] - s0[o] + s0[p + 1] - s0[o + 1]);
                u[c] = fx + sy + 0.3 + 0.01 * i, v[c] = fy - sx - 0.2 - 0.01 * j;
            }
        double ep = INFINITY, es = INFINITY, eh = INFINITY;
        if (!mk_hhd_split(&K, u, v)) {
            ep = es = eh = 0;
            for (int k = 0; k < nv1 * (ny + 1); k++) ep = fmax(ep, fabs(K.phi[k] - p0[k])), es = fmax(es, fabs(K.psi[k] - s0[k]));
            for (int j = 0; j < ny; j++)
                for (int i = 0; i < nx; i++) {
                    const int c = j * nx + i;
                    eh = fmax(eh, fmax(fabs(K.hx[c] - 0.3 - 0.01 * i), fabs(K.hy[c] + 0.2 + 0.01 * j)));
                }
        }
        CHECK(ep < 1e-12 && es < 1e-12 && eh < 1e-12, "%d x %d: a field made of known potentials and a harmonic shift + strain splits back into them (phi %.2g, psi %.2g, h %.2g)", nx, ny, ep, es,
              eh);
        free(p0), free(s0), free(u), free(v), mk_hhd_free(&K);
    }
    /* refusals */
    {
        mk_hhd K;
        CHECK(mk_hhd_init(&K, 1, 10) == -1 && mk_hhd_init(&K, 10, 1) == -1 && K.gx == NULL, "mk_hhd_init refuses a frame 1 pixel across");
        CHECK(mk_hhd_split(&K, U, V) == -1, "mk_hhd_split refuses an uninitialized frame");
        mk_hhd_init(&K, 8, 6);
        double u[48], v[48];
        for (int k = 0; k < 48; k++) u[k] = v[k] = 1;
        u[17] = NAN;
        CHECK(mk_hhd_split(&K, u, v) == -1 && isnan(K.gx[0]) && isnan(K.hy[47]) && isnan(K.phi[0]), "mk_hhd_split refuses a NaN in the field, every part NaN");
        mk_hhd_free(&K);
    }
    /* from 2 x 2 to 3840 x 2160, three fields (uniform random, the camera's forward + roll over the wall, smooth waves):
       the parts add up to the field within 4 ulps of each pixel's largest value, summed g + r + h; |cos| between any two
       parts within 1e-10; with OpenMP, the same bytes at 1, 2, 3, 7 and 16 threads */
    {
        static const int fs[][2] = {{2, 2}, {3, 2}, {2, 3}, {37, 23}, {64, 37}, {67, 71}, {256, 144}, {1920, 1031}, {1920, 1080}, {3840, 2160}};
        static const char *kn[3] = {"random", "camera", "smooth"};
        enum { NF = sizeof fs / sizeof fs[0] };
        double su = 0, cm = 0;
        int sw = 0, cw = 0, ckw = 0, okp = 1;
#ifdef _OPENMP
        int same = 1, nth = 0;
#endif
        for (int t = 0; t < NF; t++)
            for (int kind = 0; kind < 3; kind++) {
                const int nx = fs[t][0], ny = fs[t][1];
                const size_t n = (size_t)nx * ny;
                double *u = malloc(sizeof(double) * n), *v = malloc(sizeof(double) * n);
                if (kind == 1) flow(&base, nx, ny, u, v);
                else
                    for (int j = 0; j < ny; j++)
                        for (int i = 0; i < nx; i++)
                            u[(size_t)j * nx + i] = kind ? sin(0.01 * i) + 0.3 * cos(0.02 * j) : t_rand(), v[(size_t)j * nx + i] = kind ? cos(0.0091 * i + 0.01 * j) : t_rand();
                mk_hhd K;
                if (mk_hhd_init(&K, nx, ny) || mk_hhd_split(&K, u, v)) {
                    okp = 0, mk_hhd_free(&K), free(u), free(v);
                    continue;
                }
                double gg = 0, rr = 0, hh = 0, gr = 0, gh = 0, rh = 0;
                for (size_t k = 0; k < n; k++) {
                    const double p[2][4] = {{u[k], K.gx[k], K.rx[k], K.hx[k]}, {v[k], K.gy[k], K.ry[k], K.hy[k]}};
                    for (int d = 0; d < 2; d++) {
                        const double m = fmax(fmax(fabs(p[d][0]), fabs(p[d][1])), fmax(fabs(p[d][2]), fabs(p[d][3])));
                        int e;
                        frexp(m, &e);
                        const double s = m > 0 ? fabs(p[d][1] + p[d][2] + p[d][3] - p[d][0]) / ldexp(1.0, e - 53) : 0;
                        if (s > su) su = s, sw = t;
                    }
                    gg += K.gx[k] * K.gx[k] + K.gy[k] * K.gy[k], rr += K.rx[k] * K.rx[k] + K.ry[k] * K.ry[k], hh += K.hx[k] * K.hx[k] + K.hy[k] * K.hy[k];
                    gr += K.gx[k] * K.rx[k] + K.gy[k] * K.ry[k], gh += K.gx[k] * K.hx[k] + K.gy[k] * K.hy[k], rh += K.rx[k] * K.hx[k] + K.ry[k] * K.hy[k];
                }
                const double c3[3] = {gg > 0 && rr > 0 ? fabs(gr) / sqrt(gg * rr) : 0, gg > 0 && hh > 0 ? fabs(gh) / sqrt(gg * hh) : 0, rr > 0 && hh > 0 ? fabs(rh) / sqrt(rr * hh) : 0};
                for (int q = 0; q < 3; q++)
                    if (!(c3[q] <= cm)) cm = c3[q], cw = t, ckw = kind;
#ifdef _OPENMP
                if (kind == 0) { /* the bytes at other thread counts against these */
                    static const int th[4] = {2, 3, 7, 16};
                    const size_t nb = sizeof(double) * (6 * n + 2 * (size_t)(nx + 1) * (ny + 1));
                    double *ref = malloc(nb);
                    const int t0 = omp_get_max_threads();
                    memcpy(ref, K.gx, nb);
                    for (int q = 0; q < 4; q++) {
                        omp_set_num_threads(th[q]);
                        same &= !mk_hhd_split(&K, u, v) && !memcmp(ref, K.gx, nb);
                    }
                    omp_set_num_threads(t0), free(ref), nth++;
                }
#endif
                mk_hhd_free(&K), free(u), free(v);
            }
        CHECK(okp && su <= 4, "%d sizes from 2 x 2 to 3840 x 2160, three fields: g + r + h is the field within 4 ulps of each pixel's largest value (worst %.2f, at %d x %d)", NF, su, fs[sw][0],
              fs[sw][1]);
        CHECK(okp && cm <= 1e-10, "the same splits: |cos| between any two parts <= 1e-10 (worst %.2g, at %d x %d, the %s field)", cm, fs[cw][0], fs[cw][1], kn[ckw]);
#ifdef _OPENMP
        CHECK(okp && same, "the same sizes, random field: the same bytes at 1, 2, 3, 7 and 16 threads (%d sizes)", nth);
#endif
    }
    /* 1920 x 1080: random flow, orthogonal parts, the time */
    {
        enum { BX = 1920, BY = 1080 };
        mk_hhd K;
        double *u = malloc(sizeof(double) * BX * BY), *v = malloc(sizeof(double) * BX * BY), dt = INFINITY, o = INFINITY;
        for (int k = 0; k < BX * BY; k++) u[k] = t_rand(), v[k] = t_rand();
        if (!mk_hhd_init(&K, BX, BY)) {
            const double t0 = now();
            const int r = mk_hhd_split(&K, u, v);
            dt = now() - t0;
            if (!r) o = orth(&K, u, v);
            mk_hhd_free(&K);
        }
        CHECK(o < 1e-13, "1920 x 1080 random flow: parts orthogonal (%.2g); split in %.3f s, the time not checked (it depends on the machine's load)", o, dt);
        free(u), free(v);
    }
    return wrong;
}
