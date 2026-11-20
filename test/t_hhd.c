/* mk_hhd on known cases (1.5.6): a pinhole camera's motion field, from Longuet-Higgins and Prazdny 1980 (The interpretation of
   a moving retinal image, Proc. R. Soc. Lond. B 208), over a plane facing the camera, a tilted plane and a box room; then the
   split's own properties: three parts orthogonal and adding up, each part's divergence and curl, a field built from known
   potentials split back into them, the refusals, and 1920 x 1080. The camera's numbers are MADE UP but plausible: a 256 x 144
   frame, focal length 1000 px (a 14.6 x 8.2 degree view) or 200 px for the room (65 x 39 degrees), walking forward at
   2 m/s toward a wall 4 m off, heading 0.1 m/s right and 0.05 m/s up, rolling at 0.8 rad/s (a quick head tilt), panning
   (yaw) at 0.02 rad/s. The camera frame: x right, y down (image rows), z forward; the flow in pixels per second at pixel
   centres. The picture of the main row: scratch/hhd-work/hhdpic.c. */
#include "../include/mk_hhd.h"
#include "t.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    double f, U, V, W, A, B, C; /* focal length (px); translation (m/s); angular velocity about x, y, z: pitch, yaw, roll (rad/s) */
    int scene;                  /* 0 a plane facing the camera at Z0, 1 a tilted plane, 2 a box room */
    double Z0, p, q;            /* the tilted plane: 1 / Z = (1 + p x / f + q y / f) / Z0 */
} cam;
static double depth(const cam *c, double x, double y)
{
    if (c->scene == 0) return c->Z0;
    if (c->scene == 1) return c->Z0 / (1 + c->p * x / c->f + c->q * y / c->f);
    double Z = 6.0, t; /* the room (MADE UP): back wall 6 m ahead, side walls 1.8 m left and 2.2 m right, floor 1.2 m down, ceiling 1.5 m up */
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
       37 x 23 takes the transform's direct-sum radix (2 ny = 46 = 2 23) */
    static const int sz[2][2] = {{37, 23}, {NX, NY}};
    for (int t = 0; t < 2; t++) {
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
        CHECK(o < 1e-13, "1920 x 1080 random flow: parts orthogonal (%.2g); split in %.3f s, not checked (1.6.1: a loaded host failed 1.6.0 here, review 2 #7)", o, dt);
        free(u), free(v);
    }
    return wrong;
}
