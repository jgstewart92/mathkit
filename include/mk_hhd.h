/* mathkit: the Helmholtz-Hodge decomposition of a 2-D motion field on an image grid. Optic flow or motion vectors
   (u, v), one per pixel, split into three parts that add up to the field and are orthogonal in the plain pixel sum:
     gradient part    grad(phi)                        expansion and contraction; curl-free
     rotational part  rot(psi) = (dpsi/dy, -dpsi/dx)   swirl; divergence-free
     harmonic rest    h                                both divergence- and curl-free: a uniform shift, a pure strain,
                                                       and what comes in through the frame's edge
   Boundary choice: phi = psi = 0 on the frame, Dirichlet on both (the Hodge-Morrey split; "normal-parallel" in Bhatia
   et al., The Helmholtz-Hodge decomposition: a survey, IEEE TVCG 2013). The gradient part meets the frame at right
   angles, the rotational part runs along it, the harmonic part takes the rest with no condition. On a rectangle this is
   the standard choice that gives three orthogonal parts and puts a uniform shift in the harmonic one: with phi Neumann
   (dphi/dn the field's normal flow) the harmonic part could carry no flux through the frame, which on a simply
   connected frame leaves no harmonic field at all, so a shift would land in the gradient part. The cost: near the
   frame, part of any expansion or swirl (the flow crossing the frame) is harmonic. A swirl or an expansion centred on a
   16:9 frame keeps 62 % of its energy in its own part, 38 % harmonic; on a square, 84 % (test/t_hhd.c measures it). The
   rates below are exact all the same: div(grad phi) = div v and curl(rot psi) = curl v at every inner corner.
   Discrete form: pixel (i, j) centred at (i, j); the potentials at the pixel corners (i - 1/2, j - 1/2), 0 on the
   frame's corners; a pixel's gradient is the 2 x 2 box difference of its four corners. So div(rot psi) and curl(grad
   phi) vanish at every inner corner to rounding, the parts are orthogonal to rounding, and phi and psi are
   least-squares fits: A phi = -div v and A psi = curl v at the inner corners (div and curl the box differences there, A
   the box Laplacian). A's sine modes have the values 2 (1 - cos(pi k / nx) cos(pi l / ny)): the solve is a sine
   transform along y (mk_fftm, length 2 ny) and the tridiagonal (-cos(pi l / ny), 2, -cos(pi l / ny)) along x for each
   mode (mk_tridiag). The box differences are exact on quadratic fields, which is what a pinhole camera's flow over a
   plane is (Longuet-Higgins and Prazdny 1980). A pattern that alternates pixel by pixel has no box divergence or curl,
   so such noise stays in the harmonic part.
   mk_hhd_init(H, nx, ny): nx, ny >= 2 (and their product under 2^31); 0, or -1 (no memory; H is then empty).
   mk_hhd_free.
   mk_hhd_split(H, u, v): u, v nx * ny, index j nx + i (x fastest), any units (pixels per frame, per second); fills H's
   parts gx, gy, rx, ry, hx, hy (the same layout) and phi, psi ((nx + 1) (ny + 1), index j (nx + 1) + i for corner
   (i - 1/2, j - 1/2)). 0, or -1 (H not initialized, a value not finite: every part NaN). Threads (with -fopenmp): the
   transforms and solves, by columns and rows, the same bits at any thread count.
   mk_hhd_readout(H, cx, cy, R), after a split; (cx, cy) the principal point in pixels (NaN: the frame's centre):
     expand   half the mean divergence of the gradient part over the inner corners: a(x - x0) gives a; an approach at
              speed W to a wall at Z gives W / Z, 1 / the time to contact. Per the field's time unit.
     rot      half the mean curl (dv/dx - du/dy) of the rotational part: positive turns +x toward +y (clockwise on
              screen when y counts rows down). A camera's roll at rate C about its forward axis (Longuet-Higgins and
              Prazdny's C; x right, y down, z forward) gives -C. A translation over unequal depths has a shear of its
              own, which counts here too.
     trans    the harmonic part's mean, which is the field's mean: the box gradients of rim-zero potentials sum to 0.
     foe      the focus of expansion: the point x0 with the least squares of (x - x0) x v' over the pixels, v' the field
              less a rigid swirl at rate rot about (cx, cy). The flow of a translating camera points away from x0
              whatever the depth, so x0 is exact for a translation and a roll over a wall facing the camera; another
              scene's shear in rot moves it a little (the tests: 0.24 px on a tilted plane, 0.12 px in a box room). A
              pan or tilt adds a near-uniform shift that one view cannot tell from a moved heading, and the read point
              moves with it. NaN when the field has no point to fit (parallel flow, or nothing left once the swirl is
              out).
     e[3]     the gradient, rotational and harmonic parts' shares of the field's energy (sum |part|^2 / sum |v|^2; they
              add to 1).
   Header only, C23.
   Experimental: its results may change between releases.
   Speed: the sine transform of length 2 ny takes its largest prime factor by a direct sum (mk_fftm). On one thread of
   an x86-64 desktop, 1920 x 1080 (factor 5) splits in 0.18 s, x 1088 (17) in 0.29 s, x 802 (401) in 3.3 s and x 1031
   (1031) in 10.7 s; 26 common heights from 144 to 4320 (letterboxed 804 and 816 among them) have factors up to 67. */
#ifndef MK_HHD_H
#define MK_HHD_H
#include "mk_fft.h"
#include "mk_linalg.h"
#include "mk_par.h"
#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
_Static_assert(sizeof(long) == 8, "mk_hhd.h needs LP64: its grid size check multiplies in a long");

/* lanes: buffers fixed by count, not by thread; columns per transform block */
enum { MK_HHD_LANES_ = 16, MK_HHD_COLS_ = 8 };

typedef struct {
    int nx, ny;
    double *gx, *gy, *rx, *ry, *hx, *hy;
    double *phi, *psi;
    mk_fftm f_;           /* length 2 ny */
    double *zr_, *zi_;    /* (nx - 1) (ny - 1) inner corners: -div then phi; curl then psi */
    double complex *buf_; /* per lane: MK_HHD_COLS_ columns of 2 ny, then a work line of 2 ny */
    double *tri_;         /* per lane: a, b, c of nx - 1 */
} mk_hhd;

typedef struct {
    double foe[2], expand, rot, trans[2], e[3];
} mk_hhd_read;

static inline void mk_hhd_free(mk_hhd *H)
{
    free(H->gx), free(H->zr_), free(H->buf_), free(H->tri_), mk_fftm_free(&H->f_);
    memset(H, 0, sizeof *H);
}
static inline int mk_hhd_init(mk_hhd *H, int nx, int ny)
{
    memset(H, 0, sizeof *H);
    if (nx < 2 || ny < 2 || nx > (1 << 28) || ny > (1 << 27) || (long)(nx + 1) * (ny + 1) > INT_MAX) return -1;
    const size_t nc = (size_t)nx * (size_t)ny, nv = (size_t)(nx + 1) * (size_t)(ny + 1), ni = (size_t)(nx - 1) * (size_t)(ny - 1);
    H->gx = malloc(sizeof(double) * (6 * nc + 2 * nv));
    H->zr_ = malloc(sizeof(double) * 2 * ni);
    H->buf_ = malloc(sizeof(double complex) * (size_t)MK_HHD_LANES_ * (MK_HHD_COLS_ + 1) * 2 * (size_t)ny);
    H->tri_ = malloc(sizeof(double) * (size_t)MK_HHD_LANES_ * 3 * (size_t)(nx - 1));
    if (!H->gx || !H->zr_ || !H->buf_ || !H->tri_ || mk_fftm_plan(&H->f_, 2 * ny)) { mk_hhd_free(H); return -1; }
    H->gy = H->gx + nc, H->rx = H->gy + nc, H->ry = H->rx + nc, H->hx = H->ry + nc, H->hy = H->hx + nc;
    H->phi = H->hy + nc, H->psi = H->phi + nv, H->zi_ = H->zr_ + ni;
    H->nx = nx, H->ny = ny;
    return 0;
}
/* the sine transform along y (k = 1 .. ny - 1) of columns a0 .. a1 - 1 of the inner corners, real and imaginary parts
   at once; times sc */
static inline void mk_hhd_dst_(mk_hhd *H, int lane, int a0, int a1, double sc)
{
    const int Nx = H->nx - 1, ny = H->ny, L = 2 * ny;
    double complex *col = H->buf_ + (size_t)lane * (MK_HHD_COLS_ + 1) * (size_t)L, *w = col + (size_t)MK_HHD_COLS_ * L;
    for (int a = a0; a < a1; a++) col[(size_t)(a - a0) * L] = col[(size_t)(a - a0) * L + ny] = 0;
    for (int b = 1; b < ny; b++) /* each column made odd about 0 and ny */
        for (int a = a0; a < a1; a++) {
            const size_t q = (size_t)(b - 1) * Nx + a;
            double complex *y = col + (size_t)(a - a0) * L;
            y[b] = CMPLX(H->zr_[q], H->zi_[q]), y[L - b] = CMPLX(-H->zr_[q], -H->zi_[q]);
        }
    for (int a = a0; a < a1; a++) mk_fftm_c(&H->f_, col + (size_t)(a - a0) * L, w, 0);
    for (int k = 1; k < ny; k++) /* the odd sequence's transform is -2 i times the sine sums */
        for (int a = a0; a < a1; a++) {
            const size_t q = (size_t)(k - 1) * Nx + a;
            const double complex y = col[(size_t)(a - a0) * L + k];
            H->zr_[q] = -cimag(y) * (sc / 2), H->zi_[q] = creal(y) * (sc / 2);
        }
}
/* the tridiagonal solves along x for modes l0 .. l1 - 1 (rows l - 1); 0 or -1 */
static inline int mk_hhd_tri_(mk_hhd *H, int lane, int l0, int l1)
{
    const int Nx = H->nx - 1;
    double *ta = H->tri_ + (size_t)lane * 3 * Nx, *tb = ta + Nx, *tc = tb + Nx;
    int bad = 0;
    for (int l = l0; l < l1; l++) {
#if defined(MK_MUTATE) && MK_MUTATE == 241
        const double c = cos(MK_FFT_PI * l / H->ny); /* planted fault: the coupling's sign turned */
#else
        const double c = -cos(MK_FFT_PI * l / H->ny);
#endif
        for (int part = 0; part < 2; part++) {
            double *row = (part ? H->zi_ : H->zr_) + (size_t)(l - 1) * Nx;
            for (int i = 0; i < Nx; i++) ta[i] = tc[i] = c, tb[i] = 2;
            bad |= mk_tridiag(Nx, ta, tb, tc, row, row);
        }
    }
    return bad ? -1 : 0;
}
static inline int mk_hhd_split(mk_hhd *H, const double *u, const double *v)
{
    if (!H->gx) return -1;
    const int nx = H->nx, ny = H->ny, Nx = nx - 1, Ny = ny - 1, nv1 = nx + 1;
    const size_t nc = (size_t)nx * ny;
    int fin = 1;
    /* -div and curl at the inner corners: cells LL (a-1, b-1), LR (a, b-1), UL (a-1, b), UR (a, b) */
    for (int b = 1; b < ny; b++)
        for (int a = 1; a < nx; a++) {
            const size_t ll = (size_t)(b - 1) * nx + (a - 1), lr = ll + 1, ul = ll + nx, ur = ul + 1, q = (size_t)(b - 1) * Nx + (a - 1);
            const double div = 0.5 * ((u[lr] - u[ll]) + (u[ur] - u[ul])) + 0.5 * ((v[ul] - v[ll]) + (v[ur] - v[lr]));
            const double curl = 0.5 * ((v[lr] - v[ll]) + (v[ur] - v[ul])) - 0.5 * ((u[ul] - u[ll]) + (u[ur] - u[lr]));
            H->zr_[q] = -div;
#if defined(MK_MUTATE) && MK_MUTATE == 245
            H->zi_[q] = -curl; /* planted fault: the curl's sign turned */
#else
            H->zi_[q] = curl;
#endif
        }
    for (size_t k = 0; k < nc; k++) fin &= isfinite(u[k]) && isfinite(v[k]);
    if (!fin) goto bad;
    const int nblk = (Nx + MK_HHD_COLS_ - 1) / MK_HHD_COLS_;
    int fail[MK_HHD_LANES_] = {0};
    for (int pass = 0; pass < 3; pass++) {
        MK_OMP_(omp parallel for schedule(static))
        for (int lane = 0; lane < MK_HHD_LANES_; lane++) {
            if (pass == 1) { /* the modes in rows, spread over the lanes */
                const int per = (Ny + MK_HHD_LANES_ - 1) / MK_HHD_LANES_, l0 = 1 + lane * per, l1 = l0 + per < ny ? l0 + per : ny;
                if (l0 < l1 && mk_hhd_tri_(H, lane, l0, l1)) fail[lane] = 1;
                continue;
            }
#if defined(MK_MUTATE) && MK_MUTATE == 242
            const double sc = pass ? 1.0 / ny : 1; /* planted fault: the inverse scaled by 1 / ny, not 2 / ny */
#else
            const double sc = pass ? 2.0 / ny : 1; /* the sine transform twice is ny / 2 times the identity */
#endif
            for (int k = lane; k < nblk; k += MK_HHD_LANES_) {
                const int a0 = k * MK_HHD_COLS_, a1 = a0 + MK_HHD_COLS_ < Nx ? a0 + MK_HHD_COLS_ : Nx;
                mk_hhd_dst_(H, lane, a0, a1, sc);
            }
        }
    }
    for (int lane = 0; lane < MK_HHD_LANES_; lane++)
        if (fail[lane]) goto bad;
    for (int b = 0; b <= ny; b++) /* the potentials, 0 on the rim */
        for (int a = 0; a <= nx; a++) {
            const size_t o = (size_t)b * nv1 + a;
            if (a == 0 || b == 0 || a == nx || b == ny) H->phi[o] = H->psi[o] = 0;
            else {
                const size_t q = (size_t)(b - 1) * Nx + (a - 1);
                H->phi[o] = H->zr_[q], H->psi[o] = H->zi_[q];
            }
        }
    fin = 1;
    /* the parts: each pixel's box differences of its corners (i, j), (i+1, j), (i, j+1), (i+1, j+1) */
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            const size_t c = (size_t)j * nx + i, o = (size_t)j * nv1 + i, p = o + nv1;
            const double *f = H->phi, *s = H->psi;
            const double fx = 0.5 * ((f[o + 1] - f[o]) + (f[p + 1] - f[p])), fy = 0.5 * ((f[p] - f[o]) + (f[p + 1] - f[o + 1]));
            const double sx = 0.5 * ((s[o + 1] - s[o]) + (s[p + 1] - s[p])), sy = 0.5 * ((s[p] - s[o]) + (s[p + 1] - s[o + 1]));
            H->gx[c] = fx, H->gy[c] = fy;
#if defined(MK_MUTATE) && MK_MUTATE == 240
            H->rx[c] = sy, H->ry[c] = sx; /* planted fault: rot psi's second component without its sign */
#else
            H->rx[c] = sy, H->ry[c] = -sx;
#endif
            H->hx[c] = u[c] - fx - H->rx[c], H->hy[c] = v[c] - fy - H->ry[c];
            fin &= isfinite(H->hx[c]) && isfinite(H->hy[c]);
        }
    if (fin) return 0;
bad:
    for (size_t k = 0; k < nc; k++) H->gx[k] = H->gy[k] = H->rx[k] = H->ry[k] = H->hx[k] = H->hy[k] = NAN;
    for (size_t k = 0; k < (size_t)nv1 * (ny + 1); k++) H->phi[k] = H->psi[k] = NAN;
    return -1;
}
static inline void mk_hhd_readout(const mk_hhd *H, double cx, double cy, mk_hhd_read *R)
{
    *R = (mk_hhd_read){{NAN, NAN}, NAN, NAN, {NAN, NAN}, {NAN, NAN, NAN}};
    if (!H->gx) return;
    const int nx = H->nx, ny = H->ny;
    if (!isfinite(cx)) cx = (nx - 1) / 2.0;
    if (!isfinite(cy)) cy = (ny - 1) / 2.0;
    double sd = 0, sc = 0, eg = 0, er = 0, eh = 0, ev = 0, mx = 0, my = 0;
    for (int b = 1; b < ny; b++)
        for (int a = 1; a < nx; a++) {
            const size_t ll = (size_t)(b - 1) * nx + (a - 1), lr = ll + 1, ul = ll + nx, ur = ul + 1;
            const double *gx = H->gx, *gy = H->gy, *rx = H->rx, *ry = H->ry;
            sd += 0.5 * ((gx[lr] - gx[ll]) + (gx[ur] - gx[ul])) + 0.5 * ((gy[ul] - gy[ll]) + (gy[ur] - gy[lr]));
            sc += 0.5 * ((ry[lr] - ry[ll]) + (ry[ur] - ry[ul])) - 0.5 * ((rx[ul] - rx[ll]) + (rx[ur] - rx[lr]));
        }
    const double ni = (double)(nx - 1) * (ny - 1), nc = (double)nx * ny;
#if defined(MK_MUTATE) && MK_MUTATE == 244
    R->expand = sd / ni, R->rot = sc / (2 * ni); /* planted fault: the divergence's mean not halved */
#else
    R->expand = sd / (2 * ni), R->rot = sc / (2 * ni);
#endif
    for (size_t k = 0; k < (size_t)nx * ny; k++) {
        const double vx = H->gx[k] + H->rx[k] + H->hx[k], vy = H->gy[k] + H->ry[k] + H->hy[k];
        eg += H->gx[k] * H->gx[k] + H->gy[k] * H->gy[k], er += H->rx[k] * H->rx[k] + H->ry[k] * H->ry[k];
        eh += H->hx[k] * H->hx[k] + H->hy[k] * H->hy[k], ev += vx * vx + vy * vy;
        mx += H->hx[k], my += H->hy[k];
    }
    R->trans[0] = mx / nc, R->trans[1] = my / nc;
    if (ev > 0) R->e[0] = eg / ev, R->e[1] = er / ev, R->e[2] = eh / ev;
    /* the focus: min over x0 of sum ((x - x0) v'y - (y - y0) v'x)^2, a 2 x 2 normal system */
    double s11 = 0, s12 = 0, s22 = 0, t1 = 0, t2 = 0;
#if defined(MK_MUTATE) && MK_MUTATE == 243
    const double w = 0; /* planted fault: the swirl left in */
#else
    const double w = R->rot;
#endif
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            const size_t k = (size_t)j * nx + i;
            const double x = i - cx, y = j - cy;
            const double vx = H->gx[k] + H->rx[k] + H->hx[k] + w * y, vy = H->gy[k] + H->ry[k] + H->hy[k] - w * x;
            const double m = x * vy - y * vx; /* (x - x0) vy - (y - y0) vx = m - x0 vy + y0 vx */
            s11 += vy * vy, s12 -= vx * vy, s22 += vx * vx, t1 += vy * m, t2 -= vx * m;
        }
    const double det = s11 * s22 - s12 * s12;
    /* a field left with no flow (a pure swirl) has no focus either */
    if (det > 1e-12 * s11 * s22 && s11 + s22 > 1e-24 * ev)
        R->foe[0] = cx + (t1 * s22 - t2 * s12) / det, R->foe[1] = cy + (s11 * t2 - s12 * t1) / det;
}
#endif
