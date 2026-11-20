# damn, readme is unreadable on line one. i didn't test for that
# i probably never will
# mathkit 1.6.2, the freeze with the second review's fixes and their review (the maintainer, 2026-11-13; 1.0.0 to 1.0.4 landed 2026-11-14, 1.0.5 to 1.3.0 on 2026-11-15, 1.4.0 and 1.4.1 on 2026-11-16, 1.5.0 and 1.5.2 on 2026-11-16, 1.5.3 on 2026-11-16, 1.5.4 on 2026-11-18, 1.5.5 and 1.5.6 on 2026-11-19, 1.6.0 and 1.6.1 on 2026-11-20, 1.6.2 staged 2026-11-20; 1.5.1 spent)

The house's CPU math kit: header-only C23, libc and libm only (and the caller's OpenMP runtime when it builds with
-fopenmp). It is the one place to take a general method from, beside cudakit for the GPU. Decided by the maintainer (the maintainer.md
#11) on the owner's request "consolidate or simplify what we have and the workflow. prevent future development that
re-solves the same problem"; stage 2 is the owner's "build all this: a CPU CG solver, quadrature, BFGS or CMA-ES, stiff
time steppers, spectral analysis, ray-tracing trees and point-to-triangle distance". SOURCES.md says which house copies
each piece merges and what was taken from each. No model uses the kit yet; owners move to it at their next release.

    #include "mathkit/include/mk_sparse.h"  # or -I mathkit/include
    make check                                          # every test, every planted fault: "wrong: 0"
    make check OMP=1                                    # the same with every test and planted fault built with -fopenmp
    make bench                                          # CG on a 64^3 Laplacian, BVH queries a second (OMP=1: threaded)

Build with our C23 (global/C23.md): `-std=c23 -ffp-contract=off`. The headers need no feature macros. `make check`
passes at -O2 and at -O3 -march=native -fno-math-errno -fno-trapping-math -funroll-loops, with -Wall -Wextra -Wpedantic
-Werror. LP64 only: the size checks take long as 64 bits, as on every machine here (x86-64 Linux).

Threads (1.5.0, mk_par.h): a caller that builds with -fopenmp gets OpenMP threads in the hot loops; one that does not gets no
threads and the same bits. Results never depend on the thread count: work is cut into blocks fixed by the item count, and sums
add the blocks' partials in block order. On a loaded machine set OMP_WAIT_POLICY=passive and OMP_NUM_THREADS to the free
threads: with libgomp's default spinning, 12 threads on the development machine at load 17 ran CG 30 times slower than one (MEASURED, below).

The freeze (1.6.0): the numerical core keeps its calls and results from here on; a fix that moves results is named in CHANGES.
The rows marked experimental (mk_dec, mk_sepvol, mk_graft, mk_ovl, mk_hhd, and mk_fftm in mk_fft.h) are not frozen: their results
may still change between releases. An underscore-suffixed name is internal to the kit (the tests reach some of them), never API.

| header | what | checks |
|---|---|---|
| mk_version.h | MK_VERSION ("1.6.2"), MK_VERSION_MAJOR, _MINOR, _PATCH; no other header depends on it | (t_vec) |
| mk_par.h | the kit's threads (1.5.0): MK_OMP_ (OpenMP pragmas, only under -fopenmp), mk_par_blk (blocks of 4096 items, at most 256, a function of the count alone), MK_PAR_FOR_, MK_PAR_SUM_ (partials added in block order). Included by mk_sparse.h, mk_geom.h, mk_sdf.h | t_par: 17 rows at 1, 2 and 12 threads + 8 (1.5.3: counts near LONG_MAX; 1.6.1: mk_welch from two threads on one plan); 17 BITS lines, a build without OpenMP's equal |
| mk_vec.h | mk_lerp (1.6.0: exact at both ends, no overflow; 1.6.2: a at t 0 beside an infinite b; the kit's one linear interpolation); mk_v3 (add, sub, scale, mad, dot, cross, len, dist, norm, lerp); mk_quat (mul, conj, norm, axis_angle, from_rotvec, rotate, unrotate, to_mat, from_mat, slerp, integrate); mk_xf rigid transforms | 15 (with the version row) |
| mk_fft.h | mk_fft_plan/free; mk_fft_c (double complex) and mk_fft_ri (split), in place, inverse unnormalized; mk_fft2_ri, mk_fft3_c; mk_rfft (real input, n/2+1 bins); any length (1.5.6, experimental): mk_fftm_plan/free, mk_fftm_c (mixed radix, Stockham: 4, 2, 3 by hand, other primes by their direct sum; the caller's work array) | 17 |
| mk_linalg.h | mk_lu, mk_lu_solve (0, or -1 for a non-finite answer), mk_lu_det (no overflow in the product); mk_chol, mk_chol_solve (as mk_lu_solve); mk_lstsq (Householder with column pivoting, columns scaled by their largest entry, rank-deficient safe, -1 for NaN or Inf in A, n <= 64); mk_sym_jacobi (symmetric decomposition, values ascending); mk_tridiag (Thomas, in metal's elimination order and bit-identical to it; refuses a zero or rounding-small pivot); mk_solve3, mk_det3, mk_sym3 | 36 |
| mk_rng.h | mk_rng (xoshiro256**, seed and stream; equals ribbon_rng.h bit for bit), _next, _u01, _u53, _below (unbiased), _gauss; mk_splitmix64, mk_mix64; mk_sum (Neumaier; a non-finite sum gives the plain IEEE result) | 10 |
| mk_opt.h | mk_nm (Nelder-Mead); mk_lm (Levenberg-Marquardt); mk_golden; mk_lbfgs (L-BFGS, strong-Wolfe line search, gradients given); mk_cmaes (CMA-ES, full covariance, optional per-coordinate scale) | 26 |
| mk_ode.h | mk_rk4_step; mk_dopri (Dormand-Prince 5(4), adaptive); mk_verlet_step (symplectic); mk_ode23s (Rosenbrock, adaptive, L-stable, Jacobian given or numerical); mk_bdf_step (BDF1/BDF2, variable step, Newton) | 22 |
| mk_stats.h | mk_select; mk_median; mk_quantile, mk_quantile_sorted (type 7); mk_welford (shifted by the first value), mk_welford_merge (Chan) | 16 |
| mk_interp.h | mk_pchip_slopes, mk_pchip_eval (bit-identical to aeromap.c's proved pchip but where its secant product underflows); mk_spline_slopes (natural); mk_hermite_eval; mk_linear_eval; the straight-line parts by mk_lerp (1.6.0); one-point tables give y0 | 10 |
| mk_sparse.h | mk_csr (from triplets, repeated positions summed in input order with any C library (1.6.2), or from ck_spd3's 3x3 block lower triangle), mk_csr_matvec; mk_cg with MK_PC_NONE, _JACOBI, _BJACOBI3, _IC0, and MK_PC_MEANFREE for a constant null space; mk_bicgstab; both in units of b's scale, -1 for a non-finite answer; 0 and 1 with the true residual, restarting from x while it is past tol (1.6.1) and falling (1.6.2). The interface follows ck_spd3_solve. Threaded: the product, dots, updates, Jacobi and block-Jacobi (IC(0) serial); CG fuses each product, update and Jacobi step with its dot | 31 |
| mk_quad.h | mk_gauss_legendre (n 1..2^30); mk_gk15 and mk_quad_adapt (adaptive Gauss-Kronrod 7-15 with an error estimate, maxint under 2^29); mk_tri_rule (degree 1, 2, 5), mk_tet_rule (degree 1, 2, 3) | 13 |
| mk_spec.h | mk_window (rectangular, Hann, Hamming, Blackman-Harris, Kaiser), mk_bessel_i0; mk_spec_init (the plan only read, 1.6.1: one plan serves many threads); mk_stft_frame; mk_welch (one-sided density or power); mk_csd, mk_coherence, mk_transfer | 8 + 1 (t_par) |
| mk_geom.h | mk_closest_pt_tri (Ericson; the edges as well on a sliver or a repeated corner); mk_ray_make, mk_ray_tri (watertight, Kahan's a b - c d); mk_aabb; mk_bvh_build (binned SAH to depth 200, then median), _first_hit, _any_hit, _count_hits, _nearest, _nearest_lab (a label left out), _box (triangles whose boxes meet a box), _inside (parity on the rays that miss every edge, else the winding number); mk_winding_number; mk_bvh_nearest_many, mk_bvh_inside_many, mk_winding_number_many (threaded batches, each answer the single call's); internal (1.5.3): mk_bvh_count_wind_, the hits' facings summed along a ray (the winding number on a closed mesh) | 24 + 3 (t_par) + 1 (t_regr) |
| mk_vpflow.h | mk_vp (a vector potential in cubic B-splines on a grid, any axis periodic), mk_vp_ok, mk_vp_eval, mk_vp_vel (u = curl A, divergence-free everywhere), mk_vp_div; mk_vp_advect (the volume-preserving map, mk_rk4_step's RK4), mk_vp_jacdet. b-helper's | 12 |
| mk_hhd.h | (experimental) the Helmholtz-Hodge split of a 2-D motion field (1.5.6, the owner's ask for optic flow and motion vectors: the maintainer's frame generation, retina's and neuro's vision, fly's motion pathway): mk_hhd_init/free, mk_hhd_split (gradient part grad phi, rotational part rot psi = (dpsi/dy, -dpsi/dx), harmonic rest; phi = psi = 0 on the frame, so the three are orthogonal and a uniform shift is harmonic; box differences at pixel corners; a sine transform by mk_fftm and mk_tridiag; 1920 x 1080 in 0.14 s on one thread, 0.03 s on 8), mk_hhd_readout (expansion and rotation rates, mean shift, focus of expansion, each part's share of the energy) | 33 + 2 (t_par) |
| mk_snet.h | mk_snet_build (labels) and mk_snet_build_iso (a scalar field), one vertex per dual cell where the labels differ, shared by every label there; Taubin relaxation held inside each dual cell; mk_snet_faces, mk_snet_faces_iso (anticlockwise from outside); mk_snet_vertex; mk_snet_split (where a surface touches itself); mk_snet_free. Memory with the surface, not the grid. s-forms' | 21 + 7 meshcheck |
| mk_sdf.h | mk_sdf_mesh (a triangle mesh to a signed distance grid, negative inside: Batty's band and sweeps, sign by row parity with the winding number for odd rows and near nodes; MK_SDF_EXACT, MK_SDF_FILL, MK_SDF_LATTICE), mk_sdf_sample (trilinear), mk_sdf_free. s-bone's. Threaded, 1.4.0's grid bit for bit: band by tiles, sweeps as a wavefront of tiles, EXACT, the winding numbers and near nodes node by node | 17 + 2 (t_regr) |
| mk_sepvol.h | (experimental) mk_sepvol: parts of closed meshes sunk into each other moved apart, vertices, order and faces kept; one mesh's parts (both give way) or a part out of a held one (unchanged byte for byte); toward the nearest point of the other surface (edge-aware), crossing and near face pairs by selfsep's rule, never outward; mk_sepvol_report (inside, crossing pairs, rounds, moves, volume before and after, the status bits: 1.5.5 says when parts still cross or a run stopped short of crushing one). cad-helper's, with the maintainer's selfsep | 73 + 11 meshcheck |
| mk_dec.h | (experimental) mk_decimate: quadric half-edge decimation keeping the input's vertices, refusing collapses that add folds, crossings or drift past a two-sided bound or break edge, angle or turn limits; held and seam vertices, open rims, a crossing and needle clean-up, fold sums by region; mk_dec_defaults (mesher's values for a unit), mk_dec_report. mesher-fable's folddec | 15 + 3 meshcheck + 3 (t_regr) |
| mk_ovl.h | (experimental) mk_overlap_volume, mk_ovl_bvh: the volume two closed meshes share by columns along an axis, and (MK_OVL_ALL) each one's volume and what each has that the other has not; odd lines left out and counted. mesher-fable's ovl with stepcad's voldiff | 9 + 9 (t_regr) |
| mk_graft.h | (experimental) mk_mtet (marching tetrahedra on mk_sdf's grid layout: closed and outward whatever the field); mk_graft_init, mk_graft_patch (a field's surface kept on one side of a seam function), mk_graft_stitch (a mesh cut at the seam with the patch welded in, the rest triangle for triangle the mesh's own); mk_gmesh, mk_graft_report. mesher-fable's bonegraft | 15 + 3 meshcheck + 5 (t_regr) |

What the checks test:
- vec: a known 90-degree turn, composition, matrix round trip, inverse transform, slerp, constant-rate integration.
- fft: a cosine's two bins, Parseval, round trip, a direct DFT in 1-D, 2-D and 3-D, split equal to complex, real input; (1.5.6)
  mk_fftm against a direct DFT both ways at 21 lengths from 1 to 2160 (radices 4, 2, 3, 5, 7 and the prime 97), a round trip
  at 3840, equal to mk_fft_c at 1024, its refusals.
- hhd (1.5.6, t_hhd): a pinhole camera's flow by Longuet-Higgins and Prazdny 1980 (MADE UP numbers: 256 x 144, f 1000 px,
  2 m/s toward a wall 4 m off, heading (50, -25) px, roll 0.8 rad/s, yaw 0.02 rad/s): the focus at the heading to 1e-6 px,
  the expansion rate W / Z and the roll as the rotational part's rate to 1e-12, the heading's sideways flow as the harmonic
  part's mean; with yaw, the rates unchanged, the pan's -f B (1 + <x^2> / f^2) in the harmonic mean and the focus moved by
  f B Z / W (one plane cannot tell a pan from a heading); each motion alone (roll 62 % rotational and 38 % harmonic on 16:9;
  forward the same split into gradient and harmonic; yaw 99.997 %, 99.71 % and 96.73 % harmonic at 15, 46 and 81 degrees across;
  sliding sideways all harmonic); a tilted plane (the expansion rate and the slant's own shear in the rotation read to 1e-12
  against their formulas, the focus within 0.24 px); a box room at 65 degrees (the focus within 0.12 px, the rotation 0.4 %
  over the roll); the parts orthogonal, adding up, the rotational part divergence-free and the gradient part curl-free to
  rounding, the harmonic part both to 1e-10; fields made of known potentials and a shift plus strain split back to 1e-12 at
  37 x 23 (the direct-sum radix 23) and 256 x 144; refusals; 1920 x 1080 orthogonal (its time printed, not checked since 1.6.1). t_par: the same bytes at
  1, 2 and 12 threads.
- linalg: a zero first pivot, the determinant, singular and indefinite refusal, residuals, least squares with unequal
  and with repeated columns, a 6 x 6 symmetric matrix with known values, 3x3 paths; the tridiagonal solve on a known
  answer, backward error on 20 random systems of 1000, bit-equality with metal's world.c, and its refusals.
- rng: the authors' reference values, ribbon equality, Gaussian moments, uniform range, chi-squared, 1 + 1e100 + 1 - 1e100.
- opt: Rosenbrock (Nelder-Mead, Levenberg-Marquardt, L-BFGS in 2-D and 10-D, CMA-ES), a 4-D coupled quadratic, an
  exponential fit, golden section, a rotated 10-D ellipsoid of condition 1e6 (CMA-ES), 3-D Rastrigin (CMA-ES with
  doubling populations reaches the global minimum; one small-step run is shown stopping in a local basin).
- ode: RK4's order 4, Dormand-Prince over 10 periods and against tolerance, backward integration, Verlet's order 2 and
  bounded energy; Robertson's kinetics to t 40 against Hairer and Wanner's reference (ode23s and BDF2), Van der Pol at
  mu 1000 against its asymptotic period, BDF1 and BDF2 orders 1 and 2, and step counts against Dormand-Prince; the
  last step onto t1, a NaN first step under hmax, NaN times and a negative n refused, maxsteps LONG_MAX (the verification work O-E').
- stats: known medians, 1200 random arrays with ties, quantiles of 1..10, Welford, stability at an offset of 1e9.
- interp: 22000 queries bit-identical to aeromap.c, monotonicity, a natural spline on sin, straight lines; one-point
  tables, secants of 1e-170, the spline's work size at n 2^30 and locate in a table of 2^30 + 4.
- sparse: a 30 x 30 Laplacian against a known solution, iterations within CG's sqrt(kappa) bound, one mode solved in one
  step, IC(0) cutting iterations about 3x, block-Jacobi on ck_spd3's layout against a dense LU, a singular periodic
  Laplacian, BiCGSTAB on upwind convection-diffusion, and CG refusing to solve a nonsymmetric system; b at 1e-170 and
  1e170, b times 2^600 giving x times 2^600 bit for bit, block-Jacobi at 1e-101 and 1e110 and on singular blocks,
  tol 0 at three scales of b, A near 1e290 with Jacobi, an answer past the largest double, and the size refusals.
- quad: Gauss-Legendre exact to degree 2n - 1 and not past it (n 1..20), tabulated n 5, Kronrod 15 exact to degree 22,
  its error estimate, 1/sqrt(x) and log x with singular ends, an oscillatory integral, triangle and tetrahedron rules
  exact to their degree and not past it.
- spec: I0 against its tabulated values, Parseval for every window, a sine's amplitude from every window, white-noise
  level (2 s^2 / fs), coherence and transfer function of a known noisy linear system, Blackman-Harris's -92 dB
  sidelobe, STFT zero padding.
- geom: closest point against an independent method (200000 cases with slivers), ray-triangle against Moller-Trumbore,
  rays aimed exactly at every edge and vertex of a closed sphere (none slips through; Moller-Trumbore misses 1598 of
  15872), the tree's first hit, any-hit and nearest point equal to brute force on a sphere and a torus, inside or out
  by parity and by winding number against the analytic solids; the box query against every triangle's box; the hits' facings
  summed along 3000 rays against the winding number of two spheres sunk into each other (0, 1 and 2; 1.5.3); and a row
  for each of the verification work's findings on 1.2.2 (F1-F7, F8a, F9, F10, F12: the open interval at tmax, NaN and 1e308 vertices, a tree 873
  deep, a cube's centre with rays through its edges, a tie in the nearest distance, counts past INT_MAX / 2, a tiny
  triangle off the ray, slivers, a repeated corner, a 5e10-wide triangle, a thin triangle's own box).
- sepvol: twoball (single mode), ballB out of a held ballA (unchanged byte for byte, its face centres at least gap/4
  off), two cubes 1 mm into each other, each to no vertex inside and no crossing pair, judged again by each part's
  winding number; no move outward along a vertex's normal; the volume change (the cubes lose 1.2 times their overlap);
  in one round every inside vertex moves the full step, and the cubes' edge vertices go along -x, not their 45-degree
  normal; a third part beside the crease untouched and entered by nothing; the refusals; (1.5.3) each of a call's allocations
  failing in turn, V as it came in each time (the verification work F-1), a part thinner than eps inside nothing (F-A), the facings' decision
  equal to the winding number's at every vertex of three closed meshes, and an open sheet left to the winding number; (1.5.5) three
  spheres and nested shells read by coverage against the other parts' winding number (the verification work F-C), on a closed mesh and on one that
  is not, a held part wound inward refused and a hollow held solid taken, two crossed bars stopped before they are crushed and
  reported as still crossing (F-B); then the results through meshcheck, the cubes before separation reading FAULT.
- regr (1.5.5, t_regr): the verification work's findings against 1.5.3, each its exact call and failing on 1.5.3's headers: mk_ovl O-1, O-1b (hangs, run in a
  child with a limit of 5 s of CPU), O-4, O-5, O-5b, O-6, O-7, O-7b, H-end; mk_dec D-2 (the same limit), D-5, D-6; mk_graft G-1, G-2,
  G-3, G-5, R-1; mk_geom G1; mk_sdf S-1, S-2.
- dec: a unit sphere and a torus decimated within a bound (every input vertex near the result and the result near the input,
  by brute force), closed with Euler unchanged, the kept vertices the input's own; a target; held vertices; regions; a 12 x 3
  torus kept a torus by the link condition alone; a hollow ball whose walls are closer than the edges grow (no crossing pair,
  by brute force); a bumpy sphere with no fold made; an open cap whose rim stays the input's rim; refusals; the sphere, torus
  and hollow ball through meshcheck.
- ovl: two spheres' lens against its analytic volume along each axis, a sphere's own volume, a box pair on the lattice
  (exact), disjoint and nested pairs, only_a + both = vol_a, an open mesh's odd lines counted, refusals.
- graft: marching tetrahedra on a sphere's and a torus's fields (closed, Euler, volume), then a dent carved into a unit
  icosphere within a 45-degree cap: the field from mk_sdf_mesh, the patch thinned by mk_decimate, stitched back closed in one
  piece, the volume against the sphere less the ball's lens, S's own vertices and its triangles away from the seam kept as
  they were, nothing inside the ball; refusals; the results through meshcheck.
- vpflow: div u at 200 points, a ball's volume through the flow (1.8e-6 of it), the map's Jacobian (1 within 2e-7), RK4
  against half the step, a quarter turn of a rigid rotation against its exact answer at 5 to 40 steps (order 4.00), the
  rotation's field exact, a periodic axis and a uniform flow across its seam, the refusals; the verification work V-A (INT_MAX - 3 nodes) and V-B (a
  refused grid gives NaN from eval, vel and div).
- snet: a ball raw and relaxed (closed, one winding, Euler 2, each vertex in its dual cell, volume), a ball from its
  distance field (radial error 0.045 cells against 0.47 from labels), a ball split into two labels sharing their 610
  interface triangles reversed, a torus (Euler 0), two voxels meeting at an edge (split cells reported), the refusals;
  then each surface through mesher's meshcheck as a tool (closed, manifold, no flips, Euler and genus, winding,
  outward volume, no crossings), the edge-touching pair reading FAULT as it must; the verification work
  N-A (an empty net's vertex is -1).
- sdf: the unit sphere and a torus against their analytic distance (within the faces' sag near the surface, within 0.1 h
  beyond with the sweeps, within the sag everywhere with MK_SDF_EXACT), the sign at every node of three nested shells and
  the winding number's agreement, a random closed mesh against brute force (distance and sign at every node), a sphere
  with a face missing, nodes 2.1e-6 inside a prism's slanted faces, a region covered twice with and without MK_SDF_FILL,
  MK_SDF_LATTICE's shared nodes, the refusals (an axis or a node count that overflows), trilinear samples; the verification work S-A to S-D
  (a lattice past the double range, sampling an empty grid, a re-cast ray that crosses nothing, refusals that leave the
  grid empty).
- par (1.5.0, t_par): built with -fopenmp, each threaded routine at 1, 2 and 12 threads gives the same bytes (17 rows: the
  product, CG with each preconditioner and mean-free, BiCGSTAB with none and Jacobi, mk_sdf_mesh plain, exact, filled on the
  lattice, on an open sphere (73 odd rows) and on a box with 517 near nodes, the three batch queries, mk_hhd_split), every output cleared
  before each run; the batches against their single calls bit for bit; the sdf cases reaching what they test; CG's block sums
  against one serial sum (1.4.0's, written out in the test); the blocks of seven counts from 2^62 + 1 to LONG_MAX (1.5.3, the verification work
  H18). check.sh then requires the BITS lines (a hash of every result) of a build without OpenMP to equal the OpenMP build's.

Planted faults: test-only `MK_MUTATE` numbers in the headers, compiled out unless defined (as in mcphys-kc's linalg.c).
check.sh builds each one, and its test must print a FAIL row (1.5.5, test/judge.sh). A timeout is never a catch, nor is a crash or a bare
nonzero exit; a fault meant to be caught by a crash would be named in check.sh's CRASHES with its signal, and none is. The real header
passes only with at least one PASS row, no FAIL row and exit 0. A failed build or run counts as WRONG, no output of an earlier run is read,
and the exit status cannot wrap to 0. test/drvctl.sh, which check.sh runs, shows each of these rules on a known-bad input. There are 140, all
caught:
- vec 1-2: the rotation's sign, a quaternion product term.
- fft 11-12: the inverse's conjugate, the bit reversal; 238, 239 (1.5.6): mk_fftm's quarter turn always forward's, its radix-3 difference turned.
- hhd 240-245 (1.5.6): rot psi's second component unturned, the tridiagonal's coupling sign turned, the inverse transform's 2 / ny left
  at 1 / ny, the swirl left in the focus fit, the divergence's mean not halved, the curl's sign turned.
- linalg 21-27: LU pivoting, Cholesky's square root, least-squares rescaling, the Jacobi sign, 3x3 pivoting, the
  tridiagonal pivot left uneliminated, its small-pivot test reduced to exact zero.
- rng 31-34: xoshiro's rotation, Box-Muller's 2, the sum's compensation, the non-finite guard.
- opt 41-45: the centroid, LM's gradient sign, golden's comparison, L-BFGS's two-loop sign, CMA-ES without covariance
  learning.
- ode 51-55: RK4's weights, Dormand-Prince's error estimate, Verlet as explicit Euler, ode23s's second stage, BDF2's
  coefficient.
- stats 61-65: the even median, the quantile's index, Welford's n - 1, the shift, the merge's d^2 term.
- interp 71-72: pchip's harmonic slope, the natural end.
- sparse 101-106: CG as steepest descent, BiCGSTAB's omega step left out of x, IC(0) reduced to a diagonal, the
  finite-answer check removed, the solvers' power-of-two units removed, the block-Jacobi test made absolute again.
- quad 111-114: a Gauss-Legendre weight, a Gauss 7 weight, a triangle point, the tetrahedron constants.
- spec 81-84: the one-sided doubling, I0 truncated, a Blackman-Harris coefficient, the cross-spectrum's conjugate.
- geom 91-95: an Ericson region test, front faces only, the tree's pruning, the winding number's sign, the box query
  losing a box that only touches.
- vpflow 131-132: a velocity that is not a curl, forward Euler for RK4.
- snet 121-125: the dual-cell hold left out, the winding flipped, one axis's faces left out, the iso crossing put at the
  midpoint, the split test reading 0.
- sdf 141-145: no sweeps, odd rows kept on parity, the squared distance in the band, the flood left out of MK_SDF_FILL,
  the near nodes kept on parity.
- sepvol 151-155: sepvol's push along the vertex's own normal, the vertex's own patch not left out of its search, one
  round only, no fallback for a crossing with no corner over the other face, near faces left alone.
- sepvol 156-159, 191-192 (1.4.1): no crossing guard (every move whole), every held mesh taken, faces of two parts guarded in single
  mode, no smoothing, two parts crossing counted as a self-crossing, the last state returned instead of the best.
- sepvol 193-195 (1.5.2): the nearest exit always (no held exit), the held exit without the check on M's own surface,
  MK_SEPVOL_PAIRS read as on whatever it says.
- dec 201, graft 202-205, 209, 210, geom 206, sdf 207-208 (1.5.2): as its CHANGES entry says.
- sepvol 212, 213, 215 (1.5.3): V0 filled after every allocation is checked (F-1), even parity read as covered twice (F-A),
  every mesh taken as closed; geom 214: the hits' facings turned; par 211: 1.5.2's ceilings.
- sepvol 216-221 (1.5.5): odd parity read as covered once on a closed mesh, and on an open one, only the held mesh's total volume looked at,
  every part wound inward refused (a cavity too), no stop before a part is crushed, the status silent about pairs left.
- regr 222-237 (1.5.5), each undoing one patch of 1.5.4: 222, 223 mk_ovl's double lattice counters on each axis, 224 the lines' origin below
  both meshes, 225 h * h first, 226 the hit search from t = 0 only; 227 mk_dec's large triangle in every cell of its box, 228 a small bound
  squared to 0, 229 maxown counted before the crossing test; 230 mk_mtet's lo summed, 231 no check of the far corner, 232 no bisection of
  an infinite psi, 233 the map's pointers left set, 234 mk_graft_refine's coordinates unchecked; 235 mk_aabb_ray's pad without the ray's
  origin; 236 a refused box leaving axis 0's fields, 237 a cell past an int lost.
- dec 161-165: the fold refusal, the distance refusals, the crossing test, the link condition, the keep flags left out.
- the review's rows (1.6.0), one for each confirmed fault fixed, each shown failing on 1.5.6's headers and passing on 1.6.0's, no planted fault: dec D1,
  stats S1, interp I1, ode O1, opt C1, sparse B1, sdf N1, spec N2, vec Q1 (CHANGES 1.6.0).
- the second review's rows (1.6.1), the same way on 1.6.0's headers: ode #1 and #8, opt #2, par #3, vec #4, stats #5, sparse #10 and #18,
  geom #13 and #19, snet #14 (CHANGES 1.6.1).
- the rows from the reviews of 1.6.1 (1.6.2), the same way on 1.6.1's headers: ode #1 with two components and A1, opt #2's later ties and
  A5, sparse A2, N4 and #15, vec A3 and A4; drvctl's control for J1 (CHANGES 1.6.2).
- ovl 171-173: the union for the overlap, lines on the lattice itself, an odd count used as it is.
- graft 181-182: marching tetrahedra wound inward, the patch's fans left out of the stitch.
- threads (1.5.0; built with -fopenmp, caught by t_par's rows at 1, 2 and 12 threads): 107 the block partials added by
  OpenMP's reduction(+:), 108 the product's rows cut by the thread count with each piece's last row left out, 146 the sweeps'
  planes in chunks by the thread count with the last chunk first, 147 the band's tiles cut by the thread count with a triangle
  in its first tile only, 148 the exact pass cut as 108, 96 each winding number summed by reduction(+:), 97 and 98 the nearest
  and inside batches cut as 108. Between them every thread row is shown failing; the rows no header fault reaches (the
  version, the block-sum row, the sdf cases' reach, the BITS comparison, the inside batch against single calls without
  OpenMP) fail on their own known-bad inputs in test/negctl.sh, which check.sh runs.

Speed, `make bench`, one thread on the development machine at load about 15-21 (2026-11-13 22:5x; MEASURED, noisy):
- CG, 64^3 Laplacian (262144 unknowns), tol 1e-8: none 230 iterations, 0.39 s; Jacobi 230, 0.38 s (constant diagonal);
  IC(0) 74 iterations, 0.35 s (4.8 ms an iteration against 1.7 ms).
- BVH on a 262144-triangle torus: built in 0.46 s. First hit 1.36 M rays a second against 1166 by brute force (1168
  times); nearest point 52900 a second against 763 (69 times; the points fill a box around the torus, many far from
  it, which is the tree's worst case).

Threads, MEASURED 2026-11-16 (1.5.0; scratch/mathkit-150-work/cmp/timing.sh: the kit call's own time, median of 5 runs under
`jobwatch run -n 5`, OMP_WAIT_POLICY=passive; hosts loaded by other sessions throughout, so these are what a loaded house sees):

| case (real meshes from the consumers) | the development machine, 1 thread | the development machine, 6 threads | speedup | a second machine, 1 thread | a second machine, 16 threads |
|---|---|---|---|---|---|
| mk_sdf_mesh femur b67_femur_l, 0.5 mm (54.9 M nodes; sweeps 98 %) | 30.4 s | 14.3 s | 2.1 | 66.5 s | 155.6 s |
| mk_sdf_mesh external oblique r, 1 mm (3.29 M nodes) | 4.56 s | 1.83 s | 2.5 | 8.13 s | 20.1 s |
| mk_sdf_mesh body rest-plain, 5 mm (7.08 M nodes, one odd row) | 28.5 s | 6.81 s | 4.2 | 93.9 s | not reached |
| mk_sdf_mesh body, 8 mm, MK_SDF_EXACT (1.76 M nodes) | 25.7 s | 12.2 s | 2.1 | 74.5 s | not reached |
| mk_cg 64^3 Laplacian, none / Jacobi / IC(0) | 0.46 / 0.49 / 0.47 s | 0.42 / 0.45 / 0.50 s | 1.1 / 1.1 / 0.95 | 1.4 / 14.2 / 10.7 s | not reached |
| mk_bicgstab 48^3 convection-diffusion, none | 0.076 s | 0.088 s | 0.87 | 1.6 s | not reached |
| mk_bvh_nearest_many, 200000 points, 262144-triangle torus | 3.41 s | 0.51 s | 6.7 | 76.6 s | not reached |
| mk_bvh_inside_many, 50000 points | 0.227 s | 0.050 s | 4.5 | 4.30 s | not reached |
| mk_winding_number_many, 1000 points | 3.83 s | 1.01 s | 3.8 | 57.5 s | not reached |

the development machine was at load 12-22 on its 12 threads (6 used, the desk's limit); a second machine at load 78-84 on 16, where one thread got about half a
CPU and 16 more threads ran 2.3-2.5 times slower than one; that run hit its 60-minute limit before the rest. The sdf cases
and the batches gain where a thread has a large piece of work; CG and BiCGSTAB do not (memory-bound passes of a fraction
of a millisecond; IC(0)'s serial solves). With libgomp's default (spinning) at every thread of a loaded machine it is far
worse: 12 threads on the development machine at load 17 ran CG on 64^3 in 9.6 s against 0.31 s on one, and the femur at 1 mm in 21 s against
6.2 s (before the sweeps were tiled); passive, the same CG took 0.55 s. Size OMP_NUM_THREADS to the free threads.
A caller without OpenMP loses nothing (the development machine at load 11-14, one thread, 1.5.0 against 1.4.0; scratch/mathkit-150-work/cmp/serial.sh, alt.sh):
mk_sdf_mesh
4.25 s against 4.99 (external oblique 1 mm) and 22.3 s against 25.1 (body 5 mm); CG none 4-14 % faster in seven alternated
rounds (the fused passes), IC(0) within 5 % either way.

Not yet:
- a sparse direct solver, multigrid on the CPU (cudakit has AmgX on the GPU), ILU for BiCGSTAB;
- a BDF driver with error control (mk_bdf_step is the step: the caller sets h);
- a fast winding number through the tree (Barill et al. 2018; mk_winding_number is O(triangles)), for open meshes (1.5.3: on a
  closed mesh mk_sepvol sums the ray hits' facings instead);
- (done in 1.5.2: mk_sdf_mesh_box)
- quadrature on the hexahedron beyond the tensor Gauss-Legendre product;
- proofs: the verification work is to fuzz and prove the parts whose failure costs most (dense solves, FFT, LM), as its queue allows.
- threads in IC(0)'s triangular solves (level scheduling would keep the bits, the backward solve as a gather in falling row
  order, but a 64^3 grid has 190 levels: about 380 barriers an iteration), in mk_snet, mk_sepvol, mk_ovl and mk_graft (not
  hot in any profiled run), and mk_decimate's greedy loop (serial by design);
- a batch of first hits (mk_bvh_first_hit_many) for partjudge's and bonegraft's edge rays, when one of them asks;

Known gaps in the checks: Nelder-Mead's shrink step has no test of its own (mk_bdf_step's Newton-failure path has one since
1.6.1, the second review's #1); no row reaches the NaN guard of mk_ode23s's and mk_bdf_step's norms with a finite value after the NaN
(mk_lu_solve spreads a NaN to every component first: 1.6.2's two-component row passes there on 1.6.1 too); mk_dec D-3 (int marks past about 6e7 collapses) is patched and has no row: reaching it takes that many collapses;
its unsigned marks wrap past 2^32 stamps, about 9e7 collapses at 1.6.0's 48 a collapse (the review's D2, documented, no row).

Known limits (1.5.5), none patched:
- mk_sepvol: parts crossing with no vertex inside each other cannot be cleared by moving vertices: the run keeps them whole and says so
  (MK_SEPVOL_STILL_CROSSES, MK_SEPVOL_CRUSHED); a part crossing itself is not undone; held cubes near 2^540 are refused with -2 FLIPPED (the
  volume overflows; held_vol NaN).
- mk_ovl O-2: meshes spanning about 1e308 can still lose crossings (mk_ray_tri's T overflows).
- mk_dec D-4: time quadratic in a hub vertex's valence.
- mk_graft G-6: a vertex pinched between two fans is not looked for; G-7: an ear can leave a mirrored pair; W-1: a rim vertex welded to an S
  vertex can move past snap by up to 1 % of an S rim edge (documented).
- mk_sdf S-3: a box one node thick gets no sweeps (far nodes keep BIG, still an upper bound); N-1: above about 2^39 h a face-diagonal ray is
  counted twice.
- (1.6.1, the second review's, documented) mk_dopri's and mk_ode23s's step floor is 1e-14 of the span, so a fast start on a long span needs
  a call of its own (Robertson's from 0 to 1e11 in one call is -1); mk_lu refuses a pivot under 1e-13 of the matrix's largest entry, which
  turns away ode23s's steps past about 1e10 on Robertson's (the answers right); mk_solve3 refuses only an exact zero pivot (a matrix
  singular to rounding gives 0 and a huge x); mk_gauss_legendre is O(n^2) (n 65536 in 23 s); an infinite value makes mk_welford's mean NaN
  from the next value on; mk_tri_rule past degree 5 and mk_tet_rule past 3 give 0 points.
- (1.6.2, the reviews of 1.6.1, documented) hmax is met to within half an ulp of t (at t 1e17, hmax 9 steps 16); mk_cg and mk_bicgstab
  end with 1 when a restart does not lower b - A x, so a tol at the attainable accuracy can end 1 where more restarts might have dipped
  under it (9 of 736 sweep cases, true residual 1.1 to 3.5 times tol), and on -1 relres is the iteration's last; test/judge.sh reads an
  exit of 124 or 137 as a timeout even when it is the FAIL count (the largest test has 73 rows).

Stage 3 fixes (2026-11-13, MEASURED by the maintainer's cudakit cross-check):
- mk_sum turned an infinite addend into NaN. Compensation now stops once the sum is not finite, so it gives the plain
  IEEE result: inf, -inf, or NaN for inf + -inf.
- mk_welford accumulates data shifted by its first value, and mk_welford_merge adds Chan's pairwise merge, cudakit's.
  On 1e6 values near a mean of 1e6, the variance is 3.5e-14 off a long-double two-pass, against 3.6e-11 unshifted;
  near 1e9 it is 4.7e-14, against 1.3e-8. Seven uneven merged pieces give 1.8e-15.

## CHANGES

### 1.6.2 (installed 2026-11-20 10:49)

Two outside reviews of 1.6.1, read without compiling: scratch/the maintainer/api-eval/reports/mathkit-r2-opus-5-5.md (an Opus model:
#1 partial, N1 to N5) and scratch/the maintainer/api-eval/reports/mathkit-r2-fable-5-2.md (a Fable model: #1 and #2 partial, J1, A1 to
A7). Each finding was run on 1.6.1 first with a short program: the claims table scratch/mathkit-162-work/claims/claims.tsv (the
programs and outputs beside it, in claims/ and claims/out/): 13 CONFIRMED (N3 by reading, with no measured effect; A6 / N5 on a made-up
test, out of the kit's reach), and skipped as settled in 1.6.0 or 1.6.1: #9 (mk_solve3's and mk_lu's pivot rules, documented), #17
(mk_welford and an infinite value, documented), P1 (more than INT_MAX distinct triplet positions, documented). Staged from 1.6.1 in
scratch/mathkit-162; work in scratch/mathkit-162-work (STATE.md). The experimental marks stay as 1.6.0 set them.
- #1 again (mk_dopri, mk_ode23s, mk_bdf_step): 1.6.1's max, if (!(q <= err)) err = q, lost a NaN that a finite q followed, so only the last
  component's NaN was kept: dopri on (sqrt(1 - t), 1) over [0, 2] returned 0 with y[0] NaN. Now if (q > err || q != q) err = q in all
  three: -1 at y (2/3, 1). In ode23s and bdf_step the fault was masked (mk_lu_solve spreads a NaN to every component), so those parts of
  the row pass on 1.6.1 too. The same bits wherever no NaN arises.
- A1 (mk_dopri, mk_ode23s): 1.6.1's snapped step (t + s) - t could round a shrunk step back up to the one just rejected (a step of 1 to 4
  ulps of t), and that step was retried to maxsteps: y' = -0.05 y from t 1e17 over 992 at rtol 1e-2 to 1e-4 took 1e7 rejections and 6e7
  calls of f, then -1; ode23s's path for a refused LU (h halved, t + 8 tying up to t + 16) the same. A step after a rejection is now below
  the one rejected (nextafter(t + h_rejected, t) - t): rtol 1e-2 and 1e-3 finish (0), tighter ones are -1 after 3 rejections, the LU path
  after 1. Ordinary runs bit for bit (every t_ode number as 1.6.1). Chosen over the review's return -1, which would stop runs that finish.
- A2 / N1 (mk_cg, mk_bicgstab): with tol under the attainable accuracy, 1.6.1 restarted after every pass to maxit (t_sparse F9's 144
  unknowns at tol 1e-15: 1e6 iterations, 1.5 s, then 1). A restart whose b - A x is not below the pass before's now ends with 1 and that
  true relres: 42 to 81 iterations there. In a sweep of 736 cases (1-D and 2-D Laplacians, x0 0 to 1e6, tol 1e-8 to 1e-16, claims
  c_a2b.c) 513 give the same rc and iterations, 238 that ran to maxit stop early (3.7e6 iterations fewer), and 9 that 1.6.1 brought to 0
  by a random walk at the floor return 1 (the true residual 1.1 to 3.5 times tol). Restarts that work run as before (1-D 5000 from x0
  1e6: 9678 and 20458 iterations); t_par's BITS lines as 1.6.1's.
- N4 (mk_cg, mk_bicgstab): precond bits past MK_PC_MEANFREE are refused, -1 before iterating (1.6.1 ran 32 as MK_PC_NONE, 33 as Jacobi, and
  -16 as NONE | MEANFREE, which gave 0 with x 0 on a nonsingular A). N3: the comment says on -1 relres is the iteration's last (it equalled
  b - A x to 6 digits in the 4 breakdowns tried).
- #15 (mk_csr_from_triplets), decided before the freeze: repeated positions are summed in input order with any C library. The sort's ties
  go by input place (mk_trip_ keeps the triplet's index in a union with its value: still 16 bytes). 1.6.1 took qsort's order of equal keys,
  and glibc 2.44 keeps input order only while it can malloc its merge buffer: with the address space capped so that it cannot, 1.6.1
  summed 3636 of 10000 positions otherwise (1.6.2: 0). Against 1.6.1 on glibc's usual path, 224 random sets (nnz 1 to 2e6, n 3 to
  500001, values over 120 binades with signed zeros) are byte for byte in rowptr, col and val; 2e7 triplets 3.10 -> 3.12 s.
- #2 again (mk_nm): 1.6.1 skipped only a tie at the start, and one reached later still stopped the search (x^2 from -3.5, step 1: 0.25 at
  {-0.5, 0.5} after 6 calls; 15 of the 20 half-integer starts -9.5 to 9.5; 3 and 4 of 20000 random 1-D starts). The spread must now be
  under tol on two iterations running (1.6.1's it && cut): 0, 0 of 20, 0 of 20000. One iteration more per converged run (about 2.7 calls),
  so results move: t_opt's Rosenbrock from (-1.2, 1) f 1.2e-20 -> 6.6e-22 (235 -> 243 calls), its 4-D quadratic 2.8e-25 -> 1.6e-25, the
  tied starts of 1.6.1's row 32 -> 34 and 65 -> 67 calls.
- A3 (mk_q_axis_angle): a zero axis gives the identity only for a finite angle; NaN or inf gives NaN, as about any other axis (1.6.1 the
  identity).
- A4 / N2 (mk_lerp): exactly a at t 0, also beside an infinite b (1.6.1: NaN, so mk_linear_eval at a node next to an inf was NaN), and -0
  kept (1.6.1: +0). Every t_interp and t_stats number as 1.6.1.
- A5 (mk_cmaes, and mk_lbfgs the same): a refusal zeroes the counts (1.6.1 left evals, and lbfgs's iters_used, from the run before).
- J1 (test/drvctl.sh): a control for 1.6.0's order in judge.sh: 124 FAIL rows, then a hang, must read WRONG (timeout); on the order before
  1.6.0 it reads caught (scratch/mathkit-162-work/neg/j1/old.txt). A6 / N5: judge.sh's rule text says an exit of 124 or 137 is read as a
  timeout even when it is the FAIL count (a loud WRONG; the largest test has 73 rows).
- Documented, not changed: A7 (hmax is met to within half an ulp of t: at t 1e17, hmax 9 steps 16; snapping down where hmax binds would
  move ordinary hmax runs by an ulp).
- Callers in the house that can see a change (rg for every function changed over the home's code trees, rgdirs.txt there, the kit's
  and the verification work's left out, 2026-11-20): the tire's mt10.c and its 38 copies (research/tire, a data disk/the development machine; mk_linear_eval: the same bits unless
  a node's value is -0 or its neighbour infinite), the reflex's rxtrain.c and its 27 copies (grown, a data disk/the development machine; mk_cmaes: only a refusal's
  count), the maintainer's strands/pack.c (mk_lbfgs: the same), scratch/dev-eval's repro t.c and t4.c (mk_quantile_sorted; mk_cg,
  mk_csr_from_triplets, mk_dopri). No caller of mk_nm, mk_q_axis_angle, mk_bicgstab, mk_ode23s or mk_bdf_step outside the kit and the verification work
  (callers162b.txt there); grown's vendored kit/mk_*.h copies keep their own code. the verification work's fuzz units for opt, ode, sparse and vec call
  changed functions.
- Rows: t_ode 20 -> 22, t_opt 24 -> 26, t_sparse 28 -> 31, t_vec 13 -> 15 (with the version row); test/drvctl.sh one control more. Each new
  row FAILs on 1.6.1's headers (scratch/mathkit-162-work/neg/negative-controls.txt; #1's row there through its dopri part). No new
  planted faults: 140.
- Checks: `jobwatch run -t 45m make check` (a second machine, its ledger: 693 s, end=ok; done by 10:38 2026-11-20): wrong 0; 259 verdict
  lines, all 140 planted faults caught by FAIL rows (41, 52, 54 and 101 now by one row more, 102 by two), the 27 driver and negative
  controls as they must be (J1's the new one), t_par's 17 BITS lines equal without OpenMP and byte for byte 1.6.1's (log
  scratch/mathkit-162-work/check.log). `jobwatch run -t 45m make check OMP=1` (a second machine, 533 s, end=ok; done by 10:47): wrong 0; 237 verdict
  lines, 140 planted faults caught, the 27 controls, 17 BITS lines equal (check-omp.log there). Every test once on 1.6.1 and on 1.6.2
  (q161/, q162/ there): only mk_nm's numbers move (above), every other printed number the same. The four changed tests also build at -O3
  -march=native -fno-math-errno -fno-builtin-pow -fno-trapping-math -funroll-loops -Werror and print the same rows; t_ode, t_opt,
  t_sparse, t_vec, t_interp and t_stats under -fsanitize=address,undefined, without stdbuf: no ASan, UBSan or LeakSanitizer report
  (scratch/mathkit-162-work/asan/run.log).

### 1.6.1 (installed 2026-11-20 01:43)

A second outside review of 1.5.6, read without compiling (scratch/the maintainer/api-eval/reports-batch2-fable/mathkit-fable-5-2.md, a
Fable model: 2 faults, 15 risks, 11 cleanups; #1 to #28 here). Each finding new since the first review was run on 1.6.0 first with a short
program: the claims table scratch/mathkit-161-work/claims/claims.tsv (the programs and outputs beside it, in claims/ and
claims/out/): 17 CONFIRMED, 2 NOT REPRODUCED (#11, #15's order half), 4 NOT A FAULT (#22, #23, #26, #27), and skipped as the first review's
or fixed in 1.6.0: #5's quantile half (S1), #6 (O1), #12 (R1), #15's INT_MAX half (P1), #16 (C1), #20 (M1), #21 (the freeze paragraph).
Staged from 1.6.0 in scratch/mathkit-161; work in scratch/mathkit-161-work (STATE.md). The experimental marks stay as 1.6.0 set them.
- #1 (mk_dopri, mk_ode23s, mk_bdf_step): the error and Newton norms were fmax(norm, q), and fmax drops a NaN, so a step with a NaN in it was
  taken: dopri and ode23s returned 0 with y NaN, bdf_step called Newton converged. The max keeps a NaN now (if (!(q <= err)) err = q), so
  the step is refused and shrunk: dopri on y' = sqrt(1 - t) over [0, 2] is -1 at y 0.666667 (2/3), bdf_step -1 with y as it was. The
  same bits wherever no NaN arises.
- #8 (mk_dopri, mk_ode23s): a step below t's ulp moved y and not t (y' = 1 from 1.7e9 over 0.01: off by 1.2e-5 at rtol 1e-8 with 0; from
  1e17 over 992: 982.75). Each step is (t + s) - t, the one t can take, refused when 0, and a default first step is at least 16 ulps of t0
  (MATLAB's hmin): 1.7e-16 and 2.3e-16 off, 0; y' = -y at 1e17 is refused at once. Ordinary results move at rounding level (claims
  c_f8b.c: relative 4e-16 to 8e-12; t_ode: Robertson's mass 0 -> 2.9e-15 off 1, Van der Pol 5529 -> 5526 steps with the same period,
  O-A's y(1) 1e-15 -> 2.2e-16 off).
- #2 (mk_nm): the spread test starts at the second iteration, so a start simplex on a level set no longer stops at its start (x^2 from
  -0.5, step 1: 0.25 -> 0 in 32 calls). An untied start moves nothing (t_opt's other rows as 1.6.0).
- #3 (mk_spec): mk_stft_frame windowed into a buffer of the const plan's, so two threads sharing one plan mixed their frames (161 and
  198 of 200 rounds wrong). The buffer is cut: the frame is packed straight into X and transformed by mk_rfft's second half
  (mk_rfft_packed_ in mk_fft.h; mk_rfft is the packing and that). The plan is only read, mk_spec has no buf and half the allocation, and
  every mk_welch, mk_csd and mk_stft_frame output is 1.6.0's bit for bit (claims c_f3b.c).
- #4 (mk_q_norm): the identity only for a zero quaternion; a NaN comes back NaN (1.6.0's identity let mk_q_integrate, mk_q_from_mat and
  mk_q_slerp turn a diverged attitude level).
- #5 (mk_median): 0.5 (lo + hi), and 0.5 lo + 0.5 hi where that is not finite (mk_lerp's idea; mk_lerp itself would move every even median
  by its second rounding): {1e308, 1.7e308} 1.35e308 (1.6.0 inf); of 2.2e6 random medians only the 2 that overflowed move.
- #10 (mk_cg, mk_bicgstab): 0 and relres were the recurrence's residual, which falls below the true one (1-D Laplacian of 5000 from x0 1e6
  random, tol 1e-12: relres 8.9e-13, true 1.5e-10, and 0). A pass now ends by forming b - A x; past tol the solver restarts from x, as a
  caller calling again would, and relres is that true residual. A restart that cannot take its first step (p A p underflowing at A near
  2^-1015, t_sparse's F11) returns 1 with x as it was. x and the iterations stay bit for bit where the true residual already met tol
  (t_par's seven solver cases: only relres moves, by up to 2e-7 of itself); more iterations where it did not (5009 -> 8720 above). F11
  is 1 now (true 1.11e-13 against tol 1e-13, where 1.6.0 said 0); its row takes 0 or 1 with x right.
- #18: CG returns -1 for a negative or NaN r z (an indefinite preconditioner; 1.6.0: 1); BiCGSTAB refuses MK_PC_MEANFREE (1.6.0 dropped it
  and iterated).
- #13 (mk_bvh_inside): where no ray votes, |winding number| > 1/2 (as mk_sdf), so a closed mesh wound inward reads inside there too, as
  its parity says elsewhere. Outward meshes unchanged.
- #19 (mk_bvh_free): nt cleared, so a query after a free or a build that failed for memory finds no tree (1.6.0: SIGSEGV on node NULL).
- #14 (mk_snet_faces): a face needing a vertex the net has not (a grid other than its own) is -1 (1.6.0 read V[-3]: SIGSEGV, or garbage).
- #7 (test/t_hhd.c): the 1080p row's time is printed, not checked: at 5 % of one CPU 1.6.0's real header read WRONG (3.60 s against the
  row's 2 s), 1.6.1's passes (claims/out/f7-quota*.txt).
- #25: the Makefile's CC ?= gcc cut (make predefines CC, so it never acted; cc is gcc 16.2.1 on the development machine and a second machine). #24: README's three stale
  Resume here blocks cut (1.6.0's README keeps them); mk_par.h's row counts corrected.
- Documented, not changed: #9 (mk_lu's whole-matrix pivot rule turns away ode23s's steps past about 1e10 on Robertson's, the answers
  right; mk_solve3 takes rounding-small pivots, 0 with x near 9e15 on [[1 2 3] [4 5 6] [7 8 9]], no caller outside the verification work's fuzz; and a
  limit met on the way: the 1e-14 span floor makes Robertson's 0 to 1e11 in one call -1, while 0 to 1 then 1 to 1e11 gives Hairer and
  Wanner's), #11 (every n 21..2000 and 2^11..2^16 settles, weights to 1.75e-15; O(n^2), n 65536 in 23 s), #15 (glibc 2.44 sums repeated
  positions in input order), #17 (an infinite value makes mk_welford's mean NaN from the next value on; the finite shift proposed would
  change only the one-value case), #27 (LP64 only), #28 (0 points past the rules' degrees). Not taken: #22 (a conventions table), #23
  (the kit is checked on this house's machines), #26 (MK_MUTATE is test-only, documented).
- Callers in the house that can see a change (rg over the home, kit and the verification work copies left out, 2026-11-20): research/tire/v10-fast's
  launchm.c, its ref copy and launchrows.c (mk_welch on one thread: the same bits); grown/stage/bone/check/femurscore_hc.c,
  scratch/body-claude/trial-c1/src/partlib.h, scratch/the maintainer/judge/partjudge.c, scratch/glance-slices/judge/bonejudge.h and
  femurscore.c (mk_bvh_inside: only an inward-wound mesh where no ray votes reads differently); forms' snetobj.c copies (mk_snet_faces on
  their own grids: the same); scratch/dev-eval's repro t4.c (mk_cg, mk_dopri, mk_q_norm). No caller of mk_ode23s, mk_bdf_step, mk_nm,
  mk_median, mk_csd or mk_stft_frame outside the kit and the verification work. the verification work's fuzz units for opt, ode, sparse, spec, stats, vec, geom and
  snet call changed functions.
- Rows: t_ode 18 -> 20, t_opt 23 -> 24, t_par 24 -> 25 (with -fopenmp), t_vec 12 -> 13 (with the version row), t_stats 15 -> 16, t_sparse
  26 -> 28, t_geom 22 -> 24, t_snet 20 -> 21; t_sparse's F11 takes 0 or 1. Each new row FAILs on 1.6.0's headers
  (scratch/mathkit-161-work/neg/negative-controls.txt). No new planted faults: 140.
- Checks: `jobwatch run -t 45m make check` (the development machine, 2026-11-20, done by 01:33, wall 859 s): wrong 0; 258 verdict lines, all 140 planted
  faults caught by FAIL rows (41, 52, 54, 61, 92, 94 and 101 now by one or two rows more), the 26 driver and negative controls as they must
  be, t_par's 17 BITS lines equal without OpenMP (log scratch/mathkit-161-work/check.log). `jobwatch run -t 45m make check OMP=1` (a second machine,
  done by 01:41, wall 539 s): wrong 0; 236 verdict lines, 140 planted faults caught, the 26 controls, 17 BITS lines equal (check-omp.log
  there). Against 1.6.0's BITS: the 7 solver lines move (relres only: x and the iterations the same), the other 10 byte for byte. Every
  test also builds at -O3 -march=native -fno-math-errno -fno-trapping-math -funroll-loops -Werror. The 11 tests of the changed headers
  under -fsanitize=address,undefined, without stdbuf: no ASan, UBSan or LeakSanitizer report (scratch/mathkit-161-work/asan/run.log).

### 1.6.0, the freeze (installed 2026-11-20 00:03)

The owner's "finalize our most common tools" (through the maintainer), on an outside review of 1.5.6 read without compiling
(scratch/the maintainer/api-eval/reports-batch1/mathkit-opus-5-5.md: 22 findings, 3 fault, 7 risk, 12 cleanup). Each finding was
first run on 1.5.6 with a short program: the claims table scratch/mathkit-160-work/claims/claims.tsv (the programs and their
outputs beside it, in claims/ and claims/out/), 15 CONFIRMED, 3 NOT REPRODUCED (D2, P1, G1), 4 NOT A FAULT (P2, C2, G3, M1). Staged from
1.5.6 in scratch/mathkit-160; work in scratch/mathkit-160-work (STATE.md).
- The freeze: the numerical core keeps its calls and results from here on. Experimental, not frozen (a line in each header's top comment
  and a mark in the table): mk_dec, mk_sepvol, mk_graft, mk_ovl, mk_hhd, and mk_fftm in mk_fft.h.
- D1 (mk_dec): mk_dec_nbrs_ stamped MARKS++, which is exactly the MARKS + 1 that mk_dec_cheap_ gives the edge's two common vertices,
  so the neighbour list taken right after a collapse left out both vertices opposite the collapsed edge: they were not re-rated, and a
  queued entry of theirs aimed at the dead vertex was dropped. On the icosphere 40798 of 80259 neighbour lists were short (50.8 %;
  torus 49.0 %). Now MARKS += 2. mk_decimate's results move (t_dec, against 1.5.6): unit sphere at 0.005 878 -> 851 vertices, held
  vertices 1554 -> 1533, torus at 0.01 1086 -> 1045, hollow ball 67 -> 79, open cap 293 -> 274 (rim 99 -> 98), the 12 x 3 torus 11 -> 11
  with 189 -> 262 link refusals; a target of 2000 still 2000; t_regr's D-2 924 -> 670 vertices left; t_graft's patch thinned 664 -> 655
  (rim 66 -> 68), its stitched volume 4.13611 -> 4.13613 (S's less the lens: 4.13634). Every bound row still holds and meshcheck reads
  dec_sphere, dec_torus, dec_hollow and graft clean. It is no longer byte for byte folddec's:
  folddec has the same two lines (scratch/mesher-body/skel7/src/folddec.c, nbrs and cheap). Row: an octahedron through the internals,
  the end's neighbour list after a link test (1.5.6: 2 of 4).
- S1, I1: mk_vec.h's new mk_lerp(a, b, t) is the kit's one linear interpolation: a + t (b - a), b exactly at t = 1, and (1 - t) a + t b
  where b - a overflows. mk_quantile_sorted, mk_linear_eval, mk_hermite_eval and mk_pchip_eval (n 2) and mk_v3_lerp use it (mk_interp.h
  and mk_stats.h now include mk_vec.h). Bits that move (scratch/mathkit-160-work/lerp: 2e6 random queries each, values 1e-20 to 1e20):
  none inside a table (mk_linear_eval) and none in mk_quantile_sorted (its fraction is never 1); at or past the last node 55.9 % of
  queries, now y[n-1] exactly where 1.5.6 gave y[n-2] + (y[n-1] - y[n-2]): 1 or 2 ulps when the two are of a size, all of y[n-1] when it
  is far smaller (y {1, 1e-17}: 1.5.6 gave 0); and the overflows (inf before, now finite). Cost, best of 7 on the development machine: mk_linear_eval on 64
  nodes 23.35 -> 23.86 ns a call, mk_v3_lerp 1.04 -> 1.49 ns. aeromap.c's pchip stays bit for bit except a 2-point table at its end.
- O1: mk_ode23s refuses t1 < t0 (-1); 1.5.6 answered 0 with y unchanged. C1: mk_cmaes refuses lambda 1 (mu 0) before evaluating; 1.5.6
  ran a generation and returned status -1 with no reason. B1: mk_csr_from_bsr3_lower refuses a block above the diagonal; 1.5.6 mirrored it
  and summed it twice. N1: mk_sdf_sample gives NaN at a NaN coordinate (1.5.6: node 0's value). N2: mk_spec_init refuses a window outside
  MK_WIN_RECT .. MK_WIN_KAISER and an fs not positive and finite (1.5.6: rectangular, and fs 1). Q1: mk_q_axis_angle with a zero axis is
  the identity (1.5.6: (cos(a/2), 0, 0, 0), not unit). A row each.
- G2: mk_sepvol's own crossing pairs are mk_geom's mk_tri_cross (mk_bvh_crossings' and mk_decimate's test, meshcheck's), not its own
  mk_ray_tri edge test: mk_sepvol_share_, _edge_strict_ and _cross_ cut. The two agree on 2e6 generic pairs and disagree on 2.3 to 13.9 %
  of near-degenerate ones (claims/out/g2.txt). One row's count moves: the held mesh crossing itself (ballA and a copy 2 mm along x,
  near-coplanar faces) 204 -> 130 crossing pairs, mk_bvh_crossings' count (a long-double referee reads 100 there: such pairs have no
  firm answer); its status 4 as before. G3: mk_bvh_count_hits is mk_bvh_count_wind_ (the one traversal; the counts were equal) and
  mk_bvh_count_hits_ is cut (its one caller, mk_bvh_inside, calls mk_bvh_count_wind_).
- Documented, not changed: D2 (32-bit marks, about 48 a collapse: past about 9e7 collapses a stale mark can read as current), P1
  (mk_csr_from_triplets: at most INT_MAX distinct positions), B2 (mk_csr_matvec: x and y must not overlap), T1 (the ctl counters add up over
  calls), G1 (a zero ray direction hits nothing: MEASURED, not arbitrary), H1 (mk_hhd: 2 ny's largest prime factor goes by its direct sum:
  1920 x 1031 10.7 s against 0.175 s at 1080; common heights have factors up to 67), C2 (to stop mk_cmaes at f = 0, give 0x1p-1074).
- test/judge.sh (J1): the timeout test comes before the rows test, so a planted fault that prints 124 FAIL rows and hangs is WRONG, as
  the rule says, not caught (c_j1.sh).
- Not changed: R1 (stream s + 1's seed words are stream s's shifted by one; no equal draws within 64 places and no correlation over 1e6
  draws: a new seeding would move every seeded result and ribbon's bit-equality for no measured gain), P2 (L's entries are A's, under
  INT_MAX), M1 (mk_sepvol's MK_SEPVOL_PAIRS, _CRUSH and _SMOOTH stay macros while mk_sepvol is experimental).
- Callers in the house that can see a change (rg over the home, kit copies left out, 2026-11-19): research/tire's ref mt10.c and its
  staged copies (mk_linear_eval: at or past the mu table's last node now that node's value exactly); research/tire/v10-fast's launchm.c
  and launchrows.c (mk_spec_init: an fs not positive and finite now -1); grown/stage/reflex/src/rxtrain.c (RX_LAMBDA=1 now refused);
  grown/stage/joint/skel and grown/stage/bone/check (mk_sdf_sample: NaN in, NaN out); scratch/mesher-helper's sepheld.c and
  grown/stage/joint/skel/heldsep.c (mk_sepvol: own crossing counts on near-degenerate pairs). No caller of mk_bvh_count_hits_,
  mk_v3_lerp, mk_csr_from_bsr3_lower, mk_q_axis_angle or mk_ode23s was found outside the kit; metal's tire10 law uses mk_tridiag only.
- Rows: t_vec 11 -> 12, t_stats 14 -> 15, t_interp 9 -> 10, t_ode 17 -> 18, t_opt 22 -> 23, t_sparse 25 -> 26 (the table said 23), t_spec 7 -> 8,
  t_sdf 16 -> 17, t_dec 14 -> 15. Each new row FAILs when built on 1.5.6's headers (scratch/mathkit-160-work/neg/negative-controls.txt).
  No new planted faults: 140.
- Checks: `jobwatch run -t 45m make check` (the development machine, 2026-11-19, done by 23:52, wall 701 s): wrong 0; 258 verdict lines, all 140 planted
  faults caught by FAIL rows, the 26 driver and negative controls as they must be, t_par's 17 BITS lines equal without OpenMP and byte
  for byte 1.5.6's (log scratch/mathkit-160-work/check.log). `jobwatch run -t 45m make check OMP=1` (the development machine, done by 00:01 on
  2026-11-20, wall 551 s): wrong 0; 236 verdict lines (the -fopenmp builds of each header alone are the tests' own here), 140 planted
  faults caught, the 26 controls, 17 BITS lines equal (log check-omp.log there).
  The 11 changed tests under -fsanitize=address,undefined, without stdbuf: no ASan or UBSan report; LeakSanitizer finds t_sdf's own test
  meshes left unfreed (1078272 bytes in 12 blocks), the same on 1.5.6's t_sdf (scratch/mathkit-160-work/asan/run.log).

### 1.5.6 (staged and installed 2026-11-19 18:06)

The owner's ask (2026-11-18/19, through the maintainer): the Helmholtz-Hodge decomposition, applied first to motion fields (optic flow and
motion vectors) for the maintainer's video frame generation, retina's and neuro's vision and fly's motion pathway. Staged from 1.5.5 in
scratch/mathkit-156 by the maintainer's helper; work, probes and pictures in scratch/hhd-work (STATE.md).
- New mk_hhd.h. mk_hhd_init(H, nx, ny), mk_hhd_split(H, u, v) fills the gradient part grad phi, the rotational part
  rot psi = (dpsi/dy, -dpsi/dx), the harmonic rest h and both potentials; mk_hhd_readout(H, cx, cy, R) gives the expansion rate (half
  the mean divergence: W / Z for an approach), the rotation rate (half the mean curl, positive turning +x toward +y), the mean
  shift (the harmonic part's mean, which is the field's), the focus of expansion (the least-squares point the flow points away
  from, a rigid swirl at the rotation rate taken out about the principal point) and each part's share of the energy.
- Boundary choice: phi = psi = 0 on the frame (Dirichlet on both, the Hodge-Morrey "normal-parallel" split). the maintainer's brief
  suggested phi Neumann; on a rectangle that leaves no harmonic part at all (a harmonic field with no flux through a simply
  connected frame is 0), so a uniform shift would land in the gradient part. With both Dirichlet the three parts are orthogonal
  and a shift is harmonic, as asked. The price: the flow that crosses the frame is harmonic, so a swirl or an expansion centred
  on a 16:9 frame keeps 62.3 % of its energy in its own part and 37.7 % in h (84.4 % on a square; MEASURED at 256 x 144,
  1920 x 1080 and 144 x 144). The rates are exact all the same.
- Discrete form: potentials at the pixel corners, each pixel's gradient the 2 x 2 box difference, so div(rot psi) = 0 and
  curl(grad phi) = 0 to rounding and the parts are orthogonal in the plain pixel sum; both solves have one operator, whose sine
  modes have the values 2 (1 - cos(pi k / nx) cos(pi l / ny)): a sine transform along y and mk_tridiag along x for each mode.
  The box differences are exact on quadratic fields, which a camera's flow over a plane is.
- Speed (MEASURED 2026-11-19 on a second machine, Ryzen 7 9700X, median of 7, scratch/hhd-work/hhdtime.c): 1920 x 1080 split 0.139 s at
  -O2 on one thread, 0.119 s at -O3 -march=native, 0.030 s with -fopenmp on 8 threads; readout 0.007 s. 3840 x 2160: 0.60, 0.50,
  0.12 s. Memory about 170 MB at 1080p (six parts, two potentials, work), 770 MB at 4K (peak PSS).
- mk_fft.h: mk_fftm, a mixed-radix complex FFT for any length (Stockham; radices 4, 2, 3 by hand, other primes by their direct sum).
  The radix-2 plan cannot do 1080p's sine transforms (lengths 2160 and 3840); within 2.6e-15 of a direct DFT at 21 lengths up to
  2160. No house copy of a mixed-radix FFT or of a Helmholtz-Hodge split was found (rg over pix, fly, neuro, spectrum, eye,
  gauge, metal, mcphys, cudakit). signatures.tsv gains rows hhd and fft_any_length (this tree's copy; global/methods/
  signatures.tsv is the maintainer's to update at install), each matching its header here and nothing in 1.5.5.
- What the camera rows recover (t_hhd; MADE UP camera, Longuet-Higgins and Prazdny 1980's flow): on a wall facing the camera the
  focus at the heading to 1e-6 px, W / Z and the roll to 1e-12. A pan's near-uniform flow lands in h (99.997 %, 99.71 % and
  96.73 % harmonic at 15, 46 and 81 degrees across, the rest gradient as the view widens) and moves the focus read by f B Z / W
  (40.03 px against 40 to first order): one plane's flow cannot tell a pan from a turned heading, so the read point is the flow's
  zero, not the heading. A tilted plane's translation has a shear of its own, which the rotation read carries (0.00125 rad/s,
  exactly its formula) and which leaves the focus 0.24 px off; a box room's leaves it 0.12 px off, 0.4 % on the rotation.
- Rows: t_hhd 33 (new), t_fft 12 -> 17, t_par 22 -> 24 rows (built with -fopenmp) and 16 -> 17 BITS lines (mk_hhd_split the same bytes at 1, 2 and 12 threads and
  without OpenMP), t_vec's version row. Planted faults 238-245, 132 -> 140 (238 is t_fft's only: mk_hhd runs forward transforms).
  t_hhd and t_fft run clean under -fsanitize=address,undefined. `jobwatch run -t 30m make check` (on a second machine, 2026-11-19 14:32 to
  14:44, wall 686 s): wrong 0, 258 PASS lines, all 140 planted faults caught by FAIL rows, the 26 driver and negative controls as they
  must be, t_par's 17 BITS lines equal without OpenMP (log scratch/hhd-work/check.log). make check OMP=1 not run for 1.5.6.
- Pictures (viewed): scratch/hhd-work/hhd-plane.png (the main row: forward, heading, roll, yaw) and hhd-room.png (the box room),
  the input and its three parts as arrows on one scale; scratch/hhd-work/hhdpic.c draws them.
- Known limits, none patched: a pattern alternating pixel by pixel has no box divergence or curl and stays in h; the focus read
  takes the swirl about the principal point given (the frame's centre by default) and does not separate a pan from the heading;
  a scene's shear counts in the rotation read; doubles only; no consumer uses it yet.

### 1.5.5 (staged 2026-11-18, installed 2026-11-19 00:38)

Codex's evaluation of 1.5.3 (reports/mathkit-1.5.3-evaluation-2026-11-18.md: the driver behind "wrong: 0", and three public
geometry cases) and the verification work's findings F-B and F-C on mk_sepvol; the regression rows 1.5.4 owed for the verification work's 1.5.3 findings. Staged from
1.5.4 in scratch/mathkit-155; work, drivers and logs in scratch/mathkit-155-work (STATE.md).
- The driver (test/check.sh, new test/judge.sh, test/negctl.sh, new test/drvctl.sh). A planted fault counts as caught only by a FAIL row: a
  timeout never does, nor a crash or a bare nonzero exit (Codex replayed exit 86 with no row passing); a fault meant to be caught by a crash
  would be named in CRASHES with its signal, and none is. The real header needs exit 0, no FAIL row and at least one PASS row. Every run
  starts from no earlier output, so a t_par that does not build can no longer pass the BITS comparison on stale lines; a failed build or run
  there is WRONG, and the number of BITS lines is checked (16). negctl's control 4 read a build without OpenMP that printed no BITS line as
  "caught" (an empty list of lines differs); it now needs both builds and runs, and as many lines from each. `exit $wrong` wrapped at 256; it is clamped to 100. CHECK flushes each row, so a row printed
  before a crash counts. test/drvctl.sh shows each rule failing on a known-bad input (21 controls: a fault that exits 86, one that hangs, a
  crash on and off the list, a wrapped count, a t_par that does not compile with stale output present, a negctl build that fails and one
  that prints no BITS line, which 1.5.4's negctl read as caught, the clamp).
- What the stricter rule exposed: fault 201 alone, caught by 0 rows in every earlier release (its write before a block aborted the run).
  t_dec.c now gives each of the header's blocks guard bytes, checked when freed: 201 fails 2 rows and no longer corrupts the heap.
- mk_sepvol_two_ (F-C): a point covered three times has odd parity and read as covered once. On a closed mesh the clean rays' facings decide
  whatever the parity (w0 > 1). On a mesh not closed the winding number decides when a clean ray's facings exceed 1, and odd parity with
  no such ray stays "covered once", which keeps mk_winding_number's O(triangles) off every vertex of an open mesh. Three spheres: 112 ->
  159 inside, the other parts' winding number's count; nested tetrahedra of 3 and 4 shells: 4 and 8 -> 8 and 12. Cost: no new ray, and
  on the kit's rows mk_winding_number is called as often as on 1.5.4 (239 calls in t_sepvol, its 1.5.4 rows, both ways; 3.78 s against 3.78 s,
  the fastest of 7 interleaved runs at load 16-24). The only results that change are vertices covered three or more times: t_sepvol's rows
  and all 11 of its result meshes are 1.5.4's, byte for byte.
- The held mesh (Codex): only the total volume's sign was tested, so a separate part wound inward passed. Each connected part of negative
  volume must now lie inside another part (the winding number of the rest of the mesh over 1/2 at one of its vertices): a cavity, as in a
  hollow solid, is taken; a part wound inward on its own returns -2 with MK_SEPVOL_HELD_FLIPPED (rep->held_inward counts them).
- Crossed bars (F-B). No vertex move clears parts that cross with no vertex inside each other, and step 2 only squeezes them. A run now stops
  when a part of a closed M would lose more than MK_SEPVOL_CRUSH (0.35) of its volume with vertices inside or pairs crossing still left, and
  hands back the best state before that (MK_SEPVOL_CRUSHED). MK_SEPVOL_STILL_CROSSES is set whenever pairs are left at the end, so a caller
  is never told success on parts that still cross; rep->part_loss gives the largest share one part lost. Why 0.35: the largest loss in any
  run that ended with no pairs left is 13.5 % (t_sepvol's plate), and at most 0.17 % on the six held real cases; the plate squeezed with
  pairs left, a known limit, reaches 33.0 %, and stays as it was; the bars go on to 65 %, a real vessel to 48 % in one round.

  | case | version | rounds | crossing pairs | volume | status |
  |---|---|---|---|---|---|
  | two crossed bars (the verification work's) | 1.5.4 | 15 | 20 -> 12 | 0.4800 -> 0.2059 (-57.1 %) | 0 |
  | | 1.5.5 | 8 | 20 -> 12 | 0.4800 -> 0.3284 (-31.6 %) | CRUSHED, STILL_CROSSES |
  | artery, single mode, crossing itself (PRIVATE) | 1.5.4 | 9 | own 174 -> 173 | 0.01 -> -0.00 ml | SELF_CROSSES |
  | | 1.5.5 | 1, nothing moved | own 174 -> 174 | 0.01 ml | SELF_CROSSES, STILL_CROSSES, CRUSHED |
  | nerve, single mode, crossing itself (PRIVATE) | 1.5.4 | 9 | own 356 -> 345 | 0.12 -> 0.05 ml | SELF_CROSSES |
  | | 1.5.5 | 3 | own 356 -> 345 | 0.12 -> 0.08 ml | SELF_CROSSES, STILL_CROSSES, CRUSHED |
  | median nerve, single mode (PRIVATE) | both | 9 | own 167 -> 164 | 1.06 -> 0.93 ml, the same bytes | 1.5.5 adds STILL_CROSSES |
  | six held cases, heart and liver against bones and muscles (PRIVATE) | both | 3 to 30 | the same | the same bytes | 1.5.5 adds STILL_CROSSES on the four with pairs left |

  gap 0.1 mm, step 0.25 mm, 4 rings, 30 rounds, as mesher-helper's sepheld calls it (scratch/mathkit-155-work/rc155/table-final.txt). The kit's
  own rows do not change.

- Regression rows (t_regr, 20 rows, new): O-1, O-1b, O-4, O-5, O-5b, O-6, O-7, O-7b, H-end, D-2, D-5, D-6, G-1, G-2, G-3, G-5, R-1, G1, S-1, S-2.
  All 20 fail on 1.5.3's headers and pass here. D-3 has no row: reaching it takes about 6e7 collapses. O-1 and O-1b read within 7 % of the verification work's
  definition, because at those coordinates the lines' positions round to more than their spacing.
- Rows: t_sepvol 63 -> 73, t_dec 13 -> 14, t_regr 20 (new), t_vec's version row. Planted faults 216-237, all caught: 110 -> 132 (216-221 in
  t_sepvol, 222-237 in t_regr, one for each patch of 1.5.4 those rows guard; scratch/mathkit-155-work/mutmap.log lists the rows each
  fails). make check and make check OMP=1: wrong 0 each (778 s and 2382 s wall, 2026-11-18 23:32 to 11-19 00:37, the second on a loaded the development machine), every one of the 132 caught by at least one FAIL
  row, the 21 driver controls and 5 negative controls as they must be, t_par's 16 BITS lines equal to 1.5.4's. The new and changed tests
  run clean under -fsanitize=address,undefined (t_sepvol 73, t_regr 20, t_dec 14).
- API: the report grows two fields at its end (held_inward, part_loss), and two status bits (16 MK_SEPVOL_STILL_CROSSES, 32 MK_SEPVOL_CRUSHED); a
  caller that tests status == 0 now reads pairs left as not clean. mk_sepvol, mk_sepvol_mark_ and mk_sepvol_two_ keep their signatures.
  mk_sepvol_heldcheck_ keeps its signature; mk_sepvol now also allocates two more arrays per call (45 allocations in row 13's case, from 43).

### 1.5.4 (2026-11-18, installed 22:07)

the verification work's fuzz findings against 1.5.3 (the maintainer.md #68, #70; the verification work/RESULT.md Results 93-94), its patch applied as written
(mathkit-1.5.3-qc-2026-11-17.patch, sha256 8100ffafaa969f16).
- mk_ovl: O-1 and O-1b, a hang: the double loop counters stopped moving past 2^53 (now integer lattice counters). O-4 to O-7, wrong
  or NaN results at extreme scales (the lines' origin is the point of the axis span nearest 0; `x * h * h`). The hit search is closed
  at both ends as its comment says. O-2 is fixed only in part: the rest is mk_ray_tri's T overflowing near 1e308.
- mk_dec: D-2, a large face entered in every cell of its box (the verification work's demo at 128 cells a side, the maintainer's run: 28.4 s on 1.5.3, 0.005 s here; the verification work's own run read 384 s; MK_DEC_BIG). D-3, int marks
  overflowed after about 6e7 collapses (now unsigned). D-5, a positive bound under about 1.5e-162 squared to 0, meaning no bound.
  D-6, maxown counted collapses that did not happen.
- mk_graft: G-1, lo's rounded sum refused a finite lo. G-2, a far node gave infinite coordinates. G-3, psi = +inf gave NaN
  coordinates (bisect). G-5, a failed second malloc freed twice. R-1, mk_graft_refine accepted NaN or inf coordinates. Doc lines
  for G-6 and W-1.
- mk_geom: G1, mk_aabb_ray's pad ignored the ray origin's size, so a ray from far away lost hits and mk_bvh_count_wind_'s facing sum
  left the winding number (unit cube, origin about 1.7e12).
- mk_sdf: S-1, a refused mk_sdf_mesh_box left axis 0's fields in *g. S-2, the (int) cell cast of a triangle over INT_MAX cells from
  node 0 was undefined (clamped first; bit-identical inside the grid). Its demo: 125 of 125 nodes read outside on 1.5.3, 0 here at 8.6e9 cells; at 1.1e12 cells with flags 4, 15 still do (N-1 below).
- Known limits, not patched: mk_dec D-4 (time quadratic in a hub vertex's valence); mk_graft G-7 (an ear can leave a mirrored pair);
  mk_sdf S-3 (a one-node-thick box gets no sweeps; far nodes keep BIG, still an upper bound) and N-1 (above about 2^39 h a
  face-diagonal ray is counted twice).
- Regression rows for these findings follow in the next release (1.5.5: t_regr); the verification work's demos (scratch/the verification work/mk-153f/findings) are their
  negative controls until then: they fault on 1.5.3 and are clean here.

### 1.5.3 (2026-11-16, installed 19:44)

the verification work's findings on 1.5.2's mk_sepvol (the maintainer.md #65, helper H16) and H18's mk_par.h overflow. 1.5.1 is spent. A patch: the
public signatures and the report are unchanged; results move only where 1.5.2 read a point eps inside a vertex as covered twice
when it was covered by nothing (F-A), and every threaded block, sum and BITS line is 1.5.2's. Staged from 1.5.2 in
scratch/mathkit-153; work, drivers and logs in scratch/mathkit-153-work (STATE.md).
- the verification work's patch (scratch/the verification work/mathkit/findings/mathkit-1.5.2-sepvol-qc-2026-11-16.patch, sha256 30d712ec79954dec...),
  read and applied as it came, three hunks in mk_sepvol.h:
  - F-1 (contract): when an allocation other than V0's failed, `goto out` copied V0's uninitialised block over V. V0 is now
    filled right after its malloc. Row 13: each of a single-mode call's 43 allocations made to fail in turn (two cubes 1 mm into
    each other, 3 rounds; while a failure is set each fresh block is filled with 0xA5, test-only, after the verification work's demo.c): all 43
    return -1 with V exactly as it came in and the report zeroed. 1.5.2: 22 of 42 (each failure after V0's malloc and before its
    copy wrote 0xA5 bytes over V). Planted fault 212.
  - F-A: a part thinner than eps (a thousandth of the mean edge) read as inside another: the point eps inside its vertex lands
    outside every part, and even parity there read "covered twice". Even parity is now confirmed by the winding number (over 1.5).
    Row 14: the verification work's tiny cube (0.2 mm beside a 2 m cube): inside0 0, moved 0 (1.5.2: 8 inside). Planted fault 213.
- F-A's cost, and the speed-up. The patch calls mk_winding_number, O(triangles), at every even-parity point, which is nearly every
  inside vertex, every round. On the kit's single-mode rows it is lost in the noise (1.5.2's t_sepvol.c: 1611 winding numbers
  against 2; the whole test 4.9 s with it, 4.6 s on 1.5.3, 5.5 s on 1.5.2, medians of 5 interleaved runs at load 11-18). On
  orig02c (below) it made the call 3.6 times as long: 1904 s against 532 s for the same 29 rounds without it (25,948 winding
  numbers, 1.95e10 triangle terms). 1.5.3 keeps the answers and drops the cost: on a closed mesh (every edge run as often one way
  as the other, checked once a call, mk_sepvol_closed_) the winding number is the hits' facings summed along any ray with no edge
  hit, exactly, wherever no hit is rounding-close to the ray's start. The three rays mk_sepvol_two_ already casts now sum the
  facings too (mk_bvh_count_wind_; the facing is the sign of the watertight test's det, mk_ray_tri_s_, which mk_ray_tri now wraps
  with the same arithmetic). When the clean rays agree, their sum decides (2 or more: inside); when they disagree, or the mesh has
  open or one-way edges, mk_winding_number does, as in the patch. Measured against the patch's brute force: orig02c's output
  byte for byte the same and the report the same, with 0 winding numbers (532 s, 17.7 s a round against 1.5.2's 19.0 s);
  t_sepvol's outputs the same; the verification work's demo.c the same as its demo150_patched.txt; row 15 at every vertex of twoball, the two
  cubes and F-A's pair (2090 points): the same decision, 0 different, with even points summing to 2 and to 0; and orig02c with
  both at every even-parity point of every round (a work copy of the header): 25,948 points, all summed by the facings (24,276 to
  2, 1,672 to 0, none to anything else), 0 decisions unlike mk_winding_number's, the largest difference from it 2.2e-13. (That
  run, brute force and facings both, took 1173 s at load 5-10, against the patch's 1904 s at load 12-20: the times move with the
  load.) Planted
  faults 214 (the facing turned: t_geom's new row and 5 of t_sepvol's) and 215 (every mesh taken as closed: row 16, an open sheet
  through a cube, where the facings read 57 vertices inside and the winding number, at most 1.484, none).
- mk_par.h (H18): mk_par_blk's (n + 255) / 256 and mk_par_nblk's (n + b - 1) / b overflowed a long within b of LONG_MAX (UBSan at
  mk_par.h:36 and :40 for n = 9223371998217895936), and mk_par_range's (k + 1) b on the last block at n = LONG_MAX. Now
  n / k + (n % k != 0), and b < n - k b. The same answer wherever 1.5.2's did not overflow: par/parcmp.c compares the two at
  71,303,877 counts (every n from -2^20 to 2^26, every block's range up to 2^22, 2^20 counts on each side of LONG_MAX and LONG_MIN,
  every 2^j +- 300, a million log-uniform): equal at all 70,254,999 where 1.5.2 stays in a long, 0 different, and 1.5.3's blocks
  tile [0, n) in at most 256 at 5,248,429. t_par's 16 BITS lines are 1.5.2's, with and without OpenMP. Row (t_par): seven counts
  from 2^62 + 1 to LONG_MAX cut into at most 256 blocks of the exact ceiling's length, tiling [0, n): 1.5.2 fails it (2 of 7) and
  its UBSan build stops at mk_par.h:40; 1.5.3 runs clean under ASan and UBSan. Planted fault 211. Not touched: mk_geom.h's batch
  calls use the same form ((n + 63) / 64), but their n points take 24 n bytes, so no n near LONG_MAX reaches them.
- Known limits, not patched (the verification work's demos, scratch/the verification work/mk-sepvol/findings/demo.c, run on 1.5.3):
  - F-B: two long bars crossed like an X are squeezed by step 2's X fallback (two 6 x 0.2 x 0.2 bars: volume 0.48 -> 0.30 with 8
    vertices left inside after the stall on 1.5.2; on 1.5.3, with F-A's fix, 0.48 -> 0.21, none inside, 15 rounds, worst 0.51);
  - F-C: a point covered by three parts reads as covered once (parity odd), so it is not moved (three balls: 159 vertices inside
    another part by the winding number, 47 of them covered three times; inside0 112);
  - held cubes near 2^540 are refused with -2 FLIPPED: the volume check overflows (held_vol NaN).
- orig02c (cad-helper's skin, 375,000 vertices, 749,996 faces, one closed component crossing itself; PRIVATE, LAN only, not in the
  tests): one mk_sepvol call, single mode, gap 0.05 mm, step 0.25 mm, 3 rings, 40 rounds, the development machine at load 12-20:

  | | time | rounds | crossing pairs | inside | moved | mean move | worst | volume |
  |---|---|---|---|---|---|---|---|---|
  | 1.5.2 | 419 s | 21 | 1307 -> 1236 | 904 -> 837 | 5754 | 0.151 mm | 2.58 mm | -1.53 ml |
  | the verification work's patch | 1904 s | 29 | 1307 -> 1241 | 867 -> 799 | 4980 | 0.236 mm | 2.93 mm | -2.11 ml |
  | 1.5.3 | 532 s | 29 | the same, byte for byte | | | | | |
  | cad-helper (sepvol 18 rounds, meshfix 42) | | | 1307 -> 0 | 891 (its count) | 4161 | 0.38 mm | 3.97 mm | -2.1 ml |

  meshcheck reads the same crossing pairs (1236, 1241) and volume (60723.3 ml from 60725.4) in the results. 37 vertices 1.5.2
  read as inside lie where the point eps inside leaves the body (F-A). Every crossing pair here is the one part's
  own (status 8, self0 1307), and mk_sepvol never undoes those by design: the guard refuses a move that adds an own crossing
  (35,446 refused), and the run stalls. For a part crossing itself cad-helper's chain stays the tool.
- Rows: t_sepvol 59 -> 63 (13 F-1, 14 F-A, 15 the facings against the winding number, 16 an open sheet; 15 and 16 call 1.5.3's
  internals and sit under #ifndef T_OLD), t_geom 21 -> 22 (the facings summed along 3000 rays against the winding number of two
  spheres sunk into each other), t_par 5 -> 6 (16 + 6 with OpenMP), t_vec's version row. New planted faults 211-215, all caught.
  On 1.5.2's headers (the new tests with -DT_OLD) rows 13, 14 and t_par's new row fail. The new rows run clean under
  -fsanitize=address,undefined (t_sepvol 63, t_par 6, t_geom 22, t_vec 11). make check and make check OMP=1: wrong 0, 110
  planted faults caught in each, t_par's BITS equal to 1.5.2's in each.
- For the verification work: H16's sv.h (op_mark, line 737) calls mk_sepvol_mark_ with 1.5.2's arguments; 1.5.3 adds `closed` before eps
  (mk_sepvol_closed_(T, nt) in single mode, 0 in held), and mk_sepvol_two_ takes it too.
- README: line 1, the table (mk_version "1.5.3"; mk_par 16 + 6; mk_geom 22 + 3 and mk_bvh_count_wind_; mk_sepvol 63 + 11
  meshcheck: it read 25 + 5, stale since 1.4.x, 59 + 11 at 1.5.2), what the checks test, the planted-fault list (110; 1.5.2's
  201-210 line added) and count.

### 1.5.2 (2026-11-16, installed 11:39)

mk_sepvol, three findings from the maintainer's C1 race and mesher-helper (the maintainer.md #55, #56, #57; PRIVATE meshes, LAN only, not in the
tests). 1.5.1 is spent. A patch: one switch added (MK_SEPVOL_PAIRS, default on), the signature and report unchanged; results move
only in held mode where a vertex's nearest exit lies ahead. Staged from 1.5.0; only mk_sepvol.h, t_sepvol.c and check.sh's sepvol
lines changed (with the version). Work, drivers and logs: scratch/mathkit-152-work (STATE.md).
- the maintainer's runaway (#55, race-c1 H2: brachioradialis r, 2947 vertices, mean edge 2.2 mm, against its elbow's capsule and four
  ligaments and bones, gap 0.2 mm, step 3 mm, 3 rings, 8 rounds a call, 2 passes, each call on the muscle near the partner).
  MEASURED with h2carve.v1's call (repro/sepcase.c); crossing pairs with the five held parts, the muscle's own crossing pairs:

  | header | crossing pairs | own pairs | inside | worst move | volume |
  |---|---|---|---|---|---|
  | input | 748 | 0 | 37 | | 32.80 ml |
  | 1.4.0 | 14256 (the capsule's call 2504 -> 4531: the maintainer's numbers) | 24125 | 39 | 45.8 mm | -504 ml |
  | 1.4.0, 30 rounds a call | 55201 | 236091 | 19 | 138.3 mm | -36757 ml |
  | 1.4.1 = 1.5.0 (8 or 30 rounds) | 581 | 0 | 19 | 5.09 mm | 32.41 ml |
  | 1.5.2 (8 rounds; 30: 645, 10.0 mm) | 633 | 0 | 19 | 8.84 mm | 32.33 ml |
  | 1.5.2, MK_SEPVOL_PAIRS 0 | 693 | 0 | 21 | 0.80 mm | 32.77 ml |

  Cause: step 2 pushes each corner behind a held face's plane to gap in front of it, with no limit on how far behind. A face of
  M wider than the held part is thick has corners behind one side or the other, pushed through it, folding M (1.4.0 alone
  against the clean radial collateral ligament: 102 -> 224 pairs, 234 own). With step 2 off 1.4.0 ran clean (687 pairs, 0.80 mm).
  1.4.1 already bounds it: its guard refuses the folds, its best state comes back, and the capsule and annular ligament, which
  cross themselves (72 and 6 pairs, meshcheck FAULT), are refused (-2). No further fix was needed for the runaway.
- the maintainer's two asks. MK_SEPVOL_PAIRS (step 2 on or off; any int expression, read at each call): for a held part that runs
  through M with no vertex of M inside it (ligaments through a coarse muscle, ribs through the oblique sheet), which no vertex
  move clears; step 2 then squeezes M round it (the rows' cube: a third of its volume; the brachioradialis: 8.8 mm against 0.8).
  h2carve kept its own copy of the held rule for this (pushcut). A total-move bound was not added: no vertex moves more than step
  in a round, so none more than rounds times step in all (now said in the header), in the caller's unit; a second bound would
  depend on those two.
- H1's entry side (#56, race-c1 H1's proposal h1/prop/mk_sepvol.h, the same rule as H2's rayexit) taken as the rule in held mode,
  with one check added: a vertex whose nearest exit lies ahead along its own outward normal leaves along its inward normal to
  where that ray leaves H, unless M's own surface comes first on that ray (H runs through M there). Without the check the thick
  tibia sent tibialis posterior vertices back through 31-41 mm of bone. MEASURED on H1's 13 muscles (b67, gap 0.3 mm, step 1 mm,
  3 rings, 10 rounds, 2 sweeps; repro/h1, h1sep built against each header): inside vertices plus crossing pairs with the partner,
  the measure mk_sepvol keeps its best state by, 1.4.1 5873, H1's patch 4716, 1.5.2 4702; inside 643, 508, 469; volume -8.50,
  -9.69, -9.51 % (the film leaving); own crossing faces 14 -> 14 in all three (meshcheck -p). Per muscle (inside + crossing):
  subscapularis l 722, 548, 565; r 788, 516, 526; supraspinatus l 289, 93, 85; r 337, 116, 111; extensor hallucis longus l 1187,
  1075, 1078; r 1279, 1108, 1113; extensor digitorum longus l 470, 432, 434; r 511, 478, 478; tibialis anterior r 17, 0, 0; l and
  quadratus lumborum 0. Worse: tibialis posterior l 160, 180, 176 and r 113, 170, 136 (the thick tibia). Volume costs: extensor
  hallucis longus -2.8 / -3.2 % on 1.4.1, -14.0 / -11.0 % on 1.5.2, supraspinatus r -5.1 -> -16.5 %: the film now goes.
- mesher-helper's crossings (#57): not the guard, and not the definition. The liver push (diaphragm 3592676, the first liver push)
  and the heart chain to pushes 4872614 and 4242366 were replayed exactly (rounds, counts and volumes as in sep2/*.log): on the
  doubles mk_sepvol returns, 0 own crossing pairs by meshcheck's strict test and by xfaces' (which also counts an edge touching a
  face's edge), and the report's self1 0. Written at %.7f metres as sepheld writes them: the liver 1 pair (faces 1127193 and
  1127297), the heart 2 (665382-665541, 692051-692052) and 11 (16 faces), the same by both tests and by meshcheck: the log's 2, 4
  and 16 faces. The guard lets faces come within a rounding of each other; writing coordinates with %.17g keeps its result (said
  in the header). The other four undone heart pushes were not replayed.
- Rows (t_sepvol.c): 9, a cube of 2.5 mm faces against a plate 0.4 mm thick at the maintainer's settings, 30 rounds (1.4.0: -112 %, 91
  own pairs, 202 -> 610; 1.4.1 bounded, 116 left; 1.5.2 clear); 10, a cube past the middle of a 1 mm plate (1.4.1: 121 of 146
  vertices left inside; 1.5.2 clear in 8 rounds, 1.40 times the overlap); 11, a bar through a slab (the own-surface check: -9 %,
  without it -32 %); 12, the switch (off: nothing moves). Row 6b now comes clear (1.4.1 stalled and handed back the input; its
  "refused > 0" check dropped, a clear check added). meshcheck reads coarse and plate clean, coarse_before and plate_before
  FAULT. The new test fails 6 rows on 1.4.1's header and 16 on 1.4.0's. Planted faults 193 (no held exit), 194 (no own-surface
  check), 195 (the switch read as on), all caught, with 151-159, 191, 192.
- The header's "1.5.1" labels (the 1.4.1 fix's working name) read 1.4.1 now.

Part b, merged by the maintainer:

- mk_dec.h, bug: mk_decimate with clean = 1 wrote before its owner table (heap-buffer-overflow at mk_dec.h:524 under ASan, heap
  corruption without it: the maintainer's subscapularis l patch). In the clean-up's relaxed try a rim ear, every triangle of u's holding
  v, leaves no new triangle to own u and what u owned; the owner stayed -1 and the collapse read and wrote ohead[-1]. Now the
  nearest of v's triangles that stays owns them, and a collapse with no owner at all is refused. Every other path is unchanged:
  t_dec's sphere, torus and hollow-ball OBJs are byte for byte 1.5.0's. Row: a needle ear on a held grid (1.5.0 aborts, exit 134;
  under ASan the same report as the maintainer's). Planted fault 201.
- mk_graft.h, fault: mk_graft_stitch on a patch that was not thinned put an edge in twice (rc 1, dup 2-5, lone 12-30 on t_graft's
  dent; the maintainer's tibialis anterior r 14-20 dup). A dense patch rim's places on S's coarser rim edges step back. The rims are now
  taken as loops. Along a patch loop alone on its loop of S's, the longest run of places falling strictly from a start is kept,
  and each vertex out of it is merged into a kept rim neighbour (rep->merged; no move past snap). Where the places already fall,
  the result is 1.5.0's: t_graft's nine rows are digit for digit 1.5.0's. Unthinned, the patch now stitches closed, but its small
  rim triangles can still cross S's rim triangles (t_graft's dent by meshcheck: 4 crossing pairs, 34 on a coarse S; thinned, 0),
  so thin it where crossings matter. Loops without a partner (a hole in the patch by the seam, a patch island, an S island, a
  patch loop on two of S's) are counted in rep->loops_s, loops_p and unpaired and placed as in 1.5.0. Rows: the unthinned patch on
  S of 5 and 3 subdivisions; a one-triangle hole by the seam. Planted faults 202 and 205.
- mk_graft.h, fault: the stitch kept the most positive piece, so a body that came out broken or inverted, its volume negative,
  lost to a closed speck and passed with rc 0 (the maintainer's supraspinatus r became a 42-face speck: against the humerus S had 3 rim
  loops and the patch 5, and the body piece had dup 35, lone 41, volume -1.47 mL). Now the piece of most volume of either sign
  is kept. rep->inverted and inverted_vol count the pieces wound against S. 2 is returned (out still given) when the kept piece is
  wound against S, or when such pieces hold over 1 % of its volume. On H3's muscles, MEASURED: 1 to 70 % where the body broke,
  under 0.1 % for the rounding specks of a clean stitch. Rows: S wound inward with a speck wound outward (the body kept); a closed
  inward piece of radius 1.5, 0.5 and 0.05 beside the body (rc 2, 2, 0). Planted faults 203, 204 and 210.
- mk_graft.h: mk_graft_refine, the maintainer's seam refinement (roomsplit refine()). Before init it splits S's edges near the seam at
  their midpoints, with no T-junction, and leaves the surface as it was. Row: the coarse S, closed, volume off by 3e-15, every
  triangle away from the seam as it was; grafted after it, crossing pairs 31 -> 3. Planted fault 209.
- mk_geom.h: mk_seg_tri_strict (mk_dec's crossing test, moved with its arithmetic; mk_dec_segtri_ now calls it), mk_tri_cross,
  and mk_bvh_crossings: a mesh's crossing face pairs, all of them or those of marked faces (the maintainer's local crossing count).
  Row: two spheres sunk into each other, against brute force by signed volumes. On 33 of the maintainer's meshes it equals his
  crossings() in 165 of 165 counts. Planted fault 206. mk_sepvol_self_ is the same count on mk_ray_tri's edge test and can call
  mk_bvh_crossings at its next release (not touched here).
- mk_sdf.h: mk_sdf_mesh_box (the maintainer's sdfbox), which takes "a box for mk_sdf_mesh" off the not-yet list. It grids a box's nodes
  only. Every triangle still signs the rows. MK_SDF_FILL is refused. BIG is widened to the box and the mesh together.
  mk_sdf_mesh's grids are bit for bit 1.5.0's (t_sdf rows, t_par BITS). Rows: two boxes against the whole grid (with
  MK_SDF_EXACT every shared node bit for bit; plain, every sign and every band value), a box away from the mesh within the
  sag, and 5 refusals. Planted faults 207 and 208.
- Not taken: sdfbox's sign rule (classvox's winding along +z, inside where it is not 0). On a closed mesh mk_sdf's parity with
  the winding number agrees. Covering a region twice is MK_SDF_FILL's case, and a box refuses it.
- Rows: dec 13 (+3 meshcheck), graft 15 (+3), geom 21, sdf 16; new planted faults 201-210, all caught.
- README table: mk_dec 12 -> 13, mk_graft 10 -> 15 (and mk_graft_refine), mk_geom 20 -> 21 (and mk_seg_tri_strict,
  mk_tri_cross, mk_bvh_crossings), mk_sdf 14 -> 16 (and mk_sdf_mesh_box). Not yet: drop "a box for mk_sdf_mesh". Planted-fault list:
  dec 201; graft 202-205, 209, 210; geom 206; sdf 207-208.

### 1.5.0 (2026-11-16, installed 07:39)

OpenMP threads in the hot loops (the maintainer's device plan, scratch/the maintainer/devflex/PLAN.md item 7; the maintainer's helper). A minor
version: the kit's interface grows (mk_par.h, three batch calls) and two solvers' last bits move once. Staged from 1.4.0; 1.4.1's
mk_sepvol fix taken whole (mk_sepvol.h and t_sepvol.c as installed, check.sh's sepvol spec lines merged by hand).
- Profile first (scratch/mathkit-150-work/PROFILE.md). The grown tools' 9.25 CPU-h are grown's own code: they vendor the kit
  but call only mk_sdf_mesh (7-13 % of `grow seed`) and never mk_sparse; no house program calls mk_sparse. The kit's real
  time is mk_sdf_mesh's sweeps (92-98 % on closed meshes), its odd rows' winding numbers (55 % on the body mesh's hole), its
  exact pass (81 %), and the BVH queries that partjudge and bonegraft make one point at a time (84-89 % of their samples).
- mk_par.h (new): MK_OMP_ (pragmas only under -fopenmp: a caller without it starts no thread and gets the same bits), fixed
  blocks (mk_par_blk: 4096 items, at most 256 blocks, from the count alone), MK_PAR_FOR_, MK_PAR_SUM_ (partials added in
  block order). No reduction(+:) on floating point anywhere.
- mk_sparse.h: the product, dots, scaling, mean-free shift, vector updates, Jacobi and block-Jacobi on blocks of rows; CG
  fuses each product with p.q, each update with r.r and Jacobi with r.z (three or four passes an iteration, was six or
  seven); IC(0)'s solves stay serial. Bits: mk_csr_matvec unchanged; mk_cg and mk_bicgstab change once against 1.4.0 when
  n > 4096 (more than one block); a vector of at most 4096 entries keeps 1.4.0's bits.
- mk_sdf.h: band by 8 x 8 tiles of rows (each tile's triangles in index order), sweeps as a wavefront of 4 x 4 tiles, the
  exact pass, the odd rows' winding numbers and the near nodes node by node, the rows' parity row by row. 1.4.0's grid bit for
  bit at every thread count (MEASURED below).
- mk_geom.h: mk_bvh_nearest_many, mk_bvh_inside_many, mk_winding_number_many: batches on threads, each answer the single
  call's bit for bit. A consumer opts in by passing its points as one array.
- Against 1.4.0 (MEASURED; scratch/mathkit-150-work/cmp, built from mathkit/include at 03:4x, then 1.4.0; 1.4.1 has the same
  mk_sparse.h, mk_sdf.h and mk_geom.h):

  | function | 1.5.0 against 1.4.0 (no OpenMP, and OpenMP at 1, 2 and 6 or 16 threads: all one result) |
  |---|---|
  | mk_csr_matvec | identical (48^3 Laplacian, 110592 rows) |
  | mk_cg none, Jacobi, block-Jacobi, IC(0), mean-free (48^3) | same iterations (209, 209, 198, 72, 158); x differs in the last bits at 88-96 % of entries, by at most 0.92e-15 to 2.25e-15 of the largest entry |
  | mk_bicgstab none, Jacobi (48^3 convection-diffusion) | same iterations (70); x within 1.75e-15 and 2.14e-15 of the largest entry |
  | mk_sdf_mesh | identical grids and counts: femur 0.5 mm plain and MK_SDF_EXACT (54.9 M nodes), external oblique 1 mm, and 2 mm with MK_SDF_FILL and MK_SDF_LATTICE, the body 5 mm (its odd row), 8 mm MK_SDF_EXACT, 8 mm FILL and LATTICE (92 filled); at 1, 2, 6 or 16 threads and without OpenMP |
  | every other function | unchanged code |

- Rows: t_par (new) 21 with OpenMP (16 thread rows at 1, 2 and 12 threads, the batches against single calls, the sdf cases'
  reach, CG's block sums against one serial sum: 1e-15 of the largest entry, 108 iterations both), 5 without; check.sh's
  BITS comparison (16 results of a build without OpenMP equal the OpenMP build's); every header compiles alone with
  -fopenmp too. Planted faults 96, 97, 98, 107, 108, 146, 147, 148, each making a result follow the thread count: 92, all
  caught (with 1.4.1's 156-159, 191, 192). test/negctl.sh: 5 negative controls for the rows no header fault reaches, all failing as they must. `make check
  OMP=1` builds every test and fault with -fopenmp. Both gave wrong: 0, 92 planted faults caught, on the tree with 1.4.1's fix
  merged (the development machine, `jobwatch run -t 30m`, 7 minutes each; before the merge also on a second machine).
- The house flags: t_par's BITS at C23_REF, C23_FAST and C23_FAST -fopenmp are one set; t_sparse, t_sdf, t_geom, t_sepvol,
  t_ovl, t_graft, t_dec, t_snet pass at C23_FAST -fopenmp (scratch/mathkit-150-work/flags.sh).
- Speed: below ("Threads, MEASURED").
- Proofs: the verification work's QcMkSp6 (mk_sparse.h, from freeze 1.0.6) is touched: the CG, BiCGSTAB, dot, product, scaling,
  mean-free and preconditioner statements change shape (block functions, fused passes), and gen/clight_sp.sh copies only a
  header's `#include <...>` lines, so the unit must also take mk_par.h (without -fopenmp MK_OMP_ is empty and the block
  loops are plain C). No Lean unit covers mk_geom.h or mk_sdf.h; QcMk, QcMkLm, QcMkLs, QcMkOq, QcMkOs, QcMkSrv read
  headers 1.5.0 does not change.

### 1.4.1 (2026-11-16)

mk_sepvol no longer makes the mesh it moves cross itself (mesher-helper the maintainer.md #50: the heart pushed out of one rib gained 49
crossing pairs; wall parts 323 to 30,477 faces; the liver 305; fa75's tubes in single mode grew many-fold). Cause: on meshes
crumpled at their edge scale, each vertex's own nearest point and own noisy normal sent neighbours half an edge in different
directions, folding the crumples through each other. The fix is in mk_sepvol.h's header comment; in short:
- the round's moves are tried together. A face that would cross one of M's faces it didn't cross at the round's start has its
  moving corners' moves halved up to MK_SEPVOL_HALVE (4) times, then refused (single mode guards within a part only). The test
  is strict, as meshcheck's;
- each push turns toward its pushed neighbours' mean direction, MK_SEPVOL_SMOOTH (2) passes, keeping its own length;
- the best measured state is returned (fewest inside vertices plus crossing pairs), and the run stops after MK_SEPVOL_STALL (8)
  rounds without a better one;
- a held mesh that is open, non-manifold, inward facing or crossing itself gives -2 with V unchanged. In single mode, a part
  that already crosses itself is reported (status bit 8) and never made worse.
- Report fields appended: status, self0, self1, refused, shortened, held_open, held_flipped, held_cross, held_vol. `moved` now
  counts vertices whose returned position differs from the input. The signature is unchanged (sepvol_obj and sepheld.c compile).
- Results shift slightly: twoball -1.20 % (was -1.19), ballB -3.84 % (was -3.63), the cubes lose 1.30x their overlap (was
  1.20x). One 1.3.0 row is split, decided by the maintainer: smoothing turns 4 of the cubes' 11 edge vertices (those within 2 rings of
  the other cube's sides) up to 33.8 degrees off -x. So all 11 must lie nearer -x than their 45-degree normal, and the 5 farther
  than 6.38 mm from those sides within 10 degrees (worst 0.0).
- New rows: a crumpled cube against a 45-degree bar (1.4.0: own pairs 0 to 7, meshcheck FAULT; now 0); the same against a thin
  bar 0.8 mm in (1.4.0: 0 to 3741, volume -18 %; now the input back); the held mesh open, doubled, inside out, one face turned (all
  -2); a self-crossing ball in single mode (1.4.0: 12 to 166; now 12). The new test file fails 9 rows on 1.4.0's headers.
  Planted faults 156-159, 191 and 192, all caught.
- MEASURED on the real cases (PRIVATE meshes, not in the tests): heart against rib 110061, muscle 5285189, and the liver against
  muscle 5320787 are separated clean (0 own crossings, 0 inside, meshcheck clean). Heart against bones 29567 and 234277, and the
  liver against muscle 5248925, are much improved but not fully separated (deep intrusions stall; more rings or splitting the
  part would go further). The wall union and a self-crossing held muscle give -2. Heart against rib went from 0.4 s to 7 s.

### 1.4.0 (2026-11-16)

Three mesh methods from mesher-fable, offered by mesher (the maintainer.md #43), each merged with the house's other copies (SOURCES.md).
A minor version: the kit grows.
- mk_dec.h, mk_decimate: mesher-fable's folddec (from mesher's objdec), with experiment/cfd's remesh and body-claude's
  skindec read beside it. Quadric half-edge collapses that keep the input's vertices, refused when they would add a fold, a
  crossing or drift past a two-sided bound, or break the edge, angle, turn or link rules. Its millimetre constants now scale
  with the caller's unit; mk_dec_defaults(o, mm) gives mesher's values. Byte for byte folddec's output on an icosphere, a
  torus and b-helper's latissimus dorsi r (1812456 -> 37142 vertices, mesher-fable's number). Speed, one thread, perf (MEASURED): the
  latissimus on a second machine at load 31-35 on 16 threads, 836 G cycles and 165 s against folddec's 1157 G and 227 s (the crossing grid's
  cell is at most 4 mean edges, folddec's a fixed 2 mm); icospheres on the development machine, 164 k vertices 90 G against 125 G cycles, 41 k 18.8 G
  against 15.1 G. mk_ovl on two icospheres of 164 k vertices: 10.4-10.7 G cycles against ovl.c's 10.5 G.
- mk_ovl.h, mk_overlap_volume and mk_ovl_bvh: mesher-fable's ovl (columns, odd lines counted) with stepcad's voldiff's
  differences (only_a, only_b, each volume, MK_OVL_ALL) and choice of axis. The same overlap and line count as ovl.c on a
  sphere pair; the lens within 0.08 % of analytic. mesher's census was read as the reference, not merged.
- mk_graft.h: mesher-fable's bonegraft's general parts. mk_mtet (marching tetrahedra, closed whatever the field) on mk_sdf's
  grid layout, and the graft in three calls (init with the seam function, patch, stitch) so a caller can thin the patch with
  mk_decimate between, as bonegraft's two steps do. The bone-specific field stays with its owner; the test builds its field
  with mk_sdf_mesh. body-claude's tet.c (isosurface stuffing) read, not merged: a volume mesher.
- mk_decimate's report is zeroed on a refusal, as mk_sepvol's.
- Rows: dec 12 (+3 meshcheck), ovl 9, graft 10 (+3 meshcheck). Planted faults 161-165, 171-173, 181-182: 78, all caught.
- Signatures: decimate, mtet, overlap_volume (positive controls in signatures.tsv). Other house copies of marching tetrahedra
  (body-claude's seatcov.c, thighvox.c, wallimp, wallband) and of quadric decimation (atlasbake.c, remesh.c, fvremesh.c) move at
  their owners' next release.
- Done in 1.5.2: mk_sdf_mesh_box, a local field on a large mesh with no whole-mesh grid.

### 1.3.0 (2026-11-15)

mk_sepvol.h, from cad-helper's sepvol (the maintainer.md #41; the maintainer-fable's proposal, C1 films), merged with the maintainer's selfsep;
and the verification work's patch for mk_geom, mk_sdf, mk_snet, mk_vpflow and mk_ode (Results 86 and 87, the maintainer.md #42), applied first,
since mk_sepvol builds on mk_geom. A minor version: the kit grows.
- mk_sepvol (SOURCES.md): one mesh's parts sunk into each other, or a part out of a held one, moved apart with the vertices,
  their order and the faces kept. Edge-aware: a vertex inside goes toward the nearest point of the other surface, not along
  its own normal; crossing and near face pairs by selfsep's rule; nothing moves outward. Two cubes 1 mm into each other, which
  sepvol left at 182 crossing pairs, end clean in meshcheck in 6 rounds; twoball and ballB out of a held ballA too.
- mk_geom gains mk_bvh_nearest_lab (the nearest point among the triangles whose label is not one left out) and mk_bvh_box
  (the triangles whose boxes meet a box), both for mk_sepvol; mk_bvh_nearest is mk_bvh_nearest_lab with no labels.
- the verification work's patch, read before applying, as written but for two speed changes below: F1 the open interval at tmax in
  first_hit; F2 the SAH binning safe on NaN and overflow (1.2.2 indexed bb[INT_MIN] on a NaN, an infinite or a 1e308 vertex:
  UBSan, MEASURED with the verification work's demo_nan); F3 SAH to depth 200, then median splits, so the 256-entry stack never drops a
  subtree (6006 triangles 873 deep: 473 hit counts and 1417 nearest points wrong on 1.2.2); F4 inside: a ray with a hit on an
  edge or a vertex does not vote, the winding number decides when the clean rays don't agree; F5 the nearest-point prune
  widened by the rounding a distance may carry, so a tie goes to the lower index; F6 counts past (INT_MAX - 1) / 2 refused
  before 2 nt + 1 overflows (gcc warns of it on 1.2.2's code at compile time; the row can't fail there); F7, F10 mk_ray_tri's
  edge functions by Kahan's a b - c d (fma), F12 the box test's pad from the largest coordinate; F8a, F9 the closest point's
  edges as candidates on a sliver and a repeated corner. V-A, V-B mk_vp: INT_MAX - 3 nodes refused, a refused grid gives NaN
  from eval, vel and div. S-A to S-D mk_sdf: a lattice past the double range refused, an empty grid samples NaN, no qsort of
  NULL, every refusal leaves the grid empty. N-A an empty net's vertex is -1. O-E' the step cap compared before the count
  in mk_dopri and mk_ode23s.
- Our calls on the verification work's open notes: N-B `(void)ci;` in place of `[[maybe_unused]]` at mk_snet.h:170, so clightgen takes
  mk_snet_build_; S-E mk_sdf.h's header says distances are float, inf beyond FLT_MAX. Known limits, not fixed now (the verification work's
  NOTES.md in scratch/the verification work/mk-geom and scratch/the verification work/mk-new): F8b, the face weights thousands of ulps off on a
  triangle of aspect about 1e5 (about 1e-12 of its size: far under any gap mk_sepvol uses); F11, "no ray slips between two
  triangles sharing an edge" holds only while the ray's sheared frame resolves them (a ray from 2e11 away slipped through
  an edge 0.06 long); S-F, mk_sdf's sliver precision is F8b's.
- Speed changes to the patch (bits unchanged against the verification work's version wherever its edges did not win by rounding):
  the closest point's edge pass runs only on a sliver (|ab x ac|^2 under 2^-20 |ab|^2 |ac|^2), a repeated corner or a weight
  below 0 or NaN, out of line (MK_COLD_, GNU noinline and cold); inlined it took mk_closest_pt_tri past gcc's inline limit
  and cost mk_sdf_mesh 30 to 47 % at -O3 (MEASURED on the development machine, rest-uv at 10 mm). The gate still adds 20 % instructions
  there (24.1 G to 28.9 G, perf; 184 M closest points, none a sliver): within noise in cycles on the development machine (17.7-21.4 G against
  17.7-19.6 G, three runs each at load 12), 15 to 25 % slower on a second machine at load 30. First hits are 10 to 25 % slower on a second machine
  (Kahan's products, F7 and F10). UNVERIFIED on a quiet machine. The prune's slack is computed once a query from the root's box (it bounds every
  node's, so the slack is never smaller than the verification work's) and again only when the best distance changes.
- Bits on fullset-7.0.0's rest-uv.obj (1774083 faces, h 10 mm, 931392 nodes; 200000 nearest points and first hits, 20000
  inside tests; MEASURED, scratch/mathkit-130-work/cmp): mk_sdf_mesh bit-identical to 1.2.2 (twoball: one node of 54288
  differs, by 6e-21 m, no sign changes); nearest points bit-identical; first-hit triangles identical, t within 3.5e-15
  relative (from Kahan's products); inside identical. F3 and F4 change none of them.
- Rows: one for each finding (mk_geom 12, mk_sdf 4, mk_snet 1, mk_vpflow 2, mk_ode 1; F2, S-B and N-A skipped under T_OLD,
  since 1.2.2 crashes there). On a 1.2.2 copy: F1, F3, F4, F5, F7, F8a, F9, F10, F12, S-A, S-D, V-A and V-B FAIL (13); F6,
  S-C and O-E' pass there (undefined behaviour no test can observe without UBSan: an int overflow, qsort of NULL with 0,
  LONG_MAX + 1 on the LONG_MAX-th step). Planted faults 95 (the box query) and 151-155 (mk_sepvol): 68, all caught.
- The table's row counts for mk_vec, mk_rng and mk_stats were stale since 1.2.2 (9, 8, 11; they are 11, 10, 14).

### 1.2.2 (2026-11-15)

the verification work's stats, rng and vec fuzzing (Result 84), its patch as written, and O-E from Result 85:
- S1: mk_quantile_sorted gave NaN beside an infinity ({1, 2, inf} at q 0.5; {-inf, 1} at q 0). Order statistics, equal
  neighbours and one infinite neighbour are now exact. S2: mk_select refuses n <= 0 and k outside [0, n) with NaN.
- R1: mk_rotl64 by 0 shifted by 64 (undefined); now k & 63. Its row passes on 1.2.1 as well, since x86 happens to give the
  right answer; the proof, not the row, is what catches it. R2: mk_rng_below(g, 0) returns 0 instead of dividing by zero
  (its row runs on 1.2.2 only: on 1.2.1 it traps).
- V1, V2: mk_v3_norm and mk_q_norm scale by a power of two first (same bits wherever |v|^2 did not under- or overflow):
  (3, 4, 0) x 2^e normalizes right from e -1070 to 1020, and a quaternion near 1e-170 keeps its direction.
- O-E: mk_dopri and mk_ode23s test the step cap before incrementing, so maxsteps = LONG_MAX cannot overflow the count.
mk_rng is proved bit-equal to the published xoshiro256** and splitmix64 (the verification work, Result 84).

### 1.2.1 (2026-11-15)

the verification work's F11 (Result 83): 1.0.6 scaled b but not A, and the solution updates formed a * p before the multiply by b's
scale, so near 2^-1015 a * p overflowed: CG returned -1 on A = 2^-1015 tridiag(-1, 2.5, -1) with x0 -1000, where every
scale to 2^-1012 solves. The updates in mk_cg and mk_bicgstab now take (a * sc) * p. sc is a power of two, so results are
bit-identical wherever nothing over- or underflowed. Rows: the F11 system (1.2.0 FAILs: ret -1, x inf) and b = 1e308 against
the same A (an unrepresentable answer is refused), the max|b| >= 2^1022 class the verification work's generators never reached.

### 1.2.0 (2026-11-15)

mk_sdf.h, a triangle mesh to a signed distance grid, s-bone's meshsdf (Batty's SDFGen, the maintainer.md #37) merged with
body-claude's sdfgrid.c, experiment/cfd's levelset.h and metal's hairgen.c (SOURCES.md). A minor version: the kit grows.
- The distance is meshsdf's: exact in a two-cell band through mk_geom's mk_closest_pt_tri, then two rounds of eight sweeps.
  MK_SDF_EXACT makes every node exact through mk_geom's BVH (on the tests the sweeps are within 0.05 h of it beyond the band
  and equal within it).
- The sign is meshsdf's row parity, with what it lacked. A row with an odd count is cast again on two nearby rays (a ray
  through an edge), and if all three are odd (a hole) its nodes take the generalized winding number. A node nearer the
  mesh than the ray's offset gets a ray of its own. The winding number was considered for every node and kept for these:
  it is O(triangles) a node (38 ms on the 1.77 M-face body), where parity is O(crossings), and on a closed mesh they agree
  (the nested shells, the random mesh). On fullset-7.0.0's rest-uv.obj, which has a 97-edge hole, meshsdf's parity puts
  117 nodes on the wrong side at 5 mm (a stretch of one row inside the body reads outside, winding 0.99 there); the kit
  signs them by the winding number.
- MK_SDF_FILL is sdfgrid's flood; MK_SDF_LATTICE is hairgen's world-lattice placement. Grid sizes are refused when an axis
  or the node count would overflow, before anything is allocated.
- For s-bone, a wrapper with meshsdf's call (mesh_sdf, 1 or 0, g->v malloc'd) is in SOURCES.md's mk_sdf.h section.
Planted faults 141-145.

Two new headers that other sessions wrote: a minor version, since the kit grows. Both were read before use and brought
to the kit's conventions; SOURCES.md says what came from where.
1. mk_vpflow.h, b-helper's (the maintainer.md #36): a volume-preserving deformation, the flow of u = curl A with A in cubic
   B-splines, any axis periodic. Its RK4 is now mk_ode.h's mk_rk4_step, bit-identical to its own (one RK4 in the kit);
   the input checks cost 2-5 % (351 against 343 ns a point-step at -O3 -march=native). Added: mk_vp_ok, NaN for a
   non-finite coordinate, h or empty axis, a periodic coordinate past 2^30 cells reduced before it becomes an int, and
   mk_vp_advect's NaN for a refused grid, steps < 1 or T not finite. Its claims, checked here: div u 1.1e-10 of |u| / h
   (finite differences at 1e-5), the ball's volume 1.8e-6 off (bar tightened from 2e-3 to 2e-5), the Jacobian 1 within
   2e-7, 40 steps against 80 6.2e-10, the periodic axis. Added rows: a rigid rotation from a quadratic potential (exact
   field to 1e-13, a quarter turn of order 4.00 over 5 to 40 steps) and a uniform flow across the periodic seam. The
   API is unchanged for fmap-v22.c (it compiles against the kit header as it is; it needs gnu23 or _DEFAULT_SOURCE for
   its own M_PI). Planted faults 131 and 132 (were 1 and 2).
2. mk_snet.h, s-forms' through the maintainer (the maintainer.md #35): surface nets on a labelled grid, shared interfaces, Taubin
   relaxation held in each dual cell, outward winding. Merged with the house's copies: memory that goes with the surface
   (wrap.c's layers; 31 MB against 275 MB on a 400^3 grid, build 2.5 times faster), a scalar-field form with mesh.c's
   held iso crossings (mk_snet_build_iso, mk_snet_faces_iso), and wrap.c's two-group test as mk_snet_split. Results
   equal s-forms' bit for bit. Every surface goes through mesher's meshcheck in check.sh. Its own faults 121-123
   kept; 124 and 125 added. For s-forms: s.vid is gone (mk_snet_vertex(s, i, j, k) finds a cell's vertex);
   mk_snet_build, mk_snet_faces, s.V, s.nv and mk_snet_free are as they were, and snetobj.c compiles against the kit
   header unchanged.

### 1.0.6 (2026-11-15)

the verification work's mk_sparse.h fuzzing and proofs (Result 80): its patch, read and applied, with one change (F9 below). One row
per finding; 9 of the 11 that run on 1.0.5 fail there (F9's second row and F10 pass there: those faults came with the
patch's own scaling), and 3 are skipped there with T_OLD because they crash 1.0.5 (F1 aborts in malloc, F6 and nb past
INT_MAX / 3 segfault).
1. F1: preconditioner kinds outside 0..3 are refused; kind 4 wrote 9 doubles per block into an n-double buffer.
2. F2: CG and BiCGSTAB work in units of 2^e, e the exponent of b's largest entry, so |b|^2 neither under- nor
   overflows: b = 1e-170 came back "converged" as x0, b = 1e170 gave NaN and -1, A = 1e-87 I with b = 1e-121 gave -1.
   Results are bit-identical wherever nothing under- or overflowed: 396 of 400 solves (both solvers, every
   preconditioner, A from 2^-200 to 2^200, b from 1e-100 to 1e100, x0 zero and not) match 1.0.5 bit for bit, and the 4
   that differ are BiCGSTAB solves where 1.0.5's dots underflowed or overflowed and it returned -1; they now converge
   (MEASURED, scratch/mathkit-106-work/cmp.c). A is not rescaled.
3. F3: block-Jacobi inverts each block scaled by a power of two to a largest entry in [0.5, 1), so its 1e-300 test is
   relative: 1e-101 I and 1e110 I lost their inverse and CG returned -1. An exactly singular block still gets a zero
   inverse at scales 1, 2^-700 and 2^700 (a block singular only to rounding passes the test here as in 1.0.5).
   cudakit's k_dinv (spd3.cu) keeps the absolute test.
4. F5-F8: mk_csr_from_triplets refuses n < 0 and nnz < 0 (nnz -2^60 wrapped to malloc(0)); mk_csr_from_bsr3_lower
   resets *A first and refuses nb < 0 and nb > INT_MAX / 3; the work sizes of CG (4 n), BiCGSTAB (8 n), block-Jacobi
   (3 n, 9 b) and bsr3 are formed in size_t. The ck_spd3 layout is unchanged: the bsr3 rows still expand it exactly and
   match a dense LU.
5. F9: CG at tol 0 runs past attainable accuracy until r z underflows. the verification work stopped it whenever r z < DBL_MIN; that
   also ended solves where a large preconditioner makes r z small from the start (A near 1e290 with Jacobi, b near 1:
   1.0.5 converges in 96 iterations, the patch stopped at a residual of 1.9e-10 with status 1; near 1e300 at 2e-5).
   CG now stops on r z <= 0 or NaN, or on an underflow once r z is 2^-100 below its start (the residual about 1e-15 of
   its first). The tol 0 run ends with status 1 and a finite x at true residual 1.1e-14, at b 1, 2^837 and 2^-900 in
   the same 385 iterations; A near 1e290 converges again. Its moderate-scale results are the verification work's bit for bit.
6. F10: both solvers return -1 unless every x[i] is finite (an answer past the largest double came back with status 0
   or 1).
New planted faults 104 (the finite check removed), 105 (the units removed), 106 (the block test absolute again).

### 1.0.5 (2026-11-15)

the verification work's interp, quad and ode fuzzing and proofs (Result 79): its patch, read and applied, with one change (I-B below),
and the calls it left to us. A tridiagonal solve at metal's request (the maintainer.md #31). One row per finding; 8 of the 8
that run on 1.0.4 fail there, and 3 are skipped there with T_OLD (I-C and O-B are undefined behaviour on 1.0.4, the
Gauss-Legendre bound does not compile against it at -Werror).
1. mk_interp (I-A to I-D): every function now takes a one-point table (y0, slope 0) and gives NaN for n < 1, where it
   read past the table. The interior PCHIP slope compares signs: d0 * d1 underflowed for secants near 1e-162 and below,
   so y = 1e-170 x gave 6.25e-171 at 0.5 for 5e-171. Results stay bit-identical to aeromap.c except there. Change to
   the verification work's patch: its sign test also turned a NaN secant into slope 0, hiding a NaN neighbour; it gives NaN again, as
   aeromap.c and 1.0.4 do. mk_locate's midpoint is lo + (hi - lo) / 2 (n > 2^30 + 1 overflowed), and the spline's work
   size is formed in size_t (4 * n wrapped to 0 at n = 2^30, then the first loop wrote past the block).
2. mk_quad_adapt (Q-A): the running sums drifted by cancellation, so status 0 came with an error up to 3.9e5 times the
   limit (6 spikes on a floor; 0 now). The sums are re-formed on each pass. It refuses maxint >= 2^29 (status -1, NaN)
   before allocating or calling f, and mk_gauss_legendre refuses n > 2^30, so 4 i + 3, 2 k - 1 and (n + 1) / 2 stay
   inside an int: the bounds the verification work's proofs need (1 <= n <= 2^30, maxint < 2^29).
3. mk_ode (O-A, O-B, O-D): the last step takes the remainder when it is within 1e-12 of the span (hmax 1/k on y' = 1
   failed for 21 of 49 k); mk_dopri refuses n < 0 (memcpy with a negative size) and ignores a NaN or infinite *h (fmin
   dropped it, and the first step was the whole span, past hmax).
4. O-C, our call: a NaN or infinite t0 or t1 gets -1 from mk_dopri and mk_ode23s before f is called (1.0.4: 0 after one
   call, nothing integrated); mk_bdf_step gives -1 when t + h is not finite. mk_ode23s and mk_bdf_step also refuse n
   outside 0..46340, the sizes mk_lu refuses, before allocating (1.0.4 asked for the n^2 block first).
5. mk_tridiag (linalg): metal's world.c and roadline.c each carry the same Thomas elimination; mk_interp's spline a third
   (c', d' form, left in place: the verification work proved it). mk_tridiag keeps metal's order, so their answers are bit-identical (0
   of 4100 values differ), and adds the refusals neither copy had: a pivot that is NaN, zero or under 2^-52 of its row
   (singular to rounding, or the system needs pivoting), a non-finite answer, n < 1. No cyclic form: neither copy needs
   one. Planted faults 26 (the pivot left uneliminated) and 27 (the small-pivot test reduced to exact zero).

### 1.0.4 (2026-11-14)

the verification work's L-BFGS finding H (Result 78), its patch as written: y.y in gamma = s.y / y.y underflowed once |y| < 1e-154 (gamma Inf,
status -1), and the pair guard's sqrt(ss yy) under- or overflowed past |s||y| near 1e154, rejecting every pair (the run
became steepest descent). Both now scale by powers of two (exact; results bit-identical wherever nothing under- or
overflowed, the 32 MGH problems included), and a pair needs s.y >= DBL_MIN, so rho cannot overflow. New row: |x|^4 with
gtol = ftol = 1e-300 reaches f 3.0e-301 (1.0.3: status -1 at 1.02e-217). The |g| near 1e161 row now accepts status 1 (the
ftol stop) as the other scales give: its old status 0 with x exactly 0 came from every pair being rejected. mk_lm's
guard comment corrected: with size_t sizes the bound only keeps the block under 2^34 bytes (Result 77).

### 1.0.3 (2026-11-14)

the verification work's mk_lm proof (Result 74): 3 * m was an int, so m >= 715827883 wrapped it, malloc got 16 bytes and the first res()
call wrote past them (ASan and UBSan confirmed). The size is now formed in size_t, and m > INT_MAX / (3 + n) is refused
before anything is allocated or called. Its row runs on 1.0.3 only: on 1.0.2 it writes past the heap block.

### 1.0.2 (2026-11-14)

the verification work's L-BFGS and CMA-ES fuzzing (Result 71) and two size findings from its proofs (Results 70, 73), landed by the maintainer.
Every new row fails on 1.0.1 (7 of 7 run there) and passes now.
1. mk_lbfgs: a rejected (s, y) pair overwrote the oldest pair's s and y but not its rho, so directions mixed pairs; now a
   rejected pair leaves the kept ones alone. A NaN gradient reported status 0 ("converged"); now it stays large and the
   run fails (-1). A gradient over about 1.3e154 squared to Inf and stopped the run at its start; the steepest-descent
   step now uses a scaled norm (|g| near 1e161 converges). n < 0 is refused, and the work size is formed in size_t (2 m
   + 5 overflowed int from m = 2^30 - 2).
2. mk_cmaes: a NaN value at the start was never replaced; NaN broke the ranking. NaN now ranks last, and a run with a NaN
   region is bit-identical to the same run with +Inf there.
3. mk_chol and mk_chol_solve refuse n < 0 and n > 46340, as mk_lu does; mk_lstsq refuses m > INT_MAX / n (i * n + j
   must fit an int).
The rejected-pair fix has no row of its own: no small input shows the mixed directions as a wrong answer.

### 1.0.1 (2026-11-14)

the verification work's second patch (findings b, Results 67-68), landed by the maintainer with regression rows:
1. mk_sym_jacobi returned -1 on finite matrices of subnormal entries: rotations at that scale round to multiples of
   2^-1074, so the relative stop test never passed. It now scales a once by a power of two (exact) and scales the values
   back. Rows: a subnormal 3x3 and 17x17; both fail on 1.0.0.
2. mk_fft_plan accepted n = 2^30, where mk_fft_ri's int doubling overflows. It now refuses n > 2^29 before allocating.
3. mk_lu and mk_sym_jacobi did not refuse n < 0, and n * n overflows past 46340. Both now refuse those sizes. The n = -1
   rows fail on 1.0.0; the n = 46341 and 2^30 rows run only on 1.0.1, since on 1.0.0 they read past the buffer.

### 1.0.0 (2026-11-14)

The first numbered release. It lands the verification work's fuzz findings (the verification work, scratch/the verification work/mathkit/FINDINGS.md, Result
66). the verification work fuzzed against long-double references and wrote the patches, and the maintainer landed them. Each finding has a
regression row built from the verification work's failing input. Every row fails on the text before 1.0.0 and passes now.
1. mk_lstsq decided rank without column pivoting. It kept a column 2.3e-14 from the others' span and dropped one 1.1e-5
   away, so x reached 3e12 and the returned residual was not x's own. Now Householder QR with column pivoting; the rows
   use the verification work's two inputs (test/data).
2. mk_lstsq dropped a column whose squared norm overflowed or underflowed (diag(1e155, 1, 1, 1) gave x0 = 0). Each
   column is now scaled by its largest entry first.
3. mk_sym_jacobi and mk_sym3 returned 0, 0 for [[0, 1e155], [1e155, 0]], had an absolute 1e-300 skip, and accepted
   Inf. The stop test is now scaled, and a non-finite entry gives -1.
4. mk_lm stopped while lambda still damped the weak direction: a kappa 1e7 linear fit stopped short of its minimum.
   It no longer stops while lambda > 1e-9.
5. mk_lm scaled by the current diag(JtJ) and missed Brown-Dennis and Penalty I of the 32 Moré-Garbow-Hillstrom
   problems. It now keeps MINPACK's running maximum: 32 of 32 by the verification work's table, both with a Jacobian and by
   differences.
6. mk_lu_det's running product overflowed: Inf, or 0 in reverse order, where the determinant is 1. Mantissa and
   exponent are now kept apart (frexp/ldexp).
7. mk_chol accepted an infinite diagonal.
8. mk_fft3_c transformed 45 of 48 values before refusing a size; it now refuses before touching the data.

Contract changes, from the verification work's notes, decided by the maintainer:
- mk_lstsq returns -1 for NaN or Inf in A. It used to drop that column and return a finite residual.
- mk_lu_solve and mk_chol_solve now return int: 0, or -1 when the answer is not finite. mk_lu stays a relative test (a
  pivot below 1e-13 of the largest entry). So 1e-310 I, which is well conditioned, is still accepted, and the solve
  reports the overflow instead of returning Inf quietly. Callers that ignored the old void return are unchanged.
- mk_nm sets *evals to 0 when n is out of range.

One fault was found and left alone: a forward-difference Jacobian is exactly 0 when the residual's rounding exceeds
|J| h, and mk_lm then returns the start. the verification work noted it; give mk_lm a Jacobian there.

### Before 1.0.0 (unnumbered, 2026-11-13)

- Stage 1: vec, fft, linalg, rng, opt, ode, stats, interp.
- Stage 2: sparse, quad, spec, geom, L-BFGS, CMA-ES, ode23s, BDF.
- Stage 3:
  - the maintainer's cudakit cross-check found two faults. mk_sum turned an infinite addend into NaN, and now gives a plain
    sum's IEEE result. Sequential mk_welford lost 1.8e-10 near a mean of 1e6; it now accumulates data shifted by its
    first value, and mk_welford_merge (Chan's) was added.
  - mk_rfft takes one sin/cos a bin instead of two.

Totals at 1.0.0: 12 headers plus mk_version.h; 126 test rows (15 new: 12 regression, 3 contract); 45 planted faults, all
caught; `make check` wrong 0 at -O2 and at -O3 -march=native.
