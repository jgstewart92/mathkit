# mathkit

mathkit is a header-only C23 kit of numerical and geometric methods: small dense linear algebra, sparse Krylov
solvers, FFTs and spectral analysis, quadrature, ODE steppers, optimizers, statistics, interpolation, random numbers,
and triangle-mesh geometry (BVH queries, signed distance grids, surface nets). It needs only libc and libm. Built with
`-fopenmp`, its hot loops run on OpenMP threads; without it they run serially and give the same results. This is
version 1.6.4 (`MK_VERSION` in `include/mk_version.h`).

## Using it

Put `include/` on the include path. Every function is `static inline`, so there is nothing to build; link with `-lm`.
Each header includes the others it needs.

```sh
cc -std=c23 -ffp-contract=off -O2 -I mathkit/include prog.c -lm            # serial
cc -std=c23 -ffp-contract=off -O2 -fopenmp -I mathkit/include prog.c -lm   # OpenMP threads
```

```c
#include "mk_sparse.h"

mk_csr A;
if (mk_csr_from_triplets(&A, n, nnz, row, col, val)) return -1; /* repeated (row, col) entries are summed */
int iters;
double relres;
int r = mk_cg(&A, b, x, MK_PC_IC0, 1e-10, 1000, &iters, &relres); /* x: first guess in, answer out */
mk_csr_free(&A);                                                   /* r: 0 converged, 1 not within maxit, -1 failed */
```

- `-std=c23`: the headers are C23; they also compile as C11 and C17 (the tests check C11).
- `-ffp-contract=off` is required. It stops the compiler from fusing `a * b + c` into one FMA, which changes rounding
  with the target and the optimization level. The kit's results are meant to be reproducible (the same at any thread
  count, with or without OpenMP), and the tests compare many of them bit for bit. Where the kit wants an FMA it calls
  `fma()`. Give the flag explicitly: compilers differ in their default. Results can still differ between C libraries
  where they go through libm functions such as `sin`, `cos` and `exp`.
- Don't build with `-ffast-math`: it lets the compiler drop the NaN and infinity checks and reorder the compensated
  sums.
- LP64 only (64-bit `long` and pointers). The size checks rely on it: mk_linalg.h, mk_sparse.h, mk_geom.h, mk_sdf.h
  and mk_hhd.h, and the headers that include them, stop the build elsewhere with a `static_assert`.
- Each header compiles alone without warnings at `-Wall -Wextra -Wpedantic -Werror`, with and without `-fopenmp`.

## Testing

```sh
make check          # one row per check, then "wrong: N, skipped: M"; exits with N (at most 100)
make check OMP=1    # the same, with every test and planted fault also built with -fopenmp
make bench          # CG on a 64^3 Laplacian (OMP=1: threaded); BVH query rates against brute force
```

`make check` runs `test/check.sh`. It needs a POSIX shell, `timeout`, and a compiler that takes GCC's flags, including
`-fopenmp`: the thread tests are always built with it. Set `CC` to choose the compiler. Some tests use POSIX calls
(`fork`, `setitimer`, `mmap`). It checks that:

- each header compiles alone, with and without `-fopenmp`;
- each test `test/t_*.c` passes on the real headers: it exits 0 with at least one PASS row and no FAIL row, within
  `LIMIT` seconds (default 120). A SKIP row, one that cannot run on this machine, is printed and counted in
  `skipped`, never as a pass;
- every planted fault is caught. The headers hold 140 test-only faults, each under `MK_MUTATE == N`. check.sh builds
  the fault's test with `-DMK_MUTATE=N`, and the test must print a FAIL row; a timeout, a crash or a bare nonzero exit
  is not a catch (the rules are in `test/judge.sh`). Never define `MK_MUTATE` outside the tests;
- `t_par`, built with `-fopenmp`, gives the same bytes at 1, 2 and 12 threads, and its result hashes (BITS lines) equal
  those of a build without OpenMP;
- the checks discriminate: `test/negctl.sh` shows each thread row that no planted fault reaches failing on a known-bad
  input, and `test/drvctl.sh` shows each of judge.sh's verdicts on a known-bad program.

The mesh rows in check.sh use an external closed-mesh checker, `meshcheck`, which is not part of this repository.
check.sh runs `$MESHCHECK` when it is set (`make check MESHCHECK=path`), else `meshcheck` on `PATH`. With neither,
each of the 24 mesh rows prints SKIP, so an otherwise clean run ends with `wrong: 0, skipped: 24` (or more); a
`MESHCHECK` that is not an executable is WRONG. check.sh runs `$MESHCHECK -t FILE.obj` on the OBJ files that t_snet,
t_sepvol, t_dec and t_graft write, and reads the second line as tab-separated fields: the Euler characteristic in field
8, the genus in field 9, the signed volume in field 19 (positive for an outward mesh), and `clean` or `FAULT` in the
last field.

Two more rows skip when what they need is missing: t_rng's bit-for-bit comparison with an independent xoshiro256**
header (`make check RIBBON_RNG_H=path/ribbon_rng.h`, a header that defines `RbRng`, `rb_seed` and `rb_next`), and
t_interp's search of a 2^30-node table, whose 8 GiB `MAP_NORESERVE` mapping is refused under strict overcommit or an
address-space limit.

## Threads and determinism

With `-fopenmp` these run on threads: `mk_csr_matvec`, `mk_cg` and `mk_bicgstab` (all but IC(0)'s triangular solves),
`mk_sdf_mesh` and `mk_sdf_mesh_box`, the batch queries `mk_bvh_nearest_many`, `mk_bvh_inside_many` and
`mk_winding_number_many`, and `mk_hhd_split`. Without `-fopenmp` the pragmas compile out.

Results do not depend on the thread count. Work is cut into blocks whose length depends only on the item count (4096
items, or longer so there are at most 256 blocks; `mk_par.h`). Each block runs serially, and a sum adds the blocks'
partials in block order, never through OpenMP's `reduction`. So a sum over more than 4096 items is a sum of block
partials, not one running sum, at any thread count and without OpenMP. Called inside the caller's own parallel region,
the kit's loops run on the calling thread (OpenMP's default of one active level).

On a machine busy with other work, set `OMP_WAIT_POLICY=passive` and `OMP_NUM_THREADS` to the free cores. With
libgomp's default spin-waiting, threads waiting at a barrier compete with the other work: on a loaded x86-64 desktop,
CG on a 64^3 grid ran about 30 times slower on 12 threads than on one. Threads pay most in `mk_sdf_mesh` and the batch
queries; CG and BiCGSTAB are memory-bound and gain little.

## Stability

The experimental parts may change their results between releases: `mk_dec.h`, `mk_sepvol.h`, `mk_graft.h`,
`mk_ovl.h`, `mk_hhd.h`, and `mk_fftm` in `mk_fft.h`. Everything else is frozen: its calls and their results stay the
same. Names ending in an underscore (`mk_dop_`, `MK_PAR_FOR_`, `mk_bvh_count_wind_`) are internal, even where a test
calls them.

## Headers

| header | what | main calls |
|---|---|---|
| mk_version.h | the version | MK_VERSION, MK_VERSION_MAJOR, MK_VERSION_MINOR, MK_VERSION_PATCH |
| mk_par.h | the block scheme behind the threads | mk_par_blk, mk_par_nblk, mk_par_range |
| mk_vec.h | 3-vectors, unit quaternions (w, x, y, z), rigid transforms | `mk_v3_*`, `mk_q_*`, `mk_xf_*`, mk_lerp |
| mk_rng.h | `xoshiro256**` streams, splitmix64, Neumaier sum | mk_rng_seed, mk_rng_u01, mk_rng_gauss, mk_sum_add |
| mk_stats.h | selection, median, type-7 quantiles, Welford | mk_select, mk_median, mk_quantile, `mk_welford_*` |
| mk_interp.h | 1-D PCHIP, natural cubic spline, linear | mk_pchip_eval, mk_spline_slopes, mk_hermite_eval |
| mk_linalg.h | small dense matrices | mk_lu, mk_chol, mk_lstsq, mk_sym_jacobi, mk_tridiag, mk_solve3 |
| mk_fft.h | radix-2 FFT in 1-D to 3-D, real input; any length | mk_fft_plan, mk_fft_c, mk_fft3_c, mk_rfft, mk_fftm_c |
| mk_spec.h | windows, STFT, Welch spectra, coherence | mk_spec_init, mk_stft_frame, mk_welch, mk_csd, mk_coherence |
| mk_quad.h | quadrature: 1-D, triangle, tetrahedron | mk_gauss_legendre, mk_quad_adapt, mk_tri_rule, mk_tet_rule |
| mk_ode.h | RK4, Dormand-Prince, Verlet, ode23s, BDF | mk_rk4_step, mk_dopri, mk_verlet_step, mk_ode23s, mk_bdf_step |
| mk_opt.h | Nelder-Mead, Levenberg-Marquardt, golden, L-BFGS, CMA-ES | mk_nm, mk_lm, mk_golden, mk_lbfgs, mk_cmaes |
| mk_sparse.h | CSR matrices, CG, BiCGSTAB | mk_csr_from_triplets, mk_csr_from_bsr3_lower, mk_cg, mk_bicgstab |
| mk_geom.h | triangle queries, BVH, winding number | mk_closest_pt_tri, mk_ray_tri, mk_bvh_build, `mk_bvh_*` |
| mk_sdf.h | triangle mesh to a signed distance grid | mk_sdf_mesh, mk_sdf_mesh_box, mk_sdf_sample, mk_sdf_free |
| mk_snet.h | surface nets on labels or a field | mk_snet_build, mk_snet_build_iso, mk_snet_faces, mk_snet_split |
| mk_vpflow.h | volume-preserving flow of a vector potential | mk_vp_ok, mk_vp_vel, mk_vp_advect, mk_vp_jacdet |
| mk_hhd.h | experimental: Helmholtz-Hodge split of 2-D motion | mk_hhd_init, mk_hhd_split, mk_hhd_readout |
| mk_ovl.h | experimental: shared volume of two closed meshes | mk_overlap_volume, mk_ovl_bvh |
| mk_sepvol.h | experimental: interpenetrating closed meshes moved apart | mk_sepvol |
| mk_dec.h | experimental: decimation without folds or crossings | mk_decimate, mk_dec_defaults |
| mk_graft.h | experimental: marching tetrahedra; a mesh region replaced | mk_mtet, mk_graft_patch, mk_graft_stitch |

Each header's top comment describes its method, layouts and limits; each call's comment gives its arguments and
return values.

## Conventions

- Names: functions and types `mk_*`, macros and constants `MK_*`.
- Doubles throughout, except `mk_sdf`'s stored distances and the scalar fields `mk_snet_build_iso` and `mk_mtet` read,
  which are `float`.
- Flat, row-major arrays with 0-based indices: a dense matrix is `a[i * n + j]`; a mesh is `V[3 v + k]` (coordinates)
  and `T[3 t + j]` (vertex indices); a grid is indexed `(k n1 + j) n0 + i`, x fastest.
- Sizes are `int` in most calls; `long` or `size_t` where arrays can be larger (statistics arrays, triplet counts,
  point batches).
- Return codes: 0 for success and -1 for a bad argument, no memory or a non-finite result. `mk_cg` and `mk_bicgstab`
  return 1 when they stop short of the tolerance; `mk_graft_stitch` returns 1 or 2, with its result, when that result
  is not a clean closed mesh. Some calls use -2 for a second failure: no memory in `mk_ode23s` and `mk_bdf_step`, an
  unusable held mesh in `mk_sepvol`. Calls that return a value say how they fail: `mk_lstsq` and `mk_lm` return -1,
  `mk_nm` NaN; `mk_lbfgs` and `mk_cmaes` set `status` to -1 (and return NaN when they cannot start).
- Memory: a call that builds an object (`mk_fft_plan`, `mk_spec_init`, `mk_csr_from_triplets`, `mk_bvh_build`,
  `mk_sdf_mesh`, `mk_snet_build`, `mk_hhd_init`, `mk_graft_init`, `mk_mtet`) allocates it with `malloc`. Free it with
  the matching `*_free` (`mk_gmesh_free` for a mesh), which is also safe after a failed build. `mk_snet_faces` returns
  a `malloc`'d array for `free`. Calls that need scratch space either take a work array from the caller
  (`MK_RK4_WORK(n)`, `MK_DOPRI_WORK(n)`, `mk_fftm_c`'s `work`) or allocate and free it within the call.
- Borrowed pointers: a tree (`mk_bvh`), a graft (`mk_graft`) and a potential (`mk_vp`) keep pointers to the caller's
  arrays, which must outlive them.
- Threads: no call has global or static mutable state. An object a call takes as `const` (a plan, tree, matrix or grid)
  is only read, so threads may share it. Each thread needs its own `mk_rng`, `mk_sum`, `mk_welford` and work arrays.
- Non-finite values: `mk_sdf_mesh`, `mk_overlap_volume`, `mk_sepvol`, `mk_decimate` and `mk_graft_init` refuse
  non-finite coordinates (`mk_bvh_build` does not check them). The solvers return -1 for a NaN or vanishing pivot or a
  non-finite answer rather than return it. The statistics calls need NaN-free input and give NaN for an empty array.
  `mk_sum` gives the IEEE result once the sum is not finite.
- Units are the caller's. Angles are in radians, except `mk_dec_opts`' `minangle` and `fold`, in degrees.
- FFTs: forward is e^{-i}, inverse e^{+i}, and neither is normalized.
- Mesh queries break ties by the lowest triangle index, so their results do not depend on the tree's shape.

## Known limits

- No sparse direct solver, multigrid or ILU. IC(0)'s triangular solves are serial.
- `mk_bdf_step` is one step without error control: the caller chooses h.
- `mk_winding_number` costs O(triangles) per point, and so does each node of an odd row in `mk_sdf_mesh`.
- `mk_gauss_legendre` is O(n^2): n 65536 takes about 23 s.
- `mk_dopri` and `mk_ode23s` stop with -1 when the step falls below 1e-14 of the span, so a fast start on a long span
  needs a call of its own. `hmax` is met to within half an ulp of t.
- `mk_lu` refuses a pivot under 1e-13 of the matrix's largest entry, which turns away `mk_ode23s` steps past about 1e10
  on Robertson's problem. `mk_solve3` refuses only an exact zero pivot: a matrix singular to rounding gives 0 and a
  huge x.
- `mk_cg` and `mk_bicgstab` return 1 when a restart does not lower the true residual, so a tol at the attainable
  accuracy can end with 1. On -1, `relres` is the iteration's last.
- An infinite value makes `mk_welford`'s mean NaN from then on.
- `mk_tri_rule` past degree 5 and `mk_tet_rule` past degree 3 return 0 points.
- In mk_sdf.h a box one node thick gets no sweeps (its far nodes keep an upper bound), and above about 2^39 h a
  face-diagonal ray can be counted twice.
- `mk_sepvol` cannot clear parts that cross with no vertex inside each other; it keeps them whole and says so
  (`MK_SEPVOL_STILL_CROSSES`, `MK_SEPVOL_CRUSHED`). It does not undo a part that crosses itself.
- `mk_overlap_volume` can lose crossings on meshes spanning about 1e308.
- `mk_decimate` takes time quadratic in a vertex's valence. Its marks are 32-bit stamps: past about 9e7 collapses (an
  input near 1e8 vertices) a stale mark can read as current.
- `mk_graft_stitch` does not look for a vertex pinched between two fans; an ear can leave a mirrored pair; a rim vertex
  welded to one of S's vertices can move past `snap` by up to 1 % of an S rim edge.

## What the tests check

- t_vec: known turns, composition, matrix round trips, inverse transforms, slerp, constant-rate integration, `mk_lerp`.
- t_fft: a cosine's two bins, Parseval, round trips, a direct DFT in 1-D, 2-D and 3-D, split equal to complex, real
  input; `mk_fftm` against a direct DFT at 21 lengths up to 2160 (radices 4, 2, 3, 5, 7 and the prime 97).
- t_linalg: solve residuals, determinants, singular and indefinite refusals, least squares with dependent columns
  (two stored cases in test/data), a 6 x 6 symmetric matrix with known values, 3x3 paths, the tridiagonal solve's
  backward error on random systems, refusals.
- t_rng: xoshiro256** against its authors' reference output, uniform ranges, Gaussian moments, chi-squared of
  `mk_rng_below`, the sum 1 + 1e100 + 1 - 1e100 = 2; with `RIBBON_RNG_H`, seeded streams against an independent
  implementation bit for bit.
- t_stats: known medians, random arrays with ties against sorting, quantiles of 1..10, Welford against a direct sum,
  at an offset of 1e9, and merged.
- t_interp: PCHIP bit for bit against a reference implementation written out in the test (22000 queries),
  monotonicity, a natural spline on sin, straight lines, one-point tables, tiny secants, tables of 2^30 nodes.
- t_opt: Rosenbrock for every method, a coupled quadratic, an exponential fit, golden section, a rotated 10-D
  ellipsoid of condition 1e6, Rastrigin with restarts, and Brown-Dennis and Penalty I at their published minima.
- t_ode: RK4's order 4, Dormand-Prince against its tolerance and backward, Verlet's order 2 and bounded energy,
  Robertson's kinetics against Hairer and Wanner's reference values, Van der Pol at mu 1000 against its asymptotic
  period, BDF orders 1 and 2, an infinite first step, refusals.
- t_sparse: a Laplacian against a known solution, CG within its sqrt(kappa) iteration bound, IC(0) and block-Jacobi
  against a dense LU, a singular periodic Laplacian, BiCGSTAB on upwind convection-diffusion, b scaled by 2^600 giving x
  scaled bit for bit, extreme scales, size refusals.
- t_quad: Gauss-Legendre exact to degree 2n - 1 and not past it (n 1..20), Kronrod 15 exact to degree 22 and its error
  estimate, endpoint singularities, an oscillatory integral, triangle and tetrahedron rules exact to their degree.
- t_spec: I0 against tabulated values, Parseval for each window, a sine's amplitude, white-noise level 2 s^2 / fs,
  coherence and transfer function of a known noisy linear system, the Blackman-Harris -92 dB sidelobe, STFT padding.
- t_geom: closest point against an independent method (200000 cases), ray-triangle against Möller-Trumbore, rays
  aimed at every edge and vertex of a closed sphere (none slips through), tree queries against brute force, inside
  tests against analytic solids, and edge cases (NaN and 1e308 vertices, deep trees, ties, slivers).
- t_sdf: a sphere and a torus against their analytic distance, the sign at every node of nested shells, a random closed
  mesh against brute force, an open mesh, overlapping spheres with and without `MK_SDF_FILL`, refusals, sampling.
- t_snet: a ball raw and relaxed (closed, Euler 2, each vertex in its dual cell), a ball from its distance field, two
  labels sharing their interface, a torus, two voxels meeting at an edge, refusals; writes OBJ files for meshcheck.
- t_vpflow: div u, a ball's volume carried through the flow, the map's Jacobian, a rigid rotation against its exact
  answer (order 4), a periodic axis, refusals.
- t_hhd: a pinhole camera's flow (Longuet-Higgins and Prazdny 1980) over a wall, a tilted plane and a room against the
  formulas for focus, expansion and roll; the parts' orthogonality, divergence and curl; known potentials split back.
- t_sepvol: balls and cubes sunk into each other, a ball out of a held one (the held mesh unchanged byte for byte),
  crossed bars, nested shells; inside judged independently by each part's winding number; writes OBJ for meshcheck.
- t_dec: a sphere and a torus within the distance bound (by brute force both ways), topology kept, close walls without
  crossings, no new folds, held vertices, regions, an open rim, refusals; writes OBJ for meshcheck.
- t_ovl: two spheres' lens against its analytic volume on each axis, a box pair on the lattice (exact), disjoint and
  nested pairs, an open mesh's odd lines, refusals.
- t_graft: marching tetrahedra on a sphere and a torus (closed, Euler, volume), then a dent carved into a sphere and
  stitched back closed, against the analytic volume, with the untouched triangles kept; writes OBJ for meshcheck.
- t_regr: exact inputs that once gave wrong results, hangs or overflows in mk_ovl, mk_dec, mk_graft, mk_geom and mk_sdf;
  calls that could hang run in a child process under a CPU limit.
- t_par: with `-fopenmp`, every threaded call at 1, 2 and 12 threads, batches against single calls, and BITS hashes.

No test reaches Nelder-Mead's shrink step on its own, the NaN guard in the error norms of `mk_ode23s` and `mk_bdf_step`
with a finite value after the NaN, or `mk_decimate`'s mark wrap.

## References

- mk_vec.h: Shepperd 1978, "Quaternion from rotation matrix".
- mk_rng.h: Blackman and Vigna 2021, "Scrambled linear pseudorandom number generators" (xoshiro256**). Steele, Lea and
  Flood 2014, "Fast splittable pseudorandom number generators" (splitmix64). Box and Muller 1958, "A note on the
  generation of random normal deviates". Neumaier 1974, "Rundungsfehleranalyse einiger Verfahren zur Summation endlicher
  Summen".
- mk_stats.h: Hoare 1961, "Algorithm 65: Find". Hyndman and Fan 1996, "Sample quantiles in statistical packages".
  Welford 1962, "Note on a method for calculating corrected sums of squares and products". Chan, Golub and LeVeque 1983,
  "Algorithms for computing the sample variance: analysis and recommendations" (shifted data, pairwise merge).
- mk_interp.h: Fritsch and Butland 1984, "A method for constructing local monotone piecewise cubic interpolants".
- mk_fft.h: Cooley and Tukey 1965, "An algorithm for the machine calculation of complex Fourier series". `mk_fftm` uses
  Stockham's self-sorting form.
- mk_spec.h: Welch 1967, "The use of fast Fourier transform for the estimation of power spectra: a method based on time
  averaging over short, modified periodograms". Harris 1978, "On the use of windows for harmonic analysis with the
  discrete Fourier transform" (Blackman-Harris).
- mk_quad.h: Piessens, de Doncker-Kapenga, Überhuber and Kahaner 1983, "QUADPACK: a subroutine package for automatic
  integration" (the Gauss-Kronrod 7-15 constants). Radon 1948, "Zur mechanischen Kubatur" (7-point triangle rule).
  Keast 1986, "Moderate-degree tetrahedral quadrature formulas".
- mk_ode.h: Dormand and Prince 1980, "A family of embedded Runge-Kutta formulae". Shampine and Reichelt 1997, "The
  MATLAB ODE suite" (ode23s). Hairer and Wanner 1996, "Solving Ordinary Differential Equations II" (the Robertson
  reference in the tests).
- mk_opt.h: Nelder and Mead 1965, "A simplex method for function minimization". Lagarias, Reeds, Wright and Wright 1998,
  "Convergence properties of the Nelder-Mead simplex method in low dimensions". Levenberg 1944, "A method for the
  solution of certain non-linear problems in least squares". Marquardt 1963, "An algorithm for least-squares estimation
  of nonlinear parameters". Moré 1978, "The Levenberg-Marquardt algorithm: implementation and theory" (the scaling).
  Nocedal 1980, "Updating quasi-Newton matrices with limited storage". Nocedal and Wright 2006, "Numerical
  Optimization", Algorithms 3.5 and 3.6 (the strong-Wolfe line search). Hansen 2016, "The CMA evolution strategy: a
  tutorial". Moré, Garbow and Hillstrom 1981, "Testing unconstrained optimization software" (test problems).
- mk_sparse.h: Hestenes and Stiefel 1952, "Methods of conjugate gradients for solving linear systems". Meijerink and
  van der Vorst 1977, "An iterative solution method for linear systems of which the coefficient matrix is a symmetric
  M-matrix" (IC(0)). Manteuffel 1980, "An incomplete factorization technique for positive definite linear systems" (the
  diagonal shift). van der Vorst 1992, "Bi-CGSTAB: a fast and smoothly converging variant of Bi-CG for the solution of
  nonsymmetric linear systems".
- mk_geom.h: Ericson 2005, "Real-Time Collision Detection", section 5.1.5 (closest point on a triangle). Woop, Benthin
  and Wald 2013, "Watertight ray/triangle intersection". Jeannerod, Louvet and Muller 2013, "Further analysis of Kahan's
  algorithm for the accurate computation of 2 x 2 determinants" (a*b - c*d). Möller and Trumbore 1997, "Fast, minimum
  storage ray-triangle intersection" (the strict edge test). Wald 2007, "On fast construction of SAH-based bounding
  volume hierarchies" (binned SAH). Van Oosterom and Strackee 1983, "The solid angle of a plane triangle" (winding
  number).
- mk_sdf.h: Batty 2007, SDFGen (exact distances in a band, then sweeps; sign by ray parity). The winding number is
  mk_geom.h's.
- mk_snet.h: Gibson 1998, "Constrained elastic surface nets: generating smooth surfaces from binary segmented data".
  Taubin 1995, "A signal processing approach to fair surface design".
- mk_hhd.h: Bhatia, Norgard, Pascucci and Bremer 2013, "The Helmholtz-Hodge decomposition: a survey".
  Longuet-Higgins and Prazdny 1980, "The interpretation of a moving retinal image" (the tests' camera flow).
- mk_dec.h: Garland and Heckbert 1997, "Surface simplification using quadric error metrics". Dey, Edelsbrunner, Guha
  and Nekhayev 1999, "Topology preserving edge contraction" (the link condition).
- mk_graft.h: Freudenthal 1942, "Simplizialzerlegungen von beschränkter Flachheit" (six tetrahedra a cube).

## Changes

1.6.4:
- The 64-bit `long` checks are spelled `_Static_assert`, so the headers compile as C11 again.

1.6.3:
- `mk_ode23s` chooses its own first step when `*h` is infinite, as `mk_dopri` does.
- The headers whose sizes and counts need a 64-bit `long` (mk_linalg, mk_sparse, mk_geom, mk_sdf, mk_hhd) refuse other
  targets at compile time.
- `make check` runs without meshcheck or an RNG header: those rows print SKIP, and the run ends `wrong: N, skipped: M`.
- t_ode reports the absolute error of y' = -2ty.
