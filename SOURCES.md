# Where each piece came from (merge, don't pick; 2026-11-13)

The inventory these copies were found in is gauge/notes/math-tools-inventory.md. Every copy was read; none was edited.

## mk_vec.h
- body-claude src/mathx.h: the v3/quat value types and the inline style; the matrix-free rotation; axis-angle and the
  rotation-vector form; rigid transforms (apply, compose, inverse).
- metal chassis.c: the same rotation formula (it confirmed the 18-multiply form); Shepperd's matrix-to-quaternion with the
  largest-diagonal branch (quat_of_matrix); unrotate as rotate by the conjugate.
- Added: to_mat, slerp, integrate (exact for constant omega, through the exponential map). mathx.h's small-angle
  rotation vector returned an unnormalized quaternion; it is normalized here.
- Not taken: mathx.h's distal_dir (anatomy, not math), clampd and smoothstep (not vector math).

## mk_fft.h
- helium-0k src/fft.{c,h}: the plan (bit-reversal table plus twiddles), the complex interface, the power-of-2 refusal.
- metal wind.c: split re/im storage with cos and sin tables (no complex type needed), forward e^-i, inverse e^+i,
  unnormalized, and the 2-D transform with columns gathered eight at a time for the cache (fft2).
- One plan serves both forms; they share the same arithmetic, so their results are bit-identical (tested).
- Added: a 3-D complex transform along each axis, and a real-input transform through one half-length complex FFT.
- Not merged: the ~14 other hand-written radix-2 copies (mcphys sound and sound2, depthref, mhphys vision.c, body-claude
  src) are the same algorithm with no feature these lack. FFTW (eye) and cuFFT (cudakit, helium lq) stay where they are.

## mk_linalg.h
- mcphys-kc linalg.c: LU with partial pivoting and its relative singular test that also refuses a NaN pivot; Cholesky;
  Householder least squares with column scaling and dropped rank-deficient columns; cyclic Jacobi with the stable tangent.
- metal chassis.c jacobi3: the scale-relative stopping rule (confirmed the linalg form); ascending order of values with
  their vectors (from chassis.c and stepcad cloudfit.c eig3).
- mcphys2 suspension.c solve3: 3x3 with partial pivoting and a finite-result check. It replaced neuro/light/matfit.c's
  solve3 (Cramer's rule: no pivoting, no singular check) and was kept over it.
- stepcad cloudfit.c eig3: the rotation by atan2 was not taken (the tangent form needs no trig and is as accurate).
- Added: LU determinant; mk_sym3 returns a right-handed frame; Jacobi reports its sweeps or -1 if it did not converge.
- metal world.c (ground column) and roadline.c (road heat), 1.0.5: the Thomas elimination both write the same way, f =
  a[i] / b[i - 1], b[i] -= f c[i - 1], d[i] -= f d[i - 1], then back substitution; mk_tridiag keeps it so their answers
  stay bit-identical, and adds the refusals (NaN, zero or rounding-small pivot, non-finite answer). mk_interp's spline
  solve (from anatomy.c, the c', d' form) stays inside mk_spline_slopes: the verification work proved it as written.

## mk_rng.h
- dna ribbon/include/ribbon_rng.h: xoshiro256** and its splitmix64 seeding from (seed, stream). Same seeds give the same
  stream as ribbon's, bit for bit (tested over 100000 draws); its (0,1] uniform.
- mcphys sound2/dsp.h: the guard against an all-zero state; the idea of one named stream per use.
- Changed: integers below n by rejection, without the bias of ribbon's multiply; Neumaier's sum in place of Kahan's (right
  when an addend exceeds the total: 1 + 1e100 + 1 - 1e100 gives 2); Box-Muller keeps its second value (metal wind.c's
  gauss drew two uniforms per value).
- Not taken: dsp.h's sum-of-four-uniforms "Gaussian", bounded at 3.46 on purpose for the engine sound. It is not a
  Gaussian, so it stays in dsp.h. The other splitmix64 copies (eye, fly, neuro, mhphys, hair) are the same function.

## mk_opt.h
- mcphys sound2/nm.h (itself after gauge wifi3d.c): the interface (function, context, start, step, iterations, tolerance)
  and the reflect/expand/contract/shrink skeleton.
- Changed to the standard rules (Lagarias et al. 1998): a reflection is kept only if it beats the second worst (nm.h
  kept it whenever it beat the worst), with both outside and inside contraction. It also counts function calls.
- Levenberg-Marquardt, the house's five copies (stepcad cloudfit.c, body-claude fit.c and mhrbody.c, body-motion-fork
  fit.c, engine props/gen.c) agree on the method. Taken:
  - from cloudfit.c, Marquardt's diagonal scaling, the 10x jump when the damped system fails, and the 1e-12 relative stop;
  - from body-claude fit.c, a small absolute term in the damping, so a zero column cannot stall it, and Cholesky for
    the solve;
  - from gen.c, forward-difference slopes when no Jacobian is given.
  On success the damping is divided by 3, on failure multiplied by 4 (cloudfit).
- golden section: the standard form; the copies (helium lq, metal moon_check, engine cam.c, chem trans.c) are all it.

## mk_ode.h
- body-claude src/gyncycle.c: the Dormand-Prince 5(4) tableau and its error coefficients, mixed tolerance
  atol + rtol max(|y|, |y_new|), the step factor 0.9 err^-1/5 clipped to [0.2, 5], and non-finite errors rejected.
- Changed: first same as last (the last stage of an accepted step is the next step's first: 6 calls a step, not 7);
  runs backward in time as well; hmax; counts steps, rejections and calls.
- RK4: the five copies (metal indoor_check and chassis_test, motor.c, mhphys calcium.c, neuro ns.c) are the classic
  form; here the caller owns the work memory, so a step never allocates.
- Velocity Verlet: added (no house copy found), symplectic, the acceleration carried between steps.

## mk_stats.h
- jobwatch.c and depthref ctcpair.c: the median as the mean of the two middle values for an even count, and NaN for an
  empty array (ctcpair). metal sense_check.c took the upper middle value; that is not the median, and was not kept.
- Changed: selection (Hoare) in place of a full sort, O(n) on average.
- Added: type-7 quantiles (R and numpy default) and Welford's running mean and variance; the house copies had none
  in shared form.

## mk_interp.h
- aeromap.c (proved by the verification work): pchip's interior slope (Fritsch-Butland), end slope (three-point shape-preserving), the
  two-node straight line, the binary-search locate and the clamping. Each expression is written as aeromap.c writes it,
  so results are bit-identical: 0 of 22000 test queries differ.
- body-claude src/anatomy.c: the natural cubic spline slopes (Thomas algorithm), and the split into a slope table plus
  a Hermite evaluation. anatomy.c's pchip ends (the end secant) are not kept: aeromap's are shape-preserving.
- Added: a linear interpolation with the same locate.

# Stage 2 (2026-11-13, the owner: "build all this: a CPU CG solver, quadrature, BFGS or CMA-ES, stiff time steppers, spectral analysis, ray-tracing trees and point-to-triangle distance")

## mk_sparse.h
- cudakit spd3.cu: the solve interface (x the first guess and the answer, tol on |r| / |b|, iterations and relative
  residual out, 0 / 1 / -1), the stop when p.Ap is not positive (lost definiteness or NaN), the 3x3 block-Jacobi
  preconditioner with its explicit 3x3 inverse (k_dinv), and its input layout (lower-triangle block CSR), which
  mk_csr_from_bsr3_lower takes unchanged so a caller can move between the GPU and the CPU.
- bioheat.c: Jacobi PCG, warm-started, converging on the residual relative to the start's.
- gauge lattice.c: the singular periodic Poisson case. The residual is kept mean-free (MK_PC_MEANFREE), because rounding
  otherwise feeds the constant mode back and CG stalls (288 x 160 ran to the cap there).
- metal wick.c (CG on its pressure nodes, its diagonal kept for a preconditioner) and fable solve.cu (block-Jacobi CG)
  are the same method; nothing further was taken from them.
- Added (no house copy): CSR from triplets with duplicates summed; IC(0) with Manteuffel's diagonal shift on breakdown;
  BiCGSTAB, right-preconditioned.

## mk_quad.h
- body-claude bfe.c and fable elem.c: the 4-point tetrahedron rule (a 0.5854101966249685, b 0.1381966011250105, w 1/24)
  and bfe.c's 3-point triangle rule ({1/6, 1/6}, {2/3, 1/6}, {1/6, 2/3}); fable's 2x2x2 hex rule is Gauss-Legendre n 2
  in each direction.
- engine chem trans.c: a Gauss-Kronrod 7-15 for two components. Its constants were checked against QUADPACK's; the kit
  carries them to full precision (trans.c's have 15 digits).
- Added: Gauss-Legendre for any n, adaptive bisection on the largest error, Radon's 7-point triangle (degree 5), Keast's
  5-point tetrahedron (degree 3).

## mk_opt.h, stage 2
- body-motion-support skin-buckling/minimize.c (L-BFGS, 48 pairs, written inline): its curvature guard,
  s.y > 1e-12 |s||y|, before a pair is kept.
- glance-core walk/gaittrain.c (CMA-ES, Hansen's, full covariance): the same default constants (they agree with the
  kit's to the formula), and per-coordinate scaling, taken as mk_cmaes_opts.scale (the first covariance diag(scale^2)).
- Added: strong-Wolfe line search with safeguarded cubic zoom; IPOP restarts are left to the caller (the test shows how).

## mk_ode.h, stage 2
- mcphys damper.c (backward Euler, Newton on both unknowns, a nested fallback when it does not settle) and mcphys-kc
  steer.c (a locally implicit step): implicit steps solved by Newton, and a step left unchanged when Newton fails, for
  the caller to retry smaller. neuro ns.c's implicit fallback is the same idea for its long steps.
- Added: ode23s (Shampine and Reichelt 1997). Chosen over ROS3P for L-stability: ROS3P's stability function tends to
  about 0.73 at infinity, so very stiff transients ring. BDF2 with a variable step ratio; numerical Jacobians; a step
  cap per call (maxsteps, also given to Dormand-Prince) so a bad problem stops instead of running on.

## mk_spec.h
- mcphys sound2 trackstat.c: Welch-style averaging (Hann, hop N/2, |X|^2 summed per block). Its unnormalized level is
  replaced by the one-sided density and power scalings.
- orders.c and specimg.c (an MDCT and a Gaussian STFT with PGHI phase) are specialized transforms and stay where they are.
  hrv_stress.c has no PSD (time-domain statistics only).
- Added: Hamming, 4-term Blackman-Harris, Kaiser (I0 by its series); cross-spectra, transfer function, coherence; STFT
  frames that zero-pad past either end.

## mk_geom.h
- body-claude bfe.c closest_tri_w (Ericson 5.1.5): taken as written, with the closest point returned as well. Its tree
  (tv_*) answers ties with the lowest triangle index; the kit's queries do the same, so results do not depend on the tree.
- skinlight.c's bvh.inc (two copies: bodyfit, body-male fit) and hair hairpt.c: a median split on the longest centroid
  axis, leaves of 4, float boxes. Replaced by a binned SAH build (16 bins) in double, leaves of up to 4, and a leaf kept
  when a split costs more.
- stepcad cfdsurf.c (ray_hit with a skip triangle, pt_tri_dist, winding) and body-claude sdfgrid.c (ray parity "nudged
  off vertices and edges"): sdfgrid.c's problem, a parity ray through an edge or vertex counted twice, is solved here
  by a watertight ray-triangle test (Woop, Benthin and Wald 2013) and a majority of three fixed directions.
- gauge meshtop: the generalized winding number (Van Oosterom and Strackee), as mk_winding_number.
- experiment/cfd/road/compare/offbody_x.c: the same nearest-point use, by brute force; it has no tree to take from.

# Stage 3 (2026-11-13)
- mk_sum: the non-finite case, from the maintainer's cudakit cross-check (an infinite addend gave NaN).
- mk_welford: cudakit reduce.cu does Welford per thread and Chan's merge between threads (wf_merge), and is accurate
  because its merges are pairwise. Its merge is taken as mk_welford_merge, with the same formula. The sequential update
  gets the shifted-data form instead (Chan, Golub and LeVeque 1983): a running mean of x - x0 rounds at the spread's
  scale, not the mean's, which a sequential stream cannot get by merging.
- mk_rfft: its second twiddle is now -conj of the first (e^(-2 pi i (m-k)/n) = -conj(e^(-2 pi i k/n))), one sin/cos a
  bin, not two: 4.2 to 2.7 us at n 512 (MEASURED in the ggml evaluation, scratch/ggml-eval/src/fft.c). Its existing
  real-against-complex test row covers it.

# 1.0.0 (2026-11-14)
- the verification work's combined patch, mathkit-qc-2026-11-14.patch (scratch/the verification work/mathkit/findings), applied as written:
  mk_lstsq (column pivoting, column scaling by the largest entry), mk_sym_jacobi (scaled stop, non-finite refusal),
  mk_lm (no stop while lambda > 1e-9; MINPACK's running-maximum scaling), mk_lu_det (frexp/ldexp), mk_chol (Inf
  diagonal), mk_fft3_c (refuse first). Test inputs: test/data/lstsq_keeps_dependent.txt and
  lstsq_drops_independent.txt, copied from the verification work's findings (sha256 6d89b38a8f19..., e358e38b9a00...); the MGH
  problems Brown-Dennis and Penalty I from the verification work's mgh.h definitions.
- The contract changes (lstsq non-finite refusal, solves reporting non-finite answers, mk_nm's *evals) are the maintainer's.

## mk_vpflow.h (1.1.0)
- b-helper's mk_vpflow.h (s2fix hand-off 7, 2026-11-15, sha256 866ff612...): the whole method, the potential's cubic
  B-spline evaluation and its curl, the periodic axes, the fading to 0 past the grid. Its first user is fmap-v22.c in the
  same folder, whose calls (the mk_vp initializer in field order, mk_vp_advect) are unchanged.
- Kept to the kit: its RK4 step is mk_ode.h's mk_rk4_step (the same expressions, so the moved points are bit-identical:
  200000 points of the test field, 40 steps, the same hash at -O2 and -O3 -march=native). mk_vp_ok refuses a grid; eval
  gives NaN for a NaN or infinite coordinate, a zero h or an empty axis, and takes whole periods out of a periodic
  coordinate past 2^30 cells before it becomes an int (it was undefined there). The planted faults moved to 131-132.
- Speed (MEASURED on the development machine, one pinned thread, 200000 points x 40 steps): -O3 -march=native 351 ns a point-step against
  the submission's 343, -O2 864 against 826. mk_rk4_step costs nothing (the same 686-695 ns with the step written out
  in an intermediate build); the 2-5 % is the input checks. One trap: a second test on the cold path pushed mk_vp_eval
  past gcc's inline limit at -O3 and doubled the time (717 ns), so the off-grid path tests x with one sum.
- No house copy: REGISTRY's methods table had none; mesher's thighpack is an area-preserving radial map per ray, a
  different method.

## mk_snet.h (1.1.0)
- s-forms' mk_snet.h (seed, through the maintainer #35): labelled grids, one vertex per dual cell shared by every label there,
  the Taubin pair held inside each dual cell, the faces wound outward on the shorter diagonal. Results are bit-identical
  to it (vertices and triangles, labels and three-label grids, raw and relaxed; MEASURED, six cases).
- experiment/cfd wrap.c: memory with the surface. wrap.c keeps two layers of vertex ids; the kit keeps each vertex's
  dual cell (ascending) and rebuilds two or three layers when it needs them, for the relaxation's neighbours and the
  faces, in place of a 4-byte table over every dual cell. On a 400^3 grid holding one ball: peak 31 MB against 275 MB,
  build 0.72 s against 1.80 s, faces 1.08 s against 1.21 s (MEASURED, the development machine, -O2). Also wrap.c's test for a cell whose
  inside (or outside) corners fall in two groups, as mk_snet_split: wrap.c fills such cells; the kit reports them and
  leaves the labels to the caller.
- body-claude src/mesh.c and wrap.c: a scalar field's vertex at the mean of the linear iso crossings (mk_snet_build_iso),
  held to [0.05, 0.95] of each edge as mesh.c does, so a sample exactly at iso does not collapse edges; on a ball's
  distance field the raw vertices are 0.045 cells from the sphere against 0.47 for labels.
- Not taken: mesh.c's sparse 8^3 blocks over a field function, its tie flips for ambiguous faces and its projection onto
  the field's zero, and skin.c / fvskin.c's remeshing (split, collapse, flip, tangential smoothing): they belong to a field
  sampled on demand or to remeshing, not to a stored grid. seed's grown core/mesh.c (not in the inventory) is another
  field-based copy of the same kind.

## mk_sdf.h (1.2.0)
- s-bone's meshsdf.c (stage bone dev, 2026-11-15, the maintainer.md #37): the method whole, Batty's SDFGen. Exact distances on the
  nodes within two cells of each triangle, two rounds of eight sweeps passing on each node's nearest triangle, the sign by
  the parity of x-rays through the rows set off a fraction of h, the crossings bucketed by row and sorted. The kit keeps its
  grid placement (n = ceil((box + 2 margin) / h) + 1, node 0 at the box less the margin) and node order (i fastest), so its
  K_GRID users wrap the kit in a few lines. On the body mesh its sign agrees with the kit's but at 117 of 7.1 M nodes at
  5 mm (393 of 32 M at 3 mm): the mesh's 97-edge hole makes rows odd, and meshsdf's parity puts part of them on the wrong side.
- mk_geom: the band and the sweeps call mk_closest_pt_tri (Ericson) in place of meshsdf's and sdfgrid's own copies;
  MK_SDF_EXACT uses mk_bvh_nearest, bounded by the sweeps' value; the winding number is mk_winding_number.
- body-claude src/sdfgrid.c: its FILL flood (outside nodes the grid's faces cannot reach are inside) as MK_SDF_FILL. Its
  parity tested every triangle against every row, kept at most 4096 crossings a row without saying so, and set the rows off
  by a fixed 1.3e-7 m: the bucketed rows replace it. Its band-only distance (+-BAND beyond) and BLUR are not taken: the
  sweeps give a true distance continuation, and a blur belongs to the caller.
- experiment/cfd levelset.h: the flood fill is the same as sdfgrid's FILL. Its EDT (Felzenszwalb) far field and its closing
  ball (wrap_levelset) are not taken: the sweeps reach the far field at triangle accuracy, and the closing ball is a wrapper's
  job for soups with gaps, not a mesh's distance.
- metal/ref/hair hairgen.c (found by the signature, not in the inventory): its lattice placement, node 0 on the world
  lattice h Z, as MK_SDF_LATTICE (nodes computed from integer indices, so two grids share theirs bit for bit). Its sign
  from the normal at the closest point is not taken (it fails at sharp edges; the head is smooth).
- cudakit_sdf.cuh: the naming only. Its primitives are sd_sphere, sd_capsule and so on, negative inside; the kit's grid
  keeps that sign and the mk_sdf_ prefix.
- Added: odd rows cast again on two nearby rays and, if still odd, signed by mk_geom's mk_winding_number (|w| > 1/2); nodes
  nearer than the rays' offset signed by a ray of their own; the refusals.
- Speed on fullset-7.0.0's rest-uv.obj (887091 vertices, 1774083 faces, a 97-edge hole; MEASURED on the development machine at load 5-12, one
  thread unless said, -O2; scratch/mathkit-120-work/speed*.log):

  | at 5 mm (222 x 375 x 85 nodes) | time | inside |
  |---|---|---|
  | mk_sdf_mesh | 16.8-20.1 s (band 1.7, sweeps 10.0, the hole's odd row 8.5 by the winding number) | 59.495 l |
  | mk_sdf_mesh, MK_SDF_EXACT | 62.6-66.7 s | 59.495 l |
  | meshsdf | 11.1-11.9 s | 59.503 l (117 nodes on the wrong side) |
  | sdfgrid, 1 thread | 201 s (distance only within 8 mm; its OBJ read included) | 59.50 l |
  | sdfgrid, 12 threads | 46.5 s | 59.50 l |

  At 3 mm (32.4 M nodes): mk_sdf_mesh 56.8-73.3 s, EXACT 249 s, meshsdf 38.2-48.7 s. The band and sweeps run as fast as
  meshsdf's; the difference is the winding number on the hole's rows, at O(triangles) a node (Barill et al.'s fast winding
  number through the tree is in "Not yet"). sdfgrid's reader takes only "f a b c": on rest-uv.obj's "f a/b/c" faces it
  read garbage and reported 0.00 l inside, so it was timed on a plain copy of the mesh.
- The wrapper for seed core's K_GRID Grid (stage bone's call, unchanged: 1 done, 0 refused or out of memory; g->v is the
  kit's malloc'd array, freed by the caller as before). Node 0 is lo in double; g->lo holds it rounded to float.

    /* mesh_sdf on mathkit's mk_sdf_mesh, for seed core's K_GRID Grid (n, h, lo, in, v): the call and the return of stage bone's
       meshsdf.c (1 done, 0 out of memory or refused), g->v malloc'd for the caller to free as before. */
    #include "mathkit/include/mk_sdf.h"
    #include "meshsdf.h"
    int mesh_sdf(const float (*v)[3], int nv, const int (*t)[3], int nt, float h, float margin, Grid *g)
    {
    	double *V = malloc(sizeof *V * 3 * (size_t)(nv > 0 ? nv : 1));
    	if (!V) return 0;
    	for (int i = 0; i < nv; i++) for (int a = 0; a < 3; a++) V[3 * (size_t)i + a] = v[i][a];
    	mk_sdf s;
    	const int r = mk_sdf_mesh(V, nv, &t[0][0], nt, h, margin, 0, &s);
    	free(V);
    	if (r) return 0;
    	for (int a = 0; a < 3; a++) g->n[a] = s.n[a], g->lo[a] = (float)s.lo[a], g->in[a] = (float)s.in[a], g->in[3 + a] = (float)s.in[3 + a];
    	g->h = h, g->v = s.d; /* the grid's own array: g->v is freed by the caller, as meshsdf's */
    	return 1;
    }

## mk_sepvol.h (1.3.0)
- cad-helper's sepvol.c (airway/method/fix, 2026-11-15, the maintainer.md #41; for the maintainer-fable's C1 films, PROPOSAL.md item 3):
  the API's shape (one mesh's parts, or a moving part out of a held one that stays byte for byte), the inside test (a point
  just inside a vertex covered twice: even ray parity), half the depth each in single mode and the whole in held mode, plus
  GAP, at most STEP a round, the move spread over RINGS rings, rounds, the volume reported. Its push along the vertex's own
  inward normal is replaced: at a 90 degree edge that normal is the diagonal, and two cubes 1 mm into each other kept 182
  crossing pairs after ten runs (MEASURED with 1.2.2's mk_geom; 165 with 1.3.0's). The kit pushes toward the nearest point
  of the other part's surface (mk_bvh_nearest_lab), the vertex's own patch and 3 rings around it left out of the search.
- body-claude src/selfsep.c (the maintainer, 2026-10-19): its rule for crossing face pairs, taken whole: a corner behind the other
  face's plane and over it goes out along that face's normal by half the depth plus half the gap, the pushes averaged per
  vertex. The vertex test alone left the two cubes with 40 crossing pairs and no vertex inside (an edge through a face
  between two outside corners); the pair rule clears them. Added: in a crossing where no corner is over the other face (an
  X), every corner behind the other's plane moves; and faces facing each other nearer than gap/2 are pushed to gap, which
  keeps meshcheck's winding probe (1 % of an edge off each face) out of the other part (without it ballB's nearest face
  centre was 0.0046 mm from ballA, at gap 0.05 mm). Its crumpled-patch relaxation is not taken: a surface folded through
  itself has no other part to leave, and relaxing it is meshfix's job, not a volume separation.
- mesher's selfx (now inside meshcheck): its pair definition (two faces sharing no vertex, an edge of one through the
  other) is the one counted in rep.cross; on the results the kit's count and meshcheck's agree (40 and 40 on a stalled cube
  run). meshfix's keep mode (cad-helper's, the same folder): umbrella relaxation of crossing and zero-area faces' vertices,
  grown by RINGS; not taken, for the same reason as selfsep's relaxation. Its fan splitting changes the vertex list, which
  mk_sepvol must not.
- Added: no vertex moves outward along its own normal (that part of a move is removed), so a part shrinks back into
  itself and never grows into a third one; mk_bvh_box in mk_geom.h for the pair search; the refusals.
- Results (test/t_sepvol.c, gap 0.05 mm, step 0.25 mm, 3 rings; MEASURED): twoball 73 vertices inside and 88 crossing pairs
  to 0 and 0 in 3 rounds, volume -1.19 %; ballB out of ballA 38 and 88 to 0 in 6 rounds, -3.63 %, ballA unchanged; the two
  cubes 98 and 112 to 0 in 6 rounds, 0.0767 ml lost against the 0.064 ml overlap (1.20 times); every result clean in
  meshcheck.
- For cad-helper's command line, a caller to copy (compiled and run on twoball, ballB in ballA and cubes2.awk's cubes: each
  exit 0 and clean in meshcheck):

    /* sepvol on mathkit's mk_sepvol: read an OBJ (v and f lines; f a/b/c takes a), separate, write it back with every other line
       as it came. HELD NULL for one mesh's parts; a held OBJ's path for the held mode. 0 done, 1 not all separated, 2 an error. */
    #include "mathkit/include/mk_sepvol.h"
    #include <stdio.h>
    static int readobj(const char *path, double **V, int *nv, int **T, int *nt)
    {
        FILE *f = fopen(path, "r");
        if (!f) return -1;
        char l[1024];
        int cv = 4096, ct = 8192;
        *V = malloc(sizeof **V * 3 * cv), *T = malloc(sizeof **T * 3 * ct), *nv = *nt = 0;
        while (*V && *T && fgets(l, sizeof l, f)) {
            if (l[0] == 'v' && l[1] == ' ') {
                if (*nv == cv && !(*V = realloc(*V, sizeof **V * 3 * (cv *= 2)))) break;
                double *p = *V + 3 * (size_t)(*nv)++;
                if (sscanf(l + 2, "%lf %lf %lf", p, p + 1, p + 2) != 3) { fclose(f); return -1; }
            } else if (l[0] == 'f' && l[1] == ' ') {
                if (*nt == ct && !(*T = realloc(*T, sizeof **T * 3 * (ct *= 2)))) break;
                char *s = l + 2;
                for (int k = 0; k < 3; k++) {
                    (*T)[3 * (size_t)*nt + k] = (int)strtol(s, &s, 10) - 1;
                    while (*s && *s != ' ') s++;
                }
                (*nt)++;
            }
        }
        fclose(f);
        return *V && *T ? 0 : -1;
    }
    int sepvol_obj(const char *in, const char *out, const char *held, double gap, double step, int rings, int rounds)
    {
        double *V, *HV = NULL;
        int *T, *HT = NULL, nv, nt, hnv = 0, hnt = 0;
        if (readobj(in, &V, &nv, &T, &nt) || (held && readobj(held, &HV, &hnv, &HT, &hnt))) return 2;
        mk_sepvol_report r;
        if (mk_sepvol(V, nv, T, nt, HV, hnv, HT, hnt, gap, step, rings, rounds, &r)) return 2;
        printf("%s: %ld vertices inside another part -> %ld, %ld crossing face pairs -> %ld, %d rounds, %ld moved (worst %.4g), volume %.6g -> %.6g\n",
               out, r.inside0, r.inside, r.cross0, r.cross, r.rounds, r.moved, r.worst, r.vol0, r.vol1);
        FILE *f = fopen(in, "r"), *o = fopen(out, "w"); /* every line as it came, the vertices' positions moved */
        if (!f || !o) return 2;
        char l[1024];
        int i = 0;
        while (fgets(l, sizeof l, f))
            if (l[0] == 'v' && l[1] == ' ') fprintf(o, "v %.9f %.9f %.9f\n", V[3 * (size_t)i], V[3 * (size_t)i + 1], V[3 * (size_t)i + 2]), i++;
            else fputs(l, o);
        fclose(f);
        if (fclose(o)) return 2;
        free(V), free(T), free(HV), free(HT);
        return r.inside || r.cross ? 1 : 0;
    }
    #ifdef SEPVOL_EX_MAIN
    int main(int argc, char **argv) /* metres: gap 0.05 mm, step 0.25 mm, 3 rings, 30 rounds */
    {
        return argc < 3 ? 2 : sepvol_obj(argv[1], argv[2], argc > 3 ? argv[3] : NULL, 5e-5, 2.5e-4, 3, 30);
    }
    #endif

# 1.4.0 (2026-11-16): three mesh methods from mesher-fable (the maintainer.md #43; written on mk_geom.h, for mesher's C1 pass)

## mk_dec.h
- mesher-fable's folddec.c (2026-11-15), the method whole: objdec's Garland-Heckbert quadrics with half-edge collapses only (the
  kept vertices are the input's, at their input coordinates; the caller writes its own text back), the link condition, the turn,
  edge and angle refusals, and folddec's own: fold lengths that may not grow (per region through a caller's function: foldlen's
  regions about x = 0, z = CZ become a callback), the two-sided distance bound (owned vertices re-measured at every collapse, sample
  points against mk_geom's tree), the crossing refusal on an insert-only hash grid, held and seam vertices, open rims, the
  crossing and needle clean-up. Its globals became one context; its millimetre constants scale with the caller's unit: the grid
  cell is half maxedge, at most 4 mean edges (folddec's 2 mm at its 4 mm maxedge), the sample spacing dist / 0.3 (1 mm at 0.3 mm), the fold tolerance
  1e-9 mean edges, the needle edge 20 minalt, the degenerate and parallel tests relative to the triangle's own size.
  mk_dec_defaults(o, mm) gives mesher's values for a unit. MEASURED: on an icosphere (10242 vertices) and a torus (10240) the
  kit's output is byte for byte folddec's; on the latissimus dorsi r (b-helper's S2 A hand-off 8, 1812456 vertices) too:
  37142 vertices, 74456 triangles, largest owned distance 0.3000 mm, mesher-fable's own numbers (ROWS.md section 1).
- mesher's objdec.c (skel7, 2026-11-08) and body-claude's skindec.c (2026-10-18, its rules): the same collapses and refusals
  folddec started from; nothing further to take.
- experiment/cfd's remesh.c and fvremesh.c (the maintainer): quadric collapses with an RMS plane-distance stop and a fixed boundary
  (--fix-boundary). The fixed boundary is keep = 1 on the rim; the RMS stop is replaced by the distance bound, which is two-sided
  and exact where remesh's is a plane average. Their own copies move at their owners' next release.
- Faults 161-165: the fold refusal, the distance refusals, the crossing test, the link condition, the keep flags left out.

## mk_ovl.h
- mesher-fable's ovl.c (2026-11-15): columns along z through a lattice offset by (0.5371, 0.5213) of a cell, each mesh's
  crossings chained with mk_geom's first hit, the odd lines left out and counted. MEASURED: two icospheres of 10 mm 12 mm apart,
  h 0.2 mm: 0.8706 ml and 3132 lines, as ovl.c; the analytic lens 0.8713 ml.
- stepcad's voldiff.c (the maintainer): its ray-parity volume difference (what the reference has that the model has not, and the
  extra) as only_a, only_b and vol_a, vol_b over both boxes (MK_OVL_ALL), and its choice of axis. Its own crossing test on float
  STL triangles with a 1e-18 degenerate cut is replaced by mk_geom's watertight test; its pictures, blobs and slices stay with it.
- mesher's census (skel7 tools/census, a collaborator job 015): read, not merged. Its overlap volumes are an adaptive integral with
  bounds over an atlas's part pairs, capped in queries; it is the reference ovl was calibrated against (within 3 % on b-helper's
  obliques), not a method to fold in.
- Changed from ovl.c: more than 4096 crossings on a line counts the line as odd (ovl.c stopped at 256 and used what it had); the
  pad below the meshes scales with their size (ovl.c's was 1 cm).
- Faults 171-173: the union for the overlap, lines on the lattice itself, an odd count used as it is.

## mk_graft.h
- mesher-fable's bonegraft.c (2026-11-15): its general parts. Marching tetrahedra (Freudenthal's six a cube, the vertex held 10 %
  off the nodes, an edge hash) as mk_mtet on mk_sdf's grid layout (nodes at lo + h i, where bonegraft's sat at org + (i + 0.5) h);
  the cut of a mesh at a function's zero with rim edges; the patch's psi (S's linear psi near S, psi itself off it, bisected on the
  cut edges); the stitch (each patch rim vertex onto S's cut rim within snap, the rims fanned into each other, welded, the largest
  piece kept). Not taken: the bones (bitten triangles, the removed volume's skin, the opening at BG_OPEN_MM, the Gaussian, the
  weights between the cut field and the muscle's own): those make the field, which is the caller's; the KEEP.txt file becomes
  out->keep. The weld test was 1e-10 m^2 (10 microns); it is 1 % of the rim edge.
- mk_sdf: the test makes the field with mk_sdf_mesh (S's distance, then carved), as the coordinator asked. mk_sdf has no box of
  its own yet, so a local field on a large mesh means a whole-mesh grid; a box option for mk_sdf_mesh is in "Not yet".
- body-claude's tet.c (isosurface stuffing, Labelle and Shewchuk 2007): read, not merged. It meshes a volume into tetrahedra on a
  body-centred lattice with vertices warped onto the surface; it shares the edge-to-cut-vertex hash and the fixed stencils but is
  a volume mesher, not a surface method.
- mk_snet: kept beside mk_mtet. Surface nets give one vertex a cell and smoother surfaces, and can touch themselves where two
  sheets share a cell; marching tetrahedra are always a closed manifold. Other house copies of marching tetrahedra (found by the
  signature): body-claude's seat62b seatcov.c, seat62a thighvox.c, and the wallimp/wallband tools; they move at their owners'
  next release.
- Faults 181-182: every marching-tetrahedra triangle wound inward; the patch's fans left out of the stitch.

# 1.5.0 (2026-11-16): OpenMP threads (the maintainer's device plan, scratch/the maintainer/devflex/PLAN.md item 7; the maintainer's helper)

## mk_par.h
- No house copy merged: a block scheme written for the kit. Read beside it: cudakit's ck_sum_f64 and b580kit's bk_pcg_solve
  (sums in one fixed order, so a GPU solve repeats bit for bit), and grown's OpenMP stages (core and stage/*: plain
  `parallel for`, scripts setting OMP_WAIT_POLICY=passive). The kit takes the fixed order and makes it independent of the thread
  count too: blocks fixed by the item count (4096 items, at most 256 blocks), partials added in block order. `reduction(+:)`
  is never used on floating point (it is planted fault 107).
- The pragmas sit behind `#ifdef _OPENMP` (MK_OMP_), because -Wall's -Wunknown-pragmas would make a bare `#pragma omp` an
  error under -Werror for a caller without -fopenmp.

## mk_sparse.h
- The product, the dots, the scaling, the mean-free shift, the vector updates, Jacobi and block-Jacobi on blocks of rows;
  IC(0)'s triangular solves stay serial (level scheduling would change nothing in the bits but is not worth it without a
  user). CG fuses each product with p.q, the update with r.r and Jacobi with r.z, and with no preconditioner uses r for z
  (r.z is r.r bit for bit): three or four passes an iteration instead of six or seven. BiCGSTAB is threaded step by step.
- A dot over more than one block differs in rounding from 1.4.0's single sum (sizes in the README's 1.5.0 entry); a vector of
  at most 4096 entries is one block and keeps 1.4.0's bits.

## mk_sdf.h
- Every stage keeps 1.4.0's grid bit for bit. The band: a node's value depends on the order its triangles reach it (a tie, or
  a float rounding between two near-equal distances), so tiles of 8 x 8 rows each list their triangles in index order. The
  sweeps: a row reads only the rows before it in the sweep, so any order that keeps those first gives the same values; a
  wavefront of 4 x 4 tiles keeps them first with one barrier a diagonal of tiles (a wavefront of single rows had four times
  as many, and on a loaded machine each costs a time slice).

## mk_geom.h
- The batches (mk_bvh_nearest_many, mk_bvh_inside_many, mk_winding_number_many) are new entry points over the single
  queries; consumers with their own loops (partjudge, bonegraft, femurscore_hc) can call them and keep their own sums
  serial.
