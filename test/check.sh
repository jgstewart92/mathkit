#!/bin/sh
# Each header compiles alone; each test passes on the real header (exit 0); each planted fault (-DMK_MUTATE=N) makes
# its test print a FAIL row (the rules: test/judge.sh). OMP=1 builds every test and every planted fault with -fopenmp
# as well. Either way the threads' own checks run: t_par built with -fopenmp (every threaded routine the same at 1, 2
# and 12 threads), its BITS lines equal to a build without OpenMP, and the planted faults that make an order follow the
# thread count, each of which must fail it. The mesh rows need an external closed-mesh checker: MESHCHECK when set,
# else meshcheck on PATH; with neither, each of them prints SKIP. RIBBON_RNG_H, when set, is the path of an independent
# xoshiro256** header for t_rng to compare with bit for bit; without it that row prints SKIP. Prints a line a check,
# then "wrong: N, skipped: M" (a skip is counted apart and is not a pass); exits with N, at most 100.
cd "$(dirname "$0")/.." || exit 2
B=build
mkdir -p "$B"
wrong=0 skipped=0
OPENMP=-fopenmp CFLAGS0=$CFLAGS
[ -n "${RIBBON_RNG_H:-}" ] && CFLAGS="$CFLAGS -DMK_RIBBON_RNG_H=\"$RIBBON_RNG_H\"" # one word, its quotes kept
[ "${OMP:-0}" = 1 ] && CFLAGS="$CFLAGS $OPENMP"
export OMP_WAIT_POLICY=passive # threads sleep at a barrier: spinning ones crawl on a loaded machine
for h in include/mk_*.h; do
	if printf '#include "%s"\nint main(void) { return 0; }\n' "$h" | $CC $CFLAGS -x c - -o "$B/alone" -lm 2> "$B/alone.err"; then
		printf '%-20s %-28s PASS\n' "$(basename "$h")" "compiles alone"
	else
		printf '%-20s %-28s WRONG\n' "$(basename "$h")" "compiles alone"; cat "$B/alone.err"; wrong=$((wrong + 1))
	fi
done
for h in include/mk_*.h; do # C11 as well: programs built as C11 include the headers
	if printf '#include "%s"\nint main(void) { return 0; }\n' "$h" | $CC $CFLAGS -std=c11 -x c - -o "$B/alone" -lm 2> "$B/alone.err"; then
		printf '%-20s %-28s PASS\n' "$(basename "$h")" "compiles alone, -std=c11"
	else
		printf '%-20s %-28s WRONG\n' "$(basename "$h")" "compiles alone, -std=c11"; cat "$B/alone.err"; wrong=$((wrong + 1))
	fi
done
if [ "${OMP:-0}" != 1 ]; then # the pragmas as well: each header alone with -fopenmp
	for h in include/mk_*.h; do
		if printf '#include "%s"\nint main(void) { return 0; }\n' "$h" | $CC $CFLAGS $OPENMP -x c - -o "$B/alone" -lm 2> "$B/alone.err"; then
			printf '%-20s %-28s PASS\n' "$(basename "$h")" "compiles alone, -fopenmp"
		else
			printf '%-20s %-28s WRONG\n' "$(basename "$h")" "compiles alone, -fopenmp"; cat "$B/alone.err"; wrong=$((wrong + 1))
		fi
	done
fi
CRASHES="" # planted faults caught by a crash rather than a FAIL row, as fault:signal (none)
. test/judge.sh
for spec in "vec 1 2" "fft 11 12 238 239" "hhd 240 241 242 243 244 245" "linalg 21 22 23 24 25 26 27" "rng 31 32 33 34" "opt 41 42 43 44 45" "ode 51 52 53 54 55" "stats 61 62 63 64 65" "interp 71 72" "vpflow 131 132" "snet 121 122 123 124 125" "sdf 141 142 143 144 145 207 208" "sparse 101 102 103 104 105 106" "quad 111 112 113 114" "spec 81 82 83 84" "geom 91 92 93 94 95 206 214" "sepvol 151 152 153 154 155 156 157 158 159 191 192 193 194 195 212 213 215 216 217 218 219 220 221" "dec 161 162 163 164 165 201" "ovl 171 172 173" "graft 181 182 202 203 204 205 209 210" "par 211" "regr 222 223 224 225 226 227 228 229 230 231 232 233 234 235 236 237"; do
	set -- $spec
	t=$1; shift
	run "$t" 0
	for m in "$@"; do run "$t" "$m"; done
done
# the threads: t_par built with -fopenmp; its BITS lines (a hash of every threaded routine's result) equal to those of a
# build without OpenMP; then the planted faults that make an order follow the thread count, each to be caught
runo 0
bits 17
for m in 96 97 98 107 108 146 147 148; do runo $m; done
# The OBJ files t_snet, t_sepvol, t_dec and t_graft write, through the external closed-mesh checker (closed, manifold,
# oriented, Euler characteristic and genus, winding, crossings). The second line of its -t output: field 8 the Euler
# characteristic, 9 the genus, 19 positive for an outward surface, the last the verdict
if [ -n "${MESHCHECK:-}" ]; then MC=$MESHCHECK; else MC=$(command -v meshcheck) || MC=; fi
mcrows() { # the group's name, a prefix for its row names, the folder written, the tests that write it, then the rows
	# as "file chi genus verdict"
	g=$1 p=$2 d=$B/$3 ts=$4; shift 4
	mkdir -p "$d" && rm -f "$d"/*.obj
	if [ -z "$MC" ]; then
		for spec; do
			printf '%-20s %-28s SKIP (no meshcheck: set MESHCHECK or put it on PATH)\n' "meshcheck" "$p${spec%% *}"; skipped=$((skipped + 1))
		done
		return
	fi
	if [ ! -x "$MC" ]; then printf '%-20s %-28s WRONG (no %s)\n' "meshcheck" "$g" "$MC"; wrong=$((wrong + 1)); return; fi
	for t in $ts; do
		if ! "$B/t_${t}_0" "$d" > /dev/null; then printf '%-20s %-28s WRONG (t_%s did not write them)\n' "meshcheck" "$g" "$t"; wrong=$((wrong + 1)); return; fi
	done
	for spec; do
		set -- $spec
		got=$("$MC" -t "$d/$1.obj" | awk -F'\t' 'NR == 2 { print $8, $9, ($19 > 0 ? "out" : "in"), $NF }')
		if [ "$got" = "$2 $3 out $4" ]; then printf '%-20s %-28s PASS (chi %s, genus %s, outward, %s)\n' "meshcheck" "$p$1" "$2" "$3" "$4"
		else printf '%-20s %-28s WRONG (read: %s; want chi %s genus %s out %s)\n' "meshcheck" "$p$1" "$got" "$2" "$3" "$4"; wrong=$((wrong + 1)); fi
	done
}
# mk_snet's surfaces; the two voxels that share only an edge must read FAULT (the negative control)
mcrows "snet surfaces" "" snetobj snet "ball_raw 2 0 clean" "ball_relaxed 2 0 clean" "ball_iso 2 0 clean" "half1 2 0 clean" \
	"half2 2 0 clean" "torus 0 1 clean" "diagonal 3 -1 FAULT"
# mk_sepvol's results: each separated pair clean (no crossing faces, winding 1 just inside and 0 just outside every
# face); the pairs before separation must read FAULT (the negative controls)
mcrows "sepvol results" "sepvol " sepvolobj sepvol "twoball 4 0 clean" "pair 4 0 clean" "cubes 4 0 clean" "three 6 0 clean" \
	"crumpled 4 0 clean" "coarse 4 0 clean" "plate 4 0 clean" "cubes_before 4 0 FAULT" "crumpled_before 4 0 FAULT" \
	"coarse_before 4 0 FAULT" "plate_before 4 0 FAULT"
# mk_decimate's and the graft's results (the hollow ball: two shells, the inner one inward)
mcrows "dec and graft results" "" decobj "dec graft" "dec_sphere 2 0 clean" "dec_torus 0 1 clean" "dec_hollow 4 0 clean" \
	"mtet_sphere 2 0 clean" "mtet_torus 0 1 clean" "graft 2 0 clean"
# test/negctl.sh: the negative controls of the rows no planted fault in the headers reaches, each to fail on its
# known-bad input; test/drvctl.sh: the verdicts of test/judge.sh itself, each on a known-bad input
CC="$CC" CFLAGS="$CFLAGS0" sh test/negctl.sh; wrong=$((wrong + $?))
CC="$CC" CFLAGS="$CFLAGS0" sh test/drvctl.sh; wrong=$((wrong + $?))
echo "wrong: $wrong, skipped: $skipped"
exit "$(clamp "$wrong")"
