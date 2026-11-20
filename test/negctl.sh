#!/bin/sh
# Negative controls for 1.5.0's rows that no planted fault in the headers reaches: each must FAIL on its known-bad input. Each
# runs in a copy of the tree under build/neg (the tests include ../include, so a changed header needs its own tree). Prints
# one line a control and exits with the count not caught. check.sh runs it; CC and CFLAGS as there (without -fopenmp). An argument
# runs that control alone (test/drvctl.sh).
cd "$(dirname "$0")/.." || exit 2
S=$PWD W=$PWD/build/neg CC=${CC:-gcc} F=${CFLAGS:--std=c23 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic -Werror} only=${1:-}
export OMP_WAIT_POLICY=passive
rm -rf "$W" && mkdir -p "$W" || exit 2
bad=0
want() { [ -z "$only" ] || [ "$only" = "$1" ]; }
tree() { mkdir -p "$W/$1" && cp -r "$S/include" "$S/test" "$W/$1/"; }
say() { # name, 1 if the known-bad input was caught
	if [ "$2" = 1 ]; then printf '%-20s %-28s PASS (fails on: %s)\n' "negative control" "$3" "$1"
	else printf '%-20s %-28s WRONG (not caught: %s)\n' "negative control" "$3" "$1"; bad=$((bad + 1)); fi
}
row() { grep -q "^FAIL.*$2" "$1" 2> /dev/null && echo 1 || echo 0; }
# 1. t_vec's version row against 1.4.0's mk_version.h
if want 1; then
	tree v && printf '#ifndef MK_VERSION_H\n#define MK_VERSION_H\n#define MK_VERSION "1.4.0"\n#define MK_VERSION_MAJOR 1\n#define MK_VERSION_MINOR 4\n#define MK_VERSION_PATCH 0\n#endif\n' > "$W/v/include/mk_version.h" &&
		$CC $F "$W/v/test/t_vec.c" -o "$W/v/t" -lm && "$W/v/t" > "$W/v.out"
	say "1.4.0's mk_version.h" "$(row "$W/v.out" "mathkit version")" "t_vec version"
fi
# 2. CG's block sums against one serial sum, with a block longer than any vector (one block: 1.4.0's sum, the same bits)
if want 2; then
	tree b && sed -i 's/enum { MK_PAR_BLK = 4096, MK_PAR_NBLK = 256 };/enum { MK_PAR_BLK = 1 << 30, MK_PAR_NBLK = 256 };/' "$W/b/include/mk_par.h" &&
		grep -q '1 << 30' "$W/b/include/mk_par.h" && $CC $F "$W/b/test/t_par.c" -o "$W/b/t" -lm && "$W/b/t" > "$W/b.out"
	say "one block (1.4.0's sum)" "$(row "$W/b.out" "block sums against one serial sum")" "t_par block sums"
fi
# 3. the sdf cases' reach: the sphere closed (no odd row) and the box flush with the nodes (no near node)
if want 3; then
	tree r && sed -i 's/box(&bx, (10 + 1e-5) \* 0.05)/box(\&bx, 10 * 0.05)/; s/sphere(&open, 48, 24, 0.97)/sphere(\&open, 48, 24, 2)/' "$W/r/test/t_par.c" &&
		grep -q 'sphere(&open, 48, 24, 2)' "$W/r/test/t_par.c" && $CC $F "$W/r/test/t_par.c" -o "$W/r/t" -lm && "$W/r/t" > "$W/r.out"
	say "closed sphere, flush box" "$(row "$W/r.out" "the sdf cases reach what they test")" "t_par sdf reach"
fi
# 4. check.sh's BITS comparison: a build without OpenMP whose bits differ (planted fault 97 drops a point there too) against the OpenMP
# build. Caught only when both built and ran, printed as many lines, and the lines differ: a build that failed leaves no lines to compare
if want 4; then
	c4=0
	if tree p && $CC $F -DMK_MUTATE=97 "$W/p/test/t_par.c" -o "$W/p/ts" -lm && $CC $F -fopenmp "$W/p/test/t_par.c" -o "$W/p/to" -lm; then
		"$W/p/ts" > "$W/ps.out" 2>&1; "$W/p/to" > "$W/po.out" 2>&1; ro=$?
		grep '^BITS' "$W/ps.out" > "$W/ps.bits"; grep '^BITS' "$W/po.out" > "$W/po.bits"
		if [ "$ro" -eq 0 ] && [ -s "$W/ps.bits" ] && [ "$(wc -l < "$W/ps.bits")" -eq "$(wc -l < "$W/po.bits")" ] && ! cmp -s "$W/ps.bits" "$W/po.bits"; then c4=1; fi
	fi
	say "97 in the build without OpenMP" "$c4" "t_par BITS comparison"
fi
# 5. the batch inside query against its single calls without OpenMP: planted fault 98 at one thread (the last point left out)
if want 5; then
	tree i && $CC $F -DMK_MUTATE=98 "$W/i/test/t_par.c" -o "$W/i/t" -lm && "$W/i/t" > "$W/i.out"
	say "98 without OpenMP" "$(row "$W/i.out" "mk_bvh_inside_many: the single calls")" "t_par inside vs single"
fi
exit "$bad"
