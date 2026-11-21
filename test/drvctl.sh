#!/bin/sh
# Negative controls for test/judge.sh, the rules behind check.sh's "wrong: 0": each known-bad input must get the verdict
# named here, on small programs in build/drv. Prints one line a control and exits with the count that did not. check.sh
# runs it; CC and CFLAGS as there (without -fopenmp).
cd "$(dirname "$0")/.." || exit 2
S=$PWD W=$PWD/build/drv CC=${CC:-gcc} CFLAGS=${CFLAGS:--std=c23 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic -Werror}
rm -rf "$W" && mkdir -p "$W/test" "$W/build" || exit 2
miss=0
say() { # what, the verdict wanted (PASS or WRONG), the output to read
	if grep -q "$2" "$3"; then printf '%-20s %-28s PASS (%s reads %s)\n' "driver control" "$4" "$1" "$(printf '%s' "$2" | tr A-Z a-z)"
	else printf '%-20s %-28s WRONG (%s should read %s): %s\n' "driver control" "$4" "$1" "$2" "$(tail -n 1 "$3")"; miss=$((miss + 1)); fi
}
# a test program: its main's body for the real header and for fault 901
prog() { printf '#include <stdio.h>\n#include <stdlib.h>\n#include <signal.h>\n#include <unistd.h>\nint main(void)\n{\n#if defined(MK_MUTATE)\n%s\n#else\n%s\n#endif\n}\n' "$2" "$1" > "$W/test/t_x.c"; }
good='puts("PASS  a row"); return 0;'
# one judged run in the small tree, in a subshell: fault (0 = the real header), CRASHES, LIMIT
try() { (cd "$W" && B=build CFLAGS0=$CFLAGS CRASHES=$2 LIMIT=${3:-5} && . "$S/test/judge.sh" && run x "$1") > "$W/out" 2>&1; }

prog "$good" 'puts("PASS  a row"); return 86;'
try 901; say "a fault that exits 86 with no FAIL row" "WRONG (not caught: exit:86" "$W/out" "exit 86, no row"
prog "$good" 'puts("PASS  a row"); fflush(stdout); for (;;) pause();'
try 901 "" 1; say "a fault that hangs" "WRONG (timeout" "$W/out" "a hang"
prog "$good" 'puts("FAIL  a row"); fflush(stdout); for (;;) pause();'
try 901 "" 1; say "a hang after a FAIL row" "WRONG (timeout" "$W/out" "a hang after a row"
prog "$good" 'for (int i = 0; i < 124; i++) { puts("FAIL  a row"); } fflush(stdout); for (;;) pause();' # timeout's exit 124 equal to the row count
try 901 "" 1; say "124 FAIL rows, then a hang" "WRONG (timeout" "$W/out" "a hang after 124 rows"
prog "$good" 'puts("PASS  a row"); fflush(stdout); abort();'
try 901; say "a crash not on the list" "WRONG (not caught: signal:6" "$W/out" "a crash, not listed"
try 901 "901:6"; say "a crash on the list" "PASS (caught by signal 6" "$W/out" "a crash, listed"
try 901 "901:11"; say "a crash by another signal than listed" "WRONG (not caught: signal:6" "$W/out" "the wrong signal"
prog "$good" 'puts("FAIL  a row"); return 1;'
try 901; say "a fault with a FAIL row" "PASS (caught by 1 rows" "$W/out" "a FAIL row"
prog "$good" 'for (int i = 0; i < 256; i++) { puts("FAIL  a row"); } return 0;'
try 901; say "256 FAIL rows, exit 0" "PASS (caught by 256 rows" "$W/out" "the count wrapped"
prog 'for (int i = 0; i < 256; i++) { puts("FAIL  a row"); } return 0;' "$good"
try 0; say "the real header, 256 FAIL rows, exit 0" "WRONG (ok, 256 rows fail" "$W/out" "real header, wrapped"
prog 'return 0;' "$good"
try 0; say "the real header, no rows at all" "WRONG (ok, 0 rows fail, 0 pass" "$W/out" "real header, no rows"
prog 'puts("PASS  a row"); return 86;' "$good"
try 0; say "the real header, exit 86" "WRONG (exit:86" "$W/out" "real header, exit 86"
prog "$good" "$good"
try 0; say "the real header, one PASS row" "PASS (1 rows" "$W/out" "real header, clean"
prog 'puts("PASS  a row"); puts("SKIP  a row"); return 0;' "$good"
try 0; say "the real header, a PASS and a SKIP row" "PASS (1 rows, 1 skipped" "$W/out" "real header, a skip"
prog 'puts("SKIP  a row"); return 0;' "$good"
try 0; say "the real header, only a SKIP row" "WRONG (ok, 0 rows fail, 0 pass, 1 skipped" "$W/out" "real header, skips only"

# t_par's comparison of the builds, with the outputs of an earlier good run left in build/
par() { # the program without OpenMP, the BITS lines wanted
	printf '#include <stdio.h>\nint main(void)\n{\n%s\n}\n' "$1" > "$W/test/t_par.c"
	(cd "$W" && B=build OPENMP=-DOMPBUILD CFLAGS0=$CFLAGS LIMIT=5 && . "$S/test/judge.sh" && { printf 'PASS  x\nBITS  a 1\nBITS  b 2\n' > build/t_par_serial.out
		printf 'BITS  a 1\nBITS  b 2\n' > build/bits.serial; cp build/bits.serial build/bits.omp; omp0=${3:-ok}; printf 'PASS  x\nBITS  a 1\nBITS  b 2\n' > build/t_par_omp_0.out; bits "$2"; }) > "$W/out" 2>&1
}
twobits='puts("PASS  x"); puts("BITS  a 1"); puts("BITS  b 2"); return 0;'
par "$twobits" 2; say "the same lines" "PASS (2 results" "$W/out" "BITS, clean"
par 'not C' 2; say "t_par not compiling, stale outputs there" "WRONG (the build without OpenMP" "$W/out" "BITS, no build"
par 'puts("BITS  a 1"); puts("BITS  b 2"); return 3;' 2; say "t_par failing, stale outputs there" "WRONG (the run without OpenMP" "$W/out" "BITS, failed run"
par "$twobits" 2 build; say "the OpenMP build failed, stale outputs there" "WRONG (the OpenMP run: build" "$W/out" "BITS, OpenMP not built"
par "$twobits" 3; say "fewer lines than the test prints" "WRONG (2 lines, the test prints 3" "$W/out" "BITS, count"
par 'puts("PASS  x"); puts("BITS  a 1"); puts("BITS  b 3"); return 0;' 2; say "different lines" "WRONG (they differ" "$W/out" "BITS, differ"

# negctl's BITS control when the build without OpenMP cannot be made, and when it is made but prints no BITS line (an
# empty list must not read as "they differ"): neither may read as caught
printf '#!/bin/sh\nfor a; do [ "$a" = -DMK_MUTATE=97 ] && exit 1; done\nexec %s "$@"\n' "$CC" > "$W/cc97" && chmod +x "$W/cc97"
CC="$W/cc97" CFLAGS="$CFLAGS" sh "$S/test/negctl.sh" 4 > "$W/out" 2>&1
say "negctl's control 4 with its build failing" "WRONG (not caught" "$W/out" "negctl, failed build"
printf '#!/bin/sh\no= p=\nfor a; do [ "$p" = -o ] && o=$a; p=$a; done\ncase " $* " in *" -DMK_MUTATE=97 "*) printf "int main(void) { return 1; }\\n" | %s -x c - -o "$o"; exit;; esac\nexec %s "$@"\n' "$CC" "$CC" > "$W/cc97e" && chmod +x "$W/cc97e"
CC="$W/cc97e" CFLAGS="$CFLAGS" sh "$S/test/negctl.sh" 4 > "$W/out" 2>&1
say "negctl's control 4, its build printing no BITS" "WRONG (not caught" "$W/out" "negctl, no BITS lines"

# the exit status
(. "$S/test/judge.sh"; exit "$(clamp 256)"); c=$?
if [ "$c" -ne 0 ]; then printf '%-20s %-28s PASS (256 wrong exits %d)\n' "driver control" "exit status" "$c"
else printf '%-20s %-28s WRONG (256 wrong exits 0)\n' "driver control" "exit status"; miss=$((miss + 1)); fi
exit "$miss"
