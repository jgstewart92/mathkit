# What check.sh counts as passed and as caught. Sourced by check.sh and by test/drvctl.sh, which shows each verdict on a known-bad
# input. Needs B, CC, CFLAGS; adds to wrong.
#   The real header passes when its test runs to its end with at least one PASS row and no FAIL row.
#   A planted fault is caught by a FAIL row. A timeout is never a catch, nor is a crash or a bare nonzero exit, with one exception:
#   a fault meant to be caught by a crash is named in CRASHES as fault:signal, and only that signal counts. Exit 124 or 137 is read
#   as a timeout even when it is the FAIL count (124 or 137 rows: a loud WRONG, never a pass; the largest test has 73 rows, 1.6.2).
LIMIT=${LIMIT:-120} # seconds a test may run
wrong=${wrong:-0}
say_() { printf '%-20s %-28s %s\n' "$1" "$2" "$3"; }
bad_() { say_ "$1" "$2" "WRONG ($3)"; wrong=$((wrong + 1)); }
# nothing an earlier run left behind may count for this one
fresh() { rm -f "$1" "$1.out" "$1.err"; }
# run a built test: r its exit status, fails and passes its rows, ended how it ended (ok, rows: it returned its FAIL count, timeout,
# signal:N, exit:N). A test returns its FAIL count, which wraps at 256.
runexe() {
	timeout -k 5 "$LIMIT" "$1" > "$1.out" 2>&1; r=$?
	fails=$(grep -c '^FAIL' "$1.out"); passes=$(grep -c '^PASS' "$1.out")
	if [ "$r" -eq 0 ]; then ended=ok
	elif [ "$r" -eq 124 ] || [ "$r" -eq 137 ]; then ended=timeout # 1.6.0: ahead of rows, so 124 FAIL rows and a hang are a timeout (the review's J1)
	elif [ $((fails % 256)) -eq "$r" ]; then ended=rows
	elif [ "$r" -gt 128 ]; then ended=signal:$((r - 128))
	else ended=exit:$r; fi
}
# the verdict on the run just made: name, mutant (0 = the real header)
judge() {
	if [ "$2" = 0 ]; then
		if [ "$ended" = ok ] && [ "$fails" -eq 0 ] && [ "$passes" -gt 0 ]; then say_ "$1" "real header" "PASS ($passes rows)"
		else bad_ "$1" "real header" "$ended, $fails rows fail, $passes pass"; grep '^FAIL' "$exe.out"; fi
	elif [ "$ended" = timeout ]; then bad_ "$1" "planted fault $2" "timeout after $LIMIT s: not a catch"
	elif [ "$fails" -gt 0 ]; then say_ "$1" "planted fault $2" "PASS (caught by $fails rows)"
	else
		case " $CRASHES " in
		*" $2:${ended#signal:} "*) say_ "$1" "planted fault $2" "PASS (caught by signal ${ended#signal:}, listed)" ;;
		*) bad_ "$1" "planted fault $2" "not caught: $ended, 0 FAIL rows" ;;
		esac
	fi
}
# build and run: test, mutant, name for the verdict, extra flags
build_run() {
	if [ "$2" = 0 ]; then def=""; else def="-DMK_MUTATE=$2"; fi
	fresh "$exe"
	if ! $CC $CFLAGS $4 $def "test/t_$1.c" -o "$exe" -lm 2> "$exe.err"; then
		bad_ "$3" "mutant $2" "build"; cat "$exe.err"; ended=build; return
	fi
	runexe "$exe"; judge "$3" "$2"
}
run() { exe="$B/t_$1_$2"; build_run "$1" "$2" "t_$1" ""; } # test, mutant (0 = none)
runo() { exe="$B/t_par_omp_$1"; build_run par "$1" "t_par -fopenmp" "$OPENMP"; [ "$1" = 0 ] && omp0=$ended; } # mutant (0 = none)
# t_par's BITS lines without OpenMP against those of runo 0 (which ran just before): both built and run through now, as many lines as
# the test prints (want), equal
bits() { # want
	exe="$B/t_par_serial"; fresh "$exe"; rm -f "$B/bits.serial" "$B/bits.omp"
	if ! $CC $CFLAGS0 test/t_par.c -o "$exe" -lm 2> "$exe.err"; then bad_ "t_par BITS" "no OpenMP = OpenMP" "the build without OpenMP"; cat "$exe.err"; return; fi
	runexe "$exe"
	if [ "$ended" != ok ]; then bad_ "t_par BITS" "no OpenMP = OpenMP" "the run without OpenMP: $ended"; return; fi
	if [ "${omp0:-}" != ok ]; then bad_ "t_par BITS" "no OpenMP = OpenMP" "the OpenMP run: ${omp0:-not run}"; return; fi
	grep '^BITS' "$exe.out" > "$B/bits.serial"; grep '^BITS' "$B/t_par_omp_0.out" > "$B/bits.omp"
	nb=$(wc -l < "$B/bits.serial")
	if [ "$nb" -ne "$1" ]; then bad_ "t_par BITS" "no OpenMP = OpenMP" "$nb lines, the test prints $1"
	elif cmp -s "$B/bits.serial" "$B/bits.omp"; then say_ "t_par BITS" "no OpenMP = OpenMP" "PASS ($nb results)"
	else bad_ "t_par BITS" "no OpenMP = OpenMP" "they differ"; diff "$B/bits.serial" "$B/bits.omp"; fi
}
# an exit status that cannot wrap to 0
clamp() { if [ "$1" -gt 100 ]; then echo 100; else echo "$1"; fi; }
