# mathkit is header-only: nothing to build to use it. `make check` builds and runs every test, then every planted fault,
# and prints "wrong: N, skipped: M". `make check OMP=1` builds every test and planted fault with -fopenmp too; the
# threads' checks (t_par with -fopenmp) run either way. The mesh rows need an external closed-mesh checker:
# MESHCHECK=path, else meshcheck on PATH, else they print SKIP. RIBBON_RNG_H=path adds t_rng's bit-for-bit comparison
# with an independent xoshiro256** header. `make bench` times CG and the BVH queries.
CFLAGS = -std=c23 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic -Werror
OMPFLAG = $(if $(filter 1,$(OMP)),-fopenmp)

check:
	CC="$(CC)" CFLAGS="$(CFLAGS)" OMP="$(OMP)" MESHCHECK="$(MESHCHECK)" RIBBON_RNG_H="$(RIBBON_RNG_H)" sh test/check.sh

clean:
	rm -rf build
.PHONY: check clean

bench:
	mkdir -p build && $(CC) $(CFLAGS) -O3 -march=native $(OMPFLAG) test/bench.c -o build/bench -lm && build/bench
.PHONY: bench
