# mathkit: header-only, nothing to build for users. `make check` builds every test, then every planted fault, and prints "wrong: N";
# the threads' checks (t_par with -fopenmp) run either way. `make check OMP=1` builds every test and planted fault with -fopenmp too.
CFLAGS = -std=c23 -O2 -ffp-contract=off -Wall -Wextra -Wpedantic -Werror
OMPFLAG = $(if $(filter 1,$(OMP)),-fopenmp)

check:
	CC="$(CC)" CFLAGS="$(CFLAGS)" OMP="$(OMP)" sh test/check.sh

clean:
	rm -rf build
.PHONY: check clean

bench:
	mkdir -p build && $(CC) $(CFLAGS) -O3 -march=native $(OMPFLAG) test/bench.c -o build/bench -lm && build/bench
.PHONY: bench
