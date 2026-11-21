/* test harness: each CHECK prints PASS or FAIL and a description, each SKIP a row not run; main returns the FAIL
   count */
#include <math.h>
#include <stdio.h>
static int wrong;
#define CHECK(cond, ...) do { int ok_ = (cond); printf(ok_ ? "PASS  " : "FAIL  "); printf(__VA_ARGS__); printf("\n"); fflush(stdout); if (!ok_) wrong++; } while (0) /* flushed: a row printed before a crash still counts */
/* a row that could not run here: check.sh counts it apart, never as a pass */
#define SKIP(...) do { printf("SKIP  "); printf(__VA_ARGS__); printf("\n"); fflush(stdout); } while (0)
#define T_PI 3.14159265358979323846
static unsigned long long t_state = 88172645463325252ULL;
[[maybe_unused]] static double t_rand(void) /* test-only uniform [-1, 1), independent of mk_rng */
{
    t_state ^= t_state << 13, t_state ^= t_state >> 7, t_state ^= t_state << 17;
    return (double)(t_state >> 11) * 0x1.0p-52 - 1.0;
}
