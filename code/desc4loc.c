/* desc4loc.c -- (1) sanity: D == plane on an N x N plane; (2) local check: D satisfies the recurrence
 * on the band |d| <= W and the near-axis strip (x <= 40 or y <= 40) for 2 <= s < S.
 *   ./desc4loc N W S
 */
#include "onesym.h"
#include <omp.h>

int main(int argc, char **argv) {
    long N = atol(argv[1]); int W = atoi(argv[2]); long S = atol(argv[3]);
    init_tables(); init_mainboxes(S + 1000);
    fprintf(stderr, "main boxes generated: %d (last e=%ld level %d)\n", nmb, mb[nmb - 1].e, mb[nmb - 1].lev);
    /* (1) sanity against the plane */
    if (N > 0) {
        uint8_t *pl = malloc((size_t)N * N);
        for (long i = 0; i < N; i++) { pl[i * N] = 1; pl[i] = 1; }
        for (long i = 1; i < N; i++) for (long j = 1; j < N; j++) pl[i * N + j] = F[(pl[(i - 1) * N + j] * 4 + pl[(i - 1) * N + j - 1]) * 4 + pl[i * N + j - 1]];
        long bad = 0, tot = 0;
        #pragma omp parallel for reduction(+:bad,tot) schedule(dynamic, 64)
        for (long i = 0; i < N; i++) for (long j = 0; j < N; j++) { tot++; if (D(i + j, j - i) != pl[i * N + j]) { if (bad < 5) fprintf(stderr, "plane mismatch x=%ld y=%ld\n", i, j); bad++; } }
        printf("plane %ldx%ld: %ld cells, mismatches %ld\n", N, N, tot, bad); fflush(stdout);
        free(pl);
    }
    /* (2) local check */
    long bad = 0, tot = 0;
    #pragma omp parallel for reduction(+:bad,tot) schedule(dynamic, 1024)
    for (long s = 2; s < S; s++) {
        /* cells to check: band |d|<=W, and the strips x<=40 (d >= s-80) and y<=40 (d <= 80-s) */
        long ranges[3][2] = {{-W, W}, {s - 80, s - 2}, {-(s - 2), 80 - s}};
        for (int rg = 0; rg < 3; rg++)
        for (long d = ranges[rg][0]; d <= ranges[rg][1]; d++) {
            if (d < -(s - 2) || d > s - 2) continue;
            if (rg > 0 && d >= -W && d <= W) continue;     /* already covered by the band */
            if (((s - d) & 1) != 0) continue;
            int v = D(s, d), x = D(s - 1, d + 1), y = D(s - 2, d), z = D(s - 1, d - 1);
            tot++;
            if (v != F[(x * 4 + y) * 4 + z]) { if (bad < 10) fprintf(stderr, "LOCAL FAIL s=%ld d=%ld D=%d f=%d\n", s, d, v, F[(x * 4 + y) * 4 + z]); bad++; }
        }
    }
    printf("local check s<%ld, band |d|<=%d + axis strips: %ld cells, failures %ld\n", S, W, tot, bad);
    return 0;
}
