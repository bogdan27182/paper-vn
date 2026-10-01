/* verify.c
 *   ./verify grid N    : compute t(i,j) for 0<=i,j<N by the recurrence (row by row) and compare
 *                        with the explicit description D; collect statistics of argument triples.
 *   ./verify local n   : check that D satisfies the recurrence on every cell of the band up to the
 *                        end of the main ordinal S_n (no use of t at all).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "desc.h"

static uint8_t F[3][3][3];

static void grid(long N) {
    long S = 2 * N;                       /* rows s < 2N cover all (i,j) with i,j<N */
    int nmax = 1;
    while (Emain(nmax + 1) < S) nmax++;
    long W = Wbox(nmax + 1) + 8;
    Band B = paint(S, W);
    fprintf(stderr, "band painted: S=%ld W=%ld (main ordinals up to S_%d started)\n", S, W, nmax);
    uint8_t *prev = malloc(N), *cur = malloc(N);
    memset(prev, 1, N);
    unsigned long long cnt[4][27]; long first[4][27][2];
    memset(cnt, 0, sizeof cnt);
    for (int c = 0; c < 4; c++) for (int k = 0; k < 27; k++) first[c][k][0] = -1;
    unsigned long long bad = 0, checked = 0;
    /* row 0 */
    for (long j = 0; j < N; j++) { checked++; if (Dget(&B, j, j) != 1) bad++; }
    for (long i = 1; i < N; i++) {
        cur[0] = 1;
        checked++; if (Dget(&B, i, -i) != 1) bad++;
        for (long j = 1; j < N; j++) {
            int x = prev[j], y = prev[j - 1], z = cur[j - 1];
            uint8_t v = F[x][y][z];
            cur[j] = v;
            long s = i + j, d = j - i;
            checked++;
            if (Dget(&B, s, d) != v) {
                if (bad < 20) fprintf(stderr, "MISMATCH i=%ld j=%ld s=%ld d=%ld t=%d D=%d\n", i, j, s, d, v, Dget(&B, s, d));
                bad++;
            }
            int k = 9 * x + 3 * y + z;
            /* categories: 0 all cells, 1 i,j>=3 (not using rows/cols 0,1,2 as arguments... i.e. args have i-1,j-1>=2),
               2 inside band |d|<=W, 3 outside band */
            int cats[4] = {1, (i >= 3 && j >= 3), labs_(d) <= W, labs_(d) > W};
            for (int c = 0; c < 4; c++) if (cats[c]) {
                if (!cnt[c][k]) { first[c][k][0] = i; first[c][k][1] = j; }
                cnt[c][k]++;
            }
        }
        uint8_t *tmp = prev; prev = cur; cur = tmp;
    }
    printf("grid N=%ld: checked %llu cells, mismatches %llu\n", N, checked, bad);
    const char *cn[4] = {"all i,j>=1", "i,j>=3", "band |d|<=W", "outside band"};
    for (int c = 0; c < 4; c++) {
        int nt = 0;
        for (int k = 0; k < 27; k++) if (cnt[c][k]) nt++;
        printf("\n[%s] %d distinct triples (x,y,z)=(t(i-1,j),t(i-1,j-1),t(i,j-1)):\n", cn[c], nt);
        for (int k = 0; k < 27; k++) {
            int x = k / 9, y = (k / 3) % 3, z = k % 3;
            if (cnt[c][k])
                printf("  (%d,%d,%d) -> f=%d   count %15llu   first at (i,j)=(%ld,%ld)\n", x, y, z, F[x][y][z], cnt[c][k], first[c][k][0], first[c][k][1]);
            else
                printf("  (%d,%d,%d) -> f=%d   never\n", x, y, z, F[x][y][z]);
        }
    }
    free(B.a); free(prev); free(cur);
}

static void local(int nmax) {
    long S = Emain(nmax) + Lbox(nmax) + 12;
    long W = Wbox(nmax + 1) + 8;
    Band B = paint(S, W);
    fprintf(stderr, "band painted: S=%ld W=%ld\n", S, W);
    unsigned long long bad = 0, checked = 0;
    for (long s = 2; s < S; s++) {
        long h = s - 2 < W - 1 ? s - 2 : W - 1;
        for (long d = -h; d <= h; d++) {
            if (((s - d) & 1) != 0) continue;
            uint8_t v = F[Dget(&B, s - 1, d + 1)][Dget(&B, s - 2, d)][Dget(&B, s - 1, d - 1)];
            checked++;
            if (v != Dget(&B, s, d)) {
                if (bad < 20) fprintf(stderr, "LOCAL FAIL s=%ld d=%ld D=%d f(..)=%d\n", s, d, Dget(&B, s, d), v);
                bad++;
            }
        }
    }
    /* cells outside the band: all four cells are given by Dfar; check a strip just outside the band */
    for (long s = 2; s < S; s++)
        for (long d = W; d <= W + 6 && d <= s - 2; d++) {
            if (((s - d) & 1) != 0) continue;
            for (int sg = -1; sg <= 1; sg += 2) {
                long dd = sg * d;
                uint8_t v = F[Dget(&B, s - 1, dd + 1)][Dget(&B, s - 2, dd)][Dget(&B, s - 1, dd - 1)];
                checked++;
                if (v != Dget(&B, s, dd)) { if (bad < 20) fprintf(stderr, "FAR FAIL s=%ld d=%ld\n", s, dd); bad++; }
            }
        }
    printf("local check up to S_%d (s < %ld): %llu cells, failures %llu\n", nmax, S, checked, bad);
    free(B.a);
}

int main(int argc, char **argv) {
    for (int x = 0; x < 3; x++) for (int y = 0; y < 3; y++) for (int z = 0; z < 3; z++) F[x][y][z] = fF3(x, y, z);
    if (argc < 3) { fprintf(stderr, "usage: %s grid N | local n\n", argv[0]); return 1; }
    if (!strcmp(argv[1], "grid")) grid(atol(argv[2]));
    else local(atoi(argv[2]));
    return 0;
}
