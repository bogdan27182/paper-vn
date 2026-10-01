/* line1.c -- facts on the line d = 1 (y = x + 1) used for non-automaticity.
 * t (three letters): t[s,1] = 0  iff  s = 3 or s is an edge time of a box of the family B.
 * tau (four letters): tau[s,1] = 3 iff s = e_B + l_{lev B}/2 + 5 for a box B of the family C of level >= 1
 *                     (for s >= 56), where the main pairs of level 0 count through the middle block.
 * Also prints the record gaps (first occurrence of each new maximal run without the marker letter).
 */
#include "onesym.h"
#include "desc.h"
#define SMAX3 (1L << 21)
static uint8_t *edge3;                 /* edge times of the three-letter family */
static void mark3(int n, long e) {     /* box of level n opening at e */
    if (e >= SMAX3) return;
    edge3[e] = 1; if (e + Lbox(n) < SMAX3) edge3[e + Lbox(n)] = 1;
    if (n >= 3) { int lev[128]; long tim[128]; int c = children(n, e, lev, tim); for (int i = 0; i < c; i++) mark3(lev[i], tim[i]); }
}
static uint8_t *green4; static long SMAX4;
static void mark4(int n, long e) {     /* box of level n >= 1 of the family C opening at e */
    if (e >= SMAX4) return;
    long g = e + Lb[n] / 2 + 5; if (g < SMAX4) green4[g] = 1;
    for (int i = 0; i < chn[n]; i++) if (chlev[n][i] >= 1) mark4(chlev[n][i], e + choff[n][i]);
}
int main(int argc, char **argv) {
    SMAX4 = argc > 1 ? atol(argv[1]) : 50000000L;
    /* ---------- three letters: recurrence on x,y < N3 ---------- */
    long N3 = 4096; uint8_t *t = malloc(N3 * N3);
    for (long i = 0; i < N3; i++) { t[i * N3] = 1; t[i] = 1; }
    for (long x = 1; x < N3; x++) for (long y = 1; y < N3; y++) t[x * N3 + y] = fF3(t[(x - 1) * N3 + y], t[(x - 1) * N3 + y - 1], t[x * N3 + y - 1]);
    edge3 = calloc(SMAX3, 1);
    for (int n = 1; Emain(n) < SMAX3; n++) mark3(n, Emain(n));
    long bad3 = 0;
    for (long x = 1; x + 1 < N3; x++) { long s = 2 * x + 1; int white = t[x * N3 + x + 1] == 0; int pred = (s == 3) || edge3[s]; if (white != pred) { if (bad3 < 5) printf("t MISMATCH x=%ld s=%ld white=%d pred=%d\n", x, s, white, pred); bad3++; } }
    printf("three letters: line y=x+1, x < %ld: mismatches with 'white iff edge time or s=3': %ld\n", N3 - 1, bad3);
    /* record gaps of the three-letter edge set along odd s (in x units, run of non-white cells) */
    printf("  record runs of non-white (x = start of run, length in cells):");
    long last = 1, best = 0;
    for (long s = 5; s < SMAX3; s += 2) if (edge3[s]) { long run = (s - last) / 2 - 1; if (run > best) { best = run; printf(" [x=%ld: %ld]", (last + 1) / 2, run); } last = s; }
    printf("\n");
    /* ---------- four letters ---------- */
    init_tables(); init_mainboxes(SMAX4 + 1000);
    green4 = calloc(SMAX4, 1);
    for (int q = 0; q < nmb; q++) if (mb[q].lev >= 1) mark4(mb[q].lev, mb[q].e);
    for (int q = 1; q + 1 < nmb; q++) if (mb[q].lev == 0 && mb[q + 1].lev == 0 && mb[q - 1].lev == 1) { long g = mb[q - 1].e + Lb[1] + 12; if (g < SMAX4) green4[g] = 1; }
    /* against the recurrence on x,y < N4 */
    long N4 = 12000; uint8_t *u = malloc(N4 * N4);
    for (long i = 0; i < N4; i++) { u[i * N4] = 1; u[i] = 1; }
    for (long x = 1; x < N4; x++) for (long y = 1; y < N4; y++) u[x * N4 + y] = F[(u[(x - 1) * N4 + y] * 4 + u[(x - 1) * N4 + y - 1]) * 4 + u[x * N4 + y - 1]];
    long bad4 = 0, bad4D = 0;
    for (long x = 28; x + 1 < N4; x++) { long s = 2 * x + 1; int g = u[x * N4 + x + 1] == 3; if (g != green4[s]) { if (bad4 < 5) printf("tau MISMATCH x=%ld s=%ld green=%d pred=%d\n", x, s, g, green4[s]); bad4++; } }
    for (long s = 57; s < SMAX4; s += 2) { int g = D(s, 1) == 3; if (g != green4[s]) { if (bad4D < 5) printf("D MISMATCH s=%ld\n", s); bad4D++; } }
    printf("four letters: line y=x+1: mismatches with the green set: recurrence (x<%ld) %ld, description (s<%ld) %ld\n", N4 - 1, bad4, SMAX4, bad4D);
    printf("  record runs without green (x = start of run, length in cells):");
    last = 55; best = 0;
    for (long s = 57; s < SMAX4; s += 2) if (green4[s]) { long run = (s - last) / 2 - 1; if (run > best) { best = run; printf(" [x=%ld: %ld]", (last + 1) / 2, run); } last = s; }
    printf("\n");
    return 0;
}
