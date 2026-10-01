/* band.c -- simulate the 4-symbol recurrence on the band |d| <= W of the rotated plane, taking the
 * cells at |d| = W+1 from the periodic environment lattice Lambda (period 6 in s, 12 in |d|).
 * Rows s = 0..S-1.  Then list the main boxes (white-bordered, colour 0) with their levels.
 *   ./band S W [start_rows]
 * If start_rows > 0, rows 0..start_rows-1 are taken from a full-plane computation (so that the
 * transient near the origin is exact); the boundary cells |d|=W+1 are Lambda from then on.
 * Prints: main boxes as "e w L level gap".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static const char *TAB = "1023022312311302010120102300301020002120300301123102300010332020";
static uint8_t F[64];
static const int LAM[6][12] = {
    {3,9,2,9,2,9,2,9,2,9,3,9}, {9,0,9,2,9,0,9,2,9,0,9,3}, {2,9,3,9,3,9,3,9,3,9,2,9},
    {9,0,9,3,9,2,9,3,9,0,9,1}, {2,9,2,9,2,9,2,9,2,9,2,9}, {9,0,9,1,9,0,9,1,9,0,9,2}};
static inline int lam(long s, long d) { long a = d < 0 ? -d : d; return LAM[s % 6][a % 12]; }

typedef struct { long e; int w; long L; int lev, parent; } Box;
#define MAXB 30000000

int main(int argc, char **argv) {
    long S = atol(argv[1]); int W = atoi(argv[2]); long start = argc > 3 ? atol(argv[3]) : 0;
    for (int i = 0; i < 64; i++) F[i] = TAB[i] - '0';
    int Wd = 2 * (W + 1) + 1;
    uint8_t *t = malloc((size_t)S * Wd);
    #define C(s, d) t[(s) * (long)Wd + (d) + W + 1]
    /* exact start: compute the triangle |d|<=s for s < start via the full recurrence (needs a full plane
       of size start/2; fine for start up to ~40000) */
    long N0 = (start + W + 4) / 2 + 2;
    uint8_t *full = NULL;
    if (start > 0) {
        full = malloc((size_t)N0 * N0);
        for (long i = 0; i < N0; i++) { full[i * N0] = 1; full[i] = 1; }
        for (long i = 1; i < N0; i++) for (long j = 1; j < N0; j++)
            full[i * N0 + j] = F[(full[(i - 1) * N0 + j] * 4 + full[(i - 1) * N0 + j - 1]) * 4 + full[i * N0 + j - 1]];
    }
    for (long s = 0; s < S; s++) {
        for (int d = -(W + 1); d <= W + 1; d++) {
            if (((s - d) & 1) != 0) continue;
            long ad = d < 0 ? -d : d;
            if (s < start && ad <= s) { C(s, d) = full[((s - d) / 2) * N0 + (s + d) / 2]; continue; }
            if (ad > s) { C(s, d) = 9; continue; }
            if (ad == s) { C(s, d) = 1; continue; }
            if (ad == W + 1 || s < 2) { C(s, d) = lam(s, d); continue; }
            C(s, d) = F[(C(s - 1, d + 1) * 4 + C(s - 2, d)) * 4 + C(s - 1, d - 1)];
        }
    }
    /* consistency check of the boundary: cells at |d| = W-1..W should equal Lambda too (report count) */
    long incons = 0;
    for (long s = start; s < S; s++) for (int d = W - 3; d <= W; d++) if (((s - d) & 1) == 0 && C(s, d) != lam(s, d)) incons++;
    fprintf(stderr, "boundary inconsistencies (|d| in [W-3,W] != Lambda): %ld\n", incons);
    /* boxes of colour 0 */
    Box *bx = malloc(sizeof(Box) * MAXB); long nb = 0;
    for (long e = 1; e < S; e++)
        for (int w = 1; w <= W - 2 && w <= e; w++) {
            if (((e - w) & 1) != 0) continue;
            int ok = 1;
            for (int d = -w; d <= w; d += 2) if (C(e, d) != 0) { ok = 0; break; }
            if (!ok) continue;
            for (long L = 2; e + L < S; L += 2) {
                long r = e + L - 1; int has = 0;
                for (int d = -(w - 1); d <= w - 1; d += 2) if (C(r, d) != 0) { has = 1; break; }
                if (!has) break;
                if (L >= 4) { r = e + L - 2; has = 0; for (int d = -w; d <= w; d += 2) if (C(r, d) != 0) { has = 1; break; } if (!has) break; }
                if (C(e + L, w) != 0 || C(e + L, -w) != 0) break;
                int closed = 1;
                for (int d = -w + 2; d <= w - 2; d += 2) if (C(e + L, d) != 0) { closed = 0; break; }
                if (closed) { if (nb >= MAXB) { fprintf(stderr, "too many boxes\n"); return 1; } bx[nb].e = e; bx[nb].w = w; bx[nb].L = L; bx[nb].lev = -2; bx[nb].parent = -1; nb++; }
            }
        }
    /* parents: boxes sorted by (e, -w); stack of open boxes (boxes are nested or disjoint) */
    for (long i = 0; i < nb; i++) bx[i].parent = -1;
    {
        long *st = malloc(sizeof(long) * (nb + 1)); long sp = 0;
        /* find_boxes emits boxes in increasing e, and for equal e in increasing w: reorder equal-e runs by decreasing w */
        for (long i = 0; i < nb; ) {
            long j2 = i; while (j2 + 1 < nb && bx[j2 + 1].e == bx[i].e) j2++;
            for (long a = i, b = j2; a < b; a++, b--) { Box tmp = bx[a]; bx[a] = bx[b]; bx[b] = tmp; }
            i = j2 + 1;
        }
        for (long i = 0; i < nb; i++) {
            while (sp > 0 && bx[st[sp - 1]].e + bx[st[sp - 1]].L <= bx[i].e) sp--;
            if (sp > 0 && bx[st[sp - 1]].w > bx[i].w && bx[i].e + bx[i].L < bx[st[sp - 1]].e + bx[st[sp - 1]].L) bx[i].parent = st[sp - 1];
            st[sp++] = i;
        }
        free(st);
    }
    /* levels: process in decreasing e (children have larger e than parents) */
    long *kids = malloc(sizeof(long) * nb); long *kc = calloc(nb, sizeof(long));
    for (long i = 0; i < nb; i++) if (bx[i].parent >= 0) kc[bx[i].parent]++;
    long *ko = malloc(sizeof(long) * (nb + 1)); ko[0] = 0; for (long i = 0; i < nb; i++) ko[i + 1] = ko[i] + kc[i];
    long *fill = calloc(nb, sizeof(long));
    for (long i = 0; i < nb; i++) if (bx[i].parent >= 0) { long p = bx[i].parent; kids[ko[p] + fill[p]++] = i; }
    for (long i = nb - 1; i >= 0; i--) {
        long n = kc[i];
        if (n == 0) { bx[i].lev = 0; continue; }
        if (n & 1) { bx[i].lev = -1; continue; }
        long h = n / 2; int ok = 1;
        for (long m = 0; m < h; m++) { if (bx[kids[ko[i] + m]].lev != h - 1 - m || bx[kids[ko[i] + n - 1 - m]].lev != h - 1 - m) { ok = 0; break; } }
        bx[i].lev = ok ? (int)h : -1;
    }
    long prev_end = -1;
    for (long i = 0; i < nb; i++) if (bx[i].parent < 0) {
        printf("%ld %d %ld %d %ld\n", bx[i].e, bx[i].w, bx[i].L, bx[i].lev, prev_end < 0 ? -1 : bx[i].e - prev_end);
        prev_end = bx[i].e + bx[i].L;
    }
    return 0;
}
