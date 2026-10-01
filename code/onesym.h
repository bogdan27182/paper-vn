/* onesym.h -- standalone explicit description D(s,d) of the one-symbol plane (q=4 example)
 * in rotated coordinates s = x+y, d = y-x.  Tables (P0, P1, MID, axis strip, corner, initial
 * segment) are filled from a small plane (N=200) computed by the recurrence: they are finite data.
 */
#ifndef ONESYM_H
#define ONESYM_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>

static const char *TAB = "1023022312311302010120102300301020002120300301123102300010332020";
static uint8_t F[64];
static const int LAM[6][12] = {
    {3,9,2,9,2,9,2,9,2,9,3,9}, {9,0,9,2,9,0,9,2,9,0,9,3}, {2,9,3,9,3,9,3,9,3,9,2,9},
    {9,0,9,3,9,2,9,3,9,0,9,1}, {2,9,2,9,2,9,2,9,2,9,2,9}, {9,0,9,1,9,0,9,1,9,0,9,2}};
static inline int lam(long s, long d) { long a = d < 0 ? -d : d; return LAM[s % 6][a % 12]; }
static const int CORRMID[6] = {3, -1, 1, -1, 2, -1}, CORROUT[6] = {-1, 2, -1, 2, -1, 3}, WALLC[6] = {2, -1, 0, -1, 0, -1};
/* interior-like pattern at relative time r, column a>=1 (a%4==1: wall column) */
static inline int IL(long r, int a) {
    r = ((r % 6) + 6) % 6;
    if (a % 4 == 1) return WALLC[r];
    int k = (a - 2) / 4, pos = a - (4 * k + 2);
    return pos == 1 ? CORRMID[r] : CORROUT[r];
}

#define MAXLEV 40
static long Lb[MAXLEV]; static int wb[MAXLEV];
static int nlevels;
static int P0[3][3], P1[15][11], MID[15][11];
static int AX[14][12], AXY0[14] = {0,0,1,3,3,4,2,4,4,4,4,12,13,18}, CORNER[20][20]; static uint8_t INIT[56][113];
/* children of S_n: (lev, off), n>=1: S_{n-1}..S_0, S_0..S_{n-1} */
static int chn[MAXLEV]; static int chlev[MAXLEV][2 * MAXLEV]; static long choff[MAXLEV][2 * MAXLEV];
/* edges of S_n and its children: (off, halfwidth) */
static int edn[MAXLEV]; static long edoff[MAXLEV][4 * MAXLEV + 2]; static int edw[MAXLEV][4 * MAXLEV + 2];
static long mid_c[MAXLEV];   /* closing offset of the first-half S_1 inside S_n (n>=2) */

static void init_levels(void) {
    Lb[0] = 2;
    nlevels = 1;
    while (nlevels < MAXLEV && Lb[nlevels - 1] <= (LONG_MAX - 8) / 3) {
        Lb[nlevels] = 3 * Lb[nlevels - 1] + 8;
        nlevels++;
    }
    for (int n = 0; n < nlevels; n++) wb[n] = 4 * n + 1;
    for (int n = 0; n < nlevels; n++) {
        int c = 0; long r = 4;
        if (n >= 1) {
            for (int k = n - 1; k >= 0; k--) { chlev[n][c] = k; choff[n][c] = r; c++; r += Lb[k] + (k > 0 ? 4 : 2); }
            for (int k = 0; k <= n - 1; k++) { chlev[n][c] = k; choff[n][c] = r; c++; r += Lb[k] + 4; }
        }
        chn[n] = c;
        int e = 0; edoff[n][e] = 0; edw[n][e++] = wb[n]; edoff[n][e] = Lb[n]; edw[n][e++] = wb[n];
        for (int i = 0; i < c; i++) { edoff[n][e] = choff[n][i]; edw[n][e++] = wb[chlev[n][i]]; edoff[n][e] = choff[n][i] + Lb[chlev[n][i]]; edw[n][e++] = wb[chlev[n][i]]; }
        edn[n] = e;
        mid_c[n] = -1; for (int i = 0; i < c; i++) if (chlev[n][i] == 1) { mid_c[n] = choff[n][i] + Lb[1]; break; }
    }
}

/* interior pattern P_n at (r, d) */
static int P(int n, long r, int d) {
    int w = wb[n]; long L = Lb[n]; int ad = d < 0 ? -d : d;
    if (n == 0) return P0[r][d + 1];
    if (n == 1) return P1[r][d + 5];
    if (r == 0 || r == L || ad == w) return 0;
    for (int i = 0; i < chn[n]; i++) if (choff[n][i] <= r && r <= choff[n][i] + Lb[chlev[n][i]] && ad <= wb[chlev[n][i]]) return P(chlev[n][i], r - choff[n][i], d);
    long c = mid_c[n];
    if (c <= r && r <= c + 14 && ad <= 5) return MID[r - c][d + 5];
    for (int i = 0; i < n; i++) { int a = wb[chlev[n][i]]; if (r == choff[n][i] + 1 && ad == a + 1) return 2; if (r == choff[n][i] + 3 && ad == a + 1) return 1; }
    int best = -1, bestw = -1;
    for (int i = 0; i < edn[n]; i++) {
        long dr = r - edoff[n][i]; int a = edw[n][i], col = -1, wd = -100;
        if (dr == 1) { col = 1; wd = a - 1; } else if (dr == 2) { col = 1; wd = a - 2; } else if (dr == 3) { col = 0; wd = a - 3; }
        else if (dr == -1) { col = 0; wd = a - 3; } else if (dr == -2) { col = 1; wd = a - 4; } else if (dr == -3) { col = 1; wd = a - 5; }
        if (col >= 0 && ad <= wd && wd > bestw) { bestw = wd; best = col; }
    }
    if (best >= 0) return best;
    return IL(r, ad);      /* wall-continuation columns and corridors, phase relative to the box */
}

/* main boxes by the law */
typedef struct { int lev; long e; int blk; int big; } MBox;
static MBox *mb; static int nmb = 0;
static void init_mainboxes(long S) {
    mb = malloc(sizeof(MBox) * 1000000);
    mb[nmb++] = (MBox){1, 15, -1, 0}; mb[nmb++] = (MBox){1, 29, -1, 0};
    long e = 55; int j = 0;
    while (e < S + 100) {
        int big = 3 * j + 2;
        if (big + 3 >= nlevels) {
            fprintf(stderr, "requested time range exceeds the representable box levels\n");
            exit(EXIT_FAILURE);
        }
        mb[nmb++] = (MBox){big, e, j, 1};
        long p = e + Lb[big] + 4;
        for (int k = big - 1; k >= 0; k--) { mb[nmb++] = (MBox){k, p, j, 0}; p += Lb[k] + (k > 0 ? 4 : 2); }
        for (int k = 0; k <= big + 2; k++) { mb[nmb++] = (MBox){k, p, j, 0}; p += Lb[k] + (k < big + 2 ? 4 : 12); }
        e = p; j++;
    }
}
/* index of the last main box with e <= s (binary search) */
static inline int box_at(long s) {
    int lo = 0, hi = nmb - 1;
    while (lo < hi) { int m = (lo + hi + 1) / 2; if (mb[m].e <= s) lo = m; else hi = m - 1; }
    return lo;
}

static void init_tables(void) {
    for (int i = 0; i < 64; i++) F[i] = TAB[i] - '0';
    init_levels();
    long N = 200; uint8_t *pl = malloc(N * N);
    for (long i = 0; i < N; i++) { pl[i * N] = 1; pl[i] = 1; }
    for (long i = 1; i < N; i++) for (long j = 1; j < N; j++) pl[i * N + j] = F[(pl[(i - 1) * N + j] * 4 + pl[(i - 1) * N + j - 1]) * 4 + pl[i * N + j - 1]];
    #define PLL(s, d) pl[(((s) - (d)) / 2) * N + ((s) + (d)) / 2]
    for (int r = 0; r <= 2; r++) for (int d = -1; d <= 1; d++) if (((1 + r - d) & 1) == 0) P0[r][d + 1] = PLL(127 + r, d);
    for (int r = 0; r <= 14; r++) for (int d = -5; d <= 5; d++) if (((5 + r - d) & 1) == 0) P1[r][d + 5] = PLL(15 + r, d);
    for (int r = 0; r <= 14; r++) for (int d = -5; d <= 5; d++) if (((9 + 18 + r - d) & 1) == 0) MID[r][d + 5] = PLL(55 + 18 + r, d);
    for (int x = 0; x <= 13; x++) for (int k = 0; k < 12; k++) AX[x][k] = pl[x * N + 24 + k];
    for (int x = 0; x < 20; x++) for (int y = 0; y < 20; y++) CORNER[x][y] = pl[x * N + y];
    for (int s2 = 0; s2 < 56; s2++) for (int d = -s2; d <= s2; d++) if (((s2 - d) & 1) == 0) INIT[s2][d + 56] = PLL(s2, d);
    free(pl);
}

/* the description */
static int D(long s, long d) {
    long xx = (s - d) / 2, yy = (s + d) / 2; int ad = (int)(d < 0 ? -d : d);
    if (xx < 20 && yy < 20) return CORNER[xx][yy];
    if (xx <= 13 && yy >= AXY0[xx]) return AX[xx][((yy - 24) % 12 + 12) % 12];
    if (yy <= 13 && xx >= AXY0[yy]) return AX[yy][((xx - 24) % 12 + 12) % 12];
    if (s < 56) return INIT[s][d + 56];
    int qc = box_at(s);
    for (int q = qc - 1; q <= qc; q++) if (q >= 0 && mb[q].e <= s && s <= mb[q].e + Lb[mb[q].lev] && ad <= wb[mb[q].lev]) return P(mb[q].lev, s - mb[q].e, (int)d);
    int blk = mb[qc].blk;
    if (qc + 1 < nmb && mb[qc + 1].big && s >= mb[qc + 1].e - 1) blk = mb[qc + 1].blk;
    /* MID block for a main S_0 pair */
    for (int q = qc - 2; q <= qc + 1; q++) if (q >= 1 && q + 1 < nmb && mb[q].lev == 0 && mb[q + 1].lev == 0 && mb[q - 1].lev == 1) {
        long c0 = mb[q - 1].e + Lb[1];
        if (c0 <= s && s <= c0 + 14 && ad <= 5) return MID[s - c0][d + 5];
    }
    /* gap-12 junction before a big box */
    for (int q = qc; q <= qc + 1; q++) if (q >= 3 && q < nmb && mb[q].big) {
        long dr = s - mb[q].e; int a = wb[mb[q].lev], col = -1, wd = -100;
        if (dr == -10) { col = 1; wd = a - 4; } else if (dr == -9) { col = 0; wd = a - 3; } else if (dr == -8) { col = 0; wd = a; }
        else if (dr == -7) { col = 1; wd = a - 1; } else if (dr == -6) { col = 1; wd = a; } else if (dr == -5) { col = 0; wd = a + 1; }
        else if (dr == -4) { col = 0; wd = a + 4; } else if (dr == -3) { col = 1; wd = a + 3; } else if (dr == -2) { col = 1; wd = a + 2; } else if (dr == -1) { col = 0; wd = a + 1; }
        if (col >= 0 && ad <= wd) return col;
    }
    /* transition rows around main edges (widest wins) */
    int best = -1, bestw = -1;
    for (int q = qc - 2; q <= qc + 2; q++) {
        if (q < 0 || q >= nmb) continue;
        long ed[2] = {mb[q].e, mb[q].e + Lb[mb[q].lev]}; int a = wb[mb[q].lev];
        for (int u = 0; u < 2; u++) {
            long dr = s - ed[u]; int col = -1, wd = -100;
            if (dr == 1) { col = 1; wd = a - 1; } else if (dr == 2) { col = 1; wd = a - 2; } else if (dr == 3) { col = 0; wd = a - 3; }
            else if (dr == -1) { col = 0; wd = a - 3; } else if (dr == -2) { col = 1; wd = a - 4; } else if (dr == -3) { col = 1; wd = a - 5; }
            if (col >= 0 && ad <= wd && wd > bestw) { bestw = wd; best = col; }
        }
    }
    if (best >= 0) return best;
    /* extra cells after main openings with e = 1 mod 6 */
    for (int q = qc - 1; q <= qc + 1; q++) if (q >= 0 && q < nmb && mb[q].e % 6 == 1) {
        int a = wb[mb[q].lev]; long dr = s - mb[q].e;
        if (dr == 1 && ad == a + 1) return 2;
        if (dr == 3 && ad == a + 1) return 1;
    }
    /* virtual box of the block: interior-like pattern with absolute phase 3 */
    if (blk >= 0 && ad >= 1 && ad <= 12 * blk + 13) return IL(s + 3, ad);
    return lam(s, d);
}
#endif
