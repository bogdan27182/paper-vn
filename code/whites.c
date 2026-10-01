/* whites.c -- facts on white cells of the four-letter description used in the proof that the ordinals are the boxes:
 * (i)  on an even column |d| = W (2 <= W <= 200), no two cells [s,W], [s+2,W] with s >= 56 are both white;
 * (ii) on an odd row s >= 57, a white cell with |d| = 1 or |d| = 3 lies on an edge row, i.e. s is the opening or
 *      closing time of a box of C of half-width >= |d|, the opening time E_j - 4 of a virtual box V_j, or a bar E_j - 8;
 *      and such a row is white on all |d'| <= its half-width.
 */
#include "onesym.h"
int main(int argc, char **argv) {
    long S = argc > 1 ? atol(argv[1]) : 2000000;
    init_tables(); init_mainboxes(S + 1000);
    int *edgew = calloc(S + 20, sizeof(int)); for (long s = 0; s < S + 20; s++) edgew[s] = -1;
    long stk_e[4096]; int stk_l[4096];
    for (int q = 0; q < nmb; q++) { if (mb[q].e >= S) break; int sp = 0; stk_e[sp] = mb[q].e; stk_l[sp] = mb[q].lev; sp++;
        while (sp) { sp--; long e = stk_e[sp]; int n = stk_l[sp]; if (e >= S) continue;
            if (wb[n] > edgew[e]) edgew[e] = wb[n]; if (e + Lb[n] < S && wb[n] > edgew[e + Lb[n]]) edgew[e + Lb[n]] = wb[n];
            for (int i = 0; i < chn[n]; i++) { stk_e[sp] = e + choff[n][i]; stk_l[sp] = chlev[n][i]; sp++; } } }
    for (int q = 0; q < nmb; q++) if (mb[q].big && mb[q].e < S) { int j = mb[q].blk; long E = mb[q].e;
        if (12 * j + 13 > edgew[E - 4]) edgew[E - 4] = 12 * j + 13; if (12 * j + 9 > edgew[E - 8]) edgew[E - 8] = 12 * j + 9; }
    long bad1 = 0, bad2 = 0, bad3 = 0;
    for (int W = 2; W <= 200; W += 2) for (long s = 56; s + 2 < S; s++) { if (((s - W) & 1) != 0 || W > s - 34) continue;   /* x, y >= 16 for the cells and the cell two rows below */
        if (D(s, W) == 0 && D(s + 2, W) == 0) { if (bad1 < 5) printf("even column pair W=%d s=%ld\n", W, s); bad1++; } }
    for (long s = 57; s < S; s += 2) for (int d = 1; d <= 3; d += 2) {
        if (D(s, d) != 0) continue;
        if (edgew[s] < d) { if (bad2 < 5) printf("white at odd row s=%ld d=%d without edge (edgew=%d)\n", s, d, edgew[s]); bad2++; continue; }
        int a = edgew[s]; for (int dd = -a; dd <= a; dd += 2) if (D(s, dd) != 0) { if (bad3 < 5) printf("edge row s=%ld not white at d=%d\n", s, dd); bad3++; break; }
    }
    printf("s < %ld: even-column white pairs %ld, odd-row whites at |d|=1,3 off edge rows %ld, edge rows not white %ld\n", S, bad1, bad2, bad3);
    return 0;
}
