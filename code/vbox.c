/* vbox.c -- checks the virtual-box reading of the diagonal of tau (four letters):
 * (1) in block j the main boxes S_{3j+2}..S_0 S_0..S_{3j+2} open exactly at E_j-4 + (child offsets of a box of level 3j+3),
 *     and S_{3j+3} opens at E_j-4+l_{3j+3}, where the virtual box V_j would close;
 * (2) on the closed region of V_j (E_j-4 <= s <= E_j-4+l_{3j+3}, |d| <= 12j+13) the description D equals P_{3j+3}(s-E_j+4, d),
 *     except on the two walls |d| = 12j+13, where D follows the wall-column pattern I(s+3, 12j+13) = 2,0,0;
 * (3) the 12-row junction before the big box of block j+1 consists of transition rows of two white bars:
 *     half-width v_{3j+5} at E_{j+1}-8 and half-width v_{3j+6} (the opening edge of V_{j+1}) at E_{j+1}-4.
 *   ./vbox J      checks blocks 0..J-1
 */
#include "onesym.h"
int main(int argc, char **argv) {
    int J = argc > 1 ? atoi(argv[1]) : 3;
    init_tables(); init_mainboxes(4000000000L);
    long Ej[8]; int nj = 0;
    for (int q = 0; q < nmb; q++) if (mb[q].big) { if (nj < 8) Ej[nj++] = mb[q].e; }
    for (int j = 0; j < J && j + 1 < nj; j++) {
        int n = 3 * j + 3; long ev = Ej[j] - 4; int V = wb[n];
        /* (1) children positions */
        int q0 = -1; for (int q = 0; q < nmb; q++) if (mb[q].big && mb[q].blk == j) q0 = q;
        int bad1 = 0;
        for (int i = 0; i < chn[n]; i++) { MBox b = mb[q0 + i]; if (b.lev != chlev[n][i] || b.e != ev + choff[n][i]) bad1++; }
        MBox nxt = mb[q0 + chn[n]];
        printf("block %d: V opens at %ld, half-width %d; children positions %s; next main box level %d at %ld (V would close at %ld)\n",
               j, ev, V, bad1 ? "MISMATCH" : "match", nxt.lev, nxt.e, ev + Lb[n]);
        /* (2) V_j region vs P_{3j+3} */
        long diff_wall = 0, diff_other = 0, cells = 0; int wallvals[4] = {0, 0, 0, 0};
        for (long s = ev; s <= ev + Lb[n]; s++) for (int d = -V; d <= V; d++) {
            if (((s - d) & 1) != 0) continue;
            int x = D(s, d), p = P(n, s - ev, d); cells++;
            if (d == V || d == -V) { if (x != p) diff_wall++; if (s > ev && s < ev + Lb[n]) wallvals[x]++; }
            else if (x != p) { if (diff_other < 5) printf("  DIFF s=%ld d=%d D=%d P=%d\n", s, d, x, p); diff_other++; }
        }
        printf("  V_%d: %ld cells, differences off the walls: %ld, on the walls: %ld; wall values 0:%d 1:%d 2:%d 3:%d\n",
               j, cells, diff_other, diff_wall, wallvals[0], wallvals[1], wallvals[2], wallvals[3]);
        /* wall values follow I(s+3, V)? */
        long badI = 0; for (long s = ev + 1; s < ev + Lb[n]; s++) { long d = V; if (((s - d) & 1) != 0) continue; if (D(s, d) != IL(s + 3, V)) badI++; }
        printf("  wall cells not equal to I(s+3,V): %ld\n", badI);
        /* (3) junction before block j+1: rows E-12..E, which rows are white on which |d| */
        long E = Ej[j + 1];
        printf("  junction before block %d (E=%ld): white extent per row (max |d| with all cells |d'|<=|d| white, 0 = none):\n   ", j + 1, E);
        for (long s = E - 12; s <= E; s++) {
            int w = -1; for (int d = (s & 1); d < 200; d += 2) { if (D(s, d) != 0) break; w = d; }
            printf(" %ld:%d", s - E, w);
        }
        printf("\n   v_{3j+4}=%d v_{3j+5}=%d v_{3j+6}=%d\n", wb[3 * j + 4], wb[3 * j + 5], wb[3 * j + 6]);
    }
    /* the prelude: s < 56 */
    printf("prelude: white extent per odd row s=13..55:");
    for (long s = 13; s <= 55; s += 2) { int w = -1; for (int d = 1; d < 60 && d <= s; d += 2) { if (D(s, d) != 0) break; w = d; } if (w >= 1) printf(" %ld:%d", s, w); }
    printf("\n");
    return 0;
}
