/* Verify the virtual-box patterns Q_n against the four-colour description.
 * gcc -O2 -fwrapv code/check_virtual_boxes.c -o /tmp/check_virtual_boxes
 * /tmp/check_virtual_boxes
 */
#include "onesym.h"
#include <assert.h>

static int Q(int n, long r, int d) {
    if (n == 0) return 0;
    if (abs(d) == wb[n] && r > 0 && r < Lb[n]) return IL(r, wb[n]);
    return P(n, r, d);
}

int main(void) {
    init_tables();
    assert(P(0, 1, 0) == 1 && Q(0, 1, 0) == 0);
    long tested = 0;
    for (int n = 1; n <= 8; n++) {
        int k = n-1, w = wb[k];
        long ell = Lb[k];
        for (int part = 0; part < 3; part++) {
            long offset = 4 + part*ell;
            for (long r = 0; r <= ell; r++) for (int d = -w; d <= w; d++) {
                if ((w+r-d) & 1) continue;
                int expected = part == 1 ? Q(k,r,d) : P(k,r,d);
                assert(P(n,offset+r,d) == expected);
                assert(Q(n,offset+r,d) == expected);
                tested++;
            }
        }
        for (int i = 0; i < chn[n]; i++) {
            int level = chlev[n][i];
            long expected = (i < n) ? (Lb[n]+4)/2 - 3*(Lb[level]+4)/2 + 4
                                     : (Lb[n]+4)/2 + (Lb[level]+4)/2 - 4;
            assert(choff[n][i] == expected);
        }
    }
    printf("P_n and Q_n assembly, levels 1..8: %ld constituent cells, no mismatches.\n", tested);
    printf("Level-zero exception: ordinary centre 1, virtual centre 0; verified.\n");
    printf("Closed formulas for all child offsets, levels 1..8: verified.\n");
    init_mainboxes(630000);
    long cells = 0;
    for (int i = 0; i < nmb; i++) if (mb[i].big && mb[i].blk <= 2) {
        int n = mb[i].lev+1;
        long e = mb[i].e-4;
        for (long r = 0; r <= Lb[n]; r++) for (int d = -wb[n]; d <= wb[n]; d++) {
            if ((e+r-d) & 1) continue;
            assert(D(e+r,d) == Q(n,r,d));
            cells++;
        }
    }
    printf("Global virtual boxes V_0, V_1, V_2: %ld cells, no mismatches.\n", cells);
    free(mb);
    return 0;
}
