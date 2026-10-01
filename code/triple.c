/* triple.c -- (N1) the triple structure of the interior patterns P_m of the four-letter description:
 *   for m >= 2, with l = l_{m-1}, w = v_{m-1}, v = v_m:
 *   A1 = rows [4, 4+l], A2 = [4+l, 4+2l], A3 = [4+2l, 4+3l], all with |d| <= w:
 *     on A1 and A3: P_m = P_{m-1} translated;  on A2 with |d| < w: P_m = P_{m-1} translated;
 *     on the walls |d| = w of A2 strictly between its edges: P_m = I(r, w)  (dashed);
 *   frame: border white; rows 1,2,3: red |d|<=v-1, red |d|<=v-2, white |d|<=v-3; rows l_m-1, l_m-2, l_m-3: white |d|<=v-3,
 *     red |d|<=v-4, red |d|<=v-5; P_m(5, +-(w+1)) = 2, P_m(7, +-(w+1)) = 1; all other frame cells I(r, |d|).
 * (N3) ordinals: all white-bordered rectangles centred on the diagonal with a non-white cell in every interior row,
 *   in the band |d| <= WB for s < SO (description D), compared with the boxes of the family C.
 */
#include "onesym.h"
static int frame_value(int m, long r, int d) {
    long L = Lb[m]; int v = wb[m], w = wb[m - 1]; int ad = d < 0 ? -d : d;
    if (r == 0 || r == L || ad == v) return 0;
    if (r == 1 && ad <= v - 1) return 1; if (r == 2 && ad <= v - 2) return 1; if (r == 3 && ad <= v - 3) return 0;
    if (r == L - 1 && ad <= v - 3) return 0; if (r == L - 2 && ad <= v - 4) return 1; if (r == L - 3 && ad <= v - 5) return 1;
    if (r == 5 && ad == w + 1) return 2; if (r == 7 && ad == w + 1) return 1;
    return IL(r, ad);
}
int main(int argc, char **argv) {
    int MMAX = argc > 1 ? atoi(argv[1]) : 10; long SO = argc > 2 ? atol(argv[2]) : 1000000; int WB = argc > 3 ? atoi(argv[3]) : 100;
    init_tables();
    for (int m = 2; m <= MMAX; m++) {
        long L = Lb[m], l = Lb[m - 1]; int v = wb[m], w = wb[m - 1];
        long nA = 0, badA = 0, nW = 0, badW = 0, nF = 0, badF = 0;
        for (long r = 0; r <= L; r++) for (int d = -v; d <= v; d++) {
            if (((v + r - d) & 1) != 0) continue;
            int ad = d < 0 ? -d : d, p = P(m, r, d);
            int inA1 = r >= 4 && r <= 4 + l && ad <= w, inA2 = r >= 4 + l && r <= 4 + 2 * l && ad <= w, inA3 = r >= 4 + 2 * l && r <= 4 + 3 * l && ad <= w;
            if (inA1 || inA3) { long o = inA1 ? 4 : 4 + 2 * l; nA++; if (p != P(m - 1, r - o, d)) badA++; }
            else if (inA2 && ad < w) { nA++; if (p != P(m - 1, r - 4 - l, d)) badA++; }
            else if (inA2) { nW++; if (p != IL(r, w)) badW++; }
            else { nF++; if (p != frame_value(m, r, d)) { if (badF < 5) printf("  frame mismatch m=%d r=%ld d=%d P=%d desc=%d\n", m, r, d, p, frame_value(m, r, d)); badF++; } }
        }
        printf("m=%2d: sub-box cells %9ld (mismatch %ld), core-wall cells %7ld (mismatch %ld), frame cells %8ld (mismatch %ld)\n", m, nA, badA, nW, badW, nF, badF);
    }
    /* m = 2: the middle block is P_1 with dashed walls */
    long badM = 0; for (int r = 0; r <= 14; r++) for (int d = -5; d <= 5; d++) { if (((5 + r - d) & 1) != 0) continue; int ad = d < 0 ? -d : d;
        int expect = (ad == 5 && r > 0 && r < 14) ? IL(r, 5) : P1[r][d + 5]; if (MID[r][d + 5] != expect) badM++; }
    printf("middle block = P_1 with dashed walls: mismatches %ld\n", badM);
    /* ---------- (N3) ordinals ---------- */
    init_mainboxes(SO + 1000);
    /* boxes of C meeting s < SO: list (e, lev) */
    long nbox = 0, cap = 1 << 24; long *be = malloc(cap * sizeof(long)); int *bl = malloc(cap * sizeof(int));
    long stk_e[4096]; int stk_l[4096]; int sp = 0;
    for (int q = 0; q < nmb; q++) { if (mb[q].e >= SO) break; stk_e[sp] = mb[q].e; stk_l[sp] = mb[q].lev; sp++;
        while (sp) { sp--; long e = stk_e[sp]; int n = stk_l[sp]; if (e >= SO) continue; if (e + Lb[n] < SO && wb[n] <= WB) { be[nbox] = e; bl[nbox] = n; nbox++; }
            for (int i = 0; i < chn[n]; i++) { stk_e[sp] = e + choff[n][i]; stk_l[sp] = chlev[n][i]; sp++; } } }
    printf("boxes of C inside s < %ld, half-width <= %d: %ld\n", SO, WB, nbox);
    /* enumerate rectangles: for each odd half-width W and top row s0 == W (mod 2) */
    long found = 0, matched = 0, extra = 0;
    uint8_t *isbox = calloc(SO, 1);
    for (int W = 1; W <= WB; W++) {
        memset(isbox, 0, SO);
        for (long i = 0; i < nbox; i++) if (wb[bl[i]] == W) isbox[be[i]] = (uint8_t)(bl[i] + 1);
        for (long s0 = W + 2; s0 < SO; s0++) {
            if (((s0 - W) & 1) != 0) continue;
            int ok = 1; for (int d = -W; d <= W; d += 2) if (D(s0, d) != 0) { ok = 0; break; } if (!ok) continue;
            /* walls downwards */
            for (long s1 = s0 + 2; s1 < SO; s1 += 2) {
                if (D(s1, W) != 0 || D(s1, -W) != 0) break;
                /* interior row s1-1 must contain a non-white cell with |d| < W */
                int nw = 0; for (int d = -(W - 1); d <= W - 1; d += 2) if (D(s1 - 1, d) != 0) { nw = 1; break; } if (!nw) break;
                if (s1 - 2 > s0) { int nw2 = 0; for (int d = -W + 2; d <= W - 2; d += 2) if (D(s1 - 2, d) != 0) { nw2 = 1; break; } if (!nw2) break; }
                int full = 1; for (int d = -W; d <= W; d += 2) if (D(s1, d) != 0) { full = 0; break; }
                if (full) { found++; if (isbox[s0] && s1 == s0 + Lb[isbox[s0] - 1]) matched++; else { if (extra < 10) printf("  EXTRA rectangle W=%d s=%ld..%ld\n", W, s0, s1); extra++; } }
            }
        }
    }
    printf("rectangles found %ld, equal to boxes of C: %ld, others: %ld (boxes expected %ld)\n", found, matched, extra, nbox);
    return 0;
}
