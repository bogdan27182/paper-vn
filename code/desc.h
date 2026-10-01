/* Explicit description D of the double sequence
 *     t(i,0)=t(0,j)=1,  t(i,j)=f(t(i-1,j),t(i-1,j-1),t(i,j-1)),
 *     f(x,y,z)=x^2y^2+y^2z^2+y^2+xz+2y+2  over F_3,
 * written WITHOUT using the recurrence.  Rotated coordinates:
 *     s = i+j (time),  d = j-i (position),  cells with s = d (mod 2), |d| <= s.
 * Colours: 0 white, 1 red, 2 blue.
 */
#ifndef DESC_H
#define DESC_H
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static inline long Lbox(int n) { return (1L << (n + 2)) - 6; }        /* time length of S_n   */
static inline long Wbox(int n) { return 2L * n - 1; }                 /* half-width of S_n    */
static inline long Emain(int n) { return (1L << (n + 2)) - 2L * n + 3; } /* opening time of S_n */

static inline long md4(long a) { return ((a % 4) + 4) % 4; }
static inline long labs_(long a) { return a < 0 ? -a : a; }

/* odd lattice with phase eps (eps=2: global environment, eps=0: reversed) */
static inline uint8_t Lam(int eps, long s, long d) {
    return ((s & 1) && md4(s - labs_(d) - eps) == 0) ? 1 : 2;
}
/* central lattice */
static inline uint8_t Gam(long s, long d) {
    return (!(s & 1) && md4(s - d) == 0) ? 1 : 2;
}

static inline uint8_t fF3(int x, int y, int z) {
    return (uint8_t)((x * x * y * y + y * y * z * z + y * y + x * z + 2 * y + 2) % 3);
}

typedef struct {
    long S;      /* rows s = 0..S-1 */
    long W;      /* band |d| <= W   */
    uint8_t *a;  /* a[s*(2W+1) + d + W] */
} Band;

static inline uint8_t *bref(Band *B, long s, long d) { return &B->a[s * (2 * B->W + 1) + d + B->W]; }

static inline void bput(Band *B, long s, long d, uint8_t v) {
    if (s >= 0 && s < B->S && labs_(d) <= B->W && labs_(d) <= s && ((s - d) & 1) == 0)
        *bref(B, s, d) = v;
}

/* value outside the band (|d| > W): boundary lines or global environment */
static inline uint8_t Dfar(long s, long d) {
    long ad = labs_(d);
    if (ad == s) return 1;
    if (ad == s - 2 && s >= 3) return 0;
    return Lam(2, s, d);
}

static inline uint8_t Dget(Band *B, long s, long d) {
    if (labs_(d) <= B->W) return *bref(B, s, d);
    return Dfar(s, d);
}

/* children of a box of level n opening at time e: first half S_{n-2},...,S_1, then S_1,...,S_{n-2} */
static int children(int n, long e, int *lev, long *tim) {
    int c = 0;
    if (n < 3) return 0;
    long a = e + 4;
    for (int k = n - 2; k >= 1; k--) { lev[c] = k; tim[c] = a; c++; a += Lbox(k) + 4; }
    long b = e + Lbox(n) - 4;          /* closing time of the last child */
    int c0 = c;
    for (int k = n - 2; k >= 1; k--) { lev[c] = k; tim[c] = b - Lbox(k); c++; b -= Lbox(k) + 4; }
    /* reverse second half so that it is in time order */
    for (int i = c0, j = c - 1; i < j; i++, j--) {
        int tl = lev[i]; lev[i] = lev[j]; lev[j] = tl;
        long tt = tim[i]; tim[i] = tim[j]; tim[j] = tt;
    }
    return c;
}

typedef int (*InsideFn)(long s, long d, const long *ctx);

/* half-gap rows of an edge (e,a): e+-1: Gamma on |d|<=a-1 ; e+-2: red on |d|<=a-2 */
static void halfgap(Band *B, long e, long a, InsideFn inside, const long *ctx) {
    long rows[4] = {e + 1, e - 1, e + 2, e - 2};
    long hs[4] = {a - 1, a - 1, a - 2, a - 2};
    for (int q = 0; q < 4; q++) {
        long s = rows[q], h = hs[q];
        if (h < 0 || s < 0 || s >= B->S) continue;
        for (long d = -h; d <= h; d++) {
            if (((s - d) & 1) != 0) continue;
            if (!inside(s, d, ctx)) continue;
            bput(B, s, d, q < 2 ? Gam(s, d) : 1);
        }
    }
}

static int inside_box(long s, long d, const long *ctx) {
    /* ctx = {e, w, L} : open interior */
    return s > ctx[0] && s < ctx[0] + ctx[2] && labs_(d) < ctx[1];
}

static void paint_box(Band *B, int n, long e) {
    if (e >= B->S) return;
    long w = Wbox(n), L = Lbox(n), sig = e + L / 2;
    int lev[128]; long tim[128];
    int nc = children(n, e, lev, tim);
    /* 1. environment: Lambda_0 before the midpoint, Lambda_2 after */
    for (long s = e + 1; s < e + L && s < B->S; s++)
        for (long d = -(w - 1); d <= w - 1; d++)
            if (((s - d) & 1) == 0) bput(B, s, d, Lam(s < sig ? 0 : 2, s, d));
    /* 2. switching diamond, radius w-3, centred at (sig,0) */
    if (n >= 2) {
        long r = w - 3;
        for (long s = sig - r; s <= sig + r; s++) {
            long h = r - labs_(s - sig);
            for (long d = -h; d <= h; d++)
                if (((s - d) & 1) == 0) bput(B, s, d, Gam(s, d));
        }
    }
    /* 3. transition rows of own edges and children's edges */
    long ctx[3] = {e, w, L};
    halfgap(B, e, w, inside_box, ctx);
    halfgap(B, e + L, w, inside_box, ctx);
    for (int c = 0; c < nc; c++) {
        halfgap(B, tim[c], Wbox(lev[c]), inside_box, ctx);
        halfgap(B, tim[c] + Lbox(lev[c]), Wbox(lev[c]), inside_box, ctx);
    }
    /* 4. children (overwrite their closed regions) */
    for (int c = 0; c < nc; c++) paint_box(B, lev[c], tim[c]);
    /* 5. border: walls and edges are white */
    for (long s = e; s <= e + L && s < B->S; s++)
        if (((s - w) & 1) == 0) { bput(B, s, w, 0); bput(B, s, -w, 0); }
    for (int q = 0; q < 2; q++) {
        long s = q ? e + L : e;
        for (long d = -w; d <= w; d++) if (((s - d) & 1) == 0) bput(B, s, d, 0);
    }
}

static int NMAIN;
static int outside_mains(long s, long d, const long *ctx) {
    (void)ctx;
    for (int n = 1; n <= NMAIN; n++)
        if (s >= Emain(n) && s <= Emain(n) + Lbox(n) && labs_(d) <= Wbox(n)) return 0;
    return 1;
}

/* paint rows 0..S-1 of the band |d|<=W */
static Band paint(long S, long W) {
    Band B; B.S = S; B.W = W;
    B.a = (uint8_t *)malloc((size_t)S * (2 * W + 1));
    memset(B.a, 9, (size_t)S * (2 * W + 1));
    for (long s = 0; s < S; s++)
        for (long d = -W; d <= W; d++)
            if (labs_(d) <= s && ((s - d) & 1) == 0) *bref(&B, s, d) = Lam(2, s, d);
    NMAIN = 0;
    while (Emain(NMAIN + 1) < S) NMAIN++;
    /* global transition rows: main edges and the virtual edge (3,1) formed by t(1,2), t(2,1) */
    halfgap(&B, 3, 1, outside_mains, NULL);
    for (int n = 1; n <= NMAIN; n++) {
        halfgap(&B, Emain(n), Wbox(n), outside_mains, NULL);
        halfgap(&B, Emain(n) + Lbox(n), Wbox(n), outside_mains, NULL);
    }
    for (int n = 1; n <= NMAIN; n++) paint_box(&B, n, Emain(n));
    /* boundary lines */
    for (long s = 0; s < S; s++) {
        bput(&B, s, s, 1); bput(&B, s, -s, 1);
        if (s >= 3) { bput(&B, s, s - 2, 0); bput(&B, s, -(s - 2), 0); }
    }
    return B;
}
#endif
