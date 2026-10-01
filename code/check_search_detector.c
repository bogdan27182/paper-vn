/* Check the search detector against direct enumeration and a positive control.
 * gcc -O2 -fopenmp -fsanitize=undefined -fno-sanitize-recover=undefined
 *     code/check_search_detector.c -o /tmp/check_search_detector
 */
#define main search_three_letters_main
#include "search_three_letters.c"
#undef main
#include "onesym.h"

static uint32_t seed = 1234567;
static unsigned random_value(void) {
    seed = 1664525U*seed+1013904223U;
    return seed;
}

static int brute_box(const uint8_t *t, int n, int e, int w, int end, int colour) {
    for (int s = e; s <= end; ++s) {
        int nonzero = 0;
        for (int d = -w; d <= w; ++d) {
            if ((s-d)%2) continue;
            int border = s == e || s == end || d == -w || d == w;
            int v = value(t,n,s,d);
            if (border && v != colour) return 0;
            if (!border && v != colour) nonzero = 1;
        }
        if (s != e && s != end && !nonzero) return 0;
    }
    return 1;
}

int main(void) {
    uint8_t a[SMALL*SMALL], b[SMALL*SMALL], f[27];
    for (int trial = 0; trial < 100; ++trial) {
        for (int x = 0; x < 3; ++x) for (int y = 0; y < 3; ++y)
            for (int z = x; z < 3; ++z) f[9*x+3*y+z] = f[9*z+3*y+x] = random_value()%3;
        memset(b,0,sizeof(b));
        for (int x = 1; x < SMALL; ++x) for (int y = 1; y < SMALL; ++y)
            b[x*SMALL+y] = f[9*b[(x-1)*SMALL+y]+3*b[(x-1)*SMALL+y-1]+b[x*SMALL+y-1]];
        generate(a,SMALL,f);
        assert(memcmp(a,b,sizeof(a)) == 0);
        const int n = 10, width = 4;
        uint8_t small[n*n];
        for (int x = 0; x < n; ++x) for (int y = 0; y < n; ++y)
            small[x*n+y] = trial%2 ? b[x*SMALL+y] : random_value()%3;
        for (int c = 0; c < 3; ++c) {
            Box found[100]; int number = boxes(small,n,c,width,found,100), brute = 0;
            for (int e = 1; e < 2*n-2; ++e) for (int w = 1; w <= width && w <= e; ++w) {
                if ((e-w)%2) continue;
                for (int end = e+2; end+w <= 2*n-2; end += 2) {
                    if (!brute_box(small,n,e,w,end,c)) continue;
                    ++brute;
                    int seen = 0;
                    for (int i = 0; i < number; ++i)
                        if (found[i].e == e && found[i].w == w && found[i].end == end) seen = 1;
                    assert(seen);
                }
            }
            assert(brute == number);
        }
    }
    puts("100 symmetric rules: triangular and full-square generation agree.");
    puts("300 small-grid rectangle enumerations agree with the direct definition.");
    init_tables();
    uint8_t *four = malloc(LARGE*LARGE);
    Box *found = malloc(CAP_LARGE*sizeof(*found));
    assert(four && found);
    for (int x = 0; x < LARGE; ++x) for (int y = 0; y < LARGE; ++y) {
        if (!x || !y) four[x*LARGE+y] = 1;
        else four[x*LARGE+y] = F[16*four[(x-1)*LARGE+y]+4*four[(x-1)*LARGE+y-1]+four[x*LARGE+y-1]];
    }
    int number = boxes(four,LARGE,0,40,found,CAP_LARGE);
    classify(found,number);
    assert(candidate(found,number));
    int confirmed = 0;
    for (int i = 0; i < number; ++i)
        if (found[i].e == 55 && found[i].w == 9 && found[i].end == 105) {
            assert(found[i].count == 8 && found[i].rank == 2);
            ++confirmed;
        }
    assert(confirmed == 1);
    puts("Positive control: the four-letter S_2 at rows 55..105 is accepted with eight descendants.");
    free(four); free(found);
    return 0;
}
