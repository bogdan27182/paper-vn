/* Exhaustive search for the required rectangle nesting over three letters.
 * Build: gcc -O3 -fopenmp code/search_three_letters.c -o /tmp/search_three_letters
 * Run: OMP_NUM_THREADS=8 /tmp/search_three_letters [first-key last-key]
 * Colour conjugacy permits boundary 0 and g(0,0,0)=1. The remaining
 * 17 symmetric table entries run through all 3^17 possibilities.
 */
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

enum { SMALL = 40, LARGE = 192, CAP_SMALL = 306, CAP_LARGE = 6820 };
typedef struct {
    int e, w, end, count, enclosed, rank;
    int descendants[8];
} Box;

static void generate(uint8_t *t, int n, const uint8_t f[27]) {
    for (int i = 0; i < n; ++i) t[i] = t[i*n] = 0;
    for (int i = 1; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            int a = t[(i-1)*n+j], b = t[(i-1)*n+j-1], c = t[i*n+j-1];
            t[i*n+j] = t[j*n+i] = f[9*a+3*b+c];
        }
    }
}

static int value(const uint8_t *t, int n, int s, int d) {
    return t[((s-d)/2)*n+(s+d)/2];
}

static int boxes(const uint8_t *t, int n, int colour, int width, Box *out, int cap) {
    int count = 0, last = 2*n-2;
    for (int e = 1; e < last; ++e) for (int w = 1; w <= width && w <= e; ++w) {
        if ((e-w)%2 || e+w+2 > last) continue;
        int good = 1;
        for (int d = -w; d <= w; d += 2)
            if (value(t,n,e,d) != colour) { good = 0; break; }
        if (!good) continue;
        for (int end = e+2; end+w <= last; end += 2) {
            int nonzero = 0;
            for (int d = -w+1; d < w; d += 2)
                if (value(t,n,end-1,d) != colour) { nonzero = 1; break; }
            if (!nonzero) break;
            if (end > e+2) {
                nonzero = 0;
                for (int d = -w+2; d < w; d += 2)
                    if (value(t,n,end-2,d) != colour) { nonzero = 1; break; }
                if (!nonzero) break;
            }
            if (value(t,n,end,w) != colour || value(t,n,end,-w) != colour) break;
            good = 1;
            for (int d = -w+2; d < w; d += 2)
                if (value(t,n,end,d) != colour) { good = 0; break; }
            if (good) {
                if (count == cap) {
                    fprintf(stderr,"Rectangle capacity bound violated.\n");
                    abort();
                }
                out[count++] = (Box){.e=e, .w=w, .end=end, .rank=-1};
                /* A second closing edge would make this entire row an
                 * all-border-colour interior row. */
                break;
            }
        }
    }
    return count;
}

static int inside(const Box *a, const Box *b) {
    return b->e < a->e && a->end < b->end && a->w < b->w;
}

static int nested_pair(const Box *b, int n) {
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
        if (inside(b+i,b+j)) return 1;
    return 0;
}

static void classify(Box *b, int n) {
    for (int i = 0; i < n; ++i) {
        b[i].count = b[i].enclosed = 0;
        b[i].rank = -1;
    }
    for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) {
        if (!inside(b+j,b+i)) continue;
        b[j].enclosed = 1;
        if (b[i].count < 8) b[i].descendants[b[i].count] = j;
        ++b[i].count;
    }
    for (int i = 0; i < n; ++i) if (!b[i].count) b[i].rank = 0;
    for (int i = 0; i < n; ++i) if (b[i].count == 2) {
        int a = b[i].descendants[0], c = b[i].descendants[1];
        if (!b[a].count && !b[c].count) b[i].rank = 1;
    }
    for (int i = 0; i < n; ++i) if (b[i].count >= 4 && b[i].count <= 8) {
        int children[8], number = 0;
        for (int a = 0; a < b[i].count; ++a) {
            int j = b[i].descendants[a], direct = 1;
            for (int c = 0; c < b[i].count; ++c)
                if (inside(b+j,b+b[i].descendants[c])) { direct = 0; break; }
            if (direct) children[number++] = j;
        }
        if (number != 4) continue;
        for (int left = 0; left < 4; ++left) for (int right = 0; right < 4; ++right) {
            if (left == right || b[children[left]].rank != 1 || b[children[right]].rank != 1) continue;
            int good = 1;
            for (int a = 0; a < 4; ++a) if (a != left && a != right)
                if (b[children[a]].rank != 0 || b[children[a]].e < b[children[left]].e ||
                    b[children[a]].e > b[children[right]].e) good = 0;
            if (good) b[i].rank = 2;
        }
    }
}

/* The presence of S_2 alone is necessary. Its possible enclosures outside
 * the finite window cannot change its collection of contained rectangles. */
static int candidate(const Box *b, int n) {
    for (int i = 0; i < n; ++i) if (b[i].rank == 2) return 1;
    return 0;
}

int main(int argc, char **argv) {
    const long total = 129140163L;
    long first = argc > 1 ? strtol(argv[1],NULL,10) : 0;
    long last = argc > 2 ? strtol(argv[2],NULL,10) : total;
    if (first < 0 || last < first || last > total) return 2;
    long tested = 0, density = 0, pairs = 0, hits = 0, level2 = 0;
    int maxsmall = 0, maxlarge = 0;
    double start = omp_get_wtime();
    printf("Normalized symmetric rules: %ld, keys [%ld,%ld), boundary 0, g(0,0,0)=1.\n",total,first,last);
    printf("Rectangle bounds: %d in the small window, %d in the large window.\n",CAP_SMALL,CAP_LARGE);
    fflush(stdout);
    #pragma omp parallel reduction(+:tested,density,pairs,hits,level2) reduction(max:maxsmall,maxlarge)
    {
        uint8_t f[27], small[SMALL*SMALL], large[LARGE*LARGE];
        Box *bs = malloc(CAP_SMALL*sizeof(*bs));
        Box *bl = malloc(CAP_LARGE*sizeof(*bl));
        if (!bs || !bl) abort();
        #pragma omp for schedule(dynamic,1024)
        for (long key = first; key < last; ++key) {
            long digits = key;
            for (int y = 0; y < 3; ++y) for (int x = 0; x < 3; ++x)
                for (int z = x; z < 3; ++z) {
                    int v = 1;
                    if (x || y || z) { v = digits%3; digits /= 3; }
                    f[9*x+3*y+z] = f[9*z+3*y+x] = v;
                }
            generate(small,SMALL,f);
            ++tested;
            int counts[3] = {0,0,0};
            for (int i = 0; i < SMALL*SMALL; ++i) ++counts[small[i]];
            if (counts[0] > 1552 || counts[1] > 1552 || counts[2] > 1552) continue;
            ++density;
            int generated = 0;
            for (int c = 0; c < 3; ++c) {
                int ns = boxes(small,SMALL,c,9,bs,CAP_SMALL);
                if (ns > maxsmall) maxsmall = ns;
                if (!nested_pair(bs,ns)) continue;
                ++pairs;
                if (!generated) { generate(large,LARGE,f); generated = 1; }
                int nl = boxes(large,LARGE,c,40,bl,CAP_LARGE);
                if (nl > maxlarge) maxlarge = nl;
                classify(bl,nl);
                for (int i = 0; i < nl; ++i) if (bl[i].rank == 2) { ++level2; break; }
                if (!candidate(bl,nl)) continue;
                ++hits;
                #pragma omp critical
                {
                    printf("CANDIDATE key=%ld colour=%d boxes=%d table=",key,c,nl);
                    for (int i = 0; i < 27; ++i) printf("%d",f[i]);
                    putchar('\n'); fflush(stdout);
                }
            }
        }
        free(bs); free(bl);
    }
    printf("Completed: %ld rules, %ld density survivors, %ld rule/colour pairs with early nesting.\n",tested,density,pairs);
    printf("Pairs with a level-2 box: %ld. Necessary-condition candidates: %ld.\n",level2,hits);
    printf("Maximum rectangle counts: small %d, large %d. No capacity failures.\n",maxsmall,maxlarge);
    printf("Elapsed seconds: %.3f.\n",omp_get_wtime()-start);
    return 0;
}
