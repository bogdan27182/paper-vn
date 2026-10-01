/* Verify the finite region R and the axis separators in the four-letter proof.
 * gcc -O2 -fopenmp -fsanitize=undefined -fno-sanitize-recover=undefined
 *     code/check_four_letters.c -o /tmp/check_four_letters
 * OMP_NUM_THREADS=8 /tmp/check_four_letters
 */
#include "onesym.h"
#include <assert.h>
#include <omp.h>

static void record_edges(int n, long e, int *widths, long limit) {
    if (e >= limit) return;
    if (widths[e] < wb[n]) widths[e] = wb[n];
    long end = e + Lb[n];
    if (end < limit && widths[end] < wb[n]) widths[end] = wb[n];
    for (int i = 0; i < chn[n]; i++)
        record_edges(chlev[n][i], e + choff[n][i], widths, limit);
}

int main(void) {
    const long limit = 613224;
    init_tables();
    init_mainboxes(limit + 1000);
    printf("Initialized levels 0..%d without signed overflow.\n", nlevels-1);
    assert(D(857, 19) == 0);

    long failures = 0, count = 0;
    uint64_t observed = 0;
    #pragma omp parallel for reduction(+:failures,count) reduction(|:observed) schedule(dynamic,1024)
    for (long s = 2; s < limit; s++) {
        /* Check each interior cell of R once: the band, then the two margins. */
        for (int region = 0; region < 3; region++) {
            long lo = region == 0 ? -60 : (region == 1 ? s-30 : 2-s);
            long hi = region == 0 ? 60 : (region == 1 ? s-2 : 30-s);
            for (long d = lo; d <= hi; d++) {
                if (d < 2-s || d > s-2 || ((s-d)&1)) continue;
                if (region && d >= -60 && d <= 60) continue;
                int a = D(s-1,d+1), b = D(s-2,d), c = D(s-1,d-1);
                int index = 16*a + 4*b + c;
                assert(index >= 0 && index < 64);
                observed |= UINT64_C(1) << index;
                if (D(s,d) != F[index]) failures++;
                count++;
            }
        }
    }
    printf("Region R: %ld recurrence windows, %ld failures.\n", count, failures);
    int present = 0;
    printf("Absent argument triples:");
    for (int i = 0; i < 64; i++) {
        if (observed & (UINT64_C(1) << i)) present++;
        else printf(" %d%d%d", i/16, (i/4)%4, i%4);
    }
    printf("\nPresent triples: %d.\n", present);
    fflush(stdout);
    assert(failures == 0 && present == 56);

    int *widths = malloc((size_t)limit * sizeof(*widths));
    assert(widths);
    for (long s = 0; s < limit; s++) widths[s] = -1;
    for (int q = 0; q < nmb; q++) {
        record_edges(mb[q].lev, mb[q].e, widths, limit);
        if (mb[q].big && mb[q].e < limit) {
            long e = mb[q].e;
            int j = mb[q].blk;
            widths[e-4] = 12*j+13;
            widths[e-8] = 12*j+9;
        }
    }
    long edge_count = 0;
    for (long s = 58; s < limit; s++) {
        int a = D(s,s-30), b = D(s,s-32);
        assert(a == lam(s,s-30) && b == lam(s,s-32));
        assert(a != 0 || b != 0);
        if (widths[s] >= 1) {
            assert(D(s,widths[s]+2) != 0);
            edge_count++;
        }
    }
    printf("Axis separators verified for 58 <= s < %ld; %ld edge endpoints have non-white neighbours.\n",
           limit, edge_count);
    free(widths);
    free(mb);
    return 0;
}
