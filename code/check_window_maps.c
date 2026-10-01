/* Check the explicit window maps against the cell rules.
 * The proof of their validity for arbitrary levels is in the manuscript.
 * gcc -O2 -fsanitize=undefined -fno-sanitize-recover=undefined
 *     code/check_window_maps.c -o /tmp/check_window_maps
 */
#include "onesym.h"
#include <assert.h>

static long compress(long u, long old, long next) {
    assert(0 <= u && u < old && next >= 21 && (old-next)%6 == 0);
    if (u <= 9) return u;
    if (old-u <= 5) return next-(old-u);
    return 10+(u-10)%6;
}

static int frame(int m, long r, int d) {
    if (r < 2 || r > Lb[m] || d+1 > wb[m]) return 0;
    long length = Lb[m-1];
    for (int i = 0; i < 3; ++i) {
        long a = 4+i*length;
        if (r-2 >= a && r <= a+length &&
            (i == 1 ? d+1 < wb[m-1] : d+1 <= wb[m-1])) return 0;
    }
    return 1;
}

static long frame_time(int m, long r) {
    long length = Lb[m-1], next = Lb[m-2];
    if (r < 4) return r;
    if (r >= 4+3*length) return Lb[m-1]-(Lb[m]-r);
    int part = (int)((r-4)/length);
    return 4+part*next+compress(r-4-part*length,length,next);
}

static long frame_count = 0, block_count = 0;
static void check_frame(int m, long r) {
    const int dr[4] = {0,-1,-2,-1}, dd[4] = {0,1,0,-1};
    for (int d = 0; d < wb[m]; ++d) {
        if ((wb[m]+r-d)%2 || !frame(m,r,d)) continue;
        long t = frame_time(m,r);
        int a = d <= 7 ? d : d-4;
        assert(frame(m-1,t,a));
        for (int i = 0; i < 4; ++i) {
            int before = P(m,r+dr[i],d+dd[i]);
            int after = P(m-1,t+dr[i],a+dd[i]);
            if (before != after) {
                fprintf(stderr,"Frame mismatch m=%d (%ld,%d) -> (%ld,%d), component %d: %d != %d\n",
                        m,r,d,t,a,i,before,after);
                exit(1);
            }
        }
        ++frame_count;
    }
}

static long block_time(long r, long length, long next) {
    if (r < length) return compress(r,length,next);
    if (r < 2*length) return next+compress(r-length,length,next);
    if (r < 2*length+4) return 2*next+(r-2*length);
    if (r < 5*length+12)
        return 2*next+4+compress(r-2*length-4,3*length+8,3*next+8);
    return 5*next+12+(r-5*length-12);
}

static int block_boundary(int j, long r, int d) {
    long length = Lb[3*j+3];
    int width = wb[3*j+3];
    if (r >= 2 && r <= length && d+1 < width) return 0;
    if (r-2 >= length && r <= 2*length && d+1 <= width) return 0;
    if (r-2 >= 2*length+4 && r <= 5*length+12 && d+1 <= width+4) return 0;
    return 1;
}

static void check_block(int j, long start, long previous, long r) {
    const int dr[4] = {0,-1,-2,-1}, dd[4] = {0,1,0,-1};
    long length = Lb[3*j+3], next = Lb[3*j];
    long t = previous+block_time(r,length,next), s = start+r;
    for (int d = 0; d <= 12*j+27; ++d) {
        if ((s-d)%2 || !block_boundary(j,r,d)) continue;
        int a = d <= 13 ? d : d-12;
        for (int i = 0; i < 4; ++i) {
            int before = D(s+dr[i],d+dd[i]);
            int after = D(t+dr[i],a+dd[i]);
            if (before != after) {
                fprintf(stderr,"Block mismatch j=%d (%ld,%d) -> (%ld,%d), component %d: %d != %d\n",
                        j,s,d,t,a,i,before,after);
                exit(1);
            }
        }
        ++block_count;
    }
}

int main(void) {
    init_tables();
    for (int m = 6; m <= 8; ++m)
        for (long r = 2; r <= Lb[m]; ++r) check_frame(m,r);
    printf("All frame windows at levels 6..8: %ld ordered windows, no mismatches.\n",frame_count);
    for (int m = 9; m <= 30; ++m) {
        for (long r = 2; r < 4; ++r) check_frame(m,r);
        long length = Lb[m-1];
        for (int part = 0; part < 3; ++part) {
            for (long u = 0; u <= 15; ++u) check_frame(m,4+part*length+u);
            for (long u = length-5; u < length; ++u) check_frame(m,4+part*length+u);
        }
        for (long r = Lb[m]-4; r <= Lb[m]; ++r) check_frame(m,r);
    }
    printf("Including representatives at levels 9..30: %ld ordered frame windows, no mismatches.\n",frame_count);
    long starts[11]; starts[0] = 51;
    for (int j = 0; j < 10; ++j) starts[j+1] = starts[j]+5*Lb[3*j+3]+20;
    init_mainboxes(starts[10]);
    for (int j = 3; j <= 9; ++j) {
        long length = Lb[3*j+3];
        long offsets[3] = {0,length,2*length+4};
        long lengths[3] = {length,length,3*length+8};
        for (int p = 0; p < 3; ++p) {
            for (long u = 0; u <= 15; ++u) check_block(j,starts[j],starts[j-1],offsets[p]+u);
            for (long u = lengths[p]-5; u < lengths[p]; ++u)
                check_block(j,starts[j],starts[j-1],offsets[p]+u);
        }
        for (long r = 2*length; r < 2*length+4; ++r) check_block(j,starts[j],starts[j-1],r);
        for (long r = 5*length+12; r < 5*length+20; ++r) check_block(j,starts[j],starts[j-1],r);
    }
    printf("Block-boundary representatives in blocks 3..9: %ld ordered windows, no mismatches.\n",block_count);
    free(mb);
    return 0;
}
