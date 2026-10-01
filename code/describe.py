"""Independent explicit description D of the configuration (no recurrence used).
Rotated coordinates: s = i+j, d = j-i, cells with s = d (mod 2)."""

def L(n):  # duration (time extent) of S_n
    return 2**(n+2) - 6

def w(n):  # half-width of S_n
    return 2*n - 1

def e_main(n):  # opening time of main ordinal S_n
    return 2**(n+2) - 2*n + 3

def children(n, e):
    """direct subordinates of an ordinal of level n opening at time e: list of (level, opening time)"""
    if n < 3:
        return []
    res = []
    a = e + 4
    for k in range(n-2, 0, -1):          # first half: S_{n-2}, ..., S_1
        res.append((k, a))
        a += L(k) + 4
    end = e + L(n)
    b = end - 4
    second = []
    for k in range(n-2, 0, -1):          # second half, built from the end: S_{n-2} closes at end-4
        second.append((k, b - L(k)))
        b -= L(k) + 4
    return res + second[::-1]

def Lam(eps, s, d):     # odd lattice with phase eps (global environment: eps=2)
    return 1 if (s % 2 == 1 and (s - abs(d) - eps) % 4 == 0) else 2

def Gam(s, d):          # central lattice
    return 1 if (s % 2 == 0 and (s - d) % 4 == 0) else 2

def paint(N):
    """returns dict-like 2D array D[s][d+N] for s < N"""
    W = N
    D = [[None]*(2*W+1) for _ in range(N)]
    def put(s, d, v):
        if 0 <= s < N and abs(d) <= s:
            D[s][d+W] = v
    def cells_row(s, dmax):
        # cells on row s with |d| <= dmax and right parity
        start = -dmax
        if (s - start) % 2: start += 1
        return range(start, dmax+1, 2)
    # global environment
    for s in range(N):
        for d in cells_row(s, s):
            put(s, d, Lam(2, s, d))
    def halfgap(edges, inside):
        # edges: list of (time, halfwidth); paint half-gap rows restricted by predicate inside(s,d)
        for (e, a) in edges:
            for s, h, kind in ((e+1, a-1, 'G'), (e+2, a-2, 'R'), (e-1, a-1, 'G'), (e-2, a-2, 'R')):
                if h < 0: continue
                for d in cells_row(s, h):
                    if inside(s, d):
                        put(s, d, Gam(s, d) if kind == 'G' else 1)
    def ordinal(n, e):
        if e >= N: return
        ww, LL = w(n), L(n)
        kids = children(n, e)
        inside = lambda s, d: e < s < e + LL and abs(d) < ww
        # environment
        sig = e + LL//2
        for s in range(e+1, min(N, e+LL)):
            for d in cells_row(s, ww-1):
                put(s, d, Lam(0 if s < sig else 2, s, d))
        # diamond (n >= 2)
        if n >= 2:
            r = ww - 3
            for s in range(sig - r, sig + r + 1):
                for d in cells_row(s, r - abs(s - sig)):
                    put(s, d, Gam(s, d))
        # half-gap rows for own edges and children's edges
        edges = [(e, ww), (e + LL, ww)]
        for (k, ek) in kids:
            edges += [(ek, w(k)), (ek + L(k), w(k))]
        halfgap(edges, inside)
        # children
        for (k, ek) in kids:
            ordinal(k, ek)
        # border
        for s in range(e, min(N, e+LL+1)):
            if (s - ww) % 2 == 0:
                put(s, ww, 0); put(s, -ww, 0)
        for s in (e, e+LL):
            for d in cells_row(s, ww):
                put(s, d, 0)
    # main ordinals and global half-gap rows
    n = 1
    mains = []
    while e_main(n) < N:
        mains.append(n); n += 1
    edges = []
    for n in mains:
        edges += [(e_main(n), w(n)), (e_main(n) + L(n), w(n))]
    def outside_all(s, d):
        for n in mains:
            if e_main(n) <= s <= e_main(n) + L(n) and abs(d) <= w(n):
                return False
        return True
    halfgap(edges, outside_all)
    for n in mains:
        ordinal(n, e_main(n))
    # special cells near origin and the boundary
    put(2, 0, 2); put(4, 0, 1)
    for s in range(N):
        put(s, s, 1); put(s, -s, 1)
        if s >= 3:
            put(s, s-2, 0); put(s, -(s-2), 0)
    return D

if __name__ == '__main__':
    import sys
    from gen import gen
    N = int(sys.argv[1])            # compare on rows s < N
    D = paint(N)
    M = N//2 + 2
    t = gen(M+2)
    bad = 0
    for s in range(N):
        for d in range(-s, s+1, 2):
            i = (s-d)//2; j = (s+d)//2
            if i > M or j > M: continue
            if D[s][d+N] != t[i][j]:
                bad += 1
                if bad < 20: print("mismatch at s=%d d=%d (i=%d,j=%d): D=%s t=%d" % (s, d, i, j, D[s][d+N], t[i][j]))
    print("mismatches:", bad)
