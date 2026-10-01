"""Label every cell (rotated coords) by the rule of the explicit description that defines it,
and tabulate which argument triples occur in which zone (cross-check of the case analysis)."""
import sys
from describe import L, w, e_main, children, Lam, Gam
from gen import gen

def labels(N):
    lab = {}
    val = {}
    def put(s, d, v, l):
        if 0 <= s < N and abs(d) <= s:
            val[(s, d)] = v; lab[(s, d)] = l
    def row(s, h):
        st = -h
        if (s - st) % 2: st += 1
        return range(st, h+1, 2)
    for s in range(N):
        for d in row(s, s):
            put(s, d, Lam(2, s, d), 'env2')
    def halfgap(edges, inside):
        for (e, a) in edges:
            for s, h, kind in ((e+1, a-1, 'G1a'), (e-1, a-1, 'G1b'), (e+2, a-2, 'G2a'), (e-2, a-2, 'G2b')):
                if h < 0: continue
                for d in row(s, h):
                    if inside(s, d):
                        if (s, d) in lab and lab[(s, d)].startswith('G') and lab[(s, d)] != kind:
                            kind2 = lab[(s, d)] + '+' + kind
                        else:
                            kind2 = kind
                        put(s, d, Gam(s, d) if kind[1] == '1' else 1, kind2)
    def ordinal(n, e):
        if e >= N: return
        ww, LL = w(n), L(n); sig = e + LL//2
        kids = children(n, e)
        inside = lambda s, d: e < s < e + LL and abs(d) < ww
        for s in range(e+1, min(N, e+LL)):
            for d in row(s, ww-1):
                put(s, d, Lam(0 if s < sig else 2, s, d), 'env0' if s < sig else 'env2')
        if n >= 2:
            r = ww - 3
            for s in range(sig-r, sig+r+1):
                for d in row(s, r - abs(s-sig)):
                    put(s, d, Gam(s, d), 'diamond')
        edges = [(e, ww), (e+LL, ww)]
        for (k, ek) in kids: edges += [(ek, w(k)), (ek+L(k), w(k))]
        halfgap(edges, inside)
        for (k, ek) in kids: ordinal(k, ek)
        for s in range(e, min(N, e+LL+1)):
            if (s - ww) % 2 == 0:
                put(s, ww, 0, 'wall'); put(s, -ww, 0, 'wall')
        for d in row(e, ww): put(e, d, 0, 'open-corner' if abs(d) == ww else 'open-edge')
        for d in row(e+LL, ww):
            if abs(d) < ww: put(e+LL, d, 0, 'close-edge')
    mains = []; n = 1
    while e_main(n) < N: mains.append(n); n += 1
    edges = [(3, 1)]
    for n in mains: edges += [(e_main(n), w(n)), (e_main(n)+L(n), w(n))]
    def outside_all(s, d):
        return all(not (e_main(n) <= s <= e_main(n)+L(n) and abs(d) <= w(n)) for n in mains)
    halfgap(edges, outside_all)
    for n in mains: ordinal(n, e_main(n))
    for s in range(N):
        put(s, s, 1, 'line-red'); put(s, -s, 1, 'line-red')
        if s >= 3: put(s, s-2, 0, 'line-white'); put(s, -(s-2), 0, 'line-white')
    return val, lab

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    val, lab = labels(N)
    t = gen(N//2 + 3)
    # sanity: val == t
    for (s, d), v in val.items():
        i, j = (s-d)//2, (s+d)//2
        if i < len(t) and j < len(t) and t[i][j] != v: print("MISMATCH", s, d); break
    from collections import defaultdict
    tab = defaultdict(set)
    for (s, d), v in val.items():
        if abs(d) > s - 2: continue
        x = val.get((s-1, d+1)); y = val.get((s-2, d)); z = val.get((s-1, d-1))
        tab[lab[(s, d)]].add((x, y, z, v))
    for l in sorted(tab):
        print("%-12s" % l, ' '.join('%d%d%d>%d' % q for q in sorted(tab[l])))
