# Which cells violate the environment lemmas once the transition squares and switch squares are removed?
from describe import L, w, e_main, children, Lam, Gam
from gen import gen
N = 2100
t = gen(N)
def T(s, d): return t[(s-d)//2][(s+d)//2]
boxes = []   # (n, e, parent_index, half)  half: 0 first half of parent, 2 second half / global
def rec(n, e, par):
    idx = len(boxes); boxes.append((n, e, par))
    for (k, ek) in children(n, e):
        rec(k, ek, idx)
n = 1
while e_main(n) + L(n) < 2*N - 10:
    rec(n, e_main(n), -1); n += 1
nmain = n - 1
# squares (rotated diamonds): list of (center_s, radius)
sq = []
for n in range(1, nmain):        # main junctions S_n -> S_{n+1}
    c = e_main(n) + L(n); sq.append((c+3, w(n)+3, 'Tmain'))
for (n, e, par) in boxes:
    kids = children(n, e)
    for q in range(len(kids)-1):
        (k1, e1), (k2, e2) = kids[q], kids[q+1]
        c = e1 + L(k1)
        if e2 - c != 4: continue
        if k2 == k1 - 1: sq.append((c+1, w(k1)+1, 'Tdec'))
        elif k2 == k1 + 1: sq.append((c+3, w(k1)+3, 'Tinc'))
    if n >= 2:
        sig = e + L(n)//2; sq.append((sig, w(n)-3, 'switch'))
def in_sq(s, d):
    return [k for (cs, r, k) in sq if abs(s-cs) + abs(d) <= r]
def region_box(s, d):
    best = None
    for idx, (n, e, par) in enumerate(boxes):
        if e <= s <= e + L(n) and abs(d) <= w(n):
            if best is None or w(n) < w(boxes[best][0]): best = idx
    return best
viol = {}
S = 2*N - 2
for s in range(2, min(S, e_main(nmax := nmain) + L(nmax))):
    for d in range(-min(s-4, 40), min(s-4, 40)+1):
        if (s - d) % 2: continue
        b = region_box(s, d)
        if b is not None:
            n, e, par = boxes[b]
            if s in (e, e+L(n)) or abs(d) == w(n): continue   # border
            eps = 0 if s < e + L(n)//2 else 2
            key = 'inside S%d' % n
        else:
            eps = 2; key = 'global'
        if T(s, d) != Lam(eps, s, d) and not in_sq(s, d):
            viol.setdefault(key, []).append(((s-d)//2, (s+d)//2, T(s, d)))
for k, v in viol.items(): print(k, len(v), v[:12])
# structure of all squares: odd rows entirely white or entirely red, even rows = central lattice
bad = 0
for (cs, r, kind) in sq:
    if kind == 'switch': continue
    if cs + r >= S: continue
    for s in range(cs - r, cs + r + 1):
        h = r - abs(s - cs)
        vals = [T(s, d) for d in range(-h, h+1) if (s-d) % 2 == 0]
        if s % 2 == 1:
            if not (all(v == 0 for v in vals) or all(v == 1 for v in vals)): bad += 1
        else:
            if any(T(s, d) != Gam(s, d) for d in range(-h, h+1) if (s-d) % 2 == 0): bad += 1
print("squares:", len(sq), "rows violating 'odd row all white or all red, even row central lattice':", bad)
print("non-switch squares:", sum(1 for q in sq if q[2] != 'switch'), " of which fully inside computed range:", sum(1 for (cs, r, k) in sq if k != 'switch' and cs + r < S), {k: sum(1 for q in sq if q[2] == k) for k in ('Tmain', 'Tdec', 'Tinc', 'switch')})
