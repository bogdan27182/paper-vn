# Check the completeness claim of Section 6.5: every environment cell with an irregular parent, or with d=0,
# lies in one of the classes (E1)-(E6).
from describe import L, w, e_main, children, Lam, Gam
from zones import labels
N = 2400
val, lab = labels(N)
# box tree with host lookup
boxes = []
def rec(n, e, par):
    idx = len(boxes); boxes.append((n, e, par))
    for (k, ek) in children(n, e): rec(k, ek, idx)
n = 1
while e_main(n) < N:
    rec(n, e_main(n), -1); n += 1
def host(s, d):
    best = None
    for idx, (n, e, par) in enumerate(boxes):
        if e < s < e + L(n) and abs(d) < w(n):
            if best is None or w(n) < w(boxes[best][0]): best = idx
    return best
def kids_of(h):
    if h is None: return [(n, e) for (n, e, par) in boxes if par == -1]
    n, e, _ = boxes[h]; return children(n, e)
def classes(s, d, h):
    cl = []
    if d == 0: cl.append('E1')
    if abs(d) == s - 4: cl.append('E6')
    if h is not None:
        n, e, _ = boxes[h]; ww = w(n); mu = e + L(n)//2
        if abs(d) == ww - 1: cl.append('E2in')
        if s == e + 3 and abs(d) == ww - 3: cl.append('E4J1')
        if s == e + 4 and abs(d) == ww - 2: cl.append('E4J1')
        if n >= 2 and s > mu and (s - mu) + abs(d) == ww - 1: cl.append('E5')
    ks = kids_of(h)
    for q, (k, ek) in enumerate(ks):
        if abs(d) == w(k) + 1 and ek + 1 <= s <= ek + L(k) + 1: cl.append('E2out')
        if s == ek + L(k) + 2 and abs(d) == w(k): cl.append('E3')
        if q + 1 < len(ks) and ks[q+1][1] == ek + L(k) + 4 and ks[q+1][0] == k - 1 and s == ek + L(k) + 3 and abs(d) == w(k) - 1:
            cl.append('E4J2')
    return cl
unclassified = []; stats = {}; special = []
for (s, d), l in lab.items():
    if not l.startswith('env') or abs(d) > s - 2 or abs(d) > 40: continue
    h = host(s, d)
    eps = 0 if l == 'env0' else 2
    irr = [p for p in ((s-1, d+1), (s-2, d), (s-1, d-1)) if val[p] != Lam(eps, *p)]
    if not irr and d != 0: continue
    cl = classes(s, d, h); special.append((s, d))
    for c in cl: stats[c] = stats.get(c, 0) + 1
    if not cl: unclassified.append((s, d, l, irr))
print("special environment cells by class:", stats)
print("distinct special cells:", sum(1 for _ in special))
print("unclassified:", len(unclassified), unclassified[:10])
