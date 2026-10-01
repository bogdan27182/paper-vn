from gen import gen
from describe import L, w, e_main, children, Gam
t = gen(2100)
ok = True
def A(n, e): return ((e - w(n))//2, (e + w(n))//2)
def rec(n, e, depth):
    global ok
    u, v = A(n, e)
    kids = children(n, e)
    half = len(kids)//2
    for idx, (k, ek) in enumerate(kids):
        uk, vk = A(k, ek)
        if idx < half: pred = (u + 2**n - 2**(k+2) + 4, v + 2**n - 2**(k+2) - 2*(n-k) + 4)
        else:          pred = (u + 2**n + 2**(k+1) + 2*(n-k) - 4, v + 2**n + 2**(k+1) - 4)
        if pred != (uk, vk): ok = False; print("offset mismatch", n, e, k, pred, (uk, vk))
        rec(k, ek, depth+1)
    if n >= 2:
        a, b = u + 2**n, u + 2**n + 2*n - 4
        if b < 2100:
            # switching square: blue except central lattice; flanking diagonal cells red
            for x in range(a, b+1):
                for y in range(a, b+1):
                    exp = 1 if (x % 2 == 0 and y % 2 == 0) else 2
                    if t[x][y] != exp: ok = False; print("switch square bad", n, e, x, y)
            for q in (a-2, b+2):
                if t[q][q] != 1: ok = False; print("flank bad", n, e, q)
            if not (a % 2 == 0): ok = False
    # transition squares after each child that is followed by a sibling 4 rows later
    for q in range(len(kids)-1):
        (k1, e1), (k2, e2) = kids[q], kids[q+1]
        if e2 - (e1 + L(k1)) != 4: continue
        c = e1 + L(k1); uB = (c - w(k1))//2
        side = 2*k1 + 3 if k2 == k1 + 1 else 2*k1 + 1
        lo, hi = uB, uB + side - 1
        if hi >= 2100: continue
        for x in range(lo, hi+1):
            for y in range(lo, hi+1):
                s = x + y
                if s % 2 == 0:
                    if t[x][y] != (1 if x % 2 == 0 and y % 2 == 0 else 2): ok = False; print("sq even bad", x, y)
        for s in range(2*lo, 2*hi+1):
            if s % 2: 
                vals = {t[x][s-x] for x in range(lo, hi+1) if lo <= s-x <= hi}
                if len(vals) != 1 or vals == {2}: ok = False; print("sq odd row bad", s, vals)
n = 1
while e_main(n) + L(n) < 4190:
    rec(n, e_main(n), 0); n += 1
# main transition squares
for n in range(1, 10):
    lo, hi = 2**(n+2) - 2*n - 1, 2**(n+2) + 1
    if hi >= 2100: break
    for s in range(2*lo, 2*hi+1):
        if s % 2:
            vals = {t[x][s-x] for x in range(lo, hi+1) if lo <= s-x <= hi}
            if len(vals) != 1 or vals == {2}: ok = False; print("main sq bad", n, s, vals)
print("all coordinate formulas OK" if ok else "PROBLEMS")
