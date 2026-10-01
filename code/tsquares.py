# Check the transition squares between S_n and S_{n+1}.
from gen import gen
N = 4200
t = gen(N)
C = {0: 'W', 1: 'R', 2: 'B'}
for n in range(1, 11):
    a = 2**(n+2) - 2*n - 1; b = 2**(n+2) + 1
    if b >= N: break
    size = b - a + 1
    top = ''.join(C[t[a][a+m]] for m in range(size))          # top row, from top-left vertex
    left = ''.join(C[t[a+m][a]] for m in range(size))         # left column, from top-left vertex
    right = ''.join(C[t[a+m][b]] for m in range(size))        # right column, from top-right vertex down
    bottom = ''.join(C[t[b][a+m]] for m in range(size))       # bottom row, from bottom-left vertex right
    verts = [C[t[a][a]], C[t[a][b]], C[t[b][a]], C[t[b][b]]]
    # white cells inside must lie on anti-diagonal segments joining corresponding white edge cells
    ok_white = True
    for i in range(a, b+1):
        for j in range(a, b+1):
            if t[i][j] != 0: continue
            s = i + j
            m1 = s - 2*a          # initial edges: anti-diagonal through (a, a+m1) and (a+m1, a)
            m2 = 2*b - s          # final edges: through (b-m2, b) and (b, b-m2)
            on_init = 0 <= m1 < size and t[a][a+m1] == 0 and all(t[a+q][a+m1-q] == 0 for q in range(m1+1))
            on_fin = 0 <= m2 < size and t[b-m2][b] == 0 and all(t[b-m2+q][b-q] == 0 for q in range(m2+1))
            if not (on_init or on_fin): ok_white = False
    # red cells: 4-connected components must be staircases joining corresponding edge cells
    seen = set(); comps = []
    for i in range(a, b+1):
        for j in range(a, b+1):
            if t[i][j] == 1 and (i, j) not in seen:
                st = [(i, j)]; seen.add((i, j)); comp = []
                while st:
                    p = st.pop(); comp.append(p)
                    for q in ((p[0]+1, p[1]), (p[0]-1, p[1]), (p[0], p[1]+1), (p[0], p[1]-1)):
                        if a <= q[0] <= b and a <= q[1] <= b and q not in seen and t[q[0]][q[1]] == 1:
                            seen.add(q); st.append(q)
                comps.append(comp)
    stair_info = []
    ok_red = True
    for comp in comps:
        cs = set(comp)
        onedge = sorted(p for p in comp if p[0] in (a, b) or p[1] in (a, b))
        # walk the path from the end on the top row / right column
        ends = [p for p in comp if sum((q in cs) for q in ((p[0]+1, p[1]), (p[0]-1, p[1]), (p[0], p[1]+1), (p[0], p[1]-1))) <= 1]
        if len(comp) == 1:
            stair_info.append(('single', comp[0])); continue
        if len(ends) != 2: ok_red = False; stair_info.append(('notpath', len(comp))); continue
        start = min(ends, key=lambda p: (p[0], -p[1]))   # the end nearer the top row
        path = [start]; prev = None
        while True:
            p = path[-1]
            nxt = [q for q in ((p[0]+1, p[1]), (p[0]-1, p[1]), (p[0], p[1]+1), (p[0], p[1]-1)) if q in cs and q != prev]
            if not nxt: break
            prev = p; path.append(nxt[0])
        if len(path) != len(comp): ok_red = False
        runs = []; cur = 1
        dirs = [(path[k+1][0]-path[k][0], path[k+1][1]-path[k][1]) for k in range(len(path)-1)]
        for k in range(1, len(dirs)):
            if dirs[k] == dirs[k-1]: cur += 1
            else: runs.append(cur+1 if not runs else cur); cur = 1
        runs.append(cur + (1 if not runs else 0))
        e1, e2 = path[0], path[-1]
        pair = (e1[0] == a and e2[1] == a and e1[1]-a == e2[0]-a) or (e1[1] == b and e2[0] == b and b-e1[0] == b-e2[1])
        stair_info.append(('init' if e1[0] == a else 'final', e1, e2, pair, dirs[0], dirs[1] if len(dirs) > 1 else None))
        if not pair: ok_red = False
    print("n=%d square [%d,%d]^2 size %d  vertices %s" % (n, a, b, size, verts))
    print("   top   :", top); print("   left  :", left); print("   right :", right); print("   bottom:", bottom)
    print("   white cells on joining segments:", ok_white, "  red comps:", len(comps), " all staircases joining corresponding cells:", ok_red)
    for x in stair_info[:3]: print("     ", x)

# step lengths of every staircase (cells per straight run)
print("\nstep-length check")
for n in range(1, 10):
    a = 2**(n+2) - 2*n - 1; b = 2**(n+2) + 1
    red = {(i, j) for i in range(a, b+1) for j in range(a, b+1) if t[i][j] == 1}
    seen = set(); allruns = set()
    for p0 in sorted(red):
        if p0 in seen: continue
        comp = set(); st = [p0]; seen.add(p0)
        while st:
            p = st.pop(); comp.add(p)
            for q in ((p[0]+1, p[1]), (p[0]-1, p[1]), (p[0], p[1]+1), (p[0], p[1]-1)):
                if q in red and q not in seen: seen.add(q); st.append(q)
        ends = [p for p in comp if sum(q in comp for q in ((p[0]+1, p[1]), (p[0]-1, p[1]), (p[0], p[1]+1), (p[0], p[1]-1))) == 1]
        if len(comp) == 1: allruns.add((1,)); continue
        p = min(ends); path = [p]; prev = None
        while True:
            nxt = [q for q in ((p[0]+1, p[1]), (p[0]-1, p[1]), (p[0], p[1]+1), (p[0], p[1]-1)) if q in comp and q != prev]
            if not nxt: break
            prev, p = p, nxt[0]; path.append(p)
        # split into maximal straight runs of cells
        runs = []; k = 0
        while k < len(path) - 1:
            dvec = (path[k+1][0]-path[k][0], path[k+1][1]-path[k][1]); l = 1
            while k + l < len(path) - 1 and (path[k+l+1][0]-path[k+l][0], path[k+l+1][1]-path[k+l][1]) == dvec: l += 1
            runs.append(l + 1); k += l
        allruns.add(tuple(runs))
    print(n, sorted(allruns, key=len))
