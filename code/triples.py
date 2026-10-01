import sys
from gen import gen, F
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
t = gen(N)
seen = {}
for i in range(1, N):
    for j in range(1, N):
        key = (t[i-1][j], t[i-1][j-1], t[i][j-1])
        if key not in seen:
            seen[key] = (i, j)
print(len(seen), "triples (x,y,z)=(t(i-1,j),t(i-1,j-1),t(i,j-1)) occur; first occurrence, f value:")
for k in sorted(seen):
    print(k, '->', F[k[0]][k[1]][k[2]], 'first at', seen[k])
print("missing:", [(x,y,z) for x in range(3) for y in range(3) for z in range(3) if (x,y,z) not in seen])
