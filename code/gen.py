# t(i,0)=t(0,j)=1, t(i,j)=f(t(i-1,j), t(i-1,j-1), t(i,j-1)) over F_3
def f(x, y, z):
    return (x*x*y*y + y*y*z*z + y*y + x*z + 2*y + 2) % 3

F = [[[f(x, y, z) for z in range(3)] for y in range(3)] for x in range(3)]

def gen(N):
    t = [bytearray(N) for _ in range(N)]
    for i in range(N):
        t[i][0] = 1
        t[0][i] = 1
    for i in range(1, N):
        prev = t[i-1]; cur = t[i]
        for j in range(1, N):
            cur[j] = F[prev[j]][prev[j-1]][cur[j-1]]
    return t
