# fig4.py -- figures of the four-letter sequence tau in rotated coordinates (time s from left to right,
# position d from top to bottom; column p shows the row s = 2p for even d and s = 2p+1 for odd d).
from pathlib import Path
from png import write_png
ROOT = Path(__file__).resolve().parent.parent
TAB = "1023022312311302010120102300301020002120300301123102300010332020"
G = [int(c) for c in TAB]
N = 460
t = [[1] * N for _ in range(N)]
for x in range(1, N):
    for y in range(1, N):
        t[x][y] = G[(t[x - 1][y] * 4 + t[x - 1][y - 1]) * 4 + t[x][y - 1]]
PAL = {0: (255, 255, 255), 1: (220, 30, 30), 2: (40, 40, 220), 3: (30, 160, 30)}
def rotated(p0, p1, D):
    rows = []
    for d in range(-D, D + 1):
        line = []
        for p in range(p0, p1):
            s = 2 * p + (d & 1)
            line.append(t[(s - d) // 2][(s + d) // 2])
        rows.append(line)
    return rows
# the main box S_3 of block 0, opening at s = 209, closing at 367 (half-width 13)
write_png(ROOT / 'fig_four_S3.png', rotated(102, 186, 16), pal=PAL, scale=8)
# the beginning of block 0: bar at 47, virtual box V_0 from 51 to 209, S_3 from 209 to 367, S_4 opening at 371
write_png(ROOT / 'fig_four_block0.png', rotated(20, 188, 24), pal=PAL, scale=4)
print("ok")
