from pathlib import Path
from gen import gen
from png import write_png
ROOT = Path(__file__).resolve().parent.parent
t = gen(300)
# Figure 1: S_5 and the transition squares around it, original coordinates (rows = x, columns = y)
a0, a1 = 50, 132
write_png(ROOT / 'fig_S5.png', [list(t[i][a0:a1]) for i in range(a0, a1)], scale=6)
# Figure 2: rotated ("brick") view: line p shows row s=2p (even d) and s=2p+1 (odd d); d horizontal
W = 12
rows = []
for p in range(58, 126):
    line = []
    for d in range(-W, W+1):
        s = 2*p + (d & 1)
        line.append(t[(s-d)//2][(s+d)//2])
    rows.append(line)
# transpose so that time runs left to right
cols = [list(r) for r in zip(*rows)]
write_png(ROOT / 'fig_S5_rot.png', cols, scale=8)
