"""Generate overview figures from the recurrence, with exact binary masks.

Run from the repository root with: python3 code/overview_figures.py
Each grid cell is drawn as a solid block of pixels without interpolation.
"""

from pathlib import Path
import re

from gen import gen
from png import PAL, write_png
from three_initial import grid as three_initial_grid

ROOT = Path(__file__).resolve().parent.parent


def four_letter_grid(size):
    source = (ROOT / "code" / "onesym.h").read_text()
    table = list(map(int, re.search(r'TAB = "([0-3]{64})"', source).group(1)))
    grid = [bytearray([1]) * size for _ in range(size)]
    for x in range(1, size):
        previous, row = grid[x - 1], grid[x]
        for y in range(1, size):
            row[y] = table[(previous[y] * 4 + previous[y - 1]) * 4 + row[y - 1]]
    return grid


def paired_figures(stem, grid, palette, border=0, scale=5):
    write_png(str(ROOT / f"{stem}_colour.png"), grid, pal=palette, scale=scale)
    binary_palette = {value: (255, 255, 255) if value == border else (0, 0, 0)
                      for value in palette}
    write_png(str(ROOT / f"{stem}_mono.png"), grid, pal=binary_palette, scale=scale)
    print(f"{stem}: {len(grid)} x {len(grid[0])} cells; white in mask = state {border}")


if __name__ == "__main__":
    paired_figures("fig_three_corner", gen(128), PAL)
    paired_figures("fig_three_initial_corner", three_initial_grid(64, red=True), PAL, border=1,
                   scale=10)
    paired_figures("fig_four_corner", four_letter_grid(192),
                   {**PAL, 3: (30, 160, 30)})
