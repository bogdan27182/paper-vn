"""Exact checks of the arithmetic formula, red-cell gaps, and origin cutoff."""

from functools import lru_cache
from pathlib import Path
import re

from three_initial import children, grid, length, opening, pattern

ROOT = Path(__file__).resolve().parent.parent


@lru_cache(None)
def odd_edges(n):
    positions = {0, length(n)} if n % 2 else set()
    for k, offset in children(n):
        positions.update(offset + r for r in odd_edges(k))
    return tuple(sorted(positions))


def main():
    counts = [0, 1, 1, 1]
    for n in range(4, 201):
        counts.append(counts[n-1] + 2*counts[n-3])
    for base in (3, 8):
        for n in range(1, 201):
            value = (base**(n*n+2*n) // (base**(3*n)-base**(2*n)-2)) % base**n
            assert value == counts[n], (base, n, value, counts[n])
    print('Arithmetic terms, bases 3 and 8: all n=1..200 agree with the recurrence.')

    for n in range(1, 19):
        positions = tuple(r for r in range((1-n) % 2, length(n)+1, 2)
                          if pattern(n, r, 1) == 0)
        assert positions == odd_edges(n), n
        if n >= 3:
            largest = max(b-a for a, b in zip(positions, positions[1:]))
            edge_distance = 0 if n % 2 else 3
            assert (positions[0], positions[-1], largest) == (
                edge_distance, length(n)-edge_distance, 2*n+8), n
    print('Levels 1..18: red positions, endpoint distances, and maximum gaps verified.')

    size = 1024
    cells = grid(size, red=True)
    actual = {2*x+1 for x in range(size-1) if cells[x][x+1] == 1}
    expected = set()
    n = 1
    while opening(n) <= 2*size-3:
        expected.update(opening(n)+r for r in odd_edges(n) if opening(n)+r <= 2*size-3)
        n += 1
    assert actual == expected
    print(f'Direct {size}x{size} recurrence: all {len(actual)} red cells on y=x+1 agree.')

    source = (ROOT / 'code' / 'onesym.h').read_text()
    table = tuple(map(int, re.search(r'TAB = "([0-3]{64})"', source).group(1)))
    size = 160
    cells = [bytearray([1])*size for _ in range(size)]
    for x in range(1, size):
        for y in range(1, size):
            cells[x][y] = table[16*cells[x-1][y] + 4*cells[x-1][y-1] + cells[x][y-1]]

    def value(s, d):
        return cells[(s-d)//2][(s+d)//2]

    candidates = []
    for s in range(1, 58):
        for width in range(1, s+1):
            if (s-width) % 2 or any(value(s, d) for d in range(-width, width+1, 2)):
                continue
            first_break = next((r for r in range(s+2, 200, 2)
                                if value(r, width) or value(r, -width)), None)
            assert first_break is not None, (s, width)
            candidates.append((s, width, first_break))
    assert len(candidates) == 53
    assert max(r for _, _, r in candidates) == 111
    certificate = ROOT / 'data' / 'origin_wall_breaks.csv'
    certificate.write_text('opening,half_width,first_nonwhite_wall_time\n' +
                           ''.join(f'{s},{w},{r}\n' for s, w, r in candidates))
    print('Four-letter origin: 53 candidate openings; every wall breaks by time 111.')
    print(f'Certificate: {certificate.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
