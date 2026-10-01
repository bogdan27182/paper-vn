"""Three initial symbols: recurrence rules and box geometry.

State 0 is the border colour for the auxiliary rule. The paper uses the conjugate
rule under a -> 1-a (mod 3), so that its borders are red (state 1).
"""

from functools import lru_cache


def auxiliary_rule(x, y, z):
    return (2 + 2*y + 2*y*y + 2*y*y*(x*x + z*z)
            + 2*x*z*(x + z) + x*x*z*z) % 3


def red_rule(x, y, z):
    return (1 - auxiliary_rule((1-x) % 3, (1-y) % 3, (1-z) % 3)) % 3


TABLE = [auxiliary_rule(x, y, z)
         for x in range(3) for y in range(3) for z in range(3)]


def grid(size, red=False):
    table = ([red_rule(x, y, z)
              for x in range(3) for y in range(3) for z in range(3)]
             if red else TABLE)
    boundary = 0 if red else 1
    cells = [bytearray([boundary]) * size for _ in range(size)]
    for x in range(1, size):
        row, previous = cells[x], cells[x-1]
        for y in range(1, size):
            row[y] = table[9*previous[y] + 3*previous[y-1] + row[y-1]]
    return cells


@lru_cache(None)
def length(n):
    if n <= 3:
        return (0, 2, 6, 14)[n]
    return length(n-1) + 2*length(n-3) + 8


@lru_cache(None)
def opening(n):
    return 5 if n == 1 else opening(n-1) + length(n-1) + 3


@lru_cache(None)
def children(n):
    first = []
    offset = 3
    for k in range(n-3, 0, -1):
        first.append((k, offset))
        offset += length(k) + 3
    return tuple(first + [(k, length(n)-length(k)-a) for k, a in reversed(first)])


@lru_cache(None)
def edges(n):
    return ((0, n), (length(n), n)) + tuple(
        (edge, k) for k, a in children(n) for edge in (a, a+length(k)))


def pattern(n, r, d):
    """The box pattern with border colour 0."""
    span, a = length(n), abs(d)
    if r in (0, span) or a == n:
        return 0
    for k, offset in children(n):
        if offset <= r <= offset+length(k) and a <= k:
            return pattern(k, r-offset, d)
    for edge, width in edges(n):
        distance = abs(r-edge)
        if 1 <= distance <= 2 and a <= width-distance:
            return 2
    if n >= 3:
        delta = abs(r-span//2)
        if d == 0 and delta == n+1:
            return 0
        if delta+a <= n-1:
            return 1 if d % 2 == 0 else 2
    return 1 if (min(r, span-r)+a-n+2) % 4 == 0 else 2


def description(s, d):
    if (s, d) == (2, 0):
        return 0
    a = abs(d)
    n = 1
    while opening(n) <= s+2:
        e, span = opening(n), length(n)
        if e <= s <= e+span and a <= n:
            return pattern(n, s-e, d)
        for edge in (e, e+span):
            distance = abs(s-edge)
            if 1 <= distance <= 2 and a <= n-distance:
                return 2
        n += 1
    return 1 if (s-a) % 4 == 0 else 2


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--size", type=int, default=1200)
    parser.add_argument("--rows", type=int, default=0)
    args = parser.parse_args()
    cells = grid(args.size)
    bad = []
    for x, row in enumerate(cells):
        for y, value in enumerate(row):
            if description(x+y, y-x) != value:
                if len(bad) < 10:
                    bad.append((x, y, value, description(x+y, y-x)))
    if bad:
        raise AssertionError(bad)
    print(f"Description agrees with the recurrence on {args.size} x {args.size} cells.")
    if args.rows:
        for s in range(2, args.rows):
            for d in range(-(s % 2), 45, 2):
                if d > s-2:
                    continue
                got = description(s, d)
                want = auxiliary_rule(description(s-1,d+1), description(s-2,d),
                                     description(s-1,d-1))
                if got != want:
                    raise AssertionError((s,d,got,want))
        print(f"Local recurrence verified for s < {args.rows}, 0 <= d <= 44.")
