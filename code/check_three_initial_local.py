"""Finite local catalogue for the three-initial-symbol proof.

All values come from the explicit pattern definition, not a recurrence grid.
The manuscript proves that these cases cover arbitrary levels.
"""

from functools import lru_cache
from three_initial import auxiliary_rule, description, length, opening


value = lru_cache(None)(description)
OFFSETS = ((0, 0), (-1, 1), (-2, 0), (-1, -1))


def instances(n):
    """Main, first/other left child, last/other right child."""
    near, far = n + 3, n + 4
    offset = length(n + 1) + 6
    return (
        opening(n),
        opening(near) + 3,
        opening(far) + offset,
        opening(near) + length(near) - length(n) - 3,
        opening(far) + length(far) - length(n) - offset,
    )


def window(s, d):
    return tuple(value(s + ds, d + dd) for ds, dd in OFFSETS)


def correct(values):
    return values[0] == auxiliary_rule(*values[1:])


def central(n, u, d):
    a = abs(d)
    if a == n:
        return 0
    if a == 0 and abs(u) == n + 3:
        return 2
    if a == 0 and abs(u) == n + 1:
        return 0
    if abs(u) + a <= n - 1:
        return 1 if d % 2 == 0 else 2
    return 1 if (abs(u) - a - n + 1) % 4 == 0 else 2


def central_meeting(n, u, d):
    for du, dd in OFFSETS:
        a, v = abs(d + dd), u + du
        if abs(v) + a <= n - 1 or (a == 0 and abs(v) == n + 1):
            return True
    return False


def central_window(n, u, d):
    return tuple(central(n, u + du, d + dd) for du, dd in OFFSETS)


def central_cases(n):
    for u in range(-n - 1, n + 4):
        for d in range(n + 1):
            if (n + 1 + u - d) % 2 == 0 and central_meeting(n, u, d):
                yield u, d


def main():
    edge_count = wall_count = central_count = compared = 0
    catalogue = set()
    for n in range(1, 20):
        for e in instances(n):
            for edge in (e, e + length(n)):
                for offset in range(-2, 5):
                    s = edge + offset
                    for d in range(n + 2):
                        if (s - d) % 2 or d > s - 2:
                            continue
                        w = window(s, d)
                        assert correct(w), ("edge", n, e, s, d, w)
                        catalogue.add(w)
                        edge_count += 1
            rows = (range(length(n) + 1) if n <= 5 else
                    tuple(range(6, 10)) + tuple(range(length(n) - 9, length(n) - 5)))
            for r in rows:
                for d in (n - 1, n, n + 1):
                    s = e + r
                    if (s - d) % 2 or d > s - 2:
                        continue
                    w = window(s, d)
                    assert correct(w), ("wall", n, e, r, d, w)
                    catalogue.add(w)
                    wall_count += 1
    for n in range(3, 16):
        for u, d in central_cases(n):
            w = central_window(n, u, d)
            assert correct(w), ("central", n, u, d, w)
            central_count += 1
            catalogue.add(w)
            for e in instances(n):
                actual = window(e + length(n)//2 + u, d)
                assert actual == w, ("central template", n, e, u, d, actual, w)
                compared += 1
    for s, d in ((2, 0), (3, 1), (4, 0)):
        assert correct(window(s, d))
    print(f"Horizontal-edge catalogue: {edge_count} ordered windows, no failures.")
    print(f"Wall catalogue: {wall_count} ordered windows, no failures.")
    print(f"Central catalogue, levels 3..15: {central_count} ordered windows, no failures.")
    print(f"Central formula compared with the five geometric contexts: {compared} matches.")
    print(f"The catalogue contains {len(catalogue)} distinct ordered value quadruples.")
    print("Initial exceptional cell and its successors: all three nonnegative-position cases verified.")

    shrink_edges = shrink_centres = 0
    for n in range(20, 33):
        small = 16 + (n - 16) % 4
        for e, target in zip(instances(n), instances(small)):
            for end in (False, True):
                edge = e + (length(n) if end else 0)
                small_edge = target + (length(small) if end else 0)
                for offset in range(-2, 5):
                    for d in range(n + 2):
                        if (edge + offset - d) % 2:
                            continue
                        if d <= n - 9:
                            image = d if d <= 1 else 2 + (d - 2) % 4
                        else:
                            image = d - (n - small)
                        assert window(edge + offset, d) == window(small_edge + offset, image), (
                            "edge reduction", n, e, end, offset, d, small, image)
                        shrink_edges += 1
        for u, d in central_cases(n):
            points = [(u + du, d + dd) for du, dd in OFFSETS]
            if all(abs(v) + abs(a) <= n - 1 for v, a in points):
                continue
            if d >= 6:
                v, a = u, d - 4
            else:
                assert abs(u) >= 9
                v, a = u - (4 if u > 0 else -4), d
            assert central_window(n, u, d) == central_window(n - 4, v, a), (
                "central reduction", n, u, d, v, a)
            shrink_centres += 1
    print(f"Edge reductions at levels 20..32: {shrink_edges} ordered windows agree.")
    print(f"Central reductions at levels 20..32: {shrink_centres} ordered windows agree.")


if __name__ == "__main__":
    main()
