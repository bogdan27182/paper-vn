"""Enumerate rectangles and verify the three-initial-symbol box counts."""

from math import comb

from three_initial import children, grid, length, opening

def main():
    size, width = 1536, 24
    cells, limit = grid(size), 2*size-2
    def value(s, d):
        return cells[(s-d)//2][(s+d)//2]
    found = set()
    for e in range(1, limit):
        for w in range(1, min(width,e)+1):
            if (e-w) % 2 or e+w > limit:
                continue
            if any(value(e,d) for d in range(-w,w+1,2)):
                continue
            span = 2
            while e+span+w <= limit and not value(e+span,w) and not value(e+span,-w):
                if all(value(e+span,d) == 0 for d in range(-w,w+1,2)):
                    found.add((e,w,span))
                span += 2
    expected = set()
    def collect(n, e):
        if e > limit:
            return
        if n <= width and e+length(n)+n <= limit:
            expected.add((e,n,length(n)))
        for k, offset in children(n):
            collect(k,e+offset)
    n = 1
    while opening(n) <= limit:
        collect(n,opening(n))
        n += 1
    assert found == expected, (found-expected, expected-found)
    print(f'All zero-bordered rectangles with half-width <= {width} in a {size} x {size} recurrence grid: '
          f'{len(found)}; exactly the predicted boxes.')
    verified = []
    for n in range(1,20):
        e, span = opening(n), length(n)
        if e+span+n > limit:
            break
        count = sum(e <= a and a+l <= e+span and w <= n for a,w,l in found)
        formula = sum(2**j*comb(n-1-2*j,j) for j in range((n-1)//3+1))
        assert count == formula, (n,count,formula)
        verified.append((n,count))
    print('Box counts from the rectangle enumeration:', verified)


if __name__ == '__main__':
    main()
