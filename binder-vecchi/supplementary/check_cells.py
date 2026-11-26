import argparse
import itertools
import sys
import zipfile
from collections import defaultdict
from math import gcd


class SparseMatrix:
    def __init__(self, shape):
        self.columns = [dict() for _ in range(shape[1])]

    def __getitem__(self, pos):
        i, j = pos
        return self.columns[j].get(i, 0)

    def __setitem__(self, pos, v):
        i, j = pos
        if v:
            self.columns[j][i] = v
        else:
            self.columns[j].pop(i, None)


def rational_rank(A):
    pivots = {}
    for col in sorted(A.columns, key=len):
        while col:
            i = max(col)
            if i not in pivots:
                pivots[i] = col
                break
            other = pivots[i]
            a, b = other[i], col[i]
            g = gcd(a, b)
            a, b = a // g, b // g
            if a != 1:
                col = {k: a*v for k,v in col.items()}
            for k, v in other.items():
                w = col.get(k, 0) - b*v
                if w:
                    col[k] = w
                else:
                    col.pop(k, None)
            if abs(a) > 1 and col:
                d = gcd(*col.values())
                if d > 1:
                    col = {k: v//d for k,v in col.items()}
    return len(pivots)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('zipfile')
    args = parser.parse_args()
    sys.path.insert(0, args.zipfile)
    from koszul import linear, graphic
    from sympy import Poly, discriminant, symbols
    x, t = symbols('x t')
    with zipfile.ZipFile(args.zipfile) as archive:
        src = archive.read('cellular.py').decode()
    replacements = {
        'np.zeros((dims[e - 1], dims[e]), dtype=np.int64)': 'SparseMatrix((dims[e - 1], dims[e]))',
        '(A[row, col] + sgn * cf) % P': 'A[row, col] + sgn * cf',
        'ranks[e] = rank_mod(A)': 'ranks[e] = rational_rank(A)',
    }
    for old, new in replacements.items():
        assert src.count(old) == 1
        src = src.replace(old, new)
    env = {'SparseMatrix': SparseMatrix, 'rational_rank': rational_rank}
    exec(compile(src, 'cellular_over_Q.py', 'exec'), env)
    examples = [linear([(1,) + v for v in itertools.product(range(2), repeat=3)], 2, 'AG(3,2)'),
                graphic([(0,1),(0,2),(0,3),(0,4),(1,2),(2,3),(3,4),(4,1)], 'W4'),
                graphic(list(itertools.combinations(range(5), 2)), 'K5')]
    for M in examples:
        H = env['refined_cellular'](M)
        polynomial = sum(c*x**d*t**i for i,h in enumerate(H) for d,c in enumerate(h))
        delta = Poly(discriminant(polynomial, x), t)
        print(M.name, H, 'disc_degree',delta.degree(),'disc_min',min(delta.all_coeffs()),flush=True)


if __name__ == '__main__':
    main()
