"""Verify three rank-four tables and columns 9, 10 and 11 of PG(3,2) over Q.

Usage: python verify_results.py binder_vecchi_verificacion.zip
Requires Python 3.10+, NumPy and SymPy. Reads the original ZIP without extracting
or modifying it. The ZIP supplies matroids, flat chains and annihilator bases.
This file constructs integer boundary matrices, computes their ranks over Q,
and computes the discriminants symbolically.
"""

import argparse
import heapq
from hashlib import sha256
from itertools import combinations, product
from math import comb, gcd
from pathlib import Path
import sys

from sympy import Poly, discriminant, symbols


ARCHIVE_SHA256 = '114e09e9652c7f0e319574ea77d9cfe73c7be932d3652fa59ce05f989a04bcb0'


def rational_rank(columns):
    """Exact elimination in Q, choosing a sparse row or column to limit fill."""
    row_count = 1 + max((max(col, default=-1) for col in columns), default=-1)
    rows = [set() for _ in range(row_count)]
    for j, col in enumerate(columns):
        for i in col:
            rows[i].add(j)

    # Nonnegative keys denote columns. Negative keys denote rows via ~i.
    queue = [(len(col), j) for j, col in enumerate(columns) if col]
    queue.extend((len(row), ~i) for i, row in enumerate(rows) if row)
    heapq.heapify(queue)
    rank = 0
    while queue:
        size, key = heapq.heappop(queue)
        if key < 0:
            i = ~key
            if len(rows[i]) != size:
                continue
            j = min(rows[i], key=lambda k: len(columns[k]))
        else:
            j = key
            if len(columns[j]) != size:
                continue
            i = min(columns[j], key=lambda k: len(rows[k]))

        pivot = columns[j]
        for r in pivot:
            rows[r].remove(j)
        for k in list(rows[i]):
            col = columns[k]
            a, b = pivot[i], col[i]
            divisor = gcd(a, b)
            a, b = a // divisor, b // divisor
            if a != 1:
                for r in col:
                    col[r] *= a
            for r, v in pivot.items():
                w = col.get(r, 0) - b * v
                if w:
                    col[r] = w
                    rows[r].add(k)
                else:
                    col.pop(r, None)
                    rows[r].discard(k)
            if abs(a) > 1 and col:
                content = gcd(*col.values())
                if content > 1:
                    for r in col:
                        col[r] //= content
            if col:
                heapq.heappush(queue, (len(col), k))
        for r in pivot:
            if rows[r]:
                heapq.heappush(queue, (len(rows[r]), ~r))
        columns[j] = {}
        rank += 1
    return rank


def wedge_image(subset, basis_map):
    """Expand the exterior product in the target basis, including signs."""
    terms = {0: 1}
    for index in subset:
        result = {}
        for mask, coeff in terms.items():
            for j, sign in basis_map[index]:
                bit = 1 << j
                if mask & bit:
                    continue
                swaps = (mask >> (j + 1)).bit_count()
                value = -coeff * sign if swaps % 2 else coeff * sign
                key = mask | bit
                result[key] = result.get(key, 0) + value
        terms = result
    return {key: value for key, value in terms.items() if value}


def homology_entries(matroid, targets):
    from cellular import chains_of_flats, perp_basis

    n = matroid.n
    flats, chains = chains_of_flats(matroid)
    bases = {ch: perp_basis(matroid, flats, ch) for layer in chains for ch in layer}
    chain_indices = [{ch: i for i, ch in enumerate(layer)} for layer in chains]
    ranks = {}

    def boundary_rank(e, a):
        if e < 1 or e >= len(chains) or n - 1 - e < a:
            return 0
        if (e, a) in ranks:
            return ranks[e, a]
        source_wedges = list(combinations(range(n - 1 - e), a))
        target_wedges = {sum(1 << j for j in subset): i for i, subset in
                         enumerate(combinations(range(n - e), a))}
        columns = []
        for ch in chains[e]:
            local = [{} for _ in source_wedges]
            for omit in range(e):
                face = ch[:omit] + ch[omit + 1:]
                # In a block with minimum m, the basis vector at u is e_u-e_m.
                # Thus e_u-e_v is the difference of the corresponding coordinates,
                # with the minimum's coordinate interpreted as zero.
                coordinates = {u: j for j, (u, _) in enumerate(bases[face])}
                basis_map = [tuple((coordinates[k], s) for k, s in ((u, 1), (v, -1))
                                   if k in coordinates) for u, v in bases[ch]]
                offset = chain_indices[e - 1][face] * len(target_wedges)
                orientation = (-1) ** omit
                for column, subset in zip(local, source_wedges):
                    for mask, value in wedge_image(subset, basis_map).items():
                        column[offset + target_wedges[mask]] = orientation * value
            columns.extend(local)
        ranks[e, a] = rational_rank(columns)
        return ranks[e, a]

    entries = {}
    for i, e in targets:
        a = n - 1 - e - i
        dimension = len(chains[e]) * comb(n - 1 - e, a)
        entries[i, e] = dimension - boundary_rank(e, a) - boundary_rank(e + 1, a)
    return entries, ranks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    archive = parser.parse_args().archive.resolve()
    if sha256(archive.read_bytes()).hexdigest() != ARCHIVE_SHA256:
        raise SystemExit('The archive differs from the reviewed ZIP.')
    sys.path.insert(0, str(archive))
    from geoms import PG
    from koszul import graphic, linear

    wheel = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4), (4, 1)]
    cases = [
        (linear([(1,) + v for v in product(range(2), repeat=3)], 2, 'AG(3,2)'), 14, 2816000),
        (graphic(wheel, 'M(W4)'), 16, 70848),
        (graphic(list(combinations(range(5), 2)), 'M(K5)'), 22, 1562500),
    ]
    t, x = symbols('t x')
    for matroid, degree, minimum in cases:
        targets = product(range(matroid.n - matroid.r + 1), range(matroid.r))
        entries, _ = homology_entries(matroid, targets)
        polynomial = sum(h * t**i * x**e for (i, e), h in entries.items())
        delta = Poly(discriminant(polynomial, x), t)
        if (delta.degree(), min(delta.all_coeffs())) != (degree, minimum):
            raise ValueError(f'Unexpected discriminant for {matroid.name}')
        print(matroid.name, 'coefficients of x by row i:', flush=True)
        for i in range(matroid.n - matroid.r + 1):
            print(i, [entries[i, e] for e in range(matroid.r)])
        print('Discriminant degree and smallest coefficient:', degree, minimum, flush=True)

    expected = {(11, e): v for e, v in enumerate([0, 0, 0, 64])}
    expected.update({(10, e): v for e, v in enumerate([0, 0, 245, 516])})
    expected.update({(9, e): v for e, v in enumerate([0, 35, 2380, 1759])})
    entries, ranks = homology_entries(PG(4, 2), expected)
    if entries != expected or (ranks[1, 4], ranks[2, 4]) != (1001, 45439):
        raise ValueError('Unexpected PG(3,2) entries or boundary ranks')
    print('PG(3,2) entries indexed by (i, exponent):', entries)
    print('Ranks of boundaries 1 and 2 at exterior degree 4:', ranks[1, 4], ranks[2, 4])
    print('Discriminant of H^9(x)/x:', entries[9, 2]**2 - 4*entries[9, 1]*entries[9, 3])


if __name__ == '__main__':
    main()
