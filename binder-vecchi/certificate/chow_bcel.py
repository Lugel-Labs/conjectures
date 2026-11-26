"""Exact Chow polynomial of Tr^7(M_q + M_q + U_{1,1}) (Bérczi–Coron–Eur–Li, Theorem C).

Inputs are the flag f-vector of the rank-4 matroid M_q (arXiv:2610.17135,
Corollary 4.11). Everything is exact integer arithmetic.

Chow polynomial (Feichtner–Yuzvinsky basis):
    H_M(t) = sum over chains of flats  empty < F_1 < ... < F_k  (F_k may be E),
             each weighted by prod_j w(rk F_j - rk F_{j-1}),
    w(d) = t + t^2 + ... + t^(d-1).
"""
from itertools import combinations, product
from math import comb

R = 7  # rank of the truncation


def mq_flag(q):
    m = (q * q + 1) // 2
    N1 = q**4
    N2 = q**3 * (q**3 + q**2 + q + 1)
    N3 = q**2 * (q**2 + 1) * (q**2 + q + 1)
    W1 = N1 + m
    W2 = N2 + N1 * m + comb(m, 2)
    W3 = N3 + N1 * comb(m, 2) + comb(m, 3)
    a = {
        (): 1,
        (1,): W1, (2,): W2, (3,): W3,
        (1, 2): q * N2 + 2 * N1 * m + 2 * comb(m, 2),
        (1, 3): q * q * N3 + m * q * q * (q * q + 1) + 3 * N1 * comb(m, 2) + 3 * comb(m, 3),
        (2, 3): (q * q + q) * N3 + m * q**4 * (q * q + 1) + 3 * N1 * comb(m, 2) + 3 * comb(m, 3),
        (1, 2, 3): q * (q * q + q) * N3 + 2 * m * q**4 * (q * q + 1) + 6 * N1 * comb(m, 2) + 6 * comb(m, 3),
    }
    n = q**4 + m
    return a, n


def multichain_count(seq, flag, top):
    """Multichains in a geometric lattice with the given nondecreasing rank sequence."""
    ranks = tuple(sorted({x for x in seq if 0 < x < top}))
    return flag[ranks]


def nondecreasing(k, lo, hi):
    if k == 0:
        yield ()
        return
    for first in range(lo, hi + 1):
        for rest in nondecreasing(k - 1, first, hi):
            yield (first,) + rest


def sum_flag(q):
    """Flag f-vector of the truncation, for rank sets S within {1..6}."""
    a, n_mq = mq_flag(q)
    alpha = {}
    for k in range(0, R):
        for S in combinations(range(1, R), k):
            total = 0
            for sa in nondecreasing(k, 0, 4):
                for sc in nondecreasing(k, 0, 1):
                    sb = tuple(s - x - y for s, x, y in zip(S, sa, sc))
                    if any(not 0 <= v <= 4 for v in sb):
                        continue
                    if any(sb[i] > sb[i + 1] for i in range(k - 1)):
                        continue
                    total += multichain_count(sa, a, 4) * multichain_count(sb, a, 4)
            alpha[S] = total
    n = 2 * n_mq + 1
    return alpha, n


def chow(alpha):
    coeffs = [0] * R  # degree r-1 = 6
    chains = [(S, c, (0,) + S) for S, c in alpha.items()]
    chains += [(S, c, (0,) + S + (R,)) for S, c in alpha.items()]
    for S, count, ranks in chains:
        poly = [1]
        ok = True
        for lo, hi in zip(ranks, ranks[1:]):
            d = hi - lo
            if d < 2:
                ok = False
                break
            w = [0] + [1] * (d - 1)  # t + ... + t^(d-1)
            new = [0] * (len(poly) + len(w) - 1)
            for i, x in enumerate(poly):
                for j, y in enumerate(w):
                    new[i + j] += x * y
            poly = new
        if ok:
            for i, x in enumerate(poly):
                coeffs[i] += count * x
    return coeffs


if __name__ == "__main__":
    import sys
    import numpy as np
    for q in map(int, sys.argv[1:]):
        alpha, n = sum_flag(q)
        H = chow(alpha)
        roots = np.roots(H[::-1])
        nonreal = [z for z in roots if abs(z.imag) > 1e-9 * max(1, abs(z))]
        print(q, "n =", n, "H =", H, "palindromic:", H == H[::-1],
              "nonreal roots:", len(nonreal))
