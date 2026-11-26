"""Explicit counterexample to Binder–Vecchi Conjecture 1.5 in degree i = 1.

M = Tr^7(M_47 + M_47 + U_{1,1}), simple, rank 7 (Bérczi–Coron–Eur–Li, Thm C).
N = M with p parallel elements added.

Lemma 1 (Binder, Prop. 3.6) with i = 1:   H^1_N = p * H_M + H^1_M.

Rouché: if |H^1_M| < p |H_M| on a circle C, then H^1_N and p*H_M have the
same number of zeros inside C. We pick C around a non-real zero of H_M,
away from the real axis, so H^1_N has a non-real zero.

Bound on H^1_M, with no need to compute it: its coefficients are the graded
pieces of Tor_1 of the face ring of the Bergman complex over Sym N^vee.
The complex is shellable, so reducing by a linear system of parameters
gives an Artinian algebra B with dim B = number of facets F, and Tor_1 is a
subquotient of B^(n-r). Hence the coefficients of H^1_M are >= 0, sum to at
most T = (n - r) * F, and sit in degrees 0..r-1. On |z| <= R:
    |H^1_M(z)| <= T * max(1, R)^(r-1).

Every step below is exact rational arithmetic.
"""
import sys
from fractions import Fraction as Q
from math import isqrt

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
from chow_bcel import R as RANK, chow, sum_flag  # noqa: E402

q = 47
alpha, n = sum_flag(q)
H = chow(alpha)                      # coefficients, increasing degree
F = alpha[tuple(range(1, RANK))]     # facets = maximal chains of proper flats
T = (n - RANK) * F


def ev(z):  # exact complex evaluation, z = (re, im) Fractions
    re, im = Q(0), Q(0)
    for c in reversed(H):
        re, im = re * z[0] - im * z[1] + c, re * z[1] + im * z[0]
    return re, im


# Locate the non-real zeros numerically (only to choose the circle).
import mpmath  # noqa: E402
mpmath.mp.dps = 80
roots = mpmath.polyroots(H[::-1], maxsteps=500, extraprec=400)
z0 = max((r for r in roots if abs(mpmath.im(r)) > 0), key=lambda r: mpmath.im(r))
others = [r for r in roots if r != z0]
sep = min(abs(z0 - r) for r in others)
cx = Q(str(mpmath.nstr(mpmath.re(z0), 30)))
cy = Q(str(mpmath.nstr(mpmath.im(z0), 30)))
rho = Q(str(mpmath.nstr(min(sep, abs(mpmath.im(z0))) / 4, 6)))

# Lipschitz constant of H on the closed disk: |H'| <= sum k|h_k| Rb^(k-1).
Rb = abs(cx) + abs(cy) + rho          # >= |z| on the disk
L = sum(k * abs(h) * Rb ** (k - 1) for k, h in enumerate(H) if k)

# Sample the circle at exact rational points; u in [-1,1], both halves.
def circle_point(u, s):
    d = 1 + u * u
    return cx + s * rho * (1 - u * u) / d, cy + s * rho * 2 * u / d

N_SAMPLES = 4096
min_sq = None
values = []  # H at consecutive points, once around the circle
for s in (1, -1):
    for k in range(N_SAMPLES + (s == -1)):
        u = Q(2 * k, N_SAMPLES) - 1
        re, im = ev(circle_point(u, s))
        values.append((re, im))
        v = re * re + im * im
        min_sq = v if min_sq is None else min(min_sq, v)


def quadrant_steps(a, b):
    """Signed quarter-turns from a to b; exact, valid when the turn is < pi/2."""
    def quad(z):
        x, y = z
        return 0 if x > 0 and y >= 0 else 1 if x <= 0 and y > 0 else 2 if x < 0 and y <= 0 else 3
    d = (quad(b) - quad(a)) % 4
    assert d != 2, "step too large"
    return {0: 0, 1: 1, 3: -1}[d]

# Between samples, arc length <= rho * 2 * (Delta u) / 2 to the nearest sample.
gap = rho * Q(2, N_SAMPLES)
lower = Q(isqrt(int(min_sq.numerator * 10**40 // min_sq.denominator)), 10**20) - L * gap
assert lower > 0, "refine the sampling"
# Consecutive samples are <= 2*gap apart along the arc. If L*2*gap < min|H|,
# H on each arc stays in a disk around H(z_k) that excludes 0 and subtends
# less than pi/2 seen from 0, so counting quadrant crossings gives the exact
# winding number = number of zeros of H_M inside the circle.
sqrt_min = lower + L * gap
assert L * 2 * gap * 2 < sqrt_min, "refine the sampling"  # angle < pi/6 per step
turns = sum(quadrant_steps(values[i], values[(i + 1) % len(values)]) for i in range(len(values)))
assert turns % 4 == 0
zeros_inside = turns // 4
assert zeros_inside >= 1, "no zero inside the circle"
assert cy - rho > 0, "circle must stay off the real axis"

bound_H1 = T * max(Q(1), Rb) ** (RANK - 1)
p0 = bound_H1 // lower + 1

print("q =", q, "| rank 7 | n =", n)
print("H_M coefficients:", H)
print("facets F =", F)
print("Tor_1 bound T = (n - r) F =", T)
print("circle centre =", float(cx), "+", float(cy), "i | radius =", float(rho))
print("min |H_M| on circle >=", float(lower))
print("zeros of H_M inside the circle (exact winding number):", zeros_inside)
print("p0 =", p0)
print("p0 has", len(str(p0)), "digits")
print("Counterexample in degree 1: M with p >= p0 parallel elements;"
      " ground set size n + p.")
