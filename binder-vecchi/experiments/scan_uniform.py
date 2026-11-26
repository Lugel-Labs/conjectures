import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from koszul import uniform_formula
from rr import *
import numpy as np, sys
# scan P_U(t,x) = sum_i t^i H^i_{U(r,n)}(x) for t>0 on a log grid
worst = []
for r in range(3, 9):
    for n in range(r+1, r+14):
        Hs = [uniform_formula(r, n, i) for i in range(n-r+1)]
        best = (0, None)
        for lt in np.linspace(-6, 6, 1201):
            t = 10**lt
            poly = [sum(t**i * Hs[i][d] for i in range(len(Hs))) for d in range(r)]
            mi = max_imag(poly)
            # normalized imaginary part relative to root size
            rts = np.roots(list(reversed(poly)))
            rel = float(np.max(np.abs(rts.imag)/(np.abs(rts)+1e-300)))
            if rel > best[0]: best = (rel, t)
        if best[0] > 1e-6:
            worst.append((r, n, best))
            print(f"U({r},{n}): P(t,x) NOT real-rooted near t={best[1]:.4g}, rel.imag={best[0]:.3g}")
            sys.stdout.flush()
print("done; failures:", len(worst))
