from math import comb
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from koszul import uniform_formula, derangement_poly
from rr import is_real_rooted, parallel_ext
# Binder–Vecchi: H^i_{U(r,n)} = I_{r-1}( D(f_i) ),  f_i = sum_k c_{k,i} x^{r-1-k}(1+x)^k,  D(x^m)=d_m(x)
def f_i(r,n,i):
    out=[0]*r
    for k in range(r):
        c = sum(comb(n-m-1, r+i-m-1)*comb(i+k-m-1, i-1) for m in range(k+1))
        for j in range(k+1):          # x^{r-1-k}(1+x)^k
            out[r-1-k+j] += c*comb(k,j)
    return out
def D(f):
    L = max(len(f),1)+1; out=[0]*(L+1)
    for m,cf in enumerate(f):
        for e,v in enumerate(derangement_poly(m)): out[e]+=cf*v
    return out
def I(g, deg):
    g = g + [0]*(deg+1-len(g)); return list(reversed(g[:deg+1]))
ok=True
for r in range(2,7):
    for n in range(r+1, r+8):
        for i in range(1, n-r+1):
            a = I(D(f_i(r,n,i)), r-1); b = uniform_formula(r,n,i)
            if a != b: ok=False; print("mismatch", r,n,i,a,b)
print("identity H^i = I_{r-1} D(f_i) verified:", ok)
# exact check: all parallel extensions of uniform matroids, all i (including Chow column)
bad=0; tot=0
for r in range(3,8):
    for n in range(r+1, r+9):
        H=[uniform_formula(r,n,i) for i in range(n-r+1)]
        for p in range(0,41):
            for h in parallel_ext(H,p):
                tot+=1
                if not is_real_rooted(h): bad+=1; print("FAIL",r,n,p)
print(f"uniform parallel extensions: {tot} polynomials checked exactly, failures={bad}")
