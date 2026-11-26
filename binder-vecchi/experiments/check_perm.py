# Check: sum_i (-1)^i t^i H^i_M(t) = sum_{J chain of NON-flats} (-t)^{|J|} prod_{blocks} A_{|b|}(t)
import itertools
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from koszul import *
from cellular import refined_cellular
def eulerian(m):
    if m == 0: return [1]
    A=[0]*m
    for w in itertools.permutations(range(m)):
        A[sum(1 for k in range(m-1) if w[k]>w[k+1])]+=1
    return A
EUL={m:eulerian(m) for m in range(0,9)}
def pmul(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j]+=x*y
    return out
def check(M):
    n=M.n; levels=M.flats(); flats=set(F for L in levels for F in L)
    nonflats=[S for S in range(1,(1<<n)-1) if S not in flats]
    total=[0]*(2*n)
    # enumerate chains of nonflats via DFS
    def dfs(chain, last):
        # contribution
        blocks=[]; prev=0
        for S in chain: blocks.append(S&~prev); prev=S
        blocks.append(((1<<n)-1)&~prev)
        poly=[1]
        for b in blocks: poly=pmul(poly, EUL[bin(b).count('1')])
        sgn=(-1)**len(chain)
        for e,v in enumerate(poly): total[e+len(chain)]+=sgn*v
        for S in nonflats:
            if (last & S)==last and S!=last and (last==0 or True):
                if last==0 or ((last & S)==last and S!=last):
                    dfs(chain+[S], S)
    dfs([],0)
    H=refined_cellular(M)
    lhs=[0]*(2*n)
    for i,h in enumerate(H):
        for d,v in enumerate(h): lhs[i+d]+=(-1)**i*v
    while lhs and lhs[-1]==0: lhs.pop()
    while total and total[-1]==0: total.pop()
    return lhs==total, lhs, total
for M in [uniform(3,4), uniform(3,5), uniform(4,6), graphic([(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],"K4")]:
    ok,l,t = check(M); print(M.name, ok, l if not ok else "")
