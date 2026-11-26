import heapq
import sys
import time
from math import gcd

from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import verify_results as verify
sys.path.insert(0, sys.argv[1])
from geoms import PG

def rank(columns):
    t0 = time.monotonic()
    nr = max(max(c,default=-1) for c in columns) + 1
    rows = [set() for _ in range(nr)]
    for j, col in enumerate(columns):
        for i in col:
            rows[i].add(j)
    heap = [(len(c),j) for j,c in enumerate(columns) if c]
    heap.extend((len(row),~i) for i,row in enumerate(rows) if row)
    heapq.heapify(heap)
    def touch(k, size):
        if size: heapq.heappush(heap,(size,k))
    rk = 0
    print('start', nr,len(columns),flush=True)
    while heap:
        size, key = heapq.heappop(heap)
        if key < 0:
            i = ~key
            if len(rows[i]) != size: continue
            j = min(rows[i], key=lambda k:len(columns[k]))
        else:
            j = key
            if len(columns[j]) != size: continue
            i = min(columns[j], key=lambda k:len(rows[k]))
        pivot = columns[j]
        for r in pivot: rows[r].discard(j)
        for c in list(rows[i]):
            col = columns[c]
            a,b=pivot[i],col[i]
            g=gcd(a,b)
            a,b=a//g,b//g
            if a != 1:
                for r in col: col[r]*=a
            for r,v in pivot.items():
                w=col.get(r,0)-b*v
                if w:
                    col[r]=w
                    rows[r].add(c)
                else:
                    col.pop(r,None)
                    rows[r].discard(c)
            if abs(a)>1 and col:
                g=gcd(*col.values())
                if g>1:
                    for r in col: col[r]//=g
            touch(c,len(col))
        for r in pivot: touch(~r,len(rows[r]))
        columns[j]={}
        rk+=1
        if rk % 10000 == 0: print('pivots',rk,'min',size,'s',round(time.monotonic()-t0),flush=True)
    print('rank',rk,'s',round(time.monotonic()-t0),flush=True)
    return rk

verify.rational_rank=rank
print(verify.homology_entries(PG(4,2),[(9,2),(9,3),(9,0)]),flush=True)
