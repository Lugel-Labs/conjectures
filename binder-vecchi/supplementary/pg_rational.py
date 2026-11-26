import sys, zipfile, time
from check_cells import SparseMatrix, rational_rank

start = time.monotonic()
path = sys.argv[1]
sys.path.insert(0, path)
from geoms import PG

with zipfile.ZipFile(path) as z:
    source = z.read('cellular.py').decode()

# H_e at e=1 uses only the chain spaces in degrees 0, 1 and 2.
for old, new in [
    ('D = r - 1', 'D = 2'),
    ('for a in range(0, n):', 'for a in [4]:'),
    ('np.zeros((dims[e - 1], dims[e]), dtype=np.int64)', 'SparseMatrix((dims[e - 1], dims[e]))'),
    ('(A[row, col] + sgn * cf) % P', 'A[row, col] + sgn * cf'),
    ('ranks[e] = rank_mod(A)', 'print("matrix ready",e,dims,round(time.monotonic()-start),flush=True)\n            ranks[e] = rational_rank(A)\n            print("rank",e,ranks[e],round(time.monotonic()-start),flush=True)'),
    ('for e in range(D + 1):\n            hom =', 'for e in [1]:\n            hom =')]:
    assert source.count(old) == 1
    source = source.replace(old,new)
env = {'SparseMatrix':SparseMatrix, 'rational_rank':rational_rank, 'time':time, 'start':start}
exec(compile(source,'cellular_exact_target.py','exec'),env)
H = env['refined_cellular'](PG(4,2))
print('PG(3,2) rational [x] H^9',H[9][1],flush=True)
