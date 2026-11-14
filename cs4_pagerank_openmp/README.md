# Case 4 — PageRank (OpenMP)

```bash
make
./serial 65536 1000000 42 degrees.csv
OMP_NUM_THREADS=4 ./parallel 65536 1000000 42 pull static 0 verify
OMP_NUM_THREADS=4 ./parallel 65536 1000000 42 push dynamic 64 verify
OMP_NUM_THREADS=4 ./parallel 65536 1000000 42 pull guided 256 noverify
```

Serial arguments: vertices, edges, seed, optional degree CSV. Parallel adds
`pull|push`, `static|dynamic|guided`, chunk, `verify|noverify`, optional degree CSV
after the seed. Vertices must be a power of two in [2,2^24]. Zero edges is valid.
Dynamic/guided require chunk >=1; static 0 selects contiguous static partitioning.
Final input: 65,536 vertices, 1,000,000 directed edges, seed 42.

R-MAT quadrant probabilities are (0.57,0.19,0.19,0.05). Self-loops and parallel
edges are retained as a directed multigraph. Both outgoing and incoming CSR are
built. These are synthetic data. Reference: Chakrabarti, Zhan and Faloutsos,
*R-MAT: A Recursive Model for Graph Mining*, SDM 2004,
https://doi.org/10.1137/1.9781611972740.43.

Damping is 0.85. Dangling mass is redistributed uniformly. All versions stop at
L1 residual <1e-10 or fail after 500 iterations. Pull owns destinations and reads
incoming neighbors; push traverses outgoing CSR with atomic destination additions.
Both precompute source contributions. Generation/CSR building is excluded and
reported separately; the entire iterative solve and synchronization are timed.

Verification runs serial pull outside timing: maximum absolute rank error <=1e-9,
L1 error <=1e-8, rank sum within 1e-8 of one, and convergence. Tests include
all-dangling graphs and small irregular graphs. Atomic/reduction results need not
be bitwise identical.

Thread work timing starts after a team barrier and stops before the loop barrier.
Imbalance is maximum accumulated edge-loop time divided by its mean across actual
workers. It includes scheduling overhead, excludes end-barrier waiting, and is
not whole-program imbalance. Degree CSV reports incoming degree; zero-degree
vertices are omitted only from the log-log plot.

Pull/push time differences include both atomics and different traversal/access
patterns, so they do not isolate atomic instruction cost. A rough pull traffic
model is 4 bytes per neighbor index plus 8 bytes per source contribution before
cache reuse, plus O(vertices) traffic for ranks, offsets and reductions.
