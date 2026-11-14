# Case 5 — Distributed sample sort (MPI)

```bash
make
./serial 100000000 4 uniform 42
mpirun --bind-to core -np 4 ./parallel 100000000 uniform 1 42 noverify
mpirun --bind-to core -np 4 ./parallel 100000000 skew 8 42 noverify
mpirun -np 4 ./parallel 10003 skew 8 42 verify
mpirun -np 2 ./pingpong > pingpong.csv
```

Serial arguments: N, shard count, `uniform|skew`, seed. Parallel: N, distribution,
oversampling factor, seed, `verify|noverify`. Defaults: N=10^8, seed=42, uniform,
oversampling=1. N must fit an MPI int displacement (<=INT_MAX). Pass the parallel
rank count as the serial shard count to recreate the same per-rank-seeded input.
Thus the recorded serial baseline is matched separately to each p.

Each rank generates its shard, sorts locally, chooses s*(p-1) regular samples
(capped by local size), gathers samples, broadcasts quantile splitters, partitions
using upper_bound, exchanges counts with Alltoall and keys with Alltoallv, then
sorts received keys with std::sort. Only sample metadata is gathered in production.
No process is required to hold the complete production input for this algorithm.

Skew produces about 90% zero keys. Equal splitters intentionally route duplicates
to the same bucket, exposing basic sample sort's severe imbalance. Oversampling
cannot split identical keys with this partition rule. A future improvement would
break ties using global indices. An overloaded rank may receive most of this
adversarial input; this is a measured limitation, not a correctness error.

Every run checks local order, boundaries between nonempty ranks, total count,
key sum and mixed 64-bit hash sum. Hashes are probabilistic multiset checks, not
an exact proof. `verify` additionally gathers at most one million keys for an
exact serial comparison; it is rejected above that size. Generation and validation
are excluded from sorting time.

Elapsed MPI time is between barriers and maximized over ranks. The four phase
times come from the rank with greatest total elapsed time, allowing a coherent
stacked chart; ending-barrier waiting is outside these four phases.
`max_exchange_seconds` is independently maximized. Exchange timing includes counts,
receive-buffer allocation and Alltoallv, so it exceeds pure transfer time.

Ping-pong warms each size and reports mean one-way time as half the round-trip
mean; the benchmark repeats it three times. Fit T=alpha+beta*bytes, state the
fitted size regime and goodness of fit. A balanced Alltoallv approximation is
(p-1)*alpha + beta*4*N*(p-1)/p^2. Collective algorithms, shared-memory transport,
allocation and imbalance limit this model.
