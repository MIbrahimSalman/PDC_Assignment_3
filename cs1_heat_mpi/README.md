# Case 1 — Heat diffusion (MPI)

```bash
make
./serial 1024 500 50
mpirun -np 4 ./parallel 1024 500 50 1d blocking verify
mpirun -np 4 ./parallel 1024 500 50 2d overlap verify
```

Serial arguments: `N steps K`. Parallel arguments: `N steps K decomposition exchange verification`.
Defaults: `1024 500 50 2d overlap noverify`. Choices: `1d|2d`, `blocking|overlap`,
`verify|noverify`. N includes boundary cells; N >= 3, steps >= 1, K >= 1.

Initial interior temperature is zero. The top boundary is 100 (including its
corners); all other boundary cells are zero. Boundaries stay fixed. Alpha is 0.2.
Both solvers run the requested number of steps; convergence checks are diagnostic,
not an early-stopping criterion. A final check is made even if steps is not a multiple of K.

The interior `(N-2) x (N-2)` is partitioned with remainder blocks, so N need not
be divisible by rank count. 1-D uses a `p x 1` Cartesian grid; 2-D uses
`MPI_Dims_create` and `MPI_Cart_create`. Empty blocks are rejected. Ghost rows are
contiguous; ghost columns use a committed vector datatype. Non-blocking exchange
posts receives and sends, computes strict interior cells, waits, then updates border
cells. The old grid remains unmodified until every outstanding request completes.

Verification computes a serial reference on each rank outside the timed section
and compares every owned cell, with absolute tolerance `1e-12`. This intentionally
simple oracle costs O(N^2) extra memory per rank: use it on small correctness runs,
then use `noverify` for large benchmarks. No grid gathering is required for timing.

`seconds` includes time stepping, halo exchange, convergence reductions, and the
ending barrier; rank maximum is reported. Allocation, initial conditions, datatype
setup, validation, and printing are excluded. The three phase times are individually
maximized across ranks and need not add to elapsed time. `exchange_seconds` in overlap
mode is posting plus exposed wait time, **not** total asynchronous transfer duration.
Use the difference in total times to study overlap; do not interpret reduced wait
time alone as measured network time hidden.

## Required experiments

- Strong scaling: N=1024 and 4096, steps=500, workers=1,2,4,... up to machine limit.
- Compare all four decomposition/exchange combinations, and K=1 versus K=50.
- Weak scaling: choose N approximately `2 + (N_base-2)*sqrt(p)` so interior cells per
  rank stay approximately constant; record actual dimensions and any rounding.
- Repeat every point at least three times; report medians, speedup against the
  separate serial program at the same N/steps/K, and efficiency=speedup/p.
- Interior work per rank is approximately `(N-2)^2/p`. Away from physical edges,
  sent halo values per step are approximately `2*(N-2)` for 1-D and
  `4*(N-2)/sqrt(p)` for square 2-D layouts. Receive traffic is the same again.
  Account for rectangular process grids and physical edges in exact counts.
- Locate efficiency below 70% using measured data. Compare with
  `T_parallel ~ T_serial/p + T_communication`; explain the model's limits.

The automated tests include uneven sizes and local blocks only one cell wide.
