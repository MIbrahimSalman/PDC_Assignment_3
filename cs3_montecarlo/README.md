# Case 3 — Monte Carlo option pricing (OpenMP)

```bash
make
./serial 100000000 42
OMP_NUM_THREADS=4 ./parallel 100000000 42 good noverify
OMP_NUM_THREADS=4 ./parallel 1000003 42 good verify
OMP_NUM_THREADS=8 ./parallel 100000000 42 same_stream noverify
```

Serial arguments: samples, seed. Parallel adds `good|same_stream` and
`verify|noverify`. Samples must be between 2 and 10^12. Parameters are S0=100,
strike=100, risk-free rate=0.05, volatility=0.2, maturity=1 year. Both executables
print price, estimated standard error, analytic Black–Scholes price and z score.
The standard error uses the sample variance of discounted payoffs with M-1 in
the denominator. Statistical agreement is separate from numerical verification.

OpenMP fits this independent-sample workload on one shared-memory workstation.
Each path gets two distinct counters passed through a 64-bit SplitMix permutation
and transformed using Box–Muller. Assignment is by global path index, not worker
identity. Counter ranges never overlap for the allowed input sizes. This is a
deterministic pseudorandom scheme, not a claim of mathematical independence or
cryptographic quality. No shared mutable generator or rand() is used.

Workers write independent partial sums for fixed 65,536-path blocks. A final
ordered reduction combines block moments in long double. This explicitly chosen
reduction preserves the result even when the worker count changes. All sample
generation, payoff computation and reduction are timed; the small partial array
allocation is excluded. Verification repeats the serial estimator outside timing
and compares both price and standard error with absolute tolerance 1e-12.

`same_stream` deliberately repeats only ceil(M/p) counters across the samples.
The naive standard error counts duplicated observations as independent and
underestimates uncertainty by approximately sqrt(p) for even partitions. Do not
use it as a valid estimator or successful serial-equivalence test.

The benchmark uses M=10^8 with 1,2,4,6 workers; eight workers are additionally run
for the requested p=8 comparison and use SMT on this six-core CPU. A fixed-time
comparison calibrates M to a twelve-second target (at least 10^8 paths) using observed throughput and
records actual elapsed times. It is an approximate fixed-time comparison; report
the timing mismatch. Fixed blocks and reduction order make results independent
of p on this toolchain; a conventional floating-point reduction need not do so.
