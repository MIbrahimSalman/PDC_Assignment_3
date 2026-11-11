# Case 2 — k-mer counting (OpenMP)

```bash
make
export OMP_NUM_THREADS=4 OMP_PROC_BIND=close OMP_PLACES=cores
./serial 100000000 8 42 uniform
./parallel 100000000 8 42 uniform atomic verify
./parallel 100000000 8 42 uniform private verify
./parallel 1000000 8 42 skew race verify
./false_sharing 100000000 packed
./false_sharing 100000000 padded
```

Serial arguments: `L k seed distribution`. Parallel adds `mode verification`.
Defaults: `100000000 8 42 uniform private verify`.
Distributions: `uniform|skew`; modes: `race|atomic|private`; verification:
`verify|noverify`. k supports 1..12; final experiments require k=8 and k=12.
L must be at least k. Use 1,2,4,... threads up to the physical-core limit.

DNA bases are represented as one-byte codes 0..3. SplitMix64 generates reproducible
input independently of thread count. Skew uses approximately 90% A. The histogram
uses 64-bit counts and a rolling two-bit encoding; it takes O(L) time rather than
recomputing all k characters of every window. Threads own disjoint intervals of
window **start positions**, reading up to k-1 bases past their interval to count
boundary-spanning windows exactly once. Empty intervals are valid.

Version A deliberately contains a data race, which is undefined behavior in C++.
It often loses updates, especially on skewed data; it is never a correct variant
even if one run happens to match. A failed comparison returns exit status 2. Versions
B and C compare every bin exactly to the serial output, plus the total L-k+1.
Validation is excluded from timing. It can be disabled explicitly for benchmarks,
but then the count check alone does not prove histogram correctness.

Timing excludes input generation (reported separately), shared histogram allocation,
and validation. Private histogram allocation, initialization, and merging are included
in parallel time. For 64-bit counters, private storage alone is:

| k | Bins | Memory per worker |
|---|---:|---:|
| 8 | 65,536 | 512 KiB |
| 12 | 16,777,216 | 128 MiB |

Multiply by configured maximum threads; also allow for the shared histogram, input,
and a reference histogram during verification. Avoid swapping in memory benchmarks.

The separate false-sharing microbenchmark uses adjacent private counters versus
`alignas(64)` counters. It assumes a 64-byte cache line; verify the target CPU's
cache-line size. Volatile forces repeated memory accesses without introducing a
shared-counter race. Run both modes repeatedly, alternate order, pin threads, and
report actual medians; improvement is hardware-dependent, not guaranteed.

For analysis, report atomic/private times for k=8 and 12, extra private memory,
lost counts from A, and millions of k-mers/second. The rolling scan reads about one
new input byte per window after each chunk's first k bytes, but histogram traffic,
cache-line fetches, atomic coherence, zeroing, and merging add traffic. Do not treat
input bytes alone as measured total memory bandwidth. Full assignment input is
L >= 100,000,000; smaller inputs are correctness tests only.
