# Verification record

## LaTeX report formatting — 23 November 2026

The report now uses `report/Group04_report.tex`: 11-point A4 article, 2.2 cm
margins, linked contents, semantic colours, eight measured plots and five TikZ
schematics. All 109 original content blocks (54 paragraphs and 15 tables included)
were retained. No new measurements or citations were introduced. A note explains
the original whole-MiB rounding of the 512 KiB histogram allocation.

The final document was compiled with pdfLaTeX three times. It has 15 pages and
13 captioned, labelled, referenced figures. The final log has no overfull or
underfull boxes, unresolved references, or LaTeX/package warnings. All pages
were rasterized and visually inspected; unchanged page images were verified by
hash after the last local fixes. MiKTeX emitted a separate administrative notice
that its update check has not yet been run; compilation itself completed normally.
The six-slide PDF remains unchanged. Compilation logs and the content inventory
are included under `report/`.

## Full five-case verification — 23 November 2026

`make check` completed with 244 successful serial comparisons in Ubuntu WSL:
49 heat, 60 k-mer, 9 Monte Carlo, 90 PageRank and 36 sample-sort configurations.
See `verification_all.log`. This includes the updated Monte Carlo check of both
price and standard error, not just price. Invalid race demonstrations are excluded.

All five serial/parallel programs and the false-sharing/ping-pong diagnostics
build with GCC 13.3.0 and Open MPI 4.1.6. Final benchmark data are retained under
`results/`; the report and summary are generated from those measurements.
Host RAM was subsequently confirmed as 7,866,085,376 bytes through Windows CIM.
The earlier small-input-only statement below describes the first implementation
stage, not the current package.

Additional six-worker tests and million-edge PageRank checks passed after
benchmarking: 11 further serial comparisons, for 255 total, plus three independent
known-value checks (one-cell heat, Black–Scholes closed-form value, and an
all-dangling PageRank fixed point). See `verification_final.log`.

The source-only ZIP was extracted to a fresh `clean_build_check` directory and
compiled with its own Makefiles. All five case-study smoke comparisons passed;
the build emitted no compiler warnings. See `clean_build.log`. These five archive
smoke checks are additional to the 255 comparisons above.

303 workload timings are recorded across 101 experiment/worker groups, exactly
three repetitions each, plus three ping-pong sweeps. The final fixed-time Monte
Carlo experiment uses at least 100 million paths for both p=1 and p=8. Earlier
three-second exploratory runs are retained under a different label.

The generated report has 13 pages; the slide deck has exactly six. Selected dense
report pages and slides were visually inspected. Source/evidence packaging excludes
compiled executables and includes a SHA-256 checksum alongside the ZIP.

## Historical first-stage verification — 22 November 2026

Built in the existing Ubuntu WSL distribution with `make` using GCC and Open MPI
4.1.6. A follow-up build after disabling deprecated MPI C++ bindings completed
without compiler warnings. The implementation uses the MPI C API from C++.

`tools/check.sh` completed successfully; raw output is in `verification.log`:

- 49 MPI heat configurations passed against the serial reference, tolerance 1e-12.
- 60 OpenMP atomic/private k-mer configurations matched every serial histogram bin.
- Workers: 1, 2, 4. Included uneven partitions, tiny inputs, k=1/8/12, both input
  distributions, blocking/overlap, and convergence intervals 1/50.
- A further N=35, steps=50, K=1, four-rank 2-D overlap check had maximum error zero.
- Independent one-step heat sanity check: N=3, steps=1 has one interior cell at
  temperature 20, maximum change 20, and total grid checksum 320; serial output matched.
- Standalone serial k-mer smoke test: L=100003, k=8 reported 99996 windows as expected.
- Both false-sharing layouts executed successfully on four threads. These single
  short runs are not a reliable performance comparison.

The deliberately incorrect race demonstration, with four threads, L=1000000,
k=8, seed=42 and skewed input, returned total 844265 versus expected 999993 and
147 mismatched bins. This is one observed run, not a deterministic expected value.
The program correctly reported FAILED for that intentionally broken variant.

WSL-reported environment:

- AMD Ryzen 5 5500U with Radeon Graphics.
- 6 cores, 12 logical processors exposed; 2 hardware threads per core.
- Approximately 3.5 GiB RAM and 1.0 GiB swap visible to WSL. This is not a claim
  about the host's total installed RAM, which still needs to be recorded.
- L1d 192 KiB across six instances; L2 3 MiB across six; L3 4 MiB.

Final scaling results are **not** available. These checks establish small-input
correctness only. Run final experiments with at least three repetitions per point
and record actual compiler version, host RAM, affinity and WSL resource limits.
Powers-of-two worker counts up to physical cores are 1,2,4 on this machine. Eight
workers use SMT; label that explicitly rather than implying eight physical cores.
