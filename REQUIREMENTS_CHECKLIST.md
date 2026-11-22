# Requirements and evidence

This checklist separates technical artifacts from actions only the members can
complete. See `CONTRIBUTIONS.md` for the consolidated provenance disclosure.

| Requirement | Evidence / status |
|---|---|
| C/C++ only for case algorithms | Five C++17 implementations; Python is analysis/document tooling |
| At least two MPI and two OpenMP cases | MPI: heat and sort; OpenMP: k-mers, Monte Carlo, PageRank |
| Separate serial baselines, same optimization | `serial.cpp`, `parallel.cpp`, per-case Makefiles; `common.mk` uses -O3 |
| Correctness comparisons | `make check`, `verification_all.log`; additional `tools/check_final.sh` |
| Heat 1-D/2-D, blocking/nonblocking, reductions | CLI options; Cartesian topology and halo implementations |
| Heat strong/weak scaling, K=1/50 | `heat*` experiments in `results/raw.tsv`; report section 1 |
| k-mer race/atomic/private and chunk handling | All variants; full bin comparisons in verification scripts |
| k=8/12, L>=10^8, false sharing | Benchmark script, memory table and measured packed/padded comparison |
| Monte Carlo M>=10^8, price/SE/closed form | Fixed parameters and counter streams; report section 3 |
| RNG failure and reproducibility demonstrations | `same_stream`; fixed-block deterministic reduction; serial comparisons |
| p=1/8 fixed-time accuracy comparison | Calibrated twelve-second workloads, M>=10^8; actual elapsed-time mismatch disclosed |
| PageRank >=10^6 edges, CSR, dangling nodes | R-MAT generator and incoming/outgoing CSR |
| Pull/push, schedules/chunks, degree/imbalance | CLI options; schedule sweep; `degrees.csv`; work-loop timers |
| Distributed sample sort N>=10^8 | Rank-local generation, sampling, splitters, Alltoall/Alltoallv |
| Sort skew/oversampling, order/count preservation | Bucket counts, hashes, small exact serial oracle |
| Sort phase chart and ping-pong model | `phases` graph; three ping-pong CSVs; linear fit in report |
| Three repeats, median, serial speedup/efficiency | `results/raw.tsv`, `summary.csv`, labelled figures |
| Machine and toolchain | `results/environment.txt`, `host_machine.json`, source hashes |
| Report PDF (<=16 pages) | `report/Group04_report.pdf`; generated and checked before packaging |
| Exactly six slides | `slides/Group04_slides.pdf`; generated and checked before packaging |
| Portable source ZIP | `Group04_Assignment3.zip`; binaries excluded; SHA-256 alongside |
| Real member hours, authored files, commit evidence | Pending; user requested skipping for now |
| Private shared Git repository/history | Pending; cannot fabricate retrospective contributions |
| Individual peer evaluations | Members must submit privately; not generated |
| Both members' presentation/viva | Members must rehearse and demonstrate understanding |
| Compliance with Section 11 | See provenance disclosure in `CONTRIBUTIONS.md` |

Measurement caveats: the Windows host was not proven idle; measurements use a
single WSL laptop, not multiple networked nodes. Full-scale sort checks use order,
count and probabilistic fingerprints rather than gathering all keys. Exact
serial comparisons use smaller data. Eight Monte Carlo workers use SMT. Per-case
READMEs and the report state timing boundaries, tolerances and model limitations.
