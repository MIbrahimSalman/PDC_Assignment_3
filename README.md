# Group 04 — Parallel and Distributed Computing, Assignment 3

Fatima Malik — 23L-0826  
Ibrahim Salman — 23L-0723

## Status

This project includes five parallel case studies, measured results, a report and
a six-slide presentation. See [CONTRIBUTIONS.md](CONTRIBUTIONS.md) for contribution
status and the consolidated provenance disclosure.

| Case | Model | Current state |
|---|---|---|
| 1: Heat diffusion | MPI | Serial; 1-D/2-D decomposition; blocking/overlap; verification; phase timers |
| 2: k-mer counting | OpenMP | Serial; intentional race; atomic/private; full histogram verification; false-sharing experiment |
| 3: Monte Carlo | OpenMP | Reproducible counter streams; ordered reduction; standard error; Black–Scholes reference |
| 4: PageRank | OpenMP | R-MAT CSR; dangling nodes; pull/push; static/dynamic/guided; verification |
| 5: Sample sort | MPI | Distributed shards; regular samples; Alltoallv; oversampling; exact small-input oracle; ping-pong |

This allocation satisfies the requirement of at least two OpenMP and two MPI cases.
C++17 is used, which the assignment permits. Source files therefore use `.cpp`.

## Build and check in Ubuntu WSL

Open Ubuntu, then:

```bash
cd '/mnt/c/Users/Fatima/Downloads/pdc assignment 3/Group04_Assignment3'
make
make check
```

Dependencies: `g++`, `make`, and MPI development tools (`mpic++`, `mpicc`, `mpirun`).
These were found in the existing Ubuntu distribution. All baseline and parallel
programs use the same `-O3 -std=c++17` optimization settings. `make check` runs
small correctness experiments on 1, 2, and 4 workers, not assignment-scale benchmarks.
It assumes at least four MPI slots; on smaller machines, adjust the worker list.
Do not silently oversubscribe performance runs.

Each case README documents arguments and measurement boundaries. For timing, set:

```bash
export OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false
```

## Reproduce measurements and documents

```bash
make check
bash tools/benchmark.sh
python3 -m pip install -r tools/requirements-docs.txt
python3 tools/render_results.py
python3 tools/format_report.py
# Compile report/Group04_report.tex three times from the report directory.
python3 tools/package.py
```

The benchmark uses three repeats per point, 1/2/4/6 physical-core workers and an
additional 8-worker SMT Monte Carlo comparison. It includes the requested full
sizes: heat N=1024/4096 for 500 steps, 100 million DNA bases, 100 million paths,
one million graph edges and 100 million sorting keys. Expect a sustained run;
close other applications first. The recorded session did not guarantee host
idleness, so read the measurement caveats in the report.

The script resumes completed labelled repetitions, and compares source hashes to
avoid mixing code versions. Archive the results directory before a new machine or
source-version experiment. Python is only used for data analysis, plots, document
generation and packaging; all case-study algorithms are C++.

Artifacts:

- `report/Group04_report.pdf`: measured report with limitations and references.
- `slides/Group04_slides.pdf`: exactly six slides.
- `results/raw.tsv`, `summary.csv`, `figures/`: raw runs, median tables and figures.
- `verification_all.log`: all 244 pre-benchmark comparisons.
- `verification_final.log`: 11 additional comparisons and three known-value checks.
- `clean_build.log`: successful source-only archive build and five smoke checks.
- `results/environment.txt`, `host_machine.json`, `sources.sha256`: environment and provenance.
- `../Group04_Assignment3.zip`: portable source/evidence archive; binaries excluded.

The completed measurement set contains 303 workload timings (101 groups, three
repeats each), plus three ping-pong sweeps. The LaTeX report is 15 pages and the deck is
six slides. The source-only archive was built and smoke-tested from a fresh extraction.

The root Makefile builds every case and `make check` runs both verification scripts.
No performance claim is based on the deliberately incorrect demonstrations.

## Group workflow

Create a private shared Git repository and have both members commit their own work
over time. No remote repository or fabricated member commits have been created.
Record actual authorship, reviews, and hours in `CONTRIBUTIONS.md`.
