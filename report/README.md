# Measured report

`Group04_report.pdf` is compiled from the single main source `Group04_report.tex`.
It uses A4 paper, 11-point text, 2.2 cm margins, a table of contents, semantic
colours, eight measured plots and five TikZ schematics. All 109 content blocks
from the preceding report are preserved in `content.json`; `formatting_audit.json`
records their inclusion. No new experiments or citations were introduced.

From the project root, regenerate the formatting and compile:

```bash
python tools/format_report.py
cd report
pdflatex -interaction=nonstopmode -halt-on-error Group04_report.tex
pdflatex -interaction=nonstopmode -halt-on-error Group04_report.tex
pdflatex -interaction=nonstopmode -halt-on-error Group04_report.tex
```

Keep `figures/` alongside the `.tex` file when compiling. The main project ZIP
contains these assets. Compilation logs are retained as `compile-pass1.txt`
through `compile-pass3.txt` and `Group04_report.log`.

`tools/render_results.py` produces the measurement summary and slides; its original
ReportLab report writer is skipped when the LaTeX source exists. The formatter
uses the saved report content and the existing measured data. If measurements
change, refresh the report content before formatting, so prose and plots agree.

The member contribution section records pending member evidence. Proposed
responsibilities are not represented as actual work. The provenance disclosure
is maintained in `../CONTRIBUTIONS.md`.

## Required structure retained

Target `Group04_report.pdf`: suggested 10–14 pages, maximum 16 excluding appendices.
Performance claims must come from measured data.

1. Introduction: both members, machine/physical and logical cores/RAM, WSL limits,
   compiler/MPI versions, flags, affinity, repetitions and timing methodology.
2. Five case sections: problem, model justification, decomposition figure,
   decisions, correctness, median speedup/efficiency table and graph, answers to
   every critical-analysis question using actual numbers.
3. Cross-case comparison: observed scaling, dominant bottleneck, and best-scaling case.
4. Actual individual contributions and references.

Heat requires weak scaling in addition to strong scaling. Sample sort requires a
phase breakdown. Every graph needs labelled axes, units, and repetitions per point.
Retain raw measurements so every table and claim can be reproduced.
