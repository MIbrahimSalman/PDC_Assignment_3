"""Render measured TSV data into plots and six PDF slides.

Preserve an existing LaTeX report; otherwise create the original ReportLab report.

Run with Python 3 plus matplotlib, numpy and reportlab after benchmark.sh finishes.
Python is used only for analysis/plots/document production; algorithms are C++.
"""
from pathlib import Path
from collections import defaultdict
from statistics import median
import csv
import json
import math
import re
import html
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import reportlab
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.graphics.shapes import Drawing, Rect, String, Line, Polygon

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'
FIGURES = RESULTS / 'figures'
FIGURES.mkdir(exist_ok=True)
(RESULTS/'document_tools.json').write_text(json.dumps({'python':sys.version,'matplotlib':matplotlib.__version__,'numpy':np.__version__,'reportlab':reportlab.Version},indent=2))
rows = defaultdict(list)
verification_text=(ROOT/'verification_all.log').read_text()
final_verification_text=(ROOT/'verification_final.log').read_text()
if 'Final six-worker comparisons and known-value invariants passed.' not in final_verification_text:
    raise RuntimeError('Run tools/check_final.sh successfully before rendering final artifacts')
comparison_count=verification_text.count('verified=yes')+final_verification_text.count('verified=yes')
with (RESULTS / 'raw.tsv').open() as stream:
    for line in stream:
        label, p, repeat, output = line.rstrip().split('\t', 3)
        fields = dict(token.split('=', 1) for token in output.split() if '=' in token)
        fields['repeat'] = int(repeat)
        rows[(label, int(p))].append(fields)
for key, group in rows.items():
    if len(group) != 3 or sorted(r['repeat'] for r in group) != [1, 2, 3]:
        raise RuntimeError(f'{key}: expected exactly three repetitions')

def value(label, p, field='seconds'):
    return median(float(row[field]) for row in rows[(label,p)])

def first(label,p):
    return rows[(label,p)][0]

def fmt(x):
    return f'{x:.4g}'

P = [1,2,4,6]
series = {
    'heat1024_strong': ('heat1024_serial', False),
    'heat4096_strong': ('heat4096_serial', False),
    'kmer8_atomic': ('kmer8_serial',False),
    'kmer8_private': ('kmer8_serial',False),
    'kmer12_atomic': ('kmer12_serial',False),
    'kmer12_private': ('kmer12_serial',False),
    'mc_good': ('mc_serial',False),
    'pagerank_strong': ('pagerank_serial',False),
    'sort_uniform_strong': ('sort_uniform_serial',True),
    'sort_skew_strong': ('sort_skew_serial',True),
}

def baseline(label,p):
    reference, matching = series[label]
    return value(reference,p if matching else 1)

def speedup(label,p):
    return baseline(label,p)/value(label,p)

with (RESULTS/'summary.csv').open('w',newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['experiment','workers','runs','median_seconds','min_seconds','max_seconds','serial_seconds','speedup','efficiency'])
    for (label,p), group in sorted(rows.items()):
        times = [float(r['seconds']) for r in group]
        extra = [baseline(label,p),speedup(label,p),speedup(label,p)/p] if label in series else ['','','']
        writer.writerow([label,p,3,median(times),min(times),max(times),*extra])

plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
palette = ['#176B87','#D97706','#7950A0','#15946D']

def savefig(name):
    plt.tight_layout()
    plt.savefig(FIGURES/f'{name}.png',dpi=180)
    plt.savefig(FIGURES/f'{name}.pdf')
    plt.close()

def scaling(name, curves, title):
    plt.figure(figsize=(6.5,2.8))
    for i,(label,legend) in enumerate(curves):
        points = P+[8] if label == 'mc_good' else P
        speeds = [speedup(label,p) for p in points]
        lower = [s-baseline(label,p)/max(float(r['seconds']) for r in rows[(label,p)]) for s,p in zip(speeds,points)]
        upper = [baseline(label,p)/min(float(r['seconds']) for r in rows[(label,p)])-s for s,p in zip(speeds,points)]
        plt.errorbar(points,speeds,yerr=[lower,upper],marker='o',capsize=2,label=legend,color=palette[i])
    plt.plot(points,points,'--',color='#9CA3AF',label='Ideal S=p')
    plt.xlabel('Workers (8 uses SMT on 6 physical cores)' if name=='mc' else 'Workers (physical cores)')
    plt.ylabel('Speedup vs serial (x)')
    plt.title(title+' | medians of 3 runs; whiskers min-max')
    plt.xticks(points)
    plt.grid(alpha=.2)
    plt.legend(fontsize=7,ncol=2)
    savefig(name)

scaling('heat',[('heat1024_strong','N=1024'),('heat4096_strong','N=4096')],'Heat: 500 steps, 2-D overlap, K=50')
scaling('kmer',[(f'kmer{k}_{mode}',f'k={k}, {mode}') for k in [8,12] for mode in ['atomic','private']],'k-mers: 100 million bases')
scaling('mc',[('mc_good','100 million paths')],'Monte Carlo: fixed sample count')
scaling('pagerank',[('pagerank_strong','Pull, contiguous static')],'PageRank: 1 million edges')
scaling('sort',[('sort_uniform_strong','Uniform'),('sort_skew_strong','90% duplicate zeros')],'Sample sort: 100 million keys')

plt.figure(figsize=(6.5,2.7))
weak0 = value('heat_weak',1)
plt.plot(P,[value('heat_weak',p)/weak0 for p in P],'o-',color=palette[0],label='Measured T(p)/T(1)')
plt.axhline(1,linestyle='--',color='#888888',label='Ideal constant time')
plt.xticks(P); plt.xlabel('MPI ranks'); plt.ylabel('Normalized elapsed time (x)')
plt.title('Heat weak scaling: ~512² interior cells/rank | 3 runs/point')
plt.legend(); plt.grid(alpha=.2); savefig('weak')

degrees = list(csv.DictReader((RESULTS/'degrees.csv').open()))
plt.figure(figsize=(6.5,2.5))
plt.loglog([int(r['degree']) for r in degrees if int(r['degree'])>0],
           [int(r['vertices']) for r in degrees if int(r['degree'])>0],'.',color=palette[0])
plt.xlabel('Incoming degree (edges)'); plt.ylabel('Vertices with that degree (count)')
plt.title('Fixed-seed R-MAT degree distribution; zero degree excluded from log axes')
plt.grid(alpha=.2); savefig('degrees')

plt.figure(figsize=(6.5,2.8))
bottom = np.zeros(len(P))
for i,(field,label) in enumerate([('local_sort_seconds','Local sort'),('sampling_seconds','Sampling/partition'),('exchange_seconds','Counts/allocation/Alltoallv'),('final_sort_seconds','Final sort')]):
    times = np.array([value('sort_uniform_strong',p,field) for p in P])
    plt.bar([str(p) for p in P],times,bottom=bottom,label=label,color=palette[i])
    bottom += times
plt.xlabel('MPI ranks'); plt.ylabel('Phase time (seconds)')
plt.title('Uniform sample sort: component medians, 3 runs; slowest-total rank')
plt.legend(fontsize=7); savefig('phases')

ping = defaultdict(list)
for repetition in [1,2,3]:
    for row in csv.DictReader((RESULTS/f'pingpong_{repetition}.csv').open()):
        ping[int(row['bytes'])].append(float(row['one_way_seconds']))
x = np.array(sorted(k for k in ping if k>=4096),dtype=float)
y = np.array([median(ping[int(k)]) for k in x])
unconstrained_beta,unconstrained_alpha = np.polyfit(x,y,1)
# Zero-byte ping-pong supplies a physical latency estimate. Anchor the bandwidth
# fit there; a free-intercept line can return a meaningless negative latency.
alpha = median(ping[0])
beta = max(0.0,float(np.dot(x,y-alpha)/np.dot(x,x)))
r2 = 1-np.sum((y-(alpha+beta*x))**2)/np.sum((y-y.mean())**2)
(RESULTS/'pingpong_fit.json').write_text(json.dumps({'alpha_seconds':alpha,'alpha_method':'median zero-byte one-way time','beta_seconds_per_byte':beta,'r_squared':r2,'fit_min_bytes':4096,'fit_max_bytes':4194304,'runs_per_size':3,'unconstrained_alpha_seconds':unconstrained_alpha,'unconstrained_beta_seconds_per_byte':unconstrained_beta},indent=2))

# Every report page is explicit. Tables/paragraphs use real data and preserve caveats.
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(name='BodySmall',fontName='Helvetica',fontSize=9.3,leading=12.4,spaceAfter=7))
styles.add(ParagraphStyle(name='CaptionSmall',fontName='Helvetica',fontSize=8,leading=10,spaceAfter=6,textColor=colors.HexColor('#475569')))
styles['Heading1'].fontSize=18; styles['Heading1'].leading=22
styles['Heading2'].fontSize=12; styles['Heading2'].leading=15
story=[]
def para(text): story.append(Paragraph(text,styles['BodySmall']))
def heading(text): story.append(Paragraph(text,styles['Heading2']))
def page(title):
    if story: story.append(PageBreak())
    story.append(Paragraph(title,styles['Heading1']))
def table(data,widths=None):
    converted=[[Paragraph(html.escape(str(cell)),styles['CaptionSmall']) for cell in row] for row in data]
    t=Table(converted,colWidths=widths,repeatRows=1,hAlign='LEFT')
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E0EEF4')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#176B87')),('BOTTOMPADDING',(0,0),(-1,-1),3),('TOPPADDING',(0,0),(-1,-1),3)]))
    story.append(t);story.append(Spacer(1,7))
def chart(name,height=195): story.append(Image(str(FIGURES/f'{name}.png'),width=470,height=height))
def scale_table(label):
    table([['p','Serial (s)','Parallel (s)','S(p)','E(p)']]+[[p,fmt(baseline(label,p)),fmt(value(label,p)),fmt(speedup(label,p)),f'{speedup(label,p)/p:.1%}'] for p in P],[35,105,110,90,100])
def diagram(text):
    if text.startswith('1-D'):
        labels=['Rank-owned blocks','Neighbor halos','Interior / border update']
    elif text.startswith('DNA'):
        labels=['Thread window ranges','Private histograms','Merge by bin']
    elif text.startswith('Global'):
        labels=['Disjoint path counters','Worker block moments','Ordered reduction']
    elif text.startswith('Pull'):
        labels=['Worker-owned vertices','CSR neighbor traversal','Ranks / L1 reduction']
    else:
        labels=['Rank-local sorted keys','Splitters / Alltoallv','Owned output buckets']
    drawing=Drawing(470,48)
    for i,label in enumerate(labels):
        x=i*160
        drawing.add(Rect(x,8,145,32,fillColor=colors.HexColor('#E0EEF4'),strokeColor=colors.HexColor('#176B87'),rx=4,ry=4))
        drawing.add(String(x+72.5,21,label,fontName='Helvetica',fontSize=8,textAnchor='middle'))
        if i<2:
            drawing.add(Line(x+147,24,x+158,24,strokeColor=colors.HexColor('#D97706')))
            drawing.add(Polygon([x+158,24,x+154,27,x+154,21],fillColor=colors.HexColor('#D97706'),strokeColor=None))
    story.append(drawing)
    story.append(Paragraph(html.escape(text).replace('\n','<br/>'),styles['CaptionSmall']))

page('Group 04 | Parallel case studies')
para('<b>CS3006 - Assignment 3</b><br/>Fatima Malik (23L-0826) and Ibrahim Salman (23L-0723)<br/>Measurements and analysis, 23 November 2026.')
heading('Project scope')
para('Five C++17 serial/parallel implementations, correctness tests, three-repeat timings, plots and a six-slide presentation are included. Contribution status and the consolidated provenance disclosure are maintained in CONTRIBUTIONS.md.')
heading('Machine and toolchain')
para('AMD Ryzen 5 5500U: 6 physical cores, 12 logical processors. Host RAM: 7,866,085,376 bytes (7.33 GiB); WSL exposes approximately 3.5 GiB plus 1 GiB swap. Ubuntu WSL2, Linux 6.6.87.2; GCC 13.3.0 and Open MPI 4.1.6. All executables use -O3 -std=c++17 -Wall -Wextra -Wpedantic; OpenMP adds -fopenmp. This is a single-machine shared-memory MPI experiment, not a multi-node network experiment.')
heading('Measurement method')
para('Every reported point is the median of three independent process invocations. The raw TSV, source SHA-256 manifest, build log and environment dump are in results/. Generation and verification are outside timing except Monte Carlo, where generating random paths is part of the computation. MPI uses barriers around the timed region and reports maximum elapsed time across ranks. OpenMP uses OMP_PROC_BIND=close, OMP_PLACES=cores, OMP_DYNAMIC=false. MPI uses --bind-to core --map-by core. Workers are 1,2,4,6; only Monte Carlo additionally uses 8 workers (SMT). No MPI oversubscription is used.')
para('S(p)=T_serial/T_parallel and E(p)=S(p)/p. The serial baseline is a separate executable, compiled with the same optimization level. Sort baselines recreate the same p-seeded shards. Runs execute sequentially; the agent cannot guarantee the Windows host was idle. WSL scheduling, background host activity, frequency changes and thermal state can affect results. Whiskers show the observed range, not confidence intervals. Interpret these as measurements of this session, and rerun on an idle machine for assessment.')
table([['Case','Model / decomposition','Final input'],['Heat','MPI / 1-D rows and 2-D blocks','1024² and 4096²; 500 steps'],['k-mers','OpenMP / window-start intervals','100 million bases; k=8,12'],['Monte Carlo','OpenMP / fixed sample blocks','100 million paths'],['PageRank','OpenMP / destination/source vertices','65,536 vertices; 1 million edges'],['Sample sort','MPI / key buckets','100 million 32-bit keys']],[80,210,180])
para(f'All 244 small-input comparisons passed before measurement: 49 heat, 60 k-mer, 9 Monte Carlo, 90 PageRank, 36 sort. Final six-worker tests bring the total to {comparison_count} comparisons, plus three known-value invariants. Deliberately incorrect variants are excluded. These tests do not prove correctness for every possible input.')

page('1 | Heat diffusion: implementation and scaling')
para('A two-grid explicit five-point stencil updates independent next-step cells from the old grid. Alpha=0.2; the top boundary is fixed at 100, other boundaries at zero. MPI fits spatial ownership and nearest-neighbor halos. Only interior cells are partitioned; uneven blocks are supported.')
diagram('1-D: [row block 0] <halo rows> [row block 1] <halo rows> ...\n2-D: Cartesian block grid; each block exchanges north/south rows and east/west columns. Nonblocking: post halos -> compute strict interior -> wait -> compute border -> swap grids.')
para('Blocking uses MPI_Sendrecv; overlap uses Irecv/Isend and a vector datatype for columns. Every K steps and on the final step, Allreduce computes maximum change. There is no early stopping. The serial oracle checks every owned cell at absolute tolerance 1e-12, outside timing. The oracle is debug-only and replicated per rank.')
chart('heat')
scale_table('heat4096_strong')
para('Table: N=4096, 500 steps, 2-D overlap, K=50. N=1024 measurements and both input baselines are retained in summary.csv. Phase maxima are separately reduced and need not sum to maximum elapsed time.')

page('1 | Heat diffusion: critical analysis')
heading('Surface-to-volume and overlap')
para('Let L=N-2. Work is approximately L²/p cells/rank/step. Away from external boundaries, 1-D sends 2L halo values; a square 2-D layout sends 4L/sqrt(p). Receives add the same traffic again. Thus communication/work grows as 2p/L for 1-D and 4sqrt(p)/L for 2-D. For a rectangular a-by-b process grid, the interior send count is 2L/a+2L/b; physical-edge ranks communicate less. At p=1 there are no real neighbor transfers.')
heat_rows=[['Layout / exchange','K','Median (s)','Exposed exchange (s)']]
for d in ['1d','2d']:
    for mode in ['blocking','overlap']:
        for k in [1,50]:
            label=f'heat_compare_{d}_{mode}_k{k}'
            heat_rows.append([f'{d} / {mode}',k,fmt(value(label,4)),fmt(value(label,4,'exchange_seconds'))])
table(heat_rows,[150,35,115,170])
b=value('heat_compare_2d_blocking_k50',4); o=value('heat_compare_2d_overlap_k50',4)
para(f'At N=1024, p=4, K=50, overlap changes total time by {(b-o)/b:.1%} relative to blocking ({b:.4g} to {o:.4g} s). This is the observed net effect, not a direct measure of transfer time hidden. Posting/wait timing excludes background transfer; shared-memory transport, progress, interior size and cache effects confound attribution. The formula predicts a lower asymptotic surface ratio for 2-D, but it does not guarantee a faster small single-node run.')
low=next((p for p in P if speedup('heat4096_strong',p)/p<.7),None)
para(f'For N=4096, efficiency first falls below 70% at p={low if low else "none of the measured points"}. At six ranks the effective overhead T(p)-Tserial/p is {value("heat4096_strong",6)-baseline("heat4096_strong",6)/6:.4g} s. This residual includes memory-bandwidth contention and changed compute throughput, not just communication. Use Tcomp/p+Tcomm as an explanatory approximation, not a fitted network-only law. K=1 versus K=50 timings above expose the combined cost of more per-cell change work and more collectives; the change metric work is also conditional on K.')
chart('weak',165)
para('Weak scaling uses N=2+round(512 sqrt(p)), 500 steps and K=50. Interior cells per rank stay near 262,144; rounding and rectangular process layouts change the exact surface. Ideal normalized time is one. The corresponding weak efficiency is T(1)/T(p), not serial strong-scaling speedup.')

page('2 | k-mer counting: implementation and scaling')
para('Each DNA base occupies one byte and each k-mer uses a rolling two-bit index. Threads own disjoint ranges of window start positions, read k-1 bases past the range if necessary, and count boundary windows once. The serial baseline has O(L) rolling-window work.')
diagram('DNA: A C G T ... -> disjoint window-start ranges -> worker histograms\nA: shared increments (intentional race). B: shared atomic increments. C: private histograms -> merge by bin.')
para('Fixed-seed SplitMix64 input is identical for every worker count. Integer verification compares every histogram bin and checks total L-k+1. Allocation/zeroing of the shared histogram and input generation are outside timing; private histogram allocation, initialization and merge are included. Verification is excluded. Counters are 64-bit.')
chart('kmer')
scale_table('kmer8_private')
para('Table: k=8 private histograms. All four k/mode curves use L=100,000,000 and seed 42. Results for k=12 include substantially larger private allocation and merge costs; compare the measured totals rather than assuming the private method always wins.')

page('2 | k-mer counting: critical analysis')
race=first('race',4)
para(f'<b>Race:</b> in the first observed skewed-input demonstration, total={race["total"]}, expected={race["expected"]}, and {race["mismatched_bins"]} bins differed. Two threads can read the same old count and each store old+1, losing an update. This is undefined behavior in C++; apparent agreement on another run would not make Version A correct.')
table([['k','p','Atomic (s)','Private (s)','Private memory']]+[[k,p,fmt(value(f'kmer{k}_atomic',p)),fmt(value(f'kmer{k}_private',p)),f'{(4**k)*8*p/2**20:.0f} MiB'] for k in [8,12] for p in P],[30,30,115,115,180])
para('Private memory is 4^k * 8 * p bytes: 512 KiB/worker at k=8 and 128 MiB/worker at k=12. The six-core CPU has 512 KiB L2/core and a shared 4 MiB L3. k=8 can exploit locality but k=12 greatly exceeds cache, with extra zeroing and merge traffic. Atomic updates pay read-modify-write/coherence costs; private bins avoid contention but multiply memory and merge work. The table gives the actual crossover, including overhead rather than only the counting loop.')
table([['p','Packed (s)','64-byte padded (s)','Packed / padded']]+[[p,fmt(value('false_sharing_packed',p)),fmt(value('false_sharing_padded',p)),fmt(value('false_sharing_packed',p)/value('false_sharing_padded',p))] for p in P],[35,130,155,150])
para('The counter microbenchmark performs 100 million volatile increments per worker. Separate counters share cache lines when packed; alignas(64) separates them under this CPU cache-line assumption. Each worker owns a distinct object, so this is false sharing rather than a data race. Padding ratios above are measured, not a guaranteed speedup; the short microbenchmark is sensitive to scheduling and CPU state.')
para(f'Rolling scans read about one new sequence byte/window after each initial k-byte window. At p=6, k=8 private reaches {value("kmer8_private",6,"million_kmers_per_second"):.3g} million windows/s; k=12 private reaches {value("kmer12_private",6,"million_kmers_per_second"):.3g}. Histogram loads/stores, cache-line transfers, atomic coherence and merges dominate additional traffic. These input-read rates are not measured DRAM bandwidth. The finite cache and added private memory explain why additional workers can stop helping; hardware counters would be needed to separate bandwidth saturation from scheduling effects.')

page('3 | Monte Carlo: implementation and scaling')
para('European call pricing uses S0=100, strike=100, r=0.05, sigma=0.2 and T=1 year. Payoff is exp(-rT)*max(S_T-K,0), with S_T=100*exp(0.03+0.2Z). Box-Muller produces normal samples. Standard error is sqrt(sample variance/M), using M-1 degrees of freedom.')
diagram('Global path indices -> disjoint counter ranges -> fixed 65,536-path blocks\nWorkers compute (payoff sum, squared-payoff sum) -> ordered block reduction -> price and standard error.')
para('OpenMP is appropriate because paths are independent and communicate only final partial moments. Two unique counters per path are mapped through a SplitMix64 permutation, avoiding shared RNG state and duplicate stream assignments. This is a pseudorandom method, not a claim of mathematically independent random variables. Ordered block reduction in long double makes numerical output independent of scheduling and p on this toolchain. Both price and standard error are checked against serial with absolute tolerance 1e-12.')
chart('mc')
scale_table('mc_good')
mc=first('mc_good',1)
para(f'For 100 million paths, price={float(mc["price"]):.8f}, standard error={float(mc["standard_error"]):.6g}, analytic Black-Scholes={float(mc["black_scholes"]):.8f}, standardized deviation={float(mc["z_score"]):.3f}. Statistical error is distinct from serial-equivalence error. Random path generation is included in computation timing.')

page('3 | Monte Carlo: critical analysis')
heading('Why random stream ownership matters')
good=first('mc_good',8); bad=first('mc_same_stream',8)
para(f'At p=8, good-stream price={float(good["price"]):.8f}; intentional repeated-stream price={float(bad["price"]):.8f}. The bad version reports naive SE={float(bad["standard_error"]):.6g}, but repeating 12.5 million unique paths eight times does not yield 100 million independent observations. An approximate corrected SE is sqrt(8) times larger: {float(bad["standard_error"])*math.sqrt(8):.6g}. Reusing the same worker seed corrupts the uncertainty estimate even when the price looks plausible. A locked shared generator would instead add serial contention; this implementation avoids both designs.')
heading('Accuracy per second and the fixed-time experiment')
table([['Workers','Paths','Actual time (s)','Standard error']]+[[p,first('mc_fixed_time12',p)['samples'],fmt(value('mc_fixed_time12',p)),fmt(value('mc_fixed_time12',p,'standard_error'))] for p in [1,8]],[65,145,130,130])
fixed_ratio=value('mc_fixed_time12',1,'standard_error')/value('mc_fixed_time12',8,'standard_error')
para(f'Each workload was calibrated from measured throughput to a twelve-second target with at least 100 million paths, then run three times; actual times above disclose the mismatch. This is approximately fixed time, not a deadline-controlled sampler. The observed SE reduction factor is {fixed_ratio:.3g}. Since SE is proportional to 1/sqrt(M), a throughput gain g yields about sqrt(g) more precision at fixed time, or 1/g of the time for a fixed target variance. Eight workers use SMT on six physical cores, so their gain cannot be read as eight-core scaling. Earlier three-second exploratory runs are retained separately in the raw data.')
heading('Amdahl estimate and limitations')
s6=speedup('mc_good',6); fraction=(1/s6-1/6)/(1-1/6)
para(f'At p=6, speedup={s6:.4g}; the effective Amdahl fraction f=(1/S-1/p)/(1-1/p) is {fraction:.4g}. This is an effective fit from one point, not a measured proportion of serial instructions. Reduction, startup, shared execution resources, CPU turbo/thermal behavior, virtualization and host load can all enter the fitted fraction. A negative or >1 fit would signal a poor model rather than a physical fraction. This workload has little communication, but its observed scaling must be used rather than the assumption that it scales perfectly.')
heading('Reproducibility')
prices=[float(first('mc_good',p)['price']) for p in P+[8]]
para(f'The spread in recorded good-stream prices across p=1,2,4,6,8 is {max(prices)-min(prices):.3g}. Fixed path indices select the same random values, each block keeps the same internal summation order, and the final reduction follows block order. A usual worker-seeded implementation changes sample assignments with p; an OpenMP floating-point reduction may also change the reduction tree. Because addition is not associative, such implementations can change the final low-order bits even for identical samples. Reproducibility here is within this compiler/platform, not a promise of bitwise cross-platform math-library behavior.')

page('4 | PageRank: implementation and scaling')
para('A fixed-seed R-MAT directed multigraph has 65,536 vertices and one million edges. Quadrant probabilities are (0.57,0.19,0.19,0.05); self-loops and duplicate edges remain. Incoming and outgoing CSR avoid per-edge object overhead. Damping is 0.85, L1 convergence tolerance 1e-10, and the safety cap is 500 iterations.')
diagram('Pull: destination-owned CSR rows <- read in-neighbor contributions (no writes shared)\nPush: source-owned CSR rows -> atomic contributions to destination ranks\nAll workers reduce dangling mass and L1 change each iteration.')
para('Dangling vertices redistribute their mass uniformly. Each next rank begins with (1-d)/n+d*dangling/n. Serial is pull; parallel pull/push verify maximum error <=1e-9, L1 error <=1e-8, mass within 1e-8 of one, and convergence. Tests include a graph with zero edges. Graph/CSR creation is excluded; contribution preparation and all iteration overhead are included.')
chart('pagerank')
scale_table('pagerank_strong')
pr=first('pagerank_strong',4)
para(f'The four-worker static pull run converged in {pr["iterations"]} iterations with rank sum {float(pr["rank_sum"]):.12f}. Static work-loop imbalance max/mean is {value("pagerank_strong",4,"imbalance"):.3g}. Timers stop before the final loop barrier; including barrier wait would hide the imbalance.')

page('4 | PageRank: critical analysis')
chart('degrees',160)
para('Equal vertex counts do not imply equal work: high-degree destinations require longer CSR scans. The degree plot shows the actual graph rather than assuming a perfect power law. Incoming degree zero is excluded from log axes, but the solver keeps all vertices. The imbalance ratio measures accumulated edge-loop time, including scheduler overhead and excluding final wait.')
schedules=[(m,s,c) for m in ['pull','push'] for s in ['static','dynamic','guided'] for c in [1,64,256]]
table([['Schedule / chunk','Pull seconds','Push seconds']]+[[f'{s} / {c}',fmt(value(f'pagerank_pull_{s}_{c}',4)),fmt(value(f'pagerank_push_{s}_{c}',4))] for s in ['static','dynamic','guided'] for c in [1,64,256]],[170,150,150])
best=min(schedules,key=lambda entry:value(f'pagerank_{entry[0]}_{entry[1]}_{entry[2]}',4))
bestlabel=f'pagerank_{best[0]}_{best[1]}_{best[2]}'
push=value('pagerank_push_dynamic_64',4); pull=value('pagerank_pull_dynamic_64',4)
para(f'The fastest tested four-worker configuration is {best[0]}, {best[1]} chunk {best[2]} at {value(bestlabel,4):.4g} s. Small dynamic chunks can improve work balance but require more scheduling operations; larger chunks amortize scheduling while leaving some imbalance. Static cyclic chunks redistribute concentrated degrees, while contiguous static creates a different ownership pattern. Results, not labels, decide the best setting.')
para(f'At dynamic chunk 64, push/pull time ratio is {push/pull:.3g} ({push:.4g}/{pull:.4g} s). Push includes random atomic destination writes and cache-line contention; pull uses random source reads and exclusive output writes. This ratio is a combined traversal/atomic cost comparison, not an isolated atomic latency measurement.')
rate=value(bestlabel,4,'million_edges_per_second')
para(f'Best configuration throughput is {rate:.3g} million edges/s. A simple pull model counts 4-byte neighbor IDs plus 8-byte contribution reads per edge before reuse, giving an effective {rate*12/1000:.3g} GB/s, plus O(vertices) traffic. Contributions fit much better in cache than the entire edge structure, so effective traffic is not measured DRAM traffic. No independent memory-bandwidth ceiling was measured: proximity to that ceiling cannot be quantified honestly. The access pattern, schedule overhead and contention identify likely costs, but a STREAM/hardware-counter run would be needed for a roofline claim.')

page('5 | Sample sort: implementation and scaling')
para('Each MPI rank owns an independently seeded shard of the 100-million-key input. Uniform keys are 32-bit values; skew sets approximately 90% of keys to zero. Local sorting is followed by p-1 regular samples/rank, global quantile splitters, bucket partitioning, Alltoall counts, Alltoallv keys, then a final local sort.')
diagram('Rank-local keys -> local sort -> regular samples -> root selects splitters -> broadcast\nSorted local buckets -> Alltoallv redistribution -> local sort -> distributed global order')
para('Production execution never gathers full input. All runs check local order, consecutive nonempty rank boundaries, element count and two multiset fingerprints. Debug inputs <=1 million keys also compare exactly with a gathered serial result. Fingerprints are probabilistic, not proof of exact equality. Serial generation recreates the same p shards; its time excludes generation just as parallel time does.')
chart('sort')
scale_table('sort_uniform_strong')
para('Table: uniform keys, oversampling factor one. Each serial baseline is matched to that p-seeded input. Algorithms and size stay fixed, but exact input values vary with p because the assignment asks for per-rank seeds. This is a limitation of strict fixed-dataset strong-scaling interpretation.')

page('5 | Sample sort: critical analysis')
chart('phases',180)
table([['Distribution / samples','p=4 max/mean','Final bucket counts']]+[[f'{dist} / factor {factor}',fmt(value(f'sort_{dist}_'+('strong' if factor==1 else 'oversample8'),4,'imbalance')),first(f'sort_{dist}_'+('strong' if factor==1 else 'oversample8'),4)['bucket_sizes']] for dist in ['uniform','skew'] for factor in [1,8]],[140,100,230])
para('Oversampling changes quantile resolution, but equal splitters still send identical keys to the same bucket. Thus the 90%-zero input overloads one rank and can leave intermediate ranks empty. More samples cannot fix this partition rule; a secondary global index or explicit equal-key distribution would be required. Correctness checks deliberately skip empty ranks when testing global boundaries.')
phases=['local_sort_seconds','sampling_seconds','exchange_seconds','final_sort_seconds']
dominant={p:max(phases,key=lambda f:value('sort_uniform_strong',p,f)) for p in [1,6]}
para(f'For uniform data the largest measured phase is {dominant[1].replace("_seconds", "").replace("_", " ")} at p=1 and {dominant[6].replace("_seconds", "").replace("_", " ")} at p=6. Stacks use medians of phases from the slowest-total rank of each run; the sum of component medians need not equal median elapsed time. Ending-barrier wait is excluded from the stack. Exchange includes count transfer and receive allocation as well as Alltoallv.')
prediction=3*alpha+beta*4*100000000*3/16
observed=value('sort_uniform_strong',4,'max_exchange_seconds')
para(f'Ping-pong uses the median zero-byte one-way time for latency alpha={alpha:.4g} s. Anchoring there, a least-squares fit of median times from 4 KiB through 4 MiB gives beta={beta:.4g} s/byte and R²={r2:.4f} (three runs/size). A balanced p=4 model, (p-1)alpha+beta*4N(p-1)/p², predicts {prediction:.4g} s; measured maximum exchange phase is {observed:.4g} s. This extrapolates beyond calibration messages and compares ideal transfer cost to a wider phase. Shared-memory paths, collective scheduling, bucket imbalance and allocation limit agreement. The unconstrained fitted intercept was {unconstrained_alpha:.3g} s; its nonphysical value shows why a single straight-line model across transport regimes needs caution.')
para('Sub-linear speedup is expected: local sorting reduces O(N log N) work per rank, but redistribution still moves O(N) keys, sample selection serializes metadata work, and final local sorting repeats work on received keys. Cache behavior and physical memory bandwidth matter even with perfect buckets. Duplicate-heavy imbalance adds a distinct bottleneck. The implementation supports N<=INT_MAX because legacy Alltoallv counts/displacements are ints.')

page('6 | Cross-case comparison and reproducibility')
choices=[('Heat','heat4096_strong','MPI','Bandwidth / halo synchronization'),('k-mers','kmer8_private','OpenMP','Histogram locality / merge traffic'),('Monte Carlo','mc_good','OpenMP','Compute / runtime overhead'),('PageRank','pagerank_strong','OpenMP','Irregularity / memory / scheduling'),('Sample sort','sort_uniform_strong','MPI','Sorting / data redistribution')]
table([['Case','Model','S(6)','E(6)','Primary costs']]+[[name,model,fmt(speedup(label,6)),f'{speedup(label,6)/6:.1%}',cost] for name,label,model,cost in choices],[80,55,55,60,220])
winner=max(choices,key=lambda row:speedup(row[1],6))
para(f'Among the representative six-worker configurations above, {winner[0]} has the highest measured speedup ({speedup(winner[1],6):.3g}x). This compares the listed algorithms and parameters, not universal properties of the domains. Absolute elapsed times, different work sizes and overheads also matter. Poor scaling is an outcome to explain, not a correctness failure.')
heading('How to reproduce')
para('In Ubuntu WSL, enter Group04_Assignment3 and run make check. Then run bash tools/benchmark.sh on an otherwise idle machine. The script runs sequentially, resumes completed labelled repetitions, and rejects source changes against its SHA-256 manifest. Archive existing results before measuring a different source version or machine. Rebuild before timing; preserve compiler and affinity settings.')
para('Run python tools/render_results.py with matplotlib, numpy and reportlab installed to regenerate summary.csv, labelled figures, this report and exactly six PDF slides. The renderer requires three repetitions per group and all planned experiments. Raw results/raw.tsv is the source of every performance claim. The sample-sort serial baseline is indexed by both distribution and shard count.')
heading('Verification and limitations')
para('The 244 pre-benchmark checks cover uneven heat blocks, one-cell-wide domains, k-mer chunk edges, k=1/8/12, zero-edge PageRank, empty sort ranks, skew, and multiple workers. Invalid argument handling rejects unsupported sizes. Full-scale timing runs disable expensive serial oracles, but sort always checks global order/count/fingerprints. Version A k-mer races and repeated-stream Monte Carlo are explicitly incorrect demonstrations, never counted as valid speedups.')
para('Measurements were made in WSL on a laptop, not a dedicated cluster. CPU frequency, virtualization, host memory pressure and scheduler interference were not isolated. Only three samples are available at each point. Min-max whiskers show sensitivity; no narrow statistical confidence claim is justified. Profiling hardware counters, memory bandwidth and multi-node MPI would distinguish additional causes. Bonus extensions were not implemented.')
heading('Artifact map')
table([['Artifact','Location'],['Buildable implementations','cs1_heat_mpi/ through cs5_samplesort_mpi/'],['Correctness scripts/log','tools/check*.sh; verification_all.log'],['Raw data / summary / figures','results/raw.tsv; summary.csv; figures/'],['Report and six-slide deck','report/Group04_report.pdf; slides/Group04_slides.pdf'],['Contribution provenance','CONTRIBUTIONS.md; actual member evidence pending']],[190,280])

page('7 | Contributions and references')
heading('Individual contributions: awaiting member confirmation')
para('Actual member hours and authored changes remain unconfirmed. Proposed responsibilities below are not a record of completed work. The consolidated contribution and provenance record is maintained in CONTRIBUTIONS.md. Repository: <link href="https://github.com/fatimamalik2307/pdc_assignment3-" color="#176B87">github.com/fatimamalik2307/pdc_assignment3-</link>.')
table([['Member','Actual hours / authored work','Proposed lead / reviewer'],['Fatima Malik\n23L-0826','Not supplied; member confirmation required','Lead 1,3,5; review 2,4'],['Ibrahim Salman\n23L-0723','Not supplied; member confirmation required','Lead 2,4; review 1,3,5']],[120,190,160])
para('The proposed split gives each member at least two leads and two reviews, with neither leading more than three. Final integration, comparison and slides are shared responsibilities. Actual reviews, authored changes, hours and commit records should be recorded as completed. Peer evaluations and presentation rehearsal remain member responsibilities.')
heading('References')
para('[1] Supplied course document: <i>Assignment 3: Parallel Case Studies with OpenMP and MPI</i>, CS3006, PDC_assignment_3.pdf. Requirements, case definitions and integrity policy; especially Sections 3-7 and 11.')
para('[2] OpenMP Architecture Review Board. <i>OpenMP API Specification, Version 5.2</i>. Worksharing, schedules, atomics and reductions. <link href="https://www.openmp.org/spec-html/5.2/openmp.html" color="#176B87">openmp.org/spec-html/5.2/openmp.html</link>. Accessed 23 November 2026.')
para('[3] Open MPI. <i>Version 4.1 documentation</i>. Cartesian topology, point-to-point requests and collectives. <link href="https://www.open-mpi.org/doc/v4.1/" color="#176B87">open-mpi.org/doc/v4.1/</link>. Installed runtime: 4.1.6. Accessed 23 November 2026.')
para('[4] D. Chakrabarti, Y. Zhan and C. Faloutsos. <i>R-MAT: A Recursive Model for Graph Mining</i>. SIAM International Conference on Data Mining, 2004, pp. 442-446. <link href="https://doi.org/10.1137/1.9781611972740.43" color="#176B87">doi:10.1137/1.9781611972740.43</link>. Synthetic graph method; no external dataset was downloaded.')
para('[5] Reproducibility evidence generated with this package: results/raw.tsv, environment.txt, host_machine.json, sources.sha256, pingpong_1.csv through pingpong_3.csv, and verification_all.log. These are the primary sources for all numerical claims in the report.')
heading('Project records')
para('The package includes code, measured graphs, this report, six slides and a ZIP. CONTRIBUTIONS.md is the source of the provenance disclosure and contribution status; GIT_COMMIT_PLAN.md records the proposed import sequence. Performance claims are supported by the retained measurement data.')

def footer(c,doc):
    c.setFont('Helvetica',8); c.setFillColor(colors.HexColor('#64748B'))
    c.drawString(42,24,'Group 04 | Parallel case studies | 23 November 2026')
    c.drawRightString(553,24,str(doc.page))

report=ROOT/'report'/'Group04_report.pdf'
if not report.with_suffix('.tex').exists():
    SimpleDocTemplate(str(report),pagesize=(595,842),rightMargin=45,leftMargin=45,topMargin=38,bottomMargin=42,title='Group 04 - Parallel Case Studies',author='Group 04 project report').build(story,onFirstPage=footer,onLaterPages=footer)

# Two presentation editions share the same case-study content.
slides=None
def slide_title(title,subtitle,index):
    slides.setFillColor(colors.HexColor('#102C3B')); slides.rect(0,0,960,540,fill=1,stroke=0)
    slides.setFillColor(colors.white);slides.setFont('Helvetica-Bold',27);slides.drawString(36,492,title)
    slides.setFont('Helvetica',13);slides.setFillColor(colors.HexColor('#BDD6E1'));slides.drawString(36,465,subtitle)
    slides.setFont('Helvetica',9);slides.drawString(36,18,'Group 04 | Parallel case studies | Contribution and provenance record: CONTRIBUTIONS.md')
    slides.drawRightString(924,18,str(index))
def slide_text(text,x,y,width=330,size=15):
    style=ParagraphStyle('slide',fontName='Helvetica',fontSize=size,leading=size*1.35,textColor=colors.white)
    p=Paragraph(text,style);w,h=p.wrap(width,400);p.drawOn(slides,x,y-h);return y-h-16
def boxes(labels,y=380):
    width=(888-20*(len(labels)-1))/len(labels)
    for i,label in enumerate(labels):
        x=36+i*(width+20)
        slides.setFillColor(colors.HexColor('#176B87'));slides.roundRect(x,y,width,48,6,fill=1,stroke=0)
        slides.setFillColor(colors.white);slides.setFont('Helvetica-Bold',12);slides.drawCentredString(x+width/2,y+19,label)
        if i<len(labels)-1:
            slides.setStrokeColor(colors.HexColor('#E9B44C'));slides.line(x+width+2,y+24,x+width+17,y+24)
            slides.line(x+width+17,y+24,x+width+12,y+28);slides.line(x+width+17,y+24,x+width+12,y+20)

def member_card(x,y,width,height,name,roll,scope,large=False):
    slides.setFillColor(colors.HexColor('#183F50'))
    slides.roundRect(x,y,width,height,12,fill=1,stroke=0)
    slides.setFillColor(colors.HexColor('#55A5B8'))
    slides.roundRect(x,y+height-5,width,5,2,fill=1,stroke=0)
    slides.setFillColor(colors.white)
    slides.setFont('Helvetica-Bold',28 if large else 24)
    slides.drawString(x+24,y+height-46,name)
    slides.setFont('Helvetica',20 if large else 18)
    slides.setFillColor(colors.HexColor('#BDD6E1'))
    slides.drawString(x+24,y+height-77,roll)
    slide_text('<b>Proposed responsibility</b><br/>'+scope,x+24,y+height-103,width-48,15)

def model_strip(y):
    items=[('Heat diffusion','MPI'),('k-mer counting','OpenMP'),('Monte Carlo','OpenMP'),('PageRank','OpenMP'),('Sample sort','MPI')]
    for i,(case,model) in enumerate(items):
        x=36+i*180
        slides.setFillColor(colors.HexColor('#183F50'))
        slides.roundRect(x,y,168,66,6,fill=1,stroke=0)
        slides.setFillColor(colors.HexColor('#BDD6E1'));slides.setFont('Helvetica',13)
        slides.drawString(x+12,y+43,case)
        slides.setFillColor(colors.white);slides.setFont('Helvetica-Bold',17)
        slides.drawString(x+12,y+17,model)

def introduction(compact):
    slide_title('Parallel case studies','CS3006 - Assignment 3  |  Group 04',1)
    if compact:
        member_card(36,245,432,185,'Fatima Malik','23L-0826','Heat diffusion, Monte Carlo, sample sort')
        member_card(492,245,432,185,'Ibrahim Salman','23L-0723','k-mer counting, PageRank')
        slides.setFillColor(colors.HexColor('#BDD6E1'));slides.setFont('Helvetica-Bold',12)
        slides.drawString(36,214,'EXPERIMENT MACHINE')
        slides.setFillColor(colors.white);slides.setFont('Helvetica',16)
        slides.drawString(36,189,'Ryzen 5 5500U  |  6 cores / 12 threads  |  7.33 GiB RAM  |  Ubuntu WSL2')
        model_strip(87)
        slide_text('Responsibilities are proposed; actual contributions remain to be confirmed.',36,65,888,11)
    else:
        member_card(36,195,432,205,'Fatima Malik','23L-0826','Heat diffusion, Monte Carlo, sample sort',True)
        member_card(492,195,432,205,'Ibrahim Salman','23L-0723','k-mer counting, PageRank',True)
        slide_text('Five workloads. Two programming models. Measured results.',36,133,888,22)
        slide_text('Responsibilities are proposed; actual contributions remain to be confirmed.',36,86,888,12)
    slides.showPage()

def overview():
    slide_title('Experiment setup','One machine, shared inputs and separate serial baselines',2)
    slide_text('<b>Machine and toolchain</b><br/><br/>AMD Ryzen 5 5500U<br/>6 physical cores / 12 logical threads<br/>7.33 GiB host RAM; ~3.5 GiB in WSL<br/>GCC 13.3.0 and Open MPI 4.1.6',36,419,425,19)
    slide_text(f'<b>How we measured</b><br/><br/>Median of three runs per point<br/>Speedup against the serial program<br/>{comparison_count} comparisons + 3 invariants passed<br/>Eight-worker Monte Carlo uses SMT',515,419,409,19)
    model_strip(87)
    slides.showPage()

content=[
 ('Heat diffusion','MPI spatial blocks update a temperature field from neighboring cells.','heat',['Owned cells','Exchange halos','Interior + border'],f'<b>Measured:</b> N=4096, six-rank speedup {speedup("heat4096_strong",6):.2f}x.<br/><br/><b>Insight:</b> surface/volume predicts communication growth; memory throughput also limits this single-node run.<br/><br/><b>Next:</b> profile memory bandwidth and MPI progress separately.'),
 ('k-mer counting','OpenMP counts rolling DNA windows in a histogram.','kmer',['Window starts','Private bins','Merge bins'],f'<b>Measured:</b> k=8 private, six-thread speedup {speedup("kmer8_private",6):.2f}x.<br/><br/><b>Insight:</b> k=12 needs 128 MiB private histogram per worker; six copies use 768 MiB.<br/><br/><b>Next:</b> test sparse/tiled histograms and NUMA-aware initialization.'),
 ('Monte Carlo pricing','OpenMP estimates a European call from independent path counters.','mc',['Path counters','Block moments','Ordered reduction'],f'<b>Measured:</b> six-thread speedup {speedup("mc_good",6):.2f}x; SE {float(mc["standard_error"]):.4g}.<br/><br/><b>Insight:</b> fixed blocks preserve results across worker counts; accuracy improves as sqrt(samples).<br/><br/><b>Next:</b> compare antithetic variates at equal runtime.'),
 ('PageRank','OpenMP pull/push traversals rank a synthetic directed graph.','pagerank',['CSR ownership','Read / atomic add','Converge'],f'<b>Measured:</b> static p=4 imbalance {value("pagerank_strong",4,"imbalance"):.2f}; best tested p=4 time {value(bestlabel,4):.3g} s.<br/><br/><b>Insight:</b> equal vertices do not mean equal edge work.<br/><br/><b>Next:</b> partition by cumulative edge count, then compare schedules.'),
 ('Distributed sample sort','MPI redistributes sorted key buckets for global ordering.','sort',['Local sort','Splitters / Alltoallv','Final sort'],f'<b>Measured:</b> skew p=4 max/mean load {value("sort_skew_strong",4,"imbalance"):.2f}.<br/><br/><b>Insight:</b> more samples cannot divide identical keys with upper-bound partitioning.<br/><br/><b>Next:</b> balance duplicate keys using global tie-breaking indices.'),
]
for filename,extended in [('Group04_slides.pdf',False),('Group04_slides_extended.pdf',True)]:
    slides=canvas.Canvas(str(ROOT/'slides'/filename),pagesize=(960,540))
    slides.setTitle('Group 04 - Parallel Case Studies ('+('seven' if extended else 'six')+' slides)')
    slides.setAuthor('Group 04 project presentation')
    introduction(compact=not extended)
    if extended: overview()
    for index,(title,subtitle,figure,labels,insight) in enumerate(content,3 if extended else 2):
        slide_title(title,subtitle,index);boxes(labels)
        figure_path=ROOT/'report'/'figures'/f'{figure}.png'
        if not figure_path.exists(): figure_path=FIGURES/f'{figure}.png'
        slides.drawImage(str(figure_path),28,85,width=565,height=250,mask='auto')
        slides.setFillColor(colors.HexColor('#BDD6E1'));slides.setFont('Helvetica',11)
        slides.drawString(36,65,'Median of 3 runs per point; whiskers show the observed minimum and maximum.')
        slide_text(insight,620,335,305,16)
        slides.showPage()
    slides.save()
print(f'Generated six-slide and seven-slide PDFs. Report format: {"LaTeX (existing PDF preserved)" if report.with_suffix(".tex").exists() else "ReportLab"}. Groups: {len(rows)}; measured runs: {sum(map(len,rows.values()))}.')
