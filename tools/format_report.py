"""Format the preserved report content as LaTeX; use only existing measurements."""
from pathlib import Path
from collections import defaultdict
from statistics import median
from html.parser import HTMLParser
import csv
import html
import json
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/'report'
FIG=REPORT/'figures'
FIG.mkdir(exist_ok=True)
blocks=json.loads((REPORT/'content.json').read_text(encoding='utf-8'))
data=defaultdict(list)
for line in (ROOT/'results/raw.tsv').read_text().splitlines():
    label,p,repeat,output=line.split('\t',3)
    data[label,int(p)].append(dict(t.split('=',1) for t in output.split() if '=' in t))
def v(label,p,field='seconds'):
    return median(float(r[field]) for r in data[label,p])
P=[1,2,4,6]
references={'heat1024_strong':'heat1024_serial','heat4096_strong':'heat4096_serial',
            'kmer8_atomic':'kmer8_serial','kmer8_private':'kmer8_serial',
            'kmer12_atomic':'kmer12_serial','kmer12_private':'kmer12_serial',
            'mc_good':'mc_serial','pagerank_strong':'pagerank_serial',
            'sort_uniform_strong':'sort_uniform_serial','sort_skew_strong':'sort_skew_serial'}
def baseline(label,p):
    return v(references[label],p if label.startswith('sort_') else 1)
def speedup(label,p): return baseline(label,p)/v(label,p)

# Colours keep the same semantic roles in figures, diagrams and LaTeX boxes.
C={'input':'#527D83','process':'#35536B','result':'#71839B','alert':'#D65A31','neutral':'#58616D','bg':'#F1F5F8'}
plt.rcParams.update({'font.size':10.5,'axes.labelsize':10.5,'legend.fontsize':10,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'text.color':'#20252B','axes.labelcolor':'#20252B',
                     'axes.edgecolor':C['neutral'],'xtick.color':C['neutral'],
                     'ytick.color':C['neutral'],'savefig.dpi':200})
captions={}
def save(name):
    plt.tight_layout(pad=.5)
    plt.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
    plt.savefig(FIG/f'{name}.png',dpi=200,bbox_inches='tight')
    plt.close()
def scaling(name,curves):
    fig,ax=plt.subplots(figsize=(5.65,2.6))
    points=P+[8] if name=='mc' else P
    for i,(label,legend,role,marker,style) in enumerate(curves):
        y=[speedup(label,p) for p in points]
        lower=[a-baseline(label,p)/max(float(r['seconds']) for r in data[label,p]) for a,p in zip(y,points)]
        upper=[baseline(label,p)/min(float(r['seconds']) for r in data[label,p])-a for a,p in zip(y,points)]
        ax.errorbar(points,y,yerr=[lower,upper],color=C[role],marker=marker,linestyle=style,
                    label=legend,capsize=2,linewidth=1.6,markersize=4)
    ax.plot(points,points,'--',color=C['neutral'],linewidth=1,label='Ideal S=p')
    ax.set_xticks(points);ax.set_xlabel('Workers (count; 8 uses SMT)' if name=='mc' else 'Workers (count)')
    ax.set_ylabel('Speedup over serial (times)')
    ax.legend(loc='upper left',ncol=2,frameon=False,columnspacing=.8,handlelength=1.6)
    selected=curves[1] if name=='kmer' else curves[0]
    label=selected[0]
    target=points[-1]
    y=speedup(label,target)
    annotation=(selected[1]+'\n' if len(curves)>1 else '')+f'{target} workers: {y:.2f}x'
    ax.annotate(annotation,xy=(target,y),xytext=(.57,.43),textcoords='axes fraction',
                color=C['alert'],fontsize=10,arrowprops={'arrowstyle':'->','color':C['alert'],'lw':1})
    ax.set_ylim(bottom=0)
    save(name)

scaling('heat',[('heat1024_strong','N=1024','process','o','-'),('heat4096_strong','N=4096','result','s','--')])
scaling('kmer',[('kmer8_atomic','k=8 atomic','process','o','-'),('kmer8_private','k=8 private','result','o','-'),
                ('kmer12_atomic','k=12 atomic','process','s','--'),('kmer12_private','k=12 private','result','s','--')])
scaling('mc',[('mc_good','100 million paths','result','o','-')])
scaling('pagerank',[('pagerank_strong','Pull, static','result','o','-')])
scaling('sort',[('sort_uniform_strong','Uniform','process','o','-'),('sort_skew_strong','90% zeros','result','s','--')])

fig,ax=plt.subplots(figsize=(5.65,2.45))
weak=[v('heat_weak',p)/v('heat_weak',1) for p in P]
ax.plot(P,weak,'o-',color=C['result'],label='Measured T(p)/T(1)')
ax.axhline(1,ls='--',color=C['neutral'],label='Ideal constant time')
ax.set_xticks(P);ax.set_xlabel('MPI ranks (count)');ax.set_ylabel('Normalized elapsed time (times)')
ax.legend(frameon=False,loc='upper left')
ax.annotate(f'p=6: {weak[-1]:.2f}x time',xy=(6,weak[-1]),xytext=(.48,.43),textcoords='axes fraction',color=C['alert'],arrowprops={'arrowstyle':'->','color':C['alert']})
save('weak')

degrees=list(csv.DictReader((ROOT/'results/degrees.csv').open()))
dx=[int(r['degree']) for r in degrees if int(r['degree'])>0]
dy=[int(r['vertices']) for r in degrees if int(r['degree'])>0]
fig,ax=plt.subplots(figsize=(5.65,2.45))
ax.loglog(dx,dy,'o',color=C['input'],markersize=2.6)
ax.set_xlabel('Incoming degree (edges)');ax.set_ylabel('Vertices (count)')
ax.annotate(f'Maximum degree: {max(dx):,}',xy=(dx[-1],dy[-1]),xytext=(.41,.75),textcoords='axes fraction',color=C['alert'],arrowprops={'arrowstyle':'->','color':C['alert']})
save('degrees')

fig,ax=plt.subplots(figsize=(5.65,2.7))
bottom=np.zeros(len(P))
shades=['#DCE5ED','#A7BACB','#718CA6','#35536B']
for i,(field,label) in enumerate([('local_sort_seconds','Local sort'),('sampling_seconds','Sampling / partition'),('exchange_seconds','Exchange / allocation'),('final_sort_seconds','Final sort')]):
    values=np.array([v('sort_uniform_strong',p,field) for p in P])
    ax.bar([str(p) for p in P],values,bottom=bottom,color=shades[i],edgecolor=C['neutral'],linewidth=.5,hatch=['','//','..','xx'][i],label=label)
    bottom+=values
ax.set_xlabel('MPI ranks (count)');ax.set_ylabel('Phase time (seconds)')
ax.legend(frameon=False,ncol=2,loc='upper center',bbox_to_anchor=(.5,1.30),handlelength=1.2,columnspacing=.8)
ax.annotate(f'{bottom[-1]:.2f} s total phases',xy=(3,bottom[-1]),xytext=(.48,.61),textcoords='axes fraction',color=C['alert'],arrowprops={'arrowstyle':'->','color':C['alert']})
save('phases')

captions={
'heat':f'Heat speedup remains limited at both grid sizes; N=4096 reaches {speedup("heat4096_strong",6):.2f} times serial speed at six ranks. Both curves use 500 steps, 2-D overlap and K=50.',
'kmer':'The atomic and private histogram methods scale differently for k=8 and k=12. All runs count 100 million bases; circle/square markers distinguish k and solid/dashed lines remain readable in greyscale.',
'mc':f'Monte Carlo reaches {speedup("mc_good",6):.2f} times serial speed at six workers for 100 million paths. Eight workers use simultaneous multithreading on six physical cores.',
'pagerank':f'Contiguous static pull reaches {speedup("pagerank_strong",6):.2f} times serial speed at six workers for a one-million-edge graph; scheduling alternatives are compared separately.',
'sort':f'Duplicate-heavy sample sort scales worse than uniform input: six-rank speedups are {speedup("sort_skew_strong",6):.2f} and {speedup("sort_uniform_strong",6):.2f}, respectively, for 100 million keys.',
'weak':f'Heat weak-scaling time rises to {weak[-1]:.2f} times the one-rank time despite approximately 512 squared interior cells per rank. Each point is the median of three runs.',
'degrees':'The fixed-seed R-MAT graph has highly unequal incoming degrees. Zero-degree vertices remain in the solver but are excluded from these logarithmic axes.',
'phases':'Uniform sample sort divides its time among local sorting, sampling/partitioning, exchange/allocation and final sorting. Bars show component medians from three runs on each run\'s slowest-total rank; they exclude the ending barrier.'}
for name in ['heat','kmer','mc','pagerank','sort']:
    captions[name]+=' Points are medians of three runs; whiskers span the minimum and maximum.'

def escape(text):
    mapping={'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_\allowbreak{}','{':r'\{','}':r'\}',
             '~':r'\textasciitilde{}','^':r'\textasciicircum{}','²':r'\textsuperscript{2}',
             '–':'--','—':'---','≤':r'$\leq$','≥':r'$\geq$','×':r'$\times$'}
    return ''.join(mapping.get(c,c) for c in str(text)).replace('\n',r'\newline{} ').replace('=',r'=\allowbreak{}').replace('/',r'/\allowbreak{}')

class Markup(HTMLParser):
    def __init__(self): super().__init__();self.parts=[]
    def handle_starttag(self,tag,attrs):
        if tag=='b': self.parts.append(r'\textbf{')
        elif tag=='i': self.parts.append(r'\emph{')
        elif tag=='br': self.parts.append(r'\newline{} ')
        elif tag=='link':
            href=dict(attrs)['href'];self.parts.append(r'\href{'+href.replace('%',r'\%')+'}{')
    def handle_endtag(self,tag):
        if tag in ['b','i','link']:self.parts.append('}')
    def handle_data(self,text): self.parts.append(escape(text))
def tex(text):
    p=Markup();p.feed(str(text));return ''.join(p.parts)

preamble=r'''\documentclass[11pt,a4paper]{article}
\usepackage[margin=2.2cm,headheight=14pt]{geometry}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{xcolor,graphicx,amsmath,amssymb,booktabs,array,longtable,listings}
\usepackage{tcolorbox,caption,tikz,pgfplots,microtype,enumitem,fancyhdr}
\usepackage{seqsplit,placeins}
\usepackage[colorlinks=true]{hyperref}
\pgfplotsset{compat=1.18}
\usetikzlibrary{arrows.meta,positioning,shapes.geometric}
\definecolor{cInput}{HTML}{527D83} % supplied input and measured input distributions
\definecolor{cProcess}{HTML}{35536B} % computation, method and comparison baseline
\definecolor{cResult}{HTML}{71839B} % computed output and measured outcomes
\definecolor{cAlert}{HTML}{D65A31} % sole high-saturation accent: point to notice
\definecolor{cNeutral}{HTML}{58616D} % axes, arrows and structural labels
\definecolor{cBg}{HTML}{F1F5F8} % quiet panel background
\hypersetup{linkcolor=cProcess,urlcolor=cProcess,citecolor=cProcess,pdftitle={Group 04: Parallel Case Studies},pdfauthor={Group 04 project report}}
\pagestyle{fancy}\fancyhf{}
\fancyhead[L]{\small Group 04 --- Parallel case studies}
\fancyhead[R]{\small\thepage}
\renewcommand{\headrulewidth}{0.3pt}
\setlength{\parindent}{0pt}
\setlength{\parskip}{4pt plus 1pt}
\setlength{\emergencystretch}{2em}
\setlength{\textfloatsep}{9pt plus 2pt minus 2pt}
\setlength{\intextsep}{8pt plus 2pt minus 2pt}
\captionsetup{font=small,labelfont=bf,skip=4pt}
\setlist{nosep,leftmargin=*}
\tcbset{colback=cProcess!5,colframe=cProcess,boxrule=.5pt,arc=1mm,left=2mm,right=2mm,top=1mm,bottom=1mm,fonttitle=\bfseries}
\lstset{basicstyle=\small\ttfamily,breaklines=true,backgroundcolor=\color{cBg},frame=single,rulecolor=\color{cNeutral}}
\newcommand{\TODO}[1]{\textcolor{cAlert}{[TODO: #1]}}
\newcommand{\term}[2]{\textbf{#1}: #2.}
\setcounter{tocdepth}{1}
\begin{document}
\begin{center}
{\LARGE\bfseries Group 04 --- Parallel case studies}\par
{\large CS3006 --- Assignment 3}\par
Fatima Malik (23L-0826) \quad Ibrahim Salman (23L-0723)\par
Measurements and analysis, 23 November 2026
\end{center}
\tableofcontents
\vspace{4pt}
\begin{tcolorbox}[title=Reading the visual language]
Muted teal identifies inputs, dark steel identifies computation, and slate identifies results.
The warm accent marks one finding to notice. Labels, marker shapes, dashes and hatching
preserve meaning without colour. Definitions explain terms at their first relevant section.
\end{tcolorbox}
\clearpage
\section*{Introduction}\addcontentsline{toc}{section}{Introduction}
'''
out=[preamble]
definitions={
0:r'\term{MPI}{Message Passing Interface, used by separate processes to exchange data} \term{rank}{an MPI process identifier} \term{OpenMP}{directives that divide work among threads sharing memory} \term{thread}{a worker within one process} \term{baseline}{the serial program used for comparison} \term{median}{the middle measured time} \term{SMT}{simultaneous multithreading, where one physical core runs multiple hardware threads} \term{WSL}{Windows Subsystem for Linux}',
1:r'\term{stencil}{a repeated neighbor-based update rule} \term{halo}{copied neighbor cells around an owned block} \term{decomposition}{division of the grid among workers} \term{blocking}{waiting for an operation to complete before continuing} \term{nonblocking}{starting communication and completing it later} \term{Allreduce}{an MPI collective that combines values and returns the result to every rank} \term{convergence}{changes becoming sufficiently small}',
2:r'\term{k-mer}{a DNA substring of length k} \term{histogram}{a count for each possible substring} \term{atomic update}{an indivisible shared-counter update} \term{data race}{unsynchronized concurrent access with at least one write} \term{cache line}{the unit of data transferred through a processor cache} \term{false sharing}{workers modifying distinct values on one cache line, causing coherence traffic}',
3:r'\term{Monte Carlo}{estimation using pseudorandom samples} \term{variance}{how spread out values are around their average} \term{standard error}{estimated uncertainty in the sample mean} \term{Box--Muller}{a transform from uniform samples to normal samples} \term{closed form}{a direct mathematical formula} \term{reduction}{combining worker results} \term{Amdahl fraction}{the effective part of execution that does not benefit from more workers}',
4:r'\term{CSR}{compressed sparse row storage: offsets plus adjacent vertex identifiers} \term{degree}{the number of incident edges in the specified direction} \term{dangling vertex}{a vertex with no outgoing edges} \term{pull}{each destination reads its incoming contributions} \term{push}{each source adds contributions to destinations} \term{chunk}{a group of iterations scheduled together} \term{L1 residual}{the sum of absolute rank changes} \term{R-MAT}{a recursive matrix method for generating an irregular graph}',
5:r'\term{splitter}{a key value separating output buckets} \term{bucket}{the keys assigned to one destination rank} \term{oversampling}{collecting extra candidate splitters} \term{fingerprint}{a compact checksum used to detect likely data changes} \term{latency}{the fixed delay of a transfer} \term{bandwidth}{the amount of data transferred per second} \term{Alltoallv}{an MPI collective that exchanges different amounts of data between ranks}'
}
def definition(number):
    out.append(r'\begin{tcolorbox}[title=Terms used in this section]\small '+definitions[number]+r'\end{tcolorbox}')
definition(0)

section=0
table_number=0
chart_counts=defaultdict(int)
preserved=[]
for block_index,b in enumerate(blocks):
    kind=b['type'];args=b['args']
    preserved.append(block_index)
    if kind=='page':
        title=args[0]
        if title.startswith('Group 04'):continue
        number,body=title.split(' | ',1);number=int(number)
        if number!=section:
            out.append(r'\FloatBarrier')
            section=number
            if number<=5:
                base=body.split(':')[0]
                out.append(r'\section{'+escape(base)+'}')
                definition(number)
            else:out.append(r'\section{'+escape(body)+'}')
        if number<=5:
            out.append(r'\subsection{'+escape(body.split(': ',1)[1].capitalize())+'}')
    elif kind=='heading':out.append(r'\subsubsection*{'+tex(args[0])+'}')
    elif kind=='para':
        # All source paragraphs are retained, including their figures and exact values.
        paragraph=args[0]
        if paragraph.startswith(('At N=1024,','For 100 million paths,','The fastest tested','Among the representative')):
            sentence,separator,remaining=paragraph.partition('. ')
            paragraph='<b>'+sentence+('.' if separator else '')+'</b>'+(' '+remaining if separator else '')
        out.append(tex(paragraph)+'\n')
        if args[0].startswith('S(p)='):
            out.append(r'\[S(p)=\frac{T_{\mathrm{serial}}}{T_{\mathrm{parallel}}},\qquad E(p)=\frac{S(p)}{p}.\] Speedup compares serial and parallel elapsed time; efficiency divides speedup by the worker count.')
    elif kind=='table':
        table_number+=1
        cells=args[0];columns=len(cells[0]);widths=args[1] if len(args)>1 else [1]*columns
        numeric=[all(re.fullmatch(r'[\d.,eE+% /-]+(?: MiB)?',str(row[j])) for row in cells[1:]) for j in range(columns)]
        specs=[]
        for width,is_numeric in zip(widths,numeric):
            frac=width/sum(widths)
            alignment=r'\raggedleft' if is_numeric else r'\raggedright'
            specs.append('>{'+alignment+r'\arraybackslash}p{\dimexpr '+f'{frac:.5f}'+r'\linewidth-'+f'{8*(columns-1)*frac:.3f}pt'+r'\relax}')
        out.append(r'\begin{center}\small\setlength{\tabcolsep}{4pt}\renewcommand{\arraystretch}{1.12}')
        out.append(r'\begin{tabular}{@{}'+''.join(specs)+r'@{}}\toprule')
        out.append(' & '.join(r'\textbf{'+tex(c)+'}' for c in cells[0])+r'\\\midrule')
        for row in cells[1:]:
            out.append(' & '.join(tex(c).replace(',',r',\allowbreak{}') if len(str(c))>28 else tex(c) for c in row)+r'\\')
        out.append(r'\bottomrule\end{tabular}\end{center}')
        if 'Private memory' in cells[0]:
            out.append(r'{\small\textbf{Rounding note.} The supplied table rounds memory to whole MiB. Its ``0 MiB'' entry represents 512 KiB after rounding, not a zero allocation.}')
    elif kind=='chart':
        name=args[0];label='fig:'+name
        out.append(r'The measured comparison is shown in Figure~\ref{'+label+'}.')
        out.append(r'\begin{figure}[htbp]\centering\includegraphics[width=0.9\textwidth]{figures/'+name+'.pdf}')
        out.append(r'\caption{'+escape(captions[name])+r'}\label{'+label+r'}\end{figure}')
    elif kind=='diagram':
        text=args[0]
        labels={1:['Owned grid blocks','Neighbor halos','Interior / border update'],2:['Window-start ranges','Worker histograms','Merge bins'],3:['Disjoint path counters','Block moments','Ordered reduction'],4:['Owned vertices','CSR traversal','Ranks / L1 reduction'],5:['Rank-local keys','Splitters / Alltoallv','Owned sorted buckets']}[section]
        out.append(r'The ownership and data flow are summarized schematically in Figure~\ref{fig:flow'+str(section)+'}.')
        out.append(r'\begin{figure}[htbp]\centering\begin{tikzpicture}[font=\small,>=Stealth,node distance=6mm, every node/.style={align=center,text width=3.4cm,minimum height=9mm,inner sep=3pt}]')
        out.append(r'\node[draw=cInput,fill=cInput!8,rounded corners=2pt] (a) {'+tex(labels[0])+'};')
        out.append(r'\node[draw=cProcess,fill=cProcess!6,right=of a] (b) {'+tex(labels[1])+'};')
        out.append(r'\node[draw=cResult,fill=cResult!10,rounded corners=2pt,right=of b] (c) {'+tex(labels[2])+'};')
        out.append(r'\draw[->,cNeutral] (a)--(b);\draw[->,cNeutral] (b)--(c);\end{tikzpicture}')
        out.append(r'\caption{Schematic: ownership is separated from the communication or reduction that combines worker results.}\label{fig:flow'+str(section)+r'}\end{figure}')
        out.append(escape(text)+'\n')
    else:raise ValueError(kind)
out.append(r'\FloatBarrier\end{document}')
(REPORT/'Group04_report.tex').write_text('\n\n'.join(out),encoding='utf-8')
(REPORT/'formatting_audit.json').write_text(json.dumps({'source_blocks':len(blocks),'preserved_block_indices':preserved,'paragraphs':sum(b['type']=='para' for b in blocks),'tables':table_number,'measured_plots':8,'schematic_diagrams':5,'measurement_source':'../results/raw.tsv','new_measurements':False},indent=2))
print('Wrote report/Group04_report.tex; preserved all',len(blocks),'source blocks.')
