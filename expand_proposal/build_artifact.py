# Build a self-contained, theme-aware HTML report with the figures embedded as data URIs.
import base64, pathlib

FIG = pathlib.Path("figures")
def data_uri(p):
    b = base64.b64encode(pathlib.Path(p).read_bytes()).decode()
    return f"data:image/png;base64,{b}"

f1 = data_uri(FIG/"fig1_scaling_p.png")
f2 = data_uri(FIG/"fig2_scaling_n.png")
f3 = data_uri(FIG/"fig3_parallel.png")
f4 = data_uri(FIG/"fig4_su_cost.png")

HTML = f"""<title>cis-Heritability — Expanse Compute &amp; Scaling Documentation</title>
<style>
:root{{
  --ground:#fbfcfc; --surface:#ffffff; --surface-2:#f2f6f6;
  --ink:#152329; --muted:#5b6a71; --faint:#8a979d;
  --accent:#0e7c86; --accent-deep:#0b5560; --amber:#a9660e;
  --line:#dde6e8; --line-strong:#c4d1d4;
  --tealwash:#eef6f6; --amberwash:#fbf2e6;
  --maxw:940px;
}}
@media (prefers-color-scheme:dark){{
  :root{{
    --ground:#0e1618; --surface:#15211f; --surface-2:#1b2a28;
    --ink:#e7eef0; --muted:#9fb0b4; --faint:#7b8c90;
    --accent:#3fb8c2; --accent-deep:#7fd6dd; --amber:#e3a94e;
    --line:#25332f; --line-strong:#33453f;
    --tealwash:#132523; --amberwash:#241d10;
  }}
}}
:root[data-theme="light"]{{
  --ground:#fbfcfc; --surface:#ffffff; --surface-2:#f2f6f6; --ink:#152329; --muted:#5b6a71; --faint:#8a979d;
  --accent:#0e7c86; --accent-deep:#0b5560; --amber:#a9660e; --line:#dde6e8; --line-strong:#c4d1d4;
  --tealwash:#eef6f6; --amberwash:#fbf2e6;
}}
:root[data-theme="dark"]{{
  --ground:#0e1618; --surface:#15211f; --surface-2:#1b2a28; --ink:#e7eef0; --muted:#9fb0b4; --faint:#7b8c90;
  --accent:#3fb8c2; --accent-deep:#7fd6dd; --amber:#e3a94e; --line:#25332f; --line-strong:#33453f;
  --tealwash:#132523; --amberwash:#241d10;
}}
*{{box-sizing:border-box}}
body{{background:var(--ground);color:var(--ink);
  font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
  line-height:1.6;margin:0;-webkit-font-smoothing:antialiased}}
.serif{{font-family:Charter,"Bitstream Charter","Iowan Old Style",Georgia,"Times New Roman",serif}}
.mono{{font-family:ui-monospace,"SF Mono","Cascadia Code",Menlo,monospace}}
.wrap{{max-width:var(--maxw);margin:0 auto;padding:0 24px}}

header.top{{border-bottom:1px solid var(--line);background:
  linear-gradient(180deg,var(--tealwash),transparent)}}
header.top .wrap{{padding-top:52px;padding-bottom:30px}}
.eyebrow{{font-size:12.5px;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);
  font-weight:600;margin:0 0 14px}}
h1{{font-family:Charter,"Iowan Old Style",Georgia,serif;font-weight:700;
  font-size:clamp(28px,4.4vw,42px);line-height:1.12;margin:0 0 14px;text-wrap:balance;letter-spacing:-.01em}}
.sub{{font-size:17px;color:var(--muted);max-width:60ch;margin:0}}
.meta{{margin-top:22px;display:flex;flex-wrap:wrap;gap:8px 20px;font-size:13.5px;color:var(--faint)}}
.meta b{{color:var(--ink);font-weight:600}}

/* metric tiles */
.tiles{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:14px;margin:34px 0 6px}}
.tile{{background:var(--surface);border:1px solid var(--line);border-radius:11px;padding:16px 16px 14px;
  position:relative;overflow:hidden}}
.tile::before{{content:"";position:absolute;left:0;top:0;bottom:0;width:3px;background:var(--accent)}}
.tile .k{{font-family:ui-monospace,Menlo,monospace;font-size:23px;font-weight:600;color:var(--ink);
  letter-spacing:-.02em;font-variant-numeric:tabular-nums}}
.tile .l{{font-size:12.5px;color:var(--muted);margin-top:5px;line-height:1.35}}

main .wrap{{padding-top:8px;padding-bottom:64px}}
section{{padding-top:40px;scroll-margin-top:20px}}
.sec-head{{display:flex;gap:14px;align-items:baseline;border-bottom:2px solid var(--line-strong);
  padding-bottom:9px;margin-bottom:20px}}
.sec-num{{font-family:ui-monospace,Menlo,monospace;font-size:14px;color:var(--accent);font-weight:600;
  flex:none;padding-top:3px}}
h2{{font-family:Charter,"Iowan Old Style",Georgia,serif;font-size:24px;font-weight:700;margin:0;letter-spacing:-.01em}}
h3{{font-size:15px;font-weight:650;margin:26px 0 8px;color:var(--ink)}}
p{{margin:0 0 14px}}
.lead{{font-size:16.5px}}
a{{color:var(--accent-deep)}}

table{{width:100%;border-collapse:collapse;font-size:13.5px;margin:6px 0 8px}}
.scroll{{overflow-x:auto;border:1px solid var(--line);border-radius:10px;margin:14px 0 20px}}
.scroll table{{margin:0}}
th,td{{text-align:left;padding:9px 13px;border-bottom:1px solid var(--line);vertical-align:top}}
thead th{{background:var(--surface-2);font-size:12px;letter-spacing:.03em;text-transform:uppercase;
  color:var(--muted);font-weight:650;white-space:nowrap}}
tbody tr:last-child td{{border-bottom:none}}
td.num,th.num{{text-align:right;font-family:ui-monospace,Menlo,monospace;font-variant-numeric:tabular-nums;white-space:nowrap}}
tbody tr:hover td{{background:var(--surface-2)}}
.total td{{font-weight:700;border-top:2px solid var(--line-strong);background:var(--tealwash)}}

figure{{margin:20px 0;border:1px solid var(--line);border-radius:12px;overflow:hidden;background:var(--surface)}}
figure img{{width:100%;display:block;background:#fff}}
figcaption{{font-size:13px;color:var(--muted);padding:11px 16px;border-top:1px solid var(--line);background:var(--surface)}}
figcaption b{{color:var(--ink)}}

.callout{{border:1px solid var(--line);border-left:4px solid var(--amber);background:var(--amberwash);
  border-radius:9px;padding:14px 18px;margin:16px 0;font-size:14px}}
.callout .t{{font-weight:700;color:var(--amber);font-size:12.5px;letter-spacing:.05em;text-transform:uppercase;margin-bottom:5px}}
.note{{border-left-color:var(--accent);background:var(--tealwash)}}
.note .t{{color:var(--accent-deep)}}
ul{{margin:6px 0 16px;padding-left:20px}}
li{{margin:5px 0}}
.kicker{{color:var(--accent-deep);font-weight:650}}
code{{font-family:ui-monospace,Menlo,monospace;font-size:.9em;background:var(--surface-2);
  padding:1px 6px;border-radius:5px;border:1px solid var(--line)}}
footer{{border-top:1px solid var(--line);color:var(--faint);font-size:13px;padding:26px 0 60px}}
footer .wrap{{display:flex;flex-wrap:wrap;gap:6px 18px;align-items:center}}
.pill{{display:inline-block;font-size:11.5px;background:var(--surface-2);border:1px solid var(--line);
  border-radius:20px;padding:2px 11px;color:var(--muted);margin:0 4px 4px 0}}
</style>

<header class="top"><div class="wrap">
  <p class="eyebrow">ACCESS Maximize · Compute usage &amp; scaling documentation</p>
  <h1 class="serif">Fine-mapping-based estimation of <em>cis</em>-heritability</h1>
  <p class="sub">Per-project documentation of homegrown code, third-party software, storage, and measured
  performance &amp; scaling on Expanse-class AMD&nbsp;EPYC hardware — the evidence base for the Service-Unit request.</p>
  <div class="meta"><span>Amariuta Lab · <b>Sunan Xu</b></span><span>Prepared <b>July&nbsp;2026</b></span>
  <span>All figures <b>measured</b>, not estimated</span><span class="mono">env: dsc180 · SDSC EPYC&nbsp;7713</span></div>

  <div class="tiles">
    <div class="tile"><div class="k">∝&nbsp;p<span style="font-size:.7em">1.7</span></div><div class="l">runtime scales with SNPs per region (measured)</div></div>
    <div class="tile"><div class="k">∝&nbsp;p<span style="font-size:.7em">1.6</span></div><div class="l">peak memory scales with region size</div></div>
    <div class="tile"><div class="k">3&ndash;100s</div><div class="l">single-core cost per model (p&nbsp;=&nbsp;1k&ndash;8k)</div></div>
    <div class="tile"><div class="k">100%</div><div class="l">embarrassingly parallel across genes/loci</div></div>
    <div class="tile"><div class="k">~380k</div><div class="l">estimated Service Units / year</div></div>
  </div>
</div></header>

<main><div class="wrap">

<section>
  <div class="sec-head"><span class="sec-num">01</span><h2 class="serif">What the project does</h2></div>
  <p class="lead">We benchmark two estimators of the <span class="kicker"><em>cis</em>-heritability of gene expression</span> —
  GCTA/REML (individual-level, the field standard) and a <b>SuSiE fine-mapping estimator we built</b> (summary-statistic based) —
  to determine <em>when</em> fine-mapping yields more reliable <em>cis</em>-h² than GCTA.</p>
  <p>Because ground-truth heritability is never known in real data, the compute core is a <b>simulation engine</b>: it draws a
  realistic local LD structure from 1000&nbsp;Genomes, simulates expression with a <em>known</em> true h² under many genetic
  architectures, runs both estimators, and compares their error. Findings are then confirmed on real GTEx expression data. As a
  method-evaluation project, cost is dominated by sweeping a large architecture grid across thousands of replicate loci to obtain
  confidence intervals on <b>bias, power, and type-I error</b>.</p>
</section>

<section>
  <div class="sec-head"><span class="sec-num">02</span><h2 class="serif">Homegrown code &mdash; and why it was necessary</h2></div>
  <p>~<b>6,000 lines across 32 scripts</b> (Python + R + shell). Small in size, but the scientific core: no existing tool performs
  this end-to-end simulate-and-compare workflow.</p>
  <div class="scroll"><table>
    <thead><tr><th>Component</th><th>Lang</th><th>What it does</th><th>Why off-the-shelf tools don't cover it</th></tr></thead>
    <tbody>
      <tr><td><b>Simulation engine</b><br><span class="mono" style="font-size:11px">simulations_test.py&nbsp;+ variants</span></td><td>Python</td><td>Draws a real <em>cis</em> region, builds &amp; regularizes LD, Cholesky-samples genotypes, simulates effects to a <b>known</b> h², runs both estimators, records error</td><td>No public tool simulates expression from a fixed known <em>cis</em>-h² <em>and</em> benchmarks GCTA vs. fine-mapping in one loop</td></tr>
      <tr><td><b>SuSiE&rarr;h² estimator</b><br><span class="mono" style="font-size:11px">run_susie.R</span></td><td>R</td><td>Runs <code>susie_rss</code>, does allele-flip QC, converts posterior effects to h² = Σ(μ²+σ²) with an analytic SE</td><td>SuSiE returns fine-mapping posteriors, <b>not</b> a heritability estimate or its SE — that derivation is ours</td></tr>
      <tr><td><b>GCTA harness</b></td><td>Python</td><td>Builds the GRM from simulated genotypes, writes GCTA inputs, runs REML, parses output</td><td>Gluing simulated matrices into GCTA's formats at scale is custom</td></tr>
      <tr><td><b>Grid + job orchestration</b></td><td>shell</td><td>Sweeps the parameter grid, manages temp dirs, retries non-converged replicates, cleans up</td><td>Bespoke to our parameter design</td></tr>
    </tbody>
  </table></div>
  <div class="callout note"><div class="t">Core design constraint we had to enforce ourselves</div>
  The LD used to <em>simulate</em> genotypes must be the exact same ordered LD handed to the fine-mapper, or the comparison is
  uninterpretable. That alignment is baked into our code — a generic pipeline does not provide it.</div>
</section>

<section>
  <div class="sec-head"><span class="sec-num">03</span><h2 class="serif">Third-party software required</h2></div>
  <div class="scroll"><table>
    <thead><tr><th>Software</th><th>Version</th><th>Role</th></tr></thead>
    <tbody>
      <tr><td><b>GCTA</b></td><td class="mono">1.94.1</td><td>REML variance-component h² — benchmark estimator</td></tr>
      <tr><td><b>SuSiE (susieR)</b></td><td class="mono">CRAN</td><td>Summary-statistic fine-mapping (<code>susie_rss</code>, IBSS)</td></tr>
      <tr><td><b>PLINK</b></td><td class="mono">1.90</td><td>Extract <em>cis</em> regions; compute pairwise LD</td></tr>
      <tr><td><b>Python stack</b></td><td class="mono">—</td><td>NumPy/SciPy (LD, Cholesky), statsmodels (OLS), pandas-plink, scikit-learn</td></tr>
      <tr><td><b>R stack</b></td><td class="mono">—</td><td>Rfast, data.table, dplyr, optparse inside <code>run_susie.R</code></td></tr>
    </tbody>
  </table></div>
  <p>One conda environment; <b>no GPU</b>. The workload is CPU + memory-bandwidth bound (dense linear algebra on p×p LD
  matrices), which maps directly onto Expanse standard compute nodes.</p>
</section>

<section>
  <div class="sec-head"><span class="sec-num">04</span><h2 class="serif">Storage &amp; file-set profile</h2></div>
  <div class="scroll"><table>
    <thead><tr><th>Category</th><th class="num">Size</th><th>File-set shape</th></tr></thead>
    <tbody>
      <tr><td>Homegrown code (6k lines)</td><td class="num">~0.5 MB</td><td>Few small text files, version-controlled</td></tr>
      <tr><td>Reference genotypes — 1000G EUR, chr1&ndash;22</td><td class="num">14 GB</td><td><b>Few large files</b>; read-only, reused by every job</td></tr>
      <tr><td>Per-run result tables</td><td class="num">0.1&ndash;0.3 MB</td><td>Many small append-only files (one row / model)</td></tr>
      <tr><td><b>Transient LD matrices</b> (<span class="mono">*.ld</span>)</td><td class="num">0.1&ndash;12 GB <i>each</i></td><td>O(p²); written then deleted per model — dominant I/O object</td></tr>
      <tr><td>Transient GRM / sumstats / phenotype</td><td class="num">KB&ndash;MB</td><td>Many small files, created &amp; deleted per model</td></tr>
    </tbody>
  </table></div>
  <p><b>Data movement:</b> the 14&nbsp;GB reference stages once to node-local NVMe and is read many times; all heavy transient
  files live on and are deleted from <b>node-local NVMe</b>, never the parallel filesystem. Only the small final table returns to
  project storage. Persistent storage is modest (tens of GB); the pressure is on transient scratch and on <b>metadata</b>
  (file create/delete rate).</p>
</section>

<section>
  <div class="sec-head"><span class="sec-num">05</span><h2 class="serif">Performance &amp; scaling <span style="font-size:14px;color:var(--faint);font-weight:400">(measured, single-core, threads pinned)</span></h2></div>

  <h3>5a · The dominant cost driver is SNPs per region (p), not sample size</h3>
  <figure><img src="{f1}" alt="Per-stage runtime and memory versus number of SNPs in the cis-region, log-log">
    <figcaption><b>Figure 1.</b> At fixed n&nbsp;=&nbsp;500, per-model runtime rises as <b>p<sup>1.7</sup></b> and peak memory as
    <b>p<sup>1.6</sup></b>. SuSiE's dense p×p operations and the Cholesky dominate; GCTA and OLS are minor. Widening the
    <em>cis</em>-window gets expensive super-linearly.</figcaption></figure>
  <div class="scroll"><table>
    <thead><tr><th class="num">p (SNPs)</th><th class="num">Build LD</th><th class="num">SuSiE</th><th class="num">PLINK LD</th><th class="num">GCTA</th><th class="num">Total / model</th><th class="num">Peak RAM</th><th class="num">LD file</th></tr></thead>
    <tbody>
      <tr><td class="num">1,544</td><td class="num">0.2s</td><td class="num">4.9s</td><td class="num">0.9s</td><td class="num">1.4s</td><td class="num"><b>8.0s</b></td><td class="num">0.25 GB</td><td class="num">93 MB</td></tr>
      <tr><td class="num">3,567</td><td class="num">0.9s</td><td class="num">12.5s</td><td class="num">4.3s</td><td class="num">1.4s</td><td class="num"><b>20.6s</b></td><td class="num">0.81 GB</td><td class="num">496 MB</td></tr>
      <tr><td class="num">8,227</td><td class="num">6.6s</td><td class="num">59.8s</td><td class="num">23.7s</td><td class="num">1.8s</td><td class="num"><b>96.0s</b></td><td class="num">3.27 GB</td><td class="num">2.8 GB</td></tr>
      <tr><td class="num">17,287</td><td class="num">45.2s</td><td class="num">253.0s</td><td class="num">104.4s</td><td class="num">3.4s</td><td class="num"><b>417.8s</b></td><td class="num">12.9 GB</td><td class="num">12.3 GB</td></tr>
    </tbody>
  </table></div>

  <h3>5b · Sample size barely moves the cost</h3>
  <figure><img src="{f2}" alt="Runtime versus sample size at fixed region size">
    <figcaption><b>Figure 2.</b> At fixed region (p&nbsp;=&nbsp;3,567), growing n from 100&rarr;1,000 changes total per-model time
    only 19.9s&rarr;25.0s. GCTA's REML grows ~O(n²) but is small; SuSiE and the LD build are <b>n-independent</b>. Scaling
    experiments therefore vary p / number-of-loci, not n.</figcaption></figure>

  <h3>5c · From 1 to 1,000 simulations — and why 1,000</h3>
  <p>The engine is a linear loop over independent replicate loci, so wall-time is linear in the number of simulations. Measured
  ground truth for one full architecture grid (30 models) on one region: <b>3.65 min</b> (7.3s/model). So <b>1,000 loci ≈ 61
  core-hours</b> at p≈1.5k, or <b>~165 SU</b> at a representative p≈3k.</p>
  <div class="callout note"><div class="t">Why 1,000 replicates is scientifically necessary</div>
  Power, bias and type-I error are means/proportions over replicate loci, so their CIs shrink as 1/√N. Monte-Carlo SE of a
  proportion ≈ √(p(1−p)/N): <b>N&nbsp;=&nbsp;100 → ±5%</b>, <b>N&nbsp;=&nbsp;1,000 → ±1.6%</b>. Resolving whether SuSiE's type-I
  error is 5% vs. 7%, or a few-percent bias difference across the h² grid, needs ~1,000 replicates — the standard in the
  fine-mapping / TWAS-simulation literature (twas_sim, SuSiE papers).</div>

  <h3>5d · Node-level parallel scaling (1 &rarr; 128 cores)</h3>
  <figure><img src="{f3}" alt="Throughput and parallel efficiency versus cores on one node">
    <figcaption><b>Figure 3.</b> Each locus is an independent Slurm-array task with zero communication, so scaling <b>across
    nodes is linear by construction</b>. Within one node, throughput keeps rising but sub-linearly: light 0.25&nbsp;GB tasks reach
    ~13.5× at 64 cores; heavy 0.8&nbsp;GB tasks saturate near ~6× as the dense p×p algebra becomes <b>memory-bandwidth bound</b>.
    The efficient strategy — a moderate number of tasks per node, spread across <b>many</b> nodes — is exactly what a Maximize
    allocation enables.</figcaption></figure>
</section>

<section>
  <div class="sec-head"><span class="sec-num">06</span><h2 class="serif">Parameters swept &amp; their scientific meaning</h2></div>
  <div class="scroll"><table>
    <thead><tr><th>Parameter</th><th>Levels</th><th>Why it matters</th></tr></thead>
    <tbody>
      <tr><td>True <em>cis</em>-h²</td><td class="mono">0.05, 0.1, 0.2, 0.4, 0.6, 0.8</td><td>Spans realistic expression heritabilities; tests bias across the range</td></tr>
      <tr><td>Number of causal SNPs</td><td class="mono">1, 3, 5, 10, 20</td><td>Sparse vs. polygenic architecture — where fine-mapping wins or fails</td></tr>
      <tr><td>Sample size <em>n</em></td><td class="mono">100, 300, 500, 800, 1000</td><td>Matches real eQTL study sizes (e.g. GTEx tissues)</td></tr>
      <tr><td>Effect-size distribution</td><td class="mono">normal, fixed, Pareto</td><td>Robustness to the unknown true effect-size prior</td></tr>
      <tr><td>LD architecture</td><td class="mono">regular, high-LD, low-LD</td><td>LD is the confounder fine-mapping is meant to resolve</td></tr>
      <tr><td>Allele-frequency bin</td><td class="mono">common, rare</td><td>Rare-variant regime stresses both estimators</td></tr>
      <tr><td>Replicate loci</td><td class="mono">up to 1,000</td><td>Confidence intervals on bias / power / type-I error</td></tr>
    </tbody>
  </table></div>
  <p>One fully-crossed sweep per replicate ≈ <b>2,700 models per locus</b>. End-targets: <b>~20,000 genes × up to 49 GTEx
  tissues</b> for real-data validation, extensible to <b>~100 GWAS traits</b>, plus a planned third estimator (<b>HEELS</b>) that
  multiplies the method axis.</p>
</section>

<section>
  <div class="sec-head"><span class="sec-num">07</span><h2 class="serif">Efficiency work &amp; filesystem bottlenecks</h2></div>
  <h3>Already done to minimize I/O and cost</h3>
  <ul>
    <li><b>All heavy transient I/O on node-local NVMe</b>, not Lustre — thousands of array tasks never overwhelm the parallel-FS metadata server.</li>
    <li><b>Aggressive per-model cleanup</b> keeps the transient-file count bounded across a full grid.</li>
    <li><b>Cholesky factored once per region</b> and reused for all genotype replicates (O(p³) paid once, not per replicate).</li>
    <li><b>Retry-with-reseed</b> on SuSiE non-convergence, so one bad replicate never wastes a node-hour.</li>
    <li><b>Append-only single result table</b> per campaign rather than thousands of tiny output files.</li>
    <li><b>Runs past the 48-hour wall limit:</b> independent, checkpointable loci — the result table skips already-completed rows, so campaigns split into requeued Slurm-array chunks and survive the wall-clock cap.</li>
  </ul>
  <div class="callout"><div class="t">Bottleneck 1 — the LD text file is O(p²) and huge</div>
  Measured <b>93 MB → 2.8 GB → 12.3 GB</b> per model as p grows (matching the 2.4&nbsp;GB <span class="mono">temp_ld.ld</span>
  still in the repo). This single object dominates I/O and must stay on local NVMe.</div>
  <div class="callout"><div class="t">Bottleneck 2 — redundant LD recomputation (planned fix)</div>
  The LD file is currently recomputed for every one of the 30 grid cells of a region, though the region's genotypes — and thus its
  LD — are identical across cells. Caching LD <b>once per region</b> removes ~29/30 of that I/O: at p=8k it saves ~23s and 2.8&nbsp;GB
  of write <em>per cell</em>.</div>
  <div class="callout"><div class="t">Bottleneck 3 — metadata churn</div>
  Each model creates/deletes ~5 small files; at hundreds of concurrent tasks this stresses the filesystem's create/delete rate more
  than its capacity — another reason it all lives on node-local NVMe.</div>
  <p style="font-size:13.5px;color:var(--muted)"><b>Independent of NAIRR:</b> this CPU-only simulation/fine-mapping workload is
  entirely separate from any NAIRR-pilot GPU/LLM usage.</p>
</section>

<section>
  <div class="sec-head"><span class="sec-num">08</span><h2 class="serif">Service-Unit budget <span style="font-size:14px;color:var(--faint);font-weight:400">(small experiment × multiplier)</span></h2></div>
  <p>We measured the true single-core cost of <em>one</em> model, then multiply by the number of models — no wasteful full-scale
  scaling runs. On Expanse, <b>1 SU = 1 core-hour</b>; a node-hour = 128 SU.</p>
  <figure><img src="{f4}" alt="Estimated Service Units per campaign, log scale bar chart">
    <figcaption><b>Figure 4.</b> Estimated SU per representative campaign at p≈3,000 (log scale). Each bar = measured per-model cost
    × number of models.</figcaption></figure>
  <div class="scroll"><table>
    <thead><tr><th>Campaign</th><th class="num"># models</th><th class="num">SU (core-hr)</th><th class="num">node-days</th></tr></thead>
    <tbody>
      <tr><td>Baseline sim campaign (1,000 loci × 5 ncausal × 6 h²)</td><td class="num">30,000</td><td class="num">165</td><td class="num">0.1</td></tr>
      <tr><td>Full architecture sweep (1,000 loci × all axes)</td><td class="num">540,000</td><td class="num">2,971</td><td class="num">1.0</td></tr>
      <tr><td>GTEx real-data pass (49 tissues × 20k genes)</td><td class="num">980,000</td><td class="num">5,393</td><td class="num">1.8</td></tr>
      <tr><td>100 GWAS traits (~1,500 loci each)</td><td class="num">150,000</td><td class="num">825</td><td class="num">0.3</td></tr>
    </tbody>
  </table></div>
  <h3>Annualized request — transparent multipliers the PI can tune</h3>
  <div class="scroll"><table>
    <thead><tr><th>Activity</th><th>Assumption</th><th class="num">SU / yr</th></tr></thead>
    <tbody>
      <tr><td>Simulation method development</td><td>~40 full sweeps × ~3,000 SU</td><td class="num">120,000</td></tr>
      <tr><td>Real-data GTEx passes</td><td>~10 passes × ~5,400 SU × ~3 (realistic p, 3 methods)</td><td class="num">160,000</td></tr>
      <tr><td>GWAS-trait fine-mapping</td><td>~5 passes × ~825 SU × ~3</td><td class="num">12,000</td></tr>
      <tr><td>Replicate/CI + retries + safety margin</td><td>~30%</td><td class="num">90,000</td></tr>
      <tr class="total"><td colspan="2"><b>Total requested</b></td><td class="num"><b>≈ 380,000</b></td></tr>
    </tbody>
  </table></div>
  <p style="font-size:13.5px;color:var(--muted)">Every base number is measured on this hardware; only the campaign <em>counts</em>
  are policy choices to tune against the allocation ceiling.</p>
</section>

</div></main>

<footer><div class="wrap">
  <span>Reproducible:</span>
  <span class="pill mono">bench.py</span><span class="pill mono">parallel_bench.py</span>
  <span class="pill mono">fullgrid_bench.py</span><span class="pill mono">analyze.py</span>
  <span class="pill mono">results/*.json</span>
  <span style="color:var(--faint)">— env <b>dsc180</b>: GCTA 1.94.1 · susieR · PLINK 1.90 · Python 3 · SDSC 2× AMD EPYC 7713</span>
</div></footer>
"""

pathlib.Path("cis_heritability_expanse_report.html").write_text(HTML)
print("wrote cis_heritability_expanse_report.html", len(HTML), "bytes")
