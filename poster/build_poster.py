"""Build the IGVF 2026 poster (36 in wide x 48 in tall) as a .pptx.

The poster is described once, as a list of blocks per column, and rendered
twice: to PowerPoint (editable) and to a PNG proof.  Rendering both from the
same geometry means the proof cannot disagree with the deck, which matters
because there is no LibreOffice on this host to render the .pptx directly.

    python build_poster.py

Outputs
    XuS_IGVF2026_poster.pptx   the deliverable
    poster_proof.png           a flat proof for checking the layout

Every figure comes from figures/, built by make_figures.py straight off the
result tables; nothing here restates a number by hand except the callouts,
which are checked against the tables in check_numbers.py.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Inches, Pt

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"

# --- sheet geometry (inches) ------------------------------------------------
W, H = 36.0, 48.0
MARGIN = 0.85
GUTTER = 0.55
NCOL = 3
COLW = (W - 2 * MARGIN - (NCOL - 1) * GUTTER) / NCOL   # 10.87 in
HEADER_H = 5.45                 # navy band: title, authors, affiliations
CONCLUSION_Y = HEADER_H + 0.42  # one-sentence take-home, full width
CONCLUSION_H = 1.15
WORKFLOW_Y = CONCLUSION_Y + CONCLUSION_H + 0.34
WORKFLOW_W = W - 2 * MARGIN - 3.0
# Measured, not guessed: a wrong constant here silently overlaps the columns.
_WF = Image.open(HERE / "figures" / "P0_workflow.png")
WORKFLOW_H = WORKFLOW_W * _WF.size[1] / _WF.size[0]
BODY_TOP = WORKFLOW_Y + WORKFLOW_H + 0.72
BODY_BOTTOM = H - 0.60

COLX = [MARGIN + i * (COLW + GUTTER) for i in range(NCOL)]

# --- palette (matches the figures) -----------------------------------------
NAVY = "104281"
BLUE = "2a78d6"
ORANGE = "eb6834"
AQUA = "1baf7a"
INK = "0b0b0b"
INK2 = "52514e"
RULE = "d9d8d3"
TINT = "eef4fc"      # pale blue callout fill
WARM = "fdf1ea"      # pale orange callout fill

FONT = "Arial"       # present on every poster-printing machine; DejaVu is not

# --- type scale (pt) --------------------------------------------------------
# Raised throughout on review: a 36 x 48 in poster is read from ~1.5 m, so body
# text below ~24 pt is not legible from the aisle.
TITLE_PT = 72
AUTHOR_PT = 32
AFFIL_PT = 21
CONCLUSION_PT = 31
SECTION_PT = 33
BODY_PT = 25
CAPTION_PT = 21
TABLE_PT = 22
REF_PT = 20
QR_SIDE = 2.05
# TODO(Sunan): point both at the real public repo, then regenerate figures/qr.png
QR_URL = "https://github.com/"
QR_URL_TEXT = ""
TABLE_HDR_H = 0.66
TABLE_ROW_H = 0.62


def rgb(h: str) -> RGBColor:
    return RGBColor.from_string(h)


# ---------------------------------------------------------------------------
# Block model
# ---------------------------------------------------------------------------
@dataclass
class Block:
    kind: str                       # section | body | bullets | figure | table | callout
    text: str = ""
    items: Sequence = field(default_factory=list)
    path: Optional[str] = None
    caption: str = ""
    size: float = 25.0
    fill: Optional[str] = None
    width_frac: float = 1.0         # figures narrower than the column
    gap_after: float = 0.30         # inches of white space below

    height: float = 0.0             # filled in by measure()
    y: float = 0.0


# --- inline emphasis --------------------------------------------------------
# The review asked for restrained bold instead of ALL-CAPS emphasis, so body
# text carries **markers** and every renderer splits on them.
_BOLD_RE = re.compile(r"\*\*(.+?)\*\*", re.S)


def strip_markup(text: str) -> str:
    return _BOLD_RE.sub(r"\1", text)


def split_bold(text: str):
    """-> [(fragment, is_bold), ...] in order."""
    out, pos = [], 0
    for m in _BOLD_RE.finditer(text):
        if m.start() > pos:
            out.append((text[pos:m.start()], False))
        out.append((m.group(1), True))
        pos = m.end()
    if pos < len(text):
        out.append((text[pos:], False))
    return out or [(text, False)]


# --- text metrics -----------------------------------------------------------
# Arial at S pt averages ~0.52*S pt per glyph; a column is COLW inches wide.
def _wrapped_lines(text: str, size: float, width_in: float) -> int:
    text = strip_markup(text)
    chars_per_line = max(8, int(width_in * 72 / (0.52 * size)))
    lines = 0
    for para in text.split("\n"):
        lines += max(1, -(-len(para) // chars_per_line))
    return lines


def _text_height(text: str, size: float, width_in: float, pad: float = 0.20) -> float:
    return _wrapped_lines(text, size, width_in) * (size * 1.30) / 72 + pad


def measure(b: Block, colw: float) -> float:
    if b.kind == "section":
        return 1.12
    if b.kind == "body":
        return _text_height(b.text, b.size, colw - 0.15)
    if b.kind == "bullets":
        h = 0.0
        for it in b.items:
            h += _text_height(it, b.size, colw - 0.70) + 0.10
        return h
    if b.kind == "callout":
        inner = colw - 0.70
        h = 0.30
        for it in b.items:
            h += _text_height(it, b.size, inner, pad=0.10) + 0.06
        return h + 0.30
    if b.kind == "figure":
        im = Image.open(FIG / b.path)
        w = colw * b.width_frac
        h = w * im.size[1] / im.size[0]
        if b.caption:
            h += _text_height(b.caption, CAPTION_PT, colw - 0.15) + 0.10
        return h
    if b.kind == "table":
        # header row + body rows
        return TABLE_HDR_H + TABLE_ROW_H * (len(b.items) - 1)
    if b.kind == "qr":
        h = QR_SIDE
        if b.caption:
            h += _text_height(b.caption, REF_PT, colw - 0.15) + 0.12
        return h
    raise ValueError(b.kind)


# ---------------------------------------------------------------------------
# Poster content
# ---------------------------------------------------------------------------
TITLE_PLAIN = ("CHEWS improves cis-heritability precision over individual-level\n"
               "GCTA for sparse regulation in small samples")
AUTHORS = "Sunan Xu¹   ·   Stephen Golzari²   ·   Tiffany Amariuta-Bartell³⁴"
AFFIL = ("¹Bioinformatics and Systems Biology Program  ·  ²Shu Chien-Gene Lay Department of "
         "Bioengineering  ·  ³Department of Medicine  ·  ⁴Halıcıoğlu Data Science Institute"
         "   —   University of California San Diego, La Jolla, CA, USA")

# The take-home, readable from the aisle.  "Precision" is used here in the
# operational sense defined in section 1 and nowhere else on the poster.
CONCLUSION = ("Estimator performance depends jointly on **regulatory architecture, sample "
              "size and LD source**. CHEWS lowers error in small samples with concentrated "
              "effects, and underestimates h² as effects become diffuse.")


def build_blocks() -> List[Block]:
    return [
        # ------------------------------------------------------------------
        Block("section", "1  The question"),
        Block("body",
              "Cis-SNP heritability (cis-h²) quantifies how strongly a gene is genetically "
              "regulated. GCTA/GREML needs individual-level genotypes; consortium data are "
              "increasingly summary statistics only.",
              size=BODY_PT),
        Block("callout", items=[
            "**Accuracy** — closeness to truth; bias and MAE, in simulation only.   "
            "**Precision** — sampling variability; SD over repeated phenotypes at a fixed locus.   "
            "**Agreement** — similarity of values; offset, slope, RMSE.   "
            "**Reliability** — reproducibility; split-half correlation.",
        ], size=22, fill=TINT),

        # ------------------------------------------------------------------
        Block("section", "2  What CHEWS estimates"),
        Block("body",
              "Truth is the quadratic form **βᵀRβ**. SuSiE writes β = Σₗ βₗ, one "
              "single-effect vector per component, each selecting exactly one variant. "
              "Let **mₗ = E[βₗ] = αₗ⊙μₗ** and **sₗ = E[βₗᵀRβₗ]**. Then",
              size=BODY_PT),
        Block("callout", items=[
            "**E[βᵀRβ | data]  =  μᵀRμ  −  Σₗ mₗᵀR mₗ  +  Σₗ sₗ**",
        ], size=26, fill=WARM),
        Block("body",
              "The subtraction removes within-component cross products that cannot occur. "
              "The pipeline has always reported the **trace form Σⱼ(μⱼ²+σⱼ²)**, which is "
              "what this reduces to at R = I; the exact form is checked against "
              "brute-force tr(R·E[ββᵀ]) to 2e−16.",
              size=BODY_PT),
        Block("callout", items=[
            "**At n = 300 the two forms track each other on average** — r = 0.996, "
            "median |difference| 0.002, and the mean bias curves overlap at every causal "
            "count (§4). Per locus they can still differ: **95th percentile 0.034, "
            "maximum 0.368**. The exact form was evaluated at n = 300 only; §5–6 use the "
            "trace form.",
        ], size=22, fill=TINT),
        # ------------------------------------------------------------------
        Block("section", "3  Benchmark"),
        Block("table", items=[
            ["estimator", "data needed", "model assumed"],
            ["GCTA / REML", "individual genotypes", "infinitesimal"],
            ["**CHEWS** (here)", "summary stats + LD", "sparse causal set"],
            ["HEELS", "summary stats + LD", "infinitesimal"],
        ]),
        Block("bullets", items=[
            "350 random 1000 Genomes EUR cis-regions (±500 kb, median 8,200 SNPs). "
            "Effects rescaled so βᵀRβ = h², so only their **distribution** varies. Swept: "
            "n, h², causal count, effect prior, causal MAF, regional LD.",
            "**20,400 replicates attempted; 20,294 complete method triplets** "
            "(unconstrained GCTA/GREML failed in 106). Intervals resample whole regions.",
        ], size=23),

        # ------------------------------------------------------------------
        Block("section", "4  Architecture decides"),
        Block("figure", path="P2_win_landscape.png", caption=(
            "% of 2,990 paired fits at n = 300 where CHEWS is closer to truth; ~100 per "
            "cell. A win rate; §5 reports error magnitude. Median and mean absolute error "
            "agree, so rare bad fits drive neither.")),
        Block("figure", path="P1_bias_only.png", caption=(
            "Bias pooled over the h² grid, n = 300; 95% region-clustered bootstrap "
            "intervals (2,000 resamples). "
            "Dashed purple is the exact LD-aware estimator (§2); it lies on the trace "
            "form.")),
        Block("body",
              "One causal SNP: CHEWS nearly unbiased (−0.008), closer on 76% of loci. At "
              "twenty, **both** CHEWS forms underestimate h² by ≈0.13 while GCTA and "
              "HEELS stay near zero, so the LD cross terms do not explain it. "
              "**Consistent with the sparse-effect model** — but CHEWS and HEELS also "
              "differ in LD source (§5), so that attribution is provisional.",
              size=BODY_PT),

        # ------------------------------------------------------------------
        Block("section", "5  Sample size cuts both ways"),
        Block("figure", path="P3_sample_size.png", caption=(
            "(A) SD of ĥ² across loci; a within-locus design agrees to ≤0.007. (B) MAE "
            "difference by n and causal count. 1,464–1,500 replicates "
            "per point.")),
        Block("body",
              "With **25 independent phenotypes per locus** (30 regions), CHEWS's sampling "
              "SD is **0.108 against GCTA's 0.246** at n = 100, one causal SNP — 2.3×, "
              "and 2–4× across conditions. But the error advantage is **an interaction, "
              "not a rule about n**: largest at small n with concentrated effects, "
              "reversed at 20 causal SNPs at every sample size.",
              size=BODY_PT),
        Block("callout", items=[
            "HEELS is **not** a clean sparsity control as run: it received in-sample LD "
            "while CHEWS received reference-panel LD, so the two differ in model **and** "
            "LD source — the matched factorial is on the methods page. At n = 100 "
            "(p/n ≈ 70) it converges on 100% of fits but **46% of estimates exceed "
            "ĥ² = 0.9**: saturation, not non-convergence.",
        ], size=22, fill=WARM),

        # ------------------------------------------------------------------
        Block("section", "6  GTEx: agreement, not truth"),
        Block("body",
              "19,432 genes, 574 EUR donors, five nested subsets; covariates as REML fixed "
              "effects, in-sample LD; GCTA returned on 96.5–98.6%. **No ground truth "
              "here** — this measures disagreement, not accuracy.",
              size=BODY_PT),
        Block("callout", items=[
            "**Updated after abstract submission:** covariates now enter REML directly "
            "rather than pre-residualizing expression. The values here replace the "
            "submitted GTEx estimates; the degrees-of-freedom explanation is on the "
            "methods page.",
        ], size=22, fill=WARM),
        Block("figure", path="P5_agreement.png", caption=(
            "Agreement with GCTA, 18,747–19,169 complete-case genes per subset.")),
        Block("body",
              "At **n = 574, 19,169 complete-case genes**, HEELS tracks REML closely but "
              "**not exactly**: r = 0.9996, slope 0.93, offset −0.0074, RMSE 0.012. CHEWS "
              "is shrunk at n = 100 (slope 0.67), reaching unit slope by n = 574.",
              size=BODY_PT),

        # ------------------------------------------------------------------
        Block("section", "7  A conditioning artefact"),
        Block("figure", path="P6_strata_compare.png", caption=(
            "Mean CHEWS − GCTA within joint strata; 50–9,572 genes per cell, <30 blank.")),
        Block("body",
              "Using GCTA in both the stratum definition and the difference creates "
              "mathematical coupling. In (B) **both** stratifying variables come from "
              "half A and the difference from half B, sharing no donor — and most of the "
              "pattern goes with the coupling. **No Simpson's paradox is claimed.**",
              size=BODY_PT),

        # ------------------------------------------------------------------
        Block("section", "8  Robustness across donors"),
        Block("body",
              "Two disjoint halves of 287 donors, Pearson r on common genes:",
              size=BODY_PT),
        Block("table", items=[
            ["gene set", "GCTA", "CHEWS", "HEELS"],
            ["all genes (18.5–19.4k)", "0.865", "**0.912**", "0.870"],
            ["prespecified ≥1 eQTL (9.7k)", "0.850", "**0.898**", "0.857"],
            ["post-hoc ĥ²>0.05 (4.2k)", "0.848", "**0.892**", "0.855"],
        ]),
        Block("body",
              "**In this one donor split** CHEWS had the higher split-half correlation on "
              "all three gene sets; **repeated donor partitions are needed to quantify "
              "robustness to the split itself.** A paired bootstrap over genes gives "
              "Δr = +0.050, but resamples genes, not splits. Reliability is precision, "
              "not accuracy.",
              size=BODY_PT),

        # ------------------------------------------------------------------
        Block("section", "9  How to choose"),
        Block("table", items=[
            ["situation", "interpretation (current evidence)"],
            ["In-sample sufficient stats", "HEELS targets the same REML equations"],
            ["Reference-panel LD", "assess LD match; not isolated here"],
            ["Small n, concentrated", "CHEWS may lower MAE despite bias"],
            ["Diffuse signal", "CHEWS underestimates h²"],
            ["Strong MAF dependence", "stratified GRM, but validate first"],
        ]),

        # ------------------------------------------------------------------
        Block("section", "10  Limitations · References"),
        Block("bullets", items=[
            "Nominal 95% Wald intervals under-covered (68–87%). Screen fine-mapping h² "
            "against [0,1] — rare fits read 3–5 while reporting clean convergence.",
            "Simulated genotypes are multivariate normal: pairwise LD preserved, discrete "
            "counts and haplotypes not. One tissue, one donor split.",
        ], size=21),
        Block("body",
              "1. Yang J. GCTA. AJHG 2011;88:76.  2. Wang G. SuSiE. JRSS-B 2020;82:1273.  "
              "3. Zou Y. susie_rss. PLoS Genet 2022;18:e1010299.  4. Li H, Mazumder R, "
              "Lin X. HEELS. Nat Commun 2023;14:7954.  5. GTEx. Science 2020;369:1318.",
              size=REF_PT),
        # TODO(Sunan): grant numbers, and point the QR at the real repository
        # before printing (regenerate figures/qr.png in make_figures.py).
        Block("qr", text="Methods, full result tables,\ncode and sensitivity analyses",
              caption="Amariuta Lab, UC San Diego   ·   Contact: sux002@ucsd.edu"),
    ]


SECTION_LEAD = 0.26          # extra air above a section header
BREAK_PENALTY = 9.0
SAFETY = 0.7                 # inches of headroom kept at the foot of a column          # inches of "virtual imbalance" charged for
                             # cutting inside a section rather than between two


def _column_height(blocks: List[Block]) -> float:
    """Height a run of blocks occupies once laid out, ignoring the trailing gap.

    Must mirror the cursor arithmetic in pack() exactly, SECTION_LEAD included -
    one lead for every section header except one starting the column.
    """
    if not blocks:
        return 0.0
    h = sum(b.height + b.gap_after for b in blocks) - blocks[-1].gap_after
    leads = sum(1 for j, b in enumerate(blocks) if b.kind == "section" and j)
    return h + leads * SECTION_LEAD


def pack(blocks: List[Block]) -> List[List[Block]]:
    """Split the ordered block list into NCOL columns of similar height.

    Cuts are allowed between any two blocks EXCEPT immediately after a section
    header - that would orphan the header at the foot of a column.  Cutting
    inside a section is legal but charged BREAK_PENALTY, so the packer only
    does it when whole-section cuts cannot fit the sheet.
    """
    from itertools import combinations

    for b in blocks:
        b.height = measure(b, COLW)

    n = len(blocks)
    # Leave headroom: the text metric is an estimate of Arial's wrapping, so a
    # column that fits by 0.1 in on paper might not fit in PowerPoint.
    capacity = BODY_BOTTOM - BODY_TOP - SAFETY
    valid = [i for i in range(1, n) if blocks[i - 1].kind != "section"]

    best, best_key = None, None
    for cuts in combinations(valid, NCOL - 1):
        bounds = (0,) + cuts + (n,)
        spans = [blocks[bounds[i]:bounds[i + 1]] for i in range(NCOL)]
        if any(not sp for sp in spans):
            continue
        totals = [_column_height(sp) for sp in spans]
        penalty = sum(BREAK_PENALTY for c in cuts if blocks[c].kind != "section")
        key = (round(max(0.0, max(totals) - capacity), 3),
               max(totals) - min(totals) + penalty)
        if best_key is None or key < best_key:
            best_key, best = key, bounds

    if best is None:
        raise SystemExit("no valid column split")

    cols = []
    for i in range(NCOL):
        run = blocks[best[i]:best[i + 1]]
        y = BODY_TOP
        for j, b in enumerate(run):
            if b.kind == "section" and j:
                y += SECTION_LEAD
            b.y = y
            y += b.height + b.gap_after
        cols.append(run)
    return cols


# ---------------------------------------------------------------------------
# PowerPoint renderer
# ---------------------------------------------------------------------------
def add_text(slide, x, y, w, h, text, *, size, bold=False, color=INK,
             align=PP_ALIGN.LEFT, space_after=6, line=1.15, italic=False):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = Emu(0)
    tf.margin_top = tf.margin_bottom = Emu(0)
    for i, para in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = line
        p.space_after = Pt(space_after)
        for frag, is_bold in split_bold(para):
            r = p.add_run()
            r.text = frag
            r.font.size = Pt(size)
            r.font.bold = bold or is_bold
            r.font.italic = italic
            r.font.name = FONT
            r.font.color.rgb = rgb(color)
    return tb


def add_rect(slide, x, y, w, h, fill=None, line=None, line_w=1.0):
    from pptx.enum.shapes import MSO_SHAPE
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y),
                                Inches(w), Inches(h))
    if fill:
        sh.fill.solid()
        sh.fill.fore_color.rgb = rgb(fill)
    else:
        sh.fill.background()
    if line:
        sh.line.color.rgb = rgb(line)
        sh.line.width = Pt(line_w)
    else:
        sh.line.fill.background()
    sh.shadow.inherit = False
    return sh


def render_pptx(cols, out: Path):
    prs = Presentation()
    prs.slide_width = Inches(W)
    prs.slide_height = Inches(H)
    slide = prs.slides.add_slide(prs.slide_layouts[6])

    # --- header band --------------------------------------------------------
    add_rect(slide, 0, 0, W, HEADER_H, fill=NAVY)
    add_rect(slide, 0, HEADER_H, W, 0.18, fill=ORANGE)
    add_text(slide, MARGIN, 0.42, W - 2 * MARGIN, 2.6, TITLE_PLAIN,
             size=TITLE_PT, bold=True, color="ffffff", align=PP_ALIGN.CENTER, line=1.07)
    add_text(slide, MARGIN, 3.42, W - 2 * MARGIN, 0.7, AUTHORS,
             size=AUTHOR_PT, bold=True, color="dce9f8", align=PP_ALIGN.CENTER)
    add_text(slide, MARGIN, 4.28, W - 2 * MARGIN, 0.9, AFFIL,
             size=AFFIL_PT, color="b9cfeb", align=PP_ALIGN.CENTER, line=1.18)

    # --- one-sentence take-home, full width ---------------------------------
    add_rect(slide, MARGIN, CONCLUSION_Y, W - 2 * MARGIN, CONCLUSION_H, fill=TINT)
    add_rect(slide, MARGIN, CONCLUSION_Y, 0.16, CONCLUSION_H, fill=ORANGE)
    add_text(slide, MARGIN + 0.55, CONCLUSION_Y + 0.22, W - 2 * MARGIN - 1.1,
             CONCLUSION_H - 0.3, CONCLUSION, size=CONCLUSION_PT, bold=False,
             color=INK, align=PP_ALIGN.CENTER, line=1.16)

    # --- what CHEWS computes, as a strip ------------------------------------
    wf = FIG / "P0_workflow.png"
    im = Image.open(wf)
    ww, wh = WORKFLOW_W, WORKFLOW_H
    slide.shapes.add_picture(str(wf), Inches((W - ww) / 2), Inches(WORKFLOW_Y),
                             width=Inches(ww))
    add_text(slide, MARGIN, WORKFLOW_Y + wh + 0.02, W - 2 * MARGIN, 0.5,
             "CHEWS  =  Cis-Heritability Estimation With SuSiE", size=CAPTION_PT,
             color=INK2, align=PP_ALIGN.CENTER, italic=True)

    # --- columns ------------------------------------------------------------
    for ci, blocks in enumerate(cols):
        x = COLX[ci]
        for b in blocks:
            if b.kind == "section":
                add_rect(slide, x, b.y, COLW, 0.92, fill=NAVY)
                add_text(slide, x + 0.28, b.y + 0.17, COLW - 0.48, 0.7, b.text,
                         size=SECTION_PT, bold=True, color="ffffff")

            elif b.kind == "body":
                add_text(slide, x, b.y, COLW, b.height, b.text,
                         size=b.size, color=INK, line=1.22)

            elif b.kind == "bullets":
                yy = b.y
                for it in b.items:
                    hh = _text_height(it, b.size, COLW - 0.70)
                    add_text(slide, x + 0.10, yy, 0.45, 0.5, "▪",
                             size=b.size, color=ORANGE, bold=True)
                    add_text(slide, x + 0.55, yy, COLW - 0.55, hh, it,
                             size=b.size, color=INK, line=1.22)
                    yy += hh + 0.10

            elif b.kind == "callout":
                add_rect(slide, x, b.y, COLW, b.height, fill=b.fill)
                add_rect(slide, x, b.y, 0.13, b.height,
                         fill=ORANGE if b.fill == WARM else BLUE)
                yy = b.y + 0.30
                for it in b.items:
                    hh = _text_height(it, b.size, COLW - 0.70, pad=0.10)
                    add_text(slide, x + 0.42, yy, COLW - 0.72, hh, it,
                             size=b.size, color=INK, line=1.22)
                    yy += hh + 0.06

            elif b.kind == "figure":
                im = Image.open(FIG / b.path)
                fw = COLW * b.width_frac
                fh = fw * im.size[1] / im.size[0]
                fx = x + (COLW - fw) / 2
                slide.shapes.add_picture(str(FIG / b.path), Inches(fx), Inches(b.y),
                                         width=Inches(fw))
                if b.caption:
                    add_text(slide, x, b.y + fh + 0.10, COLW,
                             b.height - fh - 0.10, b.caption,
                             size=CAPTION_PT, color=INK2, line=1.18, italic=True)

            elif b.kind == "qr":
                qr = FIG / "qr.png"
                slide.shapes.add_picture(str(qr), Inches(x), Inches(b.y),
                                         width=Inches(QR_SIDE))
                add_text(slide, x + QR_SIDE + 0.40, b.y + 0.22,
                         COLW - QR_SIDE - 0.45, QR_SIDE - 0.3, b.text,
                         size=23, bold=True, color=INK, line=1.20)
                if QR_URL_TEXT:
                    add_text(slide, x + QR_SIDE + 0.40, b.y + QR_SIDE - 0.72,
                             COLW - QR_SIDE - 0.45, 0.5, QR_URL_TEXT,
                             size=REF_PT, color=BLUE, line=1.15)
                if b.caption:
                    add_text(slide, x, b.y + QR_SIDE + 0.12, COLW,
                             b.height - QR_SIDE - 0.12, b.caption,
                             size=REF_PT, color=INK2, line=1.18)

            elif b.kind == "table":
                rows, ncol = len(b.items), len(b.items[0])
                shp = slide.shapes.add_table(rows, ncol, Inches(x), Inches(b.y),
                                             Inches(COLW), Inches(b.height))
                tbl = shp.table
                tbl.first_row = True
                widths = ({2: [0.44, 0.56], 3: [0.31, 0.30, 0.39],
                           4: [0.52, 0.16, 0.16, 0.16]}[ncol])
                for j, frac in enumerate(widths):
                    tbl.columns[j].width = Inches(COLW * frac)
                tbl.rows[0].height = Inches(TABLE_HDR_H)
                for i in range(1, rows):
                    tbl.rows[i].height = Inches(TABLE_ROW_H)
                for i, row in enumerate(b.items):
                    for j, val in enumerate(row):
                        cell = tbl.cell(i, j)
                        cell.text = ""
                        cell.vertical_anchor = MSO_ANCHOR.MIDDLE
                        cell.margin_left = Inches(0.14)
                        cell.margin_right = Inches(0.08)
                        cell.margin_top = cell.margin_bottom = Inches(0.04)
                        cell.fill.solid()
                        cell.fill.fore_color.rgb = rgb(
                            NAVY if i == 0 else ("f6f5f2" if i % 2 else "ffffff"))
                        p = cell.text_frame.paragraphs[0]
                        for frag, is_bold in split_bold(val):
                            r = p.add_run()
                            r.text = frag
                            r.font.size = Pt(TABLE_PT)
                            r.font.bold = (i == 0) or (j == 0) or is_bold
                            r.font.name = FONT
                            r.font.color.rgb = rgb("ffffff" if i == 0 else INK)

    prs.save(out)
    print(f"wrote {out.name}")


# ---------------------------------------------------------------------------
# PNG proof - same geometry, so it cannot disagree with the deck
# ---------------------------------------------------------------------------
def render_proof(cols, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import matplotlib.image as mpimg

    PROOF_W = 16.0                      # proof sheet width, inches
    k = PROOF_W / W                     # poster inch -> proof inch
    fig = plt.figure(figsize=(PROOF_W, H * k), dpi=110)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), W, H, color="white"))

    def wrap(text, size, width_in):
        text = strip_markup(text)
        cpl = max(8, int(width_in * 72 / (0.52 * size)))
        out_lines = []
        for para in text.split("\n"):
            words, line = para.split(" "), ""
            for w_ in words:
                if len(line) + len(w_) + 1 <= cpl:
                    line = (line + " " + w_).strip()
                else:
                    out_lines.append(line)
                    line = w_
            out_lines.append(line)
        return "\n".join(out_lines)

    def txt(x, y, s, size, color="#0b0b0b", bold=False, w=None, ha="left", it=False):
        ax.text(x, y, wrap(s, size, w) if w else strip_markup(s), fontsize=size * k,
                color=color, ha=ha, va="top", family="DejaVu Sans",
                fontweight="bold" if bold else "normal",
                style="italic" if it else "normal", linespacing=1.22)

    ax.add_patch(Rectangle((0, 0), W, HEADER_H, color="#" + NAVY))
    ax.add_patch(Rectangle((0, HEADER_H), W, 0.18, color="#" + ORANGE))
    txt(W / 2, 0.42, TITLE_PLAIN, TITLE_PT, "white", True, ha="center")
    txt(W / 2, 3.42, AUTHORS, AUTHOR_PT, "#dce9f8", True, ha="center")
    txt(W / 2, 4.28, AFFIL, AFFIL_PT, "#b9cfeb", w=W - 2 * MARGIN, ha="center")

    ax.add_patch(Rectangle((MARGIN, CONCLUSION_Y), W - 2 * MARGIN, CONCLUSION_H,
                           color="#" + TINT))
    ax.add_patch(Rectangle((MARGIN, CONCLUSION_Y), 0.16, CONCLUSION_H,
                           color="#" + ORANGE))
    txt(W / 2, CONCLUSION_Y + 0.22, CONCLUSION, CONCLUSION_PT, "#0b0b0b",
        w=W - 2 * MARGIN - 1.1, ha="center")

    _wf = mpimg.imread(FIG / "P0_workflow.png")
    _ww, _wh = WORKFLOW_W, WORKFLOW_H
    ax.imshow(_wf, extent=((W - _ww) / 2, (W + _ww) / 2,
                           WORKFLOW_Y + _wh, WORKFLOW_Y), aspect="auto", zorder=5)
    txt(W / 2, WORKFLOW_Y + _wh + 0.02,
        "CHEWS  =  Cis-Heritability Estimation With SuSiE", CAPTION_PT,
        "#52514e", ha="center", it=True)

    for ci, blocks in enumerate(cols):
        x = COLX[ci]
        for b in blocks:
            if b.kind == "section":
                ax.add_patch(Rectangle((x, b.y), COLW, 0.92, color="#" + NAVY))
                txt(x + 0.28, b.y + 0.19, b.text, SECTION_PT, "white", True)
            elif b.kind == "body":
                txt(x, b.y, b.text, b.size, w=COLW - 0.15)
            elif b.kind == "bullets":
                yy = b.y
                for it in b.items:
                    hh = _text_height(it, b.size, COLW - 0.70)
                    txt(x + 0.08, yy, "▪", b.size, "#" + ORANGE, True)
                    txt(x + 0.55, yy, it, b.size, w=COLW - 0.55)
                    yy += hh + 0.10
            elif b.kind == "callout":
                ax.add_patch(Rectangle((x, b.y), COLW, b.height, color="#" + b.fill))
                ax.add_patch(Rectangle((x, b.y), 0.13, b.height,
                                       color="#" + (ORANGE if b.fill == WARM else BLUE)))
                yy = b.y + 0.30
                for it in b.items:
                    hh = _text_height(it, b.size, COLW - 0.70, pad=0.10)
                    txt(x + 0.42, yy, it, b.size, w=COLW - 0.72)
                    yy += hh + 0.06
            elif b.kind == "figure":
                im = mpimg.imread(FIG / b.path)
                fw = COLW * b.width_frac
                fh = fw * im.shape[0] / im.shape[1]
                fx = x + (COLW - fw) / 2
                ax.imshow(im, extent=(fx, fx + fw, b.y + fh, b.y), aspect="auto",
                          zorder=5)
                if b.caption:
                    txt(x, b.y + fh + 0.10, b.caption, CAPTION_PT, "#52514e",
                        w=COLW - 0.15, it=True)
            elif b.kind == "qr":
                q = mpimg.imread(FIG / "qr.png")
                ax.imshow(q, extent=(x, x + QR_SIDE, b.y + QR_SIDE, b.y),
                          aspect="auto", zorder=5)
                txt(x + QR_SIDE + 0.40, b.y + 0.22, b.text, 23, "#0b0b0b",
                    True, w=COLW - QR_SIDE - 0.45)
                if QR_URL_TEXT:
                    txt(x + QR_SIDE + 0.40, b.y + QR_SIDE - 0.72, QR_URL_TEXT, REF_PT,
                        "#" + BLUE, w=COLW - QR_SIDE - 0.45)
                if b.caption:
                    txt(x, b.y + QR_SIDE + 0.12, b.caption, REF_PT, "#52514e",
                        w=COLW - 0.15)

            elif b.kind == "table":
                rows = len(b.items)
                widths = ({2: [0.44, 0.56], 3: [0.31, 0.30, 0.39],
                           4: [0.52, 0.16, 0.16, 0.16]}[len(b.items[0])])
                yy = b.y
                for i, row in enumerate(b.items):
                    rh = TABLE_HDR_H if i == 0 else TABLE_ROW_H
                    fillc = "#" + NAVY if i == 0 else ("#f6f5f2" if i % 2 else "#ffffff")
                    ax.add_patch(Rectangle((x, yy), COLW, rh, color=fillc))
                    xx = x
                    for j, val in enumerate(row):
                        txt(xx + 0.14, yy + 0.19, val, TABLE_PT,
                            "white" if i == 0 else "#0b0b0b",
                            bold=(i == 0 or j == 0))
                        xx += COLW * widths[j]
                    yy += rh
                ax.add_patch(Rectangle((x, b.y), COLW, b.height, fill=False,
                                       edgecolor="#" + RULE, lw=1))

    fig.savefig(out, dpi=100)
    plt.close(fig)
    print(f"wrote {out.name}")


_TEX = re.compile(r"\$[^$]*\$")


def assert_no_tex(cols):
    """PowerPoint renders no TeX; $h^2$ would print literally. Catch it here."""
    bad = []
    for label, text in [("title", TITLE_PLAIN), ("authors", AUTHORS),
                        ("affil", AFFIL), ("conclusion", CONCLUSION)]:
        if _TEX.search(text):
            bad.append((label, text[:70]))
    for ci, blocks in enumerate(cols):
        for b in blocks:
            parts = [b.text, b.caption] + list(b.items or [])
            for part in parts:
                if isinstance(part, (list, tuple)):
                    parts.extend(part)
                    continue
                if isinstance(part, str) and _TEX.search(part):
                    bad.append((f"col{ci+1}/{b.kind}", part[:70]))
    if bad:
        for where, what in bad:
            print(f"  TeX markup in {where}: {what}")
        raise SystemExit("LaTeX markup would print literally in PowerPoint")


def warn_placeholders(cols):
    """A QR that resolves nowhere is worse than no QR - say so loudly."""
    if QR_URL.rstrip("/") in ("https://github.com", "http://github.com", ""):
        print("\n  !! QR code still points at a placeholder (" + QR_URL + ").")
        print("     Set QR_URL and QR_URL_TEXT, regenerate figures/qr.png, and")
        print("     re-run before printing. The URL line is hidden until then.\n")


def report(cols):
    print(f"\nsheet {W:.0f} x {H:.0f} in   column width {COLW:.2f} in   "
          f"body {BODY_TOP:.2f}–{BODY_BOTTOM:.2f} in")
    ok = True
    for ci, blocks in enumerate(cols):
        end = blocks[-1].y + blocks[-1].height
        slack = BODY_BOTTOM - end
        flag = "OK " if slack >= 0 else "OVERFLOW"
        if slack < 0:
            ok = False
        print(f"  column {ci+1}: {len(blocks):2d} blocks, ends at {end:5.2f} in, "
              f"slack {slack:+5.2f} in   [{flag}]")
    return ok


# ===========================================================================
# Round 4 (PI feedback) - docs/pi_feedback_round4_design.md
#
# Built ON TOP of the round-3 content: the round-3 block list is taken as is
# and only the §5 and §7 blocks are swapped.  `python build_poster.py` still
# writes the round-3 deck, byte for byte; `--round 4` writes *_round4 files.
# ===========================================================================
ROUND4_SEC5_FIG = Block("figure", path="P3_sample_size_r4.png", caption=(
    "(A) SD of ĥ² across loci; a within-locus design agrees to ≤0.007. (B) h² "
    "0.05–0.8 pooled, ~300 paired fits per cell; dotted: 95% CI spans 50%."))
ROUND4_SEC5_BODY = Block("body",
    "With **25 independent phenotypes per locus** (30 regions), CHEWS's sampling "
    "SD is **0.108 against GCTA's 0.246** at n = 100, one causal SNP — 2.3×, "
    "and 2–4× across conditions. With one causal SNP CHEWS is closer to truth on "
    "**75–82% of loci at every sample size**; what shrinks with n is the margin "
    "(MAE −73% at n = 100, −58% at n = 1,000). At twenty causal SNPs GCTA is "
    "closer on 72–77% of loci for n ≥ 300.",
    size=BODY_PT)

# §7 and the extra limitation come from round4_text.py, which reads part_g's
# pre-registered calibration verdict
ROUND4_SEC7: List[Block] = []
ROUND4_LIMITATIONS: List[str] = []


def build_blocks_round4() -> List[Block]:
    blocks = build_blocks()
    out: List[Block] = []
    i = 0
    while i < len(blocks):
        b = blocks[i]
        if b.kind == "figure" and b.path == "P3_sample_size.png":
            out.append(ROUND4_SEC5_FIG)
        elif b.kind == "body" and b.text.startswith("With **25 independent phenotypes"):
            out.append(ROUND4_SEC5_BODY)
        elif (b.kind == "bullets" and b.items and b.items[0].startswith("Nominal 95% Wald")
              and ROUND4_LIMITATIONS):
            out.append(Block("bullets", items=list(b.items) + ROUND4_LIMITATIONS, size=b.size))
        elif b.kind == "section" and b.text.startswith("7  ") and ROUND4_SEC7:
            out.extend(ROUND4_SEC7)
            # drop the round-3 §7 figure + paragraph: everything up to §8
            i += 1
            while not (blocks[i].kind == "section" and blocks[i].text.startswith("8  ")):
                i += 1
            continue
        else:
            out.append(b)
        i += 1
    return out


def word_count(blocks: List[Block]) -> int:
    """Words a reader meets in the columns: body, bullets, callouts, captions,
    table cells.  Section headers and the reference list are not counted."""
    n = 0
    for b in blocks:
        if b.kind == "body" and b.size != REF_PT:
            n += len(strip_markup(b.text).split())
        elif b.kind in ("bullets", "callout"):
            n += sum(len(strip_markup(t).split()) for t in b.items)
        elif b.kind == "table":
            n += sum(len(strip_markup(c).split()) for row in b.items for c in row)
        if b.caption and b.kind == "figure":
            n += len(b.caption.split())
    return n


if __name__ == "__main__":
    import sys

    round4 = "--round" in sys.argv and sys.argv[sys.argv.index("--round") + 1] == "4"
    if round4:
        import round4_text                      # fills ROUND4_SEC7 from the tables
        ROUND4_SEC7[:] = round4_text.sec7_blocks(Block, BODY_PT)
        ROUND4_LIMITATIONS[:] = [round4_text.LIMITATION]
        blocks = build_blocks_round4()
        print(f"round 4: {word_count(blocks)} words in the columns "
              f"(round 3: {word_count(build_blocks())})")
        deck, proof = "XuS_IGVF2026_poster_round4.pptx", "poster_proof_round4.png"
    else:
        blocks = build_blocks()
        deck, proof = "XuS_IGVF2026_poster.pptx", "poster_proof.png"
    cols = pack(blocks)
    assert_no_tex(cols)
    warn_placeholders(cols)
    good = report(cols)
    render_pptx(cols, HERE / deck)
    render_proof(cols, HERE / proof)
    if not good:
        raise SystemExit("layout overflows the sheet - shorten a column")
