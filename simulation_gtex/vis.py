# ── GTEx Heritability Results: GCTA vs SuSiE Visualization ───────────────────
#
# This script reads the output of gtex_simulations.py and produces 4 figures:
#   1. Scatter plot: GCTA h2 vs SuSiE h2 (direct agreement)
#   2. Bland-Altman plot: difference vs average (systematic bias)
#   3. Box plots stratified by n_snps bins (effect of genetic complexity)
#   4. Box plots stratified by GCTA h2 bins (behavior across h2 range)
#
# Usage:
#   python visualize_gtex_results.py --input gtex_Whole_Blood_h2_results_2026_04_22.txt
#
# Optional: if you have SuSiE per-gene sumstats, pass --sumstats_dir to
#           stratify by number of credible sets.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
import argparse
import os
import glob


# ── Style ────────────────────────────────────────────────────────────────────

GCTA_COLOR = "#4A90D9"
SUSIE_COLOR = "#5CB85C"
ACCENT_COLOR = "#D94A4A"
BG_COLOR = "#FAFAFA"
GRID_COLOR = "#E0E0E0"

def style_ax(ax, xlabel="", ylabel="", title=""):
    ax.set_facecolor(BG_COLOR)
    ax.grid(True, alpha=0.4, color=GRID_COLOR, linestyle='-', linewidth=0.5)
    ax.set_xlabel(xlabel, fontsize=12, fontweight='medium')
    ax.set_ylabel(ylabel, fontsize=12, fontweight='medium')
    ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=10)


# ── Plot 1: Scatter (Direct Agreement) ──────────────────────────────────────

def plot_scatter_agreement(df, output_dir):
    """GCTA h2 vs SuSiE h2 scatter with identity line and correlation."""
    fig, ax = plt.subplots(figsize=(8, 7))

    ax.scatter(df['gcta_hsq'], df['susie_hsq'],
               c=GCTA_COLOR, edgecolors='white', s=60, alpha=0.7,
               linewidths=0.5, zorder=3)

    # Identity line
    lims = [
        min(df['gcta_hsq'].min(), df['susie_hsq'].min()) - 0.05,
        max(df['gcta_hsq'].max(), df['susie_hsq'].max()) + 0.05
    ]
    ax.plot(lims, lims, '--', color=ACCENT_COLOR, linewidth=1.5,
            alpha=0.8, label='y = x', zorder=2)

    # Regression line
    mask = np.isfinite(df['gcta_hsq']) & np.isfinite(df['susie_hsq'])
    if mask.sum() > 2:
        slope, intercept, r, p, se = stats.linregress(
            df.loc[mask, 'gcta_hsq'], df.loc[mask, 'susie_hsq']
        )
        x_fit = np.linspace(lims[0], lims[1], 100)
        ax.plot(x_fit, slope * x_fit + intercept, '-', color='#FF8C00',
                linewidth=1.5, alpha=0.8,
                label=f'Fit: slope={slope:.2f}, r={r:.2f}', zorder=2)

    style_ax(ax,
             xlabel='GCTA Heritability Estimate',
             ylabel='SuSiE Heritability Estimate',
             title='GCTA vs SuSiE: Direct Agreement')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')
    ax.legend(fontsize=10, framealpha=0.9)

    n = len(df)
    ax.text(0.03, 0.97, f'n = {n} genes', transform=ax.transAxes,
            fontsize=10, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=GRID_COLOR, alpha=0.9))

    fig.tight_layout()
    path = os.path.join(output_dir, 'scatter_agreement.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close(fig)


# ── Plot 2: Bland-Altman ────────────────────────────────────────────────────

def plot_bland_altman(df, output_dir):
    """Difference (GCTA - SuSiE) vs Average: reveals systematic bias."""
    avg = (df['gcta_hsq'] + df['susie_hsq']) / 2
    diff = df['gcta_hsq'] - df['susie_hsq']

    mean_diff = diff.mean()
    std_diff = diff.std()
    upper_loa = mean_diff + 1.96 * std_diff
    lower_loa = mean_diff - 1.96 * std_diff

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.scatter(avg, diff, c=GCTA_COLOR, edgecolors='white', s=60,
               alpha=0.7, linewidths=0.5, zorder=3)

    # Mean difference line
    ax.axhline(mean_diff, color=ACCENT_COLOR, linestyle='-', linewidth=1.5,
               label=f'Mean diff = {mean_diff:.4f}', zorder=2)

    # Limits of agreement
    ax.axhline(upper_loa, color='gray', linestyle='--', linewidth=1,
               label=f'+1.96 SD = {upper_loa:.4f}', zorder=2)
    ax.axhline(lower_loa, color='gray', linestyle='--', linewidth=1,
               label=f'−1.96 SD = {lower_loa:.4f}', zorder=2)

    # Zero reference
    ax.axhline(0, color='black', linestyle=':', linewidth=0.8, alpha=0.5)

    # Shade limits of agreement
    xlims = ax.get_xlim()
    ax.fill_between(xlims, lower_loa, upper_loa, alpha=0.06,
                    color=GCTA_COLOR, zorder=1)
    ax.set_xlim(xlims)

    style_ax(ax,
             xlabel='Average of GCTA and SuSiE h²',
             ylabel='Difference (GCTA − SuSiE)',
             title='Bland-Altman Plot: Systematic Bias')
    ax.legend(fontsize=9, framealpha=0.9, loc='upper right')

    fig.tight_layout()
    path = os.path.join(output_dir, 'bland_altman.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close(fig)


# ── Plot 3: Box plots stratified by n_snps ──────────────────────────────────

def plot_boxplot_by_nsnps(df, output_dir):
    """Side-by-side box plots of GCTA vs SuSiE h2, binned by n_snps."""
    # Create bins
    bins = pd.qcut(df['n_snps'], q=4, duplicates='drop')
    df = df.copy()
    df['snp_bin'] = bins

    bin_labels = sorted(df['snp_bin'].unique(), key=lambda x: x.left)

    fig, ax = plt.subplots(figsize=(10, 6))
    _side_by_side_box(ax, df, 'snp_bin', bin_labels,
                      title='Heritability Estimates by SNP Count (quartiles)',
                      xlabel='Number of SNPs in cis-window')

    fig.tight_layout()
    path = os.path.join(output_dir, 'boxplot_by_nsnps.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close(fig)


# ── Plot 4: Box plots stratified by GCTA h2 bins ────────────────────────────

def plot_boxplot_by_h2(df, output_dir):
    """Side-by-side box plots binned by GCTA h2 estimate range."""
    df = df.copy()

    # Define meaningful h2 bins
    edges = [-np.inf, 0.0, 0.1, 0.2, 0.4, np.inf]
    labels_str = ['< 0', '0–0.1', '0.1–0.2', '0.2–0.4', '> 0.4']
    df['h2_bin'] = pd.cut(df['gcta_hsq'], bins=edges, labels=labels_str)

    # Drop empty bins
    nonempty = [l for l in labels_str if (df['h2_bin'] == l).sum() > 0]

    if len(nonempty) < 2:
        print("  Skipping h2-bin plot: not enough variation in GCTA h2.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    _side_by_side_box(ax, df, 'h2_bin', nonempty,
                      title='Heritability Estimates by GCTA h² Range',
                      xlabel='GCTA h² bin')

    fig.tight_layout()
    path = os.path.join(output_dir, 'boxplot_by_h2_bin.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close(fig)


# ── Shared box plot helper ───────────────────────────────────────────────────

def _side_by_side_box(ax, df, group_col, group_order, title="", xlabel=""):
    """
    Draw side-by-side GCTA vs SuSiE box plots, similar to the user's
    simulation plotting style but adapted for real-data groupings.
    """
    gcta_data = [df.loc[df[group_col] == g, 'gcta_hsq'].dropna().values
                 for g in group_order]
    susie_data = [df.loc[df[group_col] == g, 'susie_hsq'].dropna().values
                  for g in group_order]

    x = np.arange(len(group_order))
    offset = 0.18
    width = 0.32

    bp_gcta = ax.boxplot(gcta_data, positions=x - offset, widths=width,
                         patch_artist=True, showfliers=True,
                         flierprops=dict(marker='o', markersize=3, alpha=0.4))
    bp_susie = ax.boxplot(susie_data, positions=x + offset, widths=width,
                          patch_artist=True, showfliers=True,
                          flierprops=dict(marker='o', markersize=3, alpha=0.4))

    for b in bp_gcta['boxes']:
        b.set_facecolor(GCTA_COLOR)
        b.set_alpha(0.75)
    for b in bp_susie['boxes']:
        b.set_facecolor(SUSIE_COLOR)
        b.set_alpha(0.75)
    for element in ['whiskers', 'caps', 'medians']:
        for line in bp_gcta[element] + bp_susie[element]:
            line.set_color('#333333')
            line.set_linewidth(1)

    # Mean markers (red dashes, matching user's original style)
    for i in range(len(group_order)):
        gcta_mean = np.mean(gcta_data[i]) if len(gcta_data[i]) else np.nan
        susie_mean = np.mean(susie_data[i]) if len(susie_data[i]) else np.nan
        if not np.isnan(gcta_mean):
            ax.hlines(gcta_mean, x[i] - offset - width/2,
                      x[i] - offset + width/2,
                      colors=ACCENT_COLOR, linestyles='--', linewidth=2)
        if not np.isnan(susie_mean):
            ax.hlines(susie_mean, x[i] + offset - width/2,
                      x[i] + offset + width/2,
                      colors=ACCENT_COLOR, linestyles='--', linewidth=2)

    # Count annotations
    for i, g in enumerate(group_order):
        n = len(gcta_data[i])
        ax.text(x[i], ax.get_ylim()[0], f'n={n}', ha='center',
                va='bottom', fontsize=8, color='gray')

    ax.set_xticks(x)
    ax.set_xticklabels([str(g) for g in group_order], fontsize=10)

    ax.legend(handles=[
        Patch(facecolor=GCTA_COLOR, edgecolor='#333', label='GCTA', alpha=0.75),
        Patch(facecolor=SUSIE_COLOR, edgecolor='#333', label='SuSiE', alpha=0.75),
        plt.Line2D([], [], color=ACCENT_COLOR, linestyle='--', linewidth=2,
                   label='Mean'),
    ], loc='upper right', fontsize=10, framealpha=0.9)

    style_ax(ax, xlabel=xlabel, ylabel='Heritability Estimate', title=title)


# ── Plot 5 (bonus): SE comparison ───────────────────────────────────────────

def plot_se_comparison(df, output_dir):
    """Scatter of GCTA SE vs SuSiE analytical SE."""
    mask = np.isfinite(df['gcta_hsq_se']) & np.isfinite(df['susie_hsq_se'])
    d = df[mask]

    if len(d) < 5:
        print("  Skipping SE comparison: not enough valid SE values.")
        return

    fig, ax = plt.subplots(figsize=(7, 7))

    ax.scatter(d['gcta_hsq_se'], d['susie_hsq_se'],
               c=SUSIE_COLOR, edgecolors='white', s=60, alpha=0.7,
               linewidths=0.5, zorder=3)

    lims = [0, max(d['gcta_hsq_se'].max(), d['susie_hsq_se'].max()) * 1.1]
    ax.plot(lims, lims, '--', color=ACCENT_COLOR, linewidth=1.5, alpha=0.8,
            label='y = x')

    style_ax(ax,
             xlabel='GCTA SE',
             ylabel='SuSiE Analytical SE',
             title='Standard Error Comparison')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')
    ax.legend(fontsize=10, framealpha=0.9)

    fig.tight_layout()
    path = os.path.join(output_dir, 'se_comparison.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close(fig)


# ── Plot 6 (bonus): Summary statistics table ────────────────────────────────

def plot_summary_table(df, output_dir):
    """Print and save a summary statistics comparison."""
    summary = pd.DataFrame({
        'Statistic': ['Count', 'Mean', 'Median', 'Std Dev', 'Min', 'Max',
                      'Genes with h² > 0', 'Genes with h² > 0.2'],
        'GCTA h²': [
            len(df),
            f"{df['gcta_hsq'].mean():.4f}",
            f"{df['gcta_hsq'].median():.4f}",
            f"{df['gcta_hsq'].std():.4f}",
            f"{df['gcta_hsq'].min():.4f}",
            f"{df['gcta_hsq'].max():.4f}",
            (df['gcta_hsq'] > 0).sum(),
            (df['gcta_hsq'] > 0.2).sum(),
        ],
        'SuSiE h²': [
            len(df),
            f"{df['susie_hsq'].mean():.4f}",
            f"{df['susie_hsq'].median():.4f}",
            f"{df['susie_hsq'].std():.4f}",
            f"{df['susie_hsq'].min():.4f}",
            f"{df['susie_hsq'].max():.4f}",
            (df['susie_hsq'] > 0).sum(),
            (df['susie_hsq'] > 0.2).sum(),
        ],
    })

    # Correlation
    mask = np.isfinite(df['gcta_hsq']) & np.isfinite(df['susie_hsq'])
    if mask.sum() > 2:
        r, p = stats.pearsonr(df.loc[mask, 'gcta_hsq'],
                               df.loc[mask, 'susie_hsq'])
        corr_row = pd.DataFrame({
            'Statistic': [f'Pearson r (p={p:.2e})'],
            'GCTA h²': [f'{r:.4f}'],
            'SuSiE h²': [f'{r:.4f}'],
        })
        summary = pd.concat([summary, corr_row], ignore_index=True)

    print("\n" + "=" * 55)
    print("  Summary: GCTA vs SuSiE Heritability Estimates")
    print("=" * 55)
    print(summary.to_string(index=False))
    print("=" * 55 + "\n")

    path = os.path.join(output_dir, 'summary_stats.csv')
    summary.to_csv(path, index=False)
    print(f"  Saved: {path}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Visualize GTEx GCTA vs SuSiE heritability results'
    )
    parser.add_argument('--input', required=True, type=str,
                        help='Path to the GTEx h2 results file '
                             '(output of gtex_simulations.py)')
    parser.add_argument('--output_dir', default='figures', type=str,
                        help='Directory to save figures (default: figures/)')
    parser.add_argument('--sumstats_dir', default=None, type=str,
                        help='Optional: path to SuSiE sumstats directory '
                             '(for credible set stratification)')
    args = parser.parse_args()

    # Load data
    print(f"Loading results from {args.input}...")
    df = pd.read_csv(args.input, sep='\t')
    print(f"  Loaded {len(df)} genes")

    # Convert numeric columns
    for col in ['gcta_hsq', 'gcta_hsq_se', 'gcta_hsq_p',
                'susie_hsq', 'susie_hsq_se', 'n_snps', 'n_samples']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # Filter out failed runs (GCTA returns 0 with SE=10000 on failure)
    n_before = len(df)
    df = df[df['gcta_hsq_se'] < 9999].copy()
    if len(df) < n_before:
        print(f"  Removed {n_before - len(df)} failed GCTA runs")

    # If sumstats dir provided, count credible sets per gene
    if args.sumstats_dir and os.path.isdir(args.sumstats_dir):
        print("  Counting SuSiE credible sets from sumstats...")
        cs_counts = {}
        for gene_name in df['gene']:
            pattern = os.path.join(args.sumstats_dir, f"{gene_name}*.txt")
            files = glob.glob(pattern)
            if files:
                try:
                    sdf = pd.read_csv(files[0], sep=' ')
                    n_cs = sdf['CS'].nunique() - (1 if -1 in sdf['CS'].values else 0)
                    cs_counts[gene_name] = n_cs
                except Exception:
                    cs_counts[gene_name] = np.nan
            else:
                cs_counts[gene_name] = np.nan
        df['n_cs'] = df['gene'].map(cs_counts)

    if len(df) < 3:
        print("Not enough data points to plot. Need at least 3 genes.")
        return

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate all plots
    print("\nGenerating plots...")
    plot_summary_table(df, args.output_dir)
    plot_scatter_agreement(df, args.output_dir)
    plot_bland_altman(df, args.output_dir)
    plot_boxplot_by_nsnps(df, args.output_dir)
    plot_boxplot_by_h2(df, args.output_dir)
    plot_se_comparison(df, args.output_dir)

    print(f"\nAll figures saved to {args.output_dir}/")


if __name__ == "__main__":
    main()