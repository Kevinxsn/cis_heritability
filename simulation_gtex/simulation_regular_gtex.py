# This script extends the 1000G simulation pipeline to work with real GTEx data.
# Instead of simulating genotypes and phenotypes, it uses:
#   - Real genotypes from GTEx v8 (EUR, HapMap3 SNPs)
#   - Real gene expression from a chosen tissue (e.g., Whole Blood)
# It still compares GCTA vs SuSiE heritability estimates, but on real eQTL data.
#
# Key differences from the 1000G simulation:
#   1. Genes are picked from the expression BED file (real genes, not random loci)
#   2. Genotypes are extracted for the cis-window around each gene
#   3. Expression values are real (not simulated), so there's one h2 per gene (not a grid)
#   4. Samples are matched between genotype .fam and expression BED columns

import shutil
from pandas_plink import read_plink
import numpy as np
import scipy.linalg as linalg
import pandas as pd
import subprocess
import os
import argparse
import gzip
import statsmodels.api as sm

WINDOW = 500000  # cis-window: 500kb on each side of the gene TSS
SEP = "\t"
OUTPUT_COLUMNS = (
    f"gene{SEP}chr{SEP}tss{SEP}n_samples{SEP}n_snps{SEP}"
    f"gcta_hsq{SEP}gcta_hsq_se{SEP}gcta_hsq_p{SEP}"
    f"susie_hsq{SEP}susie_hsq_se\n"
)


# ── helpers (unchanged from simulation) ──────────────────────────────────────

def write_pheno(phen, out):
    """Write phenotype file for GCTA (FID, IID, pheno)."""
    with open(out + '.phen', 'w') as f:
        for i in range(len(phen)):
            f.write(f'{i}\t{i}\t{phen[i]}\n')


def write_grm(grm, nsnp, out):
    """Write GRM in GCTA format."""
    nindv = grm.shape[0]
    with gzip.open(out + '.grm.gz', 'w') as out_file, \
         open(out + '.grm.id', 'w') as out_id_file:
        for i in range(nindv):
            for j in range(0, i + 1):
                out_file.write(
                    f'{i+1}\t{j+1}\t{nsnp}\t{grm[i,j]}\n'.encode('utf-8')
                )
            out_id_file.write(f'{i}\t{i}\n')


def run_gcta(gcta, out):
    """Run GCTA REML."""
    command = (
        f'{gcta} --grm-gz {out} --pheno {out}.phen '
        f'--reml --reml-no-constrain --out {out}'
    )
    subprocess.run(command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def clean_up(out):
    for ext in ['.grm.gz', '.grm.id', '.phen']:
        path = out + ext
        if os.path.exists(path):
            os.remove(path)


def get_gcta_hsq(geno, gexpr, out, gcta):
    """
    Calculate heritability using GCTA-GREML.

    :param geno: numpy.ndarray n x p centered/scaled genotype matrix
    :param gexpr: numpy.ndarray phenotype vector (length n)
    :param out: output prefix
    :param gcta: path to gcta binary
    :returns: (hsq, hsq_se, hsq_p)
    """
    nindv, nsnp = geno.shape
    grm = np.dot(geno, geno.T) / nsnp
    write_pheno(gexpr, out)
    write_grm(grm, nsnp, out)
    run_gcta(gcta, out)
    clean_up(out)

    hsq_f = f'{out}.hsq'
    if os.path.exists(hsq_f):
        hsq_out = pd.read_table(hsq_f)
        hsq = hsq_out['Variance'][0]
        hsq_se = hsq_out['SE'][0]
        hsq_p = hsq_out['Variance'][8]
        os.remove(hsq_f)
        return hsq, hsq_se, hsq_p
    else:
        return 0, 10000, 1


def build_susie_sumstats(bim, geno, gexpr, sumstats_path, ld_path, gene_name,
                         gene_pref, plink_path, n_samples):
    """
    Build summary statistics and run SuSiE via R script.

    :param bim: DataFrame from read_plink
    :param geno: numpy.ndarray n x p centered/scaled genotype matrix
    :param gexpr: numpy.ndarray phenotype vector
    :param sumstats_path: directory for summary stats files
    :param ld_path: directory for LD files
    :param gene_name: identifier for the gene
    :param gene_pref: PLINK prefix for the extracted gene region
    :param plink_path: path to plink binary
    :param n_samples: number of samples
    :returns: filepath to SuSiE output
    """
    betas = []
    pvals = []
    std_errs = []
    for i in range(geno.shape[1]):
        x = sm.add_constant(geno[:, i].reshape(-1, 1))
        mod = sm.OLS(gexpr, x)
        results = mod.fit()
        betas.append(results.params[1])
        pvals.append(results.pvalues[1])
        std_errs.append(results.bse[1])

    sumstats = pd.DataFrame({
        'Gene': gene_name,
        'SNP': bim['snp'].values,
        'A1': bim['a0'].values,
        'A2': bim['a1'].values,
        'BETA': betas,
        'SE': std_errs,
        'P': pvals
    })
    filename = os.path.join(sumstats_path, gene_name + ".txt")
    sumstats.to_csv(filename, sep="\t", index=False, header=True)

    ld_path_gene = os.path.join(ld_path, gene_name)
    subprocess.run([
        plink_path, '--bfile', gene_pref, '--r', 'inter-chr',
        '--ld-window-r2', '0', '--keep-allele-order', '--extract',
        filename, '--out', ld_path_gene
    ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    subprocess.run([
        'Rscript', 'run_susie.R',
        '-c', filename,
        '-l', ld_path_gene + ".ld",
        '-n', str(n_samples),
        '-b', gene_pref + ".bim",
        '-o', filename
    ])
    return filename


def calculate_analytical_se(susie_output_file):
    """
    Calculate SE using the 4th Moment Formula.
    """
    try:
        df = pd.read_csv(susie_output_file, sep=" ")
        if 'POSTERIOR' not in df.columns or 'POSTERIOR_SD' not in df.columns:
            return np.nan

        mu = df['POSTERIOR']
        sigma = df['POSTERIOR_SD']
        E_b2 = mu**2 + sigma**2
        E_b4 = mu**4 + 6 * (mu**2) * (sigma**2) + 3 * (sigma**4)
        var_b2 = E_b4 - (E_b2**2)
        return np.sqrt(np.sum(var_b2))
    except Exception as e:
        print(f"Error in Analytical SE calc: {e}")
        return np.nan


def resolve_exe(exe_name_or_path):
    if os.path.sep in exe_name_or_path:
        return exe_name_or_path
    resolved = shutil.which(exe_name_or_path)
    if resolved is None:
        raise RuntimeError(
            f"Cannot find executable '{exe_name_or_path}' on PATH. "
            f"Install it or pass a full path."
        )
    return resolved


def create_directory(directory):
    os.makedirs(directory, exist_ok=True)
    print(f"Directory '{directory}' ready.")


# ── GTEx-specific functions ──────────────────────────────────────────────────

def load_expression_bed(expression_bed_path):
    """
    Load a GTEx normalized expression BED file.

    Format: #chr  start  end  gene_id  sample1  sample2  ...
    Returns: (gene_info_df, expression_matrix_df, sample_ids)
    """
    print(f"Loading expression data from {expression_bed_path}...")
    df = pd.read_csv(expression_bed_path, sep='\t', compression='gzip')

    # First 4 columns are: #chr, start, end, gene_id
    gene_info = df.iloc[:, :4].copy()
    gene_info.columns = ['chr', 'start', 'end', 'gene_id']

    # Remaining columns are sample expression values
    sample_ids = list(df.columns[4:])
    expr_matrix = df.iloc[:, 4:]

    print(f"  Loaded {len(gene_info)} genes, {len(sample_ids)} samples")
    return gene_info, expr_matrix, sample_ids


def match_samples(fam_path, expr_sample_ids):
    """
    Find overlapping samples between genotype .fam and expression data.

    GTEx .fam IIDs and expression column headers should use the same
    donor IDs (e.g., GTEX-XXXX).

    Returns: (common_ids, fam_indices, expr_indices)
    """
    fam = pd.read_csv(fam_path, sep=r'\s+', header=None,
                       names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])

    fam_ids = list(fam['IID'])
    common = [sid for sid in fam_ids if sid in set(expr_sample_ids)]

    fam_idx = [fam_ids.index(sid) for sid in common]
    expr_idx = [expr_sample_ids.index(sid) for sid in common]

    print(f"  Sample overlap: {len(common)} "
          f"(fam={len(fam_ids)}, expr={len(expr_sample_ids)})")
    return common, fam_idx, expr_idx


def extract_gene_region(plink_path, geno_prefix, gene_chr, tss, window,
                        output_prefix, keep_file=None):
    """
    Use PLINK to extract a cis-window around a gene's TSS from the
    combined genotype file.

    :param plink_path: path to plink binary
    :param geno_prefix: prefix for the combined genotype file (all chr)
    :param gene_chr: chromosome (e.g., "1", "chr1")
    :param tss: transcription start site position
    :param window: bp window on each side
    :param output_prefix: output PLINK prefix
    :param keep_file: optional file with FID/IID to keep (matched samples)
    :returns: True if successful (enough SNPs), False otherwise
    """
    chrom = str(gene_chr).replace('chr', '')
    start_pos = max(0, tss - window)
    end_pos = tss + window

    cmd = [
        plink_path,
        '--bfile', geno_prefix,
        '--chr', chrom,
        '--from-bp', str(start_pos),
        '--to-bp', str(end_pos),
        '--make-bed',
        '--out', output_prefix
    ]
    if keep_file:
        cmd.extend(['--keep', keep_file])

    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    bim_file = output_prefix + '.bim'
    if not os.path.exists(bim_file):
        return False

    n_snps = sum(1 for _ in open(bim_file))
    if n_snps < 50:
        # Clean up insufficient region
        for ext in ['.bed', '.bim', '.fam', '.log', '.nosex']:
            p = output_prefix + ext
            if os.path.exists(p):
                os.remove(p)
        return False

    return True


def load_real_genotypes(plink_prefix, sample_order=None):
    """
    Load genotypes from extracted PLINK files, center and scale.

    :param plink_prefix: PLINK prefix for the gene region
    :param sample_order: if given, reorder rows to match
    :returns: (bim_df, genotype_matrix n x p)
    """
    bim, fam, G = read_plink(plink_prefix, verbose=False)
    G = G.T  # n_samples x n_snps

    # Handle missing values: replace NaN with column mean
    G = G.compute()
    col_means = np.nanmean(G, axis=0)
    nan_mask = np.isnan(G)
    G[nan_mask] = np.take(col_means, np.where(nan_mask)[1])

    # Center and scale
    G -= np.mean(G, axis=0)
    std_devs = np.std(G, axis=0)
    std_devs[std_devs == 0] = 1
    G /= std_devs

    return bim, G


# ── Main ─────────────────────────────────────────────────────────────────────

def arg_parser():
    parser = argparse.ArgumentParser(
        prog='GTEx Heritability Analysis: SuSiE vs GCTA',
        description='''Run GCTA and SuSiE heritability estimation on real GTEx
        expression data. Uses real genotypes and real gene expression from a
        specified tissue.'''
    )
    parser.add_argument("--plink", default="plink",
                        type=str, help="Path to plink binary")
    parser.add_argument("--gcta", default="gcta64",
                        help="GCTA executable (default: gcta64)")
    parser.add_argument('--geno_prefix', required=True, type=str,
                        help="Prefix for GTEx genotype PLINK files (all chr combined)")
    parser.add_argument('--expression_bed', required=True, type=str,
                        help="Path to tissue expression BED.gz file "
                             "(e.g., Whole_Blood.v8.normalized_expression.bed.gz)")
    parser.add_argument('--covariate_file', default=None, type=str,
                        help="Path to covariate file (optional, for future use)")
    parser.add_argument('--tmpdir', required=True, type=str,
                        help="Path to temporary/working directory")
    parser.add_argument('--current_date', required=True, type=str,
                        help="Date tag for output file naming")
    parser.add_argument('--n_genes', default=100, type=int,
                        help="Number of random genes to process (default: 100)")
    parser.add_argument('--min_snps', default=50, type=int,
                        help="Minimum SNPs in cis-window to keep a gene (default: 50)")
    return parser.parse_args()


def main():
    args = arg_parser()
    plink_exe = resolve_exe(args.plink)
    gcta_exe = resolve_exe(args.gcta)

    # ── Setup directories ────────────────────────────────────────────────
    genes_folder = os.path.join(args.tmpdir, 'gene_regions')
    create_directory(genes_folder)

    sumstats_path = os.path.join(args.tmpdir, 'sumstats')
    create_directory(sumstats_path)

    ld_path_directory = os.path.join(args.tmpdir, 'ld')
    create_directory(ld_path_directory)

    gcta_path_directory = os.path.join(args.tmpdir, 'hsq')
    create_directory(gcta_path_directory)

    # ── Load expression data ─────────────────────────────────────────────
    gene_info, expr_matrix, expr_sample_ids = load_expression_bed(
        args.expression_bed
    )

    # ── Match samples ────────────────────────────────────────────────────
    fam_path = args.geno_prefix + '.fam'
    common_ids, fam_idx, expr_idx = match_samples(fam_path, expr_sample_ids)
    n_samples = len(common_ids)

    if n_samples < 50:
        raise RuntimeError(
            f"Only {n_samples} overlapping samples found. Check that the "
            f"genotype .fam IIDs match the expression BED column headers."
        )
    print(f"Using {n_samples} matched samples.")

    # Write a PLINK keep file for sample subsetting
    keep_file = os.path.join(args.tmpdir, 'keep_samples.txt')
    fam_df = pd.read_csv(fam_path, sep=r'\s+', header=None,
                          names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'])
    keep_df = fam_df.iloc[fam_idx][['FID', 'IID']]
    keep_df.to_csv(keep_file, sep='\t', index=False, header=False)

    # ── Pick random genes ────────────────────────────────────────────────
    n_total_genes = len(gene_info)
    gene_indices = np.random.choice(n_total_genes, size=min(args.n_genes * 3, n_total_genes),
                                     replace=False)
    np.random.shuffle(gene_indices)

    # ── Output file ──────────────────────────────────────────────────────
    tissue_name = os.path.basename(args.expression_bed).split('.')[0]
    output_txt = f"gtex_{tissue_name}_h2_results_{args.current_date}.txt"
    if not os.path.exists(output_txt) or os.path.getsize(output_txt) == 0:
        with open(output_txt, "a") as f:
            f.write(OUTPUT_COLUMNS)

    # ── Process genes ────────────────────────────────────────────────────
    genes_processed = 0
    for gene_idx in gene_indices:
        if genes_processed >= args.n_genes:
            break

        row = gene_info.iloc[gene_idx]
        gene_id = row['gene_id']
        gene_chr = str(row['chr']).replace('chr', '')
        gene_tss = int(row['start'])  # BED start = TSS for expression data
        gene_name = gene_id.split('.')[0]  # Remove Ensembl version suffix

        print(f"\n{'='*60}")
        print(f"Processing gene {genes_processed+1}/{args.n_genes}: "
              f"{gene_name} (chr{gene_chr}:{gene_tss})")

        # Extract cis-window genotypes
        gene_plink_prefix = os.path.join(genes_folder, gene_name)
        success = extract_gene_region(
            plink_exe, args.geno_prefix, gene_chr, gene_tss,
            WINDOW, gene_plink_prefix, keep_file
        )
        if not success:
            print(f"  Skipping {gene_name}: not enough SNPs in cis-window")
            continue

        # Load real genotypes
        try:
            bim, G = load_real_genotypes(gene_plink_prefix)
        except Exception as e:
            print(f"  Skipping {gene_name}: error loading genotypes: {e}")
            continue

        n_snps = G.shape[1]
        print(f"  Loaded {G.shape[0]} samples x {n_snps} SNPs")

        # Get real expression for matched samples
        expr_values = expr_matrix.iloc[gene_idx, expr_idx].values.astype(float)

        # The PLINK --keep may reorder samples, so we need to align
        # Read the .fam of the extracted region to get actual sample order
        extracted_fam = pd.read_csv(
            gene_plink_prefix + '.fam', sep=r'\s+', header=None,
            names=['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO']
        )
        extracted_ids = list(extracted_fam['IID'])

        # Reorder expression to match the extracted genotype sample order
        id_to_expr = dict(zip(common_ids, expr_values))
        expr_aligned = np.array([id_to_expr[sid] for sid in extracted_ids
                                  if sid in id_to_expr])

        if len(expr_aligned) != G.shape[0]:
            print(f"  Skipping {gene_name}: sample mismatch after alignment "
                  f"({len(expr_aligned)} vs {G.shape[0]})")
            continue

        # Standardize expression
        expr_aligned -= np.mean(expr_aligned)
        expr_std = np.std(expr_aligned)
        if expr_std == 0:
            print(f"  Skipping {gene_name}: zero variance in expression")
            continue
        expr_aligned /= expr_std

        # ── Run GCTA ─────────────────────────────────────────────────────
        max_retries = 3
        success = False

        for attempt in range(max_retries):
            try:
                gcta_prefix = os.path.join(gcta_path_directory, gene_name)
                h2g, hsq_se, hsq_p = get_gcta_hsq(
                    G, expr_aligned, gcta_prefix, gcta_exe
                )

                # ── Run SuSiE ───────────────────────────────────────────
                susie_result = build_susie_sumstats(
                    bim, G, expr_aligned,
                    sumstats_path, ld_path_directory,
                    gene_name, gene_plink_prefix,
                    plink_exe, len(expr_aligned)
                )

                df_susie = pd.read_csv(susie_result, sep=" ")
                if 'POSTERIOR' not in df_susie.columns or \
                   'POSTERIOR_SD' not in df_susie.columns:
                    raise KeyError("SuSiE output missing POSTERIOR columns.")

                est_h2 = np.sum(
                    df_susie['POSTERIOR']**2 + df_susie['POSTERIOR_SD']**2
                )
                analytical_se = calculate_analytical_se(susie_result)

                # Write result
                with open(output_txt, "a") as f:
                    f.write(
                        f"{gene_name}{SEP}{gene_chr}{SEP}{gene_tss}{SEP}"
                        f"{len(expr_aligned)}{SEP}{n_snps}{SEP}"
                        f"{h2g}{SEP}{hsq_se}{SEP}{hsq_p}{SEP}"
                        f"{est_h2}{SEP}{analytical_se}\n"
                    )

                success = True
                genes_processed += 1
                print(f"  GCTA h2={h2g:.4f} (SE={hsq_se:.4f}), "
                      f"SuSiE h2={est_h2:.4f} (SE={analytical_se:.4f})")
                break

            except KeyError as e:
                print(f"  Attempt {attempt+1}/{max_retries} failed for "
                      f"{gene_name}: {e}")
            except Exception as e:
                print(f"  Attempt {attempt+1}/{max_retries} failed for "
                      f"{gene_name}: {e}")

        if not success:
            print(f"  Skipping {gene_name} after {max_retries} failed attempts.")

    print(f"\nDone! Processed {genes_processed} genes. Results in {output_txt}")


if __name__ == "__main__":
    main()