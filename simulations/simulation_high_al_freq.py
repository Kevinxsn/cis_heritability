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
from sklearn import linear_model as lm
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold

NUM_PEOPLE = [300]
WINDOW = 500000
ncausal_arr = [1, 3, 5, 10] # number of causal SNPs for eQTL file
h2_arr = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8] # True heritability value
# New parameters for MAF and Effect Size architecture
maf_bins = {'rare': (0.01, 0.05), 'common': (0.05, 0.5)}
#effect_dists = ['normal', 'fixed', 'pareto']

effect_dists = ['normal']

SEP = "\t"
OUTPUT_COLUMNS = f"gene{SEP}ncausal{SEP}h2{SEP}maf_bin{SEP}effect_dist{SEP}hsq{SEP}hsq_se{SEP}hsq_p{SEP}su_hsq{SEP}su_hsq_se\n"


def ld_frac_r2_ge_quick(G_std, r2_threshold=0.8, max_snps=200, rng=None):
    """
    Quick LD summary: fraction of off-diagonal SNP pairs with r^2 >= r2_threshold,
    computed on a SNP subset (max_snps) for speed.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_samples, n_snps = G_std.shape
    m = min(max_snps, n_snps)
    if m < 2:
        return 0.0

    if n_snps > m:
        idx = rng.choice(n_snps, size=m, replace=False)
        Gs = G_std[:, idx]
    else:
        Gs = G_std

    LD_sub = (Gs.T @ Gs) / n_samples
    if hasattr(LD_sub, "compute"):
        LD_sub = LD_sub.compute()

    LD_sub = np.asarray(LD_sub)
    iu = np.triu_indices(LD_sub.shape[0], k=1)
    r2 = LD_sub[iu] ** 2
    if r2.size == 0:
        return 0.0
    return float(np.mean(r2 >= r2_threshold))


def get_ld(prefix):
    # return cholesky L
    bim, fam, G = read_plink(prefix, verbose=False)
    G = G.T # Transpose to have SNPs as columns

    # estimate LD for population from PLINK data
    n_samples, n_snps = G.shape
    mafs = (np.mean(G, axis=0) / 2).compute()

    # Center and standardize G
    G -= mafs * 2
    std_devs = np.std(G, axis=0)
    std_devs[std_devs == 0] = 1
    G /= std_devs

    # Regularize the LD matrix to make it Positive Semi Definite
    LD = np.dot(G.T, G) / n_samples
    regularization = 0.1
    LD += np.eye(n_snps) * regularization

    # re-adjust to get proper correlation matrix
    LD = LD / (1 + regularization)

    # compute cholesky decomp for faster sampling/simulation
    L = linalg.cholesky(LD, lower=True)

    return (bim, L)

def sim_geno(L, n): 
    """
    Sample genotypes from an MVN approximation.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param n: int the number of genotypes to sample

    :return: numpy.ndarray n x p centered/scaled genotype matrix
    """
    p, p = L.shape
    Z = L.dot(np.random.normal(size=(n, p)).T).T #pxn
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0) 
    return Z


def sim_beta(L, ncausal, eqtl_h2, mafs=None, maf_range=(0.0, 0.5), effect_dist='normal', rescale=True):
    n_snps = L.shape[0]
    if eqtl_h2 > 0 and ncausal > 0:
        
        # 1. Filter valid SNPs by MAF range
        if mafs is not None:
            valid_snps = np.where((mafs >= maf_range[0]) & (mafs <= maf_range[1]))[0]
            if len(valid_snps) < ncausal:
                # If a region doesn't have enough rare/common variants, expand the pool
                # print(f"Not enough SNPs in MAF range {maf_range}. Falling back to all SNPs.")
                valid_snps = np.arange(n_snps)
        else:
            valid_snps = np.arange(n_snps)

        # Pick causal SNPs from the filtered list (replace=False ensures we don't pick the same SNP twice)
        c_qtls = np.random.choice(valid_snps, ncausal, replace=False)
        b_qtls = np.zeros(n_snps)

        # 2. Apply different effect size architectures
        if effect_dist == 'normal':
            effects = np.random.normal(loc=0, scale=np.sqrt(eqtl_h2 / ncausal), size=ncausal)
        
        elif effect_dist == 'fixed':
            # Every causal SNP has identical impact
            signs = np.random.choice([-1, 1], ncausal)
            effects = signs * np.sqrt(eqtl_h2 / ncausal)
            
        elif effect_dist == 'pareto':
            # Heavy tailed: 1 or 2 massive effects, the rest are small
            effects = np.random.pareto(a=2.0, size=ncausal)
            signs = np.random.choice([-1, 1], ncausal)
            effects = signs * effects

        b_qtls[c_qtls] = effects
        
        # Rescale ensures that regardless of the architecture, the true h2 remains exactly what we want
        if rescale:
            s2g = compute_s2g(L, b_qtls)
            if s2g > 0:
                b_qtls *= np.sqrt(eqtl_h2 / s2g)
    else:
        b_qtls = np.zeros(n_snps)
        
    return b_qtls


def sim_trait(g, h2g): 
    """
    Simulate a complex trait as a function of latent genetic value and env noise.

    :param g: numpy.ndarray of latent genetic values
    :param h2g: float the heritability of the trait in the population

    :return: (numpy.ndarray, float) simulated phenotype, sd of Y
    """
    n = len(g)
    if h2g > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ( (1.0 / h2g ) - 1 )
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e #appears to add noise to the genetic values, adding the same about of noise if the seed is the same 
    else:
        e = np.random.normal(0, 1, n)
        y = e
    # standardize
    y -= np.mean(y)
    y_std = np.std(y)
    y /= y_std
    y_std = y_std.item()
    return y, y_std  


def compute_s2g(L, beta):
    """
    Compute genetic variance given betas and LD cholesky factor

    s2g := beta' V beta = beta' L L' beta

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param beta: numpy.ndarray genotype effects

    :return: float s2g
    """
    Ltb = np.dot(L.T, beta)
    s2g = np.dot(Ltb.T, Ltb)

    return s2g


def pick_rand_gene(plink_path, geno_prefix_path, random_genes_folder):
    '''
        :plink_path: path to plink binary
        :geno_prefix_path: prefix path of the genotype file (no bed, bim, fam)
        :random_genes_folder: path to folder, where random genes are stored
        :return: random gene output name
    '''
    while True:
        # choose random chromosome from 1 to 22
        chr = np.random.randint(1, 23)
        print(f"Chose chromosome {chr}")
        file = f"{geno_prefix_path}{chr}.bim"

        # Count the number of lines in the file
        result = subprocess.run(['wc', '-l', file], stdout=subprocess.PIPE, text=True)
        num_variants = int(result.stdout.split()[0])  # Extract the number of lines

        # choose random line and print the content of it
        line = np.random.randint(1, num_variants + 1)
        result = subprocess.run(['sed', '-n', f"{line}p", file],
                                stdout=subprocess.PIPE, text=True)
        random_variant = result.stdout.strip()

        # in the chosen line, extract 4th column, which is position
        random_position = int(random_variant.strip().split()[3])
        start_pos = max(0, random_position - WINDOW)
        end_pos = random_position + WINDOW

        print(f"Random Position: {random_position}")
        print(f"Start Position: {start_pos}")
        print(f"End Position: {end_pos}")
        geno_prefix = os.path.basename(geno_prefix_path)
        output = f'{random_genes_folder}/{geno_prefix}{chr}_{random_position}'
        
        # run plink to create a bed
        result = subprocess.run(
            f"ls {random_genes_folder} | grep {str(random_position)} | wc -l",
            stdout=subprocess.PIPE, text=True, shell=True
        )
        print(f"ls {random_genes_folder} | grep {str(random_position)} | wc -l")

        if int(result.stdout.strip()) != 0:
            print("Gene already exists. Rerunning...")
            # Remove existing files
            subprocess.run(f"rm {output}*", shell=True)
            continue

        result = subprocess.run([plink_path, 
                                '--bfile', f'{geno_prefix_path}{chr}', 
                                '--chr', str(chr),
                                '--from-bp', str(start_pos),
                                '--to-bp', str(end_pos),
                                '--make-bed',
                                '--out', output],
                                stdout=subprocess.PIPE, text=True)
        
        if not os.path.exists(f'{output}.bim'):
            print("No SNPs found in the region. Rerunning...")
            subprocess.run(f'rm {output}*', shell=True)
            continue

        result = subprocess.run(['wc', '-l', f'{output}.bim'],
                                stdout=subprocess.PIPE, text=True)
        vars_random_gene = int(result.stdout.split()[0])
        if vars_random_gene < 50:
            print('Not enough snps in chosen gene. Rerunning...')
            subprocess.run(f'rm {output}*', shell=True)
        else:
            print('Successfully picked a random gene')
            return output, str(chr)


def get_gcta_hsq(geno, gexpr, out, gcta):
    '''
        Calculates heritability estimates using GCTA
        :geno: numpy.ndarray n x p centered/scaled genotype matrix
        :gexpr: (numpy.ndarray, float) simulated phenotype, sd of Y
        :out: output folder
        :gcta: path to gcta binary file

        :returns: hsq, hsq_se, hsq_p - all floats
    '''
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
        command = f'rm {hsq_f}'
        subprocess.run(command.split())
        return hsq, hsq_se, hsq_p
    else:
        return 0, 10000, 1


def build_susie_sumstats(bim, geno, gexpr, sumstats_path, ld_path, gene,
                        gene_pref, plink_path, cNUM_PEOPLE):
    '''
        Builds SuSie Summary Statistics file

        :bim: bim file from plink random gene input
        :geno: numpy.ndarray n x p centered/scaled genotype matrix
        :gexpr: (numpy.ndarray, float) simulated phenotype, sd of Y
        :sumstats_path: path to where to save sumstats
        :ld_path: path to where to save ld files
        :gene: gene name from @pick_rand_gene
        :gene_pref: gene prefix with plink input

        :returns: filepath to fine mapped summary statistics

    '''
    betas = []
    pvals = []
    std_errs = []
    for i in range(geno.shape[1]):
        x = sm.add_constant(geno[:, i].reshape(-1,1))
        mod = sm.OLS(gexpr, x)
        results = mod.fit()
        betas.append(results.params[1])
        pvals.append(results.pvalues[1])
        std_errs.append(results.bse[1])
    sumstats = pd.DataFrame({
        'Gene': os.path.basename(gene),
        'SNP': bim['snp'],
        'A1': bim['a0'],
        'A2': bim['a1'],
        'BETA': betas,
        'SE': std_errs,
        'P': pvals
    })
    filename = os.path.join(sumstats_path, gene + ".txt")
    sumstats.to_csv(filename, sep="\t", index=False, header=True)
    ld_path_gene = os.path.join(ld_path, gene)

    subprocess.run([plink_path, '--bfile', gene_pref, '--r', 'inter-chr',
                    '--ld-window-r2', '0', '--keep-allele-order', '--extract',
                    filename, '--out', ld_path_gene])

    subprocess.run(['Rscript', 'run_susie.R',
                    '-c', filename,
                    '-l', ld_path_gene + ".ld",
                    '-n', str(cNUM_PEOPLE),
                    '-b', gene_pref + ".bim",
                    '-o', filename])
    return filename


def clean_up(out):
    command = f'rm {out}.grm.gz'
    subprocess.run(command.split())
    command = f'rm {out}.grm.id'
    subprocess.run(command.split())
    command = f'rm {out}.phen'
    subprocess.run(command.split())


def write_pheno(phen, out):
    out_file = open(out + '.phen', 'w')
    nindv = phen.shape[0]
    for i in range(nindv):
        out_file.write('{}\t{}\t{}\n'.format(i, i, phen[i]))
    out_file.close()


def write_grm(grm, nsnp, out):
    out_id_file = open(out+'.grm.id', 'w')
    out_file = gzip.open(out+'.grm.gz', 'w')
    nindv = grm.shape[0]
    for i in range(nindv):
        for j in range(0, i+1):
            out_file.write('{}\t{}\t{}\t{}\n'.format(i+1, j+1, nsnp,
                grm[i,j]).encode('utf-8'))
        out_id_file.write('{}\t{}\n'.format(i, i))
    out_file.close()
    out_id_file.close()


def run_gcta(gcta, out):
    command = f'{gcta} --grm-gz {out} --pheno {out + ".phen"} --reml --reml-no-constrain --out {out}'
    subprocess.run(command.split())
    # command = f'rm {out}.log'
    # subprocess.run(command.split())


def resolve_exe(exe_name_or_path: str) -> str:
    # If user passed a path, use it; otherwise search PATH
    if os.path.sep in exe_name_or_path:
        return exe_name_or_path
    resolved = shutil.which(exe_name_or_path)
    if resolved is None:
        raise RuntimeError(
            f"Cannot find executable '{exe_name_or_path}' on PATH. "
            f"Install it (e.g., conda install -c bioconda plink) or pass a full path."
        )
    return resolved


def arg_parser():
    parser = argparse.ArgumentParser(
        prog='Build Simulations for heritability analysis: SuSiE vs GCTA',
        description='''This brief command line creates a simulation for SuSiE
         vs GCTA heritability estimates. We simulate for eQTLs with different
         gene heritabilities as well as different number of causal SNPs.'''
    )
    parser.add_argument("--plink", default="plink",
                        type=str, help="Path to plink binary file")
    parser.add_argument("--gcta", default="gcta64",
                        help="gcta executable name or full path (default: gcta64 on PATH)")
    parser.add_argument('--geno_prefix', required=True,
                        type=str, help="Path to geno file")
    parser.add_argument('--tmpdir', required=True,
                        type=str, help="Path to temporary directory")
    parser.add_argument('--current_date', required=True, type=str,
                        help="Current date")

    # --- High-LD-region filtering (optional but enabled by default) ---
    parser.add_argument('--ld_mean_r2_threshold', type=float, default=0.05,
                        help=("Require regions to have mean pairwise r^2 (estimated quickly on a SNP subset) "
                              ">= this threshold before running simulations. Set to a negative value to disable."))
    parser.add_argument('--ld_score_max_snps', type=int, default=200,
                        help="Max SNPs to use for the quick LD score (smaller = faster).")
    parser.add_argument('--num_regions', type=int, default=25,
                        help="Number of accepted regions (genes) to simulate (default: 25).")
    return parser.parse_args()
def create_directory(directory):
    try:
        os.mkdir(directory)
        print(f"Directory '{directory}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory}' already exists.")
    except PermissionError:
        print(f"Permission denied: Unable to create '{directory}'.")
    except Exception as e:
        print(f"An error occurred: {e}")
        
        
def calculate_analytical_se(susie_output_file):
    """
    Calculates SE using the 4th Moment Formula:
    Var(h2) = Sum( E[b^4] - (E[b^2])^2 )
    Assumes SuSiE output has 'post_mean' and 'post_sd' columns.
    """
    try:
        df = pd.read_csv(susie_output_file, sep=" ")
        
        # NOTE: Verify these column names match your R script output!
        if 'POSTERIOR' not in df.columns or 'POSTERIOR_SD' not in df.columns:
            print("WARNING: 'post_mean' or 'post_sd' columns missing. Cannot calc Analytical SE.")
            return np.nan

        mu = df['POSTERIOR']
        sigma = df['POSTERIOR_SD']
        
        # 2nd Moment: E[X^2] = mu^2 + sigma^2
        E_b2 = mu**2 + sigma**2
        
        # 4th Moment (Normal Assumption): E[X^4] = mu^4 + 6*mu^2*sigma^2 + 3*sigma^4
        E_b4 = mu**4 + 6 * (mu**2) * (sigma**2) + 3 * (sigma**4)
        
        # Variance of the squared term
        var_b2 = E_b4 - (E_b2**2)
        
        # Total SE (Sum of Variances, assuming independence)
        se_h2 = np.sqrt(np.sum(var_b2))
        return se_h2
        
    except Exception as e:
        print(f"Error in Analytical SE calc: {e}")
        return np.nan

def get_region_info(prefix):
    """
    Reads PLINK data, computes MAFs, standardizes genotypes, and computes 
    the Cholesky decomposition of the LD matrix.
    """
    bim, fam, G = read_plink(prefix, verbose=False)
    G = G.T # Transpose to have SNPs as columns

    n_samples, n_snps = G.shape
    
    # 1. Calculate the TRUE allele frequency (AF) for math centering
    afs = (np.mean(G, axis=0) / 2).compute()
    
    # 2. Calculate the folded MAF (<= 0.5) for returning/filtering
    mafs = np.minimum(afs, 1.0 - afs)

    # Center and standardize G using the TRUE allele frequencies (afs)
    G -= afs * 2
    std_devs = np.std(G, axis=0)
    std_devs[std_devs == 0] = 1
    G /= std_devs

    # Regularize the LD matrix to make it Positive Semi Definite
    LD = np.dot(G.T, G) / n_samples
    regularization = 0.1
    LD += np.eye(n_snps) * regularization

    # re-adjust to get proper correlation matrix
    LD = LD / (1 + regularization)

    # compute cholesky decomp for faster sampling/simulation
    try:
        L = linalg.cholesky(LD, lower=True)
    except linalg.LinAlgError:
        # Fallback: If the matrix is still exceptionally stubborn, add a tiny 
        # bit more regularization to force positive-definiteness so the sim doesn't crash.
        LD += np.eye(n_snps) * 0.01
        LD = LD / 1.01
        L = linalg.cholesky(LD, lower=True)

    return (bim, L, mafs)





def main():
    args = arg_parser()
    plink_exe = resolve_exe(args.plink)
    gcta_exe  = resolve_exe(args.gcta)

    rgenes_folder = os.path.join(args.tmpdir, 'rgenes')
    create_directory(rgenes_folder)
    current_date = args.current_date
    sumstats_path = os.path.join(args.tmpdir, 'sumstats')
    create_directory(sumstats_path)
    print('Sumstats_path', sumstats_path)

    ld_path_directory = os.path.join(args.tmpdir, 'ld')
    create_directory(ld_path_directory)
    print('LD_path_directory', ld_path_directory)

    gcta_path_directory = os.path.join(args.tmpdir, 'hsq')
    create_directory(gcta_path_directory)
    print('gcta_path_directory', gcta_path_directory)

    # Output file tailored for the MAF/Effect Architecture tests
    output_txt = f"gcta_result_maf_arch_{current_date}.txt"
    if not os.path.exists(output_txt) or os.path.getsize(output_txt) == 0:
        with open(output_txt, "a") as f:
            f.write(OUTPUT_COLUMNS)

    # We can just loop normally now since we aren't rejecting regions based on LD
    for i in range(args.num_regions):
        random_gene, chr = pick_rand_gene(plink_exe, args.geno_prefix, rgenes_folder)

        # Get region info directly (no LD filtering)
        bim, L, mafs = get_region_info(random_gene)
        
        print(f"Simulating region {i+1}/{args.num_regions}: {random_gene}")

        for cNUM_PEOPLE in NUM_PEOPLE:
            Z = sim_geno(L, cNUM_PEOPLE)
            for ncausal in ncausal_arr:
                for h2 in h2_arr:
                    for maf_name, maf_range in maf_bins.items():
                        for effect_dist in effect_dists:
                            
                            max_retries = 5
                            success = False
                            
                            for attempt in range(max_retries):
                                try:
                                    # Unique ID for this specific simulation iteration
                                    curr_gene = f'{os.path.basename(random_gene)}_{ncausal}_{h2:.3f}_{maf_name}_{effect_dist}_{cNUM_PEOPLE}'
                                    
                                    b = sim_beta(L, ncausal, h2, mafs=mafs, maf_range=maf_range, effect_dist=effect_dist)
                                    gexpr = sim_trait(np.dot(Z, b), h2)[0]
                                    
                                    h2g, hsq_se, hsq_p = get_gcta_hsq(Z,
                                                                      gexpr,
                                                                      os.path.join(gcta_path_directory, curr_gene),
                                                                      gcta_exe)
                                    susie_result = build_susie_sumstats(bim,
                                                         Z, 
                                                         gexpr,
                                                         sumstats_path,
                                                         ld_path_directory,
                                                         curr_gene,
                                                         args.geno_prefix + chr, plink_exe, 
                                                         cNUM_PEOPLE)

                                    df = pd.read_csv(susie_result, sep=" ")
                                    
                                    if 'POSTERIOR' not in df.columns or 'POSTERIOR_SD' not in df.columns:
                                        raise KeyError("SuSiE output missing POSTERIOR columns.")

                                    est_h2 = np.sum(df['POSTERIOR']**2 + df['POSTERIOR_SD']**2)
                                    analytical_se_result = calculate_analytical_se(susie_result)
                                    
                                    with open(output_txt, "a") as f:
                                        f.write(f"{curr_gene}{SEP}{ncausal}{SEP}{h2}{SEP}{maf_name}{SEP}{effect_dist}{SEP}{h2g}{SEP}{hsq_se}{SEP}{hsq_p}{SEP}{est_h2}{SEP}{analytical_se_result}\n")
                                    
                                    success = True
                                    break 
                                    
                                except KeyError as e:
                                    print(f"Attempt {attempt + 1}/{max_retries} failed for {curr_gene}: {e}. Retrying...")
                                except pd.errors.EmptyDataError:
                                    print(f"Attempt {attempt + 1}/{max_retries} failed. SuSiE output empty for {curr_gene}. Retrying...")
                                except Exception as e:
                                    print(f"Attempt {attempt + 1}/{max_retries} failed with an unexpected error: {e}. Retrying...")
                            
                            if not success:
                                print(f"Skipping {curr_gene} entirely after {max_retries} failed attempts.")

if __name__ == "__main__":
    main()