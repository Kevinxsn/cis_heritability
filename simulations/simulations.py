# This script is used to simulate how GCTA and SuSiE estimate heritability in
# eQTLs with low number of causal SNPs and high number of causal SNPs
# A lot of code for this simulation was taken from:
#   https://github.com/mancusolab/twas_sim?tab=readme-ov-file

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
ncausal_arr = [1, 2, 5, 10] # number of causal SNPs for eQTL file
h2_arr = [0.05, 0.1, 0.2, 0.4, 0.6] # True heritability value
SEP = "\t"
OUTPUT_COLUMNS = f"gene{SEP}hsq{SEP}hsq_se{SEP}hsq_p\n"

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


def sim_beta(L, ncausal, eqtl_h2, rescale=True):
    """
    Sample qtl effects under a specified architecture.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param ncausal: int number of causal SNPs
    :param eqtl_h2: float the heritability of gene expression
    :param rescale: bool whether to rescale effects such that var(b V b) = h2 (default=True)

    :return: numpy.ndarray of causal effects
    """

    n_snps = L.shape[0]
    if eqtl_h2 > 0 and ncausal > 0:
        # choose how many eQTLs
        n_qtls = ncausal
        c_qtls = np.random.choice(range(int(n_snps)), n_qtls)
        b_qtls = np.zeros(int(n_snps))

        # sample effects from normal prior
        b_qtls[c_qtls] = np.random.normal(
            loc=0, scale=np.sqrt(eqtl_h2 / n_qtls), size=n_qtls
        )
        if rescale:
            s2g = compute_s2g(L, b_qtls)
            b_qtls *= np.sqrt(eqtl_h2 / s2g)
    
    else:
        b_qtls = np.zeros(int(n_snps))
    
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

    output_txt = f"gcta_result_eur_num_people_{current_date}.txt"
    if not os.path.exists(output_txt) or os.path.getsize(output_txt) == 0:
        # If the file doesn't exist or is empty, write the header
        with open(output_txt, "a") as f:
            f.write(OUTPUT_COLUMNS)

    # TODO: Change to 1000.
    for i in range(25):
        random_gene, chr = pick_rand_gene(plink_exe,
                                          args.geno_prefix, rgenes_folder)
        bim, L = get_ld(random_gene)
        for cNUM_PEOPLE in NUM_PEOPLE:
            Z = sim_geno(L, cNUM_PEOPLE)
            for ncausal in ncausal_arr:
                for h2 in h2_arr:
                    curr_gene = f'{os.path.basename(random_gene)}_{ncausal}_{h2:.3f}_{cNUM_PEOPLE}'
                    b = sim_beta(L, ncausal, h2)
                    gexpr = sim_trait(np.dot(Z, b), h2)[0]
                    h2g, hsq_se, hsq_p = get_gcta_hsq(Z,
                                                      gexpr,
                                                      os.path.join(gcta_path_directory, curr_gene),
                                                      gcta_exe)
                    build_susie_sumstats(bim,
                                         Z, 
                                         gexpr,
                                         sumstats_path,
                                         ld_path_directory,
                                         curr_gene,
                                         args.geno_prefix + chr,plink_exe, 
                                         cNUM_PEOPLE)

                    with open(output_txt, "a") as f:
                        f.write(f"{curr_gene}{SEP}{h2g}{SEP}{hsq_se}{SEP}{hsq_p}\n")


if __name__ == "__main__":
    main()
