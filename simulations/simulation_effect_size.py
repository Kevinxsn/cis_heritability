# This script simulates how GCTA and SuSiE estimate heritability 
# focusing on different EFFECT SIZE architectures, independent of LD structure.

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

NUM_PEOPLE = [300]
WINDOW = 500000
ncausal_arr = [1, 3, 5, 10, 20] 
h2_arr = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8] 

# Define the effect size distributions we want to test
EFFECT_DISTS = ['normal', 'fixed', 'pareto']

SEP = "\t"
OUTPUT_COLUMNS = f"gene{SEP}ncausal{SEP}h2{SEP}effect_dist{SEP}hsq{SEP}hsq_se{SEP}hsq_p{SEP}su_hsq{SEP}su_hsq_se\n"

def get_ld(prefix):
    """
    Computes the LD matrix and its Cholesky decomposition without any filtering.
    """
    bim, fam, G = read_plink(prefix, verbose=False)
    G = G.T  

    n_samples, n_snps = G.shape

    mafs = (np.mean(G, axis=0) / 2).compute()
    G -= mafs * 2
    std_devs = np.std(G, axis=0)
    std_devs[std_devs == 0] = 1
    G /= std_devs

    # Compute full LD and apply regularization for Cholesky stability
    LD = np.dot(G.T, G) / n_samples
    regularization = 0.1
    LD += np.eye(n_snps) * regularization
    LD = LD / (1 + regularization)

    L = linalg.cholesky(LD, lower=True)
    return bim, L

def sim_geno(L, n): 
    p, p = L.shape
    Z = L.dot(np.random.normal(size=(n, p)).T).T 
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0) 
    return Z

def sim_beta(L, ncausal, eqtl_h2, effect_dist='normal', rescale=True):
    n_snps = L.shape[0]
    if eqtl_h2 > 0 and ncausal > 0:
        c_qtls = np.random.choice(range(int(n_snps)), ncausal, replace=False)
        b_qtls = np.zeros(int(n_snps))

        # Apply different effect size architectures
        if effect_dist == 'normal':
            effects = np.random.normal(loc=0, scale=np.sqrt(eqtl_h2 / ncausal), size=ncausal)
            
        elif effect_dist == 'fixed':
            signs = np.random.choice([-1, 1], size=ncausal)
            effects = signs * np.sqrt(eqtl_h2 / ncausal)
            
        elif effect_dist == 'pareto':
            raw_effects = np.random.pareto(a=2.0, size=ncausal)
            signs = np.random.choice([-1, 1], size=ncausal)
            effects = signs * raw_effects
        else:
            raise ValueError(f"Unknown effect_dist: {effect_dist}")

        b_qtls[c_qtls] = effects
        
        # Rescale ensures true h2 remains exactly what we set
        if rescale:
            s2g = compute_s2g(L, b_qtls)
            if s2g > 0:
                b_qtls *= np.sqrt(eqtl_h2 / s2g)
    else:
        b_qtls = np.zeros(int(n_snps))
        
    return b_qtls

def sim_trait(g, h2g): 
    n = len(g)
    if h2g > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ( (1.0 / h2g ) - 1 )
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e 
    else:
        e = np.random.normal(0, 1, n)
        y = e
    y -= np.mean(y)
    y_std = np.std(y)
    y /= y_std
    y_std = y_std.item()
    return y, y_std  

def compute_s2g(L, beta):
    Ltb = np.dot(L.T, beta)
    s2g = np.dot(Ltb.T, Ltb)
    return s2g

def pick_rand_gene(plink_path, geno_prefix_path, random_genes_folder):
    while True:
        chr = np.random.randint(1, 23)
        file = f"{geno_prefix_path}{chr}.bim"

        result = subprocess.run(['wc', '-l', file], stdout=subprocess.PIPE, text=True)
        num_variants = int(result.stdout.split()[0]) 

        line = np.random.randint(1, num_variants + 1)
        result = subprocess.run(['sed', '-n', f"{line}p", file], stdout=subprocess.PIPE, text=True)
        random_variant = result.stdout.strip()

        random_position = int(random_variant.strip().split()[3])
        start_pos = max(0, random_position - WINDOW)
        end_pos = random_position + WINDOW

        geno_prefix = os.path.basename(geno_prefix_path)
        output = f'{random_genes_folder}/{geno_prefix}{chr}_{random_position}'
        
        result = subprocess.run(f"ls {random_genes_folder} | grep {str(random_position)} | wc -l", stdout=subprocess.PIPE, text=True, shell=True)

        if int(result.stdout.strip()) != 0:
            subprocess.run(f"rm {output}*", shell=True)
            continue

        result = subprocess.run([plink_path, '--bfile', f'{geno_prefix_path}{chr}', '--chr', str(chr),
                                '--from-bp', str(start_pos), '--to-bp', str(end_pos),
                                '--make-bed', '--out', output], stdout=subprocess.PIPE, text=True)
        
        if not os.path.exists(f'{output}.bim'):
            subprocess.run(f'rm {output}*', shell=True)
            continue

        result = subprocess.run(['wc', '-l', f'{output}.bim'], stdout=subprocess.PIPE, text=True)
        vars_random_gene = int(result.stdout.split()[0])
        if vars_random_gene < 50:
            subprocess.run(f'rm {output}*', shell=True)
        else:
            return output, str(chr)

def get_gcta_hsq(geno, gexpr, out, gcta):
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
        subprocess.run(f'rm {hsq_f}'.split())
        return hsq, hsq_se, hsq_p
    else:
        return 0, 10000, 1

def build_susie_sumstats(bim, geno, gexpr, sumstats_path, ld_path, gene, gene_pref, plink_path, cNUM_PEOPLE):
    betas, pvals, std_errs = [], [], []
    for i in range(geno.shape[1]):
        x = sm.add_constant(geno[:, i].reshape(-1,1))
        mod = sm.OLS(gexpr, x)
        results = mod.fit()
        betas.append(results.params[1])
        pvals.append(results.pvalues[1])
        std_errs.append(results.bse[1])
        
    sumstats = pd.DataFrame({'Gene': os.path.basename(gene), 'SNP': bim['snp'], 'A1': bim['a0'],
                             'A2': bim['a1'], 'BETA': betas, 'SE': std_errs, 'P': pvals})
    filename = os.path.join(sumstats_path, gene + ".txt")
    sumstats.to_csv(filename, sep="\t", index=False, header=True)
    
    ld_path_gene = os.path.join(ld_path, gene)
    subprocess.run([plink_path, '--bfile', gene_pref, '--r', 'inter-chr', '--ld-window-r2', '0', 
                    '--keep-allele-order', '--extract', filename, '--out', ld_path_gene])
    subprocess.run(['Rscript', 'run_susie.R', '-c', filename, '-l', ld_path_gene + ".ld",
                    '-n', str(cNUM_PEOPLE), '-b', gene_pref + ".bim", '-o', filename])
    return filename

def clean_up(out):
    subprocess.run(f'rm {out}.grm.gz'.split())
    subprocess.run(f'rm {out}.grm.id'.split())
    subprocess.run(f'rm {out}.phen'.split())

def write_pheno(phen, out):
    with open(out + '.phen', 'w') as out_file:
        for i in range(phen.shape[0]):
            out_file.write(f'{i}\t{i}\t{phen[i]}\n')

def write_grm(grm, nsnp, out):
    out_id_file = open(out+'.grm.id', 'w')
    out_file = gzip.open(out+'.grm.gz', 'w')
    nindv = grm.shape[0]
    for i in range(nindv):
        for j in range(0, i+1):
            out_file.write(f'{i+1}\t{j+1}\t{nsnp}\t{grm[i,j]}\n'.encode('utf-8'))
        out_id_file.write(f'{i}\t{i}\n')
    out_file.close()
    out_id_file.close()

def run_gcta(gcta, out):
    subprocess.run(f'{gcta} --grm-gz {out} --pheno {out + ".phen"} --reml --reml-no-constrain --out {out}'.split())

def resolve_exe(exe_name_or_path: str) -> str:
    if os.path.sep in exe_name_or_path: return exe_name_or_path
    resolved = shutil.which(exe_name_or_path)
    if resolved is None: raise RuntimeError(f"Cannot find executable '{exe_name_or_path}' on PATH.")
    return resolved

def arg_parser():
    parser = argparse.ArgumentParser(description='Simulates SuSiE vs GCTA focusing on Effect Size distributions (LD Independent).')
    parser.add_argument("--plink", default="plink", type=str)
    parser.add_argument("--gcta", default="gcta64", type=str)
    parser.add_argument('--geno_prefix', required=True, type=str)
    parser.add_argument('--tmpdir', required=True, type=str)
    parser.add_argument('--current_date', required=True, type=str)
    parser.add_argument('--num_regions', type=int, default=25)
    return parser.parse_args()

def create_directory(directory):
    try: os.mkdir(directory)
    except FileExistsError: pass

def calculate_analytical_se(susie_output_file):
    try:
        df = pd.read_csv(susie_output_file, sep=" ")
        if 'POSTERIOR' not in df.columns or 'POSTERIOR_SD' not in df.columns: return np.nan
        mu, sigma = df['POSTERIOR'], df['POSTERIOR_SD']
        E_b2 = mu**2 + sigma**2
        E_b4 = mu**4 + 6 * (mu**2) * (sigma**2) + 3 * (sigma**4)
        var_b2 = E_b4 - (E_b2**2)
        return np.sqrt(np.sum(var_b2))
    except Exception:
        return np.nan

def main():
    args = arg_parser()
    plink_exe, gcta_exe = resolve_exe(args.plink), resolve_exe(args.gcta)

    rgenes_folder = os.path.join(args.tmpdir, 'rgenes')
    sumstats_path = os.path.join(args.tmpdir, 'sumstats')
    ld_path_directory = os.path.join(args.tmpdir, 'ld')
    gcta_path_directory = os.path.join(args.tmpdir, 'hsq')
    for d in [rgenes_folder, sumstats_path, ld_path_directory, gcta_path_directory]: create_directory(d)

    output_txt = f"effect_size_gcta_result_eur_{args.current_date}.txt"
    if not os.path.exists(output_txt) or os.path.getsize(output_txt) == 0:
        with open(output_txt, "a") as f: f.write(OUTPUT_COLUMNS)

    accepted = 0

    while accepted < args.num_regions:
        random_gene, chr = pick_rand_gene(plink_exe, args.geno_prefix, rgenes_folder)

        # Simply grab the bim and the Cholesky factor L, no filtering!
        bim, L = get_ld(random_gene)

        accepted += 1
        print(f"Accepted region {accepted}/{args.num_regions}")
        
        for cNUM_PEOPLE in NUM_PEOPLE:
            Z = sim_geno(L, cNUM_PEOPLE)
            for ncausal in ncausal_arr:
                for h2 in h2_arr:
                    # Loop through your different effect size distributions
                    for dist in EFFECT_DISTS:
                        
                        max_retries = 5
                        success = False
                        
                        for attempt in range(max_retries):
                            try:
                                curr_gene = f'{os.path.basename(random_gene)}_{ncausal}_{h2:.3f}_{dist}_{cNUM_PEOPLE}'
                                
                                # Pass the effect_dist down to sim_beta
                                b = sim_beta(L, ncausal, h2, effect_dist=dist)
                                gexpr = sim_trait(np.dot(Z, b), h2)[0]
                                
                                h2g, hsq_se, hsq_p = get_gcta_hsq(Z, gexpr, os.path.join(gcta_path_directory, curr_gene), gcta_exe)
                                
                                susie_result = build_susie_sumstats(bim, Z, gexpr, sumstats_path, ld_path_directory, 
                                                                    curr_gene, args.geno_prefix + chr, plink_exe, cNUM_PEOPLE)

                                df = pd.read_csv(susie_result, sep=" ")
                                
                                if 'POSTERIOR' not in df.columns or 'POSTERIOR_SD' not in df.columns:
                                    raise KeyError("SuSiE output missing POSTERIOR columns.")

                                est_h2 = np.sum(df['POSTERIOR']**2 + df['POSTERIOR_SD']**2)
                                analytical_se_result = calculate_analytical_se(susie_result)
                                
                                with open(output_txt, "a") as f:
                                    f.write(f"{curr_gene}{SEP}{ncausal}{SEP}{h2}{SEP}{dist}{SEP}{h2g}{SEP}{hsq_se}{SEP}{hsq_p}{SEP}{est_h2}{SEP}{analytical_se_result}\n")
                                
                                success = True
                                break 
                                
                            except KeyError as e:
                                print(f"Attempt {attempt + 1}/{max_retries} failed for {curr_gene}: {e}. Retrying...")
                            except pd.errors.EmptyDataError:
                                print(f"Attempt {attempt + 1}/{max_retries}: SuSiE output empty for {curr_gene}. Retrying...")
                        
                        if not success:
                            print(f"Skipping {curr_gene} entirely after {max_retries} failed attempts.")

if __name__ == "__main__":
    main()