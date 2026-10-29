# This script is used to simulate how GCTA and SuSiE estimate heritability in
# eQTLs with low number of causal SNPs and high number of causal SNPs under LOW LD.

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
ncausal_arr = [1, 3, 5, 10, 20] # number of causal SNPs for eQTL file
h2_arr = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8] # True heritability value
SEP = "\t"
OUTPUT_COLUMNS = f"gene{SEP}parameter{SEP}hsq{SEP}hsq_se{SEP}hsq_p{SEP}su_hsq{SEP}su_hsq_se\n"

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


def get_ld(prefix, r2_pair_threshold=0.8, max_frac_pairs=0.001, score_max_snps=500, rng=None):
    """
    Reject high-LD regions BEFORE computing full LD/Cholesky.
    Accept ONLY if fraction of sampled SNP pairs with r^2 >= r2_pair_threshold is <= max_frac_pairs.
    """
    bim, fam, G = read_plink(prefix, verbose=False)
    G = G.T  # (n_samples x n_snps)

    n_samples, n_snps = G.shape

    mafs = (np.mean(G, axis=0) / 2).compute()
    G -= mafs * 2
    std_devs = np.std(G, axis=0)
    std_devs[std_devs == 0] = 1
    G /= std_devs

    frac = ld_frac_r2_ge_quick(
        G, r2_threshold=r2_pair_threshold, max_snps=score_max_snps, rng=rng
    )

    # MODIFIED: Reject if the fraction is HIGHER than the max_frac_pairs threshold (Low LD requirement)
    if max_frac_pairs is not None and frac > max_frac_pairs:
        return (None, None, frac)

    # proceed as before (full LD + regularization + Cholesky)
    LD = np.dot(G.T, G) / n_samples
    regularization = 0.1
    LD += np.eye(n_snps) * regularization
    LD = LD / (1 + regularization)

    L = linalg.cholesky(LD, lower=True)
    return (bim, L, frac)

def sim_geno(L, n): 
    p, p = L.shape
    Z = L.dot(np.random.normal(size=(n, p)).T).T #pxn
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0) 
    return Z

def sim_beta(L, ncausal, eqtl_h2, rescale=True):
    n_snps = L.shape[0]
    if eqtl_h2 > 0 and ncausal > 0:
        n_qtls = ncausal
        c_qtls = np.random.choice(range(int(n_snps)), n_qtls)
        b_qtls = np.zeros(int(n_snps))
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
        print(f"Chose chromosome {chr}")
        file = f"{geno_prefix_path}{chr}.bim"

        result = subprocess.run(['wc', '-l', file], stdout=subprocess.PIPE, text=True)
        num_variants = int(result.stdout.split()[0]) 

        line = np.random.randint(1, num_variants + 1)
        result = subprocess.run(['sed', '-n', f"{line}p", file],
                                stdout=subprocess.PIPE, text=True)
        random_variant = result.stdout.strip()

        random_position = int(random_variant.strip().split()[3])
        start_pos = max(0, random_position - WINDOW)
        end_pos = random_position + WINDOW

        print(f"Random Position: {random_position}")
        geno_prefix = os.path.basename(geno_prefix_path)
        output = f'{random_genes_folder}/{geno_prefix}{chr}_{random_position}'
        
        result = subprocess.run(
            f"ls {random_genes_folder} | grep {str(random_position)} | wc -l",
            stdout=subprocess.PIPE, text=True, shell=True
        )

        if int(result.stdout.strip()) != 0:
            print("Gene already exists. Rerunning...")
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
    subprocess.run(f'rm {out}.grm.gz'.split())
    subprocess.run(f'rm {out}.grm.id'.split())
    subprocess.run(f'rm {out}.phen'.split())

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

def resolve_exe(exe_name_or_path: str) -> str:
    if os.path.sep in exe_name_or_path:
        return exe_name_or_path
    resolved = shutil.which(exe_name_or_path)
    if resolved is None:
        raise RuntimeError(
            f"Cannot find executable '{exe_name_or_path}' on PATH. "
            f"Install it or pass a full path."
        )
    return resolved

def arg_parser():
    parser = argparse.ArgumentParser(
        prog='Build Simulations for heritability analysis: SuSiE vs GCTA (Low LD)',
        description='''Simulates SuSiE vs GCTA heritability estimates under LOW LD scenarios.'''
    )
    parser.add_argument("--plink", default="plink",
                        type=str, help="Path to plink binary file")
    parser.add_argument("--gcta", default="gcta64",
                        help="gcta executable name or full path")
    parser.add_argument('--geno_prefix', required=True,
                        type=str, help="Path to geno file")
    parser.add_argument('--tmpdir', required=True,
                        type=str, help="Path to temporary directory")
    parser.add_argument('--current_date', required=True, type=str,
                        help="Current date")

    # MODIFIED: Changed argument definition to reflect Low LD expectations
    parser.add_argument('--ld_max_frac_pairs', type=float, default=0.002,
                        help=("Require regions to have a MAXIMUM fraction of highly correlated SNPs (r^2 >= 0.8) "
                              "<= this threshold before running simulations."))
    parser.add_argument('--ld_score_max_snps', type=int, default=200,
                        help="Max SNPs to use for the quick LD score (smaller = faster).")
    parser.add_argument('--num_regions', type=int, default=25,
                        help="Number of accepted regions (genes) to simulate (default: 25).")
    return parser.parse_args()

def create_directory(directory):
    try:
        os.mkdir(directory)
    except FileExistsError:
        pass

def calculate_analytical_se(susie_output_file):
    try:
        df = pd.read_csv(susie_output_file, sep=" ")
        if 'POSTERIOR' not in df.columns or 'POSTERIOR_SD' not in df.columns:
            return np.nan

        mu = df['POSTERIOR']
        sigma = df['POSTERIOR_SD']
        
        E_b2 = mu**2 + sigma**2
        E_b4 = mu**4 + 6 * (mu**2) * (sigma**2) + 3 * (sigma**4)
        var_b2 = E_b4 - (E_b2**2)
        
        se_h2 = np.sqrt(np.sum(var_b2))
        return se_h2
    except Exception as e:
        return np.nan

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

    # Note: Using the low_ld output name from our previous setup
    output_txt = f"low_ld_gcta_result_eur_num_people_{current_date}.txt"
    if not os.path.exists(output_txt) or os.path.getsize(output_txt) == 0:
        # If the file doesn't exist or is empty, write the header
        with open(output_txt, "a") as f:
            f.write(OUTPUT_COLUMNS)

    rng = np.random.default_rng()

    # Apply the threshold for Low LD
    ld_thresh = args.ld_max_frac_pairs

    accepted = 0
    while accepted < args.num_regions:
        random_gene, chr = pick_rand_gene(plink_exe,
                                          args.geno_prefix, rgenes_folder)

        # Get LD, filtering for LOW LD
        bim, L, frac = get_ld(
            random_gene,
            r2_pair_threshold=0.8,
            max_frac_pairs=ld_thresh,
            score_max_snps=args.ld_score_max_snps,
            rng=rng
        )
        print(f"Quick LD score: frac(pairs with r^2>=0.8) = {frac:.4f}")

        if L is None:
            print("Region had TOO MUCH LD (exceeded max_frac_pairs). Removing and resampling...")
            subprocess.run(f'rm {random_gene}*', shell=True)
            continue

        accepted += 1
        for cNUM_PEOPLE in NUM_PEOPLE:
            Z = sim_geno(L, cNUM_PEOPLE)
            for ncausal in ncausal_arr:
                for h2 in h2_arr:
                    max_retries = 5
                    success = False
                    
                    for attempt in range(max_retries):
                        try:
                            curr_gene = f'{os.path.basename(random_gene)}_{ncausal}_{h2:.3f}_{cNUM_PEOPLE}'
                            
                            # Simulating inside the retry loop ensures new effects/noise are drawn if it fails
                            b = sim_beta(L, ncausal, h2)
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
                            
                            # Explicitly check for columns before calculation to trigger the except block
                            if 'POSTERIOR' not in df.columns or 'POSTERIOR_SD' not in df.columns:
                                raise KeyError("SuSiE output missing POSTERIOR columns.")

                            est_h2 = np.sum(df['POSTERIOR']**2 + df['POSTERIOR_SD']**2)
                            analytical_se_result = calculate_analytical_se(susie_result)
                            
                            with open(output_txt, "a") as f:
                                f.write(f"{curr_gene}{SEP}{ncausal}_{h2}{SEP}{h2g}{SEP}{hsq_se}{SEP}{hsq_p}{SEP}{est_h2}{SEP}{analytical_se_result}\n")
                            
                            success = True
                            

                            
                            break # Success! Break out of the retry loop
                            
                        except KeyError as e:
                            print(f"Attempt {attempt + 1}/{max_retries} failed for {curr_gene}: {e}. Retrying...")
                    
                    if not success:
                        print(f"Skipping {curr_gene} entirely after {max_retries} failed attempts.")
        #target_dir = args.tmpdir
        #command = f"cd {target_dir} && rm -rf -- *"

        #subprocess.run(command, shell=True, executable='/bin/bash')
        #print("Directory cleared.")
                    

if __name__ == "__main__":
    main()