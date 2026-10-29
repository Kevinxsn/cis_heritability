import pandas as pd
import numpy as np
import scipy.linalg as linalg
import subprocess
import os
import argparse
from pandas_plink import read_plink
import statsmodels.api as sm
import shutil
import csv
import os
from datetime import datetime


# --- CONFIGURATION ---
NUM_PEOPLE = 300
N_CAUSAL = 1
H2 = 0.2
N_REPEATS = 1  # Number of times to repeat SuSiE to find Empirical SE

def get_ld(prefix):
    """
    Reads PLINK files and generates a regularized LD matrix.
    (Simplified from original script)
    """
    print(f"Reading genotypes from {prefix}...")
    bim, fam, G = read_plink(prefix, verbose=False)
    G = G.T
    n_samples, n_snps = G.shape
    
    # Standardize Genotypes
    mafs = (np.mean(G, axis=0) / 2).compute()
    G -= mafs * 2
    std_devs = np.std(G, axis=0)
    std_devs[std_devs == 0] = 1
    G /= std_devs
    
    # LD Calculation with Regularization
    LD = np.dot(G.T, G) / n_samples
    regularization = 0.1
    LD += np.eye(n_snps) * regularization
    LD = LD / (1 + regularization)
    
    # Cholesky Decomposition
    L = linalg.cholesky(LD, lower=True)
    return bim, L

def sim_geno(L, n):
    """Generates Genotype Matrix Z from LD matrix"""
    p, p = L.shape
    Z = L.dot(np.random.normal(size=(n, p)).T).T
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0)
    return Z

def sim_trait(L, Z, ncausal, h2):
    """Simulates Phenotype Y based on Z and Heritability"""
    n_snps = L.shape[0]
    n_samples = Z.shape[0]
    
    # Simulate Betas
    b_qtls = np.zeros(n_snps)
    c_qtls = np.random.choice(range(n_snps), ncausal)
    b_qtls[c_qtls] = np.random.normal(loc=0, scale=np.sqrt(h2 / ncausal), size=ncausal)
    
    # Rescale betas to match target heritability exactly
    Ltb = np.dot(L.T, b_qtls)
    s2g = np.dot(Ltb.T, Ltb)
    b_qtls *= np.sqrt(h2 / s2g)
    
    # Simulate Phenotype
    g = np.dot(Z, b_qtls)
    s2g = np.var(g, ddof=1)
    s2e = s2g * ((1.0 / h2) - 1)
    y = g + np.random.normal(0, np.sqrt(s2e), n_samples)
    y = (y - np.mean(y)) / np.std(y)
    
    return y

def get_sumstats(geno, gexpr, bim):
    """Calculates GWAS summary statistics (Beta, SE, P)"""
    betas, std_errs, pvals = [], [], []
    for i in range(geno.shape[1]):
        x = sm.add_constant(geno[:, i])
        mod = sm.OLS(gexpr, x)
        res = mod.fit()
        betas.append(res.params[1])
        std_errs.append(res.bse[1])
        pvals.append(res.pvalues[1])
    
    return pd.DataFrame({
        'Gene': 'random',
        'SNP': bim['snp'],
        'A1': bim['a0'],
        'A2': bim['a1'],
        'BETA': betas,
        'SE': std_errs,
        'P': pvals
    })

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

def append_result(x, y, out_csv):
    file_exists = os.path.exists(out_csv)

    with open(out_csv, "a", newline="") as f:
        w = csv.writer(f)
        if not file_exists:
            w.writerow(["timestamp", "x", "y"])   # header once
        w.writerow([datetime.now().isoformat(timespec="seconds"), x, y])


def main():
    parser = argparse.ArgumentParser(description="Compare Empirical vs Analytical SE for SuSiE")
    parser.add_argument('--geno_prefix', required=True, help="Path to PLINK genotype prefix (no .bed extension)")
    parser.add_argument('--plink', default="plink", help="Path to PLINK executable")
    args = parser.parse_args()
    plink_exe = resolve_exe(args.plink)

    # --- STEP 1: Generate Data ONCE ---
    print("\n--- Generating Simulation Data ---")
    # Note: Using a single chunk of SNPs for the test. 
    # In practice, ensure this prefix points to a manageable region (e.g. 1000 SNPs)
    bim, L = get_ld(args.geno_prefix)
    Z = sim_geno(L, NUM_PEOPLE)
    y = sim_trait(L, Z, N_CAUSAL, H2)
    
    # Save Summary Stats
    sumstats = get_sumstats(Z, y, bim)
    sumstats_file = "temp_sumstats.txt"
    sumstats.to_csv(sumstats_file, sep="\t", index=False)
    print(f"Summary statistics saved to {sumstats_file}")

    # Create LD file for SuSiE (using PLINK)
    ld_file = "temp_ld"
    subprocess.run([plink_exe, '--bfile', args.geno_prefix, '--r', 'inter-chr', 
                    '--ld-window-r2', '0', '--extract', sumstats_file, 
                    '--out', ld_file], stdout=subprocess.DEVNULL)
    print(f"LD matrix saved to {ld_file}.ld")

    # --- STEP 2: Run SuSiE Repeatedly ---
    print(f"\n--- Running SuSiE {N_REPEATS} times ---")
    susie_h2_estimates = []
    analytical_se_result = None
    
    for i in range(N_REPEATS):
        output_file = f"susie_output_run_{i}.txt"
        
        # Calling the R script. Ensure 'run_susie.R' is in this folder.
        # Ensure your R script does NOT set a fixed seed, or this will be identical every time.
        cmd = ['Rscript', 'run_susie.R',
               '-c', sumstats_file,
               '-l', ld_file + ".ld",
               '-n', str(NUM_PEOPLE),
               '-b', args.geno_prefix + ".bim",
               '-o', output_file]
        
        subprocess.run(cmd, check=True)
        
        # Parse Output
        try:
            df = pd.read_csv(output_file, sep=" ")
            # Calculate h2 estimate from this run (Sum of squared effects)
            # Adjust 'post_mean' if your R script uses a different name
            est_h2 = np.sum(df['POSTERIOR']**2)
            susie_h2_estimates.append(est_h2)
            
            # Calculate Analytical SE (Only need to do this once, e.g., on first run)
            if i == 0:
                analytical_se_result = calculate_analytical_se(output_file)
                
        except Exception as e:
            print(f"Failed to parse run {i}: {e}")

    # --- STEP 3: Report Results ---
    empirical_se = np.std(susie_h2_estimates)
    mean_h2 = np.mean(susie_h2_estimates)

    
    print("\n" + "="*30)
    print("       FINAL RESULTS       ")
    print("="*30)
    print(f"True Heritability (Simulated): {H2}")
    print(f"SuSiE Mean Estimate (20 runs): {mean_h2:.5f}")
    print("-" * 30)
    print(f"1. EMPIRICAL SE (std of 20 runs): {empirical_se:.5f}")
    print(f"2. ANALYTICAL SE (formula):       {analytical_se_result:.5f}")
    print("-" * 30)
    
    if empirical_se < 1e-6:
        print(">> OBSERVATION: SuSiE appears DETERMINISTIC (SE ~ 0).")
        print(">> Bootstrapping will not work unless randomness is introduced.")
    elif analytical_se_result < empirical_se:
        print(">> OBSERVATION: Analytical SE is SMALLER than Empirical SE.")
        print(">> This implies covariance (LD) is increasing error, which the formula ignores.")
    else:
        print(">> OBSERVATION: Analytical SE is LARGER (conservative).")
    
    
    append_result(mean_h2, analytical_se_result, out_csv='resultH=0.2.csv')
if __name__ == "__main__":
    main()