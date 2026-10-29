import sys,os,time,json
sys.path.insert(0,"/projects/ps-renlab2/sux002/cis_heritability/simulations")
import numpy as np,statsmodels.api as sm,subprocess,pandas as pd
import simulations_test as sim
region="bench_regions/reg_22_kb100"; n=500
PLINK="/projects/ps-renlab2/sunan/anaconda3/envs/dsc180/bin/plink"
GCTA="/projects/ps-renlab2/sunan/anaconda3/envs/dsc180/bin/gcta64"
work="bench_work/fullgrid"; os.makedirs(work,exist_ok=True)
ncausal_arr=[1,3,5,10,20]; h2_arr=[0.05,0.1,0.2,0.4,0.6,0.8]
t0=time.time()
bim,L=sim.get_ld(region); Z=sim.sim_geno(L,n)
setup=time.time()-t0
ncell=0; tcell0=time.time()
for nc in ncausal_arr:
    for h2 in h2_arr:
        curr=f"c{ncell}"
        b=sim.sim_beta(L,nc,h2); g=sim.sim_trait(np.dot(Z,b),h2)[0]
        sim.get_gcta_hsq(Z,g,os.path.join(work,curr),GCTA)
        betas=[];ses=[];pv=[]
        for i in range(Z.shape[1]):
            x=sm.add_constant(Z[:,i].reshape(-1,1)); r=sm.OLS(g,x).fit()
            betas.append(r.params[1]);ses.append(r.bse[1]);pv.append(r.pvalues[1])
        ss=pd.DataFrame({"Gene":curr,"SNP":bim["snp"],"A1":bim["a0"],"A2":bim["a1"],"BETA":betas,"SE":ses,"P":pv})
        fn=os.path.join(work,curr+".txt"); ss.to_csv(fn,sep="\t",index=False)
        ldg=os.path.join(work,curr+"_ld")
        subprocess.run([PLINK,"--bfile",region,"--r","inter-chr","--ld-window-r2","0","--keep-allele-order","--extract",fn,"--out",ldg],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        subprocess.run(["Rscript","run_susie.R","-c",fn,"-l",ldg+".ld","-n",str(n),"-b",region+".bim","-o",fn],stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        for f in [fn,ldg+".ld",ldg+".log"]:
            try:os.remove(f)
            except OSError:pass
        ncell+=1
grid=time.time()-tcell0
tot=time.time()-t0
print(json.dumps({"p":L.shape[0],"n":n,"ncells":ncell,"setup_s":setup,"grid_s":grid,"total_s":tot,"per_cell_s":grid/ncell,"per_region_min":tot/60}))
