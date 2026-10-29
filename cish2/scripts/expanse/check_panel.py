"""Check the rebuilt 1000G.EUR panel against regions recorded by the old runs.

Every results/<campaign>/regions.tsv stores, per simulated region, the window
and the .bim index range [i0, i1) it covered in the old panel.  Region draws
pick a variant *index*, so the old regions are reproduced only if the new panel
has the same positions in the same order; identical i0/i1 for every region is
the test.

    python -I check_panel.py <geno_prefix> <results_dir>
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

prefix, results = sys.argv[1], sys.argv[2]
regions = pd.concat(
    [pd.read_csv(f, sep="\t").assign(campaign=os.path.basename(os.path.dirname(f)))
     for f in sorted(glob.glob(os.path.join(results, "*", "regions.tsv")))],
    ignore_index=True)

rows = []
for chrom, grp in regions.groupby("chrom"):
    pos = pd.read_csv(f"{prefix}{chrom}.bim", sep=r"\s+", header=None, usecols=[3])[3].to_numpy()
    i0, i1 = np.searchsorted(pos, grp["start_bp"]), np.searchsorted(pos, grp["end_bp"] + 1)
    centre_ok = np.isin(grp["center_bp"], pos)
    for r, a, b, c in zip(grp.itertuples(), i0, i1, centre_ok):
        rows.append(dict(campaign=r.campaign, region_id=r.region_id, chrom=chrom,
                         old_i0=r.i0, new_i0=int(a), old_n=r.i1 - r.i0, new_n=int(b - a),
                         centre_in_panel=bool(c)))
out = pd.DataFrame(rows)
out["same_index"] = (out["old_i0"] == out["new_i0"]) & (out["old_n"] == out["new_n"])
out["same_count"] = out["old_n"] == out["new_n"]

print(f"{len(out)} recorded regions on chromosomes {sorted(out['chrom'].unique())}")
print(f"  identical index range : {out['same_index'].sum()}/{len(out)}")
print(f"  identical SNP count   : {out['same_count'].sum()}/{len(out)}")
print(f"  centre SNP present    : {out['centre_in_panel'].sum()}/{len(out)}")
by_chr = out.groupby("chrom").agg(regions=("same_index", "size"), same_index=("same_index", "sum"),
                                  max_i0_shift=("new_i0", lambda s: 0),
                                  median_count_diff=("new_n", "median"))
by_chr["max_i0_shift"] = (out["new_i0"] - out["old_i0"]).abs().groupby(out["chrom"]).max()
by_chr["median_count_diff"] = (out["new_n"] - out["old_n"]).groupby(out["chrom"]).median()
print(by_chr.to_string())
bad = out[~out["same_index"]]
if len(bad):
    print("\nfirst mismatches:")
    print(bad.head(15).to_string(index=False))
