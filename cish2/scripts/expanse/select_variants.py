"""Write the variant IDs to keep from an unfiltered EUR fileset.

The rule mirrors the candidate filters in calibrate_1000g.py; pick the row that
reproduced the old panel there.

    python -I select_variants.py <raw prefix> <out list> --freq maf>=0.01 [--rs-only] [--drop-dup-pos] [--exclude-ids FILE]
"""
import argparse

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("raw")
ap.add_argument("out")
ap.add_argument("--freq", required=True, help="mac>=K or maf>=X")
ap.add_argument("--rs-only", action="store_true")
ap.add_argument("--drop-dup-pos", action="store_true")
ap.add_argument("--exclude-ids", help="file of IDs to drop (sites the old panel lacked)")
a = ap.parse_args()

bim = pd.read_csv(a.raw + ".bim", sep=r"\s+", header=None,
                  names=["chr", "id", "cm", "pos", "a1", "a2"])
frq = pd.read_csv(a.raw + ".frq", sep=r"\s+")
assert len(bim) == len(frq)
# v5b VCFs have no IDs ('.'), so name variants chr:pos:A1:A2 here and rewrite the
# raw .bim.  A few sites appear twice in the VCF with identical alleles; both are
# kept (as PLINK keeps duplicate rsIDs) and the repeat gets a _2, _3 suffix.
base = bim["chr"].astype(str) + ":" + bim["pos"].astype(str) + ":" + bim["a1"] + ":" + bim["a2"]
rep = base.groupby(base).cumcount()
bim["id"] = base.where(rep == 0, base + "_" + (rep + 1).astype(str))
if (rep > 0).any():
    print(f"duplicate sites (suffixed): {list(bim.loc[rep > 0, 'id'])}")
bim.to_csv(a.raw + ".bim", sep="\t", header=False, index=False)
mac = (frq["MAF"] * frq["NCHROBS"]).round()

kind, thr = a.freq.split(">=")
keep = (mac >= int(thr)) if kind == "mac" else (frq["MAF"] >= float(thr))
if a.rs_only:
    keep &= bim["id"].str.startswith("rs")
if a.drop_dup_pos:
    keep &= ~bim["pos"].duplicated(keep=False)
if a.exclude_ids:
    keep &= ~bim["id"].isin(set(open(a.exclude_ids).read().split()))
# --extract matches by ID, so an ID kept here must be unique in the fileset
ids = bim.loc[keep, "id"]
assert not ids.duplicated().any(), "kept IDs are not unique; add --rs-only/--drop-dup-pos"
ids.to_csv(a.out, index=False, header=False)
print(f"{a.raw}: kept {keep.sum():,} of {len(bim):,}")
