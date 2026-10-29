"""Find the variant filter that reproduces the old 1000G.EUR.22 panel.

The old panel is gone; what survives is a PLINK extract of chr22:17883968-19883968
(expand_proposal/bench_regions/reg_22_kb1000.bim, 17,287 SNPs) and the chr22
total (308,143 variants, from that extract's log).  Compare candidate filters
on a freshly converted, unfiltered EUR chr22 against both.

    python -I calibrate_1000g.py <raw22 prefix> <old window .bim>
"""
import sys

import pandas as pd

raw, old_bim = sys.argv[1], sys.argv[2]
LO, HI, OLD_TOTAL = 17883968, 19883968, 308143

bim = pd.read_csv(raw + ".bim", sep=r"\s+", header=None,
                  names=["chr", "id", "cm", "pos", "a1", "a2"])
frq = pd.read_csv(raw + ".frq", sep=r"\s+")
assert (bim["id"].values == frq["SNP"].values).all()
bim["maf"] = frq["MAF"].values
bim["mac"] = (frq["MAF"] * frq["NCHROBS"]).round().astype(int)
dup = bim["pos"].duplicated(keep=False)
rs = bim["id"].str.startswith("rs")  # v5b VCFs carry no rsIDs; kept for completeness

old = pd.read_csv(old_bim, sep=r"\s+", header=None,
                  names=["chr", "id", "cm", "pos", "a1", "a2"])
def keyset(df):
    return set(zip(df["pos"], [frozenset(x) for x in zip(df["a1"], df["a2"])]))


old_keys = keyset(old)
print(f"raw chr22 variants: {len(bim):,}   old window: {len(old):,}   old chr22 total: {OLD_TOTAL:,}")

rows = []
for rs_only in (False, True):
    for drop_dup in (False, True):
        for name, keep in [("mac>=1", bim["mac"] >= 1), ("mac>=2", bim["mac"] >= 2),
                           ("maf>=0.005", bim["maf"] >= 0.005), ("maf>=0.01", bim["maf"] >= 0.01),
                           ("maf>=0.05", bim["maf"] >= 0.05)]:
            m = keep.copy()
            if rs_only:
                m &= rs
            if drop_dup:
                m &= ~dup
            sel = bim[m]
            win = sel[(sel["pos"] >= LO) & (sel["pos"] <= HI)]
            keys = keyset(win)
            a1_match = (win.merge(old, on="pos", suffixes=("", "_o"))
                        .eval("a1 == a1_o").mean())
            rows.append(dict(rs_only=rs_only, drop_dup_pos=drop_dup, freq=name,
                             chr22_total=len(sel), window=len(win),
                             only_new=len(keys - old_keys), only_old=len(old_keys - keys),
                             a1_match=round(float(a1_match), 4)))
out = pd.DataFrame(rows).sort_values(["only_new", "only_old"], key=lambda s: s.abs())
out["total_diff"] = out["chr22_total"] - OLD_TOTAL
pd.set_option("display.width", 200)
print(out.to_string(index=False))

# show the disagreements for the polymorphic-in-EUR set (mac>=1)
sel = bim[bim["mac"] >= 1]
win = sel[(sel["pos"] >= LO) & (sel["pos"] <= HI)]
extra = keyset(win) - old_keys
missing = old_keys - keyset(win)
print("\nmac>=1 window, in new only:")
print(win[[ (p, frozenset((a, b))) in extra for p, a, b in zip(win["pos"], win["a1"], win["a2"])]].to_string())
print("in old only:", sorted(missing, key=lambda k: k[0])[:20])
print("\nmac>=1 chr22 positions seen twice:", int(sel["pos"].duplicated(keep=False).sum()))
print(sel[sel["pos"].duplicated(keep=False)].head(10).to_string())
