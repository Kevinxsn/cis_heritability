#!/usr/bin/env bash
# Download 1000 Genomes Phase 3 (20130502) autosome VCFs + sample panel from EBI.
# Resumable: curl -C - continues partial files; finished files are re-checked with gzip -t.
#   sbatch scripts/expanse/download_1000g.sbatch
#SBATCH --job-name=cish2_dl1kg
#SBATCH --account=sdp202
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/home/sux002/work/cis_heritability/logs/%x_%j.out
set -euo pipefail
source /home/sux002/work/cis_heritability/cish2/scripts/expanse/config.sh
cd "$VCF_DIR"

BASE=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
curl -sSf --retry 5 -o integrated_call_samples_v3.20130502.ALL.panel \
    "$BASE/integrated_call_samples_v3.20130502.ALL.panel"
awk '$3=="EUR"{print $1, $1}' integrated_call_samples_v3.20130502.ALL.panel > EUR.keep
echo "EUR samples: $(wc -l < EUR.keep)"

# Single EBI connections sometimes stall at ~10 KB/s while a fresh one runs at
# ~8 MB/s, so drop any transfer slower than 200 KB/s for 60 s and resume it.
fetch() {
    local f="ALL.chr$1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    local want have
    [ -f "$f.ok" ] && [ -f "$f" ] && { echo "ok chr$1 (verified in an earlier run)"; return 0; }
    want=$(curl -sSfI --retry 5 "$BASE/$f" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
    for attempt in $(seq 1 30); do
        have=$(stat -c %s "$f" 2>/dev/null || echo 0)
        if [ "$have" -lt "$want" ]; then
            curl -sSf --retry 3 --speed-limit 200000 --speed-time 60 -C - -o "$f" "$BASE/$f" || true
            have=$(stat -c %s "$f" 2>/dev/null || echo 0)
        fi
        if [ "$have" -eq "$want" ]; then
            if gzip -t "$f"; then touch "$f.ok"; echo "ok chr$1 $(du -h "$f" | cut -f1)"; return 0; fi
            echo "corrupt chr$1 -- restarting from zero"; rm -f "$f"
        elif [ "$have" -gt "$want" ]; then
            echo "oversized chr$1 -- restarting from zero"; rm -f "$f"
        fi
    done
    echo "FAILED chr$1 ($have of $want bytes)"; return 1
}
export -f fetch; export BASE
seq 1 22 | xargs -P 4 -I{} bash -c 'fetch {}'
