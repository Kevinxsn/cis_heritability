#!/usr/bin/env bash
# Round 4 poster rebuild (docs/pi_feedback_round4_design.md §4), in two passes:
#   1. round 3 with the UNCHANGED commands, then compare every output against
#      the checksums recorded before any round-4 edit (logs/round4/round3_output_md5.txt)
#   2. round 4 on top: new tables, new figures, *_round4 deck, round-4 pins
#   sbatch scripts/expanse/poster_round4.sbatch
#SBATCH --job-name=poster_round4
#SBATCH --account=sdp202
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/home/sux002/work/cis_heritability/logs/round4/%x_%j.out
set -euo pipefail
source /home/sux002/work/cis_heritability/cish2/scripts/expanse/config.sh
source "$PROJ/cish2/scripts/env.sh"
cd "$PROJ/poster"

echo "===== pass 1: round 3, unchanged commands"
python revision_analyses.py
python make_figures.py
python build_poster.py
python check_numbers.py
echo "===== round-3 outputs vs pre-round-4 checksums"
(cd "$PROJ" && md5sum -c logs/round4/round3_output_md5.txt) || echo "!! some round-3 outputs differ (see above)"

echo "===== pass 2: round 4"
python revision_analyses.py e f g
python make_figures.py --round 4
python build_poster.py --round 4
python check_numbers.py --round 4
