#!/usr/bin/env bash
# Rewrite the absolute paths baked into code and configs after moving servers.
#
#   scripts/relocate_paths.sh <new_project_root> <new_data_root> [--dry-run]
#
#   new_project_root  where cis_heritability/ now lives (contains PROJECT.md)
#   new_data_root     the directory holding 1000G.EUR.<chr>.{bed,bim,fam}
#                     and gtex/ (was /projects/ps-renlab2/sux002/data)
#
# Only code and configs are touched (cish2/**/*.py|*.sh|*.yaml|*.R).  The
# config.yaml snapshots and logs under results/ are records of past runs and
# are deliberately left alone.  The conda env path in scripts/env.sh is not
# rewritten -- set CISH2_ENV instead.
set -euo pipefail

NEW_ROOT="${1:?usage: relocate_paths.sh <new_project_root> <new_data_root> [--dry-run]}"
NEW_DATA="${2:?usage: relocate_paths.sh <new_project_root> <new_data_root> [--dry-run]}"
DRY="${3:-}"

OLD_ROOT="/projects/ps-renlab2/sux002/cis_heritability"
OLD_DATA="/projects/ps-renlab2/sux002/data"

NEW_ROOT="${NEW_ROOT%/}"
NEW_DATA="${NEW_DATA%/}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # .../cish2

mapfile -t FILES < <(grep -rlF -e "${OLD_ROOT}" -e "${OLD_DATA}" "${HERE}" \
    --include='*.py' --include='*.sh' --include='*.yaml' --include='*.R' \
    | grep -v '/__pycache__/' | grep -v 'relocate_paths.sh$' || true)

if [ "${#FILES[@]}" -eq 0 ]; then
    echo "nothing to rewrite -- no old paths found under ${HERE}"
    exit 0
fi

echo "rewriting ${#FILES[@]} files:"
printf '  %s\n' "${FILES[@]}"
[ "${DRY}" = "--dry-run" ] && exit 0

for f in "${FILES[@]}"; do
    # project root first: it does not overlap the data root, but order keeps it obvious
    sed -i -e "s#${OLD_ROOT}#${NEW_ROOT}#g" -e "s#${OLD_DATA}#${NEW_DATA}#g" "$f"
done

left="$(grep -rnF -e "${OLD_ROOT}" -e "${OLD_DATA}" "${HERE}" \
    --include='*.py' --include='*.sh' --include='*.yaml' --include='*.R' \
    | grep -v 'relocate_paths.sh' || true)"
if [ -n "${left}" ]; then
    echo "WARNING: old paths still present:"; echo "${left}"; exit 1
fi
echo "done."
