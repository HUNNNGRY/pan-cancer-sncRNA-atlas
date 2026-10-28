#!/usr/bin/env bash
# Pairwise Jaccard similarity between per-condition peak BED files.
# Writes ${DATA_DIR}/exSeek-dev/test/${DST}/jaccard/{beds.txt,cmb.txt,jac.txt};
# TCGA_sncRNA_cfPeak.R reads jac.txt back to build the tissue similarity matrix.
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

test_dir="${DATA_DIR}/exSeek-dev/test"
jac_dir="${test_dir}/${DST}/jaccard"
mkdir -p "${jac_dir}"

ls "${test_dir}"/TCGA-*_small/*_*.bed_optimal_cfpeak.bed > "${jac_dir}/beds.txt"
echo "peak BED files: $(wc -l < "${jac_dir}/beds.txt")"

awk -f "${EXSEEK_SCRIPTS_DIR}/combinations.awk" \
  <<< "$(tr '\n' ' ' < "${jac_dir}/beds.txt")" > "${jac_dir}/cmb.txt"
echo "combinations: $(wc -l < "${jac_dir}/cmb.txt")"

bash "${EXSEEK_SCRIPTS_DIR}/getBedJaccard2.sh" "${jac_dir}/cmb.txt" "${jac_dir}"

echo "done: ${jac_dir}/jac.txt"
