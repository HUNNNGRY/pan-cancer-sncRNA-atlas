#!/usr/bin/env bash
# Partial correlation between sncRNA expression and numeric clinical variables.
# Inputs are produced by TCGA_sncRNA_cfPeak.R:
#   output/<DST>/<cancer>_sncRNA_log2CPM_topVar.txt
#   output/<DST>/pcc_clinical/<cancer>_clinical_num.txt
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

cd "${PROJECT_DIR}"
mkdir -p "output/${DST}/pcc_clinical"

METHOD="${METHOD:-pearson}"
for cell in "${CANCERS[@]}"; do
  out_file="output/${DST}/pcc_clinical/${cell}_clinical_num_${METHOD}.txt"
  if [ -s "${out_file}" ]; then
    echo "skip ${cell}"
    continue
  fi
  submit "${cell}_${DST}" "${PROJECT_DIR}/tmp/${DST}_${cell}.outerr" "
    $(activate_cmd "${ENV_PCC}")
    ${SHARED_UTILS_DIR}/run-covarPCC.R \
      --log2CPM output/${DST}/${cell}_sncRNA_log2CPM_topVar.txt \
      --log2TPM output/${DST}/pcc_clinical/${cell}_clinical_num.txt \
      --smpTbl ${DATA_DIR}/exSeek-dev/data/${DST}/sample_table2.txt \
      --covar Null -c ${CORES} -m ${METHOD} \
      -o ${out_file}
  "
done
