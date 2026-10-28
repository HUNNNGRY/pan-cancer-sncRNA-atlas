#!/usr/bin/env bash
# Partial correlation between sncRNA expression and mRNA expression (log2TPM),
# adjusted for tumour purity, one job per cancer type.
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

cd "${PROJECT_DIR}"
mkdir -p "output/${DST}/pcc_mRNA"

METHOD="${METHOD:-pearson}"
OUT_SURFIX="${OUT_SURFIX:-_rmPurity}"
for cell in "${CANCERS[@]}"; do
  out_file="output/${DST}/pcc_mRNA/${cell}_${METHOD}${OUT_SURFIX}.txt"
  if [ -s "${out_file}" ]; then
    echo "skip ${cell}"
    continue
  fi
  submit "${cell}_${DST}" "${PROJECT_DIR}/tmp/${DST}_${cell}.outerr" "
    $(activate_cmd "${ENV_PCC}")
    ${SHARED_UTILS_DIR}/run-covarPCC.R \
      --log2CPM output/${DST}/${cell}_sncRNA_log2CPM_topVar.txt \
      --log2TPM ${TCGA_DIR}/data/${cell}_long/log2TPM_filter.txt \
      --smpTbl ${DATA_DIR}/exSeek-dev/data/${DST}/sample_table2.txt \
      --covar purity -c ${CORES} -m ${METHOD} \
      -o ${out_file}
  "
done
# 4k x 20k matrix: 16 cores, ~3 h per cancer type
