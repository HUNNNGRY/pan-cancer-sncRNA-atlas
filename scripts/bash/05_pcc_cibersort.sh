#!/usr/bin/env bash
# Partial correlation between sncRNA expression and CIBERSORTx immune cell
# fractions (22 cell types, absolute mode), adjusted for tumour purity.
# CIBERSORTx_abs.txt is written by the CIBERSORTx section of TCGA_sncRNA_cfPeak.R.
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

cd "${PROJECT_DIR}"
mkdir -p "output/${DST}/pcc_CIBERSORT"

METHOD="${METHOD:-pearson}"
OUT_SURFIX="${OUT_SURFIX:-rmPurity}"
for cell in "${CANCERS[@]}"; do
  out_file="output/${DST}/pcc_CIBERSORT/${cell}_${METHOD}_${OUT_SURFIX}.txt"
  if [ -s "${out_file}" ]; then
    echo "skip ${cell}"
    continue
  fi
  submit "${cell}_${DST}_cibersort" "${PROJECT_DIR}/tmp/${DST}_${cell}_cibersort.outerr" "
    $(activate_cmd "${ENV_PCC}")
    ${SHARED_UTILS_DIR}/run-covarPCC.R \
      --log2CPM output/${DST}/${cell}_sncRNA_log2CPM_topVar.txt \
      --log2TPM ${TCGA_DIR}/data/${cell}_long/CIBERSORTx_abs.txt \
      --smpTbl ${DATA_DIR}/exSeek-dev/data/${DST}/sample_table2.txt \
      --covar purity -c ${CORES} -m ${METHOD} \
      -o ${out_file}
  "
done
