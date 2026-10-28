#!/usr/bin/env bash
# Pathway GSEA on the sncRNA-mRNA partial correlation tables (PCC step output).
# For each cancer type this writes output/<DST>/pcc_mRNA/rmPurity_<DST>_<cancer>.gmtLong,
# which the ORA/GSEA scripts in this directory consume.
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

cd "${PROJECT_DIR}"

METHOD="${METHOD:-pearson}"
OUT_SURFIX="${OUT_SURFIX:-rmPurity}"
CORES_GSEA="${CORES_GSEA:-5}"
for cell in "${CANCERS[@]}"; do
  in_file="output/${DST}/pcc_mRNA/${cell}_${METHOD}_${OUT_SURFIX}.txt"
  if [ ! -s "${in_file}" ]; then
    echo "missing ${in_file}, skip ${cell}"
    continue
  fi
  submit "${cell}_${DST}_gsea" "${PROJECT_DIR}/tmp/${DST}_${cell}_gsea.outerr" "
    $(activate_cmd "${ENV_R}")
    ${SHARED_UTILS_DIR}/run-SncPccTab2GseaGmt.R \
      --input ${in_file} \
      --dataset ${DST} --cell ${cell} --method ${METHOD} --cores ${CORES_GSEA} \
      --PathSource CHG,ImmReg,KEGG,MSigDB --sncRES_cutoff 0.9 --p_label p.adjust --p_cutoff 0.05 \
      --outPrefix output/${DST}/pcc_mRNA/${OUT_SURFIX}
  "
done
# 4k x 20k matrix: 3 cores, ~1 h per cancer type
