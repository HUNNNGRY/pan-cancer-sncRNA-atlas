#!/usr/bin/env bash
# GO/KEGG/MSigDB over-representation analysis on the per-cancer differential
# tables (output/<DST>/diff_<cancer>_<contrast>_rmPurity.diff).
# Uses the sncRNA-target GMT built by 04_snc_pcc_gsea_gmt.sh.
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

cd "${PROJECT_DIR}"
mkdir -p "output/${DST}/pathway"

OUT_SURFIX="${OUT_SURFIX:-rmPurity}"
# contrasts: YesvsNo = metastasis vs primary, PrimaryvsNormal = tumour vs normal
PVAL1="${PVAL1:-0.01}"; PADJ1="${PADJ1:-0.01}"   # PrimaryvsNormal
PVAL2="${PVAL2:-0.01}"; PADJ2="${PADJ2:-1}"      # YesvsNo

for cell in "${CANCERS[@]}"; do
  for diff_type in YesvsNo PrimaryvsNormal; do
    if [ "${diff_type}" = "PrimaryvsNormal" ]; then
      pval="${PVAL1}"; padj="${PADJ1}"
    else
      pval="${PVAL2}"; padj="${PADJ2}"
    fi
    submit "${DST}.${cell}.${diff_type}" \
      "${PROJECT_DIR}/logs/${DST}.${cell}.${diff_type}.outerr" "
      $(activate_cmd "${ENV_R}")
      ${SHARED_SCRIPTS_DIR}/run-clusterProfiler-sncGmt.R \
        --matrix output/${DST}/diff_${cell}_${diff_type}_${OUT_SURFIX}.diff \
        --cutoffP ${pval} \
        --cutoffPadj ${padj} \
        --pValue 1 \
        --adjustPval BH \
        --qValue 1 \
        --fig TRUE \
        --PathSource CHG,ImmReg,KEGG,MSigDB \
        --GmtLong output/${DST}/pcc_mRNA/${OUT_SURFIX}_${DST}_${cell}.gmtLong \
        --outdir output/${DST}/pathway/${cell}_${diff_type}
    "
  done
done
