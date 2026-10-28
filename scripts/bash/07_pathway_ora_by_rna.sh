#!/usr/bin/env bash
# Same over-representation analysis as 06_pathway_ora.sh, but restricted to one
# sncRNA / genomic-region class at a time (features are "<tx>|<class>|..." keys).
# Writes output/<DST>/pathway/<cancer>_<contrast>_<class>/.
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

cd "${PROJECT_DIR}"
mkdir -p "output/${DST}/pathway" "${PROJECT_DIR}/tmp"

OUT_SURFIX="${OUT_SURFIX:-rmPurity}"
PVAL1="${PVAL1:-0.01}"; PADJ1="${PADJ1:-0.01}"   # PrimaryvsNormal
PVAL2="${PVAL2:-0.01}"; PADJ2="${PADJ2:-1}"      # YesvsNo
RNA_CLASSES=("rRNA" "pri_miRNA" "lncRNA" "mRNA" "piRNA" "snoRNA" "snRNA"
             "srpRNA" "tRNA" "tucpRNA" "Y_RNA"
             "intron" "promoter" "enhancer" "repeats")

for cell in "${CANCERS[@]}"; do
  for diff_type in YesvsNo PrimaryvsNormal; do
    if [ "${diff_type}" = "PrimaryvsNormal" ]; then
      pval="${PVAL1}"; padj="${PADJ1}"
    else
      pval="${PVAL2}"; padj="${PADJ2}"
    fi
    in_file="output/${DST}/diff_${cell}_${diff_type}_${OUT_SURFIX}.diff"
    if [ ! -s "${in_file}" ]; then
      echo "missing ${in_file}, skip ${cell} ${diff_type}"
      continue
    fi
    for rna in "${RNA_CLASSES[@]}"; do
      # keep the header plus the rows of this class; class is field 2 of the key
      sub_file="tmp/tmp_${DST}_${cell}_${diff_type}_${OUT_SURFIX}_${rna}.diff"
      awk -v RNA="${rna}" -F'\t' '{
                split($1, parts, "|");
                if (parts[1] == "gene_id" || parts[1] == "feature" || parts[2] == RNA) {
                  print $0;
                }
              }' "${in_file}" > "${sub_file}"
      submit "${DST}.${cell}.${diff_type}.${rna}" \
        "${PROJECT_DIR}/logs/${DST}.${cell}.${diff_type}.${rna}.outerr" "
        $(activate_cmd "${ENV_R}")
        ${SHARED_SCRIPTS_DIR}/run-clusterProfiler-sncGmt.R \
          --matrix ${sub_file} \
          --cutoffP ${pval} \
          --cutoffPadj ${padj} \
          --pValue 1 \
          --adjustPval BH \
          --qValue 1 \
          --fig TRUE \
          --PathSource CHG,ImmReg,KEGG,MSigDB \
          --GmtLong output/${DST}/pcc_mRNA/${OUT_SURFIX}_${DST}_${cell}.gmtLong \
          --outdir output/${DST}/pathway/${cell}_${diff_type}_${rna}
      "
    done
  done
done
