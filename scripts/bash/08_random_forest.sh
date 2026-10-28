#!/usr/bin/env bash
# Cross-validated random forest classification with TCGA sncRNA profiles.
#   MODE=subtype : LUAD histological subtype, within-LUAD cross-validation
#   MODE=too     : normal tissue-of-origin across 13 normal tissues
# The annotation tables are written by the random-forest section of
# TCGA_sncRNA_cfPeak.R.
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

cd "${PROJECT_DIR}"
MODE="${MODE:-too}"
CELL="${CELL:-LUAD}"
RSCRIPT="${RSCRIPT:-Rscript}"
COUNTS="${DATA_DIR}/output/${DST}/call_peak_all/count_matrix"
RF_SCRIPT="${REPO_DIR}/scripts/random_forest_tissue_classification_cv_v2.R"

case "${MODE}" in
  subtype)
    feature_matrix="${COUNTS}/cfpeakCNN_b5_d50_p1_log2CPM_rmPurity.txt"
    annotation="${DATA_DIR}/exSeek-dev/data/${DST}/sample_table_${CELL}_subtype.txt"
    output="${PROJECT_DIR}/output/${DST}/subtype_${CELL}/"
    folds=10;  feature_selection=5000; trees=500; mtry=50
    ;;
  too)
    feature_matrix="${COUNTS}/cfpeakCNN_b5_d50_p1_CPM_rmPurity.txt"
    annotation="${DATA_DIR}/exSeek-dev/data/${DST}/sample_table_NC_TOO.txt"
    output="${PROJECT_DIR}/output/${DST}/TOO/"
    folds=5;   feature_selection=1000; trees=500; mtry=20
    ;;
  *)
    echo "unknown MODE '${MODE}', use subtype or too" >&2
    exit 1
    ;;
esac

mkdir -p "${output}"
submit "RF_${MODE}_${DST}" "${PROJECT_DIR}/tmp/${DST}_${MODE}_rf.outerr" "
  ${RSCRIPT} ${RF_SCRIPT} \
    --feature_matrix ${feature_matrix} \
    --annotation ${annotation} \
    --trees ${trees} --mtry ${mtry} --folds ${folds} \
    --scaling standard \
    --feature_selection ${feature_selection} \
    --output ${output} \
    --model_output random_forest_model.rds
"
