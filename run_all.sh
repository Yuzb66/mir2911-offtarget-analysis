#!/usr/bin/env bash
#===============================================================================
# Run one complete CMV-MIR2911 pipeline in order.
#
# Usage:
#   bash run_all.sh HEK293T
#   bash run_all.sh mouse_tissue
#===============================================================================
set -euo pipefail

DATASET="${1:-}"
case "${DATASET}" in
  HEK293T)
    PIPELINE_SCRIPT="1_293T_RNA_seq_pipeline.sh"
    DESEQ_SCRIPT="2_deseq_HEK293T.R"
    ;;
  mouse_tissue)
    PIPELINE_SCRIPT="1_2911_RNA_seq_pipeline.sh"
    DESEQ_SCRIPT="2_deseq_mouse_tissue.R"
    ;;
  *)
    echo "Usage: bash run_all.sh {HEK293T|mouse_tissue}" >&2
    exit 1
    ;;
esac

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${ROOT}/config.sh"

# Export pipeline paths so the R scripts read the same locations as config.sh.
export MATRIX_FILE
export SAMPLE_INFO_FILE
export ANALYSIS_PARAMS_FILE
export TARGET_ROOT

WORKDIR="${DATA_DIR}"
cd "${WORKDIR}"

echo ">> Running ${DATASET} pipeline in ${WORKDIR}"
bash "${PIPELINE_SCRIPT}"
Rscript "${DESEQ_SCRIPT}"
bash 3_prepare_mrna.sh
bash 4_mir2911_target_analysis.sh
Rscript 5_mir2911_fisher.R
echo ">> ${DATASET} pipeline finished"
