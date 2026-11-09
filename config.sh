#!/usr/bin/env bash
#===============================================================================
# Central configuration for both CMV-MIR2911 RNA-seq pipelines.
#
# Edit ONLY this file and <dataset>/samples.tsv before running.
# Every path containing TODO_SET_ME must be replaced with a real path.
#===============================================================================

# Project root = directory containing this config.sh.
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#------------------------------------------------------------------------------
# Shared parameters
#------------------------------------------------------------------------------
THREADS=8                  # threads per fastp / HISAT2 job
MAX_JOBS=6                 # maximum number of samples processed in parallel
SORT_THREADS=4             # threads per samtools sort
FEATURECOUNT_THREADS=24    # threads for featureCounts

STRANDED=0                 # 0=unstranded, 1=forward, 2=reverse
MIN_COUNT=5                # DESeq2 low-count filter: minimum count
MIN_SAMPLES=3              # DESeq2 low-count filter: minimum number of samples
PADJ_CUT=0.05              # padj cutoff used to define down/up gene lists

# miRanda parameters.
# NOTE: miRanda parses -en with atoi(), so -2.0 is effectively -2.
# "Strict" is miRanda's default; there is no -strict command-line flag.
MIRANDA_SCORE=160
MIRANDA_ENERGY=-2.0
MIRANDA_OPTS=(-sc "${MIRANDA_SCORE}" -en "${MIRANDA_ENERGY}")

# miR2911 mature sequence.
MIR2911_NAME="mir2911"
MIR2911_SEQ="GGCCGGGGGACGGGCTGGGA"

# Ensembl cDNA release used by 3_prepare_mrna.sh.
# Use a release number for reproducibility, or "current" for current_fasta.
ENSEMBL_RELEASE=113

#------------------------------------------------------------------------------
# Helpers shared by all shell scripts
#------------------------------------------------------------------------------
die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_cmd() {
  local c
  for c in "$@"; do
    command -v "$c" >/dev/null 2>&1 || die "required command not found: $c"
  done
}

require_configured_path() {
  local name="$1" value="$2"
  if [[ "${value}" == *TODO_SET_ME* ]]; then
    die "${name} is not configured yet: ${value}
       Edit ${PROJECT_ROOT}/config.sh and the dataset samples.tsv first."
  fi
}

#------------------------------------------------------------------------------
# Dataset-specific settings
#------------------------------------------------------------------------------
: "${DATASET:?Set DATASET=HEK293T or DATASET=mouse_tissue before sourcing config.sh}"

case "${DATASET}" in
  HEK293T)
    DATA_DIR="${PROJECT_ROOT}/HEK293T"
    SAMPLES_TSV="${DATA_DIR}/samples.tsv"

    # TODO: replace with the real human HISAT2 index prefix and GTF.
    GENOME_INDEX="TODO_SET_ME/human/Homo_sapiens"
    GTF="TODO_SET_ME/human/Homo_sapiens.gtf"

    SPECIES="human"
    SPECIES_FA_DIR="homo_sapiens"
    ENSEMBL_FA_PREFIX="Homo_sapiens"
    GENOME_ASSEMBLY="GRCh38"

    CLEAN_DIR="${DATA_DIR}/clean"
    BAM_DIR="${DATA_DIR}/bam"
    COUNT_DIR="${DATA_DIR}/counts"
    MATRIX_FILE="${COUNT_DIR}/2911_matrix.csv"
    SAMPLE_INFO_FILE="${COUNT_DIR}/sample_info.tsv"
    ANALYSIS_PARAMS_FILE="${COUNT_DIR}/analysis_params.tsv"

    MRNA_DIR="${DATA_DIR}/mrna_db"
    MRNA_FA="${MRNA_DIR}/mrna.fa"
    TARGET_ROOT="${DATA_DIR}/mir2911_analysis"

    # One target analysis: 2911 vs NC.
    DEG_LABELS=(all)
    DEG_FILES=("${DATA_DIR}/deg_2911_vs_NC.csv")
    ;;

  mouse_tissue)
    DATA_DIR="${PROJECT_ROOT}/mouse_tissue"
    SAMPLES_TSV="${DATA_DIR}/samples.tsv"

    # TODO: replace with the real mouse HISAT2 index prefix and GTF.
    GENOME_INDEX="TODO_SET_ME/mouse/genome"
    GTF="TODO_SET_ME/mouse/Mus_musculus.gtf"

    SPECIES="mouse"
    SPECIES_FA_DIR="mus_musculus"
    ENSEMBL_FA_PREFIX="Mus_musculus"
    GENOME_ASSEMBLY="GRCm39"

    CLEAN_DIR="${DATA_DIR}/clean"
    BAM_DIR="${DATA_DIR}/bam"
    COUNT_DIR="${DATA_DIR}/counts"
    MATRIX_FILE="${COUNT_DIR}/2911_tissue_matrix.csv"
    SAMPLE_INFO_FILE="${COUNT_DIR}/sample_info.tsv"
    ANALYSIS_PARAMS_FILE="${COUNT_DIR}/analysis_params.tsv"

    MRNA_DIR="${DATA_DIR}/mrna_db"
    MRNA_FA="${MRNA_DIR}/mrna.fa"
    TARGET_ROOT="${DATA_DIR}/mir2911_analysis"

    # Per-tissue target analysis (three independent Fisher tests).
    DEG_LABELS=(liver kidney heart)
    DEG_FILES=(
      "${DATA_DIR}/deg_2911_vs_NC_liver.csv"
      "${DATA_DIR}/deg_2911_vs_NC_kidney.csv"
      "${DATA_DIR}/deg_2911_vs_NC_heart.csv"
    )
    ;;

  *)
    die "unknown DATASET '${DATASET}' (expected HEK293T or mouse_tissue)"
    ;;
esac

# Ensembl cDNA URL derived from the dataset settings above.
if [[ "${ENSEMBL_RELEASE}" == "current" ]]; then
  ENSEMBL_FA_BASE="https://ftp.ensembl.org/pub/current_fasta"
else
  ENSEMBL_FA_BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta"
fi
ENSEMBL_CDNA_URL="${ENSEMBL_FA_BASE}/${SPECIES_FA_DIR}/cdna/${ENSEMBL_FA_PREFIX}.${GENOME_ASSEMBLY}.cdna.all.fa.gz"
