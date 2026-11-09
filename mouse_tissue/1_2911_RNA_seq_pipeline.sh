#!/usr/bin/env bash
#===============================================================================
# mouse tissue RNA-seq pipeline
#   fastp -> HISAT2 -> samtools -> featureCounts -> expression matrix
#
# Inputs : mouse_tissue/samples.tsv, reference paths in config.sh
# Outputs: clean/, bam/, counts/2911_tissue_matrix.csv,
#          counts/sample_info.tsv, counts/analysis_params.tsv
#===============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
DATASET="mouse_tissue"
source "${SCRIPT_DIR}/../config.sh"

require_cmd fastp hisat2 samtools featureCounts awk sort uniq wc

#------------------------------------------------------------------------------
# Helpers
#------------------------------------------------------------------------------
pids=()

wait_all() {
  local p rc=0
  for p in "$@"; do
    if ! wait "$p"; then
      rc=1
    fi
  done
  return "${rc}"
}

flush_jobs() {
  if (( ${#pids[@]} > 0 )); then
    wait_all "${pids[@]}"
    pids=()
  fi
}

#------------------------------------------------------------------------------
# 1. Read and validate the sample sheet
#------------------------------------------------------------------------------
[[ -f "${SAMPLES_TSV}" ]] || die "sample sheet not found: ${SAMPLES_TSV}"
require_configured_path "GENOME_INDEX" "${GENOME_INDEX}"
require_configured_path "GTF" "${GTF}"
[[ -f "${GTF}" ]] || die "GTF not found: ${GTF}"

SAMPLES=()
SAMPLE_GROUPS=()
REPLICATES=()
TISSUES=()
R1_FILES=()
R2_FILES=()

while IFS=$'\t' read -r sample group rep tissue r1 r2 || [[ -n "${sample:-}" ]]; do
  sample="${sample%$'\r'}"
  r2="${r2%$'\r'}"
  [[ -z "${sample}" || "${sample}" == "Sample" ]] && continue
  [[ -n "${group}" && -n "${rep}" && -n "${tissue}" && -n "${r1}" && -n "${r2}" ]] \
    || die "incomplete row in ${SAMPLES_TSV}: ${sample}"
  [[ "${group}" == "2911" || "${group}" == "NC" ]] \
    || die "invalid Group '${group}' for sample ${sample}"
  [[ "${tissue}" == "liver" || "${tissue}" == "kidney" || "${tissue}" == "heart" ]] \
    || die "invalid Tissue '${tissue}' for sample ${sample}"

  # Resolve relative read paths against the dataset directory.
  [[ "${r1}" == /* ]] || r1="${DATA_DIR}/${r1}"
  [[ "${r2}" == /* ]] || r2="${DATA_DIR}/${r2}"
  require_configured_path "R1 for ${sample}" "${r1}"
  require_configured_path "R2 for ${sample}" "${r2}"
  [[ -f "${r1}" ]] || die "R1 not found for ${sample}: ${r1}"
  [[ -f "${r2}" ]] || die "R2 not found for ${sample}: ${r2}"

  SAMPLES+=("${sample}")
  SAMPLE_GROUPS+=("${group}")
  REPLICATES+=("${rep}")
  TISSUES+=("${tissue}")
  R1_FILES+=("${r1}")
  R2_FILES+=("${r2}")
done < "${SAMPLES_TSV}"

(( ${#SAMPLES[@]} > 0 )) || die "no samples found in ${SAMPLES_TSV}"
dup_samples=$(printf '%s\n' "${SAMPLES[@]}" | sort | uniq -d)
[[ -z "${dup_samples}" ]] || die "duplicate sample names: ${dup_samples}"

mkdir -p "${CLEAN_DIR}" "${BAM_DIR}" "${COUNT_DIR}"

#------------------------------------------------------------------------------
# 2. fastp: quality control and filtering (forced for this dataset)
#------------------------------------------------------------------------------
echo ">> [1/4] fastp quality control"
pids=()
for i in "${!SAMPLES[@]}"; do
  s="${SAMPLES[$i]}"
  fastp -w "${THREADS}" \
        -i "${R1_FILES[$i]}" -I "${R2_FILES[$i]}" \
        -o "${CLEAN_DIR}/${s}_R1_clean.fq.gz" \
        -O "${CLEAN_DIR}/${s}_R2_clean.fq.gz" \
        -h "${CLEAN_DIR}/${s}.fastp.html" \
        -j "${CLEAN_DIR}/${s}.fastp.json" &
  pids+=($!)
  if (( ${#pids[@]} >= MAX_JOBS )); then
    flush_jobs
  fi
done
flush_jobs

#------------------------------------------------------------------------------
# 3. HISAT2: alignment to the mouse reference genome
#------------------------------------------------------------------------------
echo ">> [2/4] HISAT2 alignment"
pids=()
for s in "${SAMPLES[@]}"; do
  hisat2 -p "${THREADS}" --dta -x "${GENOME_INDEX}" \
         -1 "${CLEAN_DIR}/${s}_R1_clean.fq.gz" \
         -2 "${CLEAN_DIR}/${s}_R2_clean.fq.gz" \
         -S "${BAM_DIR}/${s}.sam" &
  pids+=($!)
  if (( ${#pids[@]} >= MAX_JOBS )); then
    flush_jobs
  fi
done
flush_jobs

for s in "${SAMPLES[@]}"; do
  [[ -s "${BAM_DIR}/${s}.sam" ]] || die "HISAT2 produced an empty SAM file for ${s}"
done

#------------------------------------------------------------------------------
# 4. SAMtools: SAM -> BAM -> sort -> index
#------------------------------------------------------------------------------
echo ">> [3/4] SAMtools sort & index"
pids=()
for s in "${SAMPLES[@]}"; do
  (
    samtools view -b -o "${BAM_DIR}/${s}.bam" "${BAM_DIR}/${s}.sam"
    samtools sort -@ "${SORT_THREADS}" \
                  -o "${BAM_DIR}/${s}_sorted.bam" \
                  "${BAM_DIR}/${s}.bam"
    samtools index "${BAM_DIR}/${s}_sorted.bam"
    rm -f "${BAM_DIR}/${s}.sam" "${BAM_DIR}/${s}.bam"
  ) &
  pids+=($!)
  if (( ${#pids[@]} >= MAX_JOBS )); then
    flush_jobs
  fi
done
flush_jobs

for s in "${SAMPLES[@]}"; do
  [[ -s "${BAM_DIR}/${s}_sorted.bam" ]] || die "sorted BAM missing for ${s}"
  [[ -f "${BAM_DIR}/${s}_sorted.bam.bai" ]] || die "BAM index missing for ${s}"
done

#------------------------------------------------------------------------------
# 5. featureCounts: gene-level counting and final expression matrix
#------------------------------------------------------------------------------
echo ">> [4/4] featureCounts quantification"
BAMS=()
for s in "${SAMPLES[@]}"; do
  BAMS+=("${BAM_DIR}/${s}_sorted.bam")
done

featureCounts -T "${FEATURECOUNT_THREADS}" \
              -p --countReadPairs \
              -s "${STRANDED}" \
              -t exon -g gene_id \
              -a "${GTF}" \
              -o "${COUNT_DIR}/2911_tissue_counts.txt" \
              "${BAMS[@]}"

[[ -s "${COUNT_DIR}/2911_tissue_counts.txt" ]] || die "featureCounts output is empty"

# Final matrix: header = Geneid,Chr,Start,End,Strand,Length + sample columns,
# comma-separated, with Ensembl gene-ID version suffixes removed.
{
  printf 'Geneid,Chr,Start,End,Strand,Length,%s\n' "$(IFS=,; echo "${SAMPLES[*]}")"
  awk 'BEGIN { OFS="," }
       $1 !~ /^#/ && $1 != "Geneid" {
         sub(/\.[0-9]+$/, "", $1)
         print
       }' "${COUNT_DIR}/2911_tissue_counts.txt"
} > "${MATRIX_FILE}"

n_genes=$(( $(wc -l < "${MATRIX_FILE}") - 1 ))
(( n_genes > 0 )) || die "expression matrix contains no genes: ${MATRIX_FILE}"
if ! awk -F, 'NR > 1 { if (seen[$1]++ > 0) exit 1 }' "${MATRIX_FILE}"; then
  die "duplicate Geneid values after version stripping in ${MATRIX_FILE}"
fi

# Sample information table consumed by the DESeq2 scripts.
{
  printf 'Sample\tGroup\tReplicate\tTissue\n'
  for i in "${!SAMPLES[@]}"; do
    printf '%s\t%s\t%s\t%s\n' \
      "${SAMPLES[$i]}" "${SAMPLE_GROUPS[$i]}" "${REPLICATES[$i]}" "${TISSUES[$i]}"
  done
} > "${SAMPLE_INFO_FILE}"

# Analysis parameters consumed by the DESeq2 scripts.
{
  printf 'Parameter\tValue\n'
  printf 'MIN_COUNT\t%s\n' "${MIN_COUNT}"
  printf 'MIN_SAMPLES\t%s\n' "${MIN_SAMPLES}"
  printf 'PADJ_CUT\t%s\n' "${PADJ_CUT}"
} > "${ANALYSIS_PARAMS_FILE}"

echo ">> Done. Expression matrix: ${MATRIX_FILE}"
