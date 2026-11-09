#!/usr/bin/env bash
#===============================================================================
# mouse tissue: miR2911 target simulation and 2x2 counts
#
# For each per-tissue DEG table defined in config.sh:
#   down genes vs other tested genes (tested - up - down)
#   mature mRNA extraction -> miRanda -> positive/negative counts
#
# Outputs under mir2911_analysis/<label>/ (label "all" uses mir2911_analysis/):
#   down_genes.txt, up_genes.txt, other_genes.txt
#   down_mrna.fa, other_mrna.fa, mir2911.fa
#   down_miranda.out, other_miranda.out
#   miranda_counts.txt, target_analysis_qc.tsv
# and mir2911_analysis/analysis_dirs.tsv for step 5.
#===============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
DATASET="mouse_tissue"
source "${SCRIPT_DIR}/../config.sh"

require_cmd awk sort seqkit miranda grep wc tr

[[ -f "${MRNA_FA}" ]] || die "mRNA FASTA not found: ${MRNA_FA}
Run step 3 (prepare_mrna.sh) first."
[[ ${#DEG_LABELS[@]} -eq ${#DEG_FILES[@]} ]] \
  || die "DEG_LABELS and DEG_FILES length mismatch in config.sh"

mkdir -p "${TARGET_ROOT}"
MANIFEST="${TARGET_ROOT}/analysis_dirs.tsv"
printf 'Label\tDir\n' > "${MANIFEST}"

count_fasta() {
  local n
  n=$(grep -c '^>' "$1" 2>/dev/null || true)
  echo "${n:-0}"
}

count_targets() {
  local f="$1"
  [[ -s "${f}" ]] || {
    echo 0
    return 0
  }
  awk -F'\t' '/^>>/ { if ($2 != "") print $2 }' "${f}" \
    | sort -u | wc -l | tr -d ' '
}

run_one() {
  local label="$1" deg_file="$2"
  local out
  if [[ "${label}" == "all" ]]; then
    out="${TARGET_ROOT}"
  else
    out="${TARGET_ROOT}/${label}"
  fi
  mkdir -p "${out}"

  [[ -f "${deg_file}" ]] || die "DEG file not found for ${label}: ${deg_file}
Run step 2 (DESeq2) first."

  # Gene lists in one pass. padj=NA is never classified as up or down.
  : > "${out}/down_genes.txt"
  : > "${out}/up_genes.txt"
  : > "${out}/other_genes.txt"
  awk -F, -v p="${PADJ_CUT}" \
      -v down_file="${out}/down_genes.txt" \
      -v up_file="${out}/up_genes.txt" \
      -v other_file="${out}/other_genes.txt" '
    BEGIN { OFS="\t" }
    NR == 1 { next }
    {
      gsub(/\r/, "")
      if (NF < 3) next
      gsub(/"/, "", $1); gsub(/"/, "", $2); gsub(/"/, "", $3)
      sub(/\.[0-9]+$/, "", $1)
      if ($1 == "") next
      tested[$1] = 1
      if (toupper($3) != "NA" && $3 != "" &&
          toupper($2) != "NA" && $2 != "" &&
          $3 ~ /^[0-9.eE+-]+$/ && $2 ~ /^-?[0-9.eE+-]+$/) {
        if (($3 + 0) < p) {
          if (($2 + 0) < 0) {
            down[$1] = 1
          } else if (($2 + 0) > 0) {
            up[$1] = 1
          }
        }
      }
    }
    END {
      for (g in down) print g > down_file
      for (g in up) print g > up_file
      for (g in tested) {
        if (!(g in down) && !(g in up)) print g > other_file
      }
    }' "${deg_file}"

  for f in "${out}/down_genes.txt" "${out}/up_genes.txt" "${out}/other_genes.txt"; do
    [[ -f "${f}" ]] || : > "${f}"
    LC_ALL=C sort -u "${f}" -o "${f}"
  done

  : > "${out}/down_mrna.fa"
  : > "${out}/other_mrna.fa"
  if [[ -s "${out}/down_genes.txt" ]]; then
    seqkit grep -f "${out}/down_genes.txt" "${MRNA_FA}" -o "${out}/down_mrna.fa"
  fi
  if [[ -s "${out}/other_genes.txt" ]]; then
    seqkit grep -f "${out}/other_genes.txt" "${MRNA_FA}" -o "${out}/other_mrna.fa"
  fi

  local n_down n_other p_down=0 p_other=0
  n_down=$(count_fasta "${out}/down_mrna.fa")
  n_other=$(count_fasta "${out}/other_mrna.fa")

  printf '>%s\n%s\n' "${MIR2911_NAME}" "${MIR2911_SEQ}" > "${out}/mir2911.fa"

  # miRanda argument order: query (miRNA) first, reference (mRNA) second.
  if (( n_down > 0 )); then
    miranda "${out}/mir2911.fa" "${out}/down_mrna.fa" \
            "${MIRANDA_OPTS[@]}" -out "${out}/down_miranda.out"
    p_down=$(count_targets "${out}/down_miranda.out")
  else
    : > "${out}/down_miranda.out"
  fi

  if (( n_other > 0 )); then
    miranda "${out}/mir2911.fa" "${out}/other_mrna.fa" \
            "${MIRANDA_OPTS[@]}" -out "${out}/other_miranda.out"
    p_other=$(count_targets "${out}/other_miranda.out")
  else
    : > "${out}/other_miranda.out"
  fi

  (( p_down <= n_down )) \
    || die "target count (${p_down}) exceeds down gene count (${n_down}) for ${label}"
  (( p_other <= n_other )) \
    || die "target count (${p_other}) exceeds other gene count (${n_other}) for ${label}"

  local n12 n22
  n12=$(( n_down - p_down ))
  n22=$(( n_other - p_other ))

  printf 'down_target\tdown_nontarget\tother_target\tother_nontarget\n%s\t%s\t%s\t%s\n' \
         "${p_down}" "${n12}" "${p_other}" "${n22}" > "${out}/miranda_counts.txt"

  {
    printf 'Metric\tValue\n'
    printf 'DEG_file\t%s\n' "${deg_file}"
    printf 'PADJ_CUT\t%s\n' "${PADJ_CUT}"
    printf 'MIRANDA_OPTS\t%s\n' "${MIRANDA_OPTS[*]}"
    printf 'down_genes_total\t%s\n' "$(wc -l < "${out}/down_genes.txt" | tr -d ' ')"
    printf 'down_genes_with_mrna\t%s\n' "${n_down}"
    printf 'down_targets\t%s\n' "${p_down}"
    printf 'down_nontargets\t%s\n' "${n12}"
    printf 'other_genes_total\t%s\n' "$(wc -l < "${out}/other_genes.txt" | tr -d ' ')"
    printf 'other_genes_with_mrna\t%s\n' "${n_other}"
    printf 'other_targets\t%s\n' "${p_other}"
    printf 'other_nontargets\t%s\n' "${n22}"
  } > "${out}/target_analysis_qc.tsv"

  printf '%s\t%s\n' "${label}" "${out}" >> "${MANIFEST}"
  echo ">> ${label}: down ${p_down}/${n_down} targets | other ${p_other}/${n_other} targets"
}

for i in "${!DEG_LABELS[@]}"; do
  run_one "${DEG_LABELS[$i]}" "${DEG_FILES[$i]}"
done

echo ">> Counts written under ${TARGET_ROOT}"
