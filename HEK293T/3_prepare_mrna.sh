#!/usr/bin/env bash
#===============================================================================
# HEK293T: prepare mature mRNA (spliced cDNA) sequences, one per gene
# (longest transcript).
#
# Input : deg_2911_vs_NC.csv (from step 2)
# Output: mrna_db/mrna.fa (headers = gene IDs, no version suffix)
#===============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
DATASET="HEK293T"
source "${SCRIPT_DIR}/../config.sh"

require_cmd zcat gzip awk seqkit sort cut sed comm wc

mkdir -p "${MRNA_DIR}"
FA_GZ="${MRNA_DIR}/cdna.all.fa.gz"
SOURCE_FILE="${MRNA_DIR}/cdna.source.txt"
REHEADED="${MRNA_DIR}/cdna_rehead.fa"
LONGEST="${MRNA_DIR}/longest_tx.txt"

# The DEG files are produced by step 2.
for f in "${DEG_FILES[@]}"; do
  [[ -f "${f}" ]] || die "DEG file not found: ${f}
Run step 2 (DESeq2) first."
done

# Download Ensembl cDNA (pinned release unless ENSEMBL_RELEASE=current).
# A source marker prevents silently reusing a cDNA file from another release.
use_existing=0
if [[ -s "${FA_GZ}" ]] && gzip -t "${FA_GZ}" 2>/dev/null; then
  if [[ -f "${SOURCE_FILE}" ]]; then
    if [[ "$(cat "${SOURCE_FILE}")" == "${ENSEMBL_CDNA_URL}" ]]; then
      use_existing=1
    fi
  else
    # Assume a user-provided local file; record the URL for future checks.
    use_existing=1
    printf '%s\n' "${ENSEMBL_CDNA_URL}" > "${SOURCE_FILE}"
  fi
fi

if (( use_existing )); then
  echo ">> Using existing cDNA FASTA: ${FA_GZ}"
else
  require_cmd wget
  rm -f "${FA_GZ}"
  wget -c -O "${FA_GZ}" "${ENSEMBL_CDNA_URL}"
  printf '%s\n' "${ENSEMBL_CDNA_URL}" > "${SOURCE_FILE}"
fi
[[ -s "${FA_GZ}" ]] || die "failed to download cDNA FASTA: ${ENSEMBL_CDNA_URL}"

#------------------------------------------------------------------------------
# 1. Rewrite headers to >gene_id|transcript_id (versions stripped)
#    Ensembl header example:
#    >ENST00000450305.2 cdna:... gene:ENSG00000223972.1 gene_biotype:...
#------------------------------------------------------------------------------
zcat "${FA_GZ}" | awk '
  /^>/ {
    g = ""
    if (match($0, /gene:[A-Za-z0-9.]+/)) {
      g = substr($0, RSTART + 5, RLENGTH - 5); sub(/\.[0-9]+$/, "", g)
    }
    t = substr($1, 2); sub(/\.[0-9]+$/, "", t)
    print ">" (g == "" ? t : g "|" t)
    next
  }
  { print }' > "${REHEADED}"

#------------------------------------------------------------------------------
# 2. Keep the longest transcript per gene as the mature-mRNA representative
#------------------------------------------------------------------------------
seqkit fx2tab -n -l "${REHEADED}" | awk -F'\t' '
  { g = $1; sub(/\|.*/, "", g); L = $2 + 0
    if (L > best[g]) { best[g] = L; name[g] = $1 } }
  END { for (g in name) print name[g] }' | LC_ALL=C sort > "${LONGEST}"

[[ -s "${LONGEST}" ]] || die "no transcripts were parsed from ${REHEADED}"
seqkit grep -f "${LONGEST}" "${REHEADED}" \
  | seqkit replace -p '\|.*' -r '' -o "${MRNA_FA}"
[[ -s "${MRNA_FA}" ]] || die "mRNA FASTA is empty: ${MRNA_FA}"

#------------------------------------------------------------------------------
# 3. Verify ID overlap with all DEG tables
#------------------------------------------------------------------------------
: > "${MRNA_DIR}/deg_ids.txt"
for f in "${DEG_FILES[@]}"; do
  tail -n +2 "${f}" \
    | cut -d, -f1 \
    | tr -d '"' \
    | tr -d '\r' \
    | sed 's/\.[0-9][0-9]*$//' \
    | awk 'NF { print }' >> "${MRNA_DIR}/deg_ids.txt"
done
LC_ALL=C sort -u "${MRNA_DIR}/deg_ids.txt" -o "${MRNA_DIR}/deg_ids.txt"
seqkit seq -n -i "${MRNA_FA}" | LC_ALL=C sort -u > "${MRNA_DIR}/mrna_ids.txt"
comm -12 "${MRNA_DIR}/deg_ids.txt" "${MRNA_DIR}/mrna_ids.txt" \
  > "${MRNA_DIR}/matched_ids.txt"

N_DEG=$(wc -l < "${MRNA_DIR}/deg_ids.txt")
N_HIT=$(wc -l < "${MRNA_DIR}/matched_ids.txt")
[[ "${N_DEG}" -gt 0 ]] || die "no gene IDs found in the DEG tables"
[[ "${N_HIT}" -gt 0 ]] || die "no DEG gene IDs matched the cDNA FASTA.
Check the GTF/FASTA assembly, Ensembl release, and gene-ID type."
N_PCT=$(awk -v n="${N_HIT}" -v d="${N_DEG}" 'BEGIN { printf "%.1f", 100 * n / d }')
echo ">> DEG genes: ${N_DEG}, matched in mrna.fa: ${N_HIT} (${N_PCT}%)"
if (( N_HIT * 2 < N_DEG )); then
  echo "WARNING: less than 50% of DEG IDs matched the cDNA FASTA; check reference versions." >&2
fi
echo ">> Mature mRNA fasta ready: ${MRNA_FA}"
