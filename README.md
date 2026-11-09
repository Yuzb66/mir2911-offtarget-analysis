# CMV-MIR2911 RNA-seq and miRNA off-target Analysis

This repository contains two related RNA-seq analysis projects that share the same overall design and analysis logic:

| Project | Samples | Design | Main question |
|---|---|---|---|
| `HEK293T` | 6 samples | 2911 vs NC, 3 replicates per group | Differential expression and miR2911 target enrichment in HEK293T cells |
| `mouse_tissue` | 18 samples | 2911 vs NC, 3 replicates per group, 3 tissues per replicate (liver, kidney, heart) | Tissue-adjusted and per-tissue differential expression, followed by per-tissue miR2911 target enrichment |

The common workflow is:

```text
samples.tsv + config.sh
        |
        v
fastp QC/filtering
        |
        v
HISAT2 alignment
        |
        v
SAMtools SAM -> BAM -> sort -> index
        |
        v
featureCounts gene-level quantification
        |
        v
expression matrix + sample metadata + analysis parameters
        |
        v
DESeq2 differential expression + volcano plots
        |
        v
Ensembl cDNA -> longest transcript per gene -> mature mRNA FASTA
        |
        v
miR2911 target simulation with miRanda
        |
        v
2 x 2 target/non-target counts
        |
        v
Fisher exact test
```

---

## 1. Repository layout

```text
CMV-MIR2911GU5/
├── config.sh                         # central configuration for both projects
├── run_all.sh                        # run one complete project in order
├── README.md
├── HEK293T/
│   ├── samples.tsv                   # single sample/input table for HEK293T
│   ├── 1_293T_RNA_seq_pipeline.sh
│   ├── 2_deseq_HEK293T.R
│   ├── 3_prepare_mrna.sh
│   ├── 4_mir2911_target_analysis.sh
│   └── 5_mir2911_fisher.R
└── mouse_tissue/
    ├── samples.tsv                   # single sample/input table for mouse tissues
    ├── 1_2911_RNA_seq_pipeline.sh
    ├── 2_deseq_mouse_tissue.R
    ├── 3_prepare_mrna.sh
    ├── 4_mir2911_target_analysis.sh
    └── 5_mir2911_fisher.R
```

Each numbered script corresponds to one stage:

| Step | HEK293T | mouse_tissue | Purpose |
|---|---|---|---|
| 1 | `1_293T_RNA_seq_pipeline.sh` | `1_2911_RNA_seq_pipeline.sh` | QC, alignment, BAM processing, gene counts, matrix |
| 2 | `2_deseq_HEK293T.R` | `2_deseq_mouse_tissue.R` | DESeq2 differential expression and volcano plots |
| 3 | `3_prepare_mrna.sh` | `3_prepare_mrna.sh` | Build one mature mRNA sequence per gene from Ensembl cDNA |
| 4 | `4_mir2911_target_analysis.sh` | `4_mir2911_target_analysis.sh` | miRanda target simulation and 2 x 2 counts |
| 5 | `5_mir2911_fisher.R` | `5_mir2911_fisher.R` | Fisher exact test and combined summary |

---

## 2. Before running: fill in the TODO placeholders

The repository is intentionally shipped with placeholder paths marked `TODO_SET_ME`. The scripts stop with an error if a required placeholder is still present.

### 2.1 `config.sh`

Edit the following entries for each dataset:

| Variable | Meaning |
|---|---|
| `GENOME_INDEX` | HISAT2 index prefix, e.g. `/path/to/hg38/genome` |
| `GTF` | GTF annotation used by featureCounts |
| `ENSEMBL_RELEASE` | Ensembl cDNA release, default `113`; use `current` for `current_fasta` |
| `SPECIES_FA_DIR`, `ENSEMBL_FA_PREFIX`, `GENOME_ASSEMBLY` | Used to build the Ensembl cDNA URL; change only if the assembly is different |

Shared parameters that may need adjustment:

| Parameter | Default | Meaning |
|---|---|---|
| `THREADS` | `8` | Threads per fastp/HISAT2 job |
| `MAX_JOBS` | `6` | Maximum parallel sample jobs |
| `SORT_THREADS` | `4` | Threads per `samtools sort` |
| `FEATURECOUNT_THREADS` | `24` | Threads for featureCounts |
| `STRANDED` | `0` | `0` unstranded, `1` forward, `2` reverse |
| `MIN_COUNT` | `5` | DESeq2 low-count filter: minimum count |
| `MIN_SAMPLES` | `3` | DESeq2 low-count filter: minimum number of samples |
| `PADJ_CUT` | `0.05` | padj cutoff for down/up gene lists |
| `MIRANDA_SCORE` | `160` | miRanda `-sc` score threshold |
| `MIRANDA_ENERGY` | `-2.0` | miRanda `-en` energy threshold |
| `MIR2911_SEQ` | `GGCCGGGGGACGGGCTGGGA` | miR2911 query sequence |

Note: miRanda parses `-en` with `atoi()`, so `-2.0` is effectively interpreted as `-2`. miRanda has no `-strict` command-line flag; strict duplex heuristics are the default. Use `-loose` only if intentionally disabling them.

### 2.2 `<project>/samples.tsv`

Both sample sheets use the same columns:

```text
Sample  Group  Replicate  Tissue  R1  R2
```

- `Sample`: unique sample name; also used as the matrix column name.
- `Group`: `2911` or `NC`.
- `Replicate`: replicate identifier.
- `Tissue`: `cell` for HEK293T; `liver`, `kidney`, or `heart` for mouse.
- `R1`, `R2`: raw paired-end FASTQ paths.

Paths may be absolute or relative to the project directory. Replace every `TODO_SET_ME/...` entry with the real path. Mouse tissue now always runs fastp, so its `R1`/`R2` columns must point to **raw reads**, not only platform-provided clean reads.

---

## 3. Software requirements

The scripts expect the following tools on `PATH`:

- `bash` (Bash 4.4 or newer is recommended for arrays and `set -u` behavior)
- `fastp`
- `HISAT2`
- `samtools`
- `featureCounts` from Subread; `--countReadPairs` requires Subread 2.0 or newer
- `seqkit`
- `miRanda` (v3.3a or compatible)
- `wget`, `gzip`, `zcat`
- R with:
  - `DESeq2`
  - `ggplot2` (used for volcano plots; it is also a DESeq2 dependency)

The shell scripts check for required commands, and the DESeq2 scripts check for `DESeq2`; volcano plotting checks for `ggplot2` and is allowed to fail with a warning without stopping the target analysis.

---

## 4. Running the pipelines

### 4.1 Run one complete project

```bash
bash run_all.sh HEK293T
bash run_all.sh mouse_tissue
```

`run_all.sh` runs steps 1 through 5 in order. It also sources `config.sh` and exports the resolved matrix, metadata, and target-analysis paths so the R scripts use exactly the same locations as the shell scripts.

### 4.2 Run step by step

HEK293T:

```bash
cd HEK293T
bash 1_293T_RNA_seq_pipeline.sh
Rscript 2_deseq_HEK293T.R
bash 3_prepare_mrna.sh
bash 4_mir2911_target_analysis.sh
Rscript 5_mir2911_fisher.R
```

Mouse tissue:

```bash
cd mouse_tissue
bash 1_2911_RNA_seq_pipeline.sh
Rscript 2_deseq_mouse_tissue.R
bash 3_prepare_mrna.sh
bash 4_mir2911_target_analysis.sh
Rscript 5_mir2911_fisher.R
```

The shell scripts and R scripts resolve their own directory, so they can be invoked from another working directory. The R scripts also accept the following environment variables when run manually:

- Step 2: `MATRIX_FILE`, `SAMPLE_INFO_FILE`, `ANALYSIS_PARAMS_FILE`
- Step 5: `TARGET_ROOT`, or a root directory as the first command-line argument

If these variables are not set, the R scripts use the default paths under `counts/` and `mir2911_analysis/`.

---

## 5. Step 1: fastp, HISAT2, SAMtools, featureCounts

### 5.1 fastp

For every sample in `samples.tsv`:

- Input: `R1` and `R2` raw FASTQ files.
- Output under `clean/`:
  - `<sample>_R1_clean.fq.gz`
  - `<sample>_R2_clean.fq.gz`
  - `<sample>.fastp.html`
  - `<sample>.fastp.json`

Up to `MAX_JOBS` samples run in parallel, with `THREADS` threads per fastp job.

### 5.2 HISAT2

Each sample is aligned to the HISAT2 index:

```bash
hisat2 -p THREADS --dta -x GENOME_INDEX \
       -1 <sample>_R1_clean.fq.gz \
       -2 <sample>_R2_clean.fq.gz \
       -S <sample>.sam
```

`--dta` is used for downstream transcript/gene quantification. SAM files are written under `bam/`.

### 5.3 SAMtools

Each SAM file is converted, sorted, and indexed:

```bash
samtools view -b -o <sample>.bam <sample>.sam
samtools sort -@ SORT_THREADS -o <sample>_sorted.bam <sample>.bam
samtools index <sample>_sorted.bam
```

Intermediate SAM and unsorted BAM files are removed after successful indexing. The script verifies that every sorted BAM and `.bai` index exists before quantification.

### 5.4 featureCounts

Gene-level quantification is performed for all samples in one call:

```bash
featureCounts -T FEATURECOUNT_THREADS \
              -p --countReadPairs \
              -s STRANDED \
              -t exon -g gene_id \
              -a GTF \
              -o counts/gene_counts.txt \
              <sorted BAM files>
```

- `-p --countReadPairs` counts fragments rather than individual reads for paired-end data.
- `-s STRANDED` is `0` by default; change it if the library is stranded.
- `-t exon -g gene_id` assumes an Ensembl-style GTF with `gene_id` attributes.
- Only uniquely assigned reads are counted with the default featureCounts behavior.

### 5.5 Expression matrix

The final matrix is written to:

- HEK293T: `counts/2911_matrix.csv`
- mouse_tissue: `counts/2911_tissue_matrix.csv`

The matrix has the columns:

```text
Geneid,Chr,Start,End,Strand,Length,<sample 1>,<sample 2>,...
```

The featureCounts comment line and header are skipped explicitly, and Ensembl gene-ID version suffixes such as `.5` are removed from `Geneid` so IDs match the versionless Ensembl cDNA FASTA used later.

The pipeline also writes:

- `counts/sample_info.tsv`: `Sample`, `Group`, `Replicate`, `Tissue`
- `counts/analysis_params.tsv`: `MIN_COUNT`, `MIN_SAMPLES`, `PADJ_CUT`

These files are consumed by the DESeq2 scripts.

---

## 6. Step 2: DESeq2 differential expression and volcano plots

### 6.1 Common input validation

Both DESeq2 scripts:

- Read the matrix, sample metadata, and analysis parameters.
- Require the annotation columns `Geneid`, `Chr`, `Start`, `End`, `Strand`, `Length`.
- Drop the annotation columns and keep one count column per sample.
- Check for duplicated gene IDs, duplicated sample names, non-numeric counts, negative counts, and non-integer counts.
- Require the sample names in `sample_info.tsv` to match the matrix columns exactly.
- Filter low-count genes with:

```r
keep <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
```

- Stop with an explicit error if no genes pass the filter or if both groups are not present.

### 6.2 HEK293T

Design:

```r
design = ~ Group
```

Contrast:

```r
results(dds, contrast = c("Group", "2911", "NC"))
```

A positive `log2FC` means higher expression in 2911; a negative `log2FC` means lower expression in 2911.

Output:

- `deg_2911_vs_NC.csv` with columns `gene_id, log2FC, padj, baseMean, pvalue`
- `volcano_2911_vs_NC.png`

### 6.3 Mouse tissue

The mouse script performs two related analyses.

Main model, adjusting for tissue:

```r
design = ~ Tissue + Group
results(dds, contrast = c("Group", "2911", "NC"))
```

Output:

- `deg_2911_vs_NC_tissues.csv` (tissue-adjusted group effect)
- `volcano_2911_vs_NC_tissues.png`

Per-tissue analyses:

```r
design = ~ Group
```

For each of `liver`, `kidney`, and `heart`:

- `deg_2911_vs_NC_liver.csv` / `deg_2911_vs_NC_kidney.csv` / `deg_2911_vs_NC_heart.csv`
- `volcano_2911_vs_NC_liver.png` / `volcano_2911_vs_NC_kidney.png` / `volcano_2911_vs_NC_heart.png`

The per-tissue DEG tables are the ones used by the mouse target analysis. The tissue-adjusted overall table is retained as a complementary result.

### 6.4 DEG table format

All DEG tables are written without quotes so that downstream shell tools can read them directly:

```text
gene_id,log2FC,padj,baseMean,pvalue
ENSG00000223972,-1.23,0.001,152.4,0.0002
...
```

Rows with `padj = NA` are retained in the table but are not classified as up- or down-regulated downstream.

### 6.5 Volcano plots

Volcano plots are drawn with `ggplot2`:

- x-axis: `log2FC`
- y-axis: `-log10(padj)`
- colors:
  - Down: `padj < 0.05` and `log2FC < 0`
  - Up: `padj < 0.05` and `log2FC > 0`
  - NS: all other genes with finite `padj` and `log2FC`
- horizontal dashed line: `padj = 0.05`
- `padj = 0` is capped at `1e-300` for plotting.
- No `|log2FC|` cutoff is applied for coloring; this is consistent with the target-analysis definition, which uses padj and the sign of log2FC.

Plotting failure is reported as a warning and does not stop the downstream target analysis.

---

## 7. Step 3: prepare mature mRNA sequences

### 7.1 Source

The script downloads Ensembl cDNA for the configured species and release:

- Human: `Homo_sapiens.GRCh38.cdna.all.fa.gz`
- Mouse: `Mus_musculus.GRCm39.cdna.all.fa.gz`

The URL is built from `ENSEMBL_RELEASE`, `SPECIES_FA_DIR`, `ENSEMBL_FA_PREFIX`, and `GENOME_ASSEMBLY` in `config.sh`. The downloaded file is cached in `mrna_db/`.

A source marker `mrna_db/cdna.source.txt` records the URL used for the cached file. If the configured release or assembly changes, the marker does not match and the cDNA file is downloaded again instead of silently reusing stale data.

### 7.2 Header rewriting

Ensembl cDNA headers are rewritten from:

```text
>ENST00000450305.2 cdna:... gene:ENSG00000223972.1 gene_biotype:...
```

to:

```text
>ENSG00000223972|ENST00000450305
```

Gene and transcript version suffixes are stripped.

### 7.3 Longest transcript per gene

Only the longest transcript for each gene is retained as the representative mature mRNA. The final FASTA is:

```text
mrna_db/mrna.fa
```

with headers containing only the versionless gene ID, for example:

```text
>ENSG00000223972
```

### 7.4 ID consistency check

The script collects the union of gene IDs from all configured DEG tables, strips quotes and version suffixes, and compares them with the mRNA FASTA IDs. It reports:

- number of DEG gene IDs
- number matched in `mrna.fa`
- percentage matched

The pipeline stops if no IDs match, which usually indicates a species, assembly, release, or gene-ID-type mismatch. A warning is printed if fewer than 50% of DEG IDs match.

The DEG tables used here are:

- HEK293T: `deg_2911_vs_NC.csv`
- mouse_tissue: `deg_2911_vs_NC_liver.csv`, `deg_2911_vs_NC_kidney.csv`, `deg_2911_vs_NC_heart.csv`

---

## 8. Step 4: miR2911 target simulation

### 8.1 Gene sets

For each configured DEG table, the script defines:

- Down-regulated genes:
  - `padj` is not `NA`
  - `padj < PADJ_CUT`
  - `log2FC < 0`
- Up-regulated genes:
  - `padj` is not `NA`
  - `padj < PADJ_CUT`
  - `log2FC > 0`
- Other tested genes:
  - all genes present in the DEG table
  - minus down-regulated genes
  - minus up-regulated genes

Genes with `padj = NA`, missing values, non-numeric values, or a log2FC of exactly zero are not classified as up/down and therefore fall into the “other” set. This avoids the previous double-counting problem and gives a clean 2 x 2 table:

| | Target | Non-target |
|---|---|---|
| Down-regulated | down target | down non-target |
| Other tested genes | other target | other non-target |

### 8.2 Mature mRNA extraction

For each gene set, the corresponding mRNA sequences are extracted from `mrna_db/mrna.fa` with:

```bash
seqkit grep -f genes.txt mrna_db/mrna.fa -o genes_mrna.fa
```

Only genes with an available cDNA sequence can enter the Fisher test. The QC file reports both the total number of genes in the list and the number with an mRNA sequence.

### 8.3 miRanda

The miR2911 query is:

```text
>mir2911
GGCCGGGGGACGGGCTGGGA
```

miRanda is called with the query first and the reference mRNA set second:

```bash
miranda mir2911.fa down_mrna.fa \
        -sc 160 -en -2.0 \
        -out down_miranda.out
```

The same is done for the “other” set.

A gene is counted as a positive target if it has at least one miRanda summary line starting with `>>`, which indicates at least one hit passing the configured score and energy thresholds.

### 8.4 Per-run outputs

HEK293T writes its single run directly under:

```text
HEK293T/mir2911_analysis/
```

Mouse tissue writes one directory per tissue:

```text
mouse_tissue/mir2911_analysis/liver/
mouse_tissue/mir2911_analysis/kidney/
mouse_tissue/mir2911_analysis/heart/
```

Each run directory contains:

```text
down_genes.txt
up_genes.txt
other_genes.txt
down_mrna.fa
other_mrna.fa
mir2911.fa
down_miranda.out
other_miranda.out
miranda_counts.txt
target_analysis_qc.tsv
```

`miranda_counts.txt` has the columns:

```text
down_target  down_nontarget  other_target  other_nontarget
```

The run root also contains `analysis_dirs.tsv`, which lists the label and directory of each completed run for step 5.

Example for mouse:

```text
Label   Dir
liver   mouse_tissue/mir2911_analysis/liver
kidney  mouse_tissue/mir2911_analysis/kidney
heart   mouse_tissue/mir2911_analysis/heart
```

The HEK293T run uses the label `all` and writes directly under `mir2911_analysis/`.

---

## 9. Step 5: Fisher exact test

`5_mir2911_fisher.R` reads `analysis_dirs.tsv` and processes every completed run.

For each run:

1. Read `miranda_counts.txt`.
2. Build the 2 x 2 table:

   ```text
                         target   non-target
   down-regulated        n11      n12
   other tested genes    n21      n22
   ```

3. Print the table with margins.
4. Run a two-sided Fisher exact test.
5. Report:
   - odds ratio
   - 95% confidence interval
   - p-value
   - target rate in down-regulated genes
   - target rate in other tested genes
6. Write the per-run table to `mir2911_target_table.csv`.

A combined summary is written to:

```text
<project>/mir2911_analysis/mir2911_fisher_summary.tsv
```

with one row per run:

```text
Label  down_target  down_nontarget  other_target  other_nontarget
odds_ratio  ci_low  ci_high  pvalue
down_target_rate  other_target_rate  note
```

Runs with missing or invalid counts are retained in the summary with `NA` values and a note. Runs with an empty row (for example, no down-regulated genes) are not tested; zero columns are still passed to `fisher.test()`. For mouse tissue, the three per-tissue tests are reported independently; no multiple-testing correction is applied by default.

---

## 10. Output overview

### 10.1 HEK293T

```text
HEK293T/
├── clean/
│   ├── <sample>_R1_clean.fq.gz
│   ├── <sample>_R2_clean.fq.gz
│   ├── <sample>.fastp.html
│   └── <sample>.fastp.json
├── bam/
│   ├── <sample>_sorted.bam
│   └── <sample>_sorted.bam.bai
├── counts/
│   ├── 2911_matrix.csv
│   ├── gene_counts.txt
│   ├── gene_counts.txt.summary
│   ├── sample_info.tsv
│   └── analysis_params.tsv
├── deg_2911_vs_NC.csv
├── volcano_2911_vs_NC.png
├── mrna_db/
│   ├── cdna.all.fa.gz
│   ├── cdna.source.txt
│   ├── cdna_rehead.fa
│   ├── longest_tx.txt
│   ├── mrna.fa
│   ├── deg_ids.txt
│   ├── mrna_ids.txt
│   └── matched_ids.txt
└── mir2911_analysis/
    ├── analysis_dirs.tsv
    ├── down_genes.txt
    ├── up_genes.txt
    ├── other_genes.txt
    ├── down_mrna.fa
    ├── other_mrna.fa
    ├── mir2911.fa
    ├── down_miranda.out
    ├── other_miranda.out
    ├── miranda_counts.txt
    ├── target_analysis_qc.tsv
    ├── mir2911_target_table.csv
    └── mir2911_fisher_summary.tsv
```

### 10.2 mouse_tissue

```text
mouse_tissue/
├── clean/
├── bam/
├── counts/
│   ├── 2911_tissue_matrix.csv
│   ├── 2911_tissue_counts.txt
│   ├── 2911_tissue_counts.txt.summary
│   ├── sample_info.tsv
│   └── analysis_params.tsv
├── deg_2911_vs_NC_tissues.csv
├── deg_2911_vs_NC_liver.csv
├── deg_2911_vs_NC_kidney.csv
├── deg_2911_vs_NC_heart.csv
├── volcano_2911_vs_NC_tissues.png
├── volcano_2911_vs_NC_liver.png
├── volcano_2911_vs_NC_kidney.png
├── volcano_2911_vs_NC_heart.png
├── mrna_db/
└── mir2911_analysis/
    ├── analysis_dirs.tsv
    ├── liver/
    ├── kidney/
    ├── heart/
    └── mir2911_fisher_summary.tsv
```

---

## 11. Important assumptions and notes

- **Gene IDs must be Ensembl gene IDs.** The matrix builder strips version suffixes such as `.5`; the cDNA FASTA is also versionless. If the GTF uses gene symbols or another ID type, adjust the quantification/annotation step accordingly.
- **Reference consistency.** The HISAT2 index, GTF, and Ensembl cDNA assembly/release should match. A mismatch is usually detected by the ID-overlap check in step 3.
- **Strandedness.** The default is `STRANDED=0` (unstranded). Set `STRANDED=1` or `2` if the library protocol requires it.
- **Mouse fastp input.** Mouse tissue always runs fastp, so `samples.tsv` must provide raw reads. If only platform-provided clean reads are available, either point `R1`/`R2` to those files and adapt the script, or obtain the raw reads.
- **`--countReadPairs`.** This option requires Subread 2.0 or newer. It changes the counts from read-level to fragment-level, which is the intended behavior for paired-end RNA-seq.
- **miRanda thresholds.** `-en` is parsed as an integer by miRanda; `-2.0` is treated as `-2`. The default strict duplex heuristics are used; `-strict` is not a valid miRanda option.
- **Target analysis universe.** Only genes with an available cDNA sequence can be tested. The Fisher test therefore uses the number of genes with mRNA as row totals, not necessarily the full DEG table size. `target_analysis_qc.tsv` reports both numbers.
- **Background definition.** The “other” set is tested genes minus down-regulated genes minus up-regulated genes. Up-regulated genes are not included in the background, and down-regulated genes are not double-counted.
- **`padj = NA`.** These genes are not classified as up/down. They remain in the “other” set.
- **Fisher test.** The test is two-sided. Mouse tissue produces three independent per-tissue tests by default; no multiple-testing correction is applied. Add a correction if required by your analysis plan.
- **Reruns overwrite outputs.** Known output files are overwritten when the corresponding step is rerun. The cDNA cache is reused only when its source marker matches the configured URL.
- **No PCA plots.** The DESeq2 scripts produce volcano plots only. PCA was removed intentionally.

---

## 12. Troubleshooting

| Symptom | Likely cause / action |
|---|---|
| `... is not configured yet: TODO_SET_ME/...` | Edit `config.sh` and/or `<project>/samples.tsv` and replace all placeholders |
| `GTF not found` | Check the `GTF` path in `config.sh` |
| HISAT2 cannot open index files | Check `GENOME_INDEX`; it must be the index prefix, not the `.ht2` file |
| `R1 not found` / `R2 not found` | Check the `R1`/`R2` columns in `samples.tsv`; relative paths are resolved against the project directory |
| featureCounts: `unrecognized option '--countReadPairs'` | Upgrade Subread to 2.0 or newer |
| DESeq2: `Non-numeric count columns detected` | The featureCounts matrix was not built correctly; check the matrix header and rerun step 1 |
| `No genes pass the low-count filter` | Lower `MIN_COUNT`/`MIN_SAMPLES`, or inspect the count matrix and sample metadata |
| `No DEG gene IDs matched the cDNA FASTA` | Check species, assembly, Ensembl release, and gene-ID type; ensure the matrix gene IDs are Ensembl IDs |
| `mRNA FASTA is empty` | Check the cDNA download and the longest-transcript extraction |
| miRanda output has no `>>` lines | No hits passed the thresholds, or the input FASTA is empty; inspect `target_analysis_qc.tsv` and the miRanda output |
| `counts file missing` in step 5 | Run step 4 first, or check `analysis_dirs.tsv` |
| Volcano plot skipped with a ggplot2 message | Install `ggplot2`; DESeq2 normally installs it as a dependency |
| `Rscript: command not found` | Install R or activate the environment containing R |

---

## 13. Quick-start checklist

1. Edit `config.sh`:
   - human/mouse `GENOME_INDEX`
   - human/mouse `GTF`
   - confirm `ENSEMBL_RELEASE`, `STRANDED`, thresholds, and thread counts
2. Edit `HEK293T/samples.tsv` and `mouse_tissue/samples.tsv`:
   - sample names, groups, replicates, tissues
   - raw `R1`/`R2` paths
3. Run:

   ```bash
   bash run_all.sh HEK293T
   bash run_all.sh mouse_tissue
   ```

4. Inspect:
   - DEG tables and volcano plots
   - `mrna_db/matched_ids.txt` and the ID-overlap report
   - `mir2911_analysis/**/target_analysis_qc.tsv`
   - `mir2911_analysis/mir2911_fisher_summary.tsv`
