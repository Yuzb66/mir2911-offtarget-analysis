#===============================================================================
# Differential expression analysis: mouse tissues, 2911 vs NC
#   Design: 2 groups x 3 replicates x 3 tissues (liver, kidney, heart)
#
# Input : counts/2911_tissue_matrix.csv
#         counts/sample_info.tsv
#         counts/analysis_params.tsv
# Output: deg_2911_vs_NC_tissues.csv       (group effect adjusted for tissue)
#         deg_2911_vs_NC_<tissue>.csv      (per-tissue analysis)
#         volcano_2911_vs_NC_tissues.png
#         volcano_2911_vs_NC_<tissue>.png
#===============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 is not installed. Install it before running this script.")
  }
  library(DESeq2)
})

# Resolve the script directory so the script can be run from anywhere.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
setwd(script_dir)

env_or <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}
matrix_file <- env_or("MATRIX_FILE", file.path("counts", "2911_tissue_matrix.csv"))
sample_info_file <- env_or("SAMPLE_INFO_FILE", file.path("counts", "sample_info.tsv"))
params_file <- env_or("ANALYSIS_PARAMS_FILE", file.path("counts", "analysis_params.tsv"))

for (f in c(matrix_file, sample_info_file)) {
  if (!file.exists(f)) {
    stop("Input file not found: ", f, "\nRun step 1 (RNA-seq pipeline) first.")
  }
}

#------------------------------------------------------------------------------
# 1. Load count matrix
#------------------------------------------------------------------------------
raw <- read.csv(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
anno_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
missing_anno <- setdiff(anno_cols, colnames(raw))
if (length(missing_anno) > 0) {
  stop("Count matrix is missing annotation columns: ",
       paste(missing_anno, collapse = ", "))
}
if (anyDuplicated(colnames(raw))) {
  stop("Duplicated column names in ", matrix_file)
}
counts <- raw[, setdiff(colnames(raw), anno_cols), drop = FALSE]
if (ncol(counts) == 0) {
  stop("No sample columns found in ", matrix_file)
}
if (anyDuplicated(raw$Geneid)) {
  stop("Duplicated Geneid values in ", matrix_file)
}
if (anyNA(raw$Geneid) || any(raw$Geneid == "")) {
  stop("Empty Geneid values in ", matrix_file)
}
if (!all(vapply(counts, is.numeric, logical(1)))) {
  stop("Non-numeric count columns detected. ",
       "The featureCounts header row may not have been removed correctly.")
}
counts <- as.matrix(counts)
if (anyNA(counts) || any(counts < 0)) {
  stop("Count matrix contains NA or negative values")
}
if (any(counts != round(counts))) {
  stop("Count matrix contains non-integer values")
}
storage.mode(counts) <- "integer"
rownames(counts) <- raw$Geneid

#------------------------------------------------------------------------------
# 2. Sample metadata
#------------------------------------------------------------------------------
sinfo <- read.delim(sample_info_file, check.names = FALSE,
                    stringsAsFactors = FALSE)
required <- c("Sample", "Group", "Replicate", "Tissue")
if (!all(required %in% colnames(sinfo))) {
  stop("sample_info.tsv must contain columns: ",
       paste(required, collapse = ", "))
}
if (anyDuplicated(sinfo$Sample)) {
  stop("Duplicated sample names in sample_info.tsv")
}
if (ncol(counts) != nrow(sinfo) ||
    !setequal(sinfo$Sample, colnames(counts))) {
  stop("Sample names in sample_info.tsv do not match the count matrix columns")
}
if (!all(sinfo$Group %in% c("NC", "2911"))) {
  stop("sample_info.tsv Group values must be NC or 2911")
}
if (!all(sinfo$Tissue %in% c("liver", "kidney", "heart"))) {
  stop("sample_info.tsv Tissue values must be liver, kidney, or heart")
}
if (length(unique(sinfo$Group)) < 2) {
  stop("Both NC and 2911 groups are required")
}
counts <- counts[, sinfo$Sample, drop = FALSE]
tissue_levels <- unique(sinfo$Tissue)
colData <- data.frame(
  row.names = sinfo$Sample,
  Group = factor(sinfo$Group, levels = c("NC", "2911")),
  Replicate = factor(sinfo$Replicate),
  Tissue = factor(sinfo$Tissue, levels = tissue_levels)
)

#------------------------------------------------------------------------------
# 3. Analysis parameters
#------------------------------------------------------------------------------
min_count <- 5L
min_samples <- 3L
if (file.exists(params_file)) {
  params <- read.delim(params_file, check.names = FALSE,
                       stringsAsFactors = FALSE)
  if (all(c("Parameter", "Value") %in% colnames(params))) {
    get_param <- function(name, default) {
      hit <- params$Value[params$Parameter == name]
      if (length(hit) == 1) as.integer(hit) else default
    }
    min_count <- get_param("MIN_COUNT", min_count)
    min_samples <- get_param("MIN_SAMPLES", min_samples)
  }
}
if (is.na(min_count) || min_count < 0) {
  stop("Invalid MIN_COUNT in analysis_params.tsv")
}
if (is.na(min_samples) || min_samples < 1 || min_samples > ncol(counts)) {
  stop("Invalid MIN_SAMPLES in analysis_params.tsv")
}

export_deg <- function(res, file) {
  deg <- as.data.frame(res)
  deg$gene_id <- rownames(deg)
  deg <- deg[, c("gene_id", "log2FoldChange", "padj", "baseMean", "pvalue")]
  colnames(deg) <- c("gene_id", "log2FC", "padj", "baseMean", "pvalue")
  deg <- deg[order(deg$padj, na.last = TRUE), ]
  write.csv(deg, file, row.names = FALSE, quote = FALSE)
  invisible(deg)
}

#------------------------------------------------------------------------------
# Volcano plot helper
#------------------------------------------------------------------------------
plot_volcano <- function(res, file, title, padj_cut = 0.05) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to draw volcano plots.")
  }
  df <- data.frame(
    gene_id = rownames(res),
    log2FC = res$log2FoldChange,
    padj = res$padj,
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$log2FC) & is.finite(df$padj), , drop = FALSE]
  if (nrow(df) == 0) {
    warning("No finite log2FC/padj values for volcano plot: ", title)
    return(invisible(NULL))
  }
  df$neglog10padj <- -log10(pmax(df$padj, 1e-300))
  df$group <- "NS"
  df$group[df$padj < padj_cut & df$log2FC > 0] <- "Up"
  df$group[df$padj < padj_cut & df$log2FC < 0] <- "Down"
  df$group <- factor(df$group, levels = c("Down", "NS", "Up"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FC, y = neglog10padj,
                                        color = group)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.3) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut),
                        linetype = "dashed", color = "grey40") +
    ggplot2::scale_color_manual(
      values = c(Down = "#2c7bb6", NS = "grey70", Up = "#d7191c"),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title,
      x = expression(log[2]~"fold change"),
      y = expression(-log[10]~"adjusted p-value"),
      color = "Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")

  ggplot2::ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
  invisible(p)
}

#------------------------------------------------------------------------------
# 4. Main model: group effect adjusted for tissue
#------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ Tissue + Group)
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]
if (nrow(dds) == 0) {
  stop("No genes pass the low-count filter (MIN_COUNT=", min_count,
       ", MIN_SAMPLES=", min_samples, ")")
}
dds <- DESeq(dds)

res <- results(dds, contrast = c("Group", "2911", "NC"))
# Optional shrinkage:
# res <- lfcShrink(dds, coef = "Group_2911_vs_NC", res = res, type = "apeglm")
summary(res)
export_deg(res, "deg_2911_vs_NC_tissues.csv")
cat(sprintf(">> Tissue-adjusted DEGs (padj<0.05): %d\n",
            sum(res$padj < 0.05, na.rm = TRUE)))

tryCatch(
  plot_volcano(res, "volcano_2911_vs_NC_tissues.png",
               "Mouse (tissue-adjusted): 2911 vs NC"),
  error = function(e) warning("Volcano plot skipped: ", conditionMessage(e))
)

#------------------------------------------------------------------------------
# 5. Per-tissue analysis (used by the downstream target analysis)
#------------------------------------------------------------------------------
for (tis in tissue_levels) {
  keep_samples <- colData$Tissue == tis
  if (sum(keep_samples) == 0) {
    next
  }
  if (min_samples > sum(keep_samples)) {
    stop("MIN_SAMPLES=", min_samples, " exceeds sample count for tissue ", tis)
  }

  cd <- droplevels(colData[keep_samples, , drop = FALSE])
  dds_t <- DESeqDataSetFromMatrix(
    countData = counts[, keep_samples, drop = FALSE],
    colData = cd,
    design = ~ Group
  )
  keep_genes <- rowSums(counts(dds_t) >= min_count) >= min_samples
  dds_t <- dds_t[keep_genes, ]
  if (nrow(dds_t) == 0) {
    stop("No genes pass the low-count filter for tissue ", tis)
  }
  dds_t <- DESeq(dds_t)
  res_t <- results(dds_t, contrast = c("Group", "2911", "NC"))
  export_deg(res_t, sprintf("deg_2911_vs_NC_%s.csv", tis))
  cat(sprintf(">> %s: %d DEGs (padj<0.05)\n",
              tis, sum(res_t$padj < 0.05, na.rm = TRUE)))
  tryCatch(
    plot_volcano(res_t, sprintf("volcano_2911_vs_NC_%s.png", tis),
                 sprintf("Mouse %s: 2911 vs NC", tis)),
    error = function(e) warning("Volcano plot skipped: ", conditionMessage(e))
  )
}

# Volcano plots are written inside the main model and per-tissue loop above.
