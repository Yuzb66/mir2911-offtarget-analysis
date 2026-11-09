#===============================================================================
# mouse tissue: Fisher exact test for miR2911 target enrichment
#
# Input : mir2911_analysis/analysis_dirs.tsv (from step 4)
# Output: mir2911_analysis/<label>/mir2911_target_table.csv
#         mir2911_analysis/mir2911_fisher_summary.tsv
#===============================================================================

# Resolve the script directory so the script can be run from anywhere.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
setwd(script_dir)

root <- Sys.getenv("TARGET_ROOT", unset = "mir2911_analysis")
if (!nzchar(root)) root <- "mir2911_analysis"
cli <- commandArgs(trailingOnly = TRUE)
if (length(cli) >= 1 && nzchar(cli[1])) {
  root <- cli[1]
}

manifest_file <- file.path(root, "analysis_dirs.tsv")
if (!file.exists(manifest_file)) {
  stop("Manifest not found: ", manifest_file,
       "\nRun step 4 (target analysis) first.")
}

runs <- read.delim(manifest_file, check.names = FALSE,
                   stringsAsFactors = FALSE)
if (!all(c("Label", "Dir") %in% colnames(runs)) || nrow(runs) == 0) {
  stop("Invalid manifest: ", manifest_file)
}
runs <- runs[!is.na(runs$Label) & !is.na(runs$Dir) &
               runs$Label != "" & runs$Dir != "", , drop = FALSE]
if (nrow(runs) == 0) {
  stop("Manifest has no valid runs: ", manifest_file)
}

summary_rows <- vector("list", nrow(runs))

for (i in seq_len(nrow(runs))) {
  label <- runs$Label[i]
  run_dir <- runs$Dir[i]
  counts_file <- file.path(run_dir, "miranda_counts.txt")
  note <- ""
  n11 <- n12 <- n21 <- n22 <- NA_real_
  odds_ratio <- ci_low <- ci_high <- pvalue <- NA_real_
  rate_down <- rate_other <- NA_real_

  if (!file.exists(counts_file)) {
    note <- "counts file missing"
    warning(label, ": ", note)
  } else {
    cnt <- read.delim(counts_file, check.names = FALSE,
                      stringsAsFactors = FALSE)
    needed <- c("down_target", "down_nontarget",
                "other_target", "other_nontarget")
    if (!all(needed %in% colnames(cnt)) || nrow(cnt) == 0) {
      note <- "invalid counts file"
      warning(label, ": ", note)
    } else {
      n11 <- as.numeric(cnt$down_target[1])
      n12 <- as.numeric(cnt$down_nontarget[1])
      n21 <- as.numeric(cnt$other_target[1])
      n22 <- as.numeric(cnt$other_nontarget[1])
      vals <- c(n11, n12, n21, n22)

      if (anyNA(vals) || any(vals < 0) || any(vals != round(vals))) {
        note <- "non-integer or negative counts"
        warning(label, ": ", note)
        n11 <- n12 <- n21 <- n22 <- NA_real_
      } else {
        tab <- matrix(
          c(n11, n12, n21, n22),
          nrow = 2, byrow = TRUE,
          dimnames = list(
            c("下调基因", "非下调基因(除上调外)"),
            c("是miRNA靶点", "非miRNA靶点")
          )
        )
        tab_margin <- addmargins(tab)
        print(tab_margin)
        write.csv(tab_margin,
                  file.path(run_dir, "mir2911_target_table.csv"),
                  quote = FALSE)

        row_down <- n11 + n12
        row_other <- n21 + n22
        if (row_down > 0 && row_other > 0) {
          rate_down <- 100 * n11 / row_down
          rate_other <- 100 * n21 / row_other
        }

        if (row_down > 0 && row_other > 0) {
          ft <- fisher.test(tab)
          odds_ratio <- unname(ft$estimate)
          ci_low <- ft$conf.int[1]
          ci_high <- ft$conf.int[2]
          pvalue <- ft$p.value

          cat(sprintf("\n[%s] Fisher's exact test\n", label))
          cat(sprintf("odds ratio = %.4f\n", odds_ratio))
          cat(sprintf("95%% CI    = %.4f - %.4f\n", ci_low, ci_high))
          cat(sprintf("p-value   = %.4g\n", pvalue))
          cat(sprintf("target rate: down = %.2f%%, other = %.2f%%\n",
                      rate_down, rate_other))

          if (!is.finite(odds_ratio) || !is.finite(pvalue)) {
            note <- "odds ratio or p-value not finite"
            cat("结论: 统计量不可靠，无法判断。\n")
          } else if (pvalue > 0.05) {
            cat("结论: p > 0.05, 下调基因中 mir2911 靶点比例与背景无显著差异。\n")
          } else if (odds_ratio < 1) {
            cat("结论: p < 0.05 但 OR < 1, 下调基因中 mir2911 靶点显著 depletion。\n")
          } else {
            cat("结论: p < 0.05 且 OR > 1, 下调基因富集 mir2911 靶点，需检查阈值或下调定义。\n")
          }
        } else {
          note <- "empty row; Fisher test skipped"
          cat(sprintf("\n[%s] %s\n", label, note))
        }
      }
    }
  }

  summary_rows[[i]] <- data.frame(
    Label = label,
    down_target = n11,
    down_nontarget = n12,
    other_target = n21,
    other_nontarget = n22,
    odds_ratio = odds_ratio,
    ci_low = ci_low,
    ci_high = ci_high,
    pvalue = pvalue,
    down_target_rate = rate_down,
    other_target_rate = rate_other,
    note = note,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
summary_file <- file.path(root, "mir2911_fisher_summary.tsv")
write.table(summary_df, summary_file, sep = "\t", quote = FALSE,
            row.names = FALSE, na = "NA")
cat("\n>> Combined summary: ", summary_file, "\n", sep = "")
