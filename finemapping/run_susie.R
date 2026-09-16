#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(susieR)
})
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
source(file.path(script_dir, "finemapping_functions.R"))

usage <- paste(
  "Usage: Rscript run_susie.R --summary FILE --genotypes FILE",
  "--output-prefix PREFIX --n-samples N --max-signals L"
)
args <- parse_named_args(
  commandArgs(trailingOnly = TRUE),
  c("summary", "genotypes", "output-prefix", "n-samples", "max-signals"),
  usage
)
summary_path <- args$summary
genotype_path <- args$genotypes
output_prefix <- args$`output-prefix`
n_samples <- as.integer(args$`n-samples`)
max_signals <- as.integer(args$`max-signals`)
if (anyNA(c(n_samples, max_signals)) || n_samples < 1L || max_signals < 1L) {
  stop("N_SAMPLES and MAX_SIGNALS must be positive integers.")
}
dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

stats <- read_summary_statistics(summary_path)
genotypes <- readRDS(genotype_path)
aligned <- align_finemapping_inputs(stats, genotypes)
R <- compute_ld_correlation(aligned$genotypes)
validate_ld(R, aligned$summary_stats$SNP)

fit <- susie_rss(
  z = aligned$summary_stats$Z,
  R = R,
  n = n_samples,
  L = max_signals,
  max_iter = 300,
  estimate_residual_variance = FALSE
)
saveRDS(fit, paste0(output_prefix, ".rds"))
pip <- data.frame(SNP = aligned$summary_stats$SNP, PIP = fit$pip)
pip <- pip[order(pip$PIP, decreasing = TRUE), , drop = FALSE]
fwrite(pip, paste0(output_prefix, "_pip.tsv"), sep = "\t")

diagnostics <- data.frame(
  converged = isTRUE(fit$converged),
  iterations = fit$niter,
  n_variants = nrow(pip),
  n_credible_sets = length(fit$sets$cs)
)
fwrite(diagnostics, paste0(output_prefix, "_diagnostics.tsv"), sep = "\t")
message("SuSiE completed. Results prefix: ", output_prefix)
