#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
source(file.path(script_dir, "finemapping_functions.R"))

usage <- paste(
  "Usage: Rscript run_finemap.R --summary FILE --genotypes FILE",
  "--output-dir DIR --n-samples N --max-signals K --finemap EXECUTABLE"
)
args <- parse_named_args(
  commandArgs(trailingOnly = TRUE),
  c("summary", "genotypes", "output-dir", "n-samples", "max-signals", "finemap"),
  usage
)
summary_path <- args$summary
genotype_path <- args$genotypes
output_dir <- args$`output-dir`
n_samples <- as.integer(args$`n-samples`)
max_signals <- as.integer(args$`max-signals`)
finemap_executable <- args$finemap
if (anyNA(c(n_samples, max_signals)) || n_samples < 1L || max_signals < 1L) {
  stop("N_SAMPLES and MAX_SIGNALS must be positive integers.")
}
if (!file.exists(finemap_executable)) stop("FINEMAP executable not found.")

stats <- read_summary_statistics(summary_path, add_finemap_placeholders = TRUE)
genotypes <- readRDS(genotype_path)
aligned <- align_finemapping_inputs(stats, genotypes)
R <- compute_ld_correlation(aligned$genotypes)
validate_ld(R, aligned$summary_stats$SNP)
paths <- write_finemap_inputs(aligned$summary_stats, R, output_dir, n_samples)

status <- system2(
  finemap_executable,
  c("--sss", "--n-causal-snps", max_signals, "--in-files", paths$master)
)
if (status != 0L) stop("FINEMAP exited with status ", status, ".")
message("FINEMAP completed. Results: ", paths$snp)
