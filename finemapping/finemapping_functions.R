# Shared functions for FINEMAP and SuSiE-RSS analyses.

read_summary_statistics <- function(path, add_finemap_placeholders = FALSE) {
  x <- data.table::fread(path, data.table = FALSE)
  required <- c("SNP", "CHROM", "POS", "BETA", "SE")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("Summary-statistics file is missing: ", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(x$SNP)) stop("SNP identifiers must be unique.")
  x$SNP <- as.character(x$SNP)
  x$BETA <- as.numeric(x$BETA)
  x$SE <- as.numeric(x$SE)
  x$Z <- x$BETA / x$SE
  if (add_finemap_placeholders) {
    # FINEMAP requires these columns in its .z schema. They were fixed
    # placeholders in this length-dosage analysis and were not variant inputs.
    x$A1 <- "A"
    x$A2 <- "C"
    x$MAF <- 0.2
  }
  x
}

compute_ld_correlation <- function(genotypes) {
  if (!is.matrix(genotypes)) genotypes <- as.matrix(genotypes)
  storage.mode(genotypes) <- "double"
  if (is.null(rownames(genotypes)) || anyDuplicated(rownames(genotypes))) {
    stop("Genotype matrix rows must have unique variant identifiers.")
  }
  if (ncol(genotypes) < 3L) stop("At least three samples are required for LD.")

  # Mean-impute missing dosages within each variant before calculating Pearson LD.
  if (anyNA(genotypes)) {
    means <- rowMeans(genotypes, na.rm = TRUE)
    if (any(!is.finite(means))) stop("One or more variants have no observed dosages.")
    missing <- which(is.na(genotypes), arr.ind = TRUE)
    genotypes[missing] <- means[missing[, "row"]]
  }
  variances <- apply(genotypes, 1L, stats::var)
  if (any(!is.finite(variances) | variances == 0)) {
    stop("Genotype matrix contains a non-finite or invariant variant.")
  }
  R <- stats::cor(t(genotypes), method = "pearson")
  dimnames(R) <- list(rownames(genotypes), rownames(genotypes))
  R
}

align_finemapping_inputs <- function(summary_stats, genotypes) {
  valid <- is.finite(summary_stats$Z) & summary_stats$SE > 0
  summary_stats <- summary_stats[valid, , drop = FALSE]
  shared <- summary_stats$SNP[summary_stats$SNP %in% rownames(genotypes)]
  if (length(shared) < 2L) stop("Fewer than two variants are shared between inputs.")
  summary_stats <- summary_stats[match(shared, summary_stats$SNP), , drop = FALSE]
  genotypes <- genotypes[shared, , drop = FALSE]
  stopifnot(identical(summary_stats$SNP, rownames(genotypes)))
  list(summary_stats = summary_stats, genotypes = genotypes)
}

validate_ld <- function(R, ids) {
  if (!is.matrix(R) || nrow(R) != ncol(R)) stop("LD matrix must be square.")
  if (!identical(rownames(R), colnames(R))) stop("LD dimnames must be identical.")
  if (!identical(ids, rownames(R))) stop("Summary statistics and LD are misaligned.")
  if (any(!is.finite(R))) stop("LD matrix contains non-finite values.")
  if (!isTRUE(isSymmetric(R, tol = 1e-8))) stop("LD matrix is not symmetric.")
  if (max(abs(diag(R) - 1)) > 1e-8) stop("LD diagonal differs from one.")
  invisible(TRUE)
}

write_finemap_inputs <- function(summary_stats, R, output_dir, n_samples) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    z = file.path(output_dir, "finemap.z"),
    ld = file.path(output_dir, "finemap.ld"),
    snp = file.path(output_dir, "finemap.snp"),
    config = file.path(output_dir, "finemap.config"),
    cred = file.path(output_dir, "finemap.cred"),
    log = file.path(output_dir, "finemap.log"),
    master = file.path(output_dir, "finemap.master")
  )
  z <- data.frame(
    rsid = summary_stats$SNP,
    chromosome = summary_stats$CHROM,
    position = summary_stats$POS,
    allele1 = summary_stats$A1,
    allele2 = summary_stats$A2,
    maf = summary_stats$MAF,
    beta = summary_stats$BETA,
    se = summary_stats$SE
  )
  data.table::fwrite(z, paths$z, sep = " ", quote = FALSE)
  utils::write.table(R, paths$ld, row.names = FALSE, col.names = FALSE,
                     quote = FALSE, sep = " ")
  master <- data.frame(
    z = normalizePath(paths$z), ld = normalizePath(paths$ld),
    snp = normalizePath(paths$snp, mustWork = FALSE),
    config = normalizePath(paths$config, mustWork = FALSE),
    cred = normalizePath(paths$cred, mustWork = FALSE),
    log = normalizePath(paths$log, mustWork = FALSE),
    n_samples = as.integer(n_samples)
  )
  utils::write.table(master, paths$master, sep = ";", row.names = FALSE,
                     quote = FALSE)
  paths
}

parse_named_args <- function(args, required, usage) {
  if (!length(args) || any(args %in% c("-h", "--help"))) {
    cat(usage, "\n")
    quit(status = 0L)
  }
  if (length(args) %% 2L != 0L) stop(usage, call. = FALSE)
  keys <- sub("^--", "", args[seq.int(1L, length(args), 2L)])
  values <- args[seq.int(2L, length(args), 2L)]
  if (anyDuplicated(keys) || any(!keys %in% required)) stop(usage, call. = FALSE)
  out <- as.list(values)
  names(out) <- keys
  if (any(!required %in% names(out))) stop(usage, call. = FALSE)
  out
}
