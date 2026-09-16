#!/usr/bin/env Rscript

# Logistic STR length GWAS: disease_status ~ allele_sum + PC1...PC10 + sex
# allele_sum is the sum of the two repeat-copy-number alleles.

suppressPackageStartupMessages(library(data.table))

`%||%` <- function(x, y) if (is.null(x)) y else x

usage <- function() paste(
  "Usage: Rscript run_str_length_gwas.R --genotypes FILE --samples FILE",
  "--covariates FILE --output FILE [--chunk-size 100]\n",
  "Genotypes: headerless TSV; locus then two allele-length columns per sample.\n",
  "Samples: one ID per line, ordered as genotype columns.\n",
  "Covariates: headered TSV with IID, disease_status, sex and PC1-PC10."
)

parse_args <- function(x) {
  allowed <- c("--genotypes", "--samples", "--covariates", "--output", "--chunk-size")
  if (!length(x) || any(x %in% c("-h", "--help"))) {
    cat(usage(), "\n"); quit(status = 0L)
  }
  if (length(x) %% 2L || any(!x[seq(1L, length(x), 2L)] %in% allowed)) {
    stop(usage(), call. = FALSE)
  }
  keys <- x[seq(1L, length(x), 2L)]
  if (anyDuplicated(keys)) stop("Each option may be supplied only once.", call. = FALSE)
  ans <- as.list(x[seq(2L, length(x), 2L)])
  names(ans) <- sub("^--", "", keys)
  required <- c("genotypes", "samples", "covariates", "output")
  if (any(!required %in% names(ans))) stop(usage(), call. = FALSE)
  ans$chunk_size <- as.integer(ans$chunk_size %||% "100")
  if (is.na(ans$chunk_size) || ans$chunk_size < 1L) {
    stop("--chunk-size must be a positive integer.", call. = FALSE)
  }
  ans
}

fit_locus <- function(locus_id, allele1, allele2, covar, model_formula) {
  dosage <- allele1 + allele2
  keep <- !is.na(dosage)
  dat <- copy(covar[keep])
  dat[, allele_sum := dosage[keep]]
  n <- nrow(dat)
  n_cases <- sum(dat$disease_status == 1L)
  n_controls <- sum(dat$disease_status == 0L)

  empty <- function(status) data.table(
    locus = locus_id, beta = NA_real_, std_error = NA_real_,
    z_value = NA_real_, p_value = NA_real_, odds_ratio = NA_real_,
    n = n, n_cases = n_cases, n_controls = n_controls,
    call_rate = mean(keep), status = status
  )
  if (!n) return(empty("no_complete_genotypes"))
  if (uniqueN(dat$allele_sum) < 2L) return(empty("invariant"))
  if (uniqueN(dat$disease_status) < 2L) return(empty("single_outcome_class"))

  fit <- tryCatch(
    suppressWarnings(glm(model_formula, data = dat, family = binomial("logit"))),
    error = identity
  )
  if (inherits(fit, "error")) return(empty("model_error"))
  tab <- summary(fit)$coefficients
  if (!"allele_sum" %in% rownames(tab)) return(empty("coefficient_missing"))

  beta <- unname(tab["allele_sum", "Estimate"])
  se <- unname(tab["allele_sum", "Std. Error"])
  z <- unname(tab["allele_sum", "z value"])
  p <- unname(tab["allele_sum", "Pr(>|z|)"])
  status <- if (isTRUE(fit$converged) && all(is.finite(c(beta, se, z, p)))) {
    "ok"
  } else "nonconverged_or_nonfinite"

  data.table(
    locus = locus_id, beta = beta, std_error = se, z_value = z,
    p_value = p, odds_ratio = exp(beta), n = stats::nobs(fit),
    n_cases = n_cases, n_controls = n_controls,
    call_rate = mean(keep), status = status
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
for (path in unlist(args[c("genotypes", "samples", "covariates")])) {
  if (!file.exists(path)) stop("Input not found: ", path, call. = FALSE)
}

sample_ids <- fread(args$samples, header = FALSE, col.names = "IID")$IID
if (!length(sample_ids) || anyNA(sample_ids) || anyDuplicated(sample_ids)) {
  stop("Sample IDs must be unique and non-missing.", call. = FALSE)
}

covar <- fread(args$covariates)
model_covariates <- c(paste0("PC", 1:10), "sex")
required <- c("IID", "disease_status", model_covariates)
missing <- setdiff(required, names(covar))
if (length(missing)) stop("Missing covariate columns: ", paste(missing, collapse = ", "))
if (anyDuplicated(covar$IID)) stop("IID must be unique in covariates.", call. = FALSE)
if (!all(na.omit(unique(covar$disease_status)) %in% c(0, 1))) {
  stop("disease_status must be coded 0=control, 1=case.", call. = FALSE)
}

covar <- covar[match(sample_ids, IID)]
if (anyNA(covar$IID)) stop("Some samples lack covariate records.", call. = FALSE)
complete_covar <- complete.cases(covar[, ..required])
if (!all(complete_covar)) {
  message("Excluding ", sum(!complete_covar), " samples with missing covariates.")
}
covar <- covar[complete_covar]
sample_indices <- which(complete_covar)
model_formula <- reformulate(c("allele_sum", model_covariates), "disease_status")

output_columns <- c(
  "locus", "beta", "std_error", "z_value", "p_value", "odds_ratio",
  "n", "n_cases", "n_controls", "call_rate", "status"
)
empty_output <- as.data.table(setNames(replicate(length(output_columns), logical(0),
                                                 simplify = FALSE), output_columns))
fwrite(empty_output, args$output, sep = "\t")

expected_columns <- 1L + 2L * length(sample_ids)
con <- if (grepl("\\.gz$", args$genotypes)) gzfile(args$genotypes, "rt") else file(args$genotypes, "rt")
on.exit(close(con), add = TRUE)
processed <- 0L

repeat {
  lines <- readLines(con, n = args$chunk_size, warn = FALSE)
  if (!length(lines)) break
  chunk <- fread(text = lines, header = FALSE, sep = "\t", na.strings = "NA",
                 showProgress = FALSE)
  if (ncol(chunk) != expected_columns) {
    stop("Expected ", expected_columns, " columns but found ", ncol(chunk),
         " near row ", processed + 1L, ".", call. = FALSE)
  }
  results <- lapply(seq_len(nrow(chunk)), function(i) {
    values <- as.numeric(chunk[i, -1L])
    fit_locus(
      as.character(chunk[[1L]][i]),
      values[2L * sample_indices - 1L],
      values[2L * sample_indices],
      covar, model_formula
    )
  })
  fwrite(rbindlist(results), args$output, sep = "\t", append = TRUE, col.names = FALSE)
  processed <- processed + nrow(chunk)
  message("Processed ", processed, " loci")
}
message("Completed ", processed, " association tests: ", args$output)
