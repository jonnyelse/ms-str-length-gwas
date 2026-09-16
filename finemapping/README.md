# Joint SNV–STR fine-mapping

These scripts document the FINEMAP and SuSiE-RSS analyses used to fine-map
combined SNV and STR association signals.

## Method

For each locus, SNV and STR marginal association statistics were combined with
an in-sample LD correlation matrix calculated from the corresponding combined
dosage matrix. Rows are variants and columns are samples. Missing dosages are
mean-imputed within variant, after which Pearson correlations are calculated
between variants. Summary statistics and LD are aligned by variant identifier.

Both methods use `BETA` and `SE`, with `Z = BETA / SE`. FINEMAP is run using
stochastic search (`--sss`). SuSiE is run with `susie_rss`, a maximum of 300
iterations and `estimate_residual_variance = FALSE`.

## Inputs

The summary-statistics file is a headered TSV with these columns:

| Column | Meaning |
| --- | --- |
| `SNP` | Unique variant identifier matching the dosage-matrix row name |
| `CHROM` | Chromosome |
| `POS` | GRCh38 position |
| `BETA` | Marginal log-odds coefficient |
| `SE` | Standard error of `BETA` |
FINEMAP also requires allele and MAF columns in its `.z` file. These schema fields are set to fixed placeholders (`A1=A`,
`A2=C`, `MAF=0.2`). 

The genotype input is an RDS file containing a numeric matrix with variants as
rows, samples as columns, and unique variant IDs as row names. It must contain
the same samples used to estimate the association statistics.

## Sample size

`N_SAMPLES` is supplied explicitly. The effective sample size used for both
programs was calculated as:

```r
N_SAMPLES <- round(2 / (1 / N_cases + 1 / N_controls))
```



## Run FINEMAP

FINEMAP v1.4.2 was used. For a locus allowing up to four causal signals:

```bash
Rscript finemapping/run_finemap.R \
  --summary data/region_summary.tsv \
  --genotypes data/region_genotypes.rds \
  --output-dir results/finemap/region_name \
  --n-samples 10342 \
  --max-signals 4 \
  --finemap /path/to/finemap_v1.4.2
```

The output directory contains the FINEMAP `.z`, `.ld`, master, SNP, config,
credible-set and log files.

## Run SuSiE

```bash
Rscript finemapping/run_susie.R \
  --summary data/region_summary.tsv \
  --genotypes data/region_genotypes.rds \
  --output-prefix results/susie/region_name \
  --n-samples 10342 \
  --max-signals 4
```

Outputs are the complete fitted model (`.rds`), a sorted PIP table and a small
diagnostics table. Credible sets remain available in `fit$sets$cs`.

## Locus settings used in the analysis

The maximum number of causal signals is set independently for every region with
`--max-signals`. The same value is passed as FINEMAP's `--n-causal-snps` and
SuSiE's `L`; it is not hard-coded in either script.

The fine-mapping regions were defined using pre-specified LD-windows.
