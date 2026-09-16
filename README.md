# STR length association analysis

This code tests autosomal short tandem repeat (STR) length for association with
multiple sclerosis susceptibility.

## Model

For each STR and individual:

`allele_sum = allele_1 + allele_2`

The fitted logistic model is:

`disease_status ~ allele_sum + PC1 + ... + PC10 + sex`

Disease status is coded 0 for controls and 1 for cases. The coefficient is the
log odds ratio per additional repeat unit; `odds_ratio` is `exp(beta)`.
Individuals missing either allele are excluded at that locus. Individuals
missing phenotype or covariates are excluded from all loci.

## Requirements

- R 4.x and the R package `data.table`
- `bcftools`, `awk` and `gzip`
- An ExpansionHunter VCF containing `VARID`, `GT` and `REPCN`


## Inputs and execution

The covariate TSV must contain `IID`, `disease_status`, `sex`, and
`PC1` through `PC10`. Additional columns are allowed. State the coding of
`sex` in the manuscript or repository data dictionary.

```bash
chmod +x prepare_str_length_genotypes.sh

./prepare_str_length_genotypes.sh \
  data/discovery_strs.vcf.gz \
  data/keep.samples \
  data/str_alleles.tsv.gz

bcftools query --samples-file data/keep.samples \
  -l data/discovery_strs.vcf.gz > data/ordered_samples.txt

Rscript run_str_length_gwas.R \
  --genotypes data/str_alleles.tsv.gz \
  --samples data/ordered_samples.txt \
  --covariates data/covariates.tsv \
  --output results/str_length_gwas.tsv \
  --chunk-size 100
```

The sample file supplied to R must have exactly the same sample order as the
genotype columns. The two `bcftools` commands above enforce this.

## Fine-mapping

The `finemapping/` directory contains the corresponding joint SNV–STR
fine-mapping implementations for FINEMAP and SuSiE-RSS. See
[`finemapping/README.md`](finemapping/README.md) for the required inputs,
analysis settings and execution commands.

## Output

| Column | Meaning |
| --- | --- |
| `locus` | STR `VARID` |
| `beta` | Log odds ratio per repeat unit |
| `std_error` | Standard error of beta |
| `z_value` | Wald statistic |
| `p_value` | Two-sided Wald p-value |
| `odds_ratio` | `exp(beta)` |
| `n` | Individuals analysed at the locus |
| `n_cases`, `n_controls` | Included cases and controls |
| `call_rate` | Proportion with two observed alleles |
| `status` | Model diagnostic |

Invariant loci and failed fits are retained with missing statistics and a
diagnostic status rather than being assigned p=1.

## Data governance
- Remove withdrawn participants when creating `keep.samples`.

