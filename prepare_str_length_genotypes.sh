#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 INPUT.vcf.gz KEEP_SAMPLES.txt OUTPUT.tsv.gz" >&2
  exit 2
fi

input_vcf=$1
keep_samples=$2
output_tsv=$3
[[ -r "$input_vcf" ]] || { echo "Cannot read VCF: $input_vcf" >&2; exit 1; }
[[ -r "$keep_samples" ]] || { echo "Cannot read samples: $keep_samples" >&2; exit 1; }

bcftools query --samples-file "$keep_samples" \
  -f '%VARID[\t%GT\t%REPCN]\n' "$input_vcf" |
awk -F '\t' -v OFS='\t' '
  {
    printf "%s", $1
    for (i = 2; i <= NF; i += 2) {
      if ($i == "./." || $(i + 1) == "." || $(i + 1) == "") {
        printf "%sNA%sNA", OFS, OFS
      } else {
        n = split($(i + 1), allele, "/")
        if (n != 2 || allele[1] == "." || allele[2] == ".")
          printf "%sNA%sNA", OFS, OFS
        else
          printf "%s%s%s%s", OFS, allele[1], OFS, allele[2]
      }
    }
    print ""
  }
' | gzip -c > "$output_tsv"
