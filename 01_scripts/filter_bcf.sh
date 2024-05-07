#!/bin/bash

# Global variables
DATAFOLDER="05_genotyping"

# User set variables
NUM_CPU=36
INPUT_BCF="mpileup_calls.bcf"
INDEL_DIST=5      # indel surrounding window size around which to remove SNPs
MISSING_PPN=0.1   # missing data proportion (default: filter if > 10%)
QUAL=20           # minimum quality for SNP at least one indiv
MIN_TOTAL_DP=10   # minimum total depth for site
MIN_AVG_DP=10     # minimum depth when averaged across all samples
MIN_GENO_DP=4     # minimum geno depth or change to missing
MAX_GENO_DP=100   # maximum geno depth or change to missing

## Filters
# 1. Filter near indels
echo "Filtering out variants within $INDEL_DIST bp of an indel"
OUTPUT_FN_1="${INPUT_BCF%.bcf}"_noindel"$INDEL_DIST".bcf
bcftools filter --threads $NUM_CPU -g $INDEL_DIST $DATAFOLDER/$INPUT_BCF \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_1
echo "Completed"

# Reporting
echo "Output BCF now has the following number of variants: "
bcftools view $DATAFOLDER/$OUTPUT_FN_1 --threads $NUM_CPU | grep -vE '^#' - | wc -l


# 2. Filter based on missing
echo "Filtering out variants missing in more than $MISSING_PPN (proportion) individuals"
OUTPUT_FN_2="${OUTPUT_FN_1%.bcf}"_miss"$MISSING_PPN".bcf
bcftools view -i "F_missing < $MISSING_PPN" --threads $NUM_CPU $DATAFOLDER/$OUTPUT_FN_1 \
        -Ob -o $DATAFOLDER/$OUTPUT_FN_2
echo "Completed"

# Reporting
echo "Output BCF now has the following number of variants: "
bcftools view $DATAFOLDER/$OUTPUT_FN_2 --threads $NUM_CPU | grep -vE '^#' - | wc -l

