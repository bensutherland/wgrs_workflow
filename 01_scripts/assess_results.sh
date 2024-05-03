#!/bin/bash
# Run this script from within the wgrs_workflow repo

# Set variables
FQ_SUFFIX="_R1.fastq.gz"
BAM_SUFFIX=".sorted.bam"
POSTDUP_SUFFIX=".sorted_mdups.bam"

# Move into the directory containing files
cd 04_samples 

# # Clear previous runs
# rm reads_per_sample.txt 2> /dev/null
# rm reads_per_sample_table.txt 2> /dev/null
# rm mappings_per_sample.txt 2> /dev/null 
# rm mappings_per_sample_table.txt 2> /dev/null
# rm mappings_post_rmdups_sample.txt 2> /dev/null
# rm mappings_post_rmdups_sample_table.txt 2> /dev/null
# 
# 
# ## Number Reads
# # Determine number of reads per file from samples file
# for i in $(ls *$FQ_SUFFIX) ; 
#     do echo $i ;
#     gunzip -c $i | wc -l | awk 'END { print $1/4}' ;
#     done >> reads_per_sample.txt
# 
# # Separate by second line into two columns
# sed 'N;s/\n/ /' reads_per_sample.txt > reads_per_sample_table.txt
# 
# # Remove intermediate file
# rm reads_per_sample.txt
# 
# 
# ## Number Mappings
# # Determine number of reads per file from samples file
# for i in $(ls *$BAM_SUFFIX) ; 
#     do echo $i ;
#     samtools view $i | wc -l ;
#     done >> mappings_per_sample.txt
# 
# # Separate by second line into two columns
# sed 'N;s/\n/ /' mappings_per_sample.txt > mappings_per_sample_table.txt
# 
# # Remove intermediate file
# rm mappings_per_sample.txt
# 

## Number Mappings, after dups removed
# Determine number of reads per file from samples file
for i in $(ls *$POSTDUP_SUFFIX) ; 
    do echo $i ;
    samtools view $i | wc -l ;
    done >> mappings_post_rmdups_sample.txt

# Separate by second line into two columns
sed 'N;s/\n/ /' mappings_post_rmdups_sample.txt > mappings_post_rmdups_sample_table.txt

# Remove intermediate file
rm mappings_post_rmdups_sample.txt



## Reporting
echo "The assessment is complete. Please note that the alignments comprise both reads, whereas the fastq is only a single read. The pipeline assumes that the number of reads for R1 and R2 is equal."

# Plot number of reads and number of mappings in a barplot
# Will save out to main directory 
#Rscript ./../ms_oyster_popgen/01_scripts/assess_reads_mappings.R

