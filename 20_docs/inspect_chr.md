#### Other more specialized approaches: ####

Inspect coverage in bins across chromosomes:
```
# Prepare genome by unzipping, then unwrapping (E. Normandeau Scripts function)
gunzip -c 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz > 03_genome/GCF_9028066
45.1_cgigas_uk_roslin_v1_genomic.fna

fasta_unwrap.py 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap.fna

# Keep chromosomes only (replace the matching string in grep with notation for chromosomes)
grep -E -A1 '^>NC_047' 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap.fna | grep -vE '^--$' - > 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap_chr_only.fna

# Use samtools to index the unwrapped, chr-only reference
samtools faidx ./03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap_chr_only.fna

# Determine window positions per chr (adjust window size as needed)
bedtools makewindows -g ./03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap_chr_only.fna.fai -w 1000000 > 03_genome/windows.bed

# Calculate coverage per window on all sorted bam files in 04_samples
01_scripts/03_bedtools_coverage.sh
#   note: this will output 04_samples/<your_sample>_cov.txt

# Use Rscript to summarize and plot
`01_scripts/plot_chr_cov_per_sample.R`
#   note: this will output 04_samples/<your_sample>_plot_aligns_per_window.pdf

```

