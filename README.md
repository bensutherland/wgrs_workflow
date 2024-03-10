# wgrs_workflow
Workflow to support genotyping of moderate/ high depth whole-genome resequence data.       
Developed by Ben J. G. Sutherland, Ph.D. (Sutherland Bioinformatics) in the working group of Timothy J. Green (VIU).     

**Note**: this software is provided 'as is', without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in action of contract, tort or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software.       

The development of this pipeline has been supported by the following organizations: [Support and Funding page](20_docs/funding_support.md).     
Please note: this pipeline was developed and inspired by code from [Eric Normandeau's stacks_workflow](https://github.com/enormandeau/stacks_workflow).      

#### Requirements:       
- Linux or Mac operating system     
- FastQC      
- MultiQC     
- Scripts repo from Eric Normandeau     
- picard tools
- java for picard
- [vcflib](https://github.com/vcflib/vcflib)      

#### Citation ####
If you find this tool useful, please cite the original article that uses the tool:        

Please also be sure to cite the tools applied within each function.      

## Sections ##
Add here (#TODO)    

## Getting started ##
Clone this repository and change into the main directory.      
```
git clone https://github.com/bensutherland/wgrs_workflow.git
cd wgrs_workflow   

```

### Prepare data ###
Copy links to all raw data in `02_raw`.      

### Quality control and trimming ###
View raw data with fastqc and multiqc:     
```
fastqc 02_raw/*.fastq.gz -o 02_raw/fastqc_raw/ -t 5
multiqc -o 02_raw/fastqc_raw/ 02_raw/fastqc_raw
```

Trim adapters and for quality in parallel:     
`./01_scripts/01_cutadapt_PE.sh <num_cores>`     
...this will trim for quality (-q 15), minimum length (-m 50), and terminal Ns. It will also remove adapters that are present in `00_archive/adapters.fasta`.    

View trimmed data with fastqc and multiqc:      
```
fastqc 02_raw/trimmed/*.fastq.gz -o 02_raw/trimmed/fastqc_trimmed/ -t 4
multiqc -o 02_raw/trimmed/fastqc_trimmed/ 02_raw/trimmed/fastqc_trimmed
# open and view 02_raw/trimmed/fastqc_trimmed/multiqc_report.html          
```

If no combining of different samples together is needed, simply copy links from the trimmed folder to the sample folder:    
`cp -l 02_raw/trimmed/*.fastq.gz 04_samples/`


### Align reads against reference genome ###
Download and index your reference genome.
e.g., `bwa index GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz`        

Copy the link to your genome and index files into `03_genome` folder.     
e.g., `cp -l GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna* 03_genome`      

Update the genome name in the following script and launch:    
`01_scripts/bwa_mem_align_reads_pe.sh <num_cores>`       

The output of this will be bam files, sorted bam files, and index bam files in `04_samples`.    
Note: can comment out line in script to delete unsorted bams if needed to preserve space.     


### Inspect alignments ###  
Inspect coverage in bins across chromosomes:      
```
# Prepare genome by unzipping, then unwrapping (E. Normandeau Scripts function)     
gunzip -c 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz > 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna  

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

### Mark duplicates ###
```
## Mark duplicates (optical and PCR): 
01_scripts/04_markduplicates.sh 
#  note: currently only runs on single core

## Assess mark duplicates output:       
# Obtain data rows from the duplicate metric individual files
grep 'Library' 04_samples/*metrics* > 04_samples/markdups_summary.txt

# Calculate the average proportion of duplicates (optical and other) for all files
awk '{ total += $9 } END { print total/NR }' 04_samples/markdups_summary.txt

```
note: this only marks duplicates, if want to remove, need to update (#TODO)


### Additional filtering ###
(not yet implemented)

### Merge reads ###
Make a list of all of the sorted bam files with duplicates marked:      
`ls -1 04_samples/*.sorted_mdups.bam > 04_samples/bamlist.txt`      

Merge the bam files:    
`samtools merge 05_genotyping/all_merged.bam -b 04_samples/bamlist.txt --threads 24`     


### Genotyping ###
Index the reference genome:     
`samtools faidx 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna`      
Note: this will not work on the compressed (.gz) genome.        


Freebayes approach:     
```
freebayes-parallel <(fasta_generate_regions.py 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai 100000) 26 -f 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna -L 05_genotyping/all_merged.bam --haplotype-length 0 -kwVa --throw-away-complex-obs --throw-away-complex-obs > 05_genotyping/genotypes.vcf
```

Samtools approach:        
```
bcftools mpileup -D -d 11500 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -f 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna 05_genotyping/all_merged.bam --threads 36 | bcftools call -mv --annotate GQ -Ob -o 05_genotyping/mpileup_calls.bcf --threads 36
```
note: the above was done with v.1.16.1, and the -D option may have been discontinued in newer versions.    

### Filtering ###     
e.g., 
`bcftools view -i '%QUAL>=20 && FORMAT/DP>10' 05_genotyping/mpileup_calls.bcf | grep -vE '^##' | less`    

Current filter:        
```
bcftools view -i 'F_missing < 0.1 & TYPE="snp" & %QUAL>=20 & FORMAT/DP>10' --min-alleles 2 --max-alleles 2 05_genotyping/mpileup_calls.bcf -Ob -o 05_genotyping/mpileup_calls_filt.bcf

# Settings:   
# F_missing:      fraction of missing genotypes per locus
# TYPE="snp":     keep only SNPs
# QUAL:           SNP quality value
# DP:             depth (per sample; #TODO: confirm)
# --min-alleles:  at least this many alleles observed per locus
# --max-alleles:  at most this many alleles observed per locus

# still to set: max depth? Or should this be done in the alignment stage? (#TODO)

```     

Add the AF to the dataset, then filter:       
```
# Add AF info
bcftools +fill-tags 05_genotyping/mpileup_calls_filt.bcf -Ob -o 05_genotyping/mpileup_calls_filt_AF.bcf -- -t AF
# note: bcftools plugin options are after the '--' indicator     

# Filter
bcftools view -i 'INFO/AF > 0.05' 05_genotyping/mpileup_calls_filt_AF.bcf -Ob -o 05_genotyping/mpileup_calls_filt_AF_0.05.bcf

```

Filter based on linkage:     
```
bcftools +prune -m 0.5 -w 50kb 05_genotyping/mpileup_calls_filt_AF_0.05.bcf -Ob -o 05_genotyping/mpileup_calls_filt_AF_0.05_LD.0.5.50kb.bcf    
```

Optional: make random subset of variants for data exploration:      
```
bcftools view 05_genotyping/mpileup_calls_filt_AF_0.05_LD.0.5.50kb.bcf | vcflib vcfrandomsample -r 0.01 > 05_genotyping/mpileup_calls_filt_AF_0.05_LD.0.5.50kb_subset_0.01.vcf

```

Bring the VCF or the subset VCF to your next analysis program.    


