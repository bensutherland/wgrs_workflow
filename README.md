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
[00. Getting started](#00-getting-started)
[01. Prepare data](#01-prepare-data)
[02. Quality control and trimming](#03-quality-control-and-trimming)


## 00. Getting started ##
Clone this repository and change into the main directory.      
```
git clone https://github.com/bensutherland/wgrs_workflow.git
cd wgrs_workflow   

```

### 01. Prepare data ###
Copy links to all raw data in `02_raw`.      


### 02. Quality control and trimming ###
View raw data with fastqc and multiqc:     
```
fastqc 02_raw/*.fastq.gz -o 02_raw/fastqc_raw/ -t 5
multiqc -o 02_raw/fastqc_raw/ 02_raw/fastqc_raw
```

Trim adapters and for quality in parallel:     
`./01_scripts/01_cutadapt_PE.sh <num_cores>`     
...this will trim for quality (-q 15), minimum length (-m 50), and terminal Ns. It will also remove adapters that are present in `00_archive/adapters.fasta`.    

note: the trimming script  assumes that the read file suffix for R1 is `R1_001.fastq.gz`.    

View trimmed data with fastqc and multiqc:      
```
fastqc 02_raw/trimmed/*.fastq.gz -o 02_raw/trimmed/fastqc_trimmed/ -t 4
multiqc -o 02_raw/trimmed/fastqc_trimmed/ 02_raw/trimmed/fastqc_trimmed
# open and view 02_raw/trimmed/fastqc_trimmed/multiqc_report.html          
```

If no combining of different samples together is needed, simply copy links from the trimmed folder to the sample folder:    
`cp -l 02_raw/trimmed/*.fastq.gz 04_samples/`


### 03. Align reads against reference genome ###
Download and index your reference genome.
e.g., `bwa index GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz`        

Copy the link to your genome and index files into `03_genome` folder.     
e.g., `cp -l GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna* 03_genome`      

Update the genome name in the following script and launch:    
`01_scripts/bwa_mem_align_reads_pe.sh <num_cores>`       

The output of this will be sorted and indexed bam files in `04_samples`.    
Note: can assess results (see below) before or after marking duplicates.     


### 04. Mark duplicates ###
```
## Mark duplicates (optical and PCR): 
01_scripts/04_markduplicates.sh 
#  note: currently only runs on single core

## Assess mark duplicates output:       
# Obtain data rows from the duplicate metric individual files
grep 'Library' 04_samples/*metrics* > 04_samples/markdups_summary.txt

# Note: see any individual metrics file to obtain column names for the above markdups summary
# Open the output in excel or similar to view

```


### 05. Inspect alignment and duplicate removal results ###  
Compare the number of trimmed reads per sample to the number of reads aligning after filters (-q 10):      
`01_scripts/assess_results`     
...output will provide a table of per-sample number of reads and number of reads aligning with MAPQ > 10.    
Note: to compare directly, it will be necessary to multiply the number of reads by two, given that both reads (R1 and R2 are in the alignment file), and assume that the number of reads in R1 and R2 per sample are equal.     

Optional: also see `20_docs/inspect_chr.md` for a more specialized approach in viewing alignments across chromosomes.    


### 06. Merge reads ###
Make a list of all of the sorted bam files with duplicates marked:      
`ls -1 04_samples/*.sorted_mdups.bam > 04_samples/bamlist.txt`      

Merge the bam files:    
`samtools merge 05_genotyping/all_merged.bam -b 04_samples/bamlist.txt --threads 24`     


### 07. Genotyping ###
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
Use the following script to filter the bcf       
`01_scripts/filter_bcf.sh`      

Run with defaults, this script will: 
- remove variants within 5 bp of an indel
- remove variants with more than 10% missing data across individuals
- keep only variants of type 'SNP'
- keep SNPs with a quality of minimum 20 (in any one individual)
- keep SNPs with an average depth across all samples of at least 10 reads
- keep only biallelic SNPs
- remove genotypes (per individual, per site) that have fewer than 10 reads or more than 200 reads
- remove variants with more than 10% missing data across individuals (again)

At each stage, the script will report the number of variants retained.    

Current filter:        
```
# first filter out any variant that is within 5 bp of an indel:    
bcftools filter --threads 20 -g 5 05_genotyping/mpileup_calls.bcf -Ob -o 05_genotyping/mpileup_calls_noindel5.bcf

# second, filter out based on other parameters
bcftools view -i 'F_missing < 0.1 & TYPE="snp" & %QUAL>=20 & FORMAT/DP>10 & AVG(FMT/DP)>10' --min-alleles 2 --max-alleles 2 05_genotyping/mpileup_calls_noindel5.bcf -Ob -o 05_genotyping/mpileup_calls_noindel5_miss0.1_bialSNP_q20_allDP10_avgDP10.bcf --threads 42

# Settings:   
# F_missing:      fraction of missing genotypes per locus
# TYPE="snp":     keep only SNPs
# QUAL:           SNP quality value
# DP:             depth (across all samples)
# --min-alleles:  at least this many alleles observed per locus
# --max-alleles:  at most this many alleles observed per locus

# then also use vcftools to remove low coverage sites:    
vcftools --bcf 05_genotyping/mpileup_calls_noindel5_miss0.1_bialSNP_q20_allDP10_avgDP10.bcf --minDP 4 --maxDP 100 --recode-bcf --recode-INFO-all --out 05_genotyping/mpileup_calls_noindel5_miss0.1_bialSNP_q20_allDP10_avgDP10_minDP4_maxDP100

# Alternately, use bcftools
bcftools filter -S . -e "(FORMAT/AD[*:0] + FORMAT/AD[*:1]) < 4 | (FORMAT/AD[*:0] + FORMAT/AD[*:1]) > 100" 05_genotyping/mpileup_calls_noindel5_miss0.1_bialSNP_q20_allDP10_avgDP10.bcf -Ob -o 05_genotyping/mpileup_calls_noindel5_miss0.1_bialSNP_q20_allDP10_avgDP10_minDP4_maxDP100.bcf --threads 28

# Re-index
bcftools index 05_genotyping/mpileup_calls_noindel5_miss0.1_bialSNP_q20_allDP10_avgDP10_minDP4_maxDP100.bcf 

# Filter again for high missing (after removing low coverage genotypes)
bcftools view -i 'F_missing < 0.1' 05_genotyping/mpileup_calls_noindel5_miss0.1_bialSNP_q20_allDP10_avgDP10_minDP4_maxDP100_bcftools.bcf -Ob -o 05_genotyping/mpileup_calls_noindel5_miss0.1_bialSNP_q20_allDP10_avgDP10_minDP4_maxDP100_bcftools_miss0.1.bcf --threads 42

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


