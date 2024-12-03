# wgrs_workflow
Workflow to support genotyping of moderate/ high depth whole-genome resequence data.       
Developed by Ben J. G. Sutherland, Ph.D. (Sutherland Bioinformatics & VIU) in the working group of Timothy J. Green (VIU).     

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
- bcftools

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

Once fastqc has been completed, move the samples to the sample folder:      
`mv 02_raw/trimmed/*.fastq.gz 04_samples/`     
...note: mv is used instead of cp -l to avoid duplicate storage if the directory is backed up elsewhere.  



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

Genotype with bcftools:        
```
bcftools mpileup -D -d 12000 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -f 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna 05_genotyping/all_merged.bam --threads 36 | bcftools call -mv --annotate GQ -Ob -o 05_genotyping/mpileup_calls.bcf --threads 36
```
note: the above was done with v.1.16.1, and the -D option may have been discontinued in newer versions.    
note: the -d 12000 was selected based on a project with 48 samples (48 x 250x max depth each)    

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
There will also be a log file produced as `<date>_<time>_filter_bcf.log`.    

Note: all intermediate filter files will be retained, so if needed, delete to save space.    

Alternately, this can be run manually to output only a single file:            
```
bcftools filter --threads 42 -g 5 05_genotyping/mpileup_calls.bcf | bcftools view -i 'F_missing < 0.1 & TYPE="snp" & %QUAL>=20 & AVG(FMT/DP)>10' --min-alleles 2 --max-alleles 2 --threads 42 | bcftools filter -S . -e "(FORMAT/AD[*:0] + FORMAT/AD[*:1]) < 4 | (FORMAT/AD[*:0] + FORMAT/AD[*:1]) > 100" --threads 42 | bcftools view -i 'F_missing < 0.1' -Ob -o 05_genotyping/mpileup_calls_filtered.bcf --threads 42     

# View number of variants retained
bcftools view 05_genotyping/mpileup_calls_filtered.bcf | grep -vE '^#' - | wc -l

# This approach can be taken, but it isn't as good as the stepwise approach because the stepwise can allow an individual check of how many variants are retained at each stage. However, it will be much faster and produce the same output as the final output from the above script (with less space taken up).     
```     

**Additional filtering (specific to different downstream methods)**   
Allele Frequency:    
```
## Allele frequency 
# note, below replace <filtered> with your filtered filename

# Note, if want to add the AF to the dataset
bcftools +fill-tags 05_genotyping/<filtered>.bcf -Ob -o 05_genotyping/<filtered>_AF.bcf -- -t AF
# note: bcftools plugin options are after the '--' indicator     

# Filter for minor allele frequency
bcftools view -q 0.05:minor 05_genotyping/<filtered>.bcf -Ob -o 05_genotyping/<filtered>_maf0.05.bcf
# note: it is important to set the type of allele being filtered by the above (here: minor)

# How many variants remain? 
bcftools view 05_genotyping/<filtered>_maf0.05.bcf | grep -vE '^#' - | wc -l    
```

Linkage:     
```
## Linkage
# note, below replace <filtered> with your filtered filename
bcftools +prune -m 0.5 -w 50kb 05_genotyping/<filtered>.bcf -Ob -o 05_genotyping/<filtered>_LD0.5w50kb.bcf    

# How many variants remain?   
bcftools view 05_genotyping/<filtered>_LD0.5w50kb.bcf | grep -vE '^#' - | wc -l

```

**Optional subsetting to simplify troubleshooting**    
```
bcftools view 05_genotyping/<LD-filtered>.bcf | vcflib vcfrandomsample -r 0.05 > 05_genotyping/<LD-filtered.bcf>_subset0.05.vcf


```

Bring the BCF or the subset VCF to your next analysis program.    


