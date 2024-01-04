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

#### Citation ####
If you find this tool useful, please cite the original article that uses the tool:        

Please also be sure to cite the tools applied within each function.      

## Sections ##

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

If no combining is needed, simply copy links from the trimmed folder to the sample folder:    
`cp -l 02_raw/trimmed/*.fastq.gz 04_samples/`


### Align against reference genome ###
Download and index your reference genome.
e.g., `bwa index GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz`        

Copy the link to your genome and index files into `03_genome` folder.     
e.g., `cp -l GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna* 03_genome`      

Update the genome name in the following script and launch:    
`01_scripts/bwa_mem_align_reads_pe.sh <num_cores>`       

The output of this will be bam files, sorted bam files, and index bam files in `04_samples`.    
Note: can comment out line in script to delete unsorted bams if needed to preserve space.     


### Inspecting alignments ###  
Optional:     
If you still have your original .bam files, you can inspect via:    
`samtools view -Sbt 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz 04_samples/COARL2-01-43986597_S84_L004_R1.bam | samtools flagstat -`      

Interactively run `01_scripts/inspect_alignments.R` to generate a chromosome-only GFF, which will be output as `03_genome/genome_chr.gff`.    

Use bedtools to calculate per-base coverage:    
`bedtools coverage -a 03_genome/genome_chr.gff -b 04_samples/COARL2-01-43986597_S84_L004_R1.sorted.bam -sorted -d > COARL2-01-43986597_S84_L004_R1_coverage.txt`


Alternate approach to quantify coverage:      
```
# Unwrap fasta
gunzip -c 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz > 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna  

fasta_unwrap.py 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap.fna

# Keep chr only
grep -E -A1 '^>NC_047' 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap.fna | grep -vE '^--$' - > 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap_chr_only.fna

# Determine lengths of chromosomes (for interest only)
fasta_lengths.py 03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap_chr_only.fna   

# Index the unwrapped, chr-only reference
samtools faidx ./03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap_chr_only.fna

# Determine window positions per chr
bedtools makewindows -g ./03_genome/GCF_902806645.1_cgigas_uk_roslin_v1_genomic_unwrap_chr_only.fna.fai -w 1000000 > 03_genome/windows.bed

# Calculate coverage per window
bedtools coverage -a 03_genome/windows.bed -b 04_samples/COARL2-01-43986597_S84_L004_R1.sorted.bam > 04_samples/COARL2-01_coverage.txt  

# Use Rscript to summarize and plot

```






