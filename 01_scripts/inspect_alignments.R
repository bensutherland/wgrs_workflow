# Checking alignments
# Following tutorial from marineomics.github.io
# B. Sutherland 2024-01-03

# Set variables
genome_gff.FN <- "~/genomes/GCF_902806645.1_genomic.gff"
chr_str <- "NC_047" # the string that is common to all chrs

## Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path) # take main directory
setwd(current.path)
rm(current.path)

#### Identify regions to calc coverage ####
# Subset the GFF for the chr or genomic regions for which you want to eval cov
# Read in GFF
data.gff <- read.table(file = genome_gff.FN, sep = "\t", quote = "")
colnames(data.gff) <- c("seqid", "source", "type", "start", "end", "score"
                        , "strand", "phase", "attributes"
                        )
head(data.gff,n = 2)
dim(data.gff) # e.g., 1,620,046, 9

# Limit to only regions, then only chr
scaff.gff <- data.gff[data.gff$type=="region", ]
dim(scaff.gff) # 236 scaffs

chr.gff <- scaff.gff[grep(pattern = chr_str, x = scaff.gff$seqid), ]
dim(chr.gff)
chr.gff

# Write out chr definitions
write.table(x = chr.gff, file = "03_genome/genome_chr.gff", row.names = F
            , sep = "\t", col.names = F, quote = F
            )

# Next, calculate coverage with bedtools (see README)


# Come back to this script to visualize across the genome
# Adapted from marineomics (Sara M. Schaal)

### PROCESS COVERAGE DATA TO VISUALIZE ACROSS THE GENOME ###
## Sara M. Schaal (adapted by B. Sutherland, 2024-01-03)

# Load libraries
library(ggplot2)
library(data.table)

# Set user variables
samps <- c("COARL2-01-43986597_S84_L004_R1") # Input samples of interest
colors <- c("steelblue2")
n <- 10000 # INPUT SIZE OF WINDOWS TO CALCULATE AVERAGE COVERAGE

## PROCESS DATA 
# first step through each sample
#for(i in 1:length(samps)){
  i <- 1

df <- fread(paste0(samps[i], "_coverage.txt"), sep = "\t", quote = "", data.table = FALSE)
head(df)

colnames(df) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "pos", "unkn")


colnames(data.gff) <- c("seqid", "source", "type", "start", "end", "score"
                        , "strand", "phase", "attributes"
)



chroms <- unique(df[,1])

df.sample.data <- NULL
  # step through each chromosome and break into the increment chunks you set with n
  for(j in 1:length(chroms)){
    df.chrom <- df[df$chrom == chroms[j], ]
    extra <- nrow(df.chrom) %% n
    # this next part is taking our windows set by n and giving each window a number id
    # this grouped dataframe is made by binding the original df.chrom with grouping values for the window size you want 
    # then finding what the last grouping value would be using because it will most like not be an even increment of your n
    # for example if you divide a chromosome by your n and get 3000.3 then you can easily input the first 3000
    # grouping values in with rep (middle part of the following cbind function) then the last grouping value will be 3001 which
    # you get by rounding using ceiling in the last part of this cbind
    grouped <- cbind(df.chrom, c(rep(1:(nrow(df.chrom)/n), each = n), rep(ceiling(nrow(df.chrom)/n), extra)))
    colnames(grouped)[4] <- "grouping"
    # now bind this chromosomes data in the full dataframe
    df.sample.data <- rbind(df.sample.data, grouped)
  }
  # finally take your new dataframe and calculate the average coverage for your increments using
  # this new grouping variable and the chromosome
  df.covAve <- aggregate(coverage~grouping + chrom, data = df.sample.data, FUN = mean)
  
  
  #### PLOTTING ####
  pdf(paste0("figures/", samps[i], "MaxcoveragePlot", n, ".pdf"), height= 15, width=15)
  
  ## I make two plots because there will undoubtedly be some loci that have really high coverage which makes
  # it hard to see all the spread of the majority of the data. This first graph is set to the max coverage
  # found in the data frame and the second plot is setting your y limit to a more reasonable value for your
  # data. For me 100 x was good but feel free to change to what is appropriate for your data.
  
  print(ggplot(data = df.covAve, aes(x = grouping, y = coverage)) +
          geom_point(col = colors[i], alpha = 0.5) +
          facet_wrap(~chrom) +
          labs(y = "Coverage", x = paste0("location every ", n, " bases"), 
               title = paste0(samps[i], "Genome Coverage up to Max Coverage")) +
          ylim(0, max(df.covAve$coverage)) +
          xlim(0, max(df.covAve$grouping)) +
          theme_bw() + 
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"), 
                legend.position = "none"))
  dev.off()
  
  pdf(paste0("figures/", samps[i], "100XcoveragePlot", n, ".pdf"), height= 15, width=15)
  
  print(ggplot(data = df.covAve, aes(x = grouping, y = coverage)) +
          geom_point(col = colors[i], alpha = 0.5) +
          facet_wrap(~chrom) +
          labs(y = "Coverage", x = paste0("location every ", n, " bases"),
               title = paste0(samps[i], "Genome Coverage up to 100X")) +
          ylim(0, 100) +
          xlim(0, max(df.covAve$grouping)) +
          theme_bw() + 
          theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"), 
                legend.position = "none"))
  dev.off()
  
}







