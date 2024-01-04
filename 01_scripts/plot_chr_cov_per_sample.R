# Plot chr coverage per sample, requires bedtools has been already used to calc cov (see README)
# Ben J. G. Sutherland, 2024-01-04

# Clear space
#rm(list=ls())

# Load packages
#install.packages("ggplot2")
library("ggplot2")

# Set working directory
current.path <- dirname(rstudioapi::getSourceEditorContext()$path)
current.path <- gsub(pattern = "\\/01_scripts", replacement = "", x = current.path)
setwd(current.path)
rm(current.path)

## Info
# sessionInfo()

# Set variables
options(scipen = 9999999)

#### 00. Read in data ####
# Identify files
cov_files.vec <- list.files(path = "04_samples/", pattern = "_cov.txt", full.names = T)

# Read in cov files
cov_files.list <- list(); coverage.df <- NULL

for(i in 1:length(cov_files.vec)){
  
  # Identify sample name
  sample.id <- gsub(pattern = "04_samples//", replacement = "", cov_files.vec[i])
  sample.id <- gsub(pattern = "_cov.txt", replacement = "", sample.id)
  
  # Read in and format the sample coverage file
  coverage.df <- read.delim(file = cov_files.vec[i], header = F)
  colnames(coverage.df) <- c("LG", "start", "stop", "num.aligns", "num.non-zero.bp", "window.length", "fraction.non.zero")
  coverage.df$start <- as.numeric(coverage.df$start)
  coverage.df$stop  <- as.numeric(coverage.df$stop)
  
  # Sort by chr
  coverage.df <- coverage.df[with(coverage.df, order(coverage.df$LG, coverage.df$start)), ]
  
  # Save to list
  cov_files.list[[sample.id]] <- coverage.df
  
}


#### 01. Determine chr lengths ####
len.chr = NULL; len.chr.list <- list()
for(i in 1:length(unique(coverage.df$LG))){

  # What are the chromosomes? 
  chr.o.i <- unique(coverage.df$LG)[i]

  # Save the length of each chr
  len.chr.list[[chr.o.i]] <- max(coverage.df[coverage.df$LG==chr.o.i, "stop"])

}

# Retain as dataframe
len.chr.df <- as.data.frame(unlist(len.chr.list))
len.chr.df$chr <- rownames(len.chr.df)
colnames(len.chr.df) <- c("len", "chr")
len.chr.df <- len.chr.df[,c("chr", "len")]


#### 02. Calculate coverage per sample per chr ####
coverage.df <- NULL ; sample.id <- NULL; per_chr_cov.list <- list()
for(i in 1:length(cov_files.list)){
  
  # Select the sample ID
  sample.id <- names(cov_files.list)[i]
  
  # Extract the coverage dataframe
  coverage.df <- cov_files.list[[sample.id]]
  
  # Calculate the number of alignments per chromosome
  aligns_by_chr.df <- aggregate(coverage.df$num.aligns, by=list(LG=coverage.df$LG), FUN=sum)
  colnames(aligns_by_chr.df) <- c("LG", "num.aligns")
  
  per_chr_cov.list[[sample.id]] <- aligns_by_chr.df
  
}

# Note: I don't think we are using this anywhere yet? 

# # To test out the method above
# sum(coverage.df[coverage.df$LG=="NC_047559.1", "num.aligns"]) # to confirm matches above

# # What is the total coverage for the sample? 
# sum(coverage.df$num.aligns) # 47,304,597 (that's about right, remember, this is only chrs, not scaffs)

# Clean up
rm(coverage.df)
rm(aligns_by_chr.df)


# Summary
head(len.chr.df)
names(per_chr_cov.list)
names(cov_files.list)
cov_files.list[[1]][1:2,] # e.g. 


#### SKIP #####
#### 03. Combine variant counts with chromosome lengths
for(i in 1:length(cov_files.list)){
  
  # Identify sample ID
  sample.id <- names(cov_files.list)[i]
  
}


# Combine variant counts with chr
head(aligns_by_chr.df)
head(len.chr.df)
chr_length_and_cumul_align.df <- merge(x = len.chr.df, y = aligns_by_chr.df, by.x = "chr", by.y = "LG")
head(chr_length_and_cumul_align.df)
colnames(chr_length_and_cumul_align.df) <- c("chr", "length", "cumul.align")
head(chr_length_and_cumul_align.df)

# Convert bp to kbp
chr_length_and_cumul_align.df$length_per.kb <- chr_length_and_cumul_align.df$length / 1000
head(chr_length_and_cumul_align.df) # This will be used below for plotting


###### END SKIP ####


#### 03. Prepare windows (binned aligns) across chr ####
sample.id <- NULL; coverage.df <- NULL; bins <- NULL; bins.df <- NULL
per_chr_cov_with_cols.list <- list()
for(i in 1:length(cov_files.list)){
  
  # Identify sample ID
  sample.id <- names(cov_files.list)[i]
  
  # Pull coverage.df file
  coverage.df <- cov_files.list[[sample.id]]
  
  # How many bins are there per chr and in total? 
  bins <- table(coverage.df$LG)
  bins.df <- as.data.frame(names(bins))
  colnames(bins.df) <- "chr"
  bins.df$bins <- as.numeric(bins)
  bins.df$cumsum <- cumsum(bins.df$bins) # cumulate total number of bins
  
  # Per chromosome, determine the bin numbers for which the chr starts, ends
  #    , and midway for the chr label
  bins.df$start.bin <- NA; bins.df$end.bin <- NA; bins.df$label.pos <- NA
  for(j in 1:nrow(bins.df)){
    
    if(j == 1){
      
      # The starting chr starts at zero
      bins.df$start.bin[j] <- 0
      
    }else{
      
      # Any subsequent chr starts at the end of the previous cumulative bin count
      bins.df$start.bin[j] <- bins.df$cumsum[j-1]
      
    }
    
    # Per chr, the end bind is at the end of the cumulative sum per chr
    bins.df$end.bin[j] <- bins.df$cumsum[j]
    
    # Labels go in the middle of each start and stop per chr
    bins.df$label.pos[j] <- ((bins.df$end.bin[j] - bins.df$start.bin[j])/2) + bins.df$start.bin[j]
    
  }
  
  head(bins.df)
  
  # Create alternating colour vector for plotting
  bins.df$plot.colour <- "black"
  bins.df$plot.colour[c(FALSE, TRUE)] <- "darkgrey"
  head(bins.df)
  
  # Inspect df to be used
  head(coverage.df)
  head(bins.df)
  
  # Bring colour from the bins df into the coverage df
  coverage.df$plot.color <- NA
  for(c in 1:length(coverage.df$plot.color)){
    
    coverage.df$plot.color[c] <- bins.df[bins.df$chr %in% coverage.df$LG[c], "plot.colour"]
    
  }
  
  #head(coverage.df, n = 60)
  
  # Add plotting order
  coverage.df <- coverage.df[with(coverage.df, order(coverage.df$LG, coverage.df$start)), ]
  coverage.df$plot.order <- seq(1:nrow(coverage.df))
  
  # Save output
  per_chr_cov_with_cols.list[[sample.id]] <- coverage.df
  
  # # Save output to file
  # write.table(x = coverage.df, file = paste0("04_samples/", sample.id, "cov_with_header.txt"), quote = F
  #             , sep = "\t", row.names = F, col.names = T
  # ) 
  # 
  
  
}

names(per_chr_cov_with_cols.list)
head(per_chr_cov_with_cols.list[[1]], n  = 10)



#### 04. Plotting ####
coverage.df <- NULL; sample.id <- NULL
for(i in 1:length(per_chr_cov_with_cols.list)){
  
  # Pull data from list
  sample.id <- names(per_chr_cov_with_cols.list)[i]
  coverage.df <- per_chr_cov_with_cols.list[[sample.id]]
  
  # Plot alignment per bin across the chr
  pdf(file = paste0("04_samples/", sample.id, "_plot_aligns_per_window.pdf"), width = 9, height = 4)
  par(mfrow=c(1,1), mar = c(7,4.1,3, 2.1))
  plot(x = coverage.df$plot.order, y = coverage.df$num.aligns/1000
       , xaxt = "n"
       , xlab = ""
       , ylab = "Alignments (x1000) per window"
       , las = 1
       , pch = 16
       , col = coverage.df$plot.color
       , cex = 0.8
       , main = sample.id
  )
  abline(v = bins.df$end.bin, lty = 2)
  abline(v = bins.df$start.bin, lty = 2)
  axis(side = 1, at = bins.df$label.pos, labels = gsub(pattern = "LG", replacement = "", x = bins.df$chr)
       , tick = T
       , las = 2
  )
  
  dev.off()
  
}


# Finished, see result plots in "04_samples"
