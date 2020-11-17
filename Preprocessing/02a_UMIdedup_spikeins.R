#!/usr/bin/env Rscrupt

#This script will:
#read in .bam files
#convert into a data frame
#parse the last 11 characters of the read ID to dedup and collapse the UMIs

## Load needed libraries
library(forcats)
library(Rsamtools)
library(stringr)
library(plyr)
library(doParallel)
library(caroline)

PATH_PREFIX <-  "~/Projects/2018_PTB/data/2020_BatchAnalysis/Batch1_RS/"
bams_list <- list.files(path = PATH_PREFIX, pattern = "filtered_aligned")

##Characterize parallelization parameters
cl <- makeCluster(length(bams_list))
registerDoParallel(cl)

##Parallelize whole script
foreach(i = 1:seq_len(bams_list),
        .packages = c("caroline", "Rsamtools", "forcats",
                      "stringr", "plyr", "doParallel")) %dopar% {

##Specify path to all the Bamfiles I want processed
bam_file <- paste0(PATH_PREFIX, bams_list)

##Scan each bam
bam <- canBam(bam_file[[i]])

##Convert bam into a data frame
### This code was obtained from Dave Tang's blog
####[https://gist.github.com/davetang/6460320]

###Function for collapsing the list of lists into a single list
###as per the Rsamtools vignette
.unlist <- function(x) {
   ## do.call(c, ...) coerces factor to integer, which is undesired
   x1 <- x[[1L]]
   if (is.factor(x1)) {
      structure(unlist(x), class = "factor", levels = levels(x1))
   } else {
      do.call(c, x)
   }
}


###store names of BAM fields
bam_field <- names(bam[[1]])

###go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))

###store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)

### End of Dave Tang's code

#In this dataframe, each row in a read
##Assign a new column called read_count and assign each row a 1
bam_df$read_count <-  1

##Assign new column with UMI
bam_df$UMI <- str_sub(bam_df$qname, -11, -1)
bam_df$UMI <-  as.factor(paste0(bam_df$rname, bam_df$UMI))

bam_df$fragment_len <- as.factor(str_sub(bam_df$UMI, 0, 3))
bam_df$fragment_len <- revalue(bam_df$fragment_len, c("80b" = "80"))

bam_df$CpG <- as.factor(str_sub(bam_df$UMI, 5, 7))
bam_df$CpG <- revalue(bam_df$CpG, c("_16" = "16", "_2C" = "2", "_4C" = "4",
                                    "_8C" = "8", "1C_" = "1",
                                    "2C_" = "2", "4C_" = "4"))

bam_df$GC <- as.factor(str_sub(bam_df$UMI, 7, 11))
bam_df$GC <- revalue(bam_df$GC, c("_35G-" = "35", "_50G-" = "50", "_65G-" = "65",
                                  "6C_35" = "35", "6C_50" = "50",
                                  "6C_65" = "65", "C_35G" = "35",
                                  "C_50G" = "50", "C_65G" = "65"))

bam_df$methylation_status <- str_extract(bam_df$UMI, "meth")
bam_df$methylation_status <-  bam_df$methylation_status %>% fct_explicit_na(na_level = "unmeth")

bam_df$frag_grp <- as.factor(paste0(bam_df$fragment_len, sep = "_",
                                    bam_df$CpG, sep = "_", bam_df$GC))

###Collapsing the PCR duplicates
##Collapse duplicates. If UMI and rname match between samples,
##collapse the reads and add the read_count columns together
read_counts <- aggregate(read_count~UMI, data = bam_df, FUN = sum)

##I will make all read_counts of every row=1 as each is a unique fragment
read_counts$read_count <- 1

##Populating ReadCounts dataframe with other information I would want to plot
read_counts$fragment_len <- as.factor(str_sub(read_counts$UMI, 0, 3))
###Clean up the factors by renaming
read_counts$fragment_len <- revalue(read_counts$fragment_len, c("80b" = "80"))

read_counts$CpG <- as.factor(str_sub(read_counts$UMI, 5, 7))
###Clean up the factors by renaming
read_counts$CpG <- revalue(read_counts$CpG, c("_16" = "16", "_2C" = "2",
                                              "_4C" = "4", "_8C" = "8", "1C_" = "1",
                                              "2C_" = "2", "4C_" = "4"))

read_counts$GC <- as.factor(str_sub(read_counts$UMI, 7, 11))
###Clean up the factors by renaming
read_counts$GC <- revalue(read_counts$GC, c("_35G-" = "35", "_50G-" = "50",
                                          "_65G-" = "65", "6C_35" = "35",
                                          "6C_50" = "50", "6C_65" = "65",
                                          "C_35G" = "35", "C_50G" = "50",
                                          "C_65G" = "65"))

read_counts$methylation_status <- str_extract(read_counts$UMI, "meth")
read_counts$methylation_status <- read_counts$methylation_status %>% fct_explicit_na(na_level = "unmeth")

read_counts$frag_grp <- read_counts$frag_grp <- as.factor(paste0(read_counts$fragment_len, sep = "_",
                                                               read_counts$CpG, sep = "_", read_counts$GC))

dim(read_counts) ##Output the number of unique reads
read_counts <- as.data.frame(read_counts)

##Save file of all reads (with PCR duplicates) for later plotting out all read distribution
name <- paste0(bams_list[[i]], "_allreads_nodedup.csv")
write.delim(bam_df, file = paste0(bams_list[[i]], "_allreads_nodedup.csv"), sep = "\t")
head(bam_df)

##Save deduplicated data, name per sample
name <- paste0(bams_list[[i]], "_Cleaned_dedup_readcounts.csv")
write.delim(read_counts, name, sep = ",")
head(read_counts)
}

stopCluster(cl)
