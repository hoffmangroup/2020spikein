#!/usr/bin/env Rscrupt

#This script will:
##read in .bam files
##convert into a data frame
##parse the last 11 characters of the read ID to dedup and collapse the UMIs

## Load needed libraries
library(forcats)
library(Rsamtools)
library(stringr)
library(plyr)
library(doParallel)
library(caroline)

PATH_PREFIX <-  "~/Projects/2018_PTB/data/2020_BatchAnalysis/Batch1_RS/"
bams_list <- list.files(path = PATH_PREFIX,
                        pattern = "filtered_human.aligned.*.bam")

##Define parallelization parameter

cl <- makeCluster(length(BamsList))
registerDoParallel(cl)

foreach(i = 1:seq_len(bams_list),
        .packages = c("caroline", "Rsamtools", "forcats",
                                                "stringr", "plyr", "doParallel")) %dopar% {

##Specify path to all the Bamfiles I want processed
bam_file <- paste0(PATH_PREFIX, bams_list)

##Scan each bam
print("start scan")
bam <- scanBam(bam_file[[i]])
print("done scan")

##Convert bam into a data frame
### This code was obtained from Dave Tang's blog
####[https://gist.github.com/davetang/6460320]

###function for collapsing the list of lists into a single list
.unlist <- function(x) {
   ## do.call(c, ...) coerces factor to integer, which is undesired
   x1 <- x[[1L]]
   if (is.factor(x1)) {
      structure(unlist(x), class = "factor", levels = levels(x1))
   } else {
      do.call(c, x)
   }
}

print("done unlist")
###store names of BAM fields
bam_field <- names(bam[[1]])

###go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
print("done list")
###store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field

dim(bam_df)

### End of Dave Tang's code

#In this dataframe, each row in a read
##Assign a new column called read_count and assign each row a 1
bam_df$read_count <- 1

##Assign new column with UMI
bam_df$UMI <- str_sub(bam_df$qname, -11, -1)
bam_df$UMI <-  as.factor(paste0(bam_df$rname, bam_df$UMI))

###Collapsing the PCR duplicates
read_counts <- aggregate(read_count ~ UMI + seq + rname + pos, data = bam_df, FUN = sum)

#I will make all read_counts of every row=1 as each is a unique fragment
read_counts$read_count <- 1

##Save deduplicated data, name per sample
name <- paste0(bams_list[[i]], "_Cleaned_dedup_readcounts.csv")
print(name)
print("saving dedup csv")
write.delim(read_counts, file = paste0(bams_list[[i]], "_Cleaned_dedup_readcounts.csv"), sep = ",")
head(read_counts)
print("dedup csv saved")

}

stopCluster(cl)