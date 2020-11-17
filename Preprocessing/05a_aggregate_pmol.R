#/usr/bin/env Rscript

## Structure of code, main and args by Michael M. Hoffman
## Copyright 2020 Michael M. Hoffman <michael.hoffman@utoronto.ca>

#Load libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(forcats)
library(purrr)
library(stringi)
library(strict)
library(tools)
library(zonator)

main <- function(filename) {

data <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)

##Remove header lines from the file
data <- data[!data$UMI == "UMI", ]

##Structure the data
data$fragment_len <- as.numeric(data$fragment_len)

####Cleaning data
data$pos <- as.numeric(data$pos)
data$start <- data$pos
data$end <- data$start + data$fragment_len
data$pos <- NULL
data$X <- NULL
data$fmol <- as.numeric(data$fmol)

##Removing uncharacterized genes
### This is already done, but as a double check
data$UMI <- as.factor(data$UMI)
data <- as.data.frame(data[!grepl("*KI*", data$UMI), ])
data <- as.data.frame(data[!grepl("*GL*", data$UMI), ])
data <- as.data.frame(data[!grepl("*EBV*", data$UMI), ])
data$UMI <- droplevels(data$UMI)
print(head(data))

##Replace NA with 0
data[is.na(data)] <- 0
print(head(data))

data$grp <- as.factor(paste0(data$chr, "_", data$start, "_", data$end))
##put in as a check for aggregation
print(head(data))

##Now save the .bed file to compare to hg38
data <- data[, c("chr", "start", "end", "pmol")]
colnames(data) <- c("chrom", "chromStart", "chromEnd", paste0(file_path_sans_ext(basename(filename))))
data$chrom <- as.character(data$chrom)
data$chromStart <- as.integer(data$chromStart)
data$chromEnd <- as.integer(data$chromEnd)
print(head(data))

write.table(data, quote = FALSE, row.names = FALSE, sep = "\t",
            file = paste0("/cluster/projects/hoffmangroup/data_samanthawilson/",
                          file_path_sans_ext(basename(filename)), ".bed"))
}

args <-
if (length(commandArgs(TRUE))
    || commandArgs()[length(commandArgs())] == "--args") {
    as.list(commandArgs(TRUE))
} else {
    list()
}

 do.call(main, args)
