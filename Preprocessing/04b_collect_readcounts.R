#!/usr/bin/env Rscript

## Structure of code, main and args by Michael M. Hoffman
## Copyright 2020 Michael M. Hoffman <michael.hoffman@utoronto.ca>

print("load libraries")
#Load libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(stringi)
library(doParallel)
library(strict)
library(tools)

main <- function(filename) {
    data <- read.csv(filename, header = FALSE, stringsAsFactors = TRUE)
    colnames(data) <- c("X", "UMI", "seq", "rname", "pos", "read_count")
    data$X <- NULL

print(head(data))

##clean data
### create chrmosome column
### remove uncharacterized genomic regions
data$chr <- substr(data$rname, 0, 5)
data$chr <- str_replace(data$chr, "_", "")
data <- as.data.frame(data[!grepl("*KI*", data$UMI), ])
data <- as.data.frame(data[!grepl("*GL*", data$UMI), ])
data <- as.data.frame(data[!grepl("*EBV*", data$UMI), ])
data <- as.data.frame(data[!grepl("*Un*", data$UMI), ])

###Calculating fragment length, GC content and CpG number
string_counter <- function(strings, pattern) {
  counts <- NULL
  for (i in 1:seq_len(strings)) {
    counts[i] <- length(attr(gregexpr(pattern, strings[i])[[1]],
                             "match.length")[attr(gregexpr(pattern,
                              strings[i])[[1]], "match.length") > 0])
  }
return(counts)
}

###To perform DNAm correction from the control data
####calclate fragment_len, GC content, and number of CpGs in sequence

print("calculating fragment length")
print("start variable calculations")
data$seq <- as.character(data$seq)
data$fragment_len <- vapply(as.character(data$seq), nchar, numeric(1))

print("calculating GC content")
data$GC <- (string_counter(data$seq, pattern = "C") +
              string_counter(data$seq, pattern = "G")) / data$fragment_len

print("calculating CpG number")
data$CpG <- string_counter(data$seq, pattern = "CG")
data[is.na(data)] <- 0

data <- data[, c("chr", "pos", "fragment_len", "UMI", "CpG", "GC", "read_count")]

print(head(data))

    write.csv(data, file = paste0("/cluster/projects/hoffmangroup/data_samanthawilson/",
                                  file_path_sans_ext(basename(filename)), "_readcount.csv"))
}

args <-
if (length(commandArgs(TRUE))
    || commandArgs()[length(commandArgs())] == "--args") {
    as.list(commandArgs(TRUE))
} else {
    list()
}

 do.call(main, args)
