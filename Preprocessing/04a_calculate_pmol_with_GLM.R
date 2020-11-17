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

print("start variable calculations")
data$seq <- as.character(data$seq)
data$fragment_len <- vapply(as.character(data$seq), nchar, numeric(1))

print("calculating GC content")
data$GC <- (string_counter(data$seq, pattern = "C") +
              string_counter(data$seq, pattern = "G")) / data$fragment_len

print("calculating CpG number")
data$CpG <- string_counter(data$seq, pattern = "CG")
data[is.na(data)] <- 0

data_melt <- melt(data,
                  id.vars = c("GC", "UMI", "fragment_len", "CpG", "pos", "chr"),
                  measure.vars = "read_count")
data_melt$value <- as.numeric(as.character(data_melt$value))

data_melt$CpG_3 <- (data_melt$CpG) ^ (1 / 3)

data_melt$GC <- as.numeric(as.character(data_melt$GC))
data_melt$fragment_len <- as.numeric(data_melt$fragment_len)

##Source the GLM from the control GLM code
load("~/Projects/2020_Control_Project/2020_BatchAnalysis/Batch3_DT/2020_Gaussian_batch3.rda")

print("start fmol prediction")
##Predict fmol values
pred_conc <- as.data.frame(predict(gaussian_glm, data_melt))
names(pred_conc) <- "value"
##Replace value with conc
data_melt$value <- NULL
data_melt <- cbind(data_melt, pred_conc)
data_melt$variable <- NULL
colnames(data_melt) <- c("GC", "UMI", "fragment_len", "CpG", "pos", "chr", "CpG3", "fmol")
data_melt <- data_melt[, c("chr", "pos", "fragment_len", "UMI", "CpG", "CpG3", "GC", "fmol")]

print(head(data_melt))

    write.csv(data_melt, file = paste0("/cluster/projects/hoffmangroup/data_samanthawilson/",
                                       file_path_sans_ext(basename(filename)), "_fmol.csv"))
}

args <-
if (length(commandArgs(TRUE))
    || commandArgs()[length(commandArgs())] == "--args") {
    as.list(commandArgs(TRUE))
} else {
    list()
}

 do.call(main, args)
