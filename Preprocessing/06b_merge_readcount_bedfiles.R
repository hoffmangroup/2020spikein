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

##Multimerge samples into single data frame
PATH_PREFIX <- "/cluster/projects/hoffmangroup/data_samanthawilson/"
samples <- list.files(path <- PATH_PREFIX,
                      pattern = "*rc_all.bed_intersect_hg38_aggbywindow.bed")
pattern <- "*rc_all.bed_intersect_hg38_aggbywindow.bed"

##Load files and merge into one file
###Multimerge function adapted from R-bloggers:
#### https://www.r-bloggers.com/merging-multiple-data-files-into-one-data-frame/
multmerge <- function(PATH_PREFIX, pattern) {
  filenames <- list.files(path = PATH_PREFIX, pattern = pattern, full.names = TRUE)
  ##List filenames and truncate to only sample number
  labels <- stri_extract_last(filenames, regex = "\\d{4}")
  ##Read in files
  datalist <- lapply(filenames, function(x) {
      read.table(file = x, header = TRUE)})
  #label colnames
  for (i in seq(datalist)) {
        for (i in seq(labels)) {
  names(datalist[[i]])[names(datalist[[i]]) == "read_count"] <- paste(labels[[i]])
}}

Reduce(function(...) full_join(..., by = c("window")), datalist)
}

print("starting multmerge")
data <- multmerge("/cluster/projects/hoffmangroup/data_samanthawilson/", "*rc_all.bed_intersect_hg38_aggbywindow.bed")
print("merge completed")

##remove columns is .x or .y in them
data <- data[, -grep(".x", colnames(data))]
data <- data[, -grep(".y", colnames(data))]

print(head(data))

##saving this for future analysis on read counts
write.csv(data, file = "/cluster/projects/hoffmangroup/data_samanthawilson/data_rc.bed")
