#!/usr/bin/env Rscript

##Load packages
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(stringr)
library(varhandle)
library(forcats)
library(purrr)
library(stringi)
library(tidyverse)
library(tools)

main <- function(filename) {

##Load in data
data <- read.table(filename, sep = "\t" , stringsAsFactors = TRUE, fill = TRUE)
colnames(data) <- c("chr_frag", "start_frag", "end_frag", "fragment_len", "GC", "CpG", "CpG_3", "read_count", "pmol", "chrom", "chromStart", "chromEnd", "overlap_bp")

###Adjust the methylation values (overlaps/300*conc)
window_size <- 300

data$multipler <- data$overlap_bp/window_size

samples <- data[ , grepl( "pmol" , names(data))]
adj_vals <- samples*data$multipler

data_adj <- cbind(data[,c("chrom", "chromStart", "chromEnd", "read_count")],adj_vals)

#Aggregate by window
##sum pmol and read
data_adj$read_count <- as.numeric(data_adj$read_count)
data_adj_binned <- ddply(data_adj, .(data_adj$chrom, data_adj$chromStart, data_adj$chromEnd), numcolwise(sum))

print(head(data_adj_binned))

write.table(data_adj_binned, sep = "\t", row.names = FALSE, quote = FALSE, file = paste0("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/",
                                       file_path_sans_ext(basename(filename)), "_adjbinned.bed"))

}

args <-
if (length(commandArgs(TRUE))
    || commandArgs()[length(commandArgs())] == "--args") {
    as.list(commandArgs(TRUE))
} else {
    list()
}

 do.call(main, args)

