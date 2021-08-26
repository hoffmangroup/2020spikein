#!/usr/bin/env Rscript

#Load libraries
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
library(tools)

main <- function(filename) {
	data <- read.table(filename, sep = "\t", stringsAsFactors = TRUE, fill = TRUE)
	data <- data[,c(1:3, 6,13:15)]
	colnames(data) <- c("chr", "start", "end", "GC", "fragment_len", "CpG1", "CpG2")

	 ##CpGs have been read in as levels - make sure they are numeric
        data$CpG1 <- as.numeric(as.character(data$CpG1))
        data$CpG2 <- as.numeric(as.character(data$CpG2))

        #Check if CpG1 and CpG2 columns match - if they do, case sentivity in bedtools was not a problem and we only need one column
        #Otherwise, we need to add the columns together 
        if (data$CpG1 == data$CpG2) {data$CpG <- data$CpG1} else {data$CpG <- data$CpG1 + data$CpG2}

        print(head(data))
        print(str(data))
        
	##define read count as 1 per line
	data$read_count <- 1

	### remove uncharacterized genomic regions
	data <- as.data.frame(data[!grepl("*KI*", data$chr), ])
	data <- as.data.frame(data[!grepl("*GL*", data$chr), ])
	data <- as.data.frame(data[!grepl("*EBV*", data$chr), ])
	data <- as.data.frame(data[!grepl("*Un*", data$chr), ])

	data_melt <- melt(data,
                  id.vars = c("GC","fragment_len", "CpG", "chr", "start", "end"),
                  measure.vars = "read_count")
	data_melt$value <- as.numeric(as.character(data_melt$value))

	data_melt$CpG_3 <- (data_melt$CpG) ^ (1 / 3)

	data_melt$fragment_len <- as.numeric(data_melt$fragment_len)
        data_melt$GC <- as.numeric(data_melt$GC)
        data_melt$CpG <- as.numeric(data_melt$CpG)

        stopifnot(str(data_melt$fragment_len) == "num")
        stopifnot(str(data_melt$GC) == "num")
        stopifnot(str(data_melt$CpG) == "num")

	##Source the GLM from the control GLM code
	load("/mnt/work1/users/home2/wilsons/Projects/2018_PTB/data/2019_ControlAlignment/2020_Gaussian_0.01ng.rda")

	print("start pmol prediction")

	##Predict fmol values
	pred_conc <- as.data.frame(predict(gaussian_glm, data_melt))
	names(pred_conc) <- "value"
	
	 ##Replace value with conc
        data_melt$read_count<- data_melt$value 
        data_melt$value <- NULL
        data_melt <- cbind(data_melt, pred_conc)
        data_melt$variable <- NULL
        colnames(data_melt) <- c("GC","fragment_len", "CpG", "chr", "start", "end", "CpG_3", "read_count", "pmol")

	print(head(data_melt))

	#restructure to bed file format
        data_melt <- data_melt[, c("chr", "start", "end", "fragment_len", "GC", "CpG", "CpG_3", "read_count", "pmol")]
        print(head(data_melt))

	write.table(data_melt, spe = "\t", row.names = FALSE, quote = FALSE, file = paste0("~/Projects/2018_PTB/data/2019_ControlAlignment/",
                                       file_path_sans_ext(basename(filename)), "_pmol.bed"))

}

args <-
if (length(commandArgs(TRUE))
    || commandArgs()[length(commandArgs())] == "--args") {
    as.list(commandArgs(TRUE))
} else {
    list()
}

 do.call(main, args)

