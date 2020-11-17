#!/usr/bin/env Rscript

##Load libraries
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
BSgenome <- "BSgenome.Hsapiens.UCSC.hg38"
library("grDevices")
library("doParallel")
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)
library(stringr)
library(cowplot)
library(reshape2)
library(plyr)
library(data.table)

##load appris data
appris <- fread("~/Projects/2018_PTB/data/Annotations/2020_appris_5000bp_genename_overlap_300bp.bed")
colnames(appris) <- c("chr", "ROI_start", "ROI_end", "txp_start",
                      "txp_end", "ENSEMBL", "SYMBOL", "chrom",
                      "GENE_NAME", "window_start", "window_end",
                      "overlap")
appris$chrom <- NULL
appris$window <- paste0(appris$chr, "_", appris$window_start, "_", appris$window_end)


##Load in rotation data from PCA
PCs <- fread("2020_fmol_filtered_rotations.csv")
colnames(PCs) <- c("window", "PC1", "PC2", "PC3", "PC4", "PC5",
                   "PC6", "PC7", "PC8", "PC9", "PC10", "PC11",
                   "PC12", "PC13", "PC14", "PC15")

##PC1 drivers
### Take top 10%
PC1 <- PCs[, c("window", "PC1")]
PC1$PC1 <- abs(PC1$PC1)
PC1_ordered <- PC1[order(-PC1, PC1), ]
PC1_top <- PC1_ordered[1:(round(nrow(PC1_ordered) * 0.1)), ]

### merge with genenames from appris
PC1_top_genes <- merge(PC1_top, appris, by = "window")
PC1_top_genes <- PC1_top_genes[order(PC1_top_genes, -PC1), ]

##get rid of ROI/txp start and ends, so I can aggregate by window
PC1_top_genes <- PC1_top_genes[, c("window", "PC1", "ENSEMBL",
                                   "SYMBOL", "GENE_NAME", "overlap")]
PC1_top_genes <- PC1_top_genes[!duplicated(PC1_top_genes), ]
sum(PC1_top_genes$GENE_NAME == "chrX" | PC1_top_genes$GENE_NAME == "chrY") / nrow(PC1_top_genes)
# Only 2% on sex chromosomes

write.table(PC1_top_genes, file = "2020_fmol_filtered_PC1_drivers.csv")

##PC2 drivers
### Take top 10%
PC2 <- PCs[, c("window", "PC2")]
PC2$PC2 <- abs(PC2$PC2)
PC2_ordered <- PC2[order(-PC2, PC2), ]
PC2_top <- PC2_ordered[1:(round(nrow(PC2_ordered) * 0.1)), ]

### merge with genenames from appris
PC2_top_genes <- merge(PC2_top, appris, by = "window")
PC2_top_genes <- PC2_top_genes[order(PC2_top_genes, -PC2), ]

##get rid of ROI/txp start and ends, so I can aggregate by window
PC2_top_genes <- PC2_top_genes[, c("window", "PC2", "ENSEMBL",
                                  "SYMBOL", "GENE_NAME", "overlap")]
PC2_top_genes <- PC2_top_genes[!duplicated(PC2_top_genes), ]
sum(PC2_top_genes$GENE_NAME == "chrX" | PC2_top_genes$GENE_NAME == "chrY") / nrow(PC2_top_genes)
# Only 2% on sex chromosomes

write.table(PC2_top_genes, file = "2020_fmol_filtered_PC2_drivers.csv")

##How many overlap between PC1 and PC2
sum(PC1_top_genes$window %in% PC2_top_genes$window)
##18075, 66% of them overlap