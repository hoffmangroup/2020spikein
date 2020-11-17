#!/usr/bin/env Rscript

##Load needed packages

library("qsea")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
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
library(compute.es)
library(grid)

###Heat-scree plots are adapted from:
#### De Souza, Rebecca AG, et al.
###"DNA methylation profiling in human Huntington's disease brain."
###Human molecular genetics 25.10 (2016): 2013-2030.

##ggplot2 theme
theme_set(theme_bw() +
          theme(axis.title.y = element_text(size = 20),
      axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 16),
      axis.text.y = element_text(vjust = 0.5, size = 16),
      axis.title.x = element_text(vjust = 0.5, size = 20),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      strip.text.x = element_text(size = 16),
      plot.title =
        element_text(hjust = 0.5, color = "black", size = 16, face = "bold")))

#Bam files for this project
PATH_PREFIX_B1 <- "~/Projects/2018_PTB/data/2020_BatchAnalysis/Batch1_RS/"
BamsList_B1 <- list.files(path = PATH_PREFIX_B1, pattern = "filtered_human.*trimmed.*.bam")
bamFile_B1 <- data.frame(paste0(PATH_PREFIX_B1, BamsList_B1))
colnames(bamFile_B1) <- "file_name"

PATH_PREFIX_B2 <- "~/Projects/2018_PTB/data/2020_BatchAnalysis/Batch2_JB/"
BamsList_B2 <- list.files(path = PATH_PREFIX_B2, pattern = "filtered_human.*trimmed.*.bam")
bamFile_B2 <- data.frame(paste0(PATH_PREFIX_B2, BamsList_B2))
colnames(bamFile_B2) <- "file_name"

PATH_PREFIX_B3 <- "~/Projects/2018_PTB/data/2020_BatchAnalysis/Batch3_DT/"
BamsList_B3 <- list.files(path = PATH_PREFIX_B3, pattern = "filtered_human.*trimmed.*.bam")
bamFile_B3 <- data.frame(paste0(PATH_PREFIX_B3, BamsList_B3))
colnames(bamFile_B3) <- "file_name"

bamFile <- rbind(bamFile_B1, bamFile_B2, bamFile_B3)

#Create Qsea set for each sample and bind together in a list
sample_table <- data.frame(sample_name = paste0("sample_", 1:15),
                          file_name = bamFile$file_name,
                          group = c(rep("B1", 5), rep("B2", 5), rep("B3", 5)),
                          stringsAsFactors = FALSE)
	
#Defining parameters again
BSgenome <- "BSgenome.Hsapiens.UCSC.hg38"
ws <- 300
chr.select <- paste0("chr", 1:22)
	
#Create Qsea set for each sample, this will be combined with the other sample MEDIPS set into a list
qseaSet <- createQseaSet(sampleTable = sample_table, BSgenome = BSgenome,
                         window_size = ws, chr.select = chr.select)
qseaSet <- addCoverage(qseaSet, uniquePos = TRUE,
                       paired = TRUE)
qseaSet <- addPatternDensity(qseaSet, "CG", name = "CpG",
                             fragment_length = 300)
qseaSet <- addLibraryFactors(qseaSet)
qseaSet <- addOffset(qseaSet, enrichmentPattern = "CpG")

str(qseaSet)#Check that SetCreated.list is infact a list of N (number of samples in data set) MEDIPS sets

#Specify the shortened sample names here
save(qseaSet, file = "2020_allbatches_QseaSets.RData")
#save the MEDIPS sets

##PCA analysis
pca <- prcomp(t(qseaSet@count_matrix))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)
rotations_qsea <- as.data.frame(pca$rotation)

##Create in meta information
loadings$batch <- c(rep("Batch_1", 5), rep("Batch_2", 5), rep("Batch_3", 5))
loadings$sample <- c("sample 1", "sample 2", "sample 3", "sample 4", "sample 5",
                     "sample 1", "sample 2", "sample 3", "sample 4", "sample 5",
                     "sample 1", "sample 2", "sample 3", "sample 4", "sample 5")
loadings$batch <- as.factor(loadings$batch)

png("2020_PCA_qsea_allbatches.png", height = 6, width = 6, units = "in", res = 300)
QSEA <- ggplot(loadings, aes(PC1, PC2, color = batch, shape = sample)) +
  geom_point(size = 2) +
  xlab("PC1 (74.9%)") + ylab("PC2 (11.3%)") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("#67a9cf", "#ef8a62", "#b2182b"), name = "Batch") +
  scale_shape_manual(values = c(16, 17, 18, 15, 25), name = "Sample")
QSEA
dev.off()

###Load in the fmol data and merge together
batch1 <- read.table("2020_human_batch1_fmol.bed", header = TRUE)
batch1$chrom <- NULL
batch1$chromStart <- NULL
batch1$chromEnd <- NULL

batch2 <- read.table("2020_human_batch2_fmol.bed",
                     header = TRUE, sep = ",")
batch2$X <- NULL
batch2$chrom <- NULL
batch2$chromStart <- NULL
batch2$chromEnd <- NULL

batch3 <- read.table("2020_human_batch3_fmol.bed",
                     header = TRUE, sep = ",")
batch3$X <- NULL
batch3$chrom <- NULL
batch3$chromStart <- NULL
batch3$chromEnd <- NULL

data_pmol <- merge(batch1, batch2, by = "window", all = TRUE)
data_pmol <- merge(data_pmol, batch3, by = "window", all = TRUE)
colnames(data_pmol) <- c("window", "S1_B1", "S2_B1", "S3_B1", "S4_B1", "S5_B1",
                         "S1_B2", "S2_B2", "S3_B2", "S4_B2", "S5_B2",
                         "S1_B3", "S2_B3", "S3_B3", "S4_B3", "S5_B3")

##PCA analysis
pca_dat <- data_pmol
pca_dat[is.na(pca_dat)] <- 0
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)
rotations_qsea <- as.data.frame(pca$rotation)
##Create in meta information
loadings$batch <- c(rep("Batch_1", 5), rep("Batch_2", 5), rep("Batch_3", 5))
loadings$sample <- c("sample 1", "sample 2", "sample 3", "sample 4", "sample 5",
                     "sample 1", "sample 2", "sample 3", "sample 4", "sample 5",
                     "sample 1", "sample 2", "sample 3", "sample 4", "sample 5")
loadings$batch <- as.factor(loadings$batch)

png("2020_PCA_fmol_allbatches.png", height = 6, width = 6, units = "in", res = 300)
pmol <- ggplot(loadings, aes(PC1, PC2, color = batch, shape = sample)) +
  geom_point(size = 2) +
  xlab("PC1 (91.7%)") + ylab("PC2 (4.2%)") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("#67a9cf", "#ef8a62", "#b2182b"), name = "Batch") +
  scale_shape_manual(values = c(16, 17, 18, 15, 25), name = "Sample")
pmol
dev.off()

loadings_pmol <- loadings

###Load in the fmol data and merge together
batch1_raw <- read.table("2020_human_batch1_rc.bed", header = TRUE)
batch1_raw$chrom <- NULL
batch1_raw$chromStart <- NULL
batch1_raw$chromEnd <- NULL

batch2_raw <- read.table("2020_human_batch2_rc.bed",
                         header = TRUE, sep = ",")
batch2_raw$X <- NULL
batch2_raw$chrom <- NULL
batch2_raw$chromStart <- NULL
batch2_raw$chromEnd <- NULL

batch3_raw <- read.table("2020_human_batch3_rc.bed",
                         header = TRUE, sep = ",")
batch3_raw$X <- NULL
batch3_raw$chrom <- NULL
batch3_raw$chromStart <- NULL
batch3_raw$chromEnd <- NULL

data_raw <- merge(batch1_raw, batch2_raw, by = "window")
data_raw <- merge(data_raw, batch3_raw, by = "window")
colnames(data_raw) <- c("window", "S1_B1", "S2_B1", "S3_B1", "S4_B1", "S5_B1",
                        "S1_B2", "S2_B2", "S3_B2", "S4_B2", "S5_B2",
                        "S1_B3", "S2_B3", "S3_B3", "S4_B3", "S5_B3")

##PCA analysis
pca_dat <- data_raw
pca_dat[is.na(pca_dat)] <- 0
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)
##Create in meta information
loadings$batch <- c(rep("Batch_1", 5), rep("Batch_2", 5), rep("Batch_3", 5))
loadings$sample <- c("sample 1", "sample 2", "sample 3", "sample 4", "sample 5",
                     "sample 1", "sample 2", "sample 3", "sample 4", "sample 5",
                     "sample 1", "sample 2", "sample 3", "sample 4", "sample 5")
loadings$batch <- as.factor(loadings$batch)

png("2020_PCA_raw_allbatches.png", height = 6, width = 6, units = "in", res = 300)
raw <- ggplot(loadings, aes(PC1, PC2, color = batch, shape = sample)) +
  geom_point(size = 2) +
  xlab("PC1 (77.9%)") + ylab("PC2 (13.7%)") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("#67a9cf", "#ef8a62", "#b2182b"), name = "Batch") +
  scale_shape_manual(values = c(16, 17, 18, 15, 25), name = "Sample")
raw
dev.off()

Loadings_raw <- loadings
###Put all the PCA on one grid plot
###Code from https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots

grid_arrange_shared_legend <- function(..., nrow = 1,
                                       ncol = length(list(...)),
                                       position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

}

png("2020_PCA_Raw_QSEA_fmol.png", height = 6, width = 18, units = "in", res = 300)
grid_arrange_shared_legend(raw + ggtitle("Raw reads"), QSEA + ggtitle("QSEA reads"), pmol + ggtitle("Picomoles"))
dev.off()

data_pmol_raw <- data_pmol
data_pmol_raw[is.na(data_pmol_raw)] <- 0
Missing <- lapply(data_pmol_raw, function(x) {
  length(which(x == 0)) / length(x)})

##Remove repetitive elements- simple repeats only
###UCSC simple repeats
repeats <- read.table("~/Projects/2018_PTB/data/2019_EPIC_HCT116_CtlProject/2020_repeatregions_mappedto300bpwindows.bed", sep = "\t", header = FALSE)
###subset to regions where there are overlaps between our windows
repeats <- subset(repeats, repeats$V9 > 0)
repeats <- repeats[, c(1:3)]
colnames(repeats) <- c("chr", "start", "end")
repeats$window <- paste0(repeats$chr, "_", repeats$start, "_", repeats$end)
repeats <- unique(repeats)
rownames(repeats) <- repeats$window

##remove repetitive regions from the fmol data
rownames(data_pmol) <- data_pmol$window
data_pmol <- data_pmol[!rownames(data_pmol) %in% rownames(repeats), ]

##remove blacklist regions
###ENCODE blacklist regions
blacklist <- read.table("~/Projects/2018_PTB/data/Annotations/2020_new_ENCODE_blacklist.csv",
                        sep = "\t", header = FALSE)
colnames(blacklist) <- c("bl_chr", "bl_start", "bl_end", "chr", "start", "end", "overlap")
blacklist$window <- paste0(blacklist$chr, "_", blacklist$start, "_", blacklist$end)
blacklist <- unique(blacklist)

##remove blacklist regions from fmol
data_pmol <- data_pmol[!rownames(data_pmol) %in% blacklist$window, ]

## remove low mappability regions (mappability score<0.5)
mappability <- read.csv("~/Projects/2018_PTB/data/2019_ControlAlignment/2020_min_mappability_300bp_windows.csv",
                        header = TRUE)
mappability$X <- NULL

mappability_05 <- subset(mappability, mappability$map_score < 0.5)

data_pmol <- data_pmol[!rownames(data_pmol) %in% mappability_05$window, ]

##PCA analysis
pca_dat <- data_pmol
pca_dat[is.na(pca_dat)] <- 0
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)
        xlab("PC1 (74.9%)") +
          ylab("PC2 (11.3%)") +
          scale_x_continuous(labels = comma) +
          scale_y_continuous(labels = comma) +
          scale_color_manual(values = c("#67a9cf", "#ef8a62", "#b2182b"), name = "Batch") +
          scale_shape_manual(values = c(16, 17, 18, 15, 25), name = "Sample")
QSEA
dev.off()

###Load in the fmol data and merge together
batch1 <- read.table("2020_human_batch1_pmol.bed", header = TRUE)
batch1$chrom <- NULL
batch1$chromStart <- NULL
batch1$chromEnd <- NULL

batch2 <- read.table("2020_human_batch2_pmol.bed",
                     header = TRUE, sep = ",")
batch2$X <- NULL
batch2$chrom <- NULL
batch2$chromStart <- NULL
batch2$chromEnd <- NULL

batch3 <- read.table("2020_human_batch3_fmol.bed",
                     header = TRUE, sep = ",")
batch3$X <- NULL
batch3$chrom <- NULL
batch3$chromStart <- NULL
batch3$chromEnd <- NULL

data_pmol <- merge(batch1, batch2, by = "window", all = TRUE)
data_pmol <- merge(data_pmol, batch3, by = "window", all = TRUE)
colnames(data_pmol) <- c("window", "S1_B1", "S2_B1", "S3_B1", "S4_B1", "S5_B1",
                         "S1_B2", "S2_B2", "S3_B2", "S4_B2", "S5_B2",
                         "S1_B3", "S2_B3", "S3_B3", "S4_B3", "S5_B3")

##PCA analysis
pca_dat <- data_pmol
pca_dat[is.na(pca_dat)] <- 0
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)
rotations_qsea <- as.data.frame(pca$rotation)
##Create in meta information
loadings$batch <- c(rep("Batch_1", 5), rep("Batch_2", 5), rep("Batch_3", 5))
loadings$sample <- c("sample 1", "sample 2", "sample 3", "sample 4", "sample 5",
                     "sample 1", "sample 2", "sample 3", "sample 4", "sample 5",
                     "sample 1", "sample 2", "sample 3", "sample 4", "sample 5")
loadings$batch <- as.factor(loadings$batch)

png("2020_PCA_fmol_allbatches.png", height = 6, width = 6, units = "in", res = 300)
pmol <- ggplot(loadings, aes(PC1, PC2, color = batch, shape = sample)) +
  geom_point(size = 2) +
  xlab("PC1 (91.7%)") + ylab("PC2 (4.2%)") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("#67a9cf", "#ef8a62", "#b2182b"), name = "Batch") +
  scale_shape_manual(values = c(16, 17, 18, 15, 25), name = "Sample")
pmol
dev.off()

loadings_pmol <- loadings

###Load in the fmol data and merge together
batch1_raw <- read.table("2020_human_batch1_rc.bed", header = TRUE)
batch1_raw$chrom <- NULL
batch1_raw$chromStart <- NULL
batch1_raw$chromEnd <- NULL

batch2_raw <- read.table("2020_human_batch2_rc.bed",
                         header = TRUE, sep = ",")
batch2_raw$X <- NULL
batch2_raw$chrom <- NULL
batch2_raw$chromStart <- NULL
batch2_raw$chromEnd <- NULL

batch3_raw <- read.table("2020_human_batch3_rc.bed",
                         header = TRUE, sep = ",")
batch3_raw$X <- NULL
batch3_raw$chrom <- NULL
batch3_raw$chromStart <- NULL
batch3_raw$chromEnd <- NULL

##Infer sex from Y signal in raw data
chry_reads <- subset(data_raw, str_detect(data_raw$window, pattern = "chrY") == "TRUE")
chry_reads[is.na(chry_reads)] <- 0
no_chrY <- lapply(chry_reads, function(x) {
  length(which(x == 0)) / length(x)})

pca <- prcomp(t(qseaSet@count_matrix))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)

###Heat scree effect size
heat_scree_plot_es <- function(Loadings, Importance, Num, Order) {
  pca_adjusted <- importance[1:length(importance)]
  pca_df <- data.frame(adjusted_variance = pca_adjusted, PC = seq(1:length(pca_adjusted)))

  scree <- ggplot(pca_df[which(pca_df$PC < Num), ],
                  aes(PC, adjusted_variance)) +
    theme_bw() +
    geom_segment(aes(x = PC, xend = PC, y = 0, yend = adjusted_variance), color = "black", size = 1) +
    geom_point(size = 3, color = "black", fill = alpha("black", 1.0), alpha = 1.0, shape = 21, stroke = 2) +
    geom_text(label = round(pca_df[which(pca_df$PC < Num), ]$adjusted_variance, 2), vjust = -1) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          plot.margin = unit(c(1, 1.5, 0.2, 2.25), "cm")) +
    ylab("Variance") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    scale_x_discrete(limits = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")) +
    xlab("Principal component")

  #### Heat
  ## correlate meta with PCS
  ## Run anova of each PC on each meta data variable
  aov_PC_meta <- lapply(1:ncol(meta_categorical),
                        function(covar) sapply(1:ncol(Loadings),
                                               function(PC) summary(aov(Loadings[, PC] ~ meta_categorical[, covar]))[[1]]$"F value"[1]))
  names(aov_PC_meta) <- colnames(meta_categorical)
  aov_PC_meta <- do.call(rbind, aov_PC_meta)
  aov_PC_meta <- as.data.frame(aov_PC_meta)

  #adjust
  aov_PC_meta_adjust <- aov_PC_meta[, 2:ncol(aov_PC_meta)]

  ##put in the sample size for each variable
  aov_PC_meta_adjust$size <- 0
  for (i in rownames(aov_PC_meta_adjust)) {
  aov_PC_meta_adjust[i, "size"] <- length(levels(meta_categorical[[i]]))
  }

##converting F-statistic to effect size
 meandis <- array(data = NA, dim = dim(aov_PC_meta_adjust))
 rownames(meandis) <- rownames(aov_PC_meta_adjust)
 meanp <- meandis
 for (i in 1:(ncol(meandis) - 1)) {
 for (j in 1:nrow(meandis)) {
 meandis[j, i] <- as.numeric(fes(aov_PC_meta_adjust[j ,i],
                                aov_PC_meta_adjust[j, "size"],
                                length(Loadings), verbose = FALSE)["d"])
 }
 }

 for (i in 1:(ncol(meandis)-1)) {
 for (j in 1:nrow(meandis)) {
 meanp[j, i] <- as.numeric(fes(aov_PC_meta_adjust[j, i],
                              aov_PC_meta_adjust[j, "size"],
                              length(Loadings), verbose = FALSE)["pval.d"])
 }
 }

  #reshape - mean distance
  avo <- meandis[, 1:(Num - 1)]
  avo_heat_num <- apply(avo, 2, as.numeric)
  avo_heat <- as.data.frame(avo_heat_num)
  colnames(avo_heat) <- sapply(1:(Num - 1),
                             function(x) paste("PC", x, sep = ""))
  avo_heat$meta <- rownames(avo)
  avo_heat_melt <- melt(avo_heat, id = c("meta"))

  # cluster meta data
  ord <- Order
  meta_var_order <- unique(avo_heat_melt$meta)[rev(ord)]
  avo_heat_melt$meta <- factor(avo_heat_melt$meta, levels = meta_var_order)

  ##reshape - distance p value
  avo_p <- meanp[, 1:(Num - 1)]
  avo_heat_num_p <- apply(avo_p, 2, as.numeric)
  avo_heat_p <- as.data.frame(avo_heat_num_p)
  colnames(avo_heat_p) <- sapply(1:(Num - 1), function(x) paste("PC", x, sep = ""))
  avo_heat_p$meta <- rownames(avo_p)
  avo_heat_melt_p <- melt(avo_heat_p, id = c("meta"))

  # cluster meta data
  meta_var_order_p <- unique(avo_heat_melt_p$meta)[rev(ord)]
  avo_heat_melt_p$meta <- factor(avo_heat_melt_p$meta, levels = meta_var_order_p)
 avo_heat_melt_p$p_value <- avo_heat_melt_p$value
 avo_heat_melt_p$value <- NULL

avo_heat_melt <- merge(avo_heat_melt, avo_heat_melt_p, by = c("meta", "variable"))

##adjust p for multipl test correction
avo_heat_melt$adj_p <- p.adjust(avo_heat_melt$p_value, method = "holm", n = nrow(avo_heat_melt))

##labelling significance for heatmap
avo_heat_melt$stars <- cut(avo_heat_melt$adj_p,
                           breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))

  heat <- ggplot(avo_heat_melt, aes(variable, meta, fill = value)) +
  geom_tile(color = "black", size = 0.5) +
  theme_bw() +
    scale_fill_gradient(limits = c(0, 5), high = "#2c7fb8", low = "#edf8b1") +
      theme(axis.text = element_text(size = 10, color = "black"),
            axis.text.x = element_text(),
          axis.title = element_text(size = 15),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 12),
          legend.position = c(1, 0), legend.justification = c(1, 0),
          plot.margin = unit(c(0, 2.25, 1, 1), "cm")) +
        labs(fill = "Cohen's d") +
    xlab("Principal component") +
    ylab(NULL) +
    geom_text(aes(label = stars), color = "black", size = 5)

  plot_grid(scree, heat, ncol = 2, align = "h")

}

##Label to variables
##Create in meta information
Loadings_meta <- data.frame(batch = as.factor(c(rep("Batch_1", 5), rep("Batch_2", 5), rep("Batch_3", 5))),
                            sample = rep(c("sample 1", "sample 2", "sample 3", "sample 4", "sample 5"), 3),
                            sequencer = as.factor(c(rep("A", 5), rep("A", 5), rep("B", 5))),
                            adapters = as.factor(c(rep("A", 5), rep("B", 5), rep("A", 5))),
                            sex = as.factor(rep(c("F", "F", "M", "M", "F"), 3)))

# Specifiy the number of PCs you want shown
Num <- 9
# should be equal to the number of samples in your dataset; for large datasets, you can opt to just see the top PCs
# Designate what order you want the variables to appear (continuous variables rbinded to categorical variables in function)
Order <- c(1, 2, 3, 4, 5)

meta_categorical <- Loadings_meta[, c("batch", "sequencer", "adapters", "sample", "sex")]
# input column numbers in meta that contain categorical variables
colnames(meta_categorical) <- c("Batch", "Sequencing Machine", "Adapters", "Sample", "Sex")

##QSEA
pca <- prcomp(t(qseaSet@count_matrix))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)

##QSEA plot
png("2020_heatscree_meand_qsea.png", height = 6, width = 12, unit = "in", res = 300)
qsea <- heat_scree_plot_es(loadings, importance, Num, Order) +
  theme(legend.position = c(-1, 0))
qsea
dev.off()

Loadings_qsea <- loadings

##pmol filtered
pca_dat <- data_pmol
pca_dat[is.na(pca_dat)] <- 0
rownames(pca_dat) <- pca_dat$window
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)
##3% of variance associated with batch

loadings <- as.data.frame(pca$x)

##fmol filtered plot
png("2020_heatscree_meand_fmolfiltered.png", height = 6, width = 6 , unit = "in", res = 300)
pmol <- heat_scree_plot_es(loadings, importance, Num, Order) +
  theme(legend.position = c(-1, 0))
pmol
dev.off()

loadings_pmol <- loadings

##Raw data
pca_dat <- data_raw
pca_dat[is.na(pca_dat)] <- 0
rownames(pca_dat) <- pca_dat$window
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)

png("2020_heatscree_meand_raw.png", height = 6, width = 6, unit = "in", res = 300)
raw <- heat_scree_plot_es(loadings, importance, Num, Order) +
  theme(legend.position = c(-1, 0))
raw
dev.off()

Loadings_raw <- loadings

##fpmol not filtered

pca_dat <- data_pmol_raw
pca_dat[is.na(pca_dat)] <- 0
rownames(pca_dat) <- pca_dat$window
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev ^ 2
importance <- vars / sum(vars)

loadings <- as.data.frame(pca$x)

png("2020_heatscree_meand_allwindows_fmol.png", height = 6, width = 6, unit = "in", res = 300)
pmol_raw <- heat_scree_plot_es(loadings, importance, Num, Order) +
  theme(legend.position = c(-1, 0))
pmol_raw
dev.off()

Loadings_fmol_raw <- loadings

png("2020_Heatscree_allbatches_meand_pmolandpmolnoblacklist.png", height = 24, width = 12, units = "in", res = 300)
grid.arrange(arrangeGrob(raw + theme(axis.title.x = element_text(""), plot.margin = unit(c(6, 6, 6, 6)))),
        arrangeGrob(qsea + theme(axis.title.x = element_text("")), plot.margin = unit(c(6, 6, 6, 6))),
        arrangeGrob(pmol_raw  + theme(axis.title.x = element_text("")), plot.margin = unit(c(6, 6, 6, 6))),
        arrangeGrob(pmol  + theme(axis.title.x = element_text("")), plot.margin = unit(c(6, 6, 6, 6))), ncol = 1)
dev.off()