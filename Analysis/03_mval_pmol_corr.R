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
library(tidyverse)
library(scales)
library(ggpubr)
library("org.Hs.eg.db")

#Load data

##Predicted concentration metrics
conc <- read.delim("~/Projects/2018_PTB/data/2019_ControlAlignment/2020_human0.01_adj_binned.csv",
                   sep = ",", header = TRUE)
conc$X <- NULL
colnames(conc) <- c("chr", "start", "end", "sample1", "sample2")
conc$window <- paste0(conc$chr, "_", conc$start, "_", conc$end)
rownames(conc) <- conc$window

epic <- read.delim("2019_HCT116_sesame_betas_windows.csv",
                   sep = ",", header = TRUE)
epic$X <- NULL
epic$CpG_beg <- NULL
epic$CpG_end <- NULL
colnames(epic) <- c("window", "Sample1", "Sample2", "Sample3")

##Remove repetitive elements- simple repeats only
repeats <- read.table("2020_repeatregions_mappedto300bpwindows.bed",
                      sep = "\t", header = FALSE)
###subset to regions where there are overlaps between our windows
repeats <- subset(repeats, repeats$V9 > 0)
repeats <- repeats[, c(1:3)]
colnames(repeats) <- c("chr", "start", "end")
repeats$window <- paste0(repeats$chr, "_", repeats$start, "_", repeats$end)
repeats <- unique(repeats)
rownames(repeats) <- repeats$window

##remove repetitive regions from the data
conc <- conc[!rownames(conc) %in% rownames(repeats), ]

##remove blacklist regions
blacklist <- read.csv("~/Projects/2018_PTB/data/Annotations/2020_new_ENCODE_blacklist.csv",
                      header = FALSE)
colnames(blacklist) <- c("bl_chr", "bl_start", "bl_end",
                         "chr", "start", "end", "overlap")
blacklist$window <- paste0(blacklist$chr, "_", blacklist$start, "_", blacklist$end)
blacklist <- unique(blacklist)

##remove blacklist regions
conc <- conc[!rownames(conc) %in% blacklist$window, ]

## remove low mappability regions (mappability score<0.5)
mappability <- read.csv("~/Projects/2018_PTB/data/2019_ControlAlignment/2020_min_mappability_300bp_windows.csv",
                        header = TRUE)
mappability$X <- NULL

mappability_05 <- subset(mappability, mappability$map_score < 0.5)

conc <- conc[!rownames(conc) %in% mappability_05$window, ]
#4551870 windows

#subset out samples with sd >0.25
conc$sd <- apply(conc[, c("sample1", "sample2")], 1, sd)

conc <- conc[!conc$sd > 0.25, ]

##subset conc and epic to only those windows present in both
windows_of_interest <- epic$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(epic) <- epic$window
epic_subset <- subset(epic, rownames(epic) %in% rownames(conc))

##mean the samples
conc_means <- data.frame(ID = conc_subset$window,
                         Means = rowMeans(conc_subset[, c("sample1", "sample2")]))
epic_means <- data.frame(ID = epic_subset$window,
                         Means = rowMeans(epic_subset[, c("Sample1", "Sample2", "Sample3")]))

means <- merge(conc_means, epic_means, by = "ID")
colnames(means) <- c("window", "conc", "epic")

##changes betas to m-values
means$epic_m <- log2(means$epic / (1 - means$epic))

##Correlation calculation
cor.test(means$conc, means$epic_m, method = "spearman")
#0.2120731

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

##plot function
corr_plot <- function(means) {
ggplot(means, aes(x = epic_m, y = conc)) +
  geom_point(color = "#2980B9", size = 2) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, color = "#2C3E50") +
  ylab("Picomoles") +
  xlab("DNA methylation (M-value)") +
  scale_y_continuous(labels = comma)
}

##Plot
png("2020_conc_epic_sesame_mapp_corr.png",
    height = 6, width = 6, units = "in", res = 300)
all_pts <- corr_plot(means)
all_pts
dev.off()

##Correlate bins with >3 CpG represented on array
##Let's correlate those this more than 3 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg3.csv",
                     row.names = 1, header = TRUE)
#50053 windows

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID = conc_subset$window,
                              Means = rowMeans(conc_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "conc", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))

##Correlation calculation
cor.test(high_means$conc, high_means$epic_m, method = "spearman")
##0.6183018    - with 3 CpG and 50053 windows

##Full plot
png("2020_conc_epicsesame_corr_mapp_morethan3CpG.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more3cpg <- corr_plot(high_means) +
                 scale_x_continuous(limits = c(-8, 8),
                                    breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
                 scale_y_continuous(limits = c(0, 500),
                                    breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500))
corr_more3cpg
dev.off()

##Plot- cut to see the small correlation
png("2020_conc_epicsesame_corr_mapp_morethan3CpG_cut.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more3cpg_cut <- corr_plot(high_means) +
  scale_y_continuous(labels = comma, limits = c(0, 0.5)) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6))

corr_more3cpg_cut
dev.off()

##Let's correlate those this more than 5 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg5.csv",
                     row.names = 1, header = TRUE)
#12602 windows

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID = conc_subset$window,
                              Means = rowMeans(conc_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "conc", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))

##Correlation calculation
cor.test(high_means$conc, high_means$epic_m, method = "spearman")
##0.7950763    - with more than 5 CpG and 12602 windows

##Full plot
png("2020_conc_epicsesame_corr_mapp_morethan5CpG.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more5cpg <- corr_plot(high_means) +
                 scale_x_continuous(limits = c(-8, 8),
                                    breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
                 scale_y_continuous(limits = c(0, 500),
                                    breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500))
corr_more5cpg
dev.off()

##Plot- cut to see the small correlation
png("2020_conc_epic_corrsesame_mapp_morethan5CpG_cut.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more5cpg_cut <- corr_plot(high_means) +
  scale_y_continuous(labels = comma, limits = c(0, 0.5)) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6))
corr_more5cpg_cut
dev.off()

### 7 CpG
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg7.csv",
                     row.names = 1, header = TRUE)
#3567 windows

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID = conc_subset$window,
                              Means = rowMeans(conc_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "conc", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))

##Correlation calculation
cor.test(high_means$conc, high_means$epic_m, method = "spearman")
##0.8153412   - with more than 7 CpG and 3567 windows

##Full plot
png("2020_conc_epicsesame_corr_mapp_morethan7CpG.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more7cpg <- corr_plot(high_means) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
                 scale_y_continuous(limits = c(0, 500),
                                    breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500))
corr_more7cpg
dev.off()

##Plot- cut to see the small correlation
png("2020_conc_epicsesame_corr_mapp_morethan7CpG_cut.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more7cpg_cut <- corr_plot(high_means) +
  scale_y_continuous(labels = comma, limits = c(0, 0.5)) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6))
corr_more7cpg_cut
dev.off()

###10 CpG
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg10.csv",
                     row.names = 1, header = TRUE)
#294 windows

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID = conc_subset$window,
                              Means = rowMeans(conc_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "conc", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))

##Correlation calculation
cor.test(high_means$conc, high_means$epic_m, method = "spearman")
##0.8534883   - with more than 10 CpG and 294 windows

##Full plot
png("2020_conc_epicsesame_corr_mapp_morethan10CpG.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more10cpg <- corr_plot(high_means) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
                 scale_y_continuous(limits = c(0, 500),
                                    breaks = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500))
corr_more10cpg
dev.off()

##Plot- cut to see the small correlation
png("2020_conc_epicsesame_corr_mapp_morethan10CpG_cut.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more10cpg_cut <- corr_plot(high_means)  +
  scale_y_continuous(labels = comma, limits = c(0, 0.5)) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6))
corr_more10cpg_cut
dev.off()

##arrange all correlations to fmol together
png("2020_fmol_mvalsesame_corr_mapp_allplts.png",
    height = 12, width = 8, units = "in", res = 300)
all_fmol_corr_plot <- grid.arrange(arrangeGrob(corr_more3cpg + theme(axis.title.x = element_blank()),
                                               top = "More than 3 CpG"),
                            arrangeGrob(corr_more3cpg_cut + theme(axis.title.x = element_blank(),
                                                                  axis.title.y = element_blank()),
                                                                  top = "More than 3 CpG: zoom"),
                            arrangeGrob(corr_more5cpg + theme(axis.title.x = element_blank()),
                                        top = "More than 5 CpG"),
                            arrangeGrob(corr_more5cpg_cut + theme(axis.title.x = element_blank(),
                                                                  axis.title.y = element_blank()),
                                        top = "More than 5 CpG: zoom"),
                            arrangeGrob(corr_more7cpg + theme(axis.title.x = element_blank()),
                                        top = "More than 7 CpG"),
                            arrangeGrob(corr_more7cpg_cut + theme(axis.title.x = element_blank(),
                                                                  axis.title.y = element_blank()),
                                        top = "More than 7 CpG: zoom"),
                            arrangeGrob(corr_more10cpg + theme(axis.title.x = element_blank()),
                                        top = "More than 10 CpG"),
                            arrangeGrob(corr_more10cpg_cut + theme(axis.title.x = element_blank(),
                                                                   axis.title.y = element_blank()),
                                        top = "More than 10 CpG: zoom"), ncol = 2)

annotate_figure(all_fmol_corr_plot,
               bottom = text_grob("DNA methylation (M-value)", color = "black",
                                   size = 16))
dev.off()

##Let's look at what these very high fmol sites are (using the 3 CpG data)
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg3.csv",
                     row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID = conc_subset$window,
                              Means = rowMeans(conc_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "conc", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))


high_pmol <- subset(high_means, high_means$conc > 20)#39 windows
sep_info <- conc[, c("chr", "start", "end", "window")]
high_pmol <- merge(high_pmol, sep_info, by = "window")

###annotate genomic position to the geneID
#### make bed file and will use bedtools to overlay the gene name
high_fmol_regions <- high_pmol[, c("chr", "start", "end", "conc")]
write.table(high_fmol_regions, sep = "\t", quote = FALSE,
            row.names = FALSE, file = "2020_high_fmol_regions_mapp_HCT116sesame_ctrlproj.bed")

##To type into command line:
###  bedtools intersect -wao
###-a 2020_high_fmol_regions_HCT116_ctrlproj.bed
###-b gencode.v33.basic.annotation.bed
###> 2020_annotated_high_fmol.bed

gene_names <- read.table("2020_annotated_high_fmol_mapp.bed", header = TRUE)
colnames(gene_names) <- c("chr", "start", "end", "conc",
                           "gene_chr", "gene_start", "gene_end",
                           "geneID", "overlap_bp")

ensemble_ids <- as.character(gene_names$geneID)
ensemble_ids <- gsub("\\..*", "", ensemble_ids)

symbols <- mapIds(org.Hs.eg.db, keys = ensemble_ids,
                  keytype = "ENSEMBL", column = "SYMBOL")
symbols <- replace(symbols, is.na(symbols), "not _annotated")
symbols <- replace(symbols, symbols == "NULL", "not _annotated")

high_pmol_genes_sorted <- high_pmol_genes[order(high_pmol_genes$chr,
                                                high_pmol_genes$start,
                                                high_pmol_genes$end,
                                                high_pmol_genes$geneID), ]


png("2020_high_fmol_anno_table.png",
    height = 6, width = 25, unit = "in", res = 300)
p <- tableGrob(table)
grid.arrange(p)
dev.off()

#### Correlation of read counts to M-values
data_readcounts <- read.delim("~/Projects/2018_PTB/data/2019_ControlAlignment/2020_Human0.01_hg38_intersect_readcounts.csv",
                              sep = "\t", header = FALSE)
colnames(data_readcounts) <- c("frag_chr", "frag_start", "frag_end",
                               "sample1", "sample2",
                               "chr", "start", "end", "overlap")
data_readcounts <- data_readcounts[, c("chr", "start", "end",
                                       "sample1", "sample2")]

data_readcounts_agg <- ddply(data_readcounts,
                             .(data_readcounts$chr, data_readcounts$start, data_readcounts$end),
                             numcolwise(sum), .progress = "text")
colnames(data_readcounts_agg) <- c("chr", "start", "end", "sample1", "sample2")

##remove repetitive and blacklist regions from the data
data_readcounts_agg$window <- paste0(data_readcounts_agg$chr, "_",
                                     data_readcounts_agg$start, "_", data_readcounts_agg$end)
rownames(data_readcounts_agg) <- data_readcounts_agg$window

#subset to windows in the pmol analysis
data_readcounts_agg <- data_readcounts_agg[rownames(data_readcounts_agg) %in% rownames(conc), ]

##subset readcounts and epic to only those windows present in both
windows_of_interest <- epic$window

rownames(data_readcounts_agg) <- data_readcounts_agg$window
readcounts_subset <- subset(data_readcounts_agg, rownames(data_readcounts_agg) %in% windows_of_interest)

rownames(epic) <- epic$window
epic_subset <- subset(epic, rownames(epic) %in% rownames(readcounts_subset))

##mean the samples
rc_means <- data.frame(ID = readcounts_subset$window,
                       Means = rowMeans(readcounts_subset[, c("sample1", "sample2")]))
epic_means <- data.frame(ID = epic_subset$window,
                         Means = rowMeans(epic_subset[, c("Sample1", "Sample2", "Sample3")]))

means <- merge(rc_means, epic_means, by = "ID")
colnames(means) <- c("window", "readcount", "epic")

##changes betas to m-values
means$epic_m <- log2(means$epic / (1 - means$epic))
print(head(means))

###Change corr_plot function so that y is readcounts
corr_plot <- function(data) {
ggplot(data, aes(x = epic_m, y = readcounts)) +
  geom_point(color = "#2980B9", size = 2) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, color = "#2C3E50") +
  ylab("Reads") +
  xlab("DNA methylation (M-value)") +
  scale_y_continuous(labels = comma)
}


##Correlate bins with >3 CpG represented on array
##Let's correlate those this more than 3 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg3.csv",
                     row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

##Correlate bins with >3 CpG represented on array
##Let's correlate those this more than 3 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg3.csv",
                     row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(data_readcounts_agg) <- data_readcounts_agg$window
readcounts_subset <- subset(data_readcounts_agg, rownames(data_readcounts_agg) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(readcounts_subset))

##mean the samples
rc_high_means <- data.frame(ID = readcounts_subset$window,
                            Means = rowMeans(readcounts_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(rc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "readcounts", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))

##Correlation calculation
cor.test(high_means$readcounts, high_means$epic_m, method = "spearman")
##0.617944

##Full plot
png("2020_readcounts_epicsesame_corr_mapp_morethan3CpG.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more3cpg_rc <- corr_plot(high_means) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
        scale_y_continuous(limits = c(0, 2500),
                           breaks = c(0, 500, 1000, 1500, 2000, 2500))
corr_more3cpg
dev.off()

##Plot- cut to see the small correlation
png("2020_readcounts_epicsesame_corr_mapp_morethan3CpG_cut.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more3cpg_cut <- corr_plot(high_means) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(labels = comma, limits = c(0, 100))
corr_more3cpg_cut
dev.off()

##Let's correlate those this more than 5 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg5.csv",
                     row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(data_readcounts_agg) <- data_readcounts_agg$window
readcounts_subset <- subset(data_readcounts_agg, rownames(data_readcounts_agg) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(readcounts_subset))

##mean the samples
rc_high_means <- data.frame(ID = readcounts_subset$window,
                            Means = rowMeans(readcounts_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(rc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "readcounts", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))

##Correlation calculation
cor.test(high_means$readcounts, high_means$epic_m, method = "spearman")
##0.8017081

##Full plot
png("2020_readcounts_epicsesame_corr_mapp_morethan5CpG.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more5cpg_rc <- corr_plot(high_means) +
        scale_x_continuous(limits = c(-8, 8),
                           breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
        scale_y_continuous(limits = c(0, 2500),
                           breaks = c(0, 500, 1000, 1500, 2000, 2500))
corr_more5cpg
dev.off()

##Plot- cut to see the small correlation
png("2020_readcounts_epicsesame_corr_mapp_morethan5CpG_cut.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more5cpg_cut <- corr_plot(high_means) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(labels = comma, limits = c(0, 100)) +
corr_more5cpg_cut
dev.off()

###Correlate read counts with >7cpg
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg7.csv",
                     row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(data_readcounts_agg) <- data_readcounts_agg$window
readcounts_subset <- subset(data_readcounts_agg, rownames(data_readcounts_agg) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(readcounts_subset))

##mean the samples
rc_high_means <- data.frame(ID = readcounts_subset$window,
                            Means = rowMeans(readcounts_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(rc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "readcounts", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))

##Correlation calculation
cor.test(high_means$readcounts, high_means$epic_m, method = "spearman")
##0.8279825

##Full plot
png("2020_readcounts_epicsesame_corr_mapp_morethan7CpG.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more7cpg_rc <- corr_plot(high_means) +
        scale_x_continuous(limits = c(-8, 8),
                           breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
        scale_y_continuous(limits = c(0, 2500),
                           breaks = c(0, 500, 1000, 1500, 2000, 2500))
corr_more7cpg
dev.off()

##Plot- cut to see the small correlation
png("2020_readcounts_epicsesame_corr_mapp_morethan7CpG_cut.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more7cpg_cut <- corr_plot(high_means) +
  scale_x_continuous(limits = c(-8, 8),
                     breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(labels = comma, limits = c(0, 100))
corr_more7cpg_cut
dev.off()

###>10 cpgs
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg10.csv",
                     row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window", "Sample1", "Sample2", "Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(data_readcounts_agg) <- data_readcounts_agg$window
readcounts_subset <- subset(data_readcounts_agg, rownames(data_readcounts_agg) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(readcounts_subset))

##mean the samples
rc_high_means <- data.frame(ID = readcounts_subset$window,
                            Means = rowMeans(readcounts_subset[, c("sample1", "sample2")]))
epic_high_means <- data.frame(ID = high_subset$window,
                              Means = rowMeans(high_subset[, c("Sample1", "Sample2", "Sample3")]))

high_means <- merge(rc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window", "readcounts", "epic")
high_means$epic_m <- log2(high_means$epic / (1 - high_means$epic))

##Correlation calculation
cor.test(high_means$readcounts, high_means$epic_m, method = "spearman")
##0.8587268

##Full plot
png("2020_readcounts_epicsesame_corr_mapp_morethan10CpG.png",
    height = 6, width = 6, units = "in", res = 300)
corr_more10cpg_rc <- corr_plot(high_means) +
        scale_x_continuous(limits = c(-8, 8),
                           breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
        scale_y_continuous(limits = c(0,2500),
                           breaks = c(0, 500, 1000, 1500, 2000, 2500))
dev.off()

###Plot all plots together
png("2020_readcounts_sesamemval_corr_mapp_allplts.png",
    height = 12, width = 8, units = "in", res = 300)
all_rc_corr_plot <- grid.arrange(arrangeGrob(corr_more3cpg + theme(axis.title.x = element_blank()),
                                             top = "More than 3 CpG"),
                            arrangeGrob(corr_more3cpg_cut + theme(axis.title.x = element_blank(),
                                                                  axis.title.y = element_blank()),
                                        top = "More than 3 CpG: zoom"),
                            arrangeGrob(corr_more5cpg + theme(axis.title.x = element_blank()),
                                        top = "More than 5 CpG"),
                            arrangeGrob(corr_more5cpg_cut + theme(axis.title.x = element_blank(),
                                                                  axis.title.y = element_blank()),
                                        top = "More than 5 CpG: zoom"),
                            arrangeGrob(corr_more7cpg + theme(axis.title.x = element_blank()),
                                        top = "More than 7 CpG"),
                            arrangeGrob(corr_more10cpg_cut + theme(axis.title.x = element_blank(),
                                                                   axis.title.y = element_blank()),
                                        top = "More than 10 CpG: zoom"), ncol = 2)

annotate_figure(all_rc_corr_plot,
               bottom = text_grob("DNA methylation (M-value)", color = "black",
                                   size = 16))
dev.off()

##Read count and fmol plots together with no zoom
png("2020_fmol_readcounts_sesamemval_corr_mapp_allplts.png",
    height = 12, width = 8, units = "in", res = 300)
all_rc_corr_plot <- grid.arrange(arrangeGrob(corr_more3cpg + scale_y_continuous(limits = c(0, 0.5)) +
                                               theme(axis.title.x = element_blank()),
                                             top = grid::textGrob("≥ 3 CpG/300 bp window (N = 50,053 windows)",
                                                                  x = 0.15, hjust = 0)),
                                 arrangeGrob(corr_more3cpg_rc + theme(axis.title.x = element_blank()),
                                             top = ""),
                            arrangeGrob(corr_more5cpg + scale_y_continuous(limits = c(0, 0.5)) +
                                          theme(axis.title.x = element_blank()),
                                        top = grid::textGrob("≥ 5 CpG/300 bp window (N = 12,602 windows)",
                                                             x = 0.15, hjust = 0)),
                            arrangeGrob(corr_more5cpg_rc + theme(axis.title.x = element_blank()), top = ""),
                            arrangeGrob(corr_more7cpg + scale_y_continuous(limits = c(0, 0.5)) +
                                          theme(axis.title.x = element_blank()),
                                        top = grid::textGrob("≥ 7 CpG/300 bp window (N = 3,067 windows)",
                                                             x = 0.15, hjust = 0)),
                            arrangeGrob(corr_more7cpg_rc + theme(axis.title.x = element_blank()), top = ""),
                            arrangeGrob(corr_more10cpg + scale_y_continuous(limits = c(0, 0.5)) +
                                          theme(axis.title.x = element_blank()),
                                        top = grid::textGrob("≥ 10 CpG/300 bp window (N = 294 windows)",
                                                             x = 0.15, hjust = 0)),
                            arrangeGrob(corr_more10cpg_rc + theme(axis.title.x = element_blank()),
                                        top = ""), ncol = 2)

annotate_figure(all_rc_corr_plot,
               bottom = text_grob("DNA methylation (M-value)", color = "black",
                                   size = 16))
dev.off()