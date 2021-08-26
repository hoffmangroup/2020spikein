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
#library(lumi)
library(scales)
library(ggpubr)

#Load data

##Predicted concentration metrics
##load in pmol data
sample1 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/S6547_all_pmol_hg38_intersect_adjbinned.bed", sep = "\t", header = TRUE)
colnames(sample1) <- c("chr", "start", "end", "sample1_read_count", "sample1_pmol")
sample1$window <- as.factor(paste0(sample1$chr,"_", sample1$start, "_", sample1$end))
sample1 <- sample1[-c(1:3)]

sample2 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/S6548_all_pmol_hg38_intersect_adjbinned.bed", sep = "\t", header = TRUE)
colnames(sample2) <- c("chr", "start", "end", "sample2_read_count", "sample2_pmol")
sample2$window <- as.factor(paste0(sample2$chr,"_", sample2$start, "_", sample2$end))
sample2 <- sample2[-c(1:3)]

#merge all samples
conc <- merge(sample1, sample2, by = "window", all = TRUE, fill = TRUE)
rownames(conc) <- conc$window

print(head(conc))
print(dim(conc))

epic <- read.delim("2019_HCT116_sesame_betas_windows.csv", sep = ",", header = TRUE)
epic$X <- NULL
epic$start <- NULL
epic$end <- NULL
epic$probepos <- NULL
colnames(epic) <- c("window","Sample1","Sample2","Sample3")

##Remove repetitive elements- simple repeats only 
repeats <- read.table("~/Annotations/2020_repeatregions_mappedto300bpwindows.bed", sep = "\t", header = FALSE)
###subset to regions where there are overlaps between our windows
repeats <- subset(repeats, repeats$V9>0)
repeats <- repeats[,c(1:3)]
colnames(repeats) <- c("chr","start","end")
repeats$window <- paste0(repeats$chr,"_",repeats$start,"_",repeats$end)
repeats <- unique(repeats)
rownames(repeats) <-repeats$window

##remove repetitive regions from the data
conc <- conc[!rownames(conc) %in% rownames(repeats),]

##remove blacklist regions
blacklist <- read.table("~/Annotations/2020_new_ENCODE_blacklist.csv", sep = "\t", header = FALSE)
colnames(blacklist) <- c("bl_chr", "bl_start", "bl_end", "chr", "start", "end", "overlap")
blacklist$window <- paste0(blacklist$chr,"_",blacklist$start,"_",blacklist$end)
blacklist <- unique(blacklist)

##remove blacklist regions
conc <-conc[!rownames(conc) %in% blacklist$window,]

## remove low mappability regions (mappability score<0.5)
mappability <- read.csv("~/Annotations/2020_min_mappability_300bp_windows.csv", header = TRUE)
mappability$X <- NULL

mappability_0.5 <- subset(mappability, mappability$map_score<0.5)

conc <- conc[!rownames(conc) %in% mappability_0.5$window,] 
#4456768 windows

#subset out samples with sd >0.05
conc$sd <- apply(conc[, c("sample1_pmol","sample2_pmol")],1,sd)

rm <- subset(conc, conc$sd > 0.05)
conc <- conc[!rownames(conc) %in% rownames(rm),]
#4455866

##subset conc and epic to only those windows present in both
windows_of_interest <- epic$window

conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(epic) <- epic$window
epic_subset <- subset(epic, rownames(epic) %in% rownames(conc))

##mean the samples
conc_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_pmol","sample2_pmol")]))
epic_means <- data.frame(ID=epic_subset$window, Means=rowMeans(epic_subset[,c("Sample1","Sample2","Sample3")]))

means <- merge(conc_means, epic_means, by = "ID")
colnames(means) <- c("window","conc","epic")

##changes betas to m-values
means$epic_m <- log2(means$epic/(1-means$epic)) 

##Correlation calculation
cor.test(means$conc, means$epic_m, method = 'spearman')
##0.1475338 

conc[is.na(conc)] <- 0

##ggplot2 theme 
source("/cluster/home/wilsons/commonly_used_code/ggplot2_theme_bw.R")

##plot function
corr_plot <- function(means){
ggplot(means, aes(x= epic_m, y= conc)) +
  geom_point(color='#2980B9', size = 2) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50') +
  ylab("Picomoles") +
  xlab("DNA methylation (M-value)") +
  #stat_cor(method = "spearman", aes(label = paste(..r.label.., expression(paste("p < 2.2x",10^-16)), sep = "*`,`~")), digits = 2, r.digits = 2, r.accuracy = 0.01)+
  scale_y_continuous(labels = comma)
}

##Plot
png("2020_conc_epic_sesame_mapp_corr.png", height = 6, width = 6, units = "in", res =300)
all_pts <- corr_plot(means)
all_pts
dev.off()


##Correlate bins with >3 CpG represented on array
##Let's correlate those this more than 3 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg3.csv", row.names = 1, header = TRUE)
#50053 windows

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window","Sample1","Sample2","Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_pmol","sample2_pmol")]))
epic_high_means <- data.frame(ID=high_subset$window, Means=rowMeans(high_subset[,c("Sample1","Sample2","Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window","conc","epic")
high_means$epic_m <- log2(high_means$epic/(1-high_means$epic))

##Correlation calculation
cor.test(high_means$conc, high_means$epic_m, method = 'spearman')
##0.6239577  - with 3 CpG and 50053 windows

##Full plot
png("2020_conc_epicsesame_corr_mapp_morethan3CpG.png", height = 6, width = 6, units = "in", res =300)
corr_more3cpg <- corr_plot(high_means) + 
                 scale_x_continuous(limits = c(-8, 8), breaks = c(-6, -4, -2, 0, 2, 4, 6)) + 
                 scale_y_continuous(limits = c(0,10), breaks = c(0, 2, 4, 6, 8, 10))
corr_more3cpg 
dev.off()

##Let's correlate those this more than 5 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg5.csv", row.names = 1, header = TRUE)
#12602 windows

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window","Sample1","Sample2","Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_pmol","sample2_pmol")]))
epic_high_means <- data.frame(ID=high_subset$window, Means=rowMeans(high_subset[,c("Sample1","Sample2","Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window","conc","epic")
high_means$epic_m <- log2(high_means$epic/(1-high_means$epic))

##Correlation calculation
cor.test(high_means$conc, high_means$epic_m, method = 'spearman')
##0.8007208    - with more than 5 CpG and 12602 windows

##Full plot
png("2020_conc_epicsesame_corr_mapp_morethan5CpG.png", height = 6, width = 6, units = "in", res =300)
corr_more5cpg <- corr_plot(high_means) +
                 scale_x_continuous(limits = c(-8, 8), breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
scale_y_continuous(limits = c(0,10), breaks = c(0, 2, 4, 6, 8, 10))
corr_more5cpg
dev.off()

### 7 CpG
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg7.csv", row.names = 1, header = TRUE)
#3567 windows

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window","Sample1","Sample2","Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_pmol","sample2_pmol")]))
epic_high_means <- data.frame(ID=high_subset$window, Means=rowMeans(high_subset[,c("Sample1","Sample2","Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window","conc","epic")
high_means$epic_m <- log2(high_means$epic/(1-high_means$epic))

##Correlation calculation
cor.test(high_means$conc, high_means$epic_m, method = 'spearman')
##0.8206853   - with more than 7 CpG and 3567 windows

##Full plot
png("2020_conc_epicsesame_corr_mapp_morethan7CpG.png", height = 6, width = 6, units = "in", res =300)
corr_more7cpg <- corr_plot(high_means) +
  scale_x_continuous(limits = c(-8, 8), breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
scale_y_continuous(limits = c(0,10), breaks = c(0, 2, 4, 6, 8, 10))
corr_more7cpg
dev.off()

###10 CpG
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg10.csv", row.names = 1, header = TRUE)
#294 windows

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window","Sample1","Sample2","Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

rownames(conc) <- conc$window
conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc))

##mean the samples
conc_high_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_pmol","sample2_pmol")]))
epic_high_means <- data.frame(ID=high_subset$window, Means=rowMeans(high_subset[,c("Sample1","Sample2","Sample3")]))

high_means <- merge(conc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window","conc","epic")
high_means$epic_m <- log2(high_means$epic/(1-high_means$epic))

##Correlation calculation
cor.test(high_means$conc, high_means$epic_m, method = 'spearman')
##0.8703139    - with more than 10 CpG and 294 windows

##Full plot
png("2020_conc_epicsesame_corr_mapp_morethan10CpG.png", height = 6, width = 6, units = "in", res =300)
corr_more10cpg <- corr_plot(high_means) +
  scale_x_continuous(limits = c(-8, 8), breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
scale_y_continuous(limits = c(0,10), breaks = c(0, 2, 4, 6, 8, 10))
corr_more10cpg
dev.off()

#### Correlation of read counts to M-values
##mean the samples
rc_means <- data.frame(ID=conc$window, Means=rowMeans(conc[,c("sample1_read_count","sample2_read_count")]))
epic_means <- data.frame(ID=epic_subset$window, Means=rowMeans(epic_subset[,c("Sample1","Sample2","Sample3")]))

means <- merge(rc_means, epic_means, by = "ID")
colnames(means) <- c("window","readcount","epic")

##changes betas to m-values
means$epic_m <- log2(means$epic/(1-means$epic))
print(head(means))

###Change corr_plot function so that y is readcounts
corr_plot <- function(data){
ggplot(data, aes(x= epic_m, y= readcounts)) +
  geom_point(color='#2980B9', size = 2) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50') +
  ylab("Reads") +
  xlab("DNA methylation (M-value)") +
  #stat_cor(method = "spearman", aes(label = paste(..r.label.., expression(paste("p < 2.2x",10^-16)), sep = "*`,`~")), digits = 2, r.digits = 2) +
  scale_y_continuous(labels = comma)
}


##Correlate bins with >3 CpG represented on array
##Let's correlate those this more than 3 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg3.csv", row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window","Sample1","Sample2","Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc_subset))

##mean the samples
rc_high_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_read_count","sample2_read_count")]))
epic_high_means <- data.frame(ID=high_subset$window, Means=rowMeans(high_subset[,c("Sample1","Sample2","Sample3")]))

high_means <- merge(rc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window","readcounts","epic")
high_means$epic_m <- log2(high_means$epic/(1-high_means$epic))

##Correlation calculation
cor.test(high_means$readcounts, high_means$epic_m, method = 'spearman')
##0.6203779  

##Full plot
png("2020_readcounts_epicsesame_corr_mapp_morethan3CpG.png", height = 6, width = 6, units = "in", res =300)
corr_more3cpg_rc <- corr_plot(high_means) +
	scale_x_continuous(limits = c(-8, 8), breaks = c(-6, -4, -2, 0, 2, 4, 6)) + 
        scale_y_continuous(limits = c(0,500), breaks = c(0, 100, 200, 300, 400, 500))
corr_more3cpg_rc
dev.off()

##Let's correlate those this more than 5 cpg present in window
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg5.csv", row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window","Sample1","Sample2","Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc_subset))

##mean the samples
rc_high_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_read_count","sample2_read_count")]))
epic_high_means <- data.frame(ID=high_subset$window, Means=rowMeans(high_subset[,c("Sample1","Sample2","Sample3")]))

high_means <- merge(rc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window","readcounts","epic")
high_means$epic_m <- log2(high_means$epic/(1-high_means$epic))

##Correlation calculation
cor.test(high_means$readcounts, high_means$epic_m, method = 'spearman')
## 0.8138254

##Full plot
png("2020_readcounts_epicsesame_corr_mapp_morethan5CpG.png", height = 6, width = 6, units = "in", res =300)
corr_more5cpg_rc <- corr_plot(high_means) +
        scale_x_continuous(limits = c(-8, 8), breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
        scale_y_continuous(limits = c(0,500), breaks = c(0, 100, 200, 300, 400, 500))
corr_more5cpg_rc
dev.off()

###Correlate read counts with >7cpg
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg7.csv", row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window","Sample1","Sample2","Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc_subset))

##mean the samples
rc_high_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_read_count","sample2_read_count")]))
epic_high_means <- data.frame(ID=high_subset$window, Means=rowMeans(high_subset[,c("Sample1","Sample2","Sample3")]))

high_means <- merge(rc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window","readcounts","epic")
high_means$epic_m <- log2(high_means$epic/(1-high_means$epic))

##Correlation calculation
cor.test(high_means$readcounts, high_means$epic_m, method = 'spearman')
##0.8461799

##Full plot
png("2020_readcounts_epicsesame_corr_mapp_morethan7CpG.png", height = 6, width = 6, units = "in", res =300)
corr_more7cpg_rc <- corr_plot(high_means) +
        scale_x_continuous(limits = c(-8, 8), breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
scale_y_continuous(limits = c(0,500), breaks = c(0, 100, 200, 300, 400, 500))
corr_more7cpg_rc
dev.off()

###>10 cpgs
high_cpg <- read.csv("2020_HCT116_sesame_betas_windows_highcpg10.csv", row.names = 1, header = TRUE)

high_cpg$CpG_beg <- NULL
high_cpg$CpG_end <- NULL
colnames(high_cpg) <- c("window","Sample1","Sample2","Sample3")

##subset conc and epic to only those windows present in both
windows_of_interest <- high_cpg$window

conc_subset <- subset(conc, rownames(conc) %in% windows_of_interest)

rownames(high_cpg) <- high_cpg$window
high_subset <- subset(high_cpg, rownames(high_cpg) %in% rownames(conc_subset))

##mean the samples
rc_high_means <- data.frame(ID=conc_subset$window, Means=rowMeans(conc_subset[,c("sample1_read_count","sample2_read_count")]))
epic_high_means <- data.frame(ID=high_subset$window, Means=rowMeans(high_subset[,c("Sample1","Sample2","Sample3")]))

high_means <- merge(rc_high_means, epic_high_means, by = "ID")
colnames(high_means) <- c("window","readcounts","epic")
high_means$epic_m <- log2(high_means$epic/(1-high_means$epic))

##Correlation calculation
cor.test(high_means$readcounts, high_means$epic_m, method = 'spearman')
##0.877321 

##Full plot
png("2020_readcounts_epicsesame_corr_mapp_morethan10CpG.png", height = 6, width = 6, units = "in", res =300)
corr_more10cpg_rc <- corr_plot(high_means) +
        scale_x_continuous(limits = c(-8, 8), breaks = c(-6, -4, -2, 0, 2, 4, 6)) +
scale_y_continuous(limits = c(0,500), breaks = c(0, 100, 200, 300, 400, 500))
corr_more10cpg_rc
dev.off()

##Read count and pmol
##y axis cut to 500 in read counts, this does cut some of the points off, but looks cleaners. The correlations reported are for all points. 
png("2020_pmol_readcounts_sesamemval_corr_mapp_allplts.png", height = 12, width = 8, units = "in", res = 300)
all_rc_corr_plot <- grid.arrange(arrangeGrob(corr_more3cpg + scale_y_continuous(limits=c(0,0.5)) + theme(axis.title.x = element_blank()), top = grid::textGrob("≥ 3 CpG/300 bp window (N = 34,735 windows)", x = 0.6, hjust =0)),
			    arrangeGrob(corr_more3cpg_rc + scale_y_continuous(limits=c(0,500)) + theme(axis.title.x = element_blank()), top = ""),
                            arrangeGrob(corr_more5cpg + scale_y_continuous(limits=c(0,0.5)) + theme(axis.title.x = element_blank()), top = grid::textGrob("≥ 5 CpG/300 bp window (N = 7,203 windows)", x = 0.6, hjust = 0)),
			    arrangeGrob(corr_more5cpg_rc + scale_y_continuous(limits=c(0,500)) + theme(axis.title.x = element_blank()), top = ""),
                            arrangeGrob(corr_more7cpg + scale_y_continuous(limits=c(0,0.5)) + theme(axis.title.x = element_blank()), top = grid::textGrob("≥ 7 CpG/300 bp window (N = 1,898 windows)", x = 0.6, hjust = 0)),
			    arrangeGrob(corr_more7cpg_rc + scale_y_continuous(limits=c(0,500)) + theme(axis.title.x = element_blank()), top =""),
			    arrangeGrob(corr_more10cpg + scale_y_continuous(limits=c(0,0.5)) + theme(axis.title.x = element_blank()), top = grid::textGrob("≥ 10 CpG/300 bp window (N = 147 windows)", x =0.6, hjust = 0)),
                            arrangeGrob(corr_more10cpg_rc + theme(axis.title.x = element_blank()), top = ""), ncol =2)

annotate_figure(all_rc_corr_plot,
               bottom = text_grob("DNA methylation (M-value)", color = "black",
                                   size = 16))
dev.off()


