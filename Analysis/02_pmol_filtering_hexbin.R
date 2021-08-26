#)!/usr/bin/env Rscript

##load libraries
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
library(BlandAltmanLeh)
library(cowplot)
library(GGally)
library(hexbin)
library(RColorBrewer)
library(Hmisc)

##load data
#data <- read.table("2020_intersect_fmol_mapability.bed")
#colnames(data) <- c("chr_umap","start_umap","end_umap", "map_score", "chr", "start", "end", "overlap")

##aggregate by window - mean mappability score
#data$window <- as.factor(paste0(data$chr,"_", data$start, "_", data$end))
#data_agg <- ddply(data,.(data$window, data$chr, data$start, data$end),numcolwise(mean), .progress = "text")
#colnames(data_agg) <- c("window","chr", "start", "end", "start_umap", "end_umap", "map_score", "overlap")
#print(head(data_agg))
#print(dim(data_agg))
#write.table(data_agg, header = TRUE, file = "2020_mappability_300bp_windows.csv")

##aggregate by window - min mapability score
#data$window <- as.factor(paste0(data$chr,"_", data$start, "_", data$end))
#data_agg <- ddply(data,.(data$window, data$chr, data$start, data$end),numcolwise(min), .progress = "text")
#colnames(data_agg) <- c("window","chr", "start", "end", "start_umap", "end_umap", "map_score", "overlap")
#print(head(data_agg))
#print(dim(data_agg))
#write.table(data_agg, header = TRUE, file = "2020_min_mappability_300bp_windows.csv")

mapp_windows <- read.csv("~/Annotations/2020_min_mappability_300bp_windows.csv", header = TRUE)
mapp_windows$row.names <- NULL
mapp_windows$X <- NULL

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

print(head(conc))
print(dim(conc))

##merge with mappability data
data_all <- merge(conc, mapp_windows, by = "window")
data_all[is.na(data_all)] <- 0 #replace NA with 0 (this could be region was not methylated or was not picked up)

##take mean between samples
data_all$sample_mean <- (data_all$sample1_pmol + data_all$sample2_pmol) / 2
print(head(data_all))
print(dim(data_all))
#5103552

##ggplot2 theme 
source("/cluster/home/wilsons/commonly_used_code/ggplot2_theme_bw.R")

##plot function
corr_plot <- function(means){
ggplot(means, aes(x=map_score, y= sample_mean)) +
  geom_point(color='#2980B9', size = 2) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50') +
  ylab("Picomole") +
  xlab("Mappability score") +
  stat_cor(method = "spearman") +
  scale_y_continuous(labels = comma)
}

##Make plot
png("2020_conc_mapp_corr.png", height = 6, width = 6, units = "in", res =300)
corr_plot(data_all)
dev.off()
#cor -0.028

###Plots with mean mappability score
##Make ggplot hexbin plot
png("2020_conc_mapp_corr_hexbin_morecounts.png", height = 6, width = 6, units = "in", res =300)
hex <- hexbin(data_all$map_score, data_all$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count) ) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Mappability score", y = "Picomoles", fill = "Log of counts") + 
  scale_y_continuous(labels = comma)
dev.off()

##Looking into what the high fmol regions are
### subset the windows that have > 5000fmol
high_pmol <- subset(data_all, data_all$sample_mean > 1) 
rownames(high_pmol) <- high_pmol$window
#### 3 windows above 2pmol, 5 over 1 pmol

### Are any of these regions on the ENCODE blacklist?
repeats <- read.table("~/Annotations/2020_repeatregions_mappedto300bpwindows.bed", sep = "\t", header = FALSE)
###subset to regions where there are overlaps between our windows
repeats <- subset(repeats, repeats$V9>0)
repeats <- repeats[,c(1:3)]
colnames(repeats) <- c("chr","start","end")
repeats$window <- paste0(repeats$chr,"_",repeats$start,"_",repeats$end)
repeats <- unique(repeats)
rownames(repeats) <-repeats$window

blacklist <- read.table("~/Annotations/2020_new_ENCODE_blacklist.csv", sep = "\t", header = FALSE)
colnames(blacklist) <- c("bl_chr", "bl_start", "bl_end", "chr", "start", "end", "overlap")
blacklist$window <- paste0(blacklist$chr,"_",blacklist$start,"_",blacklist$end)
blacklist <- unique(blacklist)

overlap_repetitiveregions <- high_pmol[rownames(high_pmol) %in% rownames(repeats),]
##Every single one is on the UCSC repetitive regions black list

overlap_blacklist <- high_pmol[rownames(high_pmol) %in% rownames(blacklist),]

##remove the blacklist regions and remake the hexbin plot.
rm <- repeats$window
rownames(data_all) <- data_all$window
data_noblacklist <- data_all[!rownames(data_all) %in% rm,]

rm <- blacklist$window
data_noblacklist <- data_noblacklist[!rownames(data_noblacklist) %in% rm,]
#4546731 windows

png("2020_conc_mapp_corr_hexbin_noblacklist_morecounts.png", height = 6, width = 6, units = "in", res =300)
hex <- hexbin(data_noblacklist$map_score, data_noblacklist$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count) ) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Mappability score", y = "Picomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

## Follow up on the high fmol sites not in the blacklist or  STR list
high_pmol2 <- subset(data_noblacklist, data_noblacklist$sample_mean > 0.5)
write.csv(high_pmol2, file = "2020_highpmol_noblacklist.csv")
#117 windows

###Plots with minimum mappability score
##Make ggplot hexbin plot
png("2020_conc_minmapp_corr_hexbin_morecounts.png", height = 6, width = 6, units = "in", res =300)
hex <- hexbin(data_all$map_score, data_all$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count) ) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Mappability score", y = "Femtomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

##Looking into what the high fmol regions are
### subset the windows that have > 5000fmol

high_pmol <- subset(data_all, data_all$sample_mean > 1)
rownames(high_pmol) <- high_pmol$window
#### 5 windows, 3 on chromosome 16

##remove the blacklist regions and remake the hexbin plot.
rm <- repeats$window
rownames(data_all) <- data_all$window
data_noblacklist <- data_all[!rownames(data_all) %in% rm,]

rm <- blacklist$window
data_noblacklist <- data_noblacklist[!rownames(data_noblacklist) %in% rm,]

png("2020_conc_minmapp_corr_hexbin_noblacklist_morecounts.png", height = 6, width = 6, units = "in", res =300)
hex <- hexbin(data_noblacklist$map_score, data_noblacklist$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count) ) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Mappability score", y = "Picomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

## Follow up on the high fmol sites not in the blacklist or  STR list
high_pmol2 <- subset(data_noblacklist, data_noblacklist$sample_mean > 0.5)
write.csv(high_pmol2, file = "2020_highpmol_minmappability_noblacklist.csv")
# 117 windows, most have high map scores

##Hexbin plot with standard deviation
data_noblacklist$sd <- apply(data_noblacklist[, c("sample1_pmol","sample2_pmol")],1,sd)

png("2020_pmol_sd_corr.png", height = 6, width = 6, units = "in", res =300)
hex <- hexbin(data_noblacklist$sd, data_noblacklist$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count) ) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Standard deviation", y = "Picomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

rm <- subset(data_noblacklist, data_noblacklist$sd > 0.05 | data_noblacklist$map_score < 0.5)
data_filtered <- data_noblacklist[!rownames(data_noblacklist) %in% rownames(rm),]
#4446375 windows 

##Set up colour palette
library(RColorBrewer)
pal <- brewer.pal(n = 8, name = "Blues")[3:8]

png("2020_pmol_sd_filtered_pubplt.png", height = 6, width = 6, units = "in", res =300)
hex <- hexbin(data_filtered$sd, data_filtered$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
pmol_sd <- ggplot(gghex, aes(x = x, y = y, fill = count)) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_gradientn(colors = pal, breaks = c(0, 5, 10, 15,20), labels = c("0" , "5", "10", "15", "20")) +
  labs(x = "Standard deviation", y = "Picomoles", fill = "Reads") +
  scale_y_continuous(labels = comma)
pmol_sd
dev.off()


##Redo hebin mappability without the sd
png("2020_pmol_mapp_filtered_pubplt.png", height = 6, width = 6, units = "in", res =300)
hex <- hexbin(data_filtered$map_score, data_filtered$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
pmol_mapp <- ggplot(gghex, aes(x = x, y = y, fill = count) ) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_gradientn(colors = pal,  breaks = c(0, 5, 10, 15,20), labels = c("0" , "5", "10", "15", "20")) +
  labs(x = "Mappability score", y = "Picomoles", fill = "Reads") +
  scale_y_continuous(labels = comma)
pmol_mapp
dev.off()

##Extracting legend, so all plots are same size in the grid.arrange()
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a_gplot) {
  tmp <- ggplot_gtable(ggplot_build(a_gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
  }

##Arranging plots together
png("2020_pmol_sd_mapp_corr_pubplt.png", height =6, width = 12, units = "in", res =300)
mylegend <- g_legend(pmol_sd) 
grid.arrange(arrangeGrob(pmol_sd + theme(legend.position = "none")),
		arrangeGrob(pmol_mapp + theme(legend.position = "none")), nrow = 1, mylegend)
dev.off()

##What are the highest pmol predictions
high_pmol <- data_filtered[data_filtered$sample_mean > 0.25 & data_filtered$map_score < 1 | data_filtered$sample_mean > 0.7,] #11 windows
write.csv(high_pmol, file = "2020_highpmol_minmappability_noblacklist.csv")
