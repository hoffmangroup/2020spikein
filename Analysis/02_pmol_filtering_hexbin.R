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
data <- read.table("2020_intersect_fmol_mapability.bed")
colnames(data) <- c("chr_umap", "start_umap", "end_umap",
                    "map_score", "chr", "start", "end", "overlap")

##aggregate by window - min mapability score
data$window <- as.factor(paste0(data$chr, "_", data$start, "_", data$end))
data_agg <- ddply(data, .(data$window, data$chr,
                          data$start, data$end), numcolwise(min),
                  .progress = "text")
colnames(data_agg) <- c("window", "chr", "start", "end",
                        "start_umap", "end_umap", "map_score",
                        "overlap")
print(head(data_agg))
print(dim(data_agg))
write.table(data_agg, header = TRUE, file = "2020_min_mappability_300bp_windows.csv")

data_agg$row.names <- NULL
data_agg$X <- NULL

##load in fmol data
conc <- read.delim("~/Projects/2018_PTB/data/2019_ControlAlignment/2020_human0.01_adj_binned.csv",
                   sep = ",", header = TRUE)
conc$X <- NULL
colnames(conc) <- c("chr", "start", "end", "sample1", "sample2")
conc$window <- paste0(conc$chr, "_", conc$start, "_", conc$end)
print(head(conc))
print(dim(conc))

##merge with mappability data
data_all <- merge(conc, data_agg, by = "window")
data_all <- data_all[, c(1:6, 10:13)]
colnames(data_all) <- c("window", "chr", "start", "end",
                        "sample1", "sample2",
                        "start_umap", "end_umap",
                        "map_score", "overlap")

##take mean between samples
data_all$sample_mean <- (data_all$sample1 + data_all$sample2) / 2
print(head(data_all))
print(dim(data_all))

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
ggplot(means, aes(x = map_score, y = sample_mean)) +
  geom_point(color = "#2980B9", size = 2) +
  geom_smooth(method = lm, se = FALSE, fullrange = TRUE, color = "#2C3E50") +
  ylab("Picomole") +
  xlab("Mappability score") +
  stat_cor(method = "spearman") +
  scale_y_continuous(labels = comma)
}

##Make plot
png("2020_conc_mapp_corr.png", height = 6, width = 6, units = "in", res = 300)
corr_plot(data_all)
dev.off()

###Plots with mean mappability score
##Make ggplot hexbin plot
png("2020_conc_mapp_corr_hexbin_morecounts.png", height = 6, width = 6, units = "in", res = 300)
hex <- hexbin(data_all$map_score, data_all$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count)) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Mappability score", y = "Picomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

##Looking into what the high fmol regions are
### subset the windows that have > 5000fmol

high_pmol <- subset(data_all, data_all$sample_mean > 5000)
rownames(high_pmol) <- high_pmol$window
#### 13 windows have high fmol

### UCSC simple repeats
repeats <- read.table("~/Projects/2018_PTB/data/2019_EPIC_HCT116_CtlProject/2020_repeatregions_mappedto300bpwindows.bed",
                      sep = "\t", header = FALSE)
###subset to regions where there are overlaps between our windows
repeats <- subset(repeats, repeats$V9 > 0)
repeats <- repeats[, c(1:3)]
colnames(repeats) <- c("chr", "start", "end")
repeats$window <- paste0(repeats$chr, "_", repeats$start, "_", repeats$end)
repeats <- unique(repeats)
rownames(repeats) <- repeats$window

##ENCODE blacklist regions
blacklist <- read.table("~/Projects/2018_PTB/data/Annotations/2020_new_ENCODE_blacklist.csv",
                        sep = "\t", header = FALSE)
colnames(blacklist) <- c("bl_chr", "bl_start", "bl_end", "chr", "start", "end", "overlap")
blacklist$window <- paste0(blacklist$chr, "_", blacklist$start, "_", blacklist$end)
blacklist <- unique(blacklist)

overlap_repetitiveregions <- high_pmol[rownames(high_pmol) %in% rownames(repeats), ]
##Every single one is on the UCSC repetitive regions black list

overlap_blacklist <- high_pmol[rownames(high_pmol) %in% rownames(blacklist), ]

##remove the blacklist regions and remake the hexbin plot.
rm <- repeats$window
rownames(data_all) <- data_all$window
data_noblacklist <- data_all[!rownames(data_all) %in% rm, ]

rm <- blacklist$window
data_noblacklist <- data_noblacklist[!rownames(data_noblacklist) %in% rm, ]

png("2020_conc_mapp_corr_hexbin_noblacklist_morecounts.png",
    height = 6, width = 6, units = "in", res = 300)
hex <- hexbin(data_noblacklist$map_score, data_noblacklist$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count)) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Mappability score", y = "Picomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

## Follow up on the high fmol sites not in the blacklist or  STR list
high_pmol2 <- subset(data_noblacklist, data_noblacklist$sample_mean > 400)
write.csv(high_pmol2, file = "2020_highpmol_noblacklist.csv")

###Plots with minimum mappability score
##Make ggplot hexbin plot
png("2020_conc_minmapp_corr_hexbin_morecounts.png",
    height = 6, width = 6, units = "in", res =300)
hex <- hexbin(data_all$map_score, data_all$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count)) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Mappability score", y = "Picomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

##Looking into what the high fmol regions are
### subset the windows that have > 5000fmol

high_pmol <- subset(data_all, data_all$sample_mean > 5000)
rownames(high_pmol) <- high_pmol$window
#### 13 windows have high fmol - same as above

##remove the blacklist regions and remake the hexbin plot.
rm <- repeats$window
rownames(data_all) <- data_all$window
data_noblacklist <- data_all[!rownames(data_all) %in% rm, ]

rm <- blacklist$window
data_noblacklist <- data_noblacklist[!rownames(data_noblacklist) %in% rm,]

png("2020_conc_minmapp_corr_hexbin_noblacklist_morecounts.png",
    height = 6, width = 6, units = "in", res = 300)
hex <- hexbin(data_noblacklist$map_score, data_noblacklist$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count)) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Mappability score", y = "Picomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

## Follow up on the high fmol sites not in the blacklist or  STR list
high_pmol2 <- subset(data_noblacklist, data_noblacklist$sample_mean > 400)
write.csv(high_pmol2, file = "2020_highfmol_minmappability_noblacklist.csv")

##Hexbin plot with standard deviation
data_noblacklist$sd <- apply(data_noblacklist[, c("sample1", "sample2")], 1, sd)

png("2020_fmol_sd_corr.png", height = 6, width = 6, units = "in", res = 300)
hex <- hexbin(data_noblacklist$sd, data_noblacklist$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
ggplot(gghex, aes(x = x, y = y, fill = count)) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_continuous(type = "viridis") +
  labs(x = "Standard deviation", y = "Picomoles", fill = "Log of counts") +
  scale_y_continuous(labels = comma)
dev.off()

rm <- subset(data_noblacklist, data_noblacklist$sd > 0.25 | data_noblacklist$map_score < 0.5)
data_filtered <- data_noblacklist[!rownames(data_noblacklist) %in% rownames(rm), ]

##Set up colour palette
pal <- brewer.pal(n = 8, name = "Blues")[3:8]

png("2020_fmol_sd_filtered_pubplt.png", height = 6, width = 6, units = "in", res = 300)
hex <- hexbin(data_filtered$sd, data_filtered$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
pmol_sd <- ggplot(gghex, aes(x = x, y = y, fill = count)) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_gradientn(colors = pal,
                       labels = c(expression(10^0), expression(10^5),
                                  expression(10^10), expression(10^15))) +
  labs(x = "Standard deviation", y = "Picomoles", fill = "Reads") +
  scale_y_continuous(labels = comma)
pmol_sd
dev.off()


##Redo hebin mappability without the sd
png("2020_fmol_mapp_filtered_pubplt.png", height = 6, width = 6, units = "in", res = 300)
hex <- hexbin(data_filtered$map_score, data_filtered$sample_mean, 20)
   gghex <- data.frame(hcell2xy(hex), count = log(hex@count),
                        xo = hex@xcm, yo = hex@ycm)
pmol_mapp <- ggplot(gghex, aes(x = x, y = y, fill = count)) +
  geom_hex(bins = 40, stat = "identity") +
  scale_fill_gradientn(colors = pal,
                       labels = c(expression(10^0), expression(10^5),
                                  expression(10^10), expression(10^15))) +
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
png("2020_fmol_sd_mapp_corr_pubplt.png", height = 6, width = 12, units = "in", res = 300)
mylegend <- g_legend(pmol_sd)
grid.arrange(arrangeGrob(pmol_sd + theme(legend.position = "none")),
       arrangeGrob(pmol_mapp + theme(legend.position = "none")), nrow = 1, mylegend)
dev.off()

high_pmol <- data_filtered[data_filtered$sample_mean > 2, ]
