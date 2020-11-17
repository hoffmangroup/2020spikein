#!/user/bin/env Rscript

#Load libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(stringr)
library(varhandle)
library(forcats)
library(scales)

##Function of the wanted plot
####Where x is the data points and y is the dataframe of means
dist_plots <- function(x, y) {
  # New facet label names for supp variable
  facet_labs <- c("G+C: 35%", "G+C: 50%", "G+C: 65%")
  names(facet_labs) <- c("35", "50", "65")

  ggplot(x, aes(x = CpGfrac, y = n, group = variable)) +
    geom_point(aes(color = methylation_status, shape = variable), size = 5) +
    geom_point(data = y %>% filter(methylation_status == "unmethylated"),
           mapping = aes(x = CpGfrac, y = means, group = methylation_status),
           pch = "—", col = "black", size = 10) +
    geom_point(data = y %>% filter(methylation_status == "methylated"),
           mapping = aes(x = CpGfrac, y = means, group = methylation_status),
           pch = "—", col = "royalblue", size = 10) +
  facet_grid(cols = vars(GC), labeller = labeller(GC = facet_labs)) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 20),
      axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 16),
      axis.text.y = element_text(vjust = 0.5, size = 16),
      axis.title.x = element_text(vjust = 0.5, size = 20),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      strip.text.x = element_text(size = 16),
      plot.title =
        element_text(hjust = 0.5, color = "black", size = 14)) +
  ylab("Number of unique reads") +
  xlab("CpG fraction") +
  theme(legend.position = "right") +
  labs(color = "Methylation status", shape = "Sample") +
  scale_colour_manual(values = c("deepskyblue3", "darkgrey")) +
  scale_shape_manual(values = c(1, 2)) +
  scale_y_continuous(labels = comma)
}

##Formatting data to what the input for the plot function needs to be
###Specify the concentration of the data you are passing
  ###This is used to label the output plots
data_wrang <- function(x) {
###Calculate CpG fraction
x$cpg_frac <- as.factor(x$cpg / x$fragment_len)

###revalue to make labels on plot nicer to read
x$cpg_frac <- revalue(x$cpg_frac,
                      c("0.0125" = "1/80", "0.025" = "1/40", "0.05" = "1/20"))
x$methylation_status <- revalue(x$methylation_status,
                                c("meth" = "methylated",
                                  "unmeth" = "unmethylated"))

##divide data into fragment length category
data_80bp <- subset(x, x$fragment_len == "80")
data_160bp <- subset(x, x$fragment_len == "160")
data_320bp <- subset(x, x$fragment_len == "320")

##melt data for each fragment length
##each row is a read and variable lists the sample name
data_80bp_melt <- melt(data_80bp,
                       id = c("fragment_len", "GC",
                              "cpg_frac", "methylation_status"),
                       measure.vars =
                         grep("read_count_", colnames(data_80bp)))
data_160bp_melt <- melt(data_160bp,
                        id = c("fragment_len", "GC",
                               "cpg_frac", "methylation_status"),
                        measure.vars =
                          grep("read_count_", colnames(data_160bp)))
data_320bp_melt <- melt(data_320bp,
                        id = c("fragment_len", "GC",
                               "cpg_frac", "methylation_status"),
                        measure.vars =
                          grep("read_count_", colnames(data_160bp)))

###revalue to make labels on plot nicer to read for each fragment length
###relabels sample number to "Sample 1" and "Sample 2"
  names <- as.matrix(unique(data_80bp_melt$variable))
#Just pulls the unique sample names from any data frame
  relabel_samplenames <- function(sample_names, names) {
  #relabels all samples
  list_samples <- list()
  for (i in 1:seq_len(sample_names)) {
    list_samples[[i]] <- ifelse(sample_names[i] == names[1],
                                paste0("Sample 1"), paste0("Sample 2"))
  }
  return(list_samples)
}

data_80bp_melt$variable <- do.call(rbind,
                                   relabel_samplenames(data_80bp_melt$variable,
                                                       names))
data_160bp_melt$variable <- do.call(rbind,
                                   relabel_samplenames(data_160bp_melt$variable,
                                                        names))
data_320bp_melt$variable <- do.call(rbind,
                                    relabel_samplenames(data_320bp_melt$variable,
                                                        names))

##Produce total counts for each grouping of fragment length
data_80bp_totcount <- data_80bp_melt %>%
  count(value, variable, cpg_frac, GC, methylation_status)
###Produces total with all the total counts I need
data_80bp_totcount <- na.omit(data_80bp_totcount)

data_160bp_totcount <- data_160bp_melt %>%
  count(value, variable, cpg_frac, GC, methylation_status)
###Produces total with all the total counts I need
data_160bp_totcount <- na.omit(data_160bp_totcount)

data_320bp_totcount <- data_320bp_melt %>%
  count(value, variable, CpGfrac, GC, methylation_status)
###Produces total with all the total counts I need
data_320bp_totcount <- na.omit(data_320bp_totcount)

##Computing mean values for each fragment length group
data_80bp_means <- data_80bp_totcount %>%
        group_by(methylation_status, GC, cpg_frac) %>%
        summarise(means = mean(n))

data_160bp_means <- data_160bp_totcount %>%
        group_by(methylation_status, GC, cpg_frac) %>%
        summarise(means = mean(n))

data_320bp_means <- data_320bp_totcount %>%
        group_by(methylation_status, GC, cpg_frac) %>%
        summarise(means = mean(n))

return(list(data_80bp_totcount, data_80bp_means,
            data_160bp_totcount, data_160bp_means,
            data_320bp_totcount, data_320bp_means))
#Return a list of tibbles that are fragment length counts, and means
}

##Read in the files and add a column with sample name and specify read_counts column per sample
##Data from spike-ins into HCT116
data_01 <- read.table("2019_synfrag0.1_spikeincellline_dedup_lowconc.csv", header = T)
data_005 <- read.table("2019_synfrag0.05_spikeincellline_dedup_lowconc.csv", header = T)
data_001 <- read.table("2019_synfrag0.01_spikeincellline_dedup_lowconc.csv", header = T)

##Extracting legend, so all plots are same size in the grid.arrange()
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a_gplot) {
  tmp <- ggplot_gtable(ggplot_build(a_gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
  }

##0.1ng spike-in
data_01_cleaned <- data_wrang(data_01)

png("2019_0.1ng_plots.png",
    height = 10, width = 24, units = "in", res = 300)
plot80_01 <- dist_plots(as.data.frame(data_01_cleaned[1]),
                         as.data.frame(data_01_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
plot160_01 <- dist_plots(as.data.frame(data_01_cleaned[3]),
                          as.data.frame(data_01_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
plot320_01 <- dist_plots(as.data.frame(data_01_cleaned[5]),
                          as.data.frame(data_01_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
mylegend <- g_legend(plot80)
grid.arrange(arrangeGrob(plot80_01 + theme(legend.position = "none"),
                         plot160_01 + theme(legend.position = "none"),
                         plot320_01 + theme(legend.position = "none"),
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##0.05ng spike-in
data_005_cleaned <- data_wrang(data_005)

png("2019_0.05ng_plots.png",
    height = 10, width = 24, units = "in", res = 300)
plot80_005 <- dist_plots(as.data.frame(data_005_cleaned[1]),
                          as.data.frame(data_005_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
plot160_005 <- dist_plots(as.data.frame(data_005_cleaned[3]),
                           as.data.frame(data_005_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
plot320_005 <- dist_plots(as.data.frame(data_005_cleaned[5]),
                           as.data.frame(data_005_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
mylegend <- g_legend(plot80)
grid.arrange(arrangeGrob(plot80_005 + theme(legend.position = "none"),
                         plot160_005 + theme(legend.position = "none"),
                         plot320_005 + theme(legend.position = "none"),
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##0.01ng spike-in
data_001_cleaned <- data_wrang(data_001)

png("2019_0.01ng_plots.png",
    height = 10, width = 24, units = "in", res = 300)
plot80_001 <- dist_plots(as.data.frame(data_001_cleaned[1]),
                          as.data.frame(data_001_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
plot160_001 <- dist_plots(as.data.frame(data_001_cleaned[3]),
                           as.data.frame(data_001_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
plot360_001 <- dist_plots(as.data.frame(data_001_cleaned[5]),
                           as.data.frame(data_001_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000),
                     breaks = c(0, 5000, 10000, 15000, 20000))
mylegend <- g_legend(plot80_001)
grid.arrange(arrangeGrob(plot80_001 + theme(legend.position = "none"),
                         plot160_001 + theme(legend.position = "none"),
                         plot360_001 + theme(legend.position = "none"),
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##Load in the synthetic spike-in data on miseq
## This is the input and outputs
## Sequencing on only the spike-in controls
data_input <- read.table("~/Projects/2018_PTB/data/2019_TrimmedData_SynFragOly/2019_synfrag_input_dedup.csv", header = T)
data_output <- read.table("~/Projects/2018_PTB/data/2019_TrimmedData_SynFragOly/2019_synfrag_output_dedup.csv", header = T)

##input samples on miseq
data_input_cleaned <- data_wrang(data_input)

png("~/Projects/2018_PTB/data/2019_TrimmedData_SynFragOly/2019_input_miseq_plots.png",
    height = 10, width = 24, units = "in", res = 300)
plot80_input <- dist_plots(as.data.frame(data_input_cleaned[1]),
                           as.data.frame(data_input_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 10000),
                     breaks = c(0, 2000, 4000, 6000, 8000, 10000))
plot160_input <- dist_plots(as.data.frame(data_input_cleaned[3]),
                            as.data.frame(data_input_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 10000),
                     breaks = c(0, 2000, 4000, 6000, 8000, 10000))
plot320_input <- dist_plots(as.data.frame(data_input_cleaned[5]),
                            as.data.frame(data_input_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 10000),
                     breaks = c(0, 2000, 4000, 6000, 8000, 10000))
mylegend <- g_legend(plot80_input)
grid.arrange(arrangeGrob(plot80_input + theme(legend.position = "none"),
                         plot160_input + theme(legend.position = "none"),
                         plot320_input + theme(legend.position = "none"),
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##output samples on miseq
data_output_cleaned <- data_wrang(data_output)

png("~/Projects/2018_PTB/data/2019_TrimmedData_SynFragOly/2019_output_miseq_plots.png",
    height = 10, width = 24, units = "in", res = 300)
plot80_output <- dist_plots(as.data.frame(data_output_cleaned[1]),
                            as.data.frame(data_output_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 10000),
                     breaks = c(0, 2000, 4000, 6000, 8000, 10000))
plot160_output <- dist_plots(as.data.frame(data_output_cleaned[3]),
                             as.data.frame(data_output_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 10000),
                     breaks = c(0, 2000, 4000, 6000, 8000, 10000))
plot320_output <- dist_plots(as.data.frame(data_output_cleaned[5]),
                             as.data.frame(data_output_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 10000),
                     breaks = c(0, 2000, 4000, 6000, 8000, 10000))
mylegend <- g_legend(plot80_output)
grid.arrange(arrangeGrob(plot80_output + theme(legend.position = "none"),
                         plot160_output + theme(legend.position = "none"),
                         plot320_output + theme(legend.position = "none"),
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##Putting plots in order of what I want for publication
png("2020_frag_GC_CpG_pubplts.png", height = 25, width = 20, units = "in", res = 300)
mylegend <- g_legend(plot80_001)
plots_pub <- grid.arrange(arrangeGrob(plot80_input + theme(legend.position = "none", axis.title.y = element_text("Reads")) +
                                                scale_y_continuous(labels = comma, limits = c(0, 10000),
                                                                   breaks = c(0, 2500, 5000, 7500, 10000)),
                                        plot160_input + theme(axis.title.y = element_blank(), legend.position = "none") +
                                                scale_y_continuous(labels = comma, limits = c(0, 10000),
                                                                   breaks = c(0, 2500, 5000, 7500, 10000)),
                                        plot320_input + theme(axis.title.y = element_blank(), legend.position = "none") +
                                                scale_y_continuous(labels = comma, limits = c(0, 10000),
                                                                   breaks = c(0, 2500, 5000, 7500, 10000)),
                                                ncol = 3),
                                        arrangeGrob(plot80_output + theme(legend.position = "none", axis.title.y = element_text("Reads")) +
                                                scale_y_continuous(labels = comma, limits = c(0, 10000),
                                                                   breaks = c(0, 2500, 5000, 7500, 10000)),
                                        plot160_output + theme(axis.title.y = element_blank(), legend.position = "none") +
                                                scale_y_continuous(labels = comma, limits = c(0, 10000),
                                                                   breaks = c(0, 2500, 5000, 7500, 10000)),
                                        plot320_output + theme(axis.title.y = element_blank(), legend.position = "none") +
                                                scale_y_continuous(labels = comma, limits = c(0, 10000),
                                                                   breaks = c(0, 2500, 5000, 7500, 10000)),
                                                ncol = 3),
                                        arrangeGrob(plot80_001 + theme(legend.position = "none", axis.title.y = element_text("Reads")) +
            scale_y_continuous(labels = comma, limits = c(0, 10000), breaks = c(0, 2500, 5000, 7500, 10000)),
                                        plot160_001 + theme(axis.title.y = element_blank(), legend.position = "none") +
            scale_y_continuous(labels = comma, limits = c(0, 10000), breaks = c(0, 2500, 5000, 7500, 10000)),
                                        plot360_001 + theme(axis.title.y = element_blank(), legend.position = "none") +
            scale_y_continuous(labels = comma, limits = c(0, 10000), breaks = c(0, 2500, 5000, 7500, 10000)),
                                                ncol = 3),
                                                ncol = 1)
dev.off()