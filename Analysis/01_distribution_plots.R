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
library(data.table)

##Function of the wanted plot
####Where x is the data points and y is the dataframe of means
dist_plots <- function(x, y) {
# New facet label names for supp variable
facet_labs <- c("G+C: 35%", "G+C: 50%", "G+C: 65%")
names(facet_labs) <- c("35", "50", "65")

ggplot(x, aes(x = cpg_frac, y = value, group = frag_grp)) +
geom_point(aes(color = methylation_status, shape = variable), size = 5) +
geom_point(data = y %>% filter(methylation_status == "unmethylated"),
           mapping = aes(x = cpg_frac, y = means, group = methylation_status),
           pch = "—", col = "black", size = 10) +
geom_point(data = y %>% filter(methylation_status == "methylated"),
           mapping = aes(x = cpg_frac, y = means, group = methylation_status),
           pch = "—", col = "royalblue", size = 10) +
facet_grid(cols = vars(GC), scales = "free", labeller = labeller(GC = facet_labs)) +
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
ylab("Reads") +
xlab("CpG fraction") +
theme(legend.position = "right") +
labs(color = "Methylation status", shape = "Sample") +
scale_colour_manual(values = c("deepskyblue3", "darkgrey")) +
scale_shape_manual(values = c(1, 2)) +
scale_y_continuous(labels = comma) + 
scale_x_discrete(limits = levels(x$fragment_grp))
}

##Formatting data to what the input for the plot function needs to be
###Specify the concentration of the data you are passing
  ###This is used to label the output plots
data_wrang <- function(x) {
###Calculate CpG fraction
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
                       id = c("frag_grp", "fragment_len", "GC",
                              "cpg_frac", "methylation_status"),
                       measure.vars =
                         grep("read_count", colnames(data_80bp)))
data_160bp_melt <- melt(data_160bp,
                        id = c("frag_grp", "fragment_len", "GC",
                               "cpg_frac", "methylation_status"),
                        measure.vars =
                          grep("read_count", colnames(data_160bp)))
data_320bp_melt <- melt(data_320bp,
                        id = c("frag_grp", "fragment_len", "GC",
                               "cpg_frac", "methylation_status"),
                        measure.vars =
                          grep("read_count", colnames(data_320bp)))

###revalue to make labels on plot nicer to read for each fragment length
###relabels sample number to "Sample 1" and "Sample 2"
names <- as.matrix(unique(data_80bp_melt$variable))
#Just pulls the unique sample names from any data frame
relabel_samplenames <- function(sample_names,names){#relabels all samples
  list_samples <- list()
  for(i in 1:length(sample_names))
  {
    list_samples[[i]]=ifelse (sample_names[i]==names[1],paste0("Sample 1"),paste0("Sample 2"))
  }
  return(list_samples)
}

data_80bp_melt$variable <- do.call(rbind, relabel_samplenames(data_80bp_melt$variable, names))
data_160bp_melt$variable <- do.call(rbind, relabel_samplenames(data_160bp_melt$variable, names))
data_320bp_melt$variable <- do.call(rbind, relabel_samplenames(data_320bp_melt$variable, names))

##Produce total counts for each grouping of fragment length
data_80bp_totcount <- data_80bp_melt %>%
  count(value, variable, frag_grp, cpg_frac, GC, methylation_status)
###Produces total with all the total counts I need
data_80bp_totcount <- na.omit(data_80bp_totcount)

data_160bp_totcount <- data_160bp_melt %>%
  count(value, variable, frag_grp, cpg_frac, GC, methylation_status)
###Produces total with all the total counts I need
data_160bp_totcount <- na.omit(data_160bp_totcount)

data_320bp_totcount <- data_320bp_melt %>%
  count(value, variable, frag_grp, cpg_frac, GC, methylation_status)
###Produces total with all the total counts I need
data_320bp_totcount <- na.omit(data_320bp_totcount)

##Computing mean values for each fragment length group
data_80bp_means <- data_80bp_totcount %>%
        group_by(methylation_status, frag_grp, GC, cpg_frac) %>%
        summarise(means = mean(value), mad = mad(value), min = min(value), max = max(value), diff = max(value) - min(value), .groups = 'keep')

data_160bp_means <- data_160bp_totcount %>%
        group_by(methylation_status, frag_grp, GC, cpg_frac) %>%
        summarise(means = mean(value),  mad = mad(value), min = min(value), max = max(value), diff = max(value) - min(value), .groups = 'keep')

data_320bp_means <- data_320bp_totcount %>%
        group_by(methylation_status, frag_grp, GC, cpg_frac) %>%
        summarise(means = mean(value),  mad = mad(value), min = min(value), max = max(value), diff = max(value) - min(value), .groups = 'keep')

return(list(data_80bp_totcount, data_80bp_means,
            data_160bp_totcount, data_160bp_means,
            data_320bp_totcount, data_320bp_means))
#Return a list of tibbles that are fragment length counts, and means
}

##Read in the files and add a column with sample name and specify read_counts column per sample
##0.1ng spike-in
S1_01 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_aligned.trimmed.6543_S7_L_R1_001.fastqtrimmed.6543_S7_L_R2_001.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S1_01 <- S1_01[c(1:2)]
colnames(S1_01) <- c("label", "S1_01_read_count")
S1_01 <- S1_01 %>% group_by(label) %>% summarise_if(is.numeric, sum)

S2_01 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_aligned.trimmed.6544_S21_L_R1_001.fastqtrimmed.6544_S21_L_R2_001.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S2_01 <- S2_01[c(1:2)]
colnames(S2_01) <- c("label", "S2_01_read_count")
S2_01 <- S2_01 %>% group_by(label) %>% summarise_if(is.numeric, sum)

data_0.1 <- merge(S1_01, S2_01, by = "label", all = TRUE, fill = TRUE)

##0.05ng spike-in
S1_005 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_aligned.trimmed.6545_S1_L_R1_001.fastqtrimmed.6545_S1_L_R2_001.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S1_005 <- S1_005[c(1:2)]
colnames(S1_005) <- c("label", "S1_005_read_count")
S1_005 <- S1_005 %>% group_by(label) %>% summarise_if(is.numeric, sum)

S2_005 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_aligned.trimmed.6546_S2_L_R1_001.fastqtrimmed.6546_S2_L_R2_001.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S2_005 <- S2_005[c(1:2)]
colnames(S2_005) <- c("label", "S2_005_read_count")
S2_005 <- S2_005 %>% group_by(label) %>% summarise_if(is.numeric, sum)

data_0.05 <- merge(S1_005, S2_005, by = "label", all = TRUE, fill = TRUE)

##0.01ng spike-in
S1_001 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_aligned.trimmed.6547_S3_L_R1_001.fastqtrimmed.6547_S3_L_R2_001.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S1_001 <- S1_001[c(1:2)]
colnames(S1_001) <- c("label", "S1_001_read_count")
S1_001 <- S1_001 %>% group_by(label) %>% summarise_if(is.numeric, sum)

S2_001 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_aligned.trimmed.6548_S4_L_R1_001.fastqtrimmed.6548_S4_L_R2_001.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S2_001 <- S2_001[c(1:2)]
colnames(S2_001) <- c("label", "S2_001_read_count")
S2_001 <- S2_001 %>% group_by(label) %>% summarise_if(is.numeric, sum)

data_0.01 <- merge(S1_001, S2_001, by = "label", all = TRUE, fill = TRUE)
data_0.01[is.na(data_0.01)] <- 0

#create methylation status and frag_grp
data_0.1$frag_grp <- substr(data_0.1$label, 1, 10)

data_0.1$methylation_status <- for (i in 1:nrow(data_0.1)) {
				if (substr(data_0.1$label[i], (nchar(data_0.1$label)[i]+1)-4, nchar(data_0.1$label[i])) == "meth") {data_0.1[i,5] = "meth"}
				else {data_0.1[i,5] = "unmeth"}}
colnames(data_0.1) <- c("label", "S1_01_read_count", "S2_01_read_count", "frag_grp", "methylation_status")


data_0.05$frag_grp <- substr(data_0.05$label, 1, 10)

data_0.05$methylation_status <- for (i in 1:nrow(data_0.05)) {
                                if (substr(data_0.05$label[i], (nchar(data_0.05$label)[i]+1)-4, nchar(data_0.05$label[i])) == "meth") {data_0.05[i,5] = "meth"}
                                else {data_0.05[i,5] = "unmeth"}}
colnames(data_0.05) <- c("label", "S1_005_read_count", "S2_005_read_count", "frag_grp", "methylation_status")

data_0.01$frag_grp <- substr(data_0.01$label, 1, 10)

data_0.01$methylation_status <- for (i in 1:nrow(data_0.01)) {
                                if (substr(data_0.01$label[i], (nchar(data_0.01$label)[i]+1)-4, nchar(data_0.01$label[i])) == "meth") {data_0.01[i,5] = "meth"}
                                else {data_0.01[i,5] = "unmeth"}}
colnames(data_0.01) <- c("label", "S1_001_read_count", "S2_001_read_count", "frag_grp", "methylation_status")

#clean up frag_grp
data_0.1$frag_grp <- as.factor(data_0.1$frag_grp)
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("320b_16C_3" = "320b_16C_35"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("320b_16C_5" = "320b_16C_50"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("320b_16C_6" = "320b_16C_65"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_1C_35G" = "80b_1C_35"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_1C_50G" = "80b_1C_50"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_1C_65G" = "80b_1C_65"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_2C_35G" = "80b_2C_35"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_2C_50G" = "80b_2C_50"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_2C_65G" = "80b_2C_65"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_4C_35G" = "80b_4C_35"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_4C_50G" = "80b_4C_50"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("80b_4C_65G" = "80b_4C_65"))

data_0.05$frag_grp <- as.factor(data_0.05$frag_grp)
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("320b_16C_3" = "320b_16C_35"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("320b_16C_5" = "320b_16C_50"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("320b_16C_6" = "320b_16C_65"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_1C_35G" = "80b_1C_35"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_1C_50G" = "80b_1C_50"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_1C_65G" = "80b_1C_65"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_2C_35G" = "80b_2C_35"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_2C_50G" = "80b_2C_50"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_2C_65G" = "80b_2C_65"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_4C_35G" = "80b_4C_35"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_4C_50G" = "80b_4C_50"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("80b_4C_65G" = "80b_4C_65"))

data_0.01$frag_grp <- as.factor(data_0.01$frag_grp)
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("320b_16C_3" = "320b_16C_35"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("320b_16C_5" = "320b_16C_50"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("320b_16C_6" = "320b_16C_65"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_1C_35G" = "80b_1C_35"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_1C_50G" = "80b_1C_50"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_1C_65G" = "80b_1C_65"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_2C_35G" = "80b_2C_35"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_2C_50G" = "80b_2C_50"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_2C_65G" = "80b_2C_65"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_4C_35G" = "80b_4C_35"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_4C_50G" = "80b_4C_50"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("80b_4C_65G" = "80b_4C_65"))

##since initial 320bp, 65% G+C is not correct, let's correct that to correct G+C content
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("320b_4C_65" = "320b_04C_35"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("320b_8C_65" = "320b_08C_50"))
data_0.1$frag_grp <- revalue(data_0.1$frag_grp, c("320b_16C_65" = "320b_016C_50"))

data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("320b_4C_65" = "320b_04C_35"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("320b_8C_65" = "320b_08C_50"))
data_0.05$frag_grp <- revalue(data_0.05$frag_grp, c("320b_16C_65" = "320b_016C_50"))

data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("320b_4C_65" = "320b_04C_35"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("320b_8C_65" = "320b_08C_50"))
data_0.01$frag_grp <- revalue(data_0.01$frag_grp, c("320b_16C_65" = "320b_016C_50"))

#Cpg fraction
data_0.1$CpG <- as.factor(substr(data_0.1$label, 4,7))
data_0.1$CpG <- revalue(data_0.1$CpG, c("_1C_" = "1"))
data_0.1$CpG <- revalue(data_0.1$CpG, c("_2C_" = "2"))
data_0.1$CpG <- revalue(data_0.1$CpG, c("_4C_" = "4"))
data_0.1$CpG <- revalue(data_0.1$CpG, c("b_16" = "16"))
data_0.1$CpG <- revalue(data_0.1$CpG, c("b_2C" = "2"))
data_0.1$CpG <- revalue(data_0.1$CpG, c("b_4C" = "4"))
data_0.1$CpG <- revalue(data_0.1$CpG, c("b_8C" = "8"))
data_0.1$CpG <- as.numeric(as.character(data_0.1$CpG))

data_0.05$CpG <- as.factor(substr(data_0.05$label, 4,7))
data_0.05$CpG <- revalue(data_0.05$CpG, c("_1C_" = "1"))
data_0.05$CpG <- revalue(data_0.05$CpG, c("_2C_" = "2"))
data_0.05$CpG <- revalue(data_0.05$CpG, c("_4C_" = "4"))
data_0.05$CpG <- revalue(data_0.05$CpG, c("b_16" = "16"))
data_0.05$CpG <- revalue(data_0.05$CpG, c("b_2C" = "2"))
data_0.05$CpG <- revalue(data_0.05$CpG, c("b_4C" = "4"))
data_0.05$CpG <- revalue(data_0.05$CpG, c("b_8C" = "8"))
data_0.05$CpG <- as.numeric(as.character(data_0.05$CpG))

data_0.01$CpG <- as.factor(substr(data_0.01$label, 4,7))
data_0.01$CpG <- revalue(data_0.01$CpG, c("_1C_" = "1"))
data_0.01$CpG <- revalue(data_0.01$CpG, c("_2C_" = "2"))
data_0.01$CpG <- revalue(data_0.01$CpG, c("_4C_" = "4"))
data_0.01$CpG <- revalue(data_0.01$CpG, c("b_16" = "16"))
data_0.01$CpG <- revalue(data_0.01$CpG, c("b_2C" = "2"))
data_0.01$CpG <- revalue(data_0.01$CpG, c("b_4C" = "4"))
data_0.01$CpG <- revalue(data_0.01$CpG, c("b_8C" = "8"))
data_0.01$CpG <- as.numeric(as.character(data_0.01$CpG))

data_0.1$fragment_len <- as.factor(substr(data_0.1$label, 1,3))
data_0.1$fragment_len <- revalue(data_0.1$fragment_len, c("80b" = "80"))

data_0.05$fragment_len <- as.factor(substr(data_0.05$label, 1,3))
data_0.05$fragment_len <- revalue(data_0.05$fragment_len, c("80b" = "80"))

data_0.01$fragment_len <- as.factor(substr(data_0.01$label, 1,3))
data_0.01$fragment_len <- revalue(data_0.01$fragment_len, c("80b" = "80"))

data_0.1$cpg_frac <- as.factor(as.numeric(data_0.1$CpG) / as.numeric(as.character(data_0.1$fragment_len)))
data_0.05$cpg_frac <- as.factor(as.numeric(data_0.05$CpG) / as.numeric(as.character(data_0.05$fragment_len)))
data_0.01$cpg_frac <- as.factor(as.numeric(data_0.01$CpG) / as.numeric(as.character(data_0.01$fragment_len)))

###revalue to make labels on plot nicer to read
set_cpg_frac <- function(x) {
			revalue(x$cpg_frac,
                      c("0.0125" = "1/80", "0.025" = "1/40", "0.05" = "1/20"))}

data_0.1$cpg_frac <- set_cpg_frac(data_0.1)
data_0.05$cpg_frac <- set_cpg_frac(data_0.05)
data_0.01$cpg_frac <- set_cpg_frac(data_0.01)

##Add two states for CpG fractions to visualize the replicate fragments better
setDT(data_0.1)[frag_grp == "320b_04C_35" & cpg_frac == "1/80", cpg_frac := "1/81"]
setDT(data_0.05)[frag_grp == "320b_04C_35" & cpg_frac == "1/80", cpg_frac := "1/81"]
setDT(data_0.01)[frag_grp == "320b_04C_35" & cpg_frac == "1/80", cpg_frac := "1/81"]

setDT(data_0.1)[frag_grp == "320b_08C_50" & cpg_frac == "1/40", cpg_frac := "1/41"]
setDT(data_0.05)[frag_grp == "320b_08C_50" & cpg_frac == "1/40", cpg_frac := "1/41"]
setDT(data_0.01)[frag_grp == "320b_08C_50" & cpg_frac == "1/40", cpg_frac := "1/41"]

setDT(data_0.1)[frag_grp == "320b_016C_50" & cpg_frac == "1/20", cpg_frac := "1/21"]
setDT(data_0.05)[frag_grp == "320b_016C_50" & cpg_frac == "1/20", cpg_frac := "1/21"]
setDT(data_0.01)[frag_grp == "320b_016C_50" & cpg_frac == "1/20", cpg_frac := "1/21"]


##Fix to the correct GC
data_0.1$GC <- as.factor(str_sub(data_0.1$frag_grp,-2,-1))
data_0.05$GC <- as.factor(str_sub(data_0.05$frag_grp,-2,-1))
data_0.01$GC <- as.factor(str_sub(data_0.01$frag_grp,-2,-1))

##Extracting legend, so all plots are same size in the grid.arrange()
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend <- function(a_gplot) {
  tmp <- ggplot_gtable(ggplot_build(a_gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
  }

##0.1ng spike-in
levels(data_0.1$cpg_frac) <-  c("1/80", "1/40", "1/20", "1/80(2)", "1/40(2)", "1/20(2)")
data_0.1$cpg_frac <- factor(data_0.1$cpg_frac, levels = c("1/80", "1/80(2)", "1/40", "1/40(2)", "1/20", "1/20(2)"))
data_0.1_cleaned <- data_wrang(data_0.1)

png("2021_0.1ng_plots.png",
    height = 10, width = 24, units = "in", res = 300)
plot80_0.1 <- dist_plots(as.data.frame(data_0.1_cleaned[1]), as.data.frame(data_0.1_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
plot160_0.1 <- dist_plots(as.data.frame(data_0.1_cleaned[3]), as.data.frame(data_0.1_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
plot320_0.1 <- dist_plots(as.data.frame(data_0.1_cleaned[5]), as.data.frame(data_0.1_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
mylegend <- g_legend(plot80_0.1)
grid.arrange(arrangeGrob(plot80_0.1 + theme(legend.position = "none"),
                         plot160_0.1 + theme(legend.position = "none"), ncol = 2, nrow = 1),
                         arrangeGrob(plot320_0.1 + theme(legend.position = "none"), ncol = 3.35,
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##0.05ng spike-in
levels(data_0.05$cpg_frac) <-  c("1/80", "1/40", "1/20", "1/80(2)", "1/40(2)", "1/20(2)")
data_0.05$cpg_frac <- factor(data_0.05$cpg_frac, levels = c("1/80", "1/80(2)", "1/40", "1/40(2)", "1/20", "1/20(2)"))
data_0.05_cleaned <- data_wrang(data_0.05)

png("2021_0.05ng_plots.png",
    height = 10, width = 24, units = "in", res = 300)
plot80_0.05 <- dist_plots(as.data.frame(data_0.05_cleaned[1]), as.data.frame(data_0.05_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
plot160_0.05 <- dist_plots(as.data.frame(data_0.05_cleaned[3]), as.data.frame(data_0.05_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
plot320_0.05 <- dist_plots(as.data.frame(data_0.05_cleaned[5]), as.data.frame(data_0.05_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
mylegend <- g_legend(plot80_0.05)
grid.arrange(arrangeGrob(plot80_0.1 + theme(legend.position = "none"),
                         plot160_0.1 + theme(legend.position = "none"), ncol = 2, nrow = 1),
                         arrangeGrob(plot320_0.1 + theme(legend.position = "none"), ncol = 3.35,
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##0.01ng spike-in
levels(data_0.01$cpg_frac) <-  c("1/80", "1/40", "1/20", "1/80(2)", "1/40(2)", "1/20(2)")
data_0.01$cpg_frac <- factor(data_0.01$cpg_frac, levels = c("1/80", "1/80(2)", "1/40", "1/40(2)", "1/20", "1/20(2)"))
data_0.01_cleaned <- data_wrang(data_0.01)
save(data_0.01_cleaned, file = "2021_001ng_data.RData")

png("2021_0.01ng_plots.png",
    height = 10, width = 36, units = "in", res = 300)
plot80_0.01 <- dist_plots(as.data.frame(data_0.01_cleaned[1]), as.data.frame(data_0.01_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
plot160_0.01 <- dist_plots(as.data.frame(data_0.01_cleaned[3]), as.data.frame(data_0.01_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
plot320_0.01 <- dist_plots(as.data.frame(data_0.01_cleaned[5]), as.data.frame(data_0.01_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000))
mylegend <- g_legend(plot80_0.01)
grid.arrange(arrangeGrob(plot80_0.01 + theme(legend.position = "none"),
                         plot160_0.01 + theme(legend.position = "none"), ncol = 2, nrow = 1),
                         arrangeGrob(plot320_0.01 + theme(legend.position = "none"), ncol = 2.5,
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##Load in the synthetic spike-in data on miseq
S1_input <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2019_synfrag_only/filtered_aligned.trimmed.6427_S3.fastqtrimmed.6427_S3.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S1_input <- S1_input[c(1:2)]
colnames(S1_input) <- c("label", "S1_input_read_count")
S1_input <- S1_input %>% group_by(label) %>% summarise_if(is.numeric, sum)

S2_input <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2019_synfrag_only/filtered_aligned.trimmed.6428_S4.fastqtrimmed.6428_S4.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S2_input <- S2_input[c(1:2)]
colnames(S2_input) <- c("label", "S2_input_read_count")
S2_input <- S2_input %>% group_by(label) %>% summarise_if(is.numeric, sum)

data_input <- merge(S1_input, S2_input, by = "label", all = TRUE, fill = TRUE)
data_input[is.na(data_input)] <- 0

S1_output <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2019_synfrag_only/filtered_aligned.trimmed.6425_S1.fastqtrimmed.6425_S1.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S1_output <- S1_output[c(1:2)]
colnames(S1_output) <- c("label", "S1_output_read_count")
S1_output <- S1_output %>% group_by(label) %>% summarise_if(is.numeric, sum)

S2_output <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2019_synfrag_only/filtered_aligned.trimmed.6426_S2.fastqtrimmed.6426_S2.fastq_sorted_dedup_sort.bed", sep = "\t", header = FALSE)
S2_output <- S2_output[c(1:2)]
colnames(S2_output) <- c("label", "S2_output_read_count")
S2_output <- S2_output %>% group_by(label) %>% summarise_if(is.numeric, sum)

data_output <- merge(S1_output, S2_output, by = "label", all = TRUE, fill = TRUE)
data_output[is.na(data_output)] <- 0

#create methylation status and frag_grp
data_input$frag_grp <- substr(data_input$label, 1, 10)

data_input$methylation_status <- for (i in 1:nrow(data_input)) {
                                if (substr(data_input$label[i], (nchar(data_input$label)[i]+1)-4, nchar(data_input$label[i])) == "meth") {data_input[i,5] = "meth"}
                                else {data_input[i,5] = "unmeth"}}
colnames(data_input) <- c("label", "S1_input_read_count", "S2_input_read_count", "frag_grp", "methylation_status")

data_output$frag_grp <- substr(data_output$label, 1, 10)

data_output$methylation_status <- for (i in 1:nrow(data_output)) {
                                if (substr(data_output$label[i], (nchar(data_output$label)[i]+1)-4, nchar(data_output$label[i])) == "meth") {data_output[i,5] = "meth"}
                                else {data_output[i,5] = "unmeth"}}
colnames(data_output) <- c("label", "S1_output_read_count", "S2_output_read_count", "frag_grp", "methylation_status")

#Cpg fraction
data_input$frag_grp <- as.factor(data_input$frag_grp)
data_input$frag_grp <- revalue(data_input$frag_grp, c("320b_16C_3" = "320b_16C_35"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("320b_16C_5" = "320b_16C_50"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("320b_16C_6" = "320b_16C_65"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_1C_35G" = "80b_1C_35"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_1C_50G" = "80b_1C_50"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_1C_65G" = "80b_1C_65"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_2C_35G" = "80b_2C_35"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_2C_50G" = "80b_2C_50"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_2C_65G" = "80b_2C_65"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_4C_35G" = "80b_4C_35"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_4C_50G" = "80b_4C_50"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("80b_4C_65G" = "80b_4C_65"))

data_output$frag_grp <- as.factor(data_output$frag_grp)
data_output$frag_grp <- revalue(data_output$frag_grp, c("320b_16C_3" = "320b_16C_35"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("320b_16C_5" = "320b_16C_50"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("320b_16C_6" = "320b_16C_65"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_1C_35G" = "80b_1C_35"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_1C_50G" = "80b_1C_50"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_1C_65G" = "80b_1C_65"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_2C_35G" = "80b_2C_35"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_2C_50G" = "80b_2C_50"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_2C_65G" = "80b_2C_65"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_4C_35G" = "80b_4C_35"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_4C_50G" = "80b_4C_50"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("80b_4C_65G" = "80b_4C_65"))


##since initial 320bp, 65% G+C is not correct, let's correct that to correct G+C content
data_input$frag_grp <- revalue(data_input$frag_grp, c("320b_4C_65" = "320b_04C_35"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("320b_8C_65" = "320b_08C_50"))
data_input$frag_grp <- revalue(data_input$frag_grp, c("320b_16C_65" = "320b_016C_50"))

data_output$frag_grp <- revalue(data_output$frag_grp, c("320b_4C_65" = "320b_04C_35"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("320b_8C_65" = "320b_08C_50"))
data_output$frag_grp <- revalue(data_output$frag_grp, c("320b_16C_65" = "320b_016C_50"))

#Cpg fraction
data_input$CpG <- as.factor(substr(data_input$label, 4,7))
data_input$CpG <- revalue(data_input$CpG, c("_1C_" = "1"))
data_input$CpG <- revalue(data_input$CpG, c("_2C_" = "2"))
data_input$CpG <- revalue(data_input$CpG, c("_4C_" = "4"))
data_input$CpG <- revalue(data_input$CpG, c("b_16" = "16"))
data_input$CpG <- revalue(data_input$CpG, c("b_2C" = "2"))
data_input$CpG <- revalue(data_input$CpG, c("b_4C" = "4"))
data_input$CpG <- revalue(data_input$CpG, c("b_8C" = "8"))
data_input$CpG <- as.numeric(as.character(data_input$CpG))

data_output$CpG <- as.factor(substr(data_output$label, 4,7))
data_output$CpG <- revalue(data_output$CpG, c("_1C_" = "1"))
data_output$CpG <- revalue(data_output$CpG, c("_2C_" = "2"))
data_output$CpG <- revalue(data_output$CpG, c("_4C_" = "4"))
data_output$CpG <- revalue(data_output$CpG, c("b_16" = "16"))
data_output$CpG <- revalue(data_output$CpG, c("b_2C" = "2"))
data_output$CpG <- revalue(data_output$CpG, c("b_4C" = "4"))
data_output$CpG <- revalue(data_output$CpG, c("b_8C" = "8"))
data_output$CpG <- as.numeric(as.character(data_output$CpG))

data_input$fragment_len <- as.factor(substr(data_input$label, 1,3))
data_input$fragment_len <- revalue(data_input$fragment_len, c("80b" = "80"))


data_output$fragment_len <- as.factor(substr(data_output$label, 1,3))
data_output$fragment_len <- revalue(data_output$fragment_len, c("80b" = "80"))

data_input$cpg_frac <- as.factor(as.numeric(data_input$CpG) / as.numeric(as.character(data_input$fragment_len)))
data_output$cpg_frac <- as.factor(as.numeric(data_output$CpG) / as.numeric(as.character(data_output$fragment_len)))

data_input$cpg_frac <- set_cpg_frac(data_input)
data_output$cpg_frac <- set_cpg_frac(data_output)

##Add two states for CpG fractions to visualize the replicate fragments better
setDT(data_input)[frag_grp == "320b_04C_35" & cpg_frac == "1/80", cpg_frac := paste0(cpg_frac, "(2)")]
setDT(data_output)[frag_grp == "320b_04C_35" & cpg_frac == "1/80", cpg_frac := paste0(cpg_frac, "(2)")]

setDT(data_input)[frag_grp == "320b_08C_50" & cpg_frac == "1/40", cpg_frac := paste0(cpg_frac, "(2)")]
setDT(data_output)[frag_grp == "320b_08C_50" & cpg_frac == "1/40", cpg_frac := paste0(cpg_frac, "(2)")]

setDT(data_input)[frag_grp == "320b_016C_50" & cpg_frac == "1/20", cpg_frac := paste0(cpg_frac, "(2)")]
setDT(data_output)[frag_grp == "320b_016C_50" & cpg_frac == "1/20", cpg_frac := paste0(cpg_frac, "(2)")]

##Fix to the correct GC
data_input$GC <- as.factor(str_sub(data_input$frag_grp,-2,-1))
data_output$GC <- as.factor(str_sub(data_output$frag_grp,-2,-1))

##input samples on miseq
data_input$cpg_frac <- factor(data_input$cpg_frac, levels = c("1/80", "1/80(2)", "1/40", "1/40(2)", "1/20", "1/20(2)"))
data_input_cleaned <- data_wrang(data_input)

png("2019_input_miseq_plots.png", height = 10, width = 24, units = "in", res = 300)
plot80_input <- dist_plots(as.data.frame(data_input_cleaned[1]), as.data.frame(data_input_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000))
plot160_input <- dist_plots(as.data.frame(data_input_cleaned[3]), as.data.frame(data_input_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000))
plot320_input <- dist_plots(as.data.frame(data_input_cleaned[5]), as.data.frame(data_input_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000))
mylegend <- g_legend(plot80_input)
grid.arrange(arrangeGrob(plot80_input + theme(legend.position = "none"),
                         plot160_input + theme(legend.position = "none"), ncol = 2, nrow = 1),
                         arrangeGrob(plot320_input + theme(legend.position = "none"), ncol = 3.35,
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##output samples on miseq
data_output$cpg_frac <- factor(data_output$cpg_frac, levels = c("1/80", "1/80(2)", "1/40", "1/40(2)", "1/20", "1/20(2)"))
data_output_cleaned <- data_wrang(data_output)

png("2019_output_miseq_plots.png", height = 10, width = 24, units = "in", res = 300)
plot80_output <- dist_plots(as.data.frame(data_output_cleaned[1]), as.data.frame(data_output_cleaned[2])) +
  ggtitle("Fragment length: 80 bp") +
  theme(legend.position = "bottom", plot.title = element_text(size = 16, face = "bold")) +
 scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000))
plot160_output <- dist_plots(as.data.frame(data_output_cleaned[3]), as.data.frame(data_output_cleaned[4])) +
  ggtitle("Fragment length: 160 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000))
plot320_output <- dist_plots(as.data.frame(data_output_cleaned[5]), as.data.frame(data_output_cleaned[6])) +
  ggtitle("Fragment length: 320 bp") +
  theme(plot.title = element_text(size = 16, face = "bold")) +
  scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000))
mylegend <- g_legend(plot80_output)
grid.arrange(arrangeGrob(plot80_output + theme(legend.position = "none"),
                         plot160_output + theme(legend.position = "none"), ncol = 2, nrow = 1),
                         arrangeGrob(plot320_output + theme(legend.position = "none"), ncol = 3.35,
                         nrow = 1), mylegend, nrow = 2, heights = c(10, 1))
dev.off()

##Putting plots in order of what I want for publication
png("2020_frag_GC_CpG_pubplts.png", height = 25, width = 20, units = "in", res = 300)
mylegend <- g_legend(plot80_0.01)
plots_pub <- grid.arrange(arrangeGrob(plot80_input + theme(legend.position = "none", axis.title.y = element_text("Reads")) +
                                                scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000)),
                                        plot160_input + theme(axis.title.y = element_blank(), legend.position = "none") +
                                                scale_y_continuous(labels = comma, limits = c(0, 50000), c(0, 10000, 20000, 30000, 40000, 50000)), ncol =2, nrow =1),
                                        arrangeGrob(plot320_input + theme(axis.title.y = element_blank(), legend.position = "none") +
                                                scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000)), ncol = 2, nrow = 1),
                                        arrangeGrob(plot80_output + theme(legend.position = "none", axis.title.y = element_text("Reads")) +
                                                scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000)),
                                        plot160_output + theme(axis.title.y = element_blank(), legend.position = "none") +
                                                scale_y_continuous(labels = comma, limits = c(0, 50000), breaks = c(0, 10000, 20000, 30000, 40000, 50000)), ncol = 2, nrow = 1),
                                        arrangeGrob(plot320_output + theme(axis.title.y = element_blank(), legend.position = "none") +
                                                scale_y_continuous(labels = comma, limits = c(0, 50000), c(0, 10000, 20000, 30000, 40000, 50000)), ncol = 2, nrow =1),
                                        arrangeGrob(plot80_0.01 + theme(legend.position = "none", axis.title.y = element_text("Reads")) +
						scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 2500, 5000, 7500, 10000)),
                                        plot160_0.01 + theme(axis.title.y = element_blank(), legend.position = "none") +
						scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000)), ncol = 2, nrow = 1),
                                        arrangeGrob(plot320_0.01 + theme(axis.title.y = element_blank(), legend.position = "none") +
						 scale_y_continuous(labels = comma, limits = c(0, 20000), breaks = c(0, 5000, 10000, 15000, 20000)), ncol = 2, nrow = 1))
dev.off()








