#!/usr/bin/env R

#Load libraries
library(rmarkdown)
library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(boot)
library(BlandAltmanLeh)
library(scales)
library(stringr)
library(forcats)
library(varhandle)
library(purrr)
library(stringi)

##Load data
###For now I'm reading in 5  samples, in future need to figure out the multmerge that is in the methylationvaladjust script
S6654 <- read.csv("filtered_aligned.trimmed.6654_S1_L002_R1_001.fastqtrimmed.6654_S1_L002_R2_001.fastq.bam_Cleaned_dedup_readcounts.csv", header = TRUE)
S6654$S6654_read_count <- 1
table(S6654$methylation_status) ##checking binding specificity
S6654 <- subset(S6654, S6654$methylation_status == "meth")

S6655 <- read.csv("filtered_aligned.trimmed.6655_S2_L002_R1_001.fastqtrimmed.6655_S2_L002_R2_001.fastq.bam_Cleaned_dedup_readcounts.csv", header = TRUE)
S6655$S6655_read_count <- 1
table(S6655$methylation_status) ##checking binding specificity
S6655 <- subset(S6655, S6655$methylation_status == "meth")

S6656 <- read.csv("filtered_aligned.trimmed.6656_S3_L002_R1_001.fastqtrimmed.6656_S3_L002_R2_001.fastq.bam_Cleaned_dedup_readcounts.csv", header = TRUE)
S6656$S6656_read_count <- 1
table(S6656$methylation_status) ##checking binding specificity
S6656 <- subset(S6656, S6656$methylation_status == "meth")

S6657 <- read.csv("filtered_aligned.trimmed.6657_S4_L002_R1_001.fastqtrimmed.6657_S4_L002_R2_001.fastq.bam_Cleaned_dedup_readcounts.csv", header = TRUE)
S6657$S6657_read_count <- 1
table(S6657$methylation_status) ##checking binding specificity
S6657 <- subset(S6657, S6657$methylation_status == "meth")

S6658 <- read.csv("filtered_aligned.trimmed.6658_S5_L002_R1_001.fastqtrimmed.6658_S5_L002_R2_001.fastq.bam_Cleaned_dedup_readcounts.csv", header = TRUE)
S6658$S6658_read_count <- 1
table(S6658$methylation_status) ##checking binding specificity
S6658 <- subset(S6658, S6658$methylation_status == "meth")

data <- merge(S6654, S6655, by = "UMI", all = TRUE)
data <- data[,c(1:8,15)]
data <- merge(data, S6656, by = "UMI", all = TRUE)
data <- data[,c(1:9, 16)]
data <- merge(data, S6657, by = "UMI", all = TRUE)
data <- data[,c(1:10,17)]
data <- merge(data, S6658, by = "UMI", all = TRUE)
data <- data[,c(1, 3:11,18)]
colnames(data) <- c("UMI", "fragment_len", "CpG", "GC", "methylation_status", "frag_grp", "S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count")    

##fill in NAs from UMI data
data$fragment_len <- as.factor(str_sub(data$UMI, 0, 3))
data$fragment_len <- revalue(data$fragment_len, c("80b" = "80"))

data$CpG <- as.factor(str_sub(data$UMI, 5, 7))
data$CpG <- revalue(data$CpG, c("_16" = "16","_2C" = "2","_4C" = "4","_8C" = "8","1C_" = "1","2C_" = "2","4C_" = "4"))

data$GC <- as.factor(str_sub(data$UMI, 7, 11))
data$GC <- revalue(data$GC, c("_35G-" = "35", "_50G-" = "50", "_65G-" = "65", "6C_35" = "35", "6C_50" = "50", "6C_65"="65", "C_35G"="35", "C_50G" = "50", "C_65G" = "65"))

#data$methylation_status <- str_extract(data$UMI_barcode, "meth")
#data$methylation_status <- data$methylation_status %>% fct_explicit_na(na_level = 'unmeth')
##remove unmethylated
#data <- subset(data, data$methylation_status == "meth")

data$frag_grp <- as.factor(paste0(data$fragment_len, sep = '_', data$CpG,sep = '_', data$GC))

data$fragment_len <- as.factor(str_sub(data$UMI, 0, 3))
data$fragment_len <- revalue(data$fragment_len, c("80b" = "80"))

data$CpG <- as.factor(str_sub(data$UMI, 5, 7))
data$CpG <- revalue(data$CpG, c("_16" = "16","_2C" = "2","_4C" = "4","_8C" = "8","1C_" = "1","2C_" = "2","4C_" = "4"))

data$GC <- as.factor(str_sub(data$UMI, 7, 11))
data$GC <- revalue(data$GC, c("_35G-" = "35", "_50G-" = "50", "_65G-" = "65", "6C_35" = "35", "6C_50" = "50", "6C_65"="65", "C_35G"="35", "C_50G" = "50", "C_65G" = "65"))

data$frag_grp <- as.factor(paste0(data$fragment_len, sep = '_', data$CpG,sep = '_', data$GC))
data$methylation_status <- NULL

##Correcting the fragment names in the 320bp fragments
data$frag_grp <- revalue(data$frag_grp, c("320_4_65" = "320_04_35"))
data$frag_grp <- revalue(data$frag_grp, c("320_8_65" = "320_08_50"))
data$frag_grp <- revalue(data$frag_grp, c("320_16_65" = "320_016_50"))

##Fix to the correct GC
data$GC <- as.factor(str_sub(data$frag_grp,-2,-1))

data[is.na(data)] <- 0
data <- data %>%
  group_by(frag_grp, fragment_len, GC, CpG) %>%
  summarise(read_count_6547 = sum(read_count_6547),
            read_count_6548 = sum(read_count_6548))

data[is.na(data)] <- 0

data <- ddply(data,.(data$frag_grp, data$fragment_len, data$GC, data$CpG),numcolwise(sum), .progress = "text")
colnames(data) <- c("frag_grp", "fragment_len", "GC", "CpG", "S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count")  

#data <- data %>%
 # group_by(frag_grp, fragment_len, GC, CpG) %>%
  #summarise_if(is.numeric(funs(sum))

##Add in concentration information
for (i in 1:nrow(data)){
  data[i,10] <- ifelse (data[i,2] == "160", 0.002, ifelse(data[i,2] == "80", 0.004, 0.001))
}

names(data)[10] <- "conc"

##Adjusting for the 0.01ng dilution
data$conc <- data$conc*0.9


##Some exploratory plots
###Melt so that data frame 
data_melt <- melt(data, id.vars = c("conc", "frag_grp", "GC", "CpG", "fragment_len"), measure.vars = c( "S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

png("2020_concvsreads_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data_melt, aes(x = conc, y = value)) + geom_point() +
  theme_bw()+
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Concentration (fmol/ng)")+
  ylab("Read counts")+
  scale_y_continuous(labels = comma)+
  scale_fill_manual("slategray3")
dev.off()

##Colour points by GC content
data$GC <- as.factor(data$GC)

png("2020_concvsreads_dedup_GCcol_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data_melt, aes(x = conc, y = value, col = GC)) + geom_point() +
  theme_bw()+
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Concentration (fmol/ng)")+
  ylab("Read counts")+
  scale_y_continuous(labels = comma)
dev.off()

data$GC <- as.numeric(as.character(data$GC))

##Correlation between read count and concentration
cor(data_melt$conc, data_melt$value)#-0.1993395

##Looking at distribution of our covariates
###Fragment length
data_melt$fragment_len <- as.numeric(as.character(data_melt$fragment_len))

png("2020_fraglen_dist_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = fragment_len, y = value, width = 20)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Fragment length (bp)") +
  ylab("Read counts") +
  scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()

###GC content
data_melt$GC <- as.numeric(as.character(data_melt$GC))

png("2020_GC_dist_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = GC, y = value, width = 4)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12),axis.text.y = element_text(vjust = 0.5, size = 12)) +
  xlab("G + C content") +
  ylab("Read counts") +
  scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()

###CpG number
png("2020_CpG_dist_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = CpG, y = value, width = 0.5)) +
  geom_bar(stat = "identity", fill = "peachpuff") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12)) +
  xlab("Number of CpGs") +
  ylab("Read counts") +
  scale_y_continuous(labels = comma,limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()

####All show non-normality, though I only expect fragment_len to have the monotonic issue

##Data transformation
###Fragement length
###Min-Max normalization
data$fragment_len <- as.numeric(data$fragment_len)

for (i in 1:nrow(data)){
  data[i,11] = (data[i,2] - min(data$fragment_len)) / ((max(data$fragment_len) - min(data$fragment_len)))
}

names(data)[11] = "frag_len_mmnorm"

##z-score normalization
for (i in 1:nrow(data)){
  data[i,12] = (data[i,2] - mean(data$fragment_len)) / sd(data$fragment_len)
}

names(data)[12] = "frag_len_z"

##Cube root
data$frag_len_3 <- (data$fragment_len) ^ (1/3)

##Plots of fragment length transformed values
##Min-Max normalization
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "frag_len_mmnorm", "CpG", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "frag_len_mmnorm", "CpG", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$frag_len_mmnorm <- as.factor(data_melt$frag_len_mmnorm)

png("2020_fraglen_maxminnorm_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = frag_len_mmnorm, y = value, width = .1)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12),axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Fragment length (bp)") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
##left skew

##Z-score
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "frag_len_z", "CpG", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "frag_len_z", "CpG", "GC", "frag_grp"), measure.vars = c( "S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$frag_len_z <- as.factor(data_melt$frag_len_z)

png("2020_fraglen_zscore_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = frag_len_z, y = value, width = .2)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Fragment length (bp)") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma,limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#left skew and negative vals

##Cube root
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "frag_len_3", "CpG", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "frag_len_3", "CpG", "GC", "frag_grp"), measure.vars = c( "S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$frag_len_3 <-as.factor(data_melt$frag_len_3)

png("2020_fraglen_cubert_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = frag_len_3, y = value, width = .1)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle=45, vjust=0.5, size=12),axis.text.y=element_text(vjust=0.5,size=12))+
  xlab("Fragment length (bp)")+
  ylab("Read counts")
  #scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
##left skew

####################################################################

##Fragment length log
data$fragment_len <- as.numeric(as.character(data$fragment_len))
data$fragment_len_log <- (160 - data$fragment_len)^2

#data_sub <- data[, c("frag_grp",  "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len_log", "CpG", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len_log", "CpG", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

png("2020_fraglenlog_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = fragment_len_log, y = value, width = 20)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12),axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Fragment length (bp)") +
  ylab("Read counts") +
  scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#right skew

data$fragment_len_log <- as.numeric(as.character(data$fragment_len_log))
data$fragment_len_log2 <- (data$fragment_len_log)^2
data$fragment_len_log3 <- (data$fragment_len_log)^(1/3)
data$fragment_len_log_rec <- 1/data$fragment_len_log
data$fragment_len_log_log <- log(data$fragment_len_log)

##Square
#data_sub <- data[, c("frag_grp",  "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len_log2", "CpG", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len_log2", "CpG", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

png("2020_fraglenlog2_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = fragment_len_log2, y = value, width = 1000000)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12),axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Fragment length (bp)") +
  ylab("Read counts") +
  scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
##right skew

##cube root
#data_sub <- data[, c("frag_grp",  "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len_log3", "CpG", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len_log3", "CpG", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$fragment_len_log3 <- as.factor(data_melt$fragment_len_log3)

png("2020_fraglenlog3_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = fragment_len_log3, y = value, width = .1)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12),axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Fragment length (bp)") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
##still right skew

##recipricol
#data_sub <- data[, c("frag_grp",  "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len_log_rec", "CpG", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len_log_rec", "CpG", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$fragment_len_log_rec <- as.factor(data_melt$fragment_len_log_rec)

png("2020_fraglenlog_rec_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = fragment_len_log_rec, y = value)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12),axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Fragment length (bp)") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
##left skew

##log
#data_sub <- data[, c("frag_grp",  "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len_log_log", "CpG", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len_log_log", "CpG", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$fragment_len_log_log <- as.factor(data_melt$fragment_len_log_log)

png("2020_fraglenloglog_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = fragment_len_log_log, y = value, width = .1)) +
  geom_bar(stat = "identity", fill = "slategray3") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12),axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Fragment length (bp)") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
##still right skew

##Looks like no data transformation will be best
##GC content transformations- DIDN'T redo in high res, as I didn't use any of these.
####This is negatively skewed, squared, cube root and log transformations are known to fix this. Let's try

data$GC <- as.numeric(as.character(data$GC))
data$GC_sq <- (data$GC) ^ 2
data$GC_3 <- (data$GC) ^ (1/3)
data$GC_log <- log(data$GC)
data$GC_rec <- 1 / data$GC

##Plot the distributions
##Square
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len", "CpG", "GC_sq", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len", "CpG", "GC_sq", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$GC_sq <-as.factor(data_melt$GC_sq)

png("2020_GC_sq_dedup_batch1.png", height = 6, width = 6, res = 300)
p <- ggplot(data = data_melt, aes(x = GC_sq, y = value)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("G + C content") +
  ylab("Read counts")
  #scale_y_continuous(labels = comma,limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
show(p)
dev.off()
#right skew

##cube root
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len", "CpG", "GC_3", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len", "CpG", "GC_3", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

png("2020_GC3_dedup_batch1.png", height = 6, width = 6, res = 300)
ggplot(data = data_melt, aes(x = GC_3, y = value)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("G + C content") +
  ylab("Read counts")
  #scale_y_continuous(labels = comma,limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#right skew

##log
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len", "CpG", "GC_log", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len", "CpG", "GC_log", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

png("2020_GC_log_dedup_batch1.png", height = 6, width = 6, res = 300)
ggplot(data = data_melt, aes(x = GC_log, y = value)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("G + C content") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma,limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#right skew

##reciprical
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len", "CpG", "GC_rec", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len", "CpG", "GC_rec", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

png("2020_GC_rec_dedup_batch1.png", height = 6, width = 6, res = 300)
ggplot(data = data_melt, aes(x = GC_rec, y = value)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("G + C content") +
  ylab("Read counts")
  #scale_y_continuous(labels = comma,limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#left skew
##Not using any transformations here either

#CpG transformations
####I see a left skew- let's look at the transformations to fix that and what returns the data back to normality
###Min-Max normalization
data$CpG <- as.numeric(as.character(data$CpG))

for (i in 1:nrow(data)){
  data[i,23] = (data[i,4] - min(data$CpG)) / ((max(data$CpG) - min(data$CpG)))
}

names(data)[23] <- "CpG_mmnorm"

##z-score normalization
for (i in 1:nrow(data)) {
  data[i,24] = (data[i,4] - mean(data$CpG)) / sd(data$CpG)
}

names(data)[24] <- "CpG_z"

##Cube root
data$CpG <- as.integer(data$CpG)
data$CpG_3 <- (data$CpG) ^(1/3)

##Plots of CpG transformed values
##Min-Max normalization
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len", "CpG_mmnorm", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len", "CpG_mmnorm", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$CpG_mmnorm <- as.factor(data_melt$CpG_mmnorm)

png("2020_CpG_maxminnorm_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = CpG_mmnorm, y = value, width = 0.1)) +
  geom_bar(stat = "identity", fill = "peachpuff") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12)) +
  xlab("Number of CpGs") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#better, sligh left skew

##zscore plot
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len", "CpG_z", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len", "CpG_z", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$CpG_z <- as.factor(data_melt$CpG_z)

png("2020_CpG_zscore_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = CpG_z, y = value, width = 0.1)) +
  geom_bar(stat = "identity", fill = "peachpuff") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12)) +
  xlab("Number of CpGs") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#bimodel and negative

##cube root
#data_sub <- data[, c("frag_grp", "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len", "CpG_3", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len", "CpG_3", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))
data_melt$CpG_3 <- as.factor(data_melt$CpG_3)

png("2020_CpG_3_dedup_batch1.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = CpG_3, y = value, width = 0.1)) +
  geom_bar(stat = "identity", fill = "peachpuff") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12)) +
  xlab("Number of CpGs") +
  ylab("Read counts") 
  #scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#better, this will be used

###Let's try out some models now
data$fragment_len <- as.numeric(as.character(data$fragment_len))
data$fragment_len <- as.factor(str_sub(data$frag_grp, 0, 3))
data$fragment_len <- revalue(data$fragment_len, c("80_" = "80"))
data$fragment_len <-as.numeric(as.character(data$fragment_len))

data$GC <- as.numeric(as.character(data$GC))

#data_sub <- data[, c("frag_grp",  "S11_00041891_read_count", "S12_00019563_read_count", "OB1080_read_count", "fragment_len", "CpG_3", "GC", "conc")]
data_melt <- melt(data, id.vars = c("conc", "fragment_len", "CpG_3", "GC", "frag_grp"), measure.vars = c("S6654_read_count", "S6655_read_count", "S6656_read_count", "S6657_read_count", "S6658_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

##Poisson
###Poisson- Getting negatives
poisson_glm <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = data_melt, family = poisson(link = "log"))
summary(poisson_glm)#In dpois(y, mu, log = TRUE) : non-integer x = 0.000842

png("2020_GLM_fraglen_CpG3_batch1diag.png")
glm.diag.plots(poisson_glm, glmdiag = glm.diag(poisson_glm))
dev.off()

png("2020_GLM_fraglen_CpG3_poisson_batch1.png")
plot(poisson_glm)
dev.off()

png("2020_BA_poisson_fraglen_CpG3_batch1.png")
bland.altman.plot(poisson_glm$data[,1], poisson_glm$resid, conf.int=.95, pch=19) ##Why are there negative values here
dev.off()

poisson_glm$data[,8] <- as.factor(as.character(poisson_glm$data[, 2]))
names(poisson_glm$data[,8]) <- "Fragment_len"

png("2020_BA_poisson_fraglen_CpG3_batch1.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(poisson_glm$data[,1], poisson_glm$resid, conf.int = .95, col.points = poisson_glm$data[,2], pch = 19, graph.sys = "ggplot2")
  
BA + aes(color = poisson_glm$data[,8]) + 
  theme_bw()+
  scale_color_manual(values = c("black", "darkgrey", "navy"))+
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size = 14), 
        axis.title.y = element_text(size = 20), axis.text.x  = element_text(size = 16),
        axis.text.y = element_text(vjust = 0.5, size = 16),axis.title.x = element_text(size = 20)) +
  xlab("Mean of measurements") +
  ylab("Difference")
dev.off()

##Calculate R2
r2_poisson <- 1 - (poisson_glm$deviance / poisson_glm$null.deviance)
r2_poisson #0.9781715

###Gaussian
gaussian_glm <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = data_melt, family = gaussian)
summary(gaussian_glm)

png("2020_GLM_gaussian_fraglen_CpG3_batch1.png")
glm.diag.plots(gaussian_glm, glmdiag = glm.diag(gaussian_glm))
dev.off()

png("2020_GLM_fraglen_CpG3_gaussianbatch1.png")
plot(gaussian_glm)
dev.off()

png("2020_BA_gaussian_fraglen_CpG3batch1.png")
bland.altman.plot(gaussian_glm$data[,1], gaussian_glm$resid, conf.int=.95, pch=19)
dev.off()

gaussian_glm$data[,8] <- as.factor(as.character(gaussian_glm$data[, 2]))
names(gaussian_glm$data[,8]) <- "Fragment_len"

png("2020_BA_gaussian_fraglen_CpG3_batch1.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(gaussian_glm$data[,1], gaussian_glm$resid, conf.int = .95, col.points = gaussian_glm$data[,2], pch = 19, graph.sys = "ggplot2")

BA + aes(color = gaussian_glm$data[,8]) + 
  theme_bw()+
  scale_color_manual(values = c("black", "darkgrey", "navy"))+
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size = 14), 
        axis.title.y = element_text(size = 20), axis.text.x  = element_text(size = 16),
        axis.text.y = element_text(vjust = 0.5, size = 16), axis.title.x = element_text(size = 20)) +
  xlab("Mean of measurements") +
  ylab("Difference")
dev.off()
##This estimates the 160bp very nicely

##Calculate R2
r2_gaussian= 1-(gaussian_glm$deviance/gaussian_glm$null.deviance)
r2_gaussian #  0.9173609
##Correlation isn't great, but makes sense in the context of wanting 160bp to perform better, best AIC I've seen
save(gaussian_glm, file = "2020_Gaussian_batch1.rda")

###Gaussian log
gaussian_log <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = data_melt, family = gaussian(link = "log"))
summary(gaussian_log)

png("2020_GLM_gaussianlog_fraglen_CpG3_diag_batch1.png")
glm.diag.plots(gaussian_log, glmdiag = glm.diag(gaussian_log))
dev.off()

png("2020_GLM_gaussianlog_fraglen_CpG3_batch1.png")
plot(gaussian_log)
dev.off()

##Bland-Altman plot (plots the residuals against the truth and looks at the deviance)
#gaussian_log$data[,2]=as.factor(as.character(gaussian_log$data[,2]))
gaussian_log$data[,8] <- as.factor(as.character(gaussian_log$data[,2]))
names(gaussian_log$data[,8]) <- "Fragment_len"

png("2020_BA_guassianlog_fraglen_CpG3_batch1.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(exp(gaussian_log$data[,1]), gaussian_log$resid, conf.int = .95, col.points = gaussian_log$data[,2], pch = 19, graph.sys="ggplot2")
  BA + aes(color=gaussian_log$data[,8]) + 
  theme_bw() +
  scale_color_manual(values = c("black", "darkgrey", "navy"))+
  theme(legend.title = element_blank(),legend.position = "right", legend.text = element_text(size = 14),
        axis.title.y = element_text(size = 20), axis.text.x  = element_text(size = 16),
        axis.text.y = element_text(vjust = 0.5, size = 16), axis.title.x = element_text(size = 20)) +
  xlab("Mean of measurements") +
  ylab("Difference")
dev.off()

##Calculate R2
r2_gaussian_log = 1-(gaussian_log$deviance/gaussian_log$null.deviance)
r2_gaussian_log #0.9815552

##Seeing if a simple lm would perform better
lm <- lm(conc ~ value + fragment_len + GC + CpG_3, data = data_melt)
summary(lm)#Residual standard error: 0.0003349 on 125 degrees of freedom
#Multiple R-squared:  0.9174,	Adjusted R-squared:  0.9147 

png("2020_lm_fraglen_CpG3_batch1.png")
plot(lm)
dev.off()

png("2020_lm_fraglen_CpG3_diag_batch1.png")
glm.diag.plots(lm, glmdiag = glm.diag(lm))
dev.off()

png("2020_BA_lm_fraglen_CpG3_batch1.png")
bland.altman.plot(lm$model[,1], lm$resid, conf.int=.95, pch=19)
dev.off()

lm$model[,6] <- as.factor(as.character(lm$model[,3]))
names(lm$model[,6]) <- "Fragment_len"

png("2020_BA_lm_fraglen_CpG3_batch1.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(lm$model[,6], gaussian_log$resid, conf.int = .95, col.points = lm$model[,3], pch = 19, graph.sys="ggplot2")
BA + aes(color=lm$model[,6]) + 
  theme_bw() +
  scale_color_manual(values = c("black", "darkgrey", "navy"))+
  theme(legend.title = element_blank(),legend.position = "right", legend.text = element_text(size = 14),
        axis.title.y = element_text(size = 20), axis.text.x  = element_text(size = 16),
        axis.text.y = element_text(vjust = 0.5, size = 16), axis.title.x = element_text(size = 20)) +
  xlab("Mean of measurements") +
  ylab("Difference")
dev.off()

#save(lm, file = "2020_lm_0.01ng.rda")
