#!/usr/bin/env Rscript

##Load needed packages

library("qsea")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
BSgenome= "BSgenome.Hsapiens.UCSC.hg38"
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
library(tidyvers)

###Heat-scree plots are adapted from:
#### De Souza, Rebecca AG, et al. "DNA methylation profiling in human Huntington's disease brain." Human molecular genetics 25.10 (2016): 2013-2030.

##ggplot2 theme 
source("/cluster/home/wilsons/commonly_used_code/ggplot2_theme_bw.R")

#Bam files for this project
PATH_PREFIX_B1 <- "/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/"
BamsList_B1 <- list.files(path=PATH_PREFIX_B1, pattern = "filtered_human.*dedup.bam")
bamFile_B1 <- data.frame(paste0(PATH_PREFIX_B1,BamsList_B1))
colnames(bamFile_B1) <- "file_name"

PATH_PREFIX_B2 <- "/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch2_JB/"
BamsList_B2 <- list.files(path=PATH_PREFIX_B2, pattern = "filtered_human.*dedup.bam")
bamFile_B2 <- data.frame(paste0(PATH_PREFIX_B2,BamsList_B2))
colnames(bamFile_B2) <- "file_name"

PATH_PREFIX_B3 <- "/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch3_DT/"
BamsList_B3 <- list.files(path=PATH_PREFIX_B3, pattern = "filtered_human.*dedup.bam")
bamFile_B3 <- data.frame(paste0(PATH_PREFIX_B3,BamsList_B3))
colnames(bamFile_B3) <- "file_name"

bamFile <- rbind(bamFile_B1, bamFile_B2, bamFile_B3)

#Create Qsea set for each sample and bind together in a list
sample_table = data.frame(sample_name = paste0("sample_", 1:15), file_name = bamFile$file_name, group = c(rep("B1",5), rep("B2", 5), rep ("B3", 5)), stringsAsFactors = FALSE) 
	
#Defining parameters again
BSgenome="BSgenome.Hsapiens.UCSC.hg38"
ws=300
chr.select=paste0("chr",1:22)
	
#Create Qsea set for each sample, this will be combined with the other sample MEDIPS set into a list
qseaSet = createQseaSet(sampleTable = sample_table, BSgenome = BSgenome, window_size = ws, chr.select = chr.select)
qseaSet = addCoverage(qseaSet, uniquePos = TRUE, paired = TRUE)
qseaSet = addPatternDensity(qseaSet, "CG", name = "CpG", fragment_length = 300)
qseaSet=addLibraryFactors(qseaSet)
qseaSet = addOffset(qseaSet, enrichmentPattern = "CpG")

str(qseaSet)#Check that SetCreated.list is infact a list of N (number of samples in data set) MEDIPS sets

#Specify the shortened sample names here
save(qseaSet, file = "2020_allbatches_QseaSets.RData")#save the MEDIPS sets

###Load in the pmol data and merge together
###Batch 1
S1_B1 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_human.aligned.unalignedtospike.trimmed.6654_S1_L002_R1_001.fastq.1unalignedtospike.trimmed.6654_S1_L002_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S1_B1) <- c("chr", "start", "end", "S1B1_read_count", "S1B1_pmol")
S1_B1$window <- paste0(S1_B1$chr,"_",S1_B1$start,"_",S1_B1$end)
S1_B1 <- S1_B1[-c(1:3)]
S1_B1 <- S1_B1 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S2_B1 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_human.aligned.unalignedtospike.trimmed.6655_S2_L002_R1_001.fastq.1unalignedtospike.trimmed.6655_S2_L002_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S2_B1) <- c("chr", "start", "end", "S2B1_read_count", "S2B1_pmol")
S2_B1$window <- paste0(S2_B1$chr,"_",S2_B1$start,"_",S2_B1$end)
S2_B1 <- S2_B1[-c(1:3)]
S2_B1 <- S2_B1 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S3_B1 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_human.aligned.unalignedtospike.trimmed.6656_S3_L002_R1_001.fastq.1unalignedtospike.trimmed.6656_S3_L002_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S3_B1) <- c("chr", "start", "end", "S3B1_read_count", "S3B1_pmol")
S3_B1$window <- paste0(S3_B1$chr,"_",S3_B1$start,"_",S3_B1$end)
S3_B1 <- S3_B1[-c(1:3)]
S3_B1 <- S3_B1 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S4_B1 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_human.aligned.unalignedtospike.trimmed.6657_S4_L002_R1_001.fastq.1unalignedtospike.trimmed.6657_S4_L002_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S4_B1) <- c("chr", "start", "end", "S4B1_read_count", "S4B1_pmol")
S4_B1$window <- paste0(S4_B1$chr,"_",S4_B1$start,"_",S4_B1$end)
S4_B1 <- S4_B1[-c(1:3)]
S4_B1 <- S4_B1 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S5_B1 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_human.aligned.unalignedtospike.trimmed.6658_S5_L002_R1_001.fastq.1unalignedtospike.trimmed.6658_S5_L002_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S5_B1) <- c("chr", "start", "end", "S5B1_read_count", "S5B1_pmol")
S5_B1$window <- paste0(S5_B1$chr,"_",S5_B1$start,"_",S5_B1$end)
S5_B1 <- S5_B1[-c(1:3)]
S5_B1 <- S5_B1 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

#merge samples together
batch1 <- merge(S1_B1, S2_B1, by = "window", all = TRUE, fill = TRUE)
batch1 <- merge(batch1, S3_B1, by = "window", all = TRUE, fill = TRUE)
batch1 <- merge(batch1, S4_B1,  by = "window", all = TRUE, fill = TRUE)
batch1 <- merge(batch1, S5_B1, by = "window", all = TRUE, fill = TRUE)

###Batch2
S1_B2 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch2_JB/filtered_human.aligned.unalignedtospike.trimmed.JB_A29_S11_L001_R1_001.fastq.1unalignedtospike.trimmed.JB_A29_S11_L001_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S1_B2) <- c("chr", "start", "end", "S1B2_read_count", "S1B2_pmol")
S1_B2$window <- paste0(S1_B2$chr,"_",S1_B2$start,"_",S1_B2$end)
S1_B2 <- S1_B2[-c(1:3)]
S1_B2 <- S1_B2 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S2_B2 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch2_JB/filtered_human.aligned.unalignedtospike.trimmed.JB_A35_S12_L001_R1_001.fastq.1unalignedtospike.trimmed.JB_A35_S12_L001_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S2_B2) <- c("chr", "start", "end", "S2B2_read_count", "S2B2_pmol")
S2_B2$window <- paste0(S2_B2$chr,"_",S2_B2$start,"_",S2_B2$end)
S2_B2 <- S2_B2[-c(1:3)]
S2_B2 <- S2_B2 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S3_B2 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch2_JB/filtered_human.aligned.unalignedtospike.trimmed.JB_A37_S13_L001_R1_001.fastq.1unalignedtospike.trimmed.JB_A37_S13_L001_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S3_B2) <- c("chr", "start", "end", "S3B2_read_count", "S3B2_pmol")
S3_B2$window <- paste0(S3_B2$chr,"_",S3_B2$start,"_",S3_B2$end)
S3_B2 <- S3_B2[-c(1:3)]
S3_B2 <- S3_B2 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S4_B2 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch2_JB/filtered_human.aligned.unalignedtospike.trimmed.JB_A56_S14_L001_R1_001.fastq.1unalignedtospike.trimmed.JB_A56_S14_L001_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S4_B2) <- c("chr", "start", "end", "S4B2_read_count", "S4B2_pmol")
S4_B2$window <- paste0(S4_B2$chr,"_",S4_B2$start,"_",S4_B2$end)
S4_B2 <- S4_B2[-c(1:3)]
S4_B2 <- S4_B2 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S5_B2 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch2_JB/filtered_human.aligned.unalignedtospike.trimmed.JB_A64_S15_L001_R1_001.fastq.1unalignedtospike.trimmed.JB_A64_S15_L001_R1_001.fastq.2_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S5_B2) <- c("chr", "start", "end", "S5B2_read_count", "S5B2_pmol")
S5_B2$window <- paste0(S5_B2$chr,"_",S5_B2$start,"_",S5_B2$end)
S5_B2 <- S5_B2[-c(1:3)]
S5_B2 <- S5_B2 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

#merge samples together
batch2 <- merge(S1_B2, S2_B2, by = "window", all = TRUE, fill = TRUE)
batch2 <- merge(batch2, S3_B2, by = "window", all = TRUE, fill = TRUE)
batch2 <- merge(batch2, S4_B2,  by = "window", all = TRUE, fill = TRUE)
batch2 <- merge(batch2, S5_B2, by = "window", all = TRUE, fill = TRUE)

##Batch3
S1_B3 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch3_DT/filtered_human.aligned.unalignedtospike.trimmed.SWID_16386952_TGL54_0001_Ct_T_PE_317_CM_A29_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S1_B3) <- c("chr", "start", "end", "S1B3_read_count", "S1B3_pmol")
S1_B3$window <- paste0(S1_B3$chr,"_",S1_B3$start,"_",S1_B3$end)
S1_B3 <- S1_B3[-c(1:3)]
S1_B3 <- S1_B3 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S2_B3 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch3_DT/filtered_human.aligned.unalignedtospike.trimmed.SWID_16386954_TGL54_0003_Ct_T_PE_316_CM_A37_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S2_B3) <- c("chr", "start", "end", "S2B3_read_count", "S2B3_pmol")
S2_B3$window <- paste0(S2_B3$chr,"_",S2_B3$start,"_",S2_B3$end)
S2_B3 <- S2_B3[-c(1:3)]
S2_B3 <- S2_B3 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S3_B3 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch3_DT/filtered_human.aligned.unalignedtospike.trimmed.SWID_16386956_TGL54_0002_Ct_T_PE_309_CM_A35_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S3_B3) <- c("chr", "start", "end", "S3B3_read_count", "S3B3_pmol")
S3_B3$window <- paste0(S3_B3$chr,"_",S3_B3$start,"_",S3_B3$end)
S3_B3 <- S3_B3[-c(1:3)]
S3_B3 <- S3_B3 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S4_B3 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch3_DT/filtered_human.aligned.unalignedtospike.trimmed.SWID_16386958_TGL54_0004_Ct_T_PE_310_CM_A56_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S4_B3) <- c("chr", "start", "end", "S4B3_read_count", "S4B3_pmol")
S4_B3$window <- paste0(S4_B3$chr,"_",S4_B3$start,"_",S4_B3$end)
S4_B3 <- S4_B3[-c(1:3)]
S4_B3 <- S4_B3 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

S5_B3 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch3_DT/filtered_human.aligned.unalignedtospike.trimmed.SWID_16386960_TGL54_0005_Ct_T_PE_309_CM_A64_sorted_dedup_sort_cut_CG_all_pmol_hg38_intersect_adjbinned.bed", header =T)
colnames(S5_B3) <- c("chr", "start", "end", "S5B3_read_count", "S5B3_pmol")
S5_B3$window <- paste0(S5_B3$chr,"_",S5_B3$start,"_",S5_B3$end)
S5_B3 <- S5_B3[-c(1:3)]
S5_B3 <- S5_B3 %>%  group_by(window) %>% summarise_if(is.numeric, sum)

#merge samples together
batch3 <- merge(S1_B3, S2_B3, by = "window", all = TRUE, fill = TRUE)
batch3 <- merge(batch3, S3_B3, by = "window", all = TRUE, fill = TRUE)
batch3 <- merge(batch3, S4_B3,  by = "window", all = TRUE, fill = TRUE)
batch3 <- merge(batch3, S5_B3, by = "window", all = TRUE, fill = TRUE)

#merge batches together
data <- merge(batch1, batch2, by = "window", all = TRUE, fill = TRUE)
data <- merge(data, batch3, by = "window", all = TRUE, fill = TRUE)
data[is.na(data)] <- 0
data$window <- as.factor(data$window)

data_raw <- data[-c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31)]
data_pmol <- data[-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)]
data_pmol_raw <- data[-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)]

#Missing <- lapply(data_pmol_raw, function(x){ length(which(x==0))/length(x)})

##Remove repetitive elements- simple repeats only 
repeats <- read.table("~/Annotations/2020_repeatregions_mappedto300bpwindows.bed", sep = "\t", header = FALSE)
###subset to regions where there are overlaps between our windows
repeats <- subset(repeats, repeats$V9>0)
repeats <- repeats[,c(1:3)]
colnames(repeats) <- c("chr","start","end")
repeats$window <- paste0(repeats$chr,"_",repeats$start,"_",repeats$end)
repeats <- unique(repeats)
rownames(repeats) <-repeats$window

##remove repetitive regions from the fmol data
rownames(data_pmol) <- data_pmol$window
data_pmol <- data_pmol[!rownames(data_pmol) %in% rownames(repeats),]

##remove blacklist regions
blacklist <- read.table("~/Annotations/2020_new_ENCODE_blacklist.csv", sep = "\t", header = FALSE)
colnames(blacklist) <- c("bl_chr", "bl_start", "bl_end", "chr", "start", "end", "overlap")
blacklist$window <- paste0(blacklist$chr,"_",blacklist$start,"_",blacklist$end)
blacklist <- unique(blacklist)

##remove blacklist regions from fmol
data_pmol <- data_pmol[!rownames(data_pmol) %in% blacklist$window,]

## remove low mappability regions (mappability score<0.5)
mappability <- read.csv("~/Annotations/2020_min_mappability_300bp_windows.csv", header = TRUE)
mappability$X <- NULL

mappability_0.5 <- subset(mappability, mappability$map_score<0.5)

data_pmol <- data_pmol[!rownames(data_pmol) %in% mappability_0.5$window,]

##Infer sex from Y signal in raw data
chrY_reads <- subset(data_raw, str_detect(data_raw$window, pattern = "chrY") == "TRUE")
chrY_reads[is.na(chrY_reads)] <- 0
no_chrY <- lapply(chrY_reads, function(x){ length(which(x==0))/length(x)})

###Heat scree effect size
source("/cluster/home/wilsons/commonly_used_code/heatscree_cohensd_plot.R")

Loadings_meta <- data.frame(batch = as.factor(c(rep("Batch_1", 5), rep("Batch_2", 5), rep("Batch_3", 5))),
                            sample = rep(c("sample_1", "sample_2", "sample_3", "sample_4", "sample_5"),3),
                            sequencer = as.factor(c(rep("A",5), rep("A", 5), rep("B", 5))),
                            adapters = as.factor(c(rep("A",5), rep("B",5), rep("A",5))),
                            sex = as.factor(rep(c("F","F","M","M","F"),3)))
                            #NPM1 = as.factor(rep(c("pos", "pos", "pos", "neg", "neg"),3)),
                            #FLT3_ITD = as.factor(rep(c("pos", "neg", "neg", "neg", "neg"),3)),
                            #IDH1 = as.factor(rep(c("neg", "pos", "pos", "neg", "neg"),3)),
                            #IDH2 = as.factor(rep(c("pos", "neg", "neg", "neg", "neg"),3)),
                            #TET2 = as.factor(rep(c("neg", "neg", "neg", "pos", "pos"),3)))

# Specifiy the number of PCs you want shown
Num<-9 # should be equal to the number of samples in your dataset; for large datasets, you can opt to just see the top PCs
# Designate what order you want the variables to appear (continuous variables rbinded to categorical variables in function)
Order<-c(1,2,3,4,5)

meta_categorical <- Loadings_meta[,c("batch","sequencer","adapters","sample", "sex")]  # input column numbers in meta that contain categorical variables
colnames(meta_categorical) <- c("Batch","Unmethylated filler", "Adapters", "Sample", "Sex")
meta_categorical$Sample <- as.factor(meta_categorical$Sample)

##QSEA
pca <- prcomp(t(qseaSet@count_matrix))
vars <- pca$sdev^2
Importance <- vars/sum(vars)

Loadings <- as.data.frame(pca$x)

##QSEA plot
png("2020_heatscree_meand_qsea.png", height = 6, width = 12 , unit = "in", res =300)
qsea <- heat_scree_plot_es(Loadings, Importance, Num, Order) + theme(legend.position = "none")
qsea
dev.off()

Loadings_qsea <- Loadings

##pmol filtered
pca_dat <- data_pmol
pca_dat[is.na(pca_dat)] <- 0
rownames(pca_dat) <- pca_dat$window
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev^2
Importance <- vars/sum(vars) ##3% of variance associated with batch

Loadings <- as.data.frame(pca$x)

##pmol filtered plot
png("2020_heatscree_meand_pmolfiltered.png", height = 12, width =12 , unit = "in", res =300)
pmol <- heat_scree_plot_es(Loadings, Importance, Num, Order) + theme(legend.position = "none")
pmol
dev.off()

Loadings_pmol <- Loadings

##Raw data
pca_dat <- data_raw
pca_dat[is.na(pca_dat)] <- 0
rownames(pca_dat) <- pca_dat$window
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev^2
Importance <- vars/sum(vars)

Loadings <- as.data.frame(pca$x)

png("2020_heatscree_meand_raw.png", height = 12, width =12 , unit = "in", res =300)
raw <- heat_scree_plot_es(Loadings, Importance, Num, Order) + theme(legend.position = "none")
raw
dev.off()

Loadings_raw <- Loadings

##pmol not filtered
pca_dat <- data_pmol_raw
pca_dat[is.na(pca_dat)] <- 0
rownames(pca_dat) <- pca_dat$window
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev^2
Importance <- vars/sum(vars)

Loadings <- as.data.frame(pca$x)

png("2020_heatscree_meand_allwindows_fmol.png", height = 12, width =12 , unit = "in", res =300)
pmol_raw <- heat_scree_plot_es(Loadings, Importance, Num, Order) + theme(legend.position = "none")
pmol_raw
dev.off()

Loadings_pmol_raw <- Loadings

png("2020_Heatscree_allbatches_meand_pmolandpmolnoblacklist.png", height = 24, width = 12, units = "in", res =300)
grid.arrange(arrangeGrob(raw + theme(legend.position = "none", axis.title.x = element_text("")),
        arrangeGrob(qsea + theme(legend.position = "none", axis.title.x = element_text(""))),
        arrangeGrob(pmol_raw  + theme(legend.position = "none", axis.title.x = element_text(""))),
        arrangeGrob(pmol  + theme(legend.position = "none", axis.title.x = element_text(""))), ncol = 1))
dev.off()


### How much of the batch effects aka PC5 are driven by repeat regions?
repeat_masker <- read.table("~/2020_RepeatMasker_300bp.bed", sep = "\t", header = FALSE)
colnames(repeat_masker) <- c("repeat_chr", "repeat_start", "repeat_end", "element", "number", "strand", "chr", "start", "end", "overlap")
repeat_masker$window <- paste0(repeat_masker$chr,"_", repeat_masker$start,"_", repeat_masker$end)

pca_dat <- data_pmol
pca_dat[is.na(pca_dat)] <- 0
rownames(pca_dat) <- pca_dat$window
pca_dat$window <- NULL
pca <- prcomp(t(pca_dat))
vars <- pca$sdev^2
Importance <- vars/sum(vars)

##PC3
rotations <- as.data.frame(pca$rotation)
rotations_PC3 <- as.data.frame(rotations$PC3)
rownames(rotations_PC3) <- rownames(rotations)
rotations_PC3 <- abs(rotations_PC3)
colnames(rotations_PC3) <- c("PC3")
rotations_PC3$window <- rownames(rotations_PC3)
rotations_PC3 <- as.data.frame(rotations_PC3[order(-rotations_PC3$PC3),])

PC3_top <- as.data.frame(rotations_PC3[1:(round(0.1*nrow(rotations_PC3))),])

##how many of top 10% of drivers of PC5 are repeat regions?
sum(rownames(PC3_top) %in% repeat_masker$window)/(nrow(PC3_top)) #71%
driving_repeats <- repeat_masker[repeat_masker$window %in% rownames(PC3_top),]
sum(driving_repeats$element == "Alu*")
nrow(driving_repeats[grep("Alu*", driving_repeats$element),])/ nrow(driving_repeats) ##0.42 are Alu
summary(driving_repeats$element)

##PC5
rotations <- as.data.frame(pca$rotation)
rotations_PC5 <- as.data.frame(rotations$PC5)
rownames(rotations_PC5) <- rownames(rotations)
rotations_PC5 <- abs(rotations_PC5)
colnames(rotations_PC5) <- c("PC5")
rotations_PC5$window <- rownames(rotations_PC5)
rotations_PC5 <- as.data.frame(rotations_PC5[order(-rotations_PC5$PC5),])

PC5_top <- as.data.frame(rotations_PC5[1:(round(0.1*nrow(rotations_PC5))),])

##how many of top 10% of drivers of PC5 are repeat regions?
sum(rownames(PC5_top) %in% repeat_masker$window)/(nrow(PC5_top)) #0.7227874, so 72% of drivers in batch are in repetitive regions
driving_repeats <- repeat_masker[repeat_masker$window %in% rownames(PC5_top),]
sum(driving_repeats$element == "Alu*")
nrow(driving_repeats[grep("Alu*", driving_repeats$element),])/ nrow(driving_repeats) ##0.43 are Alu

