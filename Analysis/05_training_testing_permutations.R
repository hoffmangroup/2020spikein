#!/usr/bin/env R

#Load libraries
library(rmarkdown)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(boot)
library(BlandAltmanLeh)
library(scales)
library(stringr)
library(forcats)
library(ggpubr)
library(caret)
library(patchwork)
library(Metrics)
source("~/commonly_used_code/ggplot2_theme_bw.R")

#Load data
###Sample 1
sample1 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_aligned.trimmed.6654_S1_L002_R1_001.fastqtrimmed.6654_S1_L002_R2_001.fastq_sorted_dedup_sort.bed")
sample1 <- as.data.frame(sample1[,1])
colnames(sample1) <- c("frag_grp")

#Assess methylation specificity and subset to only methylated fragments for GLM
sample1$methylation_status <- str_extract(sample1$frag_grp, "meth")
sample1$methylation_status <- sample1$methylation_status %>% fct_explicit_na(na_level = 'unmeth')

print(sum(sample1$methylation_status == "meth")/nrow(sample1))
#0.9047528

sample1 <- subset(sample1, sample1$methylation_status == "meth")
sample1$methylation_status <- droplevels(sample1$methylation_status)
sample1$sample1_read_count <- 1

#Aggregate by spike-in control fragment (frag_grp)
sample1 <- sample1 %>%
        group_by(frag_grp) %>%
        summarise(sample1_read_count = sum(sample1_read_count))

###Sample2
sample2 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_aligned.trimmed.6655_S2_L002_R1_001.fastqtrimmed.6655_S2_L002_R2_001.fastq_sorted_dedup_sort.bed")
sample2 <- as.data.frame(sample2[,1])
colnames(sample2) <- c("frag_grp")

#Assess methylation specificity and subset to only methylated fragments for GLM
sample2$methylation_status <- str_extract(sample2$frag_grp, "meth")
sample2$methylation_status <- sample2$methylation_status %>% fct_explicit_na(na_level = 'unmeth')

print(sum(sample2$methylation_status == "meth")/nrow(sample2))
#0.9793557

sample2 <- subset(sample2, sample2$methylation_status == "meth")
sample2$methylation_status <- droplevels(sample2$methylation_status)
sample2$sample2_read_count <- 1

#Aggregate by spike-in control fragment (frag_grp)
sample2 <- sample2 %>%
        group_by(frag_grp) %>%
        summarise(sample2_read_count = sum(sample2_read_count))


###Sample 3
sample3 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_aligned.trimmed.6656_S3_L002_R1_001.fastqtrimmed.6656_S3_L002_R2_001.fastq_sorted_dedup_sort.bed")
sample3 <- as.data.frame(sample3[,1])
colnames(sample3) <- c("frag_grp")

#Assess methylation specificity and subset to only methylated fragments for GLM
sample3$methylation_status <- str_extract(sample3$frag_grp, "meth")
sample3$methylation_status <- sample3$methylation_status %>% fct_explicit_na(na_level = 'unmeth')

print(sum(sample3$methylation_status == "meth")/nrow(sample3))
# 0.9809234

sample3 <- subset(sample3, sample3$methylation_status == "meth")
sample3$methylation_status <- droplevels(sample3$methylation_status)
sample3$sample3_read_count <- 1

#Aggregate by spike-in control fragment (frag_grp)
sample3 <- sample3 %>%
        group_by(frag_grp) %>%
        summarise(sample3_read_count = sum(sample3_read_count))

###Sample4
sample4 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_aligned.trimmed.6657_S4_L002_R1_001.fastqtrimmed.6657_S4_L002_R2_001.fastq_sorted_dedup_sort.bed")
sample4 <- as.data.frame(sample4[,1])
colnames(sample4) <- c("frag_grp")

#Assess methylation specificity and subset to only methylated fragments for GLM
sample4$methylation_status <- str_extract(sample4$frag_grp, "meth")
sample4$methylation_status <- sample4$methylation_status %>% fct_explicit_na(na_level = 'unmeth')

print(sum(sample4$methylation_status == "meth")/nrow(sample4))
#0.9839075

sample4 <- subset(sample4, sample4$methylation_status == "meth")
sample4$methylation_status <- droplevels(sample4$methylation_status)
sample4$sample4_read_count <- 1

#Aggregate by spike-in control fragment (frag_grp)
sample4 <- sample4 %>%
        group_by(frag_grp) %>%
        summarise(sample4_read_count = sum(sample4_read_count))

###Sample5
sample5 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch1_RS/filtered_aligned.trimmed.6658_S5_L002_R1_001.fastqtrimmed.6658_S5_L002_R2_001.fastq_sorted_dedup_sort.bed")
sample5 <- as.data.frame(sample5[,1])
colnames(sample5) <- c("frag_grp")

#Assess methylation specificity and subset to only methylated fragments for GLM
sample5$methylation_status <- str_extract(sample5$frag_grp, "meth")
sample5$methylation_status <- sample5$methylation_status %>% fct_explicit_na(na_level = 'unmeth')

print(sum(sample5$methylation_status == "meth")/nrow(sample5))
#0.9829601

sample5 <- subset(sample5, sample5$methylation_status == "meth")
sample5$methylation_status <- droplevels(sample5$methylation_status)
sample5$sample5_read_count <- 1

#Aggregate by spike-in control fragment (frag_grp)
sample5 <- sample5 %>%
        group_by(frag_grp) %>%
        summarise(sample5_read_count = sum(sample5_read_count))

##merge
data <- merge(sample1, sample2, by= "frag_grp", all = TRUE)
data <- merge(data, sample3, by= "frag_grp", all = TRUE)
data <- merge(data, sample4, by= "frag_grp", all = TRUE)
data <- merge(data, sample5, by= "frag_grp", all = TRUE)


##fill in NAs from UMI data
data$fragment_len <- as.factor(str_sub(data$frag_grp, 0, 3))
data$fragment_len <- revalue(data$fragment_len, c("80b" = "80"))

data$CpG <- as.factor(str_sub(data$frag_grp, 5, 7))
data$CpG <- revalue(data$CpG, c("_16" = "16","_2C" = "2","_4C" = "4","_8C" = "8","1C_" = "1","2C_" = "2","4C_" = "4"))

data$GC <- as.factor(str_sub(data$frag_grp, 7, 11))
data$GC <- revalue(data$GC, c("_35G-" = "35", "_50G-" = "50", "_65G-" = "65", "6C_35" = "35", "6C_50" = "50", "6C_65"="65", "C_35G"="35", "C_50G" = "50", "C_65G" = "65"))

data$grp <- as.factor(paste0(data$fragment_len, sep = '_', data$CpG,sep = '_', data$GC))

##since initial 320bp, 65% G+C is not correct, let's correct that to correct G+C content
data$grp <- revalue(data$grp, c("320_4_65" = "320_04_35"))
data$grp <- revalue(data$grp, c("320_8_65" = "320_08_50"))
data$grp <- revalue(data$grp, c("320_16_65" = "320_016_50"))

##Fix to the correct GC
data$GC <- as.factor(str_sub(data$grp,-2,-1))

data <- as.data.frame(data)

##Add in concentration information
for (i in 1:nrow(data)){
  data[i,11] <- ifelse (data[i,7] == "160", 0.002, ifelse(data[i,7] == "80", 0.004, 0.001))
}

names(data)[[11]] <- "conc"
data$GC <- as.numeric(as.character(data$GC))
data$fragment_len <- as.numeric(as.character(data$fragment_len))
data$CpG <- as.numeric(as.character(data$CpG))

#Make sure the data is in correct format before modelling
stopifnot(str(data$fragment_len) == "num")
stopifnot(str(data$CpG) == "num")
stopifnot(str(data$GC) == "num")
data$CpG_3 <- (data$CpG) ^ (1/3)

set.seed(1234)
results <- list()

#Function to run 100 permutations per sampling of 1-25 spikes for training data
permutations <- function(data,iteration, samp) {
for (j in seq_along(1:iteration)) {

sample_size <- floor(samp/26 * nrow(data))

train_set <- sample(seq_len(nrow(data)), size = sample_size)

train <- data[train_set,]

test <- data[-train_set,]

train_melt <- melt(train, measure.vars = c("sample1_read_count", "sample2_read_count", "sample3_read_count", "sample4_read_count", "sample5_read_count"))

###Gaussian
gaussian_glm <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = train_melt, family = gaussian)
summary(gaussian_glm)

##Calculate R2
#r2_gaussian= 1-(gaussian_glm$deviance/gaussian_glm$null.deviance)

test_melt <-  melt(test, measure.vars = c("sample1_read_count", "sample2_read_count", "sample3_read_count", "sample4_read_count", "sample5_read_count"))

pred_conc <- as.data.frame(predict(gaussian_glm, test_melt))
names(pred_conc) <- c("predicted_pmol")
test_dat <- cbind(test_melt, pred_conc)

#test_dat$residuals <- test_dat$predicted_pmol - test_dat$conc
mean_abs_err <- mae(test_dat$conc, test_dat$predicted_pmol)
print(mean_abs_err)

results[[j]] <- mean_abs_err
}
#maes <- do.call(rbind, results)
#return(maes)
return(results)
}

output <- list()

collect_errors <- function(train_num) { 
for (i in seq_along(1:train_num)) { 
print(i) #it does go in order
errs <- permutations(data, 100, i)
errs <- as.data.frame(unlist(errs))

output[[i]] <- errs
names(output[[i]]) <- paste0("train_",i)
}
all_maes <- do.call(cbind, output)
return(all_maes)
}

data_mae <- collect_errors(25)
names(data_mae) <- c("train_1", "train_2", "train_3", "train_4", "train_5", "train_6", "train_7", "train_8", "train_9", "train_10", "train_11", "train_12", "train_13", "train_14", "train_15", "train_16", "train_17", "train_18", "train_19", "train_20", "train_21", "train_22", "train_23", "train_24", "train_25") 

data_mae <- as.data.frame(data_mae)

#save data for future reproducibility
write.table(data_mae, sep = "\t", file = "2022_data_training_permutations_mae.tsv")

data_mae_t <- t(data_mae)
data_mae_melt <- melt(data_mae_t)

#plot
png("2022_mae_train.png", height = 6, width = 10, units = "in", res = 300)
ggplot(data_mae_melt, aes(x = Var1, y = value)) +
geom_point(size =2) +
labs(x = "Number of spike-in fragments in training data", y = "Mean absolute error") + 
stat_summary(fun = "mean", geom = "point", col = "red") +
scale_x_discrete(labels = c("train_1" = "1", "train_2" = "2", "train_3" = "3",
				"train_4" = "4", "train_5" = "5", "train_6" = "6",
				"train_7" = "7", "train_8" =  "8", "train_9" = "9",
				"train_10" = "10", "train_11" = "11", "train_12" = "12",
				"train_13" = "13", "train_14" = "14", "train_15" = "15",
				"train_16" = "16", "train_17" = "17", "train_18" = "18",
				"train_19" = "19", "train_20" = "20", "train_21" =  "21",
				"train_22" = "22", "train_23" =  "23", "train_24" = "24", "train_25" = "25"))
dev.off()

#boxplot
png("2022_mae_train_boxplot.png", height = 6, width = 10, units = "in", res = 300)
ggplot(data_mae_melt, aes(x = Var1, y = value)) +
geom_boxplot() +
labs(x = "Number of spike-in fragments in training data", y = "Mean absolute error") +
scale_x_discrete(labels = c("train_1" = "1", "train_2" = "2", "train_3" = "3",
                                "train_4" = "4", "train_5" = "5", "train_6" = "6",
                                "train_7" = "7", "train_8" =  "8", "train_9" = "9",
                                "train_10" = "10", "train_11" = "11", "train_12" = "12",
                                "train_13" = "13", "train_14" = "14", "train_15" = "15",
                                "train_16" = "16", "train_17" = "17", "train_18" = "18",
                                "train_19" = "19", "train_20" = "20", "train_21" =  "21",
                                "train_22" = "22", "train_23" =  "23", "train_24" = "24", "train_25" = "25"))
dev.off()

#boxplot starting at 4 spikes in training
data_mae_6 <- data_mae[,-c(1:5)]

data_mae_t <- t(data_mae_6)
data_mae_melt <- melt(data_mae_t)

png("2022_mae_train_6boxplot.png", height = 6, width = 10, units = "in", res = 300)
ggplot(data_mae_melt, aes(x = Var1, y = value)) +
geom_boxplot() +
labs(x = "Number of spike-in fragments in training data", y = "Mean absolute error (pmol)") +
scale_x_discrete(labels = c("train_6" = "6",
                                "train_7" = "7", "train_8" =  "8", "train_9" = "9",
                                "train_10" = "10", "train_11" = "11", "train_12" = "12",
                                "train_13" = "13", "train_14" = "14", "train_15" = "15",
                                "train_16" = "16", "train_17" = "17", "train_18" = "18",
                                "train_19" = "19", "train_20" = "20", "train_21" =  "21",
                                "train_22" = "22", "train_23" =  "23", "train_24" = "24", "train_25" = "25")) +
scale_y_continuous(expand = c(0, 0), limits = c(0,0.002), breaks = c(0, 0.001, 0.002))
dev.off()

#log the last box plot
png("2022_mae_train_6logboxplot.png", height = 6, width = 10, units = "in", res = 300)
ggplot(data_mae_melt, aes(x = Var1, y = log10(value))) +
geom_boxplot() +
labs(x = "Number of spike-in fragments in training data", y = "Log10(Mean absolute error)") +
scale_x_discrete(labels = c("train_6" = "6",
                                "train_7" = "7", "train_8" =  "8", "train_9" = "9",
                                "train_10" = "10", "train_11" = "11", "train_12" = "12",
                                "train_13" = "13", "train_14" = "14", "train_15" = "15",
                                "train_16" = "16", "train_17" = "17", "train_18" = "18",
                                "train_19" = "19", "train_20" = "20", "train_21" =  "21",
                                "train_22" = "22", "train_23" =  "23", "train_24" = "24", "train_25" = "25"))
dev.off()

##Sinaplot with top 99 data points out of 100
png("2022_mae_train_6sinaplot.png", height = 6, width = 10, units = "in", res = 300)
ggplot(data_mae_melt, aes(x = Var1, y = value)) +
geom_sina() +
labs(x = "Number of spike-in fragments in training data", y = "Mean absolute error (pmol)") +
scale_x_discrete(labels = c("train_6" = "6",
                                "train_7" = "7", "train_8" =  "8", "train_9" = "9",
                                "train_10" = "10", "train_11" = "11", "train_12" = "12",
                                "train_13" = "13", "train_14" = "14", "train_15" = "15",
                                "train_16" = "16", "train_17" = "17", "train_18" = "18",
                                "train_19" = "19", "train_20" = "20", "train_21" =  "21",
                                "train_22" = "22", "train_23" =  "23", "train_24" = "24", "train_25" = "25")) +
scale_y_continuous(expand = c(0, 0), limit = c(0,0.002), breaks = c(0, 0.001, 0.002))
dev.off()



##Let's look at the model with 1 spike in training looks like - why is it better than we expect?
S6548_corr + annotate(geom = "text", x = 0.002, y = 0.004, label = "r = 0.96, p-value = 2.89e-05")
##Correlation plot between predicted and expected
##plot function
corr_plot <- function(data, x, y){
ggplot(data, aes(x= x, y= y, color = as.factor(fragment_len))) +
  geom_point(size = 2) +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color='#2C3E50') +
  ylab("Predicted picomole") +
  xlab("Actual picomole") +
 # stat_cor(method = "pearson") +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  scale_colour_manual(name = "Fragment length",
                        values = c("black", "darkgrey", "steelblue"),
                        breaks = c("80", "160", "320"),
                        labels = c("80 bp", "160 bp", "320 bp"))
}         


gaussian_glm_function <- function(train, test) {
###Gaussian
gaussian_glm <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = train, family = gaussian)
summary(gaussian_glm)

##Calculate R2
r2_gaussian= 1-(gaussian_glm$deviance/gaussian_glm$null.deviance)
print(r2_gaussian)

##predict conc in testing
stopifnot(str(test_melt$fragment_len) == "num")
stopifnot(str(test_melt$CpG) == "num")
stopifnot(str(test_melt$CpG_3) == "num")
stopifnot(str(test_melt$GC) == "num")

pred_conc <- as.data.frame(predict(gaussian_glm, test, se = TRUE))
names(pred_conc) <- c("predicted_pmol", "std_err", "resid_scale")
test_dat <- cbind(test, pred_conc)

test_dat$residuals <- test_dat$conc - test_dat$predicted_pmol
print(head(test_dat))

return(test_dat)

}

sample_size <- floor(1/26 * nrow(data))

train_set <- sample(seq_len(nrow(data)), size = sample_size)

train <- data[train_set,]

test <- data[-train_set,]

train_melt <- melt(train, measure.vars = c("sample1_read_count", "sample2_read_count", "sample3_read_count", "sample4_read_count", "sample5_read_count"))
test_melt <-  melt(test, measure.vars = c("sample1_read_count", "sample2_read_count", "sample3_read_count", "sample4_read_count", "sample5_read_count"))

spike_1_training <- gaussian_glm_function(train_melt, test_melt)
         
png("2022_actual_predict_spike1_corr.png", height = 6, width = 6, unit = "in", res = 300)
corr_plot(spike_1_training, spike_1_training$conc, spike_1_training$predicted_pmol)
dev.off()



