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

##Load data

S6547 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_aligned.trimmed.6547_S3_L_R1_001.fastqtrimmed.6547_S3_L_R2_001.fastq_sorted_dedup_sort.bed")
S6547 <- as.data.frame(S6547[,1])
colnames(S6547) <- c("frag_grp")

#Assess methylation specificity and subset to only methylated fragments for GLM
S6547$methylation_status <- str_extract(S6547$frag_grp, "meth")
S6547$methylation_status <- S6547$methylation_status %>% fct_explicit_na(na_level = 'unmeth')

print(sum(S6547$methylation_status == "meth")/nrow(S6547))
# 0.9724638

S6547 <- subset(S6547, S6547$methylation_status == "meth")
S6547$methylation_status <- droplevels(S6547$methylation_status)
S6547$S6547_read_count <- 1

#Aggregate by spike-in control fragment (frag_grp)
S6547 <- S6547 %>%
        group_by(frag_grp) %>%
        summarise(S6547_read_count = sum(S6547_read_count))

S6548 <- read.table("/cluster/projects/hoffmangroup/data_samanthawilson/2020_0.01ng_HCT116/filtered_aligned.trimmed.6548_S4_L_R1_001.fastqtrimmed.6548_S4_L_R2_001.fastq_sorted_dedup_sort.bed")
S6548 <- as.data.frame(S6548[,1])
colnames(S6548) <- c("frag_grp")

#Assess methylation specificity and subset to only methylated fragments for GLM
S6548$methylation_status <- str_extract(S6548$frag_grp, "meth")
S6548$methylation_status <- S6548$methylation_status %>% fct_explicit_na(na_level = 'unmeth')

print(sum(S6548$methylation_status == "meth")/nrow(S6548))
#0.9783167

S6548 <- subset(S6548, S6548$methylation_status == "meth")
S6548$methylation_status <- droplevels(S6548$methylation_status)
S6548$S6548_read_count <- 1

#Aggregate by spike-in control fragment (frag_grp)
S6548 <- S6548 %>%
        group_by(frag_grp) %>%
        summarise(S6548_read_count = sum(S6548_read_count))


##merge
data <- merge(S6547, S6548, by= "frag_grp", all = TRUE)

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
  data[i,8] <- ifelse (data[i,4] == "160", 0.002, ifelse(data[i,4] == "80", 0.004, 0.001))
}

names(data)[8] <- "conc"

##Adjusting for the 0.01ng dilution
data$conc <- data$conc*0.9

##Some exploratory plots
###Melt so that data frame 
data_melt <- melt(data, id.vars = c("conc", "grp", "GC", "CpG", "fragment_len"), measure.vars = c("S6547_read_count", "S6548_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

png("2020_concvsreads_dedup_0.01ng.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data_melt, aes(x = conc, y = value)) + geom_point() +
  theme_bw()+
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Concentration (pmol/ng)")+
  ylab("Read counts")+
  scale_y_continuous(labels = comma)+
  scale_fill_manual("slategray3")
dev.off()

##Colour points by GC content
data$GC <- as.factor(data$GC)

png("2020_concvsreads_dedup_GCcol_0.01ng.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data_melt, aes(x = conc, y = value, col = GC)) + geom_point() +
  theme_bw()+
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12))+
  xlab("Concentration (pmol/ng)")+
  ylab("Read counts")+
  scale_y_continuous(labels = comma)
dev.off()

data$GC <- as.numeric(as.character(data$GC))

##Correlation between read count and concentration
cor(data_melt$conc, data_melt$value)#-0.0956277

##Cube root CpG to try to return distribution to normality
data$CpG <- as.integer(as.character(data$CpG))
data$CpG_3 <- data$CpG ^(1/3)

##cube root
data_sub <- data[, c("frag_grp", "S6547_read_count", "S6548_read_count", "fragment_len", "CpG_3", "GC", "conc")]
data_melt <- melt(data_sub, id.vars = c("conc", "fragment_len", "CpG_3", "GC", "frag_grp"), measure.vars = c("S6547_read_count", "S6548_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

png("2020_CpG_3_dedup_0.01ng.png", height = 6, width = 6, units = "in", res = 300)
ggplot(data = data_melt, aes(x = CpG_3, y = value, width = 0.1)) +
  geom_bar(stat = "identity", fill = "peachpuff") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), axis.text.x  = element_text(angle = 45, vjust = 0.5, size = 12), axis.text.y = element_text(vjust = 0.5, size = 12)) +
  xlab("Number of CpGs") +
  ylab("Read counts") +
  scale_y_continuous(labels = comma, limits = c(0, 80000), breaks = c(20000, 40000, 60000, 80000))
dev.off()
#better, this will be used

###Let's try out some models now
data$fragment_len <- as.numeric(as.character(data$fragment_len))
data$GC <- as.numeric(as.character(data$GC))

data_sub <- data[, c("frag_grp", "S6547_read_count", "S6548_read_count", "fragment_len", "CpG_3", "GC", "conc")]
data_melt <- melt(data_sub, id.vars = c("conc", "fragment_len", "CpG_3", "GC", "frag_grp"), measure.vars = c("S6547_read_count", "S6548_read_count"))
data_melt$value <- as.numeric(as.character(data_melt$value))

##Set ggplot theme
theme_set(theme_bw() +
  theme(legend.title = element_blank(), legend.position = "right", legend.text = element_text(size = 14),
        axis.title.y = element_text(size = 20), axis.text.x  = element_text(size = 16),
        axis.text.y = element_text(vjust = 0.5, size = 16),axis.title.x = element_text(size = 20))
)


#Make sure the data is in correct format before modelling
stopifnot(str(data_melt$fragment_len) == "num")
stopifnot(str(data$CpG) == "num")
stopifnot(str(data_melt$CpG_3) == "num")
stopifnot(str(data_melt$GC) == "num")

##Poisson
###Poisson- Getting negatives
poisson_glm <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = data_melt, family = poisson(link = "log"))
summary(poisson_glm)#In dpois(y, mu, log = TRUE) : non-integer x = 0.000842

##Calculate R2
r2_poisson <- 1 - (poisson_glm$deviance / poisson_glm$null.deviance)
r2_poisson # 0.9829384

png("2020_GLM_fraglen_CpG3_0.01ngdiag.png")
glm.diag.plots(poisson_glm, glmdiag = glm.diag(poisson_glm))
dev.off()

png("2020_GLM_fraglen_CpG3_poisson_0.01ng.png")
plot(poisson_glm)
dev.off()

png("2020_BA_poisson_fraglen_CpG3_0.01ng.png")
bland.altman.plot(poisson_glm$data[,1], poisson_glm$resid, conf.int=.95, pch=19) ##Why are there negative values here
dev.off()

poisson_glm$data[,8] <- as.factor(as.character(poisson_glm$data[, 2]))
names(poisson_glm$data[,8]) <- "Fragment_len"

png("2020_BA_poisson_fraglen_CpG3_0.01ng.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(poisson_glm$data[,1], poisson_glm$resid, conf.int = .95, col.points = poisson_glm$data[,2], pch = 19, graph.sys = "ggplot2") 
  
BA + 
  xlab("Mean of measurements (femotomoles)") +
  ylab("Difference (picomoles)") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y = c(-0.05, -0.06, -0.07)) + 
  aes(color = poisson_glm$data[,8]) + 
  scale_color_manual(values = c("black", "darkgrey", "navy"))

dev.off()

###Gaussian
gaussian_glm <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = data_melt, family = gaussian)
summary(gaussian_glm)

png("2020_GLM_gaussian_fraglen_CpG3_0.01ngdiag.png")
glm.diag.plots(gaussian_glm, glmdiag = glm.diag(gaussian_glm))
dev.off()

png("2020_GLM_fraglen_CpG3_gaussian0.01ng.png")
plot(gaussian_glm)
dev.off()

png("2020_BA_gaussian_fraglen_CpG30.01ng.png")
bland.altman.plot(gaussian_glm$data[,1], gaussian_glm$resid, conf.int=.95, pch=19)
dev.off()

gaussian_glm$data[,8] <- as.factor(as.character(gaussian_glm$data[, 2]))
names(gaussian_glm$data[,8]) <- "Fragment_len"

png("2020_BA_gaussian_fraglen_CpG3_0.01ng.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(gaussian_glm$data[,1], gaussian_glm$resid, conf.int = .95, col.points = gaussian_glm$data[,2], pch = 19, graph.sys = "ggplot2")

BA +
  xlab("Mean of measurements (femotomoles)") +
  ylab("Difference (picomoles)") +
  stat_cor(method = "spearman", label.x = 0.001, label.y = c(-0.0015, -0.0018, -0.0021)) +
  aes(color = gaussian_glm$data[,8]) + 
  scale_color_manual(values = c("black", "darkgrey", "navy"))

dev.off()

##This estimates the 160bp very nicely

##Calculate R2
r2_gaussian= 1-(gaussian_glm$deviance/gaussian_glm$null.deviance)
r2_gaussian #0.9429097

##Correlation isn't great, but makes sense in the context of wanting 160bp to perform better, best AIC I've seen
save(gaussian_glm, file = "2020_Gaussian_0.01ng.rda")

###Gaussian log
gaussian_log <- glm(formula = conc ~ value + fragment_len + GC + CpG_3, data = data_melt, family = gaussian(link = "log"))
summary(gaussian_log)

png("2020_GLM_gaussianlog_fraglen_CpG3_diag_0.01ng.png")
glm.diag.plots(gaussian_log, glmdiag = glm.diag(gaussian_log))
dev.off()

png("2020_GLM_gaussianlog_fraglen_CpG3_0.01ng.png")
plot(gaussian_log)
dev.off()

##Bland-Altman plot (plots the residuals against the truth and looks at the deviance)
#gaussian_log$data[,2]=as.factor(as.character(gaussian_log$data[,2]))
gaussian_log$data[,8] <- as.factor(as.character(gaussian_log$data[,2]))
names(gaussian_log$data[,8]) <- "Fragment_len"

png("2020_BA_guassianlog_fraglen_CpG3_0.01ng.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(exp(gaussian_log$data[,1]), gaussian_log$resid, conf.int = .95, col.points = gaussian_log$data[,2], pch = 19, graph.sys="ggplot2")

BA +
  xlab("Mean of measurements (femotomoles)") +
  ylab("Difference (picomoles)") +
  stat_cor(method = "spearman", label.x.npc = "left", label.y = c(0.9, 0.88, 0.86)) +
  aes(color = gaussian_log$data[,8]) + 
  scale_color_manual(values = c("black", "darkgrey", "navy"))

dev.off()

##Calculate R2
r2_gaussian_log = 1-(gaussian_log$deviance/gaussian_log$null.deviance)
r2_gaussian_log #0.9848901

##Seeing if a simple lm would perform better
lm <- lm(conc ~ value + fragment_len + GC + CpG_3, data = data_melt)
summary(lm)#R2=0.9343, adjusted=0.9287 
##residual standard error -0.0005083

png("2020_lm_fraglen_CpG3_0.01ng.png")
plot(lm)
dev.off()

png("2020_lm_fraglen_CpG3_diag_0.01ng.png")
glm.diag.plots(lm, glmdiag = glm.diag(lm))
dev.off()

png("2020_BA_lm_fraglen_CpG3_0.01ng.png")
bland.altman.plot(lm$model[,1], lm$resid, conf.int=.95, pch=19)
dev.off()

lm$model[,6] <- as.factor(as.character(lm$model[,3]))
names(lm$model[,6]) <- "Fragment_len"

png("2020_BA_lm_fraglen_CpG3_0.01ng.png", height = 6, width = 6, units = "in", res = 300)
BA <- bland.altman.plot(lm$model[,1], lm$residuals, conf.int = .95, pch = 19, graph.sys="ggplot2")

BA +
  xlab("Mean of measurements (picomoles)") +
  ylab("Difference (picomoles)") +
  stat_cor(method = "spearman", label.x = 0.0012, label.y = c(0.0016, 0.0014, 0.0012)) +
  aes(color = as.factor(lm$model[, 3])) +    
  scale_color_manual(values = c("black", "darkgrey", "navy"))

dev.off()

save(lm, file = "2020_lm_0.01ng.rda")
