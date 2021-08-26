###Heat scree effect size
heat_scree_plot_es<-function(Loadings, Importance, Num, Order){
  #adjust<-1-Importance[1]
  pca_adjusted<-Importance[1:length(Importance)]
  pca_df<-data.frame(adjusted_variance=pca_adjusted, PC=seq(1:length(pca_adjusted)))

  scree<-ggplot(pca_df[which(pca_df$PC<Num),],aes(PC,adjusted_variance)) +
        #geom_bar(stat = "identity",color="black",fill="grey") + 
        theme_bw() +
        geom_segment(aes(x = PC, xend = PC, y = 0, yend = adjusted_variance), color = "black", size = 1) +
        geom_point(size = 3, color = "black", fill=alpha("black", 1.0), alpha=1.0, shape=21, stroke=2) +
        geom_text(label = round(pca_df[which(pca_df$PC<Num),]$adjusted_variance,2), vjust = -1) +
        theme(axis.text = element_text(size =12),
              axis.title = element_text(size =15),
              plot.margin=unit(c(1,1.5,0.2,2.25),"cm"))+ylab("Variance")+
              scale_y_continuous(limits = c(0,1), breaks = c(0, 0.2,0.4, 0.6, 0.8,1.0)) +
    scale_x_discrete(limits = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8")) +
        xlab("Principal component")

  #### Heat
  ## correlate meta with PCS
  ## Run anova of each PC on each meta data variable
  aov_PC_meta<-lapply(1:ncol(meta_categorical), function(covar) sapply(1:ncol(Loadings), function(PC) summary(aov(Loadings[,PC]~meta_categorical[,covar]))[[1]]$"F value"[1]))
  #cor_PC_meta<-lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Loadings), function(PC) (cor.test(Loadings[,PC],as.numeric(meta_continuous[,covar]),alternative = "two.sided", method="spearman", na.action=na.omit)$p.value)))
  names(aov_PC_meta)<-colnames(meta_categorical)
  #names(cor_PC_meta)<-colnames(meta_continuous)
  aov_PC_meta<-do.call(rbind, aov_PC_meta)
  #cor_PC_meta<-do.call(rbind, cor_PC_meta)
  #aov_PC_meta<-rbind(aov_PC_meta, cor_PC_meta)
  aov_PC_meta<-as.data.frame(aov_PC_meta)
  #adjust
  aov_PC_meta_adjust<-aov_PC_meta[,2:ncol(aov_PC_meta)]

  ##put in the sample size for each variable
  aov_PC_meta_adjust$size <- 0
  for (i in rownames(aov_PC_meta_adjust)) {
  aov_PC_meta_adjust[i, "size"] <- length(levels(meta_categorical[[i]]))
  }

  ##converting F-statistic to effect size
 meandis <- array(data = NA, dim = dim(aov_PC_meta_adjust))
 rownames(meandis) <- rownames(aov_PC_meta_adjust)
 meanp <- meandis
 for(i in 1:(ncol(meandis)-1)){
 for(j in 1:nrow(meandis)){
 meandis[j,i] <- as.numeric(fes(aov_PC_meta_adjust[j,i], aov_PC_meta_adjust[j, "size"], length(Loadings), verbose = FALSE)["d"])
 }
 }

 for(i in 1:(ncol(meandis)-1)){
 for(j in 1:nrow(meandis)){
 meanp[j,i] <- as.numeric(fes(aov_PC_meta_adjust[j,i], aov_PC_meta_adjust[j, "size"], length(Loadings), verbose = FALSE)["pval.d"])
 }
 }
#reshape - mean distance
  avo <- meandis[,1:(Num-1)]
  avo_heat_num<-apply(avo,2, as.numeric)
  avo_heat<-as.data.frame(avo_heat_num)
  colnames(avo_heat)<-sapply(1:(Num-1), function(x) paste("PC",x, sep=""))
  avo_heat$meta<-rownames(avo)
  avo_heat_melt<-melt(avo_heat, id=c("meta"))

  # cluster meta data
  ord <- Order
  meta_var_order<-unique(avo_heat_melt$meta)[rev(ord)]
  avo_heat_melt$meta <- factor(avo_heat_melt$meta, levels = meta_var_order)

  ##reshape - distance p value
  avo_p <- meanp[,1:(Num-1)]
  avo_heat_num_p<-apply(avo_p,2, as.numeric)
  avo_heat_p<-as.data.frame(avo_heat_num_p)
  colnames(avo_heat_p)<-sapply(1:(Num-1), function(x) paste("PC",x, sep=""))
  avo_heat_p$meta<-rownames(avo_p)
  avo_heat_melt_p<-melt(avo_heat_p, id=c("meta"))

  # cluster meta data
  meta_var_order_p<-unique(avo_heat_melt_p$meta)[rev(ord)]
  avo_heat_melt_p$meta <- factor(avo_heat_melt_p$meta, levels = meta_var_order_p)
 avo_heat_melt_p$p_value <- avo_heat_melt_p$value
 avo_heat_melt_p$value <- NULL

avo_heat_melt <- merge(avo_heat_melt, avo_heat_melt_p, by = c("meta", "variable"))

##adjust p for multipl test correction
avo_heat_melt$adj_p <- p.adjust(avo_heat_melt$p_value, method = "holm", n = nrow(avo_heat_melt))

##labelling significance for heatmap
avo_heat_melt$stars <- cut(avo_heat_melt$adj_p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
 heat<-ggplot(avo_heat_melt, aes(variable,meta, fill = value)) +
  geom_tile(color = "black",size=0.5) +
  theme_bw() +
        scale_fill_gradient(limits = c(0,5), high = "#2c7fb8", low = "#edf8b1") +
      theme(axis.text = element_text(size =10, color="black"),
            axis.text.x = element_text(),
          axis.title = element_text(size =15),
          legend.text = element_text(size =14),
          legend.title = element_text(size =12),
	  legend.position = "none",	
          legend.position = c(1, 0), legend.justification = c(1,0),
          plot.margin=unit(c(0,2.25,1,1),"cm")) +
        labs(fill = "Cohen's d") +
        xlab("Principal component") +
        ylab(NULL) +
        geom_text(aes(label = stars), color = "black", size = 5)

  plot_grid(scree, heat, ncol = 2, align = "h")

}

