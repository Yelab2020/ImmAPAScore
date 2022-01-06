library(tidyverse)
library(magrittr)
library(pheatmap)
library(Rtsne)
library(ggfortify)
library(mvtnorm)
library(ggfortify)
library(randomcoloR)
library(RColorBrewer)
library(ggpubr)

mypalette<-brewer.pal(12,"Paired")

Imm_APA_Score <- readRDS("/work/gywang/project/Immune_APA/Data/TCGA_Data/Imm_APA_Score.rds.gz")
TCGA_12_APA_DIFF <- readr::read_rds("/home/gywang/project/APA_CancerImmune/Data/12Normal_Tumor_APA/TCGA_12types_APA_T_N_Diff.rds.gz")
APAevents_count <- Imm_APA_Score$APAevents_count[[1]]
ImmAPA_events <- APAevents_count[APAevents_count$CanTypes10Path15 %in% "Sig",]$APA.events

ImmPathways <- readr::read_rds("/home/yye/Project/2020APA_Imm/ImmunePathwayList.rds")
ImmPathwaysD <- data.frame(term=rep(names(ImmPathways),times=unlist(lapply(ImmPathways,length))),gene=unlist(ImmPathways))

#  qvalue < 0.1
All_ImmAPA_RankScore <- Imm_APA_Score$RankScore %>% lapply(function(x){
  dplyr::filter(x,APAevents %in% ImmAPA_events)
}) %>% bind_rows()
All_ImmAPA_ESS <- All_ImmAPA_RankScore %>% split(All_ImmAPA_RankScore$Description) %>% lapply(function(x){
  x$FDR_sig <- ifelse(x$qvalues < 0.1,1,0)
  x$ESS <- ifelse(x$NES <0 ,- x$FDR_sig,x$FDR_sig)
  sub_ImmPath <- split(x,x$APAevents) %>% lapply(function(x){
    sum(x$ESS)
  }) %>% bind_rows()
  rownames(sub_ImmPath) <- unique(x$Description)
  return(sub_ImmPath)
}) %>% bind_rows()

All_ImmAPA_ESS <- All_ImmAPA_ESS %>% apply( 2, function(x){
  x <-  ifelse(is.na(x),0,x)
return(x)
  })


All_ImmAPA_ESS_n <- All_ImmAPA_ESS %>% t() %>% as.data.frame() 
All_ImmAPA_ESS_n$Sig_Path_count <- All_ImmAPA_ESS_n %>% apply(1,function(x){length(which(abs(x)>5))})
All_ImmAPA_ESS_n <- All_ImmAPA_ESS_n %>% filter(Sig_Path_count >5)

pheatmap::pheatmap(t(All_ImmAPA_ESS),scale = "none",cluster_col = T,
                   cluster_rows = T,display_numbers = F,number_format = "%.1e",
                   #cellwidth=,#  cellheight = ,
                   show_rownames= F,
                   show_colnames = T,
                   fontsize=17,
                   angle_col = 270,
                   main = "Ranked by regulation to Immune Pathway",
                   filename = "/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/Heatmap of 1991 ImmAPA  by ESS.pdf",
                   width = 11,height = 15)
pheatmap::pheatmap(t(All_ImmAPA_ESS_n[,-17]),scale = "none",cluster_col = T,
                   cluster_rows = T,display_numbers = F,number_format = "%.1e",
                   #cellwidth=,#  cellheight = ,
                   show_rownames= F,
                   show_colnames = T,
                   fontsize=17,
                   angle_col = 270,
                   main = "Ranked by regulation to Immune Pathway",
                   filename = "/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/Heatmap of 541 ImmAPA  by ESS.pdf",
                   width = 11,height = 15)


readr::write_rds(All_ImmAPA_ESS,"/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_1991_ES_score.rds")
readr::write_rds(All_ImmAPA_ESS_n,"/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_score.rds")



All_ImmAPA_ESS <- All_ImmAPA_ESS_n[,-17]
All_ImmAPA_ESS_PCA <- prcomp(All_ImmAPA_ESS)
pdf("/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/PCA of 541 ImmAPA clustering by ESS.pdf",width = 10,height = 7)
autoplot(All_ImmAPA_ESS_PCA,data=All_ImmAPA_ESS)
dev.off()

set.seed(123)
All_ImmAPA_ESS_tSNE <- Rtsne(
  All_ImmAPA_ESS,
  dims = 2,
  pca = T,
  perplexity = 10,
  #theta = 0.0,
  max_iter = 1000,
  check_duplicates = FALSE
)
plot(All_ImmAPA_ESS_tSNE$Y)
tsnes_plot=All_ImmAPA_ESS_tSNE$Y %>% as.data.frame()
tsnes_plot$APA_events <- rownames(All_ImmAPA_ESS)
meta_info_ImmAPA <- ifelse(rownames(All_ImmAPA_ESS) %in% rownames(top20_from_1991),"top20ImmAPA","noneTOP")

colnames(tsnes_plot)[1:2] <- c("tSNE1","tSNE2")
tsnes_plot$group <- meta_info_ImmAPA
ggplot(tsnes_plot, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=group))


plot_list <- list()
for (i in c(2,4,6,8,10,12)) {
  tsnes_kmean <- kmeans(tsnes_plot[,1:2],i,iter.max = 10000)
  tsnes_plot$kean_cluster <- paste("Cluster_",tsnes_kmean$cluster ,sep = "")
  p1 <- ggplot(tsnes_plot, aes(x = tSNE1, y = tSNE2)) + 
    geom_point(aes(color=kean_cluster))+
    scale_color_manual(values = mypalette)
  plot_list[[i/2]] <- p1 
  All_ImmAPA_ESS_cluster <- cbind(All_ImmAPA_ESS,tsnes_plot)
  ESS_mean_by_cluster <- All_ImmAPA_ESS_cluster %>% split(All_ImmAPA_ESS_cluster$kean_cluster) %>% lapply(function(x){
    sub <- x[,1:16] %>% apply(2,mean)
    sub <- c(cluster=unique(x$kean_cluster),sub)
    return(sub)
  }) %>% bind_rows()
  ESS_mean_by_cluster <- ESS_mean_by_cluster  %>% dplyr::select(-1) %>% apply(2,as.numeric) %>% set_rownames(ESS_mean_by_cluster$cluster)
  p2 <- pheatmap::pheatmap(t(ESS_mean_by_cluster))
  pdf(paste("/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/Heatmap of 541 ImmAPA clustering to ",i,"clusters by ESS.pdf"),width = 7,height = 8)
  print(p2)
  dev.off()
}

pdf("/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/tSNE of ImmAPA clustering by ESS.pdf",width = 21,height = 15)
cowplot::plot_grid(plotlist = plot_list,ncol = 3)
dev.off()

All_ImmAPA_ESS_cluster <- cbind(All_ImmAPA_ESS,tsnes_plot)
ESS_mean_by_cluster <- All_ImmAPA_ESS_cluster %>% split(All_ImmAPA_ESS_cluster$kean_cluster) %>% lapply(function(x){
  sub <- x[,1:16] %>% apply(2,mean)
  sub <- c(cluster=unique(x$kean_cluster),sub)
  return(sub)
}) %>% bind_rows()


ESS_mean_by_cluster <- ESS_mean_by_cluster  %>% dplyr::select(-1) %>% apply(2,as.numeric) %>% set_rownames(ESS_mean_by_cluster$cluster)
pheatmap::pheatmap(ESS_mean_by_cluster)

colnames(All_ImmAPA_ESS_cluster)[19] <- "APA_events"
readr::write_rds(All_ImmAPA_ESS_cluster,"/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_all_cluster.rds.gz",compress = "gz")












