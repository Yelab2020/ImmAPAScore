library(tidyverse)
library(magrittr)
library(pheatmap)
library(Rtsne)
library(scatterpie)
library(stats)

options(stringsAsFactors = F)

my.cor.test_spearman <- function(...) {
  obj<- try(cor.test(...,method="spearman"), silent=TRUE)
  if (is(obj, "try-error"))  return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}

my.cor.test_pearson <- function(...) {
  obj<- try(cor.test(...,method="pearson"), silent=TRUE)
  if (is(obj, "try-error"))return(c("NA","NA")) else return(obj[c("estimate","p.value")])
}

TCGA_APA <- readr::read_rds("/work/gywang/project/Immune_APA/Data/TCGA_Data/pancan32_APA.DATA.rds.gz")
tcga_infiltration <- read.csv("/work/gywang/project/Immune_APA/Data/TCGA_Data/TCGA_ICI/infiltration_estimation_for_tcga.csv")
tcga_infiltration$cell_type %>% str_sub(14,15) %>% table()
TCGA_TIMER2 <- tcga_infiltration[,1:7] %>% set_colnames(c("barcode", "B_cell","T_cell.CD4." ,"T_cell.CD8." ,"Neutrophil" ,"Macrophage"  ,           "Myeloid"))

All_ImmAPA_ESS_cluster <- readr::read_rds("/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_all_cluster.rds.gz")
CYT_PRF1_GZMA <- readr::read_rds("/work/yye/AliRstudio/Project/2018Immunotherapy/ProcessedData/CYT_PRF1_GZMA.rds")
TIFScore <- readr::read_rds("/work/yye/AliRstudio/Project/2018Immunotherapy/ProcessedData/TIFScore.rds.gz")
GEPScore <- readr::read_rds("/work/yye/AliRstudio/Project/2018Immunotherapy/ProcessedData/GEPScore.rds.gz")
All_score_APA <- merge(TCGA_APA,CYT_PRF1_GZMA,by = "cancer_types")
All_score_APA <- merge(All_score_APA,TIFScore,by = "cancer_types")
All_score_APA <- merge(All_score_APA,GEPScore,by = "cancer_types")

All_score_APA <- All_score_APA %>% mutate(All_immune_feature=map2(All_score_APA$CYT,All_score_APA$TIFscore, function(x,y){
  merge(x,y,by = "barcode")
}))
All_score_APA <- All_score_APA %>% mutate(All_immune_feature=map2(All_score_APA$All_immune_feature,All_score_APA$GEPscore, function(x,y){
  sub <- merge(x,y,by = "barcode")
}))

All_score_APA <- All_score_APA %>% mutate(All_immune_feature=lapply(All_score_APA$All_immune_feature, function(x){
  x$barcode <- gsub("\\.","-",str_sub(x$barcode,1,15))
  sub <- merge(x,TCGA_TIMER2,by= "barcode")
  return(sub)
}))


All_score_APA <- All_score_APA %>% mutate(ImmAPA_Immfeature=map2(All_score_APA$All_immune_feature,All_score_APA$APA, function(x,y){
  y <- y[y$event_id %in% rownames(All_ImmAPA_ESS_cluster),-1:-3] %>% t() %>% as.data.frame() %>% set_colnames(y[y$event_id %in% rownames(All_ImmAPA_ESS_cluster),]$event_id)
  y$barcode <- gsub("\\.","-",str_sub(rownames(y),1,15))
  sub <- merge(x,y,by = "barcode")
  return(sub)
}))

readr::write_rds(All_score_APA,"/work/gywang/project/Immune_APA/Data/ImmAPA/541_ImmAPA_PDUI_ImmFeature.rds.gz")


lapply(unique(All_ImmAPA_ESS_cluster$combined_3cluster), function(cluster){
  purrr::map2(All_score_APA$ImmAPA_Immfeature, All_score_APA$cancer_types,function(x,cans){
    #x <- All_score_APA$ImmAPA_Immfeature[[6]]
    APAevens <- intersect(colnames(x) ,rownames(All_ImmAPA_ESS_cluster[All_ImmAPA_ESS_cluster$combined_3cluster %in% cluster,]))
    x <- dplyr::select(x,1:10,APAevens)
    lapply(2:10, function(pos){
      features <- x[,pos]
      sub <-  apply(x[,-1:-10],2,function(apa){
        corr <- my.cor.test_pearson(features,as.numeric(apa))
        ifelse(as.character( corr[1]) %in%  "NA",return(c(Corr =NA ,P=NA)), return(c(Corr =corr$estimate ,P=corr$p.value )) )
      }) %>% t() %>% as.data.frame()
      sub$fdr <- p.adjust(sub$P,method="fdr")
      sub$sig <- ifelse(sub$fdr > 0.1,"NonSig",ifelse(sub$Corr.cor > 0,"PositiveCorr","NegativeCorr" ))
      sub1 <- c(ImmFeature= colnames(x)[pos],NegativeCorr=length(sub[sub$sig %in% "NegativeCorr", ]$sig),NonSig=length(sub[sub$sig %in% "NonSig", ]$sig),PositiveCorr=length(sub[sub$sig %in% "PositiveCorr", ]$sig))
      return(sub1)
    }) %>% bind_rows() -> sigCanCluster
    sigCanCluster$CancerType <- cans
    return(sigCanCluster)
  }) %>% bind_rows() -> oncluster
  oncluster$ImmAPAcluster <- cluster
  tt <- oncluster
  tt$ImmFeature <- factor(tt$ImmFeature,levels= unique(tt$ImmFeature))
  tt$x_posi <- as.numeric(factor(tt$ImmFeature))
  tt$y_posi <- as.numeric(factor(tt$CancerType))
  tt$group <- 1:279
  tt$r <- 0.4
  tt[,2:4] <- apply(tt[2:4],2,as.numeric)
  
  p <- ggplot(tt)+
    geom_scatterpie(aes(x=x_posi ,y=y_posi,group=group,r=r) ,data=tt,
                    cols=c("PositiveCorr","NegativeCorr","NonSig")) +
    coord_equal()+
    scale_y_discrete(limits=1:31,breaks=1:31,labels=unique(tt$CancerType))+
    scale_x_discrete(limits=1:9,breaks=1:9,labels=tt$ImmFeature[1:9])+
    labs(title = paste("ImmAPA event in",cluster),group="Related Types")+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_line(linetype = "longdash",color = "gray",size = 0.3),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(color="black",size=10),
            axis.text.x =element_text(color="black",size=10,angle = 90,hjust = 1),
            axis.text.y = element_text(color="black",size=10,face = "bold"),
            axis.title.y = element_blank(),
            #axis.line.y = element_line(colour = "black"),
            #axis.line.x = element_line(colour = "black"),
            panel.border = element_rect(color = "black",fill = NA))+
    scale_fill_manual(values = c("#F44336","#2196F3","#BDBDBD"))
  return(p)
  #readr::write_rds(oncluster,paste("/work/gywang/project/Immune_APA/Data/ImmAPA/541ImmAPA_PDUI_corr_with_ImmFeature_in",cluster,".rds.gz",sep = ""))
}) -> ImmAPA_3cluster_corr_PLOT

readr::write_rds(ImmAPA_3cluster_corr_PLOT,"/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_3cluster_corr_PLOT.rds")


pdf("/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/ImmAPA_541_3_cluster/3cluster APA corr WITH All immune feature.pdf",width = 12,height = 12)
cowplot::plot_grid(plotlist = ImmAPA_3cluster_corr_PLOT,ncol = 3)
dev.off()



#################### only ImmFeature

lapply(unique(All_ImmAPA_ESS_cluster$combined_3cluster), function(cluster){
  purrr::map2(All_score_APA$ImmAPA_Immfeature, All_score_APA$cancer_types,function(x,cans){
    #x <- All_score_APA$ImmAPA_Immfeature[[6]]
    APAevens <- intersect(colnames(x) ,rownames(All_ImmAPA_ESS_cluster[All_ImmAPA_ESS_cluster$combined_3cluster %in% cluster,]))
    x <- dplyr::select(x,1:10,APAevens)
    lapply(2:10, function(pos){
      features <- x[,pos]
      sub <-  apply(x[,-1:-10],2,function(apa){
        corr <- my.cor.test_pearson(features,as.numeric(apa))
        ifelse(as.character( corr[1]) %in%  "NA",return(c(Corr =NA ,P=NA)), return(c(Corr =corr$estimate ,P=corr$p.value )) )
      }) %>% t() %>% as.data.frame()
      sub$fdr <- p.adjust(sub$P,method="fdr")
      sub$sig <- ifelse(sub$fdr > 0.05,"NonSig",ifelse(sub$Corr.cor > 0.2,"PositiveCorr",ifelse(sub$Corr.cor < -0.2,  "NegativeCorr","NonSig" )))
      sub1 <- c(ImmFeature= colnames(x)[pos],NegativeCorr=length(sub[sub$sig %in% "NegativeCorr", ]$sig),NonSig=length(sub[sub$sig %in% "NonSig", ]$sig),PositiveCorr=length(sub[sub$sig %in% "PositiveCorr", ]$sig))
      return(sub1)
    }) %>% bind_rows() -> sigCanCluster
    sigCanCluster$CancerType <- cans
    return(sigCanCluster)
  }) %>% bind_rows() -> oncluster
  oncluster$ImmAPAcluster <- cluster
  tt <- oncluster
  tt$ImmFeature <- factor(tt$ImmFeature,levels= unique(tt$ImmFeature))
  tt$x_posi <- as.numeric(factor(tt$ImmFeature))
  tt$y_posi <- as.numeric(factor(tt$CancerType))
  tt$group <- 1:279
  tt$r <- 0.4
  tt[,2:4] <- apply(tt[2:4],2,as.numeric)
  tt <- dplyr::filter(tt ,ImmFeature %in% unique(tt$ImmFeature)[1:3])
  p <- ggplot(tt)+
    geom_scatterpie(aes(x=x_posi ,y=y_posi,group=group,r=r) ,data=tt,
                    cols=c("PositiveCorr","NegativeCorr","NonSig")) +
    coord_equal()+
    scale_y_discrete(limits=1:31,breaks=1:31,labels=unique(tt$CancerType))+
    scale_x_discrete(limits=1:3,breaks=1:3,labels=tt$ImmFeature[1:3])+
    labs(title = paste("ImmAPA event in",cluster),group="Related Types")+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_line(linetype = "longdash",color = "gray",size = 0.3),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(color="black",size=10),
            axis.text.x =element_text(color="black",size=10,angle = 90,hjust = 1),
            axis.text.y = element_text(color="black",size=10,face = "bold"),
            axis.title.y = element_blank(),
            #axis.line.y = element_line(colour = "black"),
            #axis.line.x = element_line(colour = "black"),
            panel.border = element_rect(color = "black",fill = NA))+
    scale_fill_manual(values = c("#F44336","#2196F3","#BDBDBD"))
  p
  return(p)
  #readr::write_rds(oncluster,paste("/work/gywang/project/Immune_APA/Data/ImmAPA/541ImmAPA_PDUI_corr_with_ImmFeature_in",cluster,".rds.gz",sep = ""))
}) -> ImmAPA_ImmFeature_3cluster_corr_PLOT
lapply(unique(All_ImmAPA_ESS_cluster$combined_3cluster), function(cluster){
  purrr::map2(All_score_APA$ImmAPA_Immfeature, All_score_APA$cancer_types,function(x,cans){
    #x <- All_score_APA$ImmAPA_Immfeature[[6]]
    APAevens <- intersect(colnames(x) ,rownames(All_ImmAPA_ESS_cluster[All_ImmAPA_ESS_cluster$combined_3cluster %in% cluster,]))
    x <- dplyr::select(x,1:10,APAevens)
    lapply(2:10, function(pos){
      features <- x[,pos]
      sub <-  apply(x[,-1:-10],2,function(apa){
        corr <- my.cor.test_pearson(features,as.numeric(apa))
        ifelse(as.character( corr[1]) %in%  "NA",return(c(Corr =NA ,P=NA)), return(c(Corr =corr$estimate ,P=corr$p.value )) )
      }) %>% t() %>% as.data.frame()
      sub$fdr <- p.adjust(sub$P,method="fdr")
      sub$sig <- ifelse(sub$fdr > 0.1,"NonSig",ifelse(sub$Corr.cor > 0.2,"PositiveCorr",ifelse(sub$Corr.cor < -0.2,  "NegativeCorr","NonSig" )))
      sub1 <- c(ImmFeature= colnames(x)[pos],NegativeCorr=length(sub[sub$sig %in% "NegativeCorr", ]$sig),NonSig=length(sub[sub$sig %in% "NonSig", ]$sig),PositiveCorr=length(sub[sub$sig %in% "PositiveCorr", ]$sig))
      return(sub1)
    }) %>% bind_rows() -> sigCanCluster
    sigCanCluster$CancerType <- cans
    return(sigCanCluster)
  }) %>% bind_rows() -> oncluster
  oncluster$ImmAPAcluster <- cluster
  tt <- oncluster
  tt$ImmFeature <- factor(tt$ImmFeature,levels= unique(tt$ImmFeature))
  tt$x_posi <- as.numeric(factor(tt$ImmFeature))
  tt$y_posi <- as.numeric(factor(tt$CancerType))
  tt$group <- 1:279
  tt$r <- 0.4
  tt[,2:4] <- apply(tt[2:4],2,as.numeric)
  tt <- dplyr::filter(tt ,ImmFeature %in% unique(tt$ImmFeature)[1:3])
  p <- ggplot(tt)+
    geom_scatterpie(aes(x=x_posi ,y=y_posi,group=group,r=r) ,data=tt,
                    cols=c("PositiveCorr","NegativeCorr","NonSig")) +
    coord_equal()+
    scale_y_discrete(limits=1:31,breaks=1:31,labels=unique(tt$CancerType))+
    scale_x_discrete(limits=1:3,breaks=1:3,labels=tt$ImmFeature[1:3])+
    labs(title = paste("ImmAPA event in",cluster),group="Related Types")+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_line(linetype = "longdash",color = "gray",size = 0.3),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(color="black",size=10),
            axis.text.x =element_text(color="black",size=10,angle = 90,hjust = 1),
            axis.text.y = element_text(color="black",size=10,face = "bold"),
            axis.title.y = element_blank(),
            #axis.line.y = element_line(colour = "black"),
            #axis.line.x = element_line(colour = "black"),
            panel.border = element_rect(color = "black",fill = NA))+
    scale_fill_manual(values = c("#F44336","#2196F3","#BDBDBD"))
  p
  return(p)
  #readr::write_rds(oncluster,paste("/work/gywang/project/Immune_APA/Data/ImmAPA/541ImmAPA_PDUI_corr_with_ImmFeature_in",cluster,".rds.gz",sep = ""))
}) -> ImmAPA_ImmFeatureFDR0.1_3cluster_corr_PLOT
pdf("/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/ImmAPA_541_3_cluster/3cluster APA Rs0.2 FDR0.05 corr WITH ImmFeature.pdf",width = 12,height = 12)
cowplot::plot_grid(plotlist =ImmAPA_ImmFeature_3cluster_corr_PLOT,ncol = 3)
dev.off()
pdf("/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/ImmAPA_541_3_cluster/3cluster APA Rs0.2 FDR0.1 corr WITH ImmFeature.pdf",width = 12,height = 12)
cowplot::plot_grid(plotlist =ImmAPA_ImmFeatureFDR0.1_3cluster_corr_PLOT,ncol = 3)
dev.off()


#################### only ICI

lapply(unique(All_ImmAPA_ESS_cluster$combined_3cluster), function(cluster){
  purrr::map2(All_score_APA$ImmAPA_Immfeature, All_score_APA$cancer_types,function(x,cans){
    #x <- All_score_APA$ImmAPA_Immfeature[[6]]
    APAevens <- intersect(colnames(x) ,rownames(All_ImmAPA_ESS_cluster[All_ImmAPA_ESS_cluster$combined_3cluster %in% cluster,]))
    x <- dplyr::select(x,1:10,APAevens)
    lapply(2:10, function(pos){
      features <- x[,pos]
      sub <-  apply(x[,-1:-10],2,function(apa){
        corr <- my.cor.test_pearson(features,as.numeric(apa))
        ifelse(as.character( corr[1]) %in%  "NA",return(c(Corr =NA ,P=NA)), return(c(Corr =corr$estimate ,P=corr$p.value )) )
      }) %>% t() %>% as.data.frame()
      sub$fdr <- p.adjust(sub$P,method="fdr")
      sub$sig <- ifelse(sub$fdr > 0.05,"NonSig",ifelse(sub$Corr.cor > 0.2,"PositiveCorr",ifelse(sub$Corr.cor < -0.2,  "NegativeCorr","NonSig" )))
      sub1 <- c(ImmFeature= colnames(x)[pos],NegativeCorr=length(sub[sub$sig %in% "NegativeCorr", ]$sig),NonSig=length(sub[sub$sig %in% "NonSig", ]$sig),PositiveCorr=length(sub[sub$sig %in% "PositiveCorr", ]$sig))
      return(sub1)
    }) %>% bind_rows() -> sigCanCluster
    sigCanCluster$CancerType <- cans
    return(sigCanCluster)
  }) %>% bind_rows() -> oncluster
  oncluster$ImmAPAcluster <- cluster
  tt <- oncluster
  tt$ImmFeature <- factor(tt$ImmFeature,levels= unique(tt$ImmFeature))
  tt$x_posi <- as.numeric(factor(tt$ImmFeature))
  tt$y_posi <- as.numeric(factor(tt$CancerType))
  tt$group <- 1:279
  tt$r <- 0.4
  tt[,2:4] <- apply(tt[2:4],2,as.numeric)
  tt <- dplyr::filter(tt ,ImmFeature %in% unique(tt$ImmFeature)[4:9])
  p <- ggplot(tt)+
    geom_scatterpie(aes(x=x_posi ,y=y_posi,group=group,r=r) ,data=tt,
                    cols=c("PositiveCorr","NegativeCorr","NonSig")) +
    coord_equal()+
    scale_y_discrete(limits=1:31,breaks=1:31,labels=unique(tt$CancerType))+
    scale_x_discrete(limits=4:9,breaks=4:9,labels=tt$ImmFeature[1:6])+
    labs(title = paste("ImmAPA event in",cluster),group="Related Types")+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_line(linetype = "longdash",color = "gray",size = 0.3),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(color="black",size=10),
            axis.text.x =element_text(color="black",size=10,angle = 90,hjust = 1),
            axis.text.y = element_text(color="black",size=10,face = "bold"),
            axis.title.y = element_blank(),
            #axis.line.y = element_line(colour = "black"),
            #axis.line.x = element_line(colour = "black"),
            panel.border = element_rect(color = "black",fill = NA))+
    scale_fill_manual(values = c("#F44336","#2196F3","#BDBDBD"))
  p
  return(p)
  #readr::write_rds(oncluster,paste("/work/gywang/project/Immune_APA/Data/ImmAPA/541ImmAPA_PDUI_corr_with_ImmFeature_in",cluster,".rds.gz",sep = ""))
}) -> ImmAPA_ICI_3cluster_corr_PLOT
pdf("/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/ImmAPA_541_3_cluster/3cluster APA corr Rs0.2 FDR0.05 WITH ICI.pdf",width = 12,height = 12)
cowplot::plot_grid(plotlist =ImmAPA_ICI_3cluster_corr_PLOT,ncol = 3)
dev.off()

lapply(unique(All_ImmAPA_ESS_cluster$combined_3cluster), function(cluster){
  purrr::map2(All_score_APA$ImmAPA_Immfeature, All_score_APA$cancer_types,function(x,cans){
    #x <- All_score_APA$ImmAPA_Immfeature[[6]]
    APAevens <- intersect(colnames(x) ,rownames(All_ImmAPA_ESS_cluster[All_ImmAPA_ESS_cluster$combined_3cluster %in% cluster,]))
    x <- dplyr::select(x,1:10,APAevens)
    lapply(2:10, function(pos){
      features <- x[,pos]
      sub <-  apply(x[,-1:-10],2,function(apa){
        corr <- my.cor.test_pearson(features,as.numeric(apa))
        ifelse(as.character( corr[1]) %in%  "NA",return(c(Corr =NA ,P=NA)), return(c(Corr =corr$estimate ,P=corr$p.value )) )
      }) %>% t() %>% as.data.frame()
      sub$fdr <- p.adjust(sub$P,method="fdr")
      sub$sig <- ifelse(sub$fdr > 0.1,"NonSig",ifelse(sub$Corr.cor > 0.2,"PositiveCorr",ifelse(sub$Corr.cor < -0.2,  "NegativeCorr","NonSig" )))
      sub1 <- c(ImmFeature= colnames(x)[pos],NegativeCorr=length(sub[sub$sig %in% "NegativeCorr", ]$sig),NonSig=length(sub[sub$sig %in% "NonSig", ]$sig),PositiveCorr=length(sub[sub$sig %in% "PositiveCorr", ]$sig))
      return(sub1)
    }) %>% bind_rows() -> sigCanCluster
    sigCanCluster$CancerType <- cans
    return(sigCanCluster)
  }) %>% bind_rows() -> oncluster
  oncluster$ImmAPAcluster <- cluster
  tt <- oncluster
  tt$ImmFeature <- factor(tt$ImmFeature,levels= unique(tt$ImmFeature))
  tt$x_posi <- as.numeric(factor(tt$ImmFeature))
  tt$y_posi <- as.numeric(factor(tt$CancerType))
  tt$group <- 1:279
  tt$r <- 0.4
  tt[,2:4] <- apply(tt[2:4],2,as.numeric)
  tt <- dplyr::filter(tt ,ImmFeature %in% unique(tt$ImmFeature)[4:9])
  p <- ggplot(tt)+
    geom_scatterpie(aes(x=x_posi ,y=y_posi,group=group,r=r) ,data=tt,
                    cols=c("PositiveCorr","NegativeCorr","NonSig")) +
    coord_equal()+
    scale_y_discrete(limits=1:31,breaks=1:31,labels=unique(tt$CancerType))+
    scale_x_discrete(limits=4:9,breaks=4:9,labels=tt$ImmFeature[1:6])+
    labs(title = paste("ImmAPA event in",cluster),group="Related Types")+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_line(linetype = "longdash",color = "gray",size = 0.3),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(color="black",size=10),
            axis.text.x =element_text(color="black",size=10,angle = 90,hjust = 1),
            axis.text.y = element_text(color="black",size=10,face = "bold"),
            axis.title.y = element_blank(),
            #axis.line.y = element_line(colour = "black"),
            #axis.line.x = element_line(colour = "black"),
            panel.border = element_rect(color = "black",fill = NA))+
    scale_fill_manual(values = c("#F44336","#2196F3","#BDBDBD"))
  p
  return(p)
  #readr::write_rds(oncluster,paste("/work/gywang/project/Immune_APA/Data/ImmAPA/541ImmAPA_PDUI_corr_with_ImmFeature_in",cluster,".rds.gz",sep = ""))
}) -> ImmAPA_ICI_3cluster_FDR0.1_corr_PLOT



pdf("/work/gywang/project/Immune_APA/Figure/ImmAPA_Clustering/ImmAPA_541_3_cluster/3cluster APA corr Rs0.2 FDR0.1 WITH ICI.pdf",width = 12,height = 12)
cowplot::plot_grid(plotlist =ImmAPA_ICI_3cluster_FDR0.1_corr_PLOT,ncol = 3)
dev.off()



