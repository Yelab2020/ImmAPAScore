library(ggpubr)

ImmAPA_541_ImmRanking <- readr::read_rds("/work/gywang/project/Immune_APA/Data/Immune_pathway_data/ImmAPA_ranking/ImmAPA_541_ImmAPA_ranking_method_NC.rds")
Top20ImmAPA <- ImmAPA_541_ImmRanking %>% top_n(20,Rank) %>% rownames()# %>% lapply(as.list)


ImmPathways <- readr::read_rds("/home/yye/Project/2020APA_Imm/ImmunePathwayList.rds")
ImmPathwaysD <- data.frame(term=rep(names(ImmPathways),times=unlist(lapply(ImmPathways,length))),gene=unlist(ImmPathways))


#TC3a <- readr::read_rds("/data/share_data/TCGA_APA/TCGA_APA_tc3a.rds.gz")

TCGA_APA_YUxiang <- readr::read_rds("/data/share_data/TCGA_APA/TCGA_APA_Yuxiang.rds.gz")
TCGA_APA_YUxiang <- TCGA_APA_YUxiang[c(1:3,5:8,10:12,14,16,17),]
TCGA_APA_YUxiang$cancer_types %>% length()


Plist <- list()
for (i in 1:20) {
SigAPA <- Top20ImmAPA[[i]]
sub_APA <- purrr::map2(TCGA_APA_YUxiang$APAexpr,TCGA_APA_YUxiang$cancer_types,function(x,y){
  sub <- x[SigAPA,] %>% dplyr::select(ends_with("PDUI")) %>% t() %>% as.data.frame()
  sub$Tissue <- str_sub(rownames(sub),1,5)
  sub$CancerType <- gsub("APA_","",y)
  colnames(sub)[1] <- "PDUI"
  rownames(sub) <- NULL
  return(sub)
})

sub_APA1 <- sub_APA %>% dplyr::bind_rows()
sub_APA1$Tissue  <- ifelse(sub_APA1$Tissue %in% "Norma","Normal","Tumor")
p1 <- ggplot(sub_APA1, aes(x = CancerType , y = PDUI))+ 
  geom_boxplot(aes(fill = Tissue),position=position_dodge(0.5),width=0.40)+
  labs(title = paste("PDUI of",str_split(SigAPA,"\\|")[[1]][2]))+
  scale_fill_manual(values = c("gray", "red"))+
  theme(axis.text.x = element_text(size = 8,angle = 90, hjust = 0.5,vjust =0.5),
        axis.text.y = element_text(size = 8),
        axis.title.y = element_blank(),                     #设定y轴的标题
        axis.title.x = element_blank(),
        title = element_text(size = 8),
        #axis.ticks.x = element_blank(),                     #设定x轴的轴须
        #axis.ticks.y = element_blank(),                     #设定y轴的轴须
        #legend.key.size = 13,
        #legend.position='none',
        legend.key.size=unit(1,'cm'),
        legend.key.width=unit(0.5,'cm'),
        legend.text = element_text(size=8,colour = 'black', hjust = 0.5, vjust = 0.5, face = 'bold'),
        legend.title = element_text(size = 8),
        plot.title = element_text(hjust = 0,size=8),       #设定图片的标题
        panel.grid.major = element_blank(),                 #设定主要网格线
        panel.grid.minor = element_blank(),                 #设定次要网格线
        plot.background = element_rect(fill = "white",colour = "white"),   #设定图像的背景  
        panel.background = element_rect(fill = "white",colour = "white"),  #设定数据图的背景
        panel.border = element_rect(fill="transparent",colour  ="black"),                      #图片边框
  )+
  stat_compare_means(aes(group=Tissue),label = "p.signif")  
p1
Plist[[i]] <- p1
}
pdf("/work/gywang/project/Immune_APA/Figure/Top20_ImmAPA/SigAPA_541_ImmAPA(Top20 NC'method) PDUI in 12 tumor normal pairs.pdf",width = 25,height = 15)
cowplot::plot_grid(plotlist = Plist,ncol =5)
dev.off()




lapply(Top20ImmAPA, function(x){
  #x <- Top20ImmAPA[13]
  SigAPA <- x
  APA_events <- SigAPA %>% str_split("\\|")
  APA_events <- APA_events[[1]][2]
  sub_APA <- purrr::map2(TCGA_APA_YUxiang$APAexpr,TCGA_APA_YUxiang$cancer_types,function(x,y){
    sub <- x[SigAPA,] %>% dplyr::select(ends_with("PDUI")) %>% t() %>% as.data.frame()
    sub$Tissue <- str_sub(rownames(sub),1,5)
    sub$CancerType <- gsub("APA_","",y)
    colnames(sub)[1] <- "PDUI"
    rownames(sub) <- NULL
    sub_APA1 <- sub
    sub_APA1$Tissue  <- ifelse(sub_APA1$Tissue %in% "Norma","Normal","Tumor")
    Tumor_mean =mean(na.omit(sub_APA1[sub_APA1$Tissue %in% "Tumor",]$PDUI))
    if(length(na.omit(sub_APA1[sub_APA1$Tissue %in% "Normal",]$PDUI)) < 1){
      Normal_mean = NA
      tetsP <- NA
    }else{
      Normal_mean = mean(na.omit(sub_APA1[sub_APA1$Tissue %in% "Normal",]$PDUI))
      tets= t.test(sub_APA1[sub_APA1$Tissue %in% "Tumor",]$PDUI, sub_APA1[sub_APA1$Tissue %in% "Normal",]$PDUI)
      tetsP <-  tets$p.value
    }
    sub_diff <- c(APA_events=APA_events, cancerType=unique(sub_APA1$CancerType),Normal_mean=Normal_mean,Tumor_mean=Tumor_mean,diff=Tumor_mean - Normal_mean,P=tetsP,log10P=log10(tetsP))
    return(sub_diff)
  })  %>% bind_rows()
  
  sub_APA$mean_dif <- mean(na.omit(as.numeric(sub_APA$diff)))
  return(sub_APA)
}) %>% bind_rows() -> top20APA_diff_N_T

top20APA_diff_N_T[3:8] <- apply(top20APA_diff_N_T[3:8],2,as.numeric)
top20APA_diff_N_T$log10P_n <- ifelse(log10(top20APA_diff_N_T$P)< -10,10,-log10(top20APA_diff_N_T$P))
top20APA_diff_N_T$sig <- ifelse(log10(top20APA_diff_N_T$P) %in% NA ,"NonSig",ifelse(log10(top20APA_diff_N_T$P) > log10(0.05),"NonSig","Sig"))
 level <- unique(top20APA_diff_N_T[order(top20APA_diff_N_T$mean_dif),]$APA_events)
top20APA_diff_N_T$APA_events <- factor(top20APA_diff_N_T$APA_events,levels=c(level))
pdf("/work/gywang/project/Immune_APA/Figure/Top20_ImmAPA/PDUI of top20 immAPA pdui in TUMOR AND NORAML.pdf",width = 8,height = 5)
ggplot(top20APA_diff_N_T,aes(x=APA_events,y=cancerType))+
  geom_point(aes(color=diff,size=log10P_n,shape=sig))+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          #panel.grid.minor=element_blank(),
          axis.title.x=element_text(color="black",size=12),
          axis.text.x =element_text(color="black",size=10,angle = 90,hjust = 1,vjust = 0.5),
          axis.text.y = element_text(color="black",size=10,face = "bold"),
          #axis.line.y = element_line(colour = "black"),
          #axis.line.x = element_line(colour = "black"),
          panel.border = element_rect(color = "black",fill = NA))+
  scale_size_continuous(limit=c(0,10),range = c(1,6),breaks=c(-log10(0.05),5,10),labels=c("0.05","1e-5","<1e-10"),name="P")+
  scale_shape_manual(values = c(1,19))+
  scale_color_gradient2(low="blue",high="red",mid = "gray")
dev.off()


