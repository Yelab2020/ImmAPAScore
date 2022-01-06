library(tidyverse)
library(stats)

options(stringsAsFactors = F)

my.cor.test_pearson <- function(...) {
  obj<- try(cor.test(...,method="pearson"), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}
my.cor.test_spearman <- function(...) {
  obj<- try(cor.test(...,method="spearman"), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}

GTEx_data <- readr::read_rds("/work/gywang/project/APA_LiLab/APA_change_AGING_GTEx/data/GTEx_FPKM_PDUI_all_data.rds")
TCGA_33canFPKM <- readr::read_rds("/data/share_data/TCGA_preliminary/TCGA_33CanFPKM_AllSamples.rds.gz")
TCGA_APA_YUXIANG <- readr::read_rds("/data/share_data/TCGA_APA/TCGA_APA_Yuxiang.rds.gz")
TCGA_APA_YUxiang <- TCGA_APA_YUXIANG[c(1:3,5:8,10:12,14,16,17),]
TCGA_APA_YUxiang$cancer_types <- TCGA_APA_YUxiang$cancer_types %>% str_sub(5)
Sub_TCGA_FPKM <- TCGA_33canFPKM %>% filter( cancertypes %in% TCGA_APA_YUxiang$cancer_types)

ImmAPA_541Ranking <- readr::read_rds("/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_ranking_n.rds.gz")
Top10ImmAPA <- ImmAPA_541Ranking %>% top_n(10,Rank) %>% rownames()

CandidateAPAevent <- Top10ImmAPA[2]
CandidateGeneName <- str_split(CandidateAPAevent,"\\|")[[1]][2]
Candidate_APA_PDUI <- filter(GTEx_PDUI, gene_name %in% CandidateGeneName) %>% t() %>% as.data.frame()
Candidate_APA_PDUI$SampleName <- rownames(Candidate_APA_PDUI)
Candidate_APA_PDUI <- Candidate_APA_PDUI[-1:-3,]
colnames(Candidate_APA_PDUI)[1] <- "PDUI"


GTEx_TCGA <- tibble(GTEx_tissue=list("Bladder","Breast",c("Cervix_ecto","Cervix_endo"),NULL,"Kidney","Kidney","Kidney","Liver","Lung","Lung","Prostate" ,"Stomach"   ,"Thyroid"),
                    TCGA_cancer=list("BLCA", "BRCA", "CESC","HNSC","KICH" ,"KIRC" ,"KIRP" ,"LIHC" ,"LUAD" ,"LUSC", "PRAD" ,"STAD", "THCA"))



GTEx_TCGA <- GTEx_TCGA %>% mutate(Candidatate_GTEx= lapply(GTEx_TCGA$GTEx_tissue, function(x){
  sub <- GTEx_data %>% filter(tissue %in% x)
  sub_APA <- sub$PDUI
  sub_FPKM <- sub$FPKM_log2
  PDUI_FPKM <-  purrr::map2(sub_APA,sub_FPKM,function(x1,y1){
    sub <- x1 %>% filter(gene_name %in% CandidateGeneName)  %>% dplyr::select(-1:-3) %>% t() %>% as.data.frame()
    sub$Tissue <- unlist(str_split(x,"_") )[1]
    sub$SampleName <- rownames(sub)
    colnames(sub)[1] <- "PDUI"
    #rownames(sub) <- NULL
    sub_FPKM <- filter(y1,hgnc_symbol %in% CandidateGeneName)  %>% dplyr::select(-1:-3) %>% t() %>% as.data.frame()
    sub_FPKM$SampleName  <- sub_FPKM %>% rownames() 
    #sub_FPKM <- sub_FPKM[-1,]
    sub_FPKM$log2FPKM <- sub_FPKM$V1 %>% as.character() %>% as.numeric()
    sub_FPKM <- sub_FPKM %>% dplyr::select(-V1)
    sub_PDUI_FPKM <- merge(sub,sub_FPKM,by.x="SampleName",by.y="SampleName")
    sub_PDUI_FPKM$source <- "GTEx"
    return(sub_PDUI_FPKM)})   %>% bind_rows()
    PDUI_FPKM$GeneName <- CandidateGeneName
    return(PDUI_FPKM)
  }))


GTEx_TCGA <- GTEx_TCGA %>% mutate(Candidatate_TCGA=lapply(GTEx_TCGA$TCGA_cancer, function(x){
  
  #x <- GTEx_TCGA$TCGA_cancer[[1]]
  sub_APA <- TCGA_APA_YUxiang %>% filter(cancer_types %in% x)
  sub_FPKM <- Sub_TCGA_FPKM %>% filter(cancertypes %in% x)
  PDUI_FPKM <-  purrr::map2(sub_APA$APAexpr,sub_FPKM$Expr,function(x1,y1){
    #x1 <- sub_APA$APAexpr[[1]]
    #y1 <- sub_FPKM$Expr[[1]]
    sub <- x1[CandidateAPAevent,] %>% dplyr::select(ends_with("PDUI")) %>% t() %>% as.data.frame()
    sub$Is_tumor <- str_sub(rownames(sub),1,5)
    sub$SampleName <- str_sub(rownames(sub),7,-6)
    sub$SampleName <- ifelse(str_sub(sub$SampleName,1,1) %in% ".",str_sub(sub$SampleName,2),str_sub(sub$SampleName,1))
    sub$SampleName <- ifelse(sub$Is_tumor %in% "Tumor",paste(sub$SampleName,"_Tumor",sep = ""),paste(sub$SampleName,"_Normal",sep = ""))
    sub <- dplyr::select(sub,- Is_tumor)
    colnames(sub)[1] <- "PDUI"
    sub_FPKM <- filter(y1,Gene %in% CandidateGeneName) %>% t() %>% as.data.frame()
    sub_FPKM$SampleName_Orig  <- sub_FPKM %>% rownames()
    sub_FPKM <- sub_FPKM[-1,]
    sub_FPKM$SampleName  <- sub_FPKM$SampleName_Orig %>% str_sub(6,12)
    # sub_FPKM$SampleName <- gsub("\\-",".",sub_FPKM$SampleName)
    sub_FPKM$Is_tumor <- ifelse(str_sub(sub_FPKM$SampleName_Orig,14,14) %in% "1", "Normal","Tumor")
    sub_FPKM$SampleName <- ifelse(sub_FPKM$Is_tumor %in% "Tumor",paste(sub_FPKM$SampleName,"_Tumor",sep = ""),paste(sub_FPKM$SampleName,"_Normal",sep = ""))
    sub_FPKM  <-  split(sub_FPKM,sub_FPKM$SampleName) %>% lapply(function(x){
      if(nrow(x) >1){
        x$V1 <- mean(as.numeric(as.character(x$V1)))
        return(x[1,-2])
      }else{ x[1] <- as.numeric(as.character(x[1]))
        return(x[-2])}
    }) %>% bind_rows()
    sub_FPKM$log2FPKM <-  log2(as.numeric(as.character(sub_FPKM$V1)))
    sub_FPKM <- sub_FPKM %>% dplyr::select(-V1)
    sub_PDUI_FPKM <- merge(sub,sub_FPKM,by.x="SampleName",by.y="SampleName")  
    sub_PDUI_FPKM$source <- "TCGA"
    return(sub_PDUI_FPKM)})  
  PDUI_FPKM <- purrr::map2(PDUI_FPKM,sub_APA$cancer_types,function(x,y){
    x$TCGA_cancer <- y
    x$GeneName <- CandidateGeneName
    return(x)
  }) %>% bind_rows()
  return(PDUI_FPKM)
  }))

GTEx_TCGA  <- GTEx_TCGA %>% mutate(Candidatate_GTEx_TCGA = purrr::map2(Candidatate_GTEx,Candidatate_TCGA,function(x,y){
  x$Is_tumor <- "Normal"
  colnames(y)[6] <- "Tissue"
  sub <- rbind(x,y)
  sub$group <- ifelse(sub$source %in% "GTEx","GTEx_N",ifelse(sub$Is_tumor %in% "Tumor","TCGA_T","TCGA_N"))
  return(sub)
  }))

GTEx_TCGA <- GTEx_TCGA  %>% mutate(PDUI_FPKM_Corr= lapply(Candidatate_GTEx_TCGA, function(x){
  corr_pearson <- split(x,x$group) %>% lapply(function(x1){
    sub <- my.cor.test_pearson(x1$PDUI,x1$log2FPKM)
    sub<- c(sub,SampleNumber=nrow(x1))
    return(sub)
  }) %>% unlist()
  corr_spearman <- split(x,x$group) %>% lapply(function(x1){
    sub <- my.cor.test_spearman(x1$PDUI,x1$log2FPKM)
    sub<- c(sub,SampleNumber=nrow(x1))
    return(sub)
  }) %>% unlist()
  corr <- data.frame(corr_pearson,corr_spearman)
  return(corr)
}))

#GTEx_TCGA <- GTEx_TCGA_COL1A1
GTEx_TCGA <- GTEx_TCGA %>% mutate(PDUI_FPKM_Corr_pearson_plot =  purrr::map2(GTEx_TCGA$Candidatate_GTEx_TCGA,GTEx_TCGA$PDUI_FPKM_Corr,function(data,corr){
  tissue <- unique(data[data$source %in% "GTEx",]$Tissue)
  p <-  ggplot(data,aes(x=PDUI,y=log2FPKM,group=group))+
    geom_point(aes(color=group))+
    ylim(c(0,15))+
    labs(title = paste("Correlation betwwen PDUI and FPKM of", CandidateGeneName , "in",tissue,"and",data[data$source %in% "TCGA",]$Tissue),x="PDUI",y="log2(FPKM)",
         subtitle = paste( "   GTEx                              ","Estimate Corr:",round(as.numeric(corr$corr_pearson[1]),3),"  ","P value:",round(as.numeric(corr$corr_pearson[2]),4),"  ", "Samples size:",corr$corr_pearson[3],"\n",
                           "Normal tissue in TCGA    ","Estimate Corr:",round(as.numeric(corr$corr_pearson[4]),3),"    ","P value:",round(as.numeric(corr$corr_pearson[5]),4),"  ", "Samples size:",corr$corr_pearson[6],"\n",
                           "Tumor tissue in TCGA     ","Estimate Corr:",round(as.numeric(corr$corr_pearson[7]),3),"   ","P value:",round(as.numeric(corr$corr_pearson[8]),4),"     ", "Samples size:",corr$corr_pearson[9]))+
    geom_smooth(method = 'lm', formula = y ~ x,aes(color=group)) +
    #stat_compare_means(method = "anova")+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(color="black",size=12),
            axis.text.x =element_text(color="black",size=10,angle = 0,hjust = 1),
            axis.text.y = element_text(color="black",size=10,face = "bold"),
            axis.line.y = element_line(colour = "black"),
            axis.line.x = element_line(colour = "black"),
            panel.border = element_rect(color = "black",fill = NA))+
    scale_color_manual(limit=c("GTEx_N","TCGA_N","TCGA_T"),values =c("green","blue","red"))
}))

GTEx_TCGA <- GTEx_TCGA %>% mutate(PDUI_FPKM_Corr_spearman_plot =  purrr::map2(GTEx_TCGA$Candidatate_GTEx_TCGA,GTEx_TCGA$PDUI_FPKM_Corr,function(data,corr){
  tissue <- unique(data[data$source %in% "GTEx",]$Tissue)
  p <-  ggplot(data,aes(x=PDUI,y=log2FPKM,group=group))+
    geom_point(aes(color=group))+
    ylim(c(0,15))+
    labs(title = paste("Correlation betwwen PDUI and FPKM of", CandidateGeneName , "in",tissue,"and",data[data$source %in% "TCGA",]$Tissue),x="PDUI",y="log2(FPKM)",
         subtitle = paste( "   GTEx                              ","Estimate Corr:",round(as.numeric(corr$corr_spearman[1]),3),"  ","P value:",round(as.numeric(corr$corr_spearman[2]),4),"  ", "Samples size:",corr$corr_spearman[3],"\n",
                           "Normal tissue in TCGA    ","Estimate Corr:",round(as.numeric(corr$corr_spearman[4]),3),"    ","P value:",round(as.numeric(corr$corr_spearman[5]),4),"  ", "Samples size:",corr$corr_spearman[6],"\n",
                           "Tumor tissue in TCGA     ","Estimate Corr:",round(as.numeric(corr$corr_spearman[7]),3),"   ","P value:",round(as.numeric(corr$corr_spearman[8]),4),"     ", "Samples size:",corr$corr_spearman[9]))+
    geom_smooth(method = 'lm', formula = y ~ x,aes(color=group)) +
    #stat_compare_means(method = "anova")+
    theme(  panel.background=element_rect(colour=NA,fill="white"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(color="black",size=12),
            axis.text.x =element_text(color="black",size=10,angle = 0,hjust = 1),
            axis.text.y = element_text(color="black",size=10,face = "bold"),
            axis.line.y = element_line(colour = "black"),
            axis.line.x = element_line(colour = "black"),
            panel.border = element_rect(color = "black",fill = NA))+
    scale_color_manual(limit=c("GTEx_N","TCGA_N","TCGA_T"),values =c("green","blue","red"))
}))

readr::write_rds(GTEx_TCGA,paste("/work/gywang/project/Immune_APA/Data/GTEx_TCGA",CandidateGeneName,"_GTEx_TCGA_FPKM_PDUI.rds.gz",sep = ""),compress = "gz")


GTEx_TCGA <- readr::read_rds("/work/gywang/project/Immune_APA/Data/GTEx_TCGA/COL1A1_GTEx_TCGA_FPKM_PDUI.rds.gz")
###################### PLOT  

All_tissue_corr <-  purrr::map2(GTEx_TCGA$Candidatate_GTEx_TCGA,GTEx_TCGA$PDUI_FPKM_Corr,function(data,corr){
  tissue <- unique(data[data$source %in% "TCGA",]$Tissue)
  corr_N <- corr[c(1,4,7),]
  corr_N$corr_pearson_Pvalue <- corr[c(2,5,8),1]
  corr_N$corr_spearman_Pvalue <- corr[c(2,5,8),2]
  corr_N$group <- rownames(corr_N) %>% str_sub(1,6)
  corr_N$tissue <- tissue
  corr_N
  }) %>% bind_rows()
CandidateGeneName <-"COL1A1"
  
All_tissue_corr$corr_pearson_sig <- ifelse(All_tissue_corr$corr_pearson_Pvalue < 0.05,"Sig","NonSig")
All_tissue_corr$corr_spearman_sig <- ifelse(All_tissue_corr$corr_spearman_Pvalue < 0.05,"Sig","NonSig")
All_tissue_corr <- All_tissue_corr %>% na.omit()
ggplot(All_tissue_corr)+
  geom_point(aes(x=tissue,y=corr_pearson,group=group,color=group,size=-log10(corr_pearson_Pvalue+10^-10),shape=corr_pearson_sig))+
  geom_hline(yintercept = 0,size=0.5,,color="gray",linetype = c("longdash"))+
  labs(title = paste("Pearson Correlation between PDUI and FPKM of", CandidateGeneName,"in TCGA and GTEx"),size= "-log10(Pvalue)",y=paste("Correlation of PDUI and FPKM of",CandidateGeneName),shape="Significance")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_text(color="black",size=12),
          axis.text.x =element_text(color="black",size=10,angle = 0,hjust = 0.5),
          axis.text.y = element_text(color="black",size=10,face = "bold"),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          panel.border = element_rect(color = "black",fill = NA))+
  scale_color_manual(limit=c("GTEx_N","TCGA_N","TCGA_T"),values =c("green","blue","red"))+
  scale_size_continuous(limit=c(0,10),range = c(1,6),breaks=c(-log10(0.05),5,10),labels=c("0.05","1e-5","<1e-10"),name="Pvalue")+
  scale_shape_manual(values = c(1,19))+
  scale_linetype_manual(values = c("longdash"))


ggplot(All_tissue_corr)+
  geom_point(aes(x=tissue,y=corr_spearman,group=group,color=group,size=-log10(corr_spearman_Pvalue+10^-10),shape=corr_spearman_sig))+
  geom_hline(yintercept = 0,size=0.5,,color="gray",linetype = c("longdash"))+
  labs(title = paste("Spearman Correlation between PDUI and FPKM of", CandidateGeneName,"in TCGA and GTEx"),size= "-log10(Pvalue)",y=paste("Correlation of PDUI and FPKM of",CandidateGeneName),shape="Significance")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_text(color="black",size=12),
          axis.text.x =element_text(color="black",size=10,angle = 0,hjust = 0.5),
          axis.text.y = element_text(color="black",size=10,face = "bold"),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          panel.border = element_rect(color = "black",fill = NA))+
  scale_color_manual(limit=c("GTEx_N","TCGA_N","TCGA_T"),values =c("green","blue","red"))+
  scale_size_continuous(limit=c(0,10),range = c(1,6),breaks=c(-log10(0.05),5,10),labels=c("0.05","1e-5","<1e-10"),name="Pvalue")+
  scale_shape_manual(values = c(1,19))+
  scale_linetype_manual(values = c("longdash"))
p2
pdf(paste("/work/gywang/project/Immune_APA/Figure/COL1A2/Correlation score of PDUI and FPKM of", CandidateGeneName,"in TCGA and GTEx.pdf"),width = 20,height = 7)
cowplot::plot_grid(plotlist = list(p1,p2))
dev.off()



GTEx_N <- All_tissue_corr %>% dplyr::filter(group %in% "GTEx_N")
TCGA_N <- All_tissue_corr %>% dplyr::filter(group %in% "TCGA_N")
TCGA_T <- All_tissue_corr %>% dplyr::filter(group %in% "TCGA_T")


All_tissue_corr %>% split(All_tissue_corr$group)  %>% lapply(function(x){
  median(x$corr_pearson)
}) %>% bind_rows()

openxlsx::write.xlsx(All_tissue_corr,"/work/gywang/project/Immune_APA/supplementary_table/COL1A1_PDUI_FPKM_corr.xlsx")


