library(tidyverse)
library(magrittr)
library(pheatmap)
library(Rtsne)
library(scatterpie)
library(stats)
library(survival)
library(ggpubr)
options(stringsAsFactors = F)

my.wilcox.test <- function(...) {
  obj<-  try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.cor.test_spearman <- function(...) {
  obj<- try(cor.test(...,method="spearman"), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}
ImmAPA_541_ImmRanking <- readr::read_rds("/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_ranking_n.rds.gz")
Top10ImmAPA <- ImmAPA_541_ImmRanking %>% top_n(10,Rank) %>% rownames()# %>% lapply(as.list)

CYT_PRF1_GZMA <- readr::read_rds("/work/yye/AliRstudio/Project/2018Immunotherapy/ProcessedData/CYT_PRF1_GZMA.rds")
TIFScore <- readr::read_rds("/work/yye/AliRstudio/Project/2018Immunotherapy/ProcessedData/TIFScore.rds.gz")
GEPScore <- readr::read_rds("/work/yye/AliRstudio/Project/2018Immunotherapy/ProcessedData/GEPScore.rds.gz")



#ExprFPKM <- readr::read_rds("/data/share_data/TCGA_preliminary/TCGA_33CanFPKM.rds.gz")
TCGA_APA <- readRDS("/work/gywang/project/Immune_APA/Data/TCGA_Data/pancan32_APA.DATA.rds.gz")
TCGA_APA$APA <- lapply(TCGA_APA$APA,function(x){
  x <- x %>% dplyr::select(-2:-3) #%>% set_colnames(c("APA.events",str_sub(colnames(x)[-1:-3],1,12)))
})

Candidate <- Top10ImmAPA[1]
ClinicalData <- read.delim("/data/share_data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
ClinicalData$PFI.time <- as.numeric(ClinicalData$PFI.time)
ClinicalData$OS.time <- as.numeric(ClinicalData$OS.time)
ClinicalData <- filter(ClinicalData,type!="READ")
ClinicalData$PFI <- as.numeric(ClinicalData$PFI)
ClinicalData$OS <- as.numeric(ClinicalData$OS)

TCGA_APA %>% dplyr::mutate(APA=purrr::map(.x=APA,function(.x){
  .x <- .x %>% dplyr::filter(event_id %in% Candidate) %>% tidyr::gather(key="barcode",value="PDUI",-event_id)
  return(.x)
})) %>% dplyr::mutate(APA = purrr::map2(.x=APA,.y=cancer_types,function(.x,.y){
  if(.y=="SKCM"){
    .x <- .x %>% dplyr::filter(substr(barcode,14,15) %in% "06") %>% dplyr::mutate(barcode=substr(gsub("\\.","-",barcode),1,12))
  }else{
    .x <- .x %>% dplyr::filter(substr(barcode,14,15) %in% c("01","03")) %>% dplyr::mutate(barcode=substr(gsub("\\.","-",barcode),1,12))
  }
  return(.x)}
)) -> CandidateAPA

CandidateAPA <- CandidateAPA[unlist(lapply(CandidateAPA$APA,nrow))>0,] %>% dplyr::mutate(Candidate_sur=purrr::map(.x=APA,function(.x){
  #.x <- CandidateAPA$APA[[1]]
  .x <- .x %>% #dplyr::filter(substr(barcode,14,15) %in% c("01","03","06")) %>%
    dplyr::mutate(barcode = gsub("\\.","-",substr(barcode,1,12))) %>% 
    dplyr::inner_join(ClinicalData,by=c("barcode"="bcr_patient_barcode"))
  .x <- .x %>% dplyr::mutate(OS.time=as.numeric(OS.time)/30,PFI.time=as.numeric(PFI.time)/30) %>% dplyr::filter(!is.na(OS.time))
  .x <- .x %>% dplyr::filter(!is.na(PDUI)) %>%
    dplyr::mutate(PDUIGroup = ifelse(PDUI >= median(PDUI),"Lengthen","Shortening"),
                  OS = as.numeric(as.character(OS)),
                  PFI = as.numeric(as.character(PFI)))
  aa_try <- try(aa <- maxstat.test(survival::Surv(OS.time, OS) ~ PDUI, data= .x,smethod="LogRank"),silent=TRUE)
  if (is(aa_try , "try-error")) {.x$PDUIGroupB <- NA}
  else { .x <- .x %>% dplyr::mutate(PDUIGroupB = ifelse(PDUI > aa_try$estimate,"Lengthen","Shortening"))
    print(nrow(.x))}
  return(.x)
  }))

################## CYT
COL1A1_CYT <- merge(CandidateAPA, CYT_PRF1_GZMA,by = "cancer_types")

COL1A1_CYT <- COL1A1_CYT %>% mutate(COL1A1_pdui_CYT= purrr::map2(COL1A1_CYT$CYT,COL1A1_CYT$Candidate_sur,function(x,y){
  x$new_id <-  gsub("\\.","\\-",str_sub(x$barcode,1,12))
  sub <- merge(x,y,by.x = "new_id",by.y = "barcode")
  sub <- dplyr::select(sub,new_id,CYT, barcode,event_id,PDUI,type,PDUIGroup,PDUIGroupB)
  return(sub)
}))

bind_cancer_CYT <- COL1A1_CYT$COL1A1_pdui_CYT %>% bind_rows()

pdf("/work/gywang/project/Immune_APA/Figure/COL1A1/COL1A1_PDUI GROUP CYT(no line).pdf",width = 10,height = 4)
ggplot(bind_cancer_CYT,aes(x=type,y=CYT,fill=PDUIGroup))+
  geom_boxplot(width=0.7)+
  stat_compare_means(aes(group=PDUIGroup),label="p.signif")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          #panel.grid.major=element_line(linetype = "longdash",color = "gray",size = 0.3),
          panel.grid.minor=element_blank(),
          axis.title.x=element_text(color="black",size=13),
          axis.text.x =element_text(color="black",size=12,angle =90,hjust = 1),
          axis.text.y = element_text(color="black",size=12,face = "bold"),
          #axis.title.y = element_blank(),
          #axis.line.y = element_line(colour = "black"),
          #axis.line.x = element_line(colour = "black"),
          panel.border = element_rect(color = "black",fill = NA))
 
dev.off()
##################TIF
COL1A1_TIF <- merge(CandidateAPA, TIFScore,by = "cancer_types")
COL1A1_TIF <- COL1A1_TIF %>% mutate(COL1A1_pdui_TIF= purrr::map2(COL1A1_TIF$TIFscore,COL1A1_TIF$Candidate_sur,function(x,y){
  x$new_id <-  gsub("\\.","\\-",str_sub(x$barcode,1,12))
  sub <- merge(x,y,by.x = "new_id",by.y = "barcode")
  sub <- dplyr::select(sub,new_id,TIF, barcode,event_id,PDUI,type,PDUIGroup,PDUIGroupB)
  return(sub)
}))
bind_cancer_TIF <- COL1A1_TIF$COL1A1_pdui_TIF %>% bind_rows()
pdf("/work/gywang/project/Immune_APA/Figure/COL1A1/COL1A1_PDUI GROUP TIF(no line).pdf",width = 10,height = 4)
ggplot(bind_cancer_TIF,aes(x=type,y=TIF,fill=PDUIGroup))+
  geom_boxplot(width=0.7)+
  stat_compare_means(aes(group=PDUIGroup),label="p.signif")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          #panel.grid.major=element_line(linetype = "longdash",color = "gray",size = 0.3),
          panel.grid.minor=element_blank(),
          axis.title.x=element_text(color="black",size=13),
          axis.text.x =element_text(color="black",size=12,angle =90,hjust = 1),
          axis.text.y = element_text(color="black",size=12,face = "bold"),
          #axis.title.y = element_blank(),
          #axis.line.y = element_line(colour = "black"),
          #axis.line.x = element_line(colour = "black"),
          panel.border = element_rect(color = "black",fill = NA))
dev.off()

##################GEP
COL1A1_GEP <- merge(CandidateAPA, GEPScore,by = "cancer_types")
COL1A1_GEP <- COL1A1_GEP %>% mutate(COL1A1_pdui_GEP= purrr::map2(COL1A1_GEP$GEPscore,COL1A1_GEP$Candidate_sur,function(x,y){
  x$new_id <-  gsub("\\.","\\-",str_sub(x$barcode,1,12))
  sub <- merge(x,y,by.x = "new_id",by.y = "barcode")
  sub <- dplyr::select(sub,new_id,GEP, barcode,event_id,PDUI,type,PDUIGroup,PDUIGroupB)
  return(sub)
}))
bind_cancer_GEP <- COL1A1_GEP$COL1A1_pdui_GEP %>% bind_rows()
pdf("/work/gywang/project/Immune_APA/Figure/COL1A1/COL1A1_PDUI GROUP GEP(no line).pdf",width = 10,height = 4)
ggplot(bind_cancer_GEP,aes(x=type,y=GEP,fill=PDUIGroup))+
  geom_boxplot(width=0.7)+
  stat_compare_means(aes(group=PDUIGroup),label="p.signif")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          #panel.grid.major=element_line(linetype = "longdash",color = "gray",size = 0.3),
          panel.grid.minor=element_blank(),
          axis.title.x=element_text(color="black",size=13),
          axis.text.x =element_text(color="black",size=12,angle =90,hjust = 1),
          axis.text.y = element_text(color="black",size=12,face = "bold"),
          #axis.title.y = element_blank(),
          #axis.line.y = element_line(colour = "black"),
          #axis.line.x = element_line(colour = "black"),
          panel.border = element_rect(color = "black",fill = NA))
dev.off()



all_three_feature <- merge(bind_cancer_CYT,bind_cancer_GEP)
all_three_feature <- merge(all_three_feature,bind_cancer_TIF)
all_three_feature$type <- as.character(all_three_feature$type)
all_three_feature %>% split(all_three_feature$type)-> tt
tt %>% lapply(function(x){
  CYT_diff <- mean(x[x$PDUIGroup %in% "Lengthen",]$CYT) - mean(x[x$PDUIGroup %in% "Shortening",]$CYT)
  CYT_p <- t.test(x[x$PDUIGroup %in% "Lengthen",]$CYT,x[x$PDUIGroup %in% "Shortening",]$CYT)
  CYT <- c(Feature="CYT",L_S_diff=CYT_diff,Pvalue=CYT_p$p.value)
  GEP_diff <- mean(x[x$PDUIGroup %in% "Lengthen",]$GEP) - mean(x[x$PDUIGroup %in% "Shortening",]$GEP)
  GEP_p <- t.test(x[x$PDUIGroup %in% "Lengthen",]$GEP,x[x$PDUIGroup %in% "Shortening",]$GEP)
  GEP <- c(Feature="GEP",L_S_diff=GEP_diff,Pvalue=GEP_p$p.value)
  TIF_diff <- mean(x[x$PDUIGroup %in% "Lengthen",]$TIF) - mean(x[x$PDUIGroup %in% "Shortening",]$TIF)
  TIF_p <- t.test(x[x$PDUIGroup %in% "Lengthen",]$TIF,x[x$PDUIGroup %in% "Shortening",]$TIF)
  TIF <- c(Feature="TIF",L_S_diff=TIF_diff,Pvalue=TIF_p$p.value)
  sub_diff <- rbind(CYT,GEP,TIF) %>% as.data.frame()
  sub_diff$L_S_diff <- as.numeric(sub_diff$L_S_diff)
  sub_diff$Pvalue <- as.numeric(sub_diff$Pvalue)
  sub_diff$cancer <- unique(x$type)
  return(sub_diff)
}) %>% bind_rows() -> all_feature_diff


all_feature_diff$Sig <- ifelse(all_feature_diff$Pvalue < 0.05,"Sig","NonSig")
pdf("/work/gywang/project/Immune_APA/Figure/COL1A1/Immune feature in diff COL1A1_PDUI group.pdf",width = 8,height = 5)
ggplot(all_feature_diff,aes(x=cancer,y=Feature ,shape=Sig ))+
  geom_point(aes(color=L_S_diff,size=-log10(Pvalue)))+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          #panel.grid.minor=element_blank(),
          axis.title.x=element_text(color="black",size=12),
          axis.text.x =element_text(color="black",size=10,angle = 90,hjust = 1,vjust = 0.5),
          axis.text.y = element_text(color="black",size=10,face = "bold"),
          #axis.line.y = element_line(colour = "black"),
          #axis.line.x = element_line(colour = "black"),
          panel.border = element_rect(color = "black",fill = NA))+
  scale_size_continuous(limit=c(0,10),range = c(1,6),breaks=c(-log10(0.05),3,5),labels=c("0.05","0.001","<1e-5"),name="P")+
  scale_shape_manual(values = c(1,19))+
  scale_color_gradient2(low="#303F9F",high="#b71c1c",mid = "#81d4fa")
dev.off()  
scale_color_gradient2(colours=c("#303F9F","#03A9F4","gray","red"))
  
  
dev.off()             


