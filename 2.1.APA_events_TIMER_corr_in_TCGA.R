library(magrittr)
library(ggplot2)
options(stringsAsFactors = F)
library(doParallel)
library(doMC)
registerDoMC()
library(foreach)
library(ppcor)
my.cor.test <- function(...) {
  obj<-try(cor.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}
my.pcor.test <- function(...) {
  obj<-try(pcor.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}
Purity <- readr::read_rds("/home/gywang/project/20200831ye/PanCan_tumorPurity.rds.gz")
PurityAll <- dplyr::bind_rows(Purity$Purity) %>% dplyr::filter(substr(barcode,14,15) %in% c("01","03")) %>%
  dplyr::mutate(barcode=gsub("\\-","\\.",substr(barcode,1,15)))

GeneAnn <- readr::read_rds("/home/gywang/project/20200831ye/GeneAnnotation_GeneTypes.rds")

TCGA_APA <- readr::read_rds("/home/gywang/project/20200831ye/pancan32_APA.DATA.rds.gz")
TCGA_Timer <- readr::read_rds("/home/gywang/project/APA_CancerImmune//immune_pathway_data/TCGA_TIMER.rds.gz")
TCGA_APA_Timer <- dplyr::inner_join(TCGA_APA,TCGA_Timer,by=c("cancer_types"))
TCGA_APA_Timer <- TCGA_APA_Timer %>% dplyr::mutate(Corr=purrr::map2(.x=AI,.y=APA,function(.x,.y){
  .x <-  .x %>% data.frame
  .x$barcode <- gsub("-",".",rownames(.x))
  .x <- .x[substr(.x$barcode,14,15) %in%c("01","03"),]
  .x$barcode <- substr(.x$barcode,1,15)
  .y <- .y[,c("event_id","loci", "Proximal_Site",colnames(.y)[substr(colnames(.y),14,15) %in% c("01","03")])]
  colnames(.y) <- c("event_id","loci", "Proximal_Site",substr(colnames(.y)[4:ncol(.y)],1,15))
  CommonSamples <- intersect(.x$barcode,grep("TCGA",colnames(.y),value=T))
  CommonSamples <- intersect(CommonSamples,PurityAll$barcode)
  TumorPurity <- PurityAll[match(CommonSamples,PurityAll$barcode),]$Purity
  
  .x <- .x[!duplicated(.x$barcode),]
  .x <- .x %>% tidyr::gather(key=ImmuneCell,value=Infiltration,-barcode) %>%
    tidyr::spread(barcode,Infiltration)
  .x <- .x[,c("ImmuneCell",CommonSamples)]
  .y <- .y[,c("event_id","loci", "Proximal_Site",CommonSamples)]
  lapply(split(.x,.x$ImmuneCell), function(y){
    as.numeric(y[,CommonSamples]) -> dy
    Corr = data.frame(APAevents = .y$event_id,ImmuneCell =rep(unique(y$ImmuneCell),times=length(.y$event_id)))
    Corr[,c("estimate.rho","p.value")] <- t(apply(.y[,CommonSamples],1,function(x){
      xval <- as.numeric(x)
      unlist(my.pcor.test( xval[!is.na(xval)],dy[!is.na(xval)], TumorPurity[!is.na(xval)]))
    }))
    Corr$FDR <- signif(p.adjust(Corr$p.value,method="fdr"),digits = 2)
    return(Corr)
  }) %>% dplyr::bind_rows() -> AllCorr
  print(nrow(.x))
  return(AllCorr)
}))
readr::write_rds(TCGA_APA_Timer,path="/home/yye/Project/2020APA_Imm/ImmAPA/ImmInfiltration/TCGA_APA_Timer_Ppearson.rds.gz")


###
#TCGA_APA_Timer  <- readr::read_rds("/home/yye/Project/2020APA_Imm/ImmAPA/ImmInfiltration/TCGA_APA_Timer.rds.gz")
TCGA_APA_Timer <- TCGA_APA_Timer %>% dplyr::mutate(CorrSig=purrr::map(.x=Corr,function(.x){
 .x <-  .x[abs(.x$estimate.rho) > 0.2 & .x$p.value < 0.05,]
  .x[!is.na(.x$APAevents),]
}))
RankScore <- readr::read_rds("/home/yye/Project/2020APA_Imm/ImmAPA/ImmAPAscore/APA_ImmScore.rds.gz")
RankScore$cancer_types <- gsub("_ImmAPA","",RankScore$cancer_types)

TCGA_APA_Timer_RankScore <- dplyr::inner_join(RankScore,TCGA_APA_Timer,by="cancer_types")
TCGA_APA_Timer_RankScore <- TCGA_APA_Timer_RankScore %>% dplyr::mutate(Overlap = purrr::map2(.x=CorrSig,.y=RankScoreSig,function(.x,.y){
  #.y <- .y[.y$Class %in% "Pos",]
  #.x <- .x[.x$estimate.rho>0,]
  PathwaysAPA <- names(table(.y$APAevents)[table(.y$APAevents) >= 1])
  aa <- lapply(split(.x,.x$ImmuneCell),function(x){
  data.frame(ImmuneCells=length(setdiff(x$APAevents,PathwaysAPA)),ImmPathways=length(setdiff(PathwaysAPA,x$APAevents)),Overlap=length(intersect(x$APAevents,PathwaysAPA)))
   }) %>% dplyr::bind_rows(.)
  aa$OverPer <- aa$Overlap/length(PathwaysAPA) *100
  return(aa)
}))


TCGA_APA_CIBERSORT  <- readr::read_rds("/home/yye/Project/2020APA_Imm/ImmAPA/ImmInfiltration/TCGA_APA_CIBERSORT.rds.gz")
TCGA_APA_CIBERSORT <- TCGA_APA_CIBERSORT %>% dplyr::mutate(CorrSig=purrr::map(.x=Corr,function(.x){
  .x <-  .x[abs(.x$estimate.rho) > 0.2 & .x$FDR < 0.05,]
  .x[!is.na(.x$APAevents),]
}))
TCGA_APA_CIBERSORT_RankScore <- dplyr::inner_join(RankScore,TCGA_APA_CIBERSORT,by="cancer_types")
TCGA_APA_CIBERSORT_RankScore <- TCGA_APA_CIBERSORT_RankScore %>% dplyr::mutate(Overlap = purrr::map2(.x=CorrSig,.y=RankScoreSig,function(.x,.y){
  #.y <- .y[.y$Class %in% "Pos",]
  #.x <- .x[.x$estimate.rho>0,]
  PathwaysAPA <- names(table(.y$APAevents)[table(.y$APAevents) >= 1])
  aa <- lapply(split(.x,.x$ImmuneCell),function(x){
    data.frame(ImmuneCells=length(setdiff(x$APAevents,PathwaysAPA)),ImmPathways=length(setdiff(PathwaysAPA,x$APAevents)),Overlap=length(intersect(x$APAevents,PathwaysAPA)))
  }) %>% dplyr::bind_rows(.)
  aa$OverPer <- aa$Overlap/length(PathwaysAPA) *100
  return(aa)
}))
TCGA_APA_CIBERSORT_RankScore <- TCGA_APA_CIBERSORT_RankScore %>% dplyr::mutate(OverlapAll = purrr::map2(.x=CorrSig,.y=RankScoreSig,function(.x,.y){
  #.y <- .y[.y$Class %in% "Pos",]
  #.x <- .x[.x$estimate.rho>0,]
  PathwaysAPA <- names(table(.y$APAevents)[table(.y$APAevents) >= 1])
  data.frame(ImmuneCells=length(setdiff(unique(.x$APAevents),PathwaysAPA)),ImmPathways=length(setdiff(PathwaysAPA,unique(.x$APAevents))),
             Overlap=length(intersect(unique(.x$APAevents),PathwaysAPA)),OverPer=length(intersect(unique(.x$APAevents),PathwaysAPA))/length(PathwaysAPA)*100)
  
}))

