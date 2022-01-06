library(magrittr)
library(tidyverse)
library(maxstat)
library(ggplot2)
library(survival)
options(stringsAsFactors = F)
my.wilcox.test <- function(...) {
  obj<-  try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

TCGA_APA <- readRDS("/work/gywang/project/Immune_APA/Data/TCGA_Data/pancan32_APA.DATA.rds.gz")
TCGA_APA$APA <- lapply(TCGA_APA$APA,function(x){
  x <- x %>% dplyr::select(-2:-3) #%>% set_colnames(c("APA.events",str_sub(colnames(x)[-1:-3],1,12)))
})

ClinicalData <- read.delim("/data/share_data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
ClinicalData$PFI.time <- as.numeric(ClinicalData$PFI.time)
ClinicalData$OS.time <- as.numeric(ClinicalData$OS.time)
ClinicalData <- filter(ClinicalData,type!="READ")
ClinicalData$PFI <- as.numeric(ClinicalData$PFI)
ClinicalData$OS <- as.numeric(ClinicalData$OS)


OneAPA_OS_SurvialAcrossTCGA <- function(Candidate){
  TCGA_APA %>% dplyr::mutate(APA=purrr::map(.x=APA,function(.x){
    .x <- .x %>% dplyr::filter(event_id %in% Candidate) %>% tidyr::gather(key="barcode",value="PDUI",-event_id)
    .x$PDUI <- .x$PDUI*100
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
    .x <- .x %>% 
      dplyr::mutate(barcode = gsub("\\.","-",substr(barcode,1,12))) %>% 
      dplyr::inner_join(ClinicalData,by=c("barcode"="bcr_patient_barcode"))
    .x <- .x %>% dplyr::mutate(OS.time=as.numeric(OS.time)/30,PFI.time=as.numeric(PFI.time)/30) %>% dplyr::filter(!is.na(OS.time))
    .x <- .x %>% dplyr::filter(!is.na(PDUI)) %>%
      dplyr::mutate(PDUIGroup = ifelse(PDUI >= median(PDUI),"Lengthen","Shortening"),
                    OS = as.numeric(as.character(OS)),
                    PFI = as.numeric(as.character(PFI)))
    aa_try <- try(aa <- maxstat.test(survival::Surv(OS.time, OS) ~ PDUI, data= .x,smethod="LogRank"),silent=TRUE)
    if (is(aa_try , "try-error")) {return(NA)} else {
      .x <- .x %>% dplyr::mutate(PDUIGroupB = ifelse(PDUI > aa_try$estimate,"Lengthen","Shortening"))
      print(nrow(.x))
      return(.x)
    }}))
  
  ## HR
  .x <- CandidateAPA$Candidate_sur[[1]]
  CandidateAPA <- CandidateAPA[!is.na(CandidateAPA$Candidate_sur),] %>% dplyr::mutate(Candidate_HR=purrr::map(.x=Candidate_sur,function(.x){
    model1 <- coxph(survival::Surv(OS.time, OS) ~ PDUI, data=.x,  na.action=na.exclude)
    HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
    Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
    HR_detail = summary(model1)
    CIShortening =  HR_detail$conf.int[,"lower .95"]
    CILengthen =  HR_detail$conf.int[,"upper .95"]
    model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ PDUIGroup,data=.x, na.action=na.exclude)
    KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$PDUIGroup)))-1)
    model2 <- survival::survdiff(survival::Surv(OS.time, OS) ~ PDUIGroupB,data=.x, na.action=na.exclude)
    KMP2 <- 1-pchisq(model2$chisq, df=length(levels(factor(.x$PDUIGroupB)))-1)
    print(nrow(.x))
    data.frame(cancertypes=unique(.x$type),HR, Coxp,CIShortening,CILengthen,KMP,BestKMP=KMP2)
  }))
  
  ## surv_plot
  CandidateAPA <- CandidateAPA[!is.na(CandidateAPA$Candidate_sur),]  %>% dplyr::mutate(Candidate_plot=purrr::map(.x=CandidateAPA[!is.na(CandidateAPA$Candidate_sur),]$Candidate_sur,function(.x){
    model1 <- survival::coxph(survival::Surv(OS.time, OS) ~ PDUI, data=.x,  na.action=na.exclude)
    HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
    Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
    HR_detail = summary(model1)
    CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
    fit<- survival::survfit(survival::Surv(OS.time, OS) ~ PDUIGroup,data=.x)
    model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ PDUIGroup,data=.x, na.action=na.exclude)
    KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$PDUIGroup)))-1)
    p <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                               # risk.table = TRUE, risk.table.col = "strata",
                               pval = F, title=paste(unique(.x$type),": p = ",signif(as.numeric(KMP),digits = 2),
                                                     "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                               fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                               legend.labs = c(paste("Lengthen (n=",table(.x$PDUIGroup)[[1]],")",sep=""),
                                               paste("Shortening (n=",table(.x$PDUIGroup)[[2]],")",sep="")))
    p$plot <- p$plot + labs(x="Months")+ theme(legend.title = element_blank(),legend.background = element_blank())
    print(nrow(.x))
    return(p$plot)
  }))  
  
  ####Best pvalue
  CandidateAPA <- CandidateAPA %>% dplyr::mutate(Candidate_plotBestCurv=purrr::map(.x=Candidate_sur,function(.x){
    model1 <- survival::coxph(survival::Surv(OS.time, OS) ~ PDUI, data=.x,  na.action=na.exclude)
    HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
    Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
    HR_detail = summary(model1)
    CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
    fit<- survival::survfit(survival::Surv(OS.time, OS) ~ PDUIGroupB,data=.x)
    model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ PDUIGroupB,data=.x, na.action=na.exclude)
    KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$PDUIGroupB)))-1)
    p <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                               # risk.table = TRUE, risk.table.col = "strata",
                               pval = F, title=paste(unique(.x$type),": p = ",signif(as.numeric(KMP),digits = 2),
                                                     "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                               fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                               legend.labs = c(paste("Lengthen (n=",table(.x$PDUIGroupB)[[1]],")",sep=""),
                                               paste("Shortening (n=",table(.x$PDUIGroupB)[[2]],")",sep="")))
    p$plot <- p$plot + labs(x="Months")+theme(legend.title = element_blank(),legend.background = element_blank())
    return(p$plot)
  }))
  
  # plot HR
  Candidate_HR <- dplyr::bind_rows(CandidateAPA$Candidate_HR)
  Candidate_HRSig <- Candidate_HR
  Candidate_HRSig <- Candidate_HRSig[order(Candidate_HRSig$HR),]
  Candidate_HRSig$Group <- ifelse(Candidate_HRSig$HR >1,"High","Low")
  Candidate_HRSig$Psig <- ifelse(Candidate_HRSig$Coxp < 0.05,"Sig","nonSig")
p_HR <- ggplot(Candidate_HRSig,aes(y=cancertypes,x= HR ,size=-log10(Coxp+10^-4)))+
    geom_segment(data= Candidate_HRSig,aes(x=CIShortening ,y=cancertypes,xend=CILengthen,yend=cancertypes),color="black",size=1,guide=F)+
    geom_point(aes(color=Psig ))+
    scale_color_manual(limits=c("Sig","nonSig"),values=c("red","gray"))+
    scale_size_continuous(limit=c(0,4),range = c(1,6),breaks=c(-log10(0.05),4),labels=c("0.05","1e-4"))+
    scale_x_log10()+
    scale_y_discrete(limit=Candidate_HRSig$cancertypes)+
    theme(panel.background=element_blank(),
          panel.grid=element_blank(),
          #axis.title=element_blank(),
          axis.title.y = element_blank(),
          axis.text=element_text(size=12,color = "black"),
          axis.ticks.y=element_blank(),
          axis.line.x = element_line(colour = "black"),
          legend.text=element_text(size=10),
          legend.title=element_text(size=12),
          legend.key = element_rect(fill="white",colour = "black"))+
    geom_vline(xintercept = 1,size=0.5)+
    geom_point(data=Candidate_HRSig[Candidate_HRSig$Coxp < 0.05,],aes(y=cancertypes,x= HR ,size=-log10(Coxp+10^-4)),shape=1,color="black")+
    labs(x=paste("HR of",Candidate,"APA PDUI"))
#p_SigBestCut <- cowplot::plot_grid(plotlist = CandidateAPA[Candidate_HR$BestKMP< 0.05,]$Candidate_plotBestCurv, rel_heights = c(3,3),rel_widths = c(3,3),ncol = 6)
p_AllBestCut <- cowplot::plot_grid(plotlist = CandidateAPA$Candidate_plotBestCurv,
                                   rel_heights = c(3,3),rel_widths = c(3,3),ncol = 6)
All_plot <- list(HR=p_HR,AllBestCut=p_AllBestCut)
result <- list(CandidateAPA,All_plot)
return(result)
}

