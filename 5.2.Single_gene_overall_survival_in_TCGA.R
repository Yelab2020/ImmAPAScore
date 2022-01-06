library(magrittr)
library(maxstat)
library(ggplot2)
library(survival)
options(stringsAsFactors = F)
my.wilcox.test <- function(...) {
  obj<-  try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

ExprFPKM <- readr::read_rds("/data/share_data/TCGA_preliminary/TCGA_33CanFPKM.rds.gz")
ExprFPKM <- ExprFPKM[ExprFPKM$cancertypes !="LAML",]
ClinicalData <- read.delim("/data/share_data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
ClinicalData$PFI.time <- as.numeric(ClinicalData$PFI.time)
ClinicalData$OS.time <- as.numeric(ClinicalData$OS.time)

if(.y=="SKCM"){
  .x <- .x %>% dplyr::filter(substr(barcode,14,15) %in% "06") %>% dplyr::mutate(barcode=substr(gsub("\\.","-",barcode),1,12))
}else{
  .x <- .x %>% dplyr::filter(substr(barcode,14,15) %in% c("01","03")) %>% dplyr::mutate(barcode=substr(gsub("\\.","-",barcode),1,12))
}
return(.x)


###################### return all plot and survival object
OneGeneSurvialAcrossTCGA <- function(Candidate){
  BestKMP_T=T
  KMP_T=F
  Coxp_T=F
  HR_T=F
  ExprFPKM %>% dplyr::mutate(Expr=purrr::map(.x=Expr,function(.x){
    .x %>% dplyr::filter(Gene %in% Candidate) %>% tidyr::gather(key="barcode",value="Expression",-Gene) %>%
      dplyr::mutate(Expression=log2(Expression+1))})) %>%
    dplyr::mutate(Expr = purrr::map(.x=Expr,function(.x){.x %>% 
        dplyr::mutate(types = data.frame(do.call(rbind,strsplit(barcode,"-")))$X1,
                      sample_class = data.frame(do.call(rbind,strsplit(barcode,"-")))$X2)})) -> Candidateexp 
  ##Candidate_sur
  Candidateexp <- Candidateexp %>% dplyr::mutate(Candidate_sur=purrr::map(.x=Expr,function(.x){
    .x <- .x %>%  dplyr::mutate(barcode = gsub("\\.","-",substr(barcode,1,12))) %>% 
      dplyr::inner_join(ClinicalData,by=c("barcode"="bcr_patient_barcode")) ##按barcode添加
    .x <- .x %>% dplyr::mutate(OS.time=as.numeric(OS.time)/30,PFI.time=as.numeric(PFI.time)/30) %>% dplyr::filter(!is.na(OS.time))
    .x <- .x %>%    dplyr::mutate(ExpressionGroup = ifelse(Expression > median(Expression),"High","Low"),##与中位数相比
                    OS = as.numeric(as.character(OS)),
                    PFI = as.numeric(as.character(PFI)))
    aa_try <- try(aa <- maxstat.test(survival::Surv(OS.time, OS) ~ Expression, data= .x,smethod="LogRank"),silent=TRUE)
    .x <- .x %>% dplyr::mutate(ExpressionGroupB = ifelse(Expression >aa_try$estimate,"High","Low"))
    return(.x)
  }))
  
  ##Candidate_HR
  Candidateexp <- Candidateexp %>% dplyr::mutate(Candidate_HR=purrr::map(.x=Candidate_sur,function(.x){
    model1 <- coxph(survival::Surv(OS.time, OS) ~ Expression, data=.x,  na.action=na.exclude)
    HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
    Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
    HR_detail = summary(model1)
    CILow =  HR_detail$conf.int[,"lower .95"]
    CIHigh =  HR_detail$conf.int[,"upper .95"]
    model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ ExpressionGroup,data=.x, na.action=na.exclude)
    KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$ExpressionGroup)))-1)
    model2 <- survival::survdiff(survival::Surv(OS.time, OS) ~ ExpressionGroupB,data=.x, na.action=na.exclude)
    KMP2 <- 1-pchisq(model2$chisq, df=length(levels(factor(.x$ExpressionGroupB)))-1)
    print(nrow(.x))
    data.frame(cancertypes=unique(.x$type),HR, Coxp,CILow,CIHigh,KMP,BestKMP=KMP2)
  }))
  
  ##Candidate_plot
  .x <- Candidateexp$Candidate_sur[[1]]
  Candidateexp <- Candidateexp %>% dplyr::mutate(Candidate_plot=purrr::map(.x=Candidateexp$Candidate_sur,function(.x){
    model1 <- survival::coxph(survival::Surv(OS.time, OS) ~ Expression, data=.x,  na.action=na.exclude)
    HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
    Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
    HR_detail = summary(model1)
    CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
    fit <- survival::survfit(survival::Surv(OS.time, OS) ~ ExpressionGroup,data=.x)
    model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ ExpressionGroup,data=.x, na.action=na.exclude)
    KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$ExpressionGroup)))-1)
    p <-  survminer::ggsurvplot(fit,palette = c("red", "blue"),
                                  pval = F, title=paste(unique(.x$type),": p = ",signif(as.numeric(KMP),digits = 2), "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                                  fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                                  legend.labs = c(paste("High (n=",table(.x$ExpressionGroup)[[1]],")",sep=""),
                                      paste("Low (n=",table(.x$ExpressionGroup)[[2]],")",sep="")))
    p$plot <- p$plot + labs(x="Months")+ theme(legend.title = element_blank(),legend.background = element_blank())
    return(p$plot)
  }))
  
  ####Best pvalue
  Candidateexp <- Candidateexp %>% dplyr::mutate(Candidate_plotBestCurv=purrr::map(.x=Candidate_sur,function(.x){
    model1 <- survival::coxph(survival::Surv(OS.time, OS) ~ Expression, data=.x,  na.action=na.exclude)
    HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
    Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
    HR_detail = summary(model1)
    CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
    fit<- survival::survfit(survival::Surv(OS.time, OS) ~ ExpressionGroupB,data=.x)
    model1 <- survival::survdiff(survival::Surv(OS.time, OS) ~ ExpressionGroupB,data=.x, na.action=na.exclude)
    KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$ExpressionGroupB)))-1)
    p <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                               # risk.table = TRUE, risk.table.col = "strata",
                               pval = F, title=paste(unique(.x$type),": p = ",signif(as.numeric(KMP),digits = 2),
                                                     "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                               fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                               legend.labs = c(paste("High (n=",table(.x$ExpressionGroupB)[[1]],")",sep=""),
                                               paste("Low (n=",table(.x$ExpressionGroupB)[[2]],")",sep="")))
    p$plot <- p$plot + labs(x="Months")+theme(legend.title = element_blank(),legend.background = element_blank())
    return(p$plot)
  }))
  
  Candidate_HR <- dplyr::bind_rows(Candidateexp$Candidate_HR)
  if(HR_T==T){
    Candidate_HRSig <- Candidate_HR[Candidate_HR$Coxp < 0.05,]
    Candidate_HRSig <- Candidate_HRSig[order(Candidate_HRSig$HR),]
    Candidate_HRSig$Group <- ifelse(Candidate_HRSig$HR >1,"High","Low")
  }else{
    Candidate_HRSig <- Candidate_HR
    Candidate_HRSig <- Candidate_HRSig[order(Candidate_HRSig$HR),]
    Candidate_HRSig$Group <- ifelse(Candidate_HRSig$HR >1,"High","Low")
  }
  Candidate_HRSig$Psig <- ifelse(Candidate_HRSig$Coxp < 0.05,"Sig","nonSig")
  
  ####### plot
  p_HR <- ggplot(Candidate_HRSig,aes(y=cancertypes,x= HR ,size=-log10(Coxp+10^-4)))+
    geom_segment(data= Candidate_HRSig,aes(x=CILow ,y=cancertypes,xend=CIHigh,yend=cancertypes),color="black",size=1,guide=F)+
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
    labs(x=paste("HR of",Candidate,"gene expression"))
  p_SigBestCut <- cowplot::plot_grid(plotlist = Candidateexp[Candidate_HR$BestKMP< 0.05,]$Candidate_plotBestCurv,
                           rel_heights = c(3,3),rel_widths = c(3,3),ncol = 6)
  p_AllBestCut <- cowplot::plot_grid(plotlist = Candidateexp$Candidate_plotBestCurv,
                           rel_heights = c(3,3),rel_widths = c(3,3),ncol = 6)
  
  All_plot <- list(HR=p_HR,SigBestCut=p_SigBestCut,AllBestCut=p_AllBestCut)
  
  result <- list(Candidateexp,All_plot)
  return(result)
  }