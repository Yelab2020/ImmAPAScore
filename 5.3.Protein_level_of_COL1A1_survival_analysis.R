library(magrittr)
library(maxstat)
library(ggplot2)
library(survival)
options(stringsAsFactors = F)
my.wilcox.test <- function(...) {
  obj<-  try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
#MelanomaProteinExpession <- read.delim("/work/HEYi/Public/Protein/Melanoma_protein_2021_RawData.txt",sep = "\t")
MelanomaProteinExpession_60filter <- read.delim("/work/HEYi/Public/Protein/Melanoma_protein_2021_FilterData60%.txt",sep = "\t")
MelanomaProteinClinical <- read.delim("/work/HEYi/Public/Protein/Melanoma_protein_2021_ClinicalData.txt")
MelanomaProteinClinical <- dplyr::select(MelanomaProteinClinical,-5,-6,-10:-12,-14,-16)
colnames(MelanomaProteinClinical)[6:9] <- c("OST","OS","Treatment","Technical_set")



Candidate_GeneName <- "COL1A1"

COL1A1_Protein_Melanoma <- OneGene_Survival_MelanomaProtein_Data(Candidate_GeneName ="COL1A1")

pdf("/work/gywang/project/Immune_APA/Figure/COL1A1/protein_level_OF_COL1A1_survival_IN_melanoma.pdf",width = 6,height = 4)
output$Candidate_Plot_Best
COL1A1_Protein_Melanoma$Candidate_Plot_Best
dev.off()

pdf("/work/gywang/project/Immune_APA/Figure/COL1A1/protein_level_OF_COL1A1_and_response_IN_melanoma.pdf")
ggplot(COL1A1_Protein_Melanoma$ProExp_Clinical,aes(x=Treatment,y=ProExp,color=Response))+
  geom_boxplot(aes(x=Treatment,y=ProExp,color=Response))+
  stat_compare_means(aes(group=Response))
ggplot(COL1A1_Protein_Melanoma$ProExp_Clinical,aes(y=ProExp,x=Response))+
  geom_boxplot(aes(y=ProExp,color=Response))+
  stat_compare_means(aes(group=Response))
dev.off()











OneGene_Survival_MelanomaProtein_Data <- function(Candidate_GeneName){
  Candidate_GENEexp <- MelanomaProteinExpession_60filter %>% filter(Gene.names %in% Candidate_GeneName) %>% t() %>% as.data.frame()  %>% set_colnames("ProExp")
  Candidate_GENEexp$SampleName <- Candidate_GENEexp %>% rownames()
  MelanomaProteinClinical_Candidate <- merge(MelanomaProteinClinical,Candidate_GENEexp,by.x="Sample.name",by.y="SampleName") 
  MelanomaProteinClinical_Candidate$ProExp <- as.numeric(MelanomaProteinClinical_Candidate$ProExp)
  MelanomaProteinClinical_Candidate <- MelanomaProteinClinical_Candidate %>% dplyr::filter(ProExp>0 )
  MelanomaProteinClinical_Candidate$CandidateGene <- Candidate_GeneName
#  MelanomaProteinClinical_Candidate
  
  # divide all sample into diff cluster
 MelanomaProteinClinical_Candidate <- MelanomaProteinClinical_Candidate %>% dplyr::mutate(ProteinExpGroup = ifelse(ProExp >= median(ProExp),"High","Low"),
                                         OS = ifelse(OS %in% "N",0,1),
                                         OST = as.numeric(as.character(OST)))
  aa_try <- try(aa <- maxstat.test(survival::Surv(OST, OS) ~ ProExp, data=MelanomaProteinClinical_Candidate,smethod="LogRank"),silent=TRUE)
  if (is(aa_try , "try-error")) {NA} else {MelanomaProteinClinical_Candidate <-  MelanomaProteinClinical_Candidate %>% dplyr::mutate(PDUIGroupB = ifelse(ProExp > aa_try$estimate,"High","Low"))}
  
  ##########  
  # caculate HR
  .x <- MelanomaProteinClinical_Candidate
  model1 <- coxph(survival::Surv(OST, OS) ~ ProExp, data=.x,  na.action=na.exclude)
  HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
  Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
  HR_detail = summary(model1)
  CILow =  HR_detail$conf.int[,"lower .95"]
  CIHigh =  HR_detail$conf.int[,"upper .95"]
  model1 <- survival::survdiff(survival::Surv(OST, OS) ~ ProteinExpGroup,data=.x, na.action=na.exclude)
  KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(.x$PDUIGroup)))-1)
  model2 <- survival::survdiff(survival::Surv(OST, OS) ~ PDUIGroupB,data=.x, na.action=na.exclude)
  KMP2 <- 1-pchisq(model2$chisq, df=length(levels(factor(.x$PDUIGroupB)))-1)
  Candidate_HR <-  data.frame(APAevents=Candidate_GeneName,HR, Coxp,CILow,CIHigh,KMP,BestKMP=KMP2)
  
  ############
  #plot the survival curve median cutoff
  x1 <-MelanomaProteinClinical_Candidate
  model1 <- survival::coxph(survival::Surv(OST, OS) ~ ProExp, data=x1,  na.action=na.exclude)
  HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
  Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
  HR_detail = summary(model1)
  CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
  fit<- survival::survfit(survival::Surv(OST, OS) ~ ProteinExpGroup,data=x1)
  model1 <- survival::survdiff(survival::Surv(OST, OS) ~ ProteinExpGroup,data=x1, na.action=na.exclude)
  KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(x1$ProteinExpGroup)))-1)
  p <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                             pval = F,
                             title=paste("Protein level of ",x1$CandidateGene,"of Melanoma patients","P value = ",signif(as.numeric(KMP),digits = 2),
                                         "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                             fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                             legend.labs = c(paste("High (n=",table(x1$ProteinExpGroup)[[1]],")",sep=""),paste("Low (n=",table(x1$ProteinExpGroup)[[2]],")",sep=""))  )
  p$plot <- p$plot + labs(x="Months")+ theme(legend.title = element_blank(),legend.background = element_blank())
  Candidate_Plot_Median <- p$plot
  
  
  #############
  #plot the survival curve Best cutoff
  x2 <-MelanomaProteinClinical_Candidate
  model1 <- survival::coxph(survival::Surv(OST, OS) ~ ProExp, data=x2,  na.action=na.exclude)
  HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
  Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
  HR_detail = summary(model1)
  CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
  fit<- survival::survfit(survival::Surv(OST, OS) ~ PDUIGroupB,data=x2)
  model1 <- survival::survdiff(survival::Surv(OST, OS) ~ PDUIGroupB,data=x2, na.action=na.exclude)
  KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(x2$PDUIGroupB)))-1)
  p <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                             pval = F,
                             title=paste("Protein level of ",x1$CandidateGene,"of Melanoma patients","P value = ",signif(as.numeric(KMP),digits = 2),
                                         "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                             fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                             legend.labs = c(paste("High (n=",table(x2$PDUIGroupB)[[1]],")",sep=""),paste("Low (n=",table(x2$PDUIGroupB)[[2]],")",sep="")) )
  p$plot <- p$plot + labs(x="Months")+ theme(legend.title = element_blank(),legend.background = element_blank())
  Candidate_Plot_Best <- p$plot
  
  #############
  # return the result
  output <- list(MelanomaProteinClinical_Candidate,Candidate_HR,Candidate_Plot_Median,Candidate_Plot_Best)
  names(output) <- c("ProExp_Clinical","Candidate_HR","Candidate_Plot_Median","Candidate_Plot_Best")
  return(output)
}
