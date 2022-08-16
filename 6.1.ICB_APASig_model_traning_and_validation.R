library(glmnet)
library(tidyverse)
library(magrittr)
library(maxstat)
library(survival) 
library(survminer)
library(timeROC)
library(MASS)
options(stringsAsFactors = F)
theme1 <-     theme(panel.background=element_rect(colour=NA,fill="white"),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    #axis.title.x = element_blank(),
                    axis.text.x = element_text(color="black",size=10,angle = 0,hjust = 0.5),
                    axis.text.y = element_text(color="black",size=10,face = "bold"),
                    axis.line.y = element_line(colour = "black"),
                    axis.line.x = element_line(colour = "black"),
                    panel.border = element_rect(color = "black",fill = NA))

my.wilcox.test <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.t.test <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)}
my.cor.test_spearman <- function(...) {
  obj<- try(cor.test(...,method="spearman"), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])
}


# 541 ImmAPA events
ImmAPA541 <- readr::read_rds("/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_ranking_n.rds.gz")
ImmAPA541_APAevents <- rownames(ImmAPA541)
ImmAPA541_GeneName <- data.frame(do.call(rbind,str_split(ImmAPA541_APAevents,"\\|")))[,2]
TCGA_APA <- readRDS("/work/gywang/project/Immune_APA/Data/TCGA_Data/pancan32_APA.DATA.rds.gz")
TCGA_APA$APA <- lapply(TCGA_APA$APA,function(x){
  x <- x %>% dplyr::select(-2:-3) # %>% set_colnames(c("APA.events",str_sub(colnames(x)[-1:-3],1,12)))
})

# TCGA meta infomation
ClinicalData <- read.delim("/data/share_data/TCGA_preliminary/TCGA_ClinicalData_20180420.txt")
ClinicalData$PFI.time <- as.numeric(ClinicalData$PFI.time)
ClinicalData$OS.time <- as.numeric(ClinicalData$OS.time)
ClinicalData <- filter(ClinicalData,type!="READ")
ClinicalData$PFI <- as.numeric(ClinicalData$PFI)
ClinicalData$OS <- as.numeric(ClinicalData$OS)
TCGA_SKCM_APA <- TCGA_APA[TCGA_APA$cancer_types %in% "SKCM",]$APA[[1]]
TCGA_SKCM_MetaInfo <- ClinicalData %>% dplyr::filter(type %in% "SKCM")
TCGA_SKCM_barcode <- data.frame(Run=colnames(TCGA_SKCM_APA)[-1:-3],
                                tumor_state=ifelse(str_sub(colnames(TCGA_SKCM_APA)[c(-1:-3)],14,15) %in% "01","primary","metastasis"),
                                MetaInfo_ID= substr(gsub("\\.","-",str_sub(colnames(TCGA_SKCM_APA)[-1:-3])),1,12))
TCGA_SKCM_MetaInfo <- dplyr::select(TCGA_SKCM_MetaInfo,1:6,14,OS,OS.time,DSS,DSS.time,PFI,PFI.time)
TCGA_SKCM_MetaInfo <- merge(TCGA_SKCM_MetaInfo,TCGA_SKCM_barcode,by.y = "MetaInfo_ID",by.x = "bcr_patient_barcode")  

# TCGA ICB predicted information
TCGA_SKCM_Response_Prediction <- read.delim("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/TCGA_SKCM_Response_Prediction.txt")
TCGA_SKCM_Response_Prediction$MetaInfo_ID <- gsub("\\.","-",TCGA_SKCM_Response_Prediction$X)
TCGA_SKCM_MetaInfo <- merge(TCGA_SKCM_MetaInfo,TCGA_SKCM_Response_Prediction,by.y = "MetaInfo_ID",by.x = "bcr_patient_barcode")  
colnames(TCGA_SKCM_MetaInfo)[c(9,12,13)] <- c("OST","PFS","PFST")


TCGA_SKCM_Metastasis_MetaInfo <- TCGA_SKCM_MetaInfo %>% dplyr::filter(tumor_state %in% "metastasis")
colnames(TCGA_SKCM_APA)[1] <- "APAevent"
TCGA_SKCM_MetaInfo$Tumor_stage <- ifelse(TCGA_SKCM_MetaInfo$ajcc_pathologic_tumor_stage %in% 
                                           c("I/II NOS","Stage I"    ,    "Stage IA" ,"Stage IB", "Stage II","Stage IIA",  "Stage IIB" ,"Stage IIC"),"StageI/II",
                                         ifelse(TCGA_SKCM_MetaInfo$ajcc_pathologic_tumor_stage %in% c("Stage III","Stage IIIA","Stage IIIB","Stage IIIC","Stage IV"),"StageIII/IV",NA))

# TCGA metastasis ImmAPA
TCGA_SKCM_APA_metastasis <- dplyr::select(TCGA_SKCM_APA,APAevent,TCGA_SKCM_MetaInfo[TCGA_SKCM_MetaInfo$tumor_state %in% "metastasis",]$Run)
TCGA_SKCM_metastasis_541ImmAPA <- TCGA_SKCM_APA_metastasis %>% filter(APAevent %in% ImmAPA541_APAevents)
TCGA_SKCM_metastasis_541ImmAPA <- TCGA_SKCM_metastasis_541ImmAPA[,-1] %>% apply(1,function(x){x <- as.numeric(x)*100;return(x)}) %>%
  set_colnames(TCGA_SKCM_metastasis_541ImmAPA$APAevent) %>% set_rownames(colnames(TCGA_SKCM_metastasis_541ImmAPA)[-1]) %>% as.data.frame()
TCGA_SKCM_metastasis_541ImmAPA_sur <- merge(TCGA_SKCM_metastasis_541ImmAPA,TCGA_SKCM_MetaInfo,by.x="row.names",by.y="Run")
TCGA_SKCM_metastasis_541ImmAPA_sur <- set_rownames(TCGA_SKCM_metastasis_541ImmAPA_sur,TCGA_SKCM_metastasis_541ImmAPA_sur$Run)
#readr::write_rds(TCGA_SKCM_metastasis_541ImmAPA_sur,"/work/gywang/project/Immune_APA/ImmAPA_ICB/data/TCGA_SKCM_metastasis_541ImmAPA_and_survival.rds")

BenefitID <- TCGA_SKCM_Metastasis_MetaInfo[TCGA_SKCM_Metastasis_MetaInfo$Responder %in% "True",]$Run
NonBenefitID <- TCGA_SKCM_Metastasis_MetaInfo[TCGA_SKCM_Metastasis_MetaInfo$Responder %in% "False",]$Run

# find different expressed ImmAPA
DiffAll <- apply(TCGA_SKCM_APA_metastasis[,TCGA_SKCM_Metastasis_MetaInfo$Run],1,function(x){                 # apply 
  c(NonBenefit_Mean = mean(as.numeric(x[names(x) %in% NonBenefitID]) ,na.rm = T),
    Benefit_Mean = mean(as.numeric(x[names(x) %in% BenefitID]) ,na.rm = T),
    Diff = mean(as.numeric(x[names(x) %in% BenefitID]),na.rm = T)- mean(as.numeric(x[names(x) %in% NonBenefitID]) ,na.rm = T),
    Pval = unlist(my.wilcox.test(as.numeric(x[names(x) %in% BenefitID]),as.numeric(x[names(x) %in% NonBenefitID]),na.rm = T)))
})  %>% t %>% data.frame %>% dplyr::mutate(APAevent=TCGA_SKCM_APA_metastasis$APAevent)
DiffAll$sig <- ifelse(DiffAll$Pval < 0.05 ,"Sig","Nonsig")

Benefit_DiffAPAevents <- DiffAll[DiffAll$sig %in% "Sig",]$APAevent
intersect(Benefit_DiffAPAevents,ImmAPA541_APAevents) -> Benifit_Diff_ImmAPA

# find ImmAPA that effect the survival
TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore <- lapply(Benifit_Diff_ImmAPA,function(x){
  sub <- TCGA_SKCM_APA_metastasis[TCGA_SKCM_APA_metastasis$APAevent %in% x,] %>% 
    tidyr::gather(key=Run,value=PDUI,-APAevent) %>% 
    dplyr::inner_join(TCGA_SKCM_MetaInfo,by="Run") %>%
    dplyr::mutate(PDUI=as.numeric(PDUI)*100) %>%
    dplyr::filter(!is.na(OST)) # sub 选出每个APA.events 对应的Pre的样本的表达值，以及各样本的生存信息
  sub <- sub %>% dplyr::mutate(PDUIGroup = ifelse(PDUI > median(PDUI),"High","Low")) ##  Median cutpoint
  sub$OS <- as.numeric(sub$OS)
  ## Cutpoint for best statistically significant p value
  aa_try <- try(aa <- maxstat.test(survival::Surv(OST, OS) ~ PDUI, data= sub,smethod="LogRank"),silent=TRUE)
  if (is(aa_try, "try-error")) {return(NA)}
  else{
    sub <- sub %>% dplyr::mutate(PDUIGroupB = ifelse(PDUI >= aa_try$estimate,"High","Low"))
    if(length(table(sub$PDUIGroup)) ==2 & length(table(sub$PDUIGroupB)) ==2){
      model1 <- survival::coxph(survival::Surv(OST,OS) ~ PDUI, data=sub,  na.action=na.exclude)
      HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
      Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
      HR_detail = summary(model1)
      CILow =  HR_detail$conf.int[,"lower .95"]
      CIHigh =  HR_detail$conf.int[,"upper .95"]
      CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") # Result[i1, c('CI_95%_for_HR')]
      model1 <- survival::survdiff(survival::Surv(OST, OS) ~ PDUIGroup,data=sub, na.action=na.exclude)
      KMP <- signif((1-pchisq(model1$chisq, df=length(levels(factor(sub$PDUIGroup)))-1)),digits = 3)
      model1 <- survival::survdiff(survival::Surv(OST, OS) ~ PDUIGroupB,data=sub, na.action=na.exclude)
      KMP_max <- signif((1-pchisq(model1$chisq, df=length(levels(factor(sub$PDUIGroupB)))-1)),digits = 3)
      return(data.frame(HR, CI,Coxp,KMP,KMP_max,maxCut=aa_try$estimate,maxHighNum=table(sub$PDUIGroupB)[1],maxLowNum=table(sub$PDUIGroupB)[2]))
    }else{return(NA)}
  }})
names(TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore) <- Benifit_Diff_ImmAPA

TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore <- TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore[!is.na(TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore)]
TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore <- dplyr::bind_rows(TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore ) %>% dplyr::mutate(APAevent=rep(names(TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore )))
TCGA_SKCM_APA_metastasis_ImmAPA_SurvScoreSig <- TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore[TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore$KMP_max < 0.1,]
TCGA_SKCM_APA_metastasis_ImmAPA_SurvScoreCoxpKMP <- TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore[TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore$KMP_max < 0.05 & TCGA_SKCM_APA_metastasis_ImmAPA_SurvScore$Coxp < 0.05,] 

# get the effective APA
sig1_APA_exp <- TCGA_SKCM_APA_metastasis[TCGA_SKCM_APA_metastasis$APAevent %in% TCGA_SKCM_APA_metastasis_ImmAPA_SurvScoreCoxpKMP$APAevent,]  

rownames(sig1_APA_exp) <- sig1_APA_exp$APAevent
sig1_APA_exp <- sig1_APA_exp[,-1] %>% apply(1,function(x){x <- as.numeric(x)*100;return(x)})  %>% set_rownames(colnames(sig1_APA_exp)[-1])  %>% as.data.frame()
sig1_APA_exp_sur <- merge(sig1_APA_exp,TCGA_SKCM_Metastasis_MetaInfo,by.x ="row.names",by.y ="Run" )
sig1_APA_exp_sur <- dplyr::filter(sig1_APA_exp_sur, !is.na(OST))
sig1_APA_exp <- sig1_APA_exp_sur %>% dplyr::select(colnames(sig1_APA_exp))

x <- as.matrix(sig1_APA_exp)
x <- ifelse(is.na(x),0,x)
y <- as.matrix(survival::Surv( sig1_APA_exp_sur$OST,sig1_APA_exp_sur$OS ))
# Lasso 回归建模 
lasso <- glmnet( x , y , family = "cox" , alpha = 1)  # 最后一行可以看出非0系数，17
print(lasso)
# 它表明，当λ值减小时，压缩参数随之减小，而系数绝对值随之增大。
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/lasso model sselect.pdf",width = 5,height = 4)
p2 <- plot(lasso, xvar = "lambda", label = TRUE)
dev.off()

#  4.2.2 挑选合适的λ值
set.seed(123)
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/lasso model select best lambda.pdf",width = 5,height = 4)
plot(fitCV)
dev.off()

fitCV <- cv.glmnet(x, y, family = "cox",
                   type.measure = "deviance",
                   nfolds = 20)

# === 输入这两个系数用来建模
fitCV$lambda.1se
coef( fitCV , s = "lambda.1se") #%>% table()
fitCV$lambda.min
coef( fitCV , s = "lambda.min")

coefficient  <-  coef(fitCV, s= fitCV$lambda.min) # 这个也可以显示哪些基因可以建模 即用来后续分析
Active.Index <-  which(as.numeric(coefficient) != 0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox  <- colnames(sig1_APA_exp)[Active.Index]
length(sig_gene_multi_cox) # 3个这里的基因输入后续多因素cox

sig_gene_multi_cox_geneName  <-  gsub("N.*.[0-9]\\|" , "" , gsub("\\|chr.*|*","",sig_gene_multi_cox))
sig1_APA_exp_sur_geneName    <-  sig1_APA_exp_sur
colnames(sig1_APA_exp_sur_geneName) <- gsub("N.*.[0-9]\\|","",gsub("\\|chr.*|*","",colnames(sig1_APA_exp_sur)))
colnames(sig1_APA_exp_sur_geneName)[c(7:9,11)] <- c("Age" , "Gender", "Race", "TumorStatus")
sig1_APA_exp_sur_geneName$Age <- as.numeric(sig1_APA_exp_sur_geneName$Age)
sig1_APA_exp_sur_geneName$Benifit <- ifelse(sig1_APA_exp_sur_geneName$Responder %in% "True","Benifit","NonBenifit")
sig1_APA_exp_sur_geneName
# =====================  4.3 利用筛选出来的少数 candidates 建立多因素 cox 回归模型
Mulcox <- coxph(Surv(OST,OS) ~ GPNMB+IFITM1+COL1A1+CSK, data =  sig1_APA_exp_sur_geneName )
Mulcox_coef <- coef(Mulcox)
Mulcox_Coxp <- signif(as.numeric(summary(Mulcox)$coefficients[,c( "Pr(>|z|)")]),digit=4)
Mulcox_coef_used = Mulcox_coef[Mulcox_Coxp < 0.05 ] # 剩下3个
Mulcox_coef_used

pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/Hazard ratio of GPNMB IFITM1 COL1A1 CSK APA event from 541 ImmAPA .pdf",width = 10,height = 8,onefile = FALSE)
ggforest( model = Mulcox                    ,
          data  = sig1_APA_exp_sur_geneName ,
          main  = "Hazard ratio"            ,
          cpositions = c(0.10, 0.22, 0.4)   ,
          fontsize = 1.0                    ,
          refLabel = "1"                    ,
          noDigits =  4 )
dev.off()
sig1_APA_exp_sur_geneName$OST_days <-  sig1_APA_exp_sur_geneName$OST
sig1_APA_exp_sur_geneName$OST_month <- sig1_APA_exp_sur_geneName$OST_days/30
sig1_APA_exp_sur_geneName$OST_year <- sig1_APA_exp_sur_geneName$OST_days/365

Mulcox <- coxph(Surv(OST_days,OS) ~ GPNMB+COL1A1, data =  sig1_APA_exp_sur_geneName)
Mulcox_coef <- coef(Mulcox)
Mulcox_Coxp <- signif(as.numeric(summary(Mulcox)$coefficients[,c( "Pr(>|z|)")]),digit=4)
Mulcox_coef_used = Mulcox_coef[Mulcox_Coxp < 0.05 ] # 剩下3个
Mulcox_coef_used
# OST 单位 天 月 年 对model 的系数没有影响
Mulcox_coef_used <- Mulcox_coef_used %>% as.data.frame() %>% set_colnames("weight")
Mulcox_coef_used$GeneName  <-  rownames(  Mulcox_coef_used  )
Mulcox_coef_used$APAevent  <-  sig_gene_multi_cox[ c(1 , 3 )]
Mulcox_coef_used %>% readr::write_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/TCGA_SKCM_APA_metastasis_OS_survival_model_GPNMB_COL1A1_APA.rds")
#rownames(Mulcox_coef_used)  <-  Mulcox_coef_used$APAevent
sig1_APA_exp_sur_geneName$ImmAPA_SurvScore  <- sig1_APA_exp_sur_geneName %>% dplyr::select(Mulcox_coef_used$GeneName) %>% apply(1,function(x){
  x <- as.data.frame(x) %>% t() %>% as.data.frame()
  apply(x,1,function(x1){
    Mulcox_coef_used[names(x1),]$weight*as.numeric(x1)
  }) %>% sum -> ImmAPA_SurvScore
  return(ImmAPA_SurvScore)
})


Mulcox1 <- coxph(Surv(OST_days,OS) ~ ImmAPA_SurvScore+Gender+Age+Gender+Tumor_stage, data =  sig1_APA_exp_sur_geneName)

sig1_APA_exp_sur_geneName <- sig1_APA_exp_sur_geneName %>% mutate(Response=ifelse(Benifit %in% "NonBenifit","NonResponse",NA))

# compare model and other effected factor
# OST 单位 天 月 年 对model 的系数没有影响
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/Hazard ratio of GPNMB+COL1A1 model ImmAPA_SurvScore+Gender+Age+Gender+Tumor_stage.pdf",width = 10,height = 6,onefile = FALSE)
ggforest( model = Mulcox1                    ,
          data  = sig1_APA_exp_sur_geneName ,
          main  = "Hazard ratio"            ,
          cpositions = c(0.10, 0.22, 0.4)   ,
          fontsize = 1.0                    ,
          refLabel = "1"                    ,
          noDigits =  4 )+labs(title = "OST_months")
dev.off()


# model validation in the training dataset
aa_try <- try(aa <- maxstat.test(survival::Surv(OST, OS) ~ ImmAPA_SurvScore, data= sig1_APA_exp_sur_geneName,smethod="LogRank"),silent=TRUE)
if (is(aa_try , "try-error")) {return(NA)} else { sig1_APA_exp_sur_geneName <- sig1_APA_exp_sur_geneName  %>% 
  dplyr::mutate(ImmAPA_SurvScore_bestCut = ifelse(ImmAPA_SurvScore > aa_try$estimate,"High","Low"))}
model1 <- survival::coxph(survival::Surv(OST,OS) ~ ImmAPA_SurvScore, data=sig1_APA_exp_sur_geneName,na.action=na.exclude)
HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
HR_detail = summary(model1)
CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
fit<- survival::survfit(survival::Surv(OST, OS) ~ ImmAPA_SurvScore_bestCut,data=sig1_APA_exp_sur_geneName)
model1 <- survival::survdiff(survival::Surv(OST, OS) ~ ImmAPA_SurvScore_bestCut,data=sig1_APA_exp_sur_geneName, na.action=na.exclude)
KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(sig1_APA_exp_sur_geneName$ImmAPA_SurvScore_bestCut)))-1)

p5 <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                            pval = F, title=paste("TCGA SKCM metastasis : p = ",signif(as.numeric(KMP),digits = 2),
                                                  "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                            fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                            legend.labs = c(paste("High (n=",nrow(sig1_APA_exp_sur_geneName[sig1_APA_exp_sur_geneName$ImmAPA_SurvScore_bestCut %in% "High",]),")",sep=""),
                                            paste("Low (n=",nrow(sig1_APA_exp_sur_geneName[sig1_APA_exp_sur_geneName$ImmAPA_SurvScore_bestCut %in% "Low",]),")",sep="")))
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/ImmAPA_SurvScore os survival TCGA_SKCM_APA_metastasis.pdf",width = 7,height = 6,onefile = FALSE)
p5
dev.off()

# validation for FPS
sig1_APA_exp_sur_geneName$PFST <- sig1_APA_exp_sur_geneName$PFST/30
aa_try <- try(aa <- maxstat.test(survival::Surv(PFST, PFS) ~ ImmAPA_SurvScore, data= sig1_APA_exp_sur_geneName,smethod="LogRank"),silent=TRUE)
if (is(aa_try , "try-error")) {return(NA)} else { sig1_APA_exp_sur_geneName <- sig1_APA_exp_sur_geneName  %>% 
  dplyr::mutate(ImmAPA_SurvScore_bestCut = ifelse(ImmAPA_SurvScore > aa_try$estimate,"High","Low"))}
model1 <- survival::coxph(survival::Surv(PFST,PFS) ~ ImmAPA_SurvScore, data=sig1_APA_exp_sur_geneName,na.action=na.exclude)
HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
HR_detail = summary(model1)
CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") #Result[i1, c('CI_95%_for_HR')]
fit<- survival::survfit(survival::Surv(PFST, PFS) ~ ImmAPA_SurvScore_bestCut,data=sig1_APA_exp_sur_geneName)
model1 <- survival::survdiff(survival::Surv(PFST, PFS) ~ ImmAPA_SurvScore_bestCut,data=sig1_APA_exp_sur_geneName, na.action=na.exclude)
KMP <- 1-pchisq(model1$chisq, df=length(levels(factor(sig1_APA_exp_sur_geneName$ImmAPA_SurvScore_bestCut)))-1)

p5 <- survminer::ggsurvplot(fit,palette = c("red", "blue"),
                            pval = F, title=paste("TCGA SKCM metastasis : p = ",signif(as.numeric(KMP),digits = 2),
                                                  "\nHR ",HR, " [95% CI:",CI,"]",sep=""),
                            fun = function(y) y*100,legend = c(0.2, 0.2), legend.title=F,font.legend=10,font.title=10,
                            legend.labs = c(paste("High (n=",nrow(sig1_APA_exp_sur_geneName[sig1_APA_exp_sur_geneName$ImmAPA_SurvScore_bestCut %in% "High",]),")",sep=""),
                                            paste("Low (n=",nrow(sig1_APA_exp_sur_geneName[sig1_APA_exp_sur_geneName$ImmAPA_SurvScore_bestCut %in% "Low",]),")",sep="")))
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/ImmAPA_SurvScore PFS survival TCGA_SKCM_APA_metastasis.pdf",width = 7,height = 6,onefile = FALSE)
p5
dev.off()


# predicted benifit validation
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/ImmAPA_SurvScore benifit diff onside test in metastasis SKCM.pdf",width=5,height=5)
lapply(split(sig1_APA_exp_sur_geneName,sig1_APA_exp_sur_geneName$Benifit), function(x){x$ImmAPA_SurvScore}) -> Benifit_ImmAPA_SurvScore
t.test(Benifit_ImmAPA_SurvScore$Benifit,Benifit_ImmAPA_SurvScore$NonBenifit,alternative = c( "less")) -> t_test_less
wilcox.test(Benifit_ImmAPA_SurvScore$Benifit,Benifit_ImmAPA_SurvScore$NonBenifit,alternative = c( "less")) -> wilcox_test_less

ggplot(sig1_APA_exp_sur_geneName,aes(y=ImmAPA_SurvScore,x=Benifit))+
  geom_boxplot()+geom_jitter(width = 0.2)+
  labs(title = paste("P value=",signif(t_test_less$p.value,digits = 4)))+
  theme1
dev.off()

# analysis the model and other ICB related signature
TIFScore <- readr::read_rds("/work/yye/AliRstudio/Project/2018Immunotherapy/ProcessedData/TIFScore.rds.gz")
GEPScore <- readr::read_rds("/work/yye/AliRstudio/Project/2018Immunotherapy/ProcessedData/GEPScore.rds.gz")
SKCM_TIF <- TIFScore[TIFScore$cancer_types %in% "SKCM",]$TIFscore[[1]]
SKCM_GEP <- GEPScore[TIFScore$cancer_types %in% "SKCM",]$GEPscore[[1]]

SKCM_Metastasis_ImmAPA_SURV <- merge(sig1_APA_exp_sur_geneName,SKCM_GEP,by.x="Row.names",by.y="barcode")
SKCM_Metastasis_ImmAPA_SURV <- merge(SKCM_Metastasis_ImmAPA_SURV,SKCM_TIF,by.x="Row.names",by.y="barcode")

corr <- my.cor.test_spearman(SKCM_Metastasis_ImmAPA_SURV$ImmAPA_SurvScore,SKCM_Metastasis_ImmAPA_SURV$GEP)
p1 <- ggplot(SKCM_Metastasis_ImmAPA_SURV,aes(x=ImmAPA_SurvScore,y=GEP))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(title = paste("R =",round(as.numeric(corr$estimate),2),"  ","P value=",signif(as.numeric(corr$p.value),3)))+
  theme1

corr <- my.cor.test_spearman(SKCM_Metastasis_ImmAPA_SURV$ImmAPA_SurvScore,SKCM_Metastasis_ImmAPA_SURV$TIF)
p2 <- ggplot(SKCM_Metastasis_ImmAPA_SURV,aes(x=ImmAPA_SurvScore,y=TIF))+
  geom_point()+
  geom_smooth(method = "lm")+
  labs(title = paste("R =",round(as.numeric(corr$estimate),2),"  ","P value=",signif(as.numeric(corr$p.value),3)))+
  theme1
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/ImmAPA_SurvScore corr with GEP and TIF in metastasis SKCM.pdf",width=10,height=5)
cowplot::plot_grid(p1,p2)
dev.off()
SKCM_Metastasis_ImmAPA_SURV$Age 
