library(glmnet)
library(tidyverse)
library(magrittr)
library(maxstat)
library(survival) 
library(survminer)
library(timeROC)
options(stringsAsFactors = F)

theme1 <-   theme(  panel.background=element_rect(colour=NA,fill="white"),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    axis.title.x=element_blank(),
                    axis.text.x =element_text(color="black",size=10,angle = 0,hjust = 0.5),
                    axis.text.y = element_text(color="black",size=10,face = "bold"),
                    axis.line.y = element_line(colour = "black"),
                    axis.line.x = element_line(colour = "black"),
                    panel.border = element_rect(color = "black",fill = NA))
ICB_APA_data <- readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/ICB_APA_data.rds.gz")

Mulcox_coef_used <- readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/TCGA_SKCM_APA_metastasis_ImmAPA_survival_model_2_APA.rds")

GSE91061_OnTreatment_541ImmAPA_sur <- readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/GSE91061_on_treatment_541ImmAPA_and_survival.rds")
GSE91061_PreTreatment_541ImmAPA_sur <- readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/GSE91061_541ImmAPA_and_survival.rds")
ICB_melanoma_Xiangya_541ImmAPA_sur <- readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/ICB_melanoma_Xiangya_541ImmAPA_and_survival.rds")
TCGA_SKCM_metastasis_541ImmAPA_sur <- readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/TCGA_SKCM_metastasis_541ImmAPA_and_survival.rds")


sig_APA_exp <- GSE91061_OnTreatment_541ImmAPA_sur
sig_APA_exp_sur <- sig_APA_exp %>% dplyr::select(-setdiff(intersect(ImmAPA541_APAevents,colnames(sig_APA_exp)),Mulcox_coef_used$APAevent))
sig_APA_exp_sur$ImmAPA_survScore  <- sig_APA_exp_sur %>% dplyr::select(Mulcox_coef_used$APAevent) %>% apply(1,function(x){
  x <- as.data.frame(x) %>% t() %>% as.data.frame()
  apply(x,1,function(x1){
    Mulcox_coef_used[names(x1),]$weight*as.numeric(x1)
  }) %>% sum -> ImmAPA_surv
  return(ImmAPA_surv)
})
sig_APA_exp_sur$Benifit_N <- ifelse(sig_APA_exp_sur$Response %in% c("CR","PR"),"CR/PR",as.character( sig_APA_exp_sur$Response))
sig_APA_exp_sur_CR_PR_PD <- dplyr::filter(sig_APA_exp_sur,Benifit_N %in% c("CR/PR","PD"))
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/TCGA_SKCM_APA_metastasis_OS_model GPNMB COL1A1 validate in GSE91061 on tretment benifit.pdf",width = 10,height = 5,onefile = FALSE)
p1 <- ggplot(sig_APA_exp_sur,aes(x=Benifit,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 on-treatment")+
  theme1
p2 <- ggplot(sig_APA_exp_sur_CR_PR_PD,aes(x=Benifit_N,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 on-treatment")+
  theme1
p3 <- ggplot(sig_APA_exp_sur,aes(x=Response,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 on-treatment")+
  theme1
cowplot::plot_grid(plotlist = list(p1,p2,p3),ncol = 3)
dev.off()

sig_APA_exp <- GSE91061_PreTreatment_541ImmAPA_sur
sig_APA_exp_sur <- sig_APA_exp %>% dplyr::select(-setdiff(intersect(ImmAPA541_APAevents,colnames(sig_APA_exp)),Mulcox_coef_used$APAevent))
sig_APA_exp_sur$ImmAPA_survScore  <- sig_APA_exp_sur %>% dplyr::select(Mulcox_coef_used$APAevent) %>% apply(1,function(x){
  x <- as.data.frame(x) %>% t() %>% as.data.frame()
  apply(x,1,function(x1){
    Mulcox_coef_used[names(x1),]$weight*as.numeric(x1)
  }) %>% sum -> ImmAPA_surv
  return(ImmAPA_surv)
})
sig_APA_exp_sur$Benifit_N <- ifelse(sig_APA_exp_sur$Response %in% c("CR","PR"),"CR/PR",as.character( sig_APA_exp_sur$Response))
sig_APA_exp_sur_CR_PR_PD <- dplyr::filter(sig_APA_exp_sur,Benifit_N %in% c("CR/PR","PD"))
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/TCGA_SKCM_APA_metastasis_OS_model GPNMB COL1A1 validate in ICB_melanoma_Xiangya benifit.pdf",width = 10,height = 5,onefile = FALSE)
p1 <- ggplot(sig_APA_exp_sur,aes(x=Benifit,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "ICB_melanoma_Xiangya")+
  theme1
p2 <- ggplot(sig_APA_exp_sur_CR_PR_PD,aes(x=Benifit_N,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "ICB_melanoma_Xiangya")+
  theme1
p3 <- ggplot(sig_APA_exp_sur,aes(x=Response,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "ICB_melanoma_Xiangya")+
  theme1
cowplot::plot_grid(plotlist = list(p1,p2,p3),ncol = 3)
dev.off()


sig_APA_exp <- ICB_melanoma_Xiangya_541ImmAPA_sur
sig_APA_exp_sur <- sig_APA_exp %>% dplyr::select(-setdiff(intersect(ImmAPA541_APAevents,colnames(sig_APA_exp)),Mulcox_coef_used$APAevent))
sig_APA_exp_sur$ImmAPA_survScore  <- sig_APA_exp_sur %>% dplyr::select(Mulcox_coef_used$APAevent) %>% apply(1,function(x){
  x <- as.data.frame(x) %>% t() %>% as.data.frame()
  apply(x,1,function(x1){
    Mulcox_coef_used[names(x1),]$weight*as.numeric(x1)
  }) %>% sum -> ImmAPA_surv
  return(ImmAPA_surv)
})
sig_APA_exp_sur$Benifit_N <- ifelse(sig_APA_exp_sur$Response %in% c("CR","PR"),"CR/PR",as.character( sig_APA_exp_sur$Response))
sig_APA_exp_sur_CR_PR_PD <- dplyr::filter(sig_APA_exp_sur,Benifit_N %in% c("CR/PR","PD"))
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/TCGA_SKCM_APA_metastasis_OS_model GPNMB COL1A1 validate in GSE91061 Pre tretment benifit.pdf",width = 10,height = 5,onefile = FALSE)
p1 <- ggplot(sig_APA_exp_sur,aes(x=Benifit,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 Pre-treatment")+
  theme1
p2 <- ggplot(sig_APA_exp_sur_CR_PR_PD,aes(x=Benifit_N,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 Pre-treatment")+
  theme1
p3 <- ggplot(sig_APA_exp_sur,aes(x=Response,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 Pre-treatment")+
  theme1
cowplot::plot_grid(plotlist = list(p1,p2,p3),ncol = 3)
dev.off()


Mulcox_coef_used <- TCGA_SKCM_metastasis_OS_model_2_APA <- readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/TCGA_SKCM_APA_metastasis_OS_survival_model_2_APA.rds")

sig_APA_exp <- TCGA_SKCM_metastasis_541ImmAPA_sur
sig_APA_exp_sur <- sig_APA_exp %>% dplyr::select(-setdiff(intersect(ImmAPA541_APAevents,colnames(sig_APA_exp)),Mulcox_coef_used$APAevent))
sig_APA_exp_sur$ImmAPA_survScore  <- sig_APA_exp_sur %>% dplyr::select(Mulcox_coef_used$APAevent) %>% apply(1,function(x){
  x <- as.data.frame(x) %>% t() %>% as.data.frame()
  apply(x,1,function(x1){
    Mulcox_coef_used[names(x1),]$weight*as.numeric(x1)
  }) %>% sum -> ImmAPA_surv
  return(ImmAPA_surv)
})
sig_APA_exp_sur$Benifit_N <- ifelse(sig_APA_exp_sur$Response %in% c("CR","PR"),"CR/PR",as.character( sig_APA_exp_sur$Response))
sig_APA_exp_sur_CR_PR_PD <- dplyr::filter(sig_APA_exp_sur,Benifit_N %in% c("CR/PR","PD"))
pdf("/work/gywang/project/Immune_APA/ImmAPA_ICB/figure/TCGA_SKCM_APA_metastasis_model/TCGA_SKCM_APA_metastasis_OS_model GPNMB COL1A1 validate in GSE91061 Pre tretment benifit.pdf",width = 10,height = 5,onefile = FALSE)
p1 <- ggplot(sig_APA_exp_sur,aes(x=Benifit,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 Pre-treatment")+
  theme1
p2 <- ggplot(sig_APA_exp_sur_CR_PR_PD,aes(x=Benifit_N,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 Pre-treatment")+
  theme1
p3 <- ggplot(sig_APA_exp_sur,aes(x=Response,y=ImmAPA_survScore))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = "GSE91061 Pre-treatment")+
  theme1
cowplot::plot_grid(plotlist = list(p1,p2,p3),ncol = 3)
dev.off()
