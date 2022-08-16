library(survivalROC)
library(survival)
library(reshape)
library(patchwork)
#计算 score的 AUC  c-index  HR response vs nonresponse 差异p值
#计算 score的 AUC  c-index  HR 差异p值
my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
my.wilcox.test.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#===1.读入 signature score
mRNA_exp=readRDS('"/work/gywang/project/Immune_APA/ImmAPA_ICB/data/GSE91061_metastasis_SKCM_24_ICB_related_signature.rds"')
mRNA_exp$clinical_ImmAPA[[1]] <- mRNA_exp$clinical_ImmAPA[[1]] %>% mutate(Benifit=.$Benifit_PRCR_PDSD)
mRNA_exp$clinical_ImmAPA[[2]] <- mRNA_exp$clinical_ImmAPA[[2]] %>% mutate(Benifit=.$Benefit)
mRNA_exp$clinical_ImmAPA[[2]] <- mRNA_exp$clinical_ImmAPA[[2]] %>% mutate(OST=.$OST/30)

#计算 score的 AUC  c-index  HR response vs nonresponse 差异p值
#计算 score的 AUC  c-index  HR 差异p值
my.t.test.p.value <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)}

my.wilcox.test.p.value <- function(...) {
  obj<-try(wilcox.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)}

# 整合几个score 在一起 注意 提取 pre 样本 否则样本数目不一致 会有问题
mRNA_exp$all_score <- lapply(c(1:2),function(xx){
  temp_score=t(do.call(rbind,	lapply(as.list(mRNA_exp[xx,c(3:26)]),function(xx1){xx1[[1]]})))%>%data.frame()
  temp_score
})
.x=mRNA_exp$all_score[[1]]
.y=mRNA_exp$clinical_ImmAPA[[1]]
x=1

# ===2.计算 OS的 AUC c-index HR
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  AUC_HR_c_index_OS = purrr::pmap(list(.x=all_score,.y=clinical_ImmAPA),function(.x,.y){
    try(  {.x=data.frame(.x)
    tempexp=.x
    temp_sur=.y
    #3.对于每个score
    lapply(c(1:ncol(tempexp)),function(x){
      #score
      gene <- tempexp[,x]
      gene <- unlist(gene)
      group <- ifelse(gene > median(gene,na.rm = T),'high', 'low')
      #=======================================1.统计score 和 各临床信息
      survival_dat <- data.frame(group = group,  #画 MK用这组 信息
                                 status = temp_sur$OS,
                                 time = temp_sur$OST,
                                 gene = gene,
                                 benefit=temp_sur$Benifit,
                                 stringsAsFactors = F)
      survival_dat$samples=temp_sur$sample_use# rownames(survival_dat)				 
      
      #1.AUC
      cutoff=12*1 #12个 月 #31年
      nobs <- length(survival_dat)
      PFS_sur_roc_year1= survivalROC(Stime=survival_dat$time,##生存时间
                                     status=survival_dat$status,## 终止事件    
                                     marker = survival_dat$gene, ## marker value    
                                     predict.time = 12,## 预测时间截点
                                     method="KM")	                                     
      PFS_sur_roc_year2= survivalROC(Stime=survival_dat$time,##生存时间
                                     status=survival_dat$status,## 终止事件    
                                     marker = survival_dat$gene, ## marker value    
                                     predict.time = 24,## 预测时间截点
                                     method="KM")	 									   #span = 0.25*nobs^(-0.20))##span,NNE法的namda
      PFS_sur_roc_year3= survivalROC(Stime=survival_dat$time,##生存时间
                                     status=survival_dat$status,## 终止事件    
                                     marker = survival_dat$gene, ## marker value    
                                     predict.time = 36,## 预测时间截点
                                     method="KM")	
      ROC=data.frame(AUC=c(PFS_sur_roc_year1$AUC,PFS_sur_roc_year2$AUC,PFS_sur_roc_year3$AUC))
      
      if(F){
        #PFS_sur_roc_year1$AUC
        if(nrow(temp_sur)!=70){ #该套数据时间最大为32 不到 36月
          ROC<- timeROC(T=survival_dat$time,#结局时间 
                        delta=survival_dat$status,#生存结局 
                        marker=survival_dat$gene,#预测变量 
                        cause=1,#阳性结局赋值，比如死亡与否
                        weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
                        times=c(12,24,36),#时间点，选取5年(60个月)和8年生存率  # SKCM_metastasis 三年用 30
                        ROC = TRUE,
                        iid = TRUE)        
          #ROC$AUC
        }else if(nrow(temp_sur)==70){
          ROC<- timeROC(T=survival_dat$time,#结局时间 
                        delta=survival_dat$status,#生存结局 
                        marker=survival_dat$gene,#预测变量 
                        cause=1,#阳性结局赋值，比如死亡与否
                        weighting="marginal",#权重计算方法，marginal是默认值，采用km计算删失分布
                        times=c(12,24,30),#时间点，选取5年(60个月)和8年生存率  # SKCM_metastasis 三年用 30
                        ROC = TRUE,
                        iid = TRUE) }}
      
      #2.c-index
      fit <- coxph(Surv(time,status)~gene,data=survival_dat)
      sum.surv <- summary(fit)
      c_index <-sum.surv$concordance
      
      #3.HR coxp
      model1 <- survival::coxph(Surv(time, status) ~ gene, data = survival_dat,  na.action=na.exclude)
      coef = signif(as.numeric(summary(model1)$coefficients[1,c( "coef")]),digit=4)		
      HR = signif(as.numeric(summary(model1)$coefficients[1,c( "exp(coef)")]),digit=4)
      Coxp <- signif(as.numeric(summary(model1)$coefficients[1,c( "Pr(>|z|)")]),digit=4)
      HR_detail = summary(model1)
      CILow =  HR_detail$conf.int[,"lower .95"]
      CIHigh =  HR_detail$conf.int[,"upper .95"]
      CI =  paste(signif(HR_detail$conf.int[,"lower .95"],2),'-',signif(HR_detail$conf.int[,"upper .95"],2),sep=" ") 			
      HR_se=sqrt(model1$var)
      
      #4.得分在 CR_PR-response 和 PD-nonresponse 的差异 p值
      response_P=c( my.t.test.p.value(survival_dat[survival_dat$benefit %in% c("Benifit"),]$gene,survival_dat[survival_dat$benefit %in% c("NonBenifit"),]$gene),
                    my.wilcox.test.p.value(survival_dat[survival_dat$benefit %in% c("Benifit"),]$gene,survival_dat[survival_dat$benefit %in% c("NonBenifit"),]$gene))
      
      #4.2 得分在 CR_PR-response 和 PD-nonresponse 的差异 logistic p值
      fit<-glm(as.factor(benefit) ~ gene ,data=survival_dat,family=binomial(link="logit"))
      response_P2=summary(fit)
      rocobj <- roc(survival_dat$benefit, survival_dat$gene)
      Signature_Benifit_AUC <- rocobj$auc[[1]]
      
      
      
      
      c(ROC$AUC, # ROC
        ROC$AUC, # ROC$inference$vect_sd_1,#ROC的SE
        c_index, # c-index 及 se
        HR,HR_se,Coxp,coef,CILow,CIHigh,response_P,response_P2$coefficients[-1,][4],Signature_Benifit_AUC) %>% 
        set_names(c('year1','year2','year3','year1_Se','year2_Se','year3_Se','C_index','C-se','HR','HR_se','Coxp',
                    'coef','CILow','CIHigh','benifit_diff_testP','benifit_diff_wilcoxP','logisticP',"Signature_Benifit_AUC"))
    })%>% bind_rows()%>% data.frame()
    },silent = T)
  }))


#===3.可视化 OS AUC c-index HR 差异p值
AUC_HR_index_all_OS=rbind(mRNA_exp$AUC_HR_c_index_OS[[1]],mRNA_exp$AUC_HR_c_index_OS[[2]])
AUC_HR_index_all_OS$dataset=c(rep('GSE91061',nrow(mRNA_exp$AUC_HR_c_index_OS[[1]])),rep('SKCM_metastasis',nrow(mRNA_exp$AUC_HR_c_index_OS[[2]])))
AUC_HR_index_all_OS$sig=rep(colnames(mRNA_exp$all_score[[1]]),2)
AUC_HR_index_all_OS$AUC_mean=apply(AUC_HR_index_all_OS[,c('year1','year2','year3')],1,mean)
write.table(AUC_HR_index_all_OS,file="/work/gywang/project/Immune_APA/ImmAPA_ICB/data/multiple_signature_compear/GSE91061_metastasis_SKCM_24_ICB_related_signature_AUC_HR_index_all_OS.tab",sep="\t",row.names = F,quote = F)

#去掉 blood
AUC_HR_index_all_OS=subset(AUC_HR_index_all_OS,!(sig%in%c('Blood','CGAs','IRGs')))
#设置三套数据展示顺序
AUC_HR_index_all_OS$dataset=factor(AUC_HR_index_all_OS$dataset, levels=c('GSE91061','SKCM_metastasis')) 


#3.2 HR barplot
# 由小到大排序 reorder(Samples,Freq)
x_order_HR <- lapply(split(AUC_HR_index_all_OS[AUC_HR_index_all_OS$dataset %in% "GSE91061",]$HR,AUC_HR_index_all_OS$sig), mean) %>% unlist()
x_order_HR <- names(x_order_HR )[order(x_order_HR ,decreasing =T)]

ggHR_OS=ggplot(arrange(AUC_HR_index_all_OS,HR),aes(x=sig,y=HR,fill=sig))+
  geom_bar(stat = "identity",position="stack",color=NA,width = 0.6)+
  facet_grid(dataset~.,scale="free_y")+
  scale_x_discrete(limit=x_order_HR,label=x_order_HR)+
  theme(panel.background = element_rect(fill=NA,color="black"),axis.title.x = element_blank(),panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=90,hjust = 1),axis.text.y = element_text(color="black"),axis.ticks.x = element_blank(),
        legend.position='none')+
  labs(ylab='HR') +
  geom_hline(aes(yintercept=1), colour="#990000", linetype="dashed")

ggsave(paste('/work/gywang/project/Immune_APA/ImmAPA_ICB/multiple_signature_compear/','',"ggggHR_OS_compare_barplot.pdf",sep=''),ggHR_OS,width=1.5,height=4)


x_order_Signature_Benifit_AUC <- lapply(split(AUC_HR_index_all_OS[AUC_HR_index_all_OS$dataset %in% "GSE91061",]$Signature_Benifit_AUC,AUC_HR_index_all_OS$sig), mean) %>% unlist()
x_order_Signature_Benifit_AUC  <- names(x_order_Signature_Benifit_AUC  )[order(x_order_Signature_Benifit_AUC,decreasing =F)]

ggSignature_Benifit_AUC_OS=ggplot(AUC_HR_index_all_OS,aes(x=sig,y=Signature_Benifit_AUC,fill=sig))+
  facet_grid(.~dataset,scale="free_y")+
  coord_flip()+
  geom_bar(stat = "identity",position="stack",color=NA,width = 0.6)+
  scale_x_discrete(limit=x_order_Signature_Benifit_AUC,label=x_order_Signature_Benifit_AUC)+
  theme(panel.background = element_rect(fill=NA,color="black"),axis.title.x = element_blank(),panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=0,hjust = 1),axis.text.y = element_text(color="black"),axis.ticks.x = element_blank(),
        legend.position='none')+
  labs(ylab='HR') +
  geom_hline(aes(yintercept=1), colour="#990000", linetype="dashed")
ggSignature_Benifit_AUC_OS
ggsave(paste('/work/gywang/project/Immune_APA/ImmAPA_ICB/multiple_signature_compear/','',"ggSignature_Benifit_AUC_barplot.pdf",sep=''),ggSignature_Benifit_AUC_OS,width=4,height=4)




#==3.3 AUC
AUC_HR_index_all_OS_long=melt(AUC_HR_index_all_OS, id.vars=c('dataset','sig'),measure.vars=c('year1','year2','year3'))
AUC_HR_index_all_OS_long$variable=as.character(AUC_HR_index_all_OS_long$variable)
#设置三套数据展示顺序
AUC_HR_index_all_OS$dataset=factor(AUC_HR_index_all_OS$dataset, levels=c('GSE91061', 'SKCM_metastasis')) 
AUC_HR_index_all_OS$AUC_mean=apply(AUC_HR_index_all_OS[,c('year1','year2','year3')],1,mean)

x_order_AUC <- lapply(split(AUC_HR_index_all_OS[AUC_HR_index_all_OS$dataset %in% "GSE91061",]$AUC_mean,AUC_HR_index_all_OS$sig), mean) %>% unlist()
x_order_AUC <- names(x_order_AUC)[order(x_order_AUC,decreasing =T)]

#===3.3.4    3个 AUC boxplot + 散点
#先画 bat

#boxplot 横过来
ggAUC_OS_box2=ggplot(AUC_HR_index_all_OS_long,aes(x=sig,y=value))+
  facet_grid(.~dataset,scale="free_y")+ #coord_flip()+
  coord_flip()+
  geom_boxplot(width=0.6,color =NA,outlier.shape = NA,fill="gray")+ #fill=NA,
  ggbeeswarm::geom_beeswarm(alpha=0.5,size=0.5)+  
  scale_x_discrete(limit=rev(x_order_AUC),label=rev(x_order_AUC))+
  theme(panel.background = element_rect(fill=NA,color="black"),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=0),axis.text.y = element_text(color="black"),axis.ticks.x = element_blank(),
        legend.position='none')+ 
  labs(ylab='AUC')  

#再改点颜色
ggAUC_OS_box2=ggAUC_OS_box2 + geom_point(data=AUC_HR_index_all_OS,aes( x= sig, y =  year1),color = 'blue',height = 0.4,size= 0.8)+ #画 多个点
  geom_point(data=AUC_HR_index_all_OS,aes( x= sig, y =  year2),color = 'red',height = 0.4,size= 0.8)+ #画 多个点
  geom_point(data=AUC_HR_index_all_OS,aes( x= sig, y =  year3),color = 'green',height = 0.4,size= 0.8) #画 多个点
ggAUC_OS_box2
ggsave(paste('/work/gywang/project/Immune_APA/ImmAPA_ICB/multiple_signature_compear/','',"gg_OS_AUC_boxplot2.pdf",sep=''),ggAUC_OS_box2,width=4,height=4)


#==3.4.1 wilcox 差异 p值  barplot
ggttest_OS=ggplot(AUC_HR_index_all_OS,aes(x=sig,y=-log2(benifit_diff_wilcoxP),fill=sig))+
  geom_bar(stat = "identity",position="stack",color=NA,width = 0.6)+
  facet_grid(dataset~.,scale="free_y")+
  scale_x_discrete(limit=x_order_AUC,label=x_order_AUC)+
  theme(panel.background = element_rect(fill=NA,color="black"),axis.title.x = element_blank(),panel.grid.major = element_blank(),
        axis.text.x = element_text(angle=45),axis.text.y = element_text(color="black"),axis.ticks.x = element_blank(),
        legend.position='none')+
  labs(ylab='HR') +
  geom_hline(aes(yintercept=4.32),colour="#990000", linetype="dashed")
ggsave(paste('/work/gywang/project/Immune_APA/ImmAPA_ICB/multiple_signature_compear/','',"ggttest_OS_benifit_compare_barplot.pdf",sep=''),ggttest_OS,width=8,height=4)

#==3.4.2 CR - PD 之间差异 p 值 wilcox 差异 p值  散点图
ggttest_OS_ponit=ggplot(AUC_HR_index_all_OS,aes(x=-log10(benifit_diff_wilcoxP),y=sig))+
  xlab("-log10P-Value")+ylab("sig")+
  geom_point(aes(shape=dataset,color=sig),alpha=0.8)+ #size=-log10(CRwilcoxP),
  scale_y_discrete(limit=rev(x_order_AUC),label=rev(x_order_AUC))+ theme_bw()+
  geom_vline(xintercept = 1.30103, colour=colors()[114], linetype="dashed")+
  theme(legend.position ="none",panel.background = element_blank())#theme_linedraw,panel.grid=element_blank(),
ggsave(paste('/work/gywang/project/Immune_APA/ImmAPA_ICB/multiple_signature_compear/','',"signature diffence significance wilcoxP test.pdf",sep=''),ggttest_OS_ponit,width=4,height=4)

ggttest_OS_ponit=ggplot(AUC_HR_index_all_OS,aes(x=-log10(benifit_diff_testP),y=sig))+
  xlab("-log10P-Value")+ylab("sig")+
  geom_point(aes(shape=dataset,color=sig),alpha=0.8)+ #size=-log10(CRwilcoxP),
  scale_y_discrete(limit=rev(x_order_AUC),label=rev(x_order_AUC))+ theme_bw()+
  geom_vline(xintercept = 1.30103, colour=colors()[114], linetype="dashed")+
  theme(legend.position ="none",panel.background = element_blank())#theme_linedraw,panel.grid=element_blank(),
ggsave(paste('/work/gywang/project/Immune_APA/ImmAPA_ICB/multiple_signature_compear/','',"signature diffence significance T test.pdf",sep=''),ggttest_OS_ponit,width=4,height=4)

ggttest_OS_ponit=ggplot(AUC_HR_index_all_OS[AUC_HR_index_all_OS$dataset %in% "GSE91061",],aes(x=-log10(benifit_diff_testP),y=sig))+
  xlab("-log10P-Value")+ylab("sig")+
  geom_point(aes(shape=dataset,color=sig),alpha=0.8)+ #size=-log10(CRwilcoxP),
  scale_y_discrete(limit=rev(x_order_AUC),label=rev(x_order_AUC))+ theme_bw()+
  geom_vline(xintercept = 1.30103, colour=colors()[114], linetype="dashed")+
  theme(legend.position ="none",panel.background = element_blank())#theme_linedraw,panel.grid=element_blank(),
ggttest_OS_ponit
ggttest_OS_ponit=ggplot(AUC_HR_index_all_OS[AUC_HR_index_all_OS$dataset != "GSE91061",],aes(x=-log10(benifit_diff_testP),y=sig))+
  xlab("-log10P-Value")+ylab("sig")+
  geom_point(aes(shape=dataset,color=sig),alpha=0.8)+ #size=-log10(CRwilcoxP),
  scale_y_discrete(limit=rev(x_order_AUC),label=rev(x_order_AUC))+ theme_bw()+
  geom_vline(xintercept = 1.30103, colour=colors()[114], linetype="dashed")+
  theme(legend.position ="none",panel.background = element_blank())#theme_linedraw,panel.grid=element_blank(),
ggttest_OS_ponit

#==3.4.3 CR - PD 之间差异 logistic p值  散点图
ggLogisticP_OS_ponit=ggplot(AUC_HR_index_all_OS,aes(x=-log10(logisticP),y=sig))+
  xlab("-log10P-Value")+ylab("sig")+#coord_fixed()+
  geom_point(aes(shape=dataset,color=sig),alpha=0.8)+ #size=-log10(logisticP),
  scale_y_discrete(limit=rev(x_order_AUC),label=rev(x_order_AUC))+ theme_bw()+
  #p 值 阈值
  geom_vline(xintercept = 1.30103, colour=colors()[114], linetype="dashed")+
  theme(legend.position ="none",panel.background = element_blank())#theme_linedraw,panel.grid=element_blank(),
ggsave(paste('/work/gywang/project/Immune_APA/ImmAPA_ICB/multiple_signature_compear/','',"ggLogisticP_OS_ponit.pdf",sep=''),ggLogisticP_OS_ponit,width=4,height=4)

# 先在12.3.compute_other_signature_score_HR_metaAnalysis_OS.r 中 计算出 HR 的 点图 HR_geom_point_plot
# AUC boxplot + response 差异 p值 + HR 的 点图 
ggAUC_OS_box2_ttest_point=ggAUC_OS_box2+ggLogisticP_OS_ponit+HR_geom_point_plot + plot_layout(ncol = 3,widths=c(1.5,1,1)) 
ggsave(paste('/work/gywang/project/Immune_APA/ImmAPA_ICB/multiple_signature_compear/','',"ggLogistic_OS_ggttest_OS_ponit_AUC_box.pdf",sep=''),ggAUC_OS_box2_ttest_point,width=8,height=4)


