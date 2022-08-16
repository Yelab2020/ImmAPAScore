library(magrittr)

# mRNA expression
TCGA_FPKM <- readRDS("/work/gywang/project/Immune_APA/Data/TCGA_Data/TCGA_33CanFPKM.rds.gz")
TCGA_SKCM_metastasis_FPKM <- TCGA_FPKM[TCGA_FPKM$cancertypes %in% "SKCM",]$Expr[[1]] %>% 
  dplyr::select(Gene,gsub("-",".",TCGA_SKCM_metastasis_ImmAPASig$bcr_patient_barcode))
colnames(TCGA_SKCM_metastasis_FPKM) <- gsub("\\.","-",colnames(TCGA_SKCM_metastasis_FPKM))
TCGA_SKCM_metastasis_FPKM <- TCGA_SKCM_metastasis_FPKM %>% set_rownames(TCGA_SKCM_metastasis_FPKM$Gene) %>%  dplyr::select(-Gene)

GSE91061_AllData <- readr::read_rds("/work/gywang/project/Immune_APA/Data/ICB/GSE91061/GSE91061_AllData.rds.gz")
GSE91061_GeneExp <- GSE91061_AllData$GeneExp[[1]] %>% dplyr::select(Gene.Name,GSE91061_ImmAPASig$Row.names)
GSE91061_GeneExp <- GSE91061_GeneExp %>% set_rownames(GSE91061_GeneExp$Gene.Name) %>%  dplyr::select(-Gene.Name)

# ImmAPASigScore
COL1A1_GPNMB_model <- readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/TCGA_SKCM_APA_metastasis_OS_survival_model_GPNMB_COL1A1_APA.rds")
Validation_data_set <-  readr::read_rds("/work/gywang/project/Immune_APA/ImmAPA_ICB/data/Validation_data_set.rds")
GSE91061_ImmAPASig <- Validation_data_set$COL1A1_GPNMB_model_data[[1]]
GSE91061_ImmAPASig <- dplyr::filter(GSE91061_ImmAPASig,Subtype != "MUCOSAL")
TCGA_SKCM_metastasis_ImmAPASig <- Validation_data_set$COL1A1_GPNMB_model_data[[3]]

GSE91061_ImmAPASigScore <- GSE91061_ImmAPASig %>% dplyr::filter(Row.names %in% colnames(GSE91061_GeneExp)) %>% 
  dplyr::select(ImmAPA_survScore) %>% t() %>% as.numeric() %>% set_names(colnames(GSE91061_GeneExp))
TCGA_SKCM_metastasis_ImmAPASigScore <- TCGA_SKCM_metastasis_ImmAPASig %>% dplyr::filter(bcr_patient_barcode %in% colnames(TCGA_SKCM_metastasis_FPKM)) %>% 
  dplyr::select(ImmAPA_survScore) %>% t() %>% as.numeric() %>% set_names(colnames(TCGA_SKCM_metastasis_FPKM))


mRNA_exp <- tibble(dataset_name=c("GSE91061_without_MUCOSAL","TCGA_SKCM_metastasis"),
                   mRNA_FPKM=list(GSE91061_GeneExp,TCGA_SKCM_metastasis_FPKM),
                   ImmAPASigScore=list(GSE91061_ImmAPASigScore,TCGA_SKCM_metastasis_ImmAPASigScore))


#=====计算 其他signature的score

#====读入 signature 信息
ICB_signature <- openxlsx::read.xlsx("/work/dy/pub_data/signature/ICB/ICB_signature-input.xlsx")

#1.lncRNA_melanoma  ENSG00000279873 ENSG00000277767 没有表达
#保证 lncRNA 样本顺序 与 mRNA 相同

#2.checkpoint  加权求和  #输入表达谱 和 genelist
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  checkpoint = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        weight=c(-0.147,-0.135,0.283,-0.168,-0.560,0.773,-0.399)
        geneList=strsplit(ICB_signature$genes[[2]],',')[[1]];
        names(weight)=geneList
        interGenes=intersect(rownames(.x),names(weight));
        #交集的genelist和系数
        geneListInter=weight[interGenes]
        apply(.x[match(names(geneListInter),rownames(.x)),],2,function(x){sum(as.numeric(x)*as.numeric(geneListInter),na.rm=T)})		
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)

#3.TIR 4 gene  没看懂
#Thus, the TIR risk score of each melanoma patient was calculated using the “predict” function of R 
#software based on the expression levels (X) and the regression coefficients (β) of the four genes

#4.IMS=IFN-γ/IMS score两个得分比值 python-code 没看懂
#1.先用house gene 进行标准化 2.每个得分是标准化后的平均值   
#For the IFN-γ signature in this paper, we used the arithmetic mean of the log2-transformed, 
#house-keeping gene-normalized expression levels of the 10-gene “preliminary” IFN-γ signature 
#(IFNG, STAT1, CCR5, CXCL9, CXCL10, CXCL11, IDO1, PRF1, GZMA, and HLA-DRA)
IMS_set = c('CCL8', 'VCAN', 'CCL2', 'BCAT1', 'ISG15', 'CD163', 'AXL', 'CCL13', 'COL6A3', 
            'SIGLEC1', 'PDGFRB', 'IL10', 'STC1', 'ADAM12', 'OLFML2B', 'FAP', 'TWIST2', 'INHBA')
IFN_GAMMA_SET=c('IFNG', 'STAT1', 'CXCL9', 'CXCL10', 'IDO1', 'HLA-DRA')
House_keeping=c('ABCF1','DNAJC14','ERCC3','G6PD','GUSB','MRPL19','OAZ1','POLR2A','PSMC4',
                'PUM1', 'SDHA', 'SF3A1', 'STK11IP', 'TBC1D10B', 'TBP', 'TFRC', 'TLK2', 'TMUB2', 'UBB')
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  IMS = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        #log2
        .x=log2(.x+1)
        #house 每个样本计算 house 平均值
        housevalue=colMeans(.x[match(House_keeping,rownames(.x)),])
        #IFN_GAMMA_SET  每列数据除以 每个样本的housekeeping 值
        IFN_GAMMA_score=colMeans(apply(.x[match(IFN_GAMMA_SET,rownames(.x)),],2,function(x){x-housevalue}))
        #IMS_set  每列数据除以 每个样本的housekeeping 值
        IMS_set_score=colMeans(apply(.x[match(IMS_set,rownames(.x)),],2,function(x){x-housevalue}))
        IFN_GAMMA_score-IMS_set_score
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#5.IPRES summing the log2 Z scores of genes
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  IPRES = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[5]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){sum(log2(scale(as.numeric(x))),na.rm=T)})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#6.IFN-γ mean
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  IFN = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[6]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){mean(as.numeric(x),na.rm=T)})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#7.Exp. Immu. mean
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  Exp_Immu = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[7]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){mean(as.numeric(x),na.rm=T)})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#8.Roh Immu  geometric mean
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  Roh_Immu = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[2]]
        geneList=strsplit(ICB_signature$genes[[8]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){psych::geometric.mean(log2(x+0.1))})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#9.MessinaPrincipal component 1 score from PCA of expression levels of 12 chemokine signature genes

#10.IMPRES  Sum of ratios of 15 checkpoint pairs 两两基因对 加和
IMPRES_signature <- openxlsx::read.xlsx("/work/dy/pub_data/signature/ICB/IMPRES-genepairs.xlsx",sheet='input')
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  IMPRES = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        #log2转化 防止0 比值后为inf
        .x=log2(.x+0.1)
        exp1=.x[match(IMPRES_signature$Gene.1,rownames(.x)),]
        exp2=.x[match(IMPRES_signature$Gene.2,rownames(.x)),]
        #计算每个样本 gene1/gene2 比例 然后加和
        temp_ratios=lapply(c(1:ncol(.x)),function(x){
          sum(exp1[,x]/exp2[,x],na.rm=T)
        })%>%unlist
        names(temp_ratios)=colnames(.x)
        temp_ratios
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#11.NRS  ??????????
# mRNA_exp <- mRNA_exp %>% dplyr::mutate(
#   NRS = purrr::map(.x=mRNA_FPKM,function(.x){ 
#     try( #报错 不停止运行程序
#       {	
#         #.x=mRNA_exp$mRNA_FPKM[[1]]
#         geneList=strsplit(ICB_signature$genes[[11]],',')[[1]];
#         interGenes=intersect(rownames(.x),geneList);
#         #交集的genelist
#         apply(.x[match(interGenes,rownames(.x)),],2,function(x){psych::geometric.mean(x)})	
#       },silent = F)	
#   })) #%>% dplyr::select(-mRNA_FPKM)
# #saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")


#12.T eff. mean
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  T_eff = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[12]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){mean(as.numeric(x),na.rm=T)})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#13.Davoli mean
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  Davoli = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[13]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){mean(as.numeric(x),na.rm=T)})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#14.GEP  ssgsea  
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  GEP = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[14]],',')[[1]];
        xm <- apply(.x,2,function(x)(log2(x+1)))
        xm <- t(apply(xm, 1,function(x){(x-mean(x))/sd(x)}))
        #rownames(xm) <- .x$symbol
        xm_ES <- GSVA::gsva(xm,list(GEP=geneList) ,method= "ssgsea" , verbose=TRUE)
        xm_ES [1,]
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)

#16.CYT  
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  CYT4gene = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[13]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        CYT4gene=apply(.x[match(c("GZMB","PRF1","GZMA","GNLY"),rownames(.x)),],2,function(x){psych::geometric.mean(log2(x+0.1))})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#17.IPS   Error in if (x <= 0) { : missing value where TRUE/FALSE needed
source('/work/dy/pub_data/immune/Immunophenogram/IPS_dy.R') 
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  IPS = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #if(length(which(is.na(.x)))==0){
        #	.x=.x
        #}else if(length(which(is.na(.x)))>0){
        #	.x[which(is.na(.x))]=0
        #}			
        IPS=IPS_dy(data.frame(.x))
        score=IPS$IPS
        names(score)=IPS$SAMPLE
        score
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#18.TIDE  ????????????????????? 

#19.NLRP3  ssgsea  
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  NLRP3 = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[19]],',')[[1]];
        #rownames(xm) <- .x$symbol
        xm_ES <- GSVA::gsva(as.matrix(.x),list(NLRP3=geneList) ,method= "ssgsea" , verbose=TRUE)
        xm_ES[1,]
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)

#20.MPS 33上调基因加和-12下调基因加和
MPS_signature <- openxlsx::read.xlsx("/work/dy/pub_data/signature/ICB/MPS-suppl.xlsx",sheet='input')
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  MPS = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        weight=MPS_signature$Sign.in.the.signature
        names(weight)=MPS_signature$Gene.Symbol
        interGenes=intersect(rownames(.x),names(weight));
        #交集的genelist和系数
        geneListInter=weight[interGenes]
        apply(.x[match(names(geneListInter),rownames(.x)),],2,function(x){sum(as.numeric(x)*as.numeric(geneListInter),na.rm=T)})
        
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#21.CSIS  基因与系数乘积求和
CSIS_signature <- openxlsx::read.xlsx("/work/dy/pub_data/signature/ICB/Cancer-Specific Immune Prognostic Signature-supply.xlsx",sheet='input')
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  CSIS = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        weight=CSIS_signature$xishu
        names(weight)=CSIS_signature$gene
        interGenes=intersect(rownames(.x),names(weight));
        #交集的genelist和系数
        geneListInter=weight[interGenes]
        apply(.x[match(names(geneListInter),rownames(.x)),],2,function(x){sum(as.numeric(x)*as.numeric(geneListInter),na.rm=T)})
        
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#22.ImmuneCells 原本是监督分类 没有计算score  我们这里计算 summing the log2 Z scores of genes
ImSig.genes <- readRDS('/work/dy/pub_data/signature/ICB/Immune_cells_analysis-master/ImSig.rds')
#aa=readRDS('/work/dy/pub_data/signature/ICB/Immune_cells_analysis-master/NatMed_103samples_pData.rds')
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  ImSig = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {
        geneList=ImSig.genes
        interGenes=intersect(rownames(.x),geneList);
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){sum(log2(scale(as.numeric(x))),na.rm=T)})				
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#23.CGAs mean
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  CGAs = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[23]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){mean(as.numeric(x),na.rm=T)})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#24.IRGs  加权求和
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  IRGs = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        weight=c(0.32196,-0.64921,-0.32677,0.23573,0.39005,0.02975,0.39830,0.14607,-0.68625,0.38166,-0.03522)
        geneList=c('LEPR','PRLHR','NR2F2','PRL','NRP1','TNFRSF10B','TNFRSF10A','PLAU','IFI30','ANGPTL5','IGF1')
        names(weight)=geneList
        interGenes=intersect(rownames(.x),names(weight));
        #交集的genelist和系数
        geneListInter=weight[interGenes]
        apply(.x[match(names(geneListInter),rownames(.x)),],2,function(x){sum(as.numeric(x)*as.numeric(geneListInter),na.rm=T)})		
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")


#25. LRRC15.CAF PCA系数加权  系数？？？？？？？

#26.Inflammatory    sum(log2(scale(c(1,2,3,5,4,8,8,5))))
# summing the log2 Z scores of genes
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  Inflammatory = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList=strsplit(ICB_signature$genes[[26]],',')[[1]];
        interGenes=intersect(rownames(.x),geneList);
        #交集的genelist
        apply(.x[match(interGenes,rownames(.x)),],2,function(x){sum(log2(scale(as.numeric(x))),na.rm=T)})	
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#27.EMT adding the sum of the log2 Z scores of 6 established mesenchymal genes 
#(AGER, FN1, MMP2, SNAI2, VIM, ZEB2) and 
#subtracting the sum of the log2 Z scores of 6 established epithelial genes (CDH1, CDH3, CLDN4, EPCAM, MAL2, and ST14)
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  EMT = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        geneList1=c('AGER', 'FN1', 'MMP2', 'SNAI2', 'VIM', 'ZEB2')
        geneList2=c('CDH1', 'CDH3', 'CLDN4', 'EPCAM', 'MAL2', 'ST14')
        #交集的genelist
        apply(.x[match(geneList1,rownames(.x),nomatch=0),],2,function(x){sum(log2(scale(as.numeric(x))),na.rm=T)})-
          apply(.x[match(geneList2,rownames(.x),nomatch=0),],2,function(x){sum(log2(scale(as.numeric(x))),na.rm=T)})
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#28.Blood
#大于得分有效 24.009-(0.8697 x ADAM17) + (0.7486x CDK2)-(0.5885x CDKN2A) + (0.3462x DPP4)-(0.2401x ERBB2) + 
#(1.7427x HLA-DRA) + (0.2481x ICOS)-(1.1975x ITGA4)-(1.0184x LARGE) + (1.1721x MYC)-(0.6531x NAB2)-(1.1491x NRAS) + 
#(0.7377x RHOC)-(1.0585x TGFB1) + (0.8328x TIMP1)
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  Blood = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        weight=c(0.8697,0.7486,-0.5885,0.3462,-0.2401,1.7427,0.2481,-1.1975,-1.0184,1.1721,-0.6531,-1.1491,0.7377,-1.0585,0.8328)
        geneList=c('ADAM17','CDK2','CDKN2A','DPP4','ERBB2','HLA-DRA','ICOS','ITGA4','LARGE','MYC','NAB2','NRAS','RHOC','TGFB1','TIMP1')
        names(weight)=geneList
        interGenes=intersect(rownames(.x),names(weight));
        #交集的genelist和系数
        geneListInter=weight[interGenes]
        24.009-apply(.x[match(names(geneListInter),rownames(.x)),],2,function(x){sum(as.numeric(x)*as.numeric(geneListInter),na.rm=T)})		
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)
#saveRDS(mRNA_exp_score,"/work/dy/project/2020circRNA/result/train_test_quantiseq.rds.gz",compress = "gzip")

#29.Traf3
Traf3_signature <- openxlsx::read.xlsx("/work/dy/pub_data/signature/ICB/Traf3-KO.xlsx",sheet='input')
mRNA_exp <- mRNA_exp %>% dplyr::mutate(
  Traf3 = purrr::map(.x=mRNA_FPKM,function(.x){ 
    try( #报错 不停止运行程序
      {	
        #.x=mRNA_exp$mRNA_FPKM[[1]]
        weight=Traf3_signature$xishu
        names(weight)=Traf3_signature$genes
        interGenes=intersect(rownames(.x),names(weight));
        #交集的genelist和系数
        geneListInter=weight[interGenes]
        apply(.x[match(names(geneListInter),rownames(.x)),],2,function(x){sum(as.numeric(x)*as.numeric(geneListInter),na.rm=T)})		
      },silent = F)	
  })) #%>% dplyr::select(-mRNA_FPKM)



mRNA_exp$clinical_ImmAPA <- list(GSE91061_ImmAPASig ,TCGA_SKCM_metastasis_ImmAPASig)
readr::write_rds(mRNA_exp,"/work/gywang/project/Immune_APA/ImmAPA_ICB/data/GSE91061_metastasis_SKCM_24_ICB_related_signature.rds")
