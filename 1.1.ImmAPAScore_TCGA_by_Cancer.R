library(magrittr)
library(ggplot2)
options(stringsAsFactors = F)
library(doParallel)
library(doMC)
registerDoMC()
library(foreach)
library(ppcor)

my.pcor.test <- function(...) {
  obj<- try(pcor.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(c(NA,NA)) else return(obj[c("estimate","p.value")])}

ImmPathways <- read.table("1793_genes_in_17_Immune_Pathway.txt") ### provided in data
Purity <- read.table("All_tumor_purity_32_cancer_type.txt") # provided in data
TCGA_APA <- readr::read_rds("pancan32_APA.DATA.rds.gz")   # TC3A: The Cancer 3' UTR Atlas (http://tc3a.org/)
TCGA_FPKM <- readr::read_rds("TCGA_33CanFPKM.rds.gz")   # TCGA data portal (http://gdac.broadinstitute.org/)
GeneAnn <- readr::read_rds("/home/gywang/project/20200831ye/GeneAnnotation_GeneTypes.rds")

TCGA_FPKMCoding <- TCGA_FPKM %>% dplyr::mutate(Expr=purrr::map(.x=Expr,function(.x){
  .x <- .x[.x$Gene %in% GeneAnn[GeneAnn$GeneTypes %in% "protein_coding","GeneName"],]
  return(.x)
}))

PurityAll <- dplyr::bind_rows(Purity$Purity) %>% dplyr::filter(substr(barcode,14,15) %in% c("01","03")) %>%
  dplyr::mutate(barcode=gsub("\\-","\\.",substr(barcode,1,12)))

TCGA_APA_mRNA <- dplyr::inner_join(TCGA_APA,TCGA_FPKMCoding,by=c("cancer_types"="cancertypes"))

foreach (i = unlist(TCGA_APA_mRNA$cancer_types)) %dopar% {
  .x <-  TCGA_APA_mRNA[TCGA_APA_mRNA$cancer_types %in% i,]$Expr[[1]]
  .x[,2:ncol(.x)] <- apply(.x[,2:ncol(.x)],2,function(x){log2(x+1)})
  .y <- TCGA_APA_mRNA[TCGA_APA_mRNA$cancer_types %in% i,]$APA[[1]]
  .y <- .y[,c("event_id","loci", "Proximal_Site",colnames(.y)[substr(colnames(.y),14,15) %in% c("01","03")])]
  colnames(.y) <- c("event_id","loci", "Proximal_Site",substr(colnames(.y)[4:ncol(.y)],1,12))
  CommonSamples <- intersect(grep("TCGA",colnames(.x),value=T),grep("TCGA",colnames(.y),value=T))
  CommonSamples <- intersect(CommonSamples,PurityAll$barcode)
  TumorPurity <- PurityAll[match(CommonSamples,PurityAll$barcode),]$Purity
  .x <- .x[,c("Gene",CommonSamples)]
  .y <- .y[,c("event_id","loci", "Proximal_Site",CommonSamples)]
  .x %>% tidyr::gather(key=Gene,value=expression,-Gene) %>% data.frame() -> .dx
  lapply(split(.x,.x$Gene), function(y){
    as.numeric(y[,CommonSamples]) -> dy
    Corr = data.frame(APAevents = .y$event_id,Genes =rep(unique(y$Gene),times=length(.y$event_id)))
    Corr[,c("estimate.rho","p.value")] <- t(apply(.y[,CommonSamples],1,function(x){
      xval <- as.numeric(x)
      unlist(my.pcor.test( xval[!is.na(xval)],dy[!is.na(xval)], TumorPurity[!is.na(xval)]))
    }))
    Corr$RankScore <- -log10(Corr$p.value) * sign(Corr$estimate.rho)
    return(Corr)
  }) %>% dplyr::bind_rows() -> AllCorr
  AllCorr[,3:5] <- apply(AllCorr[,3:5],2,function(x)signif(x,digits = 3))
  readr::write_rds(AllCorr,path=paste("/home/yye/Project/2020APA_Imm/ImmAPA/",i,"_Corr.rds.gz",sep=""),compress = "gz")
}

FileIns <- list.files(path="/home/yye/Project/2020APA_Imm/ImmAPA/",pattern = "*Corr.rds.gz")
RankScore  <- RankScore %>% dplyr::mutate(RankScore= lapply(FileIns,function(x){
  aa <- readr::read_rds(paste("/home/yye/Project/2020APA_Imm/ImmAPA/",x,sep=""))
  aa <- aa[!is.na(aa$estimate.rho) & !is.na(aa$p.value),]
  aa$RankScore <- -log10(aa$p.value+10^-20) * sign(aa$estimate.rho)
  APA_RES <- lapply(split(aa,aa$APAevents),function(sub){
    subRankList <- sub$RankScore
    names(subRankList) <- sub$Genes
    subRankList  <- rev(sort(subRankList))
    fgseaRes_sub <- clusterProfiler::GSEA(geneList = subRankList,
                                          nPerm=1000,TERM2GENE=ImmPathwaysD,
                                          minGSSize=5,
                                          maxGSSize=1000,
                                          pvalueCutoff = 1)
    fgseaRes_subD <- data.frame(fgseaRes_sub)[c(2:10)]
    fgseaRes_subD$APA_RES <- ifelse(fgseaRes_subD$enrichmentScore > 0,1-2*fgseaRes_subD$pvalue,2*fgseaRes_subD$pvalue-1)
    fgseaRes_subD$APA_events <- rep(unique(sub$APA_events),times=nrow( fgseaRes_sub))
    print(unique(sub$APA_events))
    return(fgseaRes_subD)
  }) %>% dplyr::bind_rows()
  return(APA_RES)
}))

readr::write_rds(RankScore,path = "/home/yye/Project/2020APA_Imm/ImmAPA/ImmAPAscore/Corr_ImmAPA_Score.rds.gz",compress = "gz")
