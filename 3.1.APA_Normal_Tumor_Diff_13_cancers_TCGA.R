library(tidyverse)
library(ggplot2)



TCGA_APA_DIF_xlsx <- list.files("/home/gywang/project/APA_CancerImmune/Data/TCGA_APA_DIF_new_xlsx/",pattern = ".xlsx")
TCGA_12_APA_DIFF <- tibble::tibble(cancertype=gsub("_mean&pvalue.xlsx","",TCGA_APA_DIF_xlsx))
TCGA_12_APA_DIFF  <- TCGA_12_APA_DIFF %>% dplyr::mutate(APA_DIF = lapply(TCGA_APA_DIF_xlsx,function(x){
  openxlsx::read.xlsx (paste("/home/gywang/project/APA_CancerImmune/Data/TCGA_APA_DIF_new_xlsx/" , x,sep="")) #  %>% as.matrix()
}))


TCGA_12_APA_DIFF <- readr::read_rds("/home/gywang/project/APA_CancerImmune/Data/TCGA_12types_APA_T_N_Diff.rds.gz")
TCGA_12_APA_DIFF <- TCGA_12_APA_DIFF %>% dplyr::mutate(SigDiff=lapply(TCGA_12_APA_DIFF$APA_DIF,function(x){
  dplyr::filter(x, pvalue<0.05 ) %>% dplyr::select(gene,MeanDiff,pvalue)
}))



Sig_DiffAPA <-lapply(TCGA_12_APA_DIFF$APA_DIF,function(x){
  dplyr::filter(x, pvalue<0.05 ) %>% dplyr::select(gene,MeanDiff,pvalue)
})  %>% dplyr::bind_rows()

Sig_DiffAPA_n <-  table(Sig_DiffAPA$gene) %>% t() %>% as.data.frame() %>% select(2,3) %>% magrittr::set_colnames(c("APA.events","CanNumber"))
Sig_DiffAPA_n2 <- lapply(split(Sig_DiffAPA,Sig_DiffAPA$gene), function(x1){
  c(as.numeric(mean(x1$MeanDiff)),abs(as.numeric(mean(x1$MeanDiff))))
}) %>% dplyr::bind_rows() %>% t() %>% as.data.frame(StringsAsFactor=F)
Sig_DiffAPA_n2$APA.events <- rownames(Sig_DiffAPA_n2)
SigDiff_APA_12Cancer <- merge(Sig_DiffAPA_n,Sig_DiffAPA_n2,by.x="APA.events",by.y = )  %>% set_colnames(c("APA.events","Cancer.Number","MeanDiff","ABS.Diff"))
SigDiff_APA_12Cancer <- dplyr::arrange(SigDiff_APA_12Cancer,desc(Cancer.Number),desc(ABS.Diff))
SigDiff_APA_12Cancer$DiffRANK <- as.numeric(rev(rownames(SigDiff_APA_12Cancer)))/8347

TCGA_12_APA_DIFF$Sig_Diff_rank <- list(SigDiff_APA_12Cancer)

readr::write_rds(TCGA_12_APA_DIFF,"/home/gywang/project/APA_CancerImmune/Data/TCGA_12types_APA_T_N_Diff.rds.gz",compress = "gz")


