library(tidyverse)

Imm_APA_Score <- readr::read_rds("/home/gywang/project/APA_CancerImmune/Data/Imm_APA_Score.rds.gz")

for (i in 1:31){
  sub <- Imm_APA_Score$RankScoreSig[[i]]$APAevents %>% as.data.frame()
  ifelse(i==1,All.sig.apa <- sub, All.sig.apa <- rbind(All.sig.apa,sub) %>% unique())
}
colnames(All.sig.apa) <- "APA.events"

for (i in 1:31) {
  apa.pair <-   Imm_APA_Score$RankScoreSig[[i]]$APAevents %>% table()  %>% as.data.frame()
  apa.cancertype <-  Imm_APA_Score$RankScoreSig[[i]]$APAevents %>% unique() %>% table()  %>% as.data.frame() 
  apa.count <- merge(apa.pair,apa.cancertype,by = ".")
  colnames(apa.count) <- c("APA.events","apa.pair","apa.cancertype")
  apa.count <- left_join(All.sig.apa,apa.count,by=c("APA.events"))
  apa.count[is.na(apa.count)] <- 0
  ifelse(i==1, 
         APA_counts <- apa.count,
         APA_counts[,c("apa.pair","apa.cancertype")] <- c(APA_counts$apa.pair + apa.count$apa.pair , APA_counts$apa.cancertype + apa.count$apa.cancertype)
  ) }
APA_counts[,c("transcript_ID","GeneName","chr","strand")] <- data.frame(do.call(rbind,strsplit(as.character(APA_counts$APA.events),"\\|")),stringsAsFactors=FALSE)

for (i in 1:22) { 
  for (b in c(i:(i+10))) {
    PanAPA_cutoff <- paste("CanTypes",i,"Path",b,sep = "") 
    APA_counts[ ,c(PanAPA_cutoff)] <- ifelse(APA_counts$apa.cancertype >=i & APA_counts$apa.pair >=b,"Sig","nonSig")
  }}

Imm_APA_Score$APAevents_count <- list(APA_counts)
readr::write_rds(Imm_APA_Score,"/home/gywang/project/APA_CancerImmune/Data/Imm_APA_Score.rds.gz",compress = "gz")
