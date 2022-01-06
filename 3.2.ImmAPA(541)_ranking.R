


Imm_APA_Score <- readRDS("/work/gywang/project/Immune_APA/Data/TCGA_Data/Imm_APA_Score.rds.gz")
TCGA_12_APA_DIFF <- readr::read_rds("/work/gywang/project/Immune_APA/Data/N_T_12cancer_APA/TCGA_12types_APA_T_N_Diff.rds.gz")
ImmAPA_541_ES <- readr::read_rds("/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_all_cluster.rds.gz")

ImmPathways <- readr::read_rds("/work/gywang/project/Immune_APA/2020APA_Imm_by_yye/ImmunePathwayList.rds")
ImmPathwaysD <- data.frame(term=rep(names(ImmPathways),times=unlist(lapply(ImmPathways,length))),gene=unlist(ImmPathways))

ALL.sig.APA <- lapply(Imm_APA_Score$RankScoreSig, function(x){
  dplyr::select(x,Description,NES,APAevents)
}) %>% dplyr::bind_rows()
ImmAPA_541_rank_score <-   ALL.sig.APA %>% dplyr::filter(APAevents %in% ImmAPA_541_ES$APA_events)


ALL.APA.name <- unique(ImmAPA_541_rank_score$APAevents)  %>% as.data.frame()
colnames(ALL.APA.name) <- "APAevents"

APA_PathSum <- lapply(split(ImmAPA_541_rank_score,ImmAPA_541_rank_score$Description), function(x){
  #x <- filter(ImmAPA_541_rank_score , Description=="TNF_Family_Members")
  sub <- group_by(x,APAevents) %>% summarise(sum(NES),count=n())
  colnames(sub)[2:3] <- c("NES.sum","CancerNumber")
  sub <- arrange(sub,desc(CancerNumber),desc(abs(NES.sum)))
  sub$Rank <- (length(sub$APAevents)+1 - as.numeric(rownames(sub)))/length(sub$APAevents)
  colnames(sub)[4] <- unique(x$Description)
  sub <- left_join(ALL.APA.name,sub,by="APAevents")
  rownames(sub) <- sub$APAevents
  sub[,4] <- ifelse(is.na(sub[,4]),0,as.numeric(as.character(sub[,4])))
  sub1 <- dplyr::select(sub,4)
  return(sub1)
}) %>% dplyr::bind_cols() #%>% set_rownames(rownames(a)) %>% t()
APA_PathSum$Interferons_Receptors <- 0
APA_PathSum <- dplyr::select(APA_PathSum,1:7,17,everything())


TCGA_APA_diff <- TCGA_12_APA_DIFF$Sig_Diff_rank[[1]]

APA_PathSum$APA.events <- rownames(APA_PathSum)
APA_PathSum <- left_join(APA_PathSum,TCGA_APA_diff)
APA_PathSum$APA_Diff_Rank <- ifelse(is.na(APA_PathSum$DiffRANK),0,APA_PathSum$DiffRANK)
rownames(APA_PathSum) <- APA_PathSum$APA.events
APA_PathSum <- dplyr::select(APA_PathSum,-18:-22)

APA_PathSum$Rank <- apply(APA_PathSum,1,function(x){
  sum(x[1:17] , 17*x[18])/36
}) %>% unlist()
APA_PathSum  <- arrange(APA_PathSum,desc(Rank))

APA_PathSum <- APA_PathSum %>% dplyr::select( -Interferons_Receptors)
readr::write_rds(APA_PathSum,"/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_ranking_n.rds.gz")

APA_PathSum <- readr::read_rds("/work/gywang/project/Immune_APA/Data/ImmAPA/ImmAPA_541_ES_ranking_n.rds.gz")
APA_PathSum$APA_events <- rownames(APA_PathSum)
APA_PathSum <- dplyr::select(APA_PathSum,APA_events,everything())
openxlsx::write.xlsx(APA_PathSum,"/work/gywang/project/Immune_APA/supplementary_table/541_ImmAPA_ranking.xlsx")


ImmAPA_541_ImmRanking <- APA_PathSum
top20_ImmAPA <- ImmAPA_541_ImmRanking %>% top_n(20,Rank) %>% arrange(desc(Rank))
rownames(top20_ImmAPA) <- data.frame(do.call(rbind,str_split(rownames(top20_ImmAPA),"\\|")))[,2]
colnames(top20_ImmAPA) <- gsub("_"," ",colnames(top20_ImmAPA))
pheatmap::pheatmap(top20_ImmAPA,scale = "none",cluster_col = FALSE,
                   cluster_rows = F,display_numbers = T,number_format = "%.2f",
                   cellwidth=30, cellheight = 12,
                   angle_col = 270,
                   main = "Ranked by regulation to Immune Pathway",
                   filename = "/work/gywang/project/Immune_APA/Figure/ImmAPA_541_ranking/Immune APA Enrichment Score Top20 heatmap.pdf",
                   width = 15,height = 8)
pheatmap::pheatmap(ImmAPA_541_ImmRanking[,1:18],scale = "none",cluster_col = FALSE,
                   cluster_rows = F,display_numbers = F,number_format = "%.1e",
                   cellwidth=30,#  cellheight = ,
                   show_rownames= F,
                   show_colnames = T,
                   fontsize=17,
                   angle_col = 270,
                   main = "Ranked by regulation to Immune Pathway",
                   filename = "/work/gywang/project/Immune_APA/Figure/ImmAPA_541_ranking/Immune APA Enrichment Score all 541 ImmAPA heatmap.pdf",
                   width = 11,height = 15)
