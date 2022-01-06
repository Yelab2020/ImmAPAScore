library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(magrittr)

Imm_APA_Score <- readRDS("/work/gywang/project/Immune_APA/Data/TCGA_Data/Imm_APA_Score.rds.gz")
TCGA_APA <- readr::read_rds("/work/gywang/project/Immune_APA/Data/TCGA_Data/pancan32_APA.DATA.rds.gz")
Immune_17pathway_gene <- readRDS("/work/gywang/project/Immune_APA/Data/Immune_pathway_data/ImmunePathwayList.rds")
ImmPort_17_gene_pathways <- read.table("/work/gywang/project/Immune_APA/Data/Immune_pathway_data/ImmPort_17_gene_pathways_GeneList.txt",sep = "\t")

ImmGene <- unlist(Immune_17pathway_gene) 
Immune_17pathway_gene <- Immune_17pathway_gene

for (i in 1:length(Immune_17pathway_gene)){
 sub <- data.frame(GeneName=Immune_17pathway_gene[i],ImmPath=rep(names(Immune_17pathway_gene)[i],times=length(Immune_17pathway_gene[i])))
 colnames(sub)[1] <- "GeneName"
 ifelse(i==1,ImmPathGene <- sub,ImmPathGene <- rbind(ImmPathGene,sub))
}

#ImmPathGene$GeneName[duplicated(ImmPathGene$GeneName)] %>% unique() -> duplicate_gene

ImmScoreAllGene <- lapply(Imm_APA_Score$RankScore,function(x){
  #x <- Imm_APA_Score$cancer_types[[1]]
  x$GeneName <- as.data.frame(do.call(rbind,str_split(x$APAevents,"\\|")))[,2]
  x %>% group_by(GeneName) %>% summarise(mean(abs(lncRES)))
})


#merge(unlist(ImmScoreAllGene))
ImmGeneAPA <- filter(x,GeneName %in% ImmGene)
# plot each immoune pathway gene number
Imm_gene_count  <- lapply(Immune_17pathway_gene, function(x){
  pathway_gene_number  <- length(x)
})  %>% as.data.frame()

Imm_gene_count <- t(Imm_gene_count) %>% as.data.frame()
colnames(Imm_gene_count)="PathwayGeneCount"
Imm_gene_count$log2_ImmGeneCount <- log2(Imm_gene_count$PathwayGeneCount)
Imm_gene_count$Pathway <- str_replace_all(rownames(Imm_gene_count),"_"," ")


pdf("/home/gywang/project/APA_CancerImmune/Figure/Figure1/Number of gene in each Immune pathway.pdf",width =5.5,height = 7)
ggplot(Imm_gene_count,aes(x=log2_ImmGeneCount,y=Pathway))+ #
  geom_col(fill="#D2B48C")+                                              #
  scale_x_continuous(breaks = seq(2,8,2),expand = c(0,0))+                 #
  scale_y_discrete(labels = element_text(size = 14,margin=margin(r=30,l=30)))+ # 
  labs(x=NULL,y=NULL,title="Log2(Number of immune-related genes)")+       # 
  theme_bw()   +                                                    # 
  theme(plot.title = element_text(hjust = 0,size=12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.ticks.y = element_blank())
dev.off()


############### get all APA PDUI
all_pdui <- lapply(TCGA_APA$APA, function(x){
  x$PDUI  <- apply(x,1,function(y){
    y[-1:-3] %>% as.numeric() %>% na.omit() %>% mean() })
  sub <- dplyr::select(x,event_id,PDUI)
  return(sub)
}) 


for (i in 1:length(all_pdui)) {
if (i==1) {sub <- all_pdui[[1]]
}else{sub <- full_join(sub,all_pdui[[i]],by="event_id")}}

colnames(sub) <- c("APA.events",TCGA_APA$cancer_types)
AllCancerPDUI <- dplyr::select(sub,-15)
AllCancerPDUI$mean.pdui <- AllCancerPDUI[,-1] %>%  apply(1,function(x){
  x  %>% na.omit() %>% mean()
})

#### replace the NA by the mean.pdui
#AllCancerPDUI  <- AllCancerPDUI %>% apply(1,function(x){ x[is.na(x)] <- x[33]  return(x)}) %>% t() %>% as.data.frame()


AllCancerPDUI$Gene.Name <- data.frame(do.call(rbind,str_split(AllCancerPDUI$APA.events,"\\|")))[,2]

gene_pudi <- AllCancerPDUI[,c(-1,-34)] %>% apply(2,as.numeric) %>% set_rownames(AllCancerPDUI$APA.event)
Imm_gene_pudi <- AllCancerPDUI[AllCancerPDUI$Gene.Name %in% ImmGene,]
Imm_gene_pudi_path <- merge(ImmPathGene,Imm_gene_pudi,by.x ="GeneName" ,by.y = "Gene.Name")


split(Imm_gene_pudi_path,Imm_gene_pudi_path$ImmPath)%>% lapply(function(x){
  c(ImmPath=unique(x$ImmPath),maxPDUI=max(na.omit(x$mean.pdui)),minPUDI=min(na.omit(x$mean.pdui)),range=max(na.omit(x$mean.pdui))-min(na.omit(x$mean.pdui))  )
}) %>% bind_rows() -> ImmGene_path_PDUI

Imm_gene_pudi_matrix <- Imm_gene_pudi[,-1] %>% apply(2,as.numeric) %>% set_rownames(Imm_gene_pudi$APA.event)

tt <- merge(Imm_gene_pudi,ImmPathGene,by.x="Gene.Name",by.y="GeneName")
tt1 <- tt[order(tt$ImmPath,tt$mean.pdui),]
tt1 <- tt1 %>% dplyr::select(ImmPath,everything())
colnames(tt1)[1:3] <- c("immune_pathway","gene_name","transcript") 
openxlsx::write.xlsx(tt1,"/work/gywang/project/Immune_APA/supplementary_table/Gene_in_all_immune_pathway_PDUI_31_cacner_types.xlsx")

tt2 <- filter(tt1,!is.na(mean.pdui))
ha  <- ComplexHeatmap::rowAnnotation(type=tt1$ImmPath)
pdf("/home/gywang/project/APA_CancerImmune/Figure/Figure1/Immune gene PDUI in 31 cancer.pdf")
ComplexHeatmap::Heatmap(tt1[,3:34],show_row_names = F,cluster_rows = F,na_col = "white",cluster_columns =F,left_annotation = ha,col=colorRampPalette(c("blue","red"))(50))
dev.off()


# cluster na to -1

tt2 <- tt1
colnames(tt2)[34] <- "Mean PDUI" 
ImmGeneCount <- tt2[,3:34] %>% apply(2,function(x){
 sub <- cbind(tt2$APA.events,as.numeric(x)) %>% as.data.frame()
 sub$V2 <- as.numeric(x)
 sub<-  filter(sub,V2 > -1)
 y <- length(unique(sub$V1))
 names(y) <- names(x)
 y
})


tt2[,3:34]  <- tt2[,3:34] %>% apply(1,function(x){
  x[is.na(x)] <- -1
  as.numeric(x)
  x
}) %>% t()
tt2$ImmPath <- gsub("_"," ",tt2$ImmPath)
ha  <- ComplexHeatmap::rowAnnotation(type=tt2$ImmPath)
top_anno = ComplexHeatmap::HeatmapAnnotation(name = "Number of Immune related genes",foo = anno_barplot(ImmGeneCount),border = F)


pdf("/home/gywang/project/APA_CancerImmune/Figure/Figure1/Immune gene PDUI in 31 cancer(convert na to mean).pdf",width = 8,height = 12)
ComplexHeatmap::Heatmap(tt2[,3:34],show_row_names = F,cluster_rows = T,
                        column_title_side = "top",
                        na_col = "white",cluster_columns =T,
                        left_annotation = ha,
                        top_annotation = top_anno,
                        show_row_dend = F,
                        show_column_dend = F,
                        heatmap_legend_param = list(
                          title = "PDUI", at = c(1, 0), 
                          title_gp = gpar(col = "black"),
                          labels = c("1", "0")
                          #border = "red"
                        ),
                        col=colorRampPalette(c("gray","blue","red"))(50))
dev.off()

pdf("/work/gywang/project/Immune_APA/Figure/F1_Algorithms_DataInfo/Immune Gene PDUI in TCGA heatmap plot.pdf",width = 10,height = 12)
ComplexHeatmap::Heatmap(tt2[,3:34],show_row_names = F,cluster_rows = T,
                        col = colorRampPalette(colors = c("#e0e0e0","#e0e0e0","#e0e0e0","#e0e0e0","#7C4DFF","#03a9f4","#fff9c4","#e57373","#ff3f2f"))(1000),
                        #col=colorRampPalette(c("gray","blue","red"))(50)
                        column_title_side = "top",
                        na_col = "white",
                        cluster_columns =T,
                        left_annotation = ha,
                        top_annotation = top_anno,
                        show_row_dend = F,
                        show_column_dend = F,
                        heatmap_legend_param = list(
                          title = "PDUI", at = c(1, 0), 
                          title_gp = gpar(col = "black"),
                          labels = c("1", "0")
                          #border = "red"
                        )
                        )
dev.off()


