library(tidyverse)
library(corrplot)
library(epitools)
library(reshape2)

options(stringsAsFactors = F)

Imm_APA_Score <- readRDS("/work/gywang/project/Immune_APA/Data/TCGA_Data/Imm_APA_Score.rds.gz")
TCGA_APA_Timer <- readr::read_rds("/work/gywang/project/Immune_APA/Data/TCGA_Data/TCGA_ICI/TCGA_APA_Timer.rds.gz")
APA_counts <- Imm_APA_Score$APAevents_count[[1]]


col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))



TCGA_APA_Timer
c=10 # number of cancer types APA events appeared
p=15 # number of apa-immune pair
Cutoff=paste("CanTypes",c,"Path",p,sep = "")


for (i in 1:31) {
  x <- TCGA_APA_Timer$ICI0.2_ImmAPA_cutoff[[i]]
  II_APA <- x[x[,Cutoff] %in% "II",] %>% group_by(ImmuneCell) %>% summarise(counts=n())
  colnames(II_APA)[2] <- "II.APA"
  ICI_APA <- x %>% group_by(ImmuneCell,) %>% summarise(counts=n())
  colnames(ICI_APA)[2] <- "ICI.APA"
  Imm_APA <- Imm_APA_Score$RankScoreSig[[i]][Imm_APA_Score$RankScoreSig[[i]]$APAevents %in% APA_counts[APA_counts[,Cutoff] %in% "Sig",]$APA.events,]$APAevents %>% unique() %>% length()
  Total.APA <- TCGA_APA_Timer$APA[[1]]$event_id  %>% unique() %>% length()
  sub <- data.frame(II_APA,ICI_APA,Imm_APA,Total.APA)
  sub$NI.APA <- sub$ICI.APA - sub$II.APA
  sub$IN.APA <- sub$Imm_APA - sub$II.APA
  sub$NN.APA <- sub$Total.APA - sub$NI.APA - sub$Imm_APA
  sub$II_Imm.Proportion <- sub$II.APA/sub$Imm_APA
  sub <- sub[,-3]
  oddsratio <- apply(sub,1,function(x){
   #sum <- c(x$NN.APA,x$IN.APA,x$NI.APA,x$II.APA)
   sum <- x[c(8,7,6,2)] %>% as.numeric()
   or <- oddsratio(sum)
   ICI.APA.OR <- c(or$measure[-1,],or$p.value[2,2],Imm_APA_Score$cancer_types[[i]]) #%>% as.data.frame(stringsAsFactor=F)
   #rownames(ICI.APA.OR) <- c(colnames(or$measure)[1:3],"p.value","CancerType")
   #colnames(ICI.APA.OR) <- x[1]
   return(ICI.APA.OR)
  }) %>% as.data.frame()
  colnames(oddsratio) <- sub$ImmuneCell
  rownames(oddsratio)[c(4:5)] <- c("Pvalue","CancerType")
  proportion <- sub$II_Imm.Proportion  %>% as.data.frame() %>% apply(1,as.factor) %>% as.data.frame() %>% t()
  colnames(proportion) <- colnames(oddsratio)
  ICI.data <- rbind(oddsratio,proportion)
  rownames(ICI.data)[6] <- "proportion"
  ifelse(i==1,Total.ICI <- ICI.data,Total.ICI <- cbind(Total.ICI,ICI.data))
}
  
Total.ICI_n <- Total.ICI %>% t() %>% as.data.frame()
Total.ICI_n$cell.type <- rownames(Total.ICI_n)[1:6]
ICI_Imm_proportion <- Total.ICI_n[,c(5,6,7)]
ICI.proportion <- ICI_Imm_proportion %>% tidyr::spread(cell.type,proportion)
rownames(ICI.proportion) <- ICI.proportion$CancerType
ICI.proportion <- ICI.proportion[,-1]
ICI.proportion_n <-  apply(ICI.proportion,2, function(x)as.numeric(as.character(x)))# as.matrix()
rownames(ICI.proportion_n) <- rownames(ICI.proportion)
rownames(ICI.proportion_n) <- rownames(ICI.proportion_n) %>% str_sub(end = -8)

ICI.proportion_n <- dplyr::select(ICI.proportion_n,cancer_type,everything())
openxlsx::write.xlsx(ICI.proportion_n,"/work/gywang/project/Immune_APA/supplementary_table/fig2a_proportion.xlsx")
pdf(paste("/home/gywang/project/APA_CancerImmune/Figure/Figure2_ICI/TIMER ICI ImmAPA Proportion (",Cutoff,")_N.pdf"),width = 10,height = 20)
corrplot(ICI.proportion_n,
               method = "pie",
               tl.cex= 0.7,cl.pos="b",
               cl.lim=c(0,0.7),
               cl.ratio=0.07,
               cl.length=6,
               col=col2(50),
               tl.srt = 90, is.corr = FALSE,
               title=paste("Cutoff (",Cutoff,")",sep=""),
         mar=c(0, 0, 1, 0))
dev.off()  


plist <- list()
for (i in 1:6) {
  CellType <- Total.ICI_n$cell.type[1:6][i]
  oddsratio_cell <- dplyr::filter(Total.ICI_n,cell.type==CellType)
  data <- oddsratio_cell
  data[,c(1:4)] <- oddsratio_cell[,c(1:4)] %>% apply(2, function(x)as.numeric(as.character(x)))# 
  data$CancerType <- gsub("_ImmAPA","",data$CancerType)
  data$class <- ifelse(data$estimate > 1 & data$Pvalue < 0.05,"risk","nonrisk")
  risk.count <- data %>% group_by(class) %>% summarise(counts=n())
  colnames(risk.count)[2] <- CellType
  data$significance <- ifelse(data$Pvalue < 0.05,"P < 0.05","P >= 0.05")
p <- ggplot(data,aes(x=estimate,y=CancerType))+
    geom_point(aes(x=estimate,y=CancerType,color=significance),size=2)+  
    #scale_color_manual(values=c("red","orange"))+
    geom_errorbar(aes(xmin=lower, xmax=upper), width=.2,size=0.5)+                # 在点上加errorbar        
    geom_vline(aes(xintercept=1), colour="black", linetype="dashed",size=0.3)+    # vline 加垂直线
    # scale_x_continuous("ICI ImmAPA Odds Ratios and 95% CIs")+
    # scale_fill_discrete(title="")+
    #labs(title=colnames(data),x="CancerType",y="CancerType")+
    xlab(paste(CellType,"Odds Ratio"))+
    ylab(NULL)+
    ggtitle("")+
   labs(color="P value")+
    theme(panel.background=element_rect(fill="white"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.title.x=element_text(color="black",size=12),
            axis.text.x =element_text(color="black",size=10),
            axis.text.y = element_text(color="black",size=10,face = "bold"),
            # axis.title.y = element_blank(),
            # axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.y = element_line(colour = "black"),
            axis.line.x = element_line(colour = "black"),
            #strip.background = element_blank(),
            legend.position = c(0.9,0.4),
            legend.justification = c(0.9,0.4)
            )
plist[[i]] <- p 
ifelse(i==1,risk.counts <- risk.count,risk.counts <- merge(risk.counts,risk.count,by="class") )
}

pdf(paste("/home/gywang/project/APA_CancerImmune/Figure/Figure2_ICI/OddsRatio under",Cutoff,"by celltype.pdf"),width = 15,height = 9)
cowplot::plot_grid(plotlist =plist,ncol=3)
dev.off()

OR_Risk_per <- t(risk.counts) %>% as.data.frame()
OR_Risk_per$celltype <- rownames(OR_Risk_per)
colnames(OR_Risk_per) <- c("nonerisk","risk","celltype")
OR_Risk_per <- OR_Risk_per[-1,]
OR_Risk_per[,c(1,2)] <- apply(OR_Risk_per[,c(1,2)],2,as.numeric)
#OR_Risk_per <- t(OR_Risk_per) %>% as.data.frame()
t1 <- melt(OR_Risk_per,id.vars = "celltype")
t2 <- t1
t2 <- t1[c(7:12,1:6),]

pdf(paste("/home/gywang/project/APA_CancerImmune/Figure/Figure2_ICI/OR proportion",Cutoff,"by celltype.pdf"))
ggplot(t2,aes(x=celltype,y=value,fill=factor(variable,levels = c('nonerisk','risk'))))+
  geom_bar(stat = "identity",width=0.6)+
  labs(x = '') + labs(fill="")+
  theme(  panel.background=element_rect(colour=NA,fill="white"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title.x=element_text(color="black",size=12),
          axis.text.x =element_text(color="black",size=10),
          axis.text.y = element_text(color="black",size=10,face = "bold"),
          # axis.title.y = element_blank(),
          # axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),legend.position = "right",
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          strip.background = element_blank())
dev.off()



  
