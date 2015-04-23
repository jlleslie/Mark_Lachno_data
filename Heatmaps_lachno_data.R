setwd("~/Desktop/Mark_Lachno_data")

library(reshape2)
library(ggplot2)
library(plyr)
library(scales)

data<-read.table(file="Mark_Epithelial_Norm_D4_group_change_untreated_D4_CIP.txt", header=T)
data$VPI_un<-NULL
data$VPI_un<-data$VPI_Mean/data$Untreated_Mean
data$D4_un<-NULL
data$D4_un<-data$D4_Mean/data$Untreated_Mean
data$CIP_un<-NULL
data$CIP_un<-data$CIP_Mean/data$Untreated_Mean
data$D4_CIP<-NULL
data$D4_CIP<-data$D4_Mean/data$CIP_Mean
data$CIP_D4<-NULL
data$CIP_D4<-data$CIP_Mean/data$D4_Mean
row.names(data)<-data$genes
datas<-data[,15:16]
datas$genes<-rownames(datas)
datas.m <- melt(datas)
datas.m <- ddply(datas.m, .(variable), transform, rescale = rescale(value))
p <- ggplot(datas.m, aes(variable, genes)) + 
  geom_tile(aes(fill = rescale),colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")

base_size <- 9
 p + theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
 theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))

data.1<-data[,14:18]
data.ma<-as.matrix(data.1)
library(gplots)
heatmap.2(data.ma, cexCol=.8, scale = "column", col=bluered, srtCol= 0,
          density.info="none", trace="none", key.xlab="Fold change",  xlab = "Comparisons",ylab ="epithelial genes",  main="Epithelial gene expression")




#######

cyto<-read.table(file="Mark_Cytokine_Norm_withRepeats_Untreated_D4_CIP.txt", header=T)
cyto$VPI_un<-NULL
cyto$VPI_un<-cyto$VPI_Mean/cyto$Untreated_Mean
cyto$D4_un<-NULL
cyto$D4_un<-cyto$D4_Mean/cyto$Untreated_Mean
cyto$CIP_un<-NULL
cyto$CIP_un<-cyto$CIP_Mean/cyto$Untreated_Mean
row.names(cyto)<-cyto$Genes
cytos<-cyto[,15:16]
cytos$genes<-rownames(cytos)
cytos.m <- melt(cytos)
cytos.m <- ddply(cytos.m, .(variable), transform, rescale = rescale(value))
p <- ggplot(cytos.m, aes(variable, genes)) + 
  geom_tile(aes(fill = rescale),colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue")

base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))

cyto.1<-cyto[,15:16]
cyto.ma<-as.matrix(cyto.1)
library(gplots)
heatmap.2(cyto.ma, cexCol=.8, scale = "column", col=bluered, srtCol= 0,
          density.info="none", trace="none", key.xlab="Fold change",  xlab = "Comparisons",ylab =" cytokine genes",  main="Cytokine expression across monocolonized mice")
##########################

allmice<-read.table(file="Mark_Epithelial_Norm_dCT_individmice.txt", header=T)
row.names(allmice)<-allmice$Genes
allmice<-allmice[,2:20]
allmice$mean_untreated<-rowMeans(allmice[,1:5])
allmice[14,20]<-mean(13.10, 13.80)
allmice<-allmice[6:20]

for (i in 1:ncol(allmice)){
  allmice[,i]<- allmice[,i]/allmice$mean_untreated
}
library(gplots)
allmice<-allmice[,1:14]
cipd4<-allmice[,1:10]
cipd4.ma<-as.matrix(cipd4)
heatmap.2(cipd4.ma, cexCol=.8, scale = "column", col=bluered, srtCol= 45,
          density.info="none", trace="none", key.xlab="Fold change",  xlab = "Comparisons",ylab ="epithelial genes",  main="Epithelial gene expression")


allmice.ma<-as.matrix(allmice)
heatmap.2(allmice.ma, cexCol=.8, scale = "column", col=bluered, srtCol= 45,
          density.info="none", trace="none", key.xlab="Fold change",  xlab = "Comparisons",ylab ="epithelial genes",  main="Epithelial gene expression")

