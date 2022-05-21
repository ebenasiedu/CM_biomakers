#This code generates plots for the DEG analysis

library(ggplot2)
library(rstatix)
library(tidyverse)
library(dplyr)
library(reshape2)    
library(gtools)   
library(ComplexHeatmap)
library(ggvenn)
library(patchwork)
library(ggsci)
library(ggrepel)
library(readxl)
library(ggvenn)
library(ggpubr)


# Read data
workDir <- getwd()
files <- list.files(path=workDir, pattern = "*.xlsx", full.names=TRUE)
filelist <- lapply(files, read_excel)
names(filelist) <- gsub(".xlsx","", list.files(workDir, pattern = "*.xlsx", full.names = FALSE), fixed = TRUE)
invisible(lapply(names(filelist), function(x) assign(x,filelist[[x]],envir=.GlobalEnv)))   # unlist from list to d/n dfs

#Deduplicate
CM.HC.GSE1124.GPL96 = CM.HC.GEO2R.GSE1124.GPL96[!duplicated(CM.HC.GEO2R.GSE1124.GPL96[c("Gene")]),]
SM.HC.GSE1124.GPL97 = SM.HC.GEO2R.GSE1124.GPL97[!duplicated(SM.HC.GEO2R.GSE1124.GPL97[c("Gene")]),]
CM.HC.GSE1124.GPL97 = CM.HC.GEO2R.GSE1124.GPL97[!duplicated(CM.HC.GEO2R.GSE1124.GPL97[c("Gene")]),]
SM.HC.GSE1124.GPL96 = SM.HC.GEO2R.GSE1124.GPL96[!duplicated(SM.HC.GEO2R.GSE1124.GPL96[c("Gene")]),]
CM.HC.GSE117613 = CM.HC.GEO2R.GSE117613[!duplicated(CM.HC.GEO2R.GSE117613[c("Gene")]),]
SM.HC.GSE117613 = SM.HC.GEO2R.GSE117613[!duplicated(SM.HC.GEO2R.GSE117613[c("Gene")]),]

#Volcano plots
#define differential expression  criteria
df <-CM.HC.GSE117613
df$Change <- "No change"
df$Change[df$logFC > 1.5 & df$adj.P.Val < 0.05] <- "Increased"
df$Change[df$logFC < -1.5 &  df$adj.P.Val  < 0.05] <- "Decreased"
df$delabel <- NA
df$delabel[df$Change != "No change"] <- df$Gene[df$Change != "No change"]
CMHC.GSE117613 <- ggplot(df, aes(x=logFC, y=-log10(adj.P.Val), col=Change, label=delabel)) + geom_point() +  theme_pubr(legend = "right") +  geom_text_repel() + 
                                          xlab(bquote('Fold Change'(Log[2]))) + ylab(bquote(Significance(-Log[10])))  + scale_color_aaas() 

#Plot volcano plots
volcano <- (CMHC.GSE1124.GPL96 | SMHC.GSE1124.GPL96) / (CMHC.GSE1124.GPL97 | SMHC.GSE1124.GPL97) /  (CMHC.GSE117613 | SMHC.GSE117613)
volcano +  plot_annotation(tag_levels = list(c('CM_HC.GPL96', 'SM_HC.GPL96', 'CM_HC.GPL97', 'SM_HC.GPL97', 'CM_HC.GSE117613', 'SM_HC.GSE117613')))

##Getting DEGs
df <- SM.HC.GSE1124.GPL96
dfOrdered <- df[order(df$adj.P.Val),]
dfSig1 <- subset(dfOrdered, adj.P.Val < 0.05)
dfSig <- subset(dfSig1, logFC > 1.5 | logFC < -1.5)
dfSig$adj.P.Val = NULL
write_csv(dfSig, "SM.HC.GSE1124.GPL96.DEGs.csv")

##  Getting CM-specific genes
workDir <- getwd()
files <- list.files(path=workDir, pattern = "*.DEGs.csv", full.names=TRUE)
filelist <- lapply(files, read_csv)
names(filelist) <- gsub(".DEGs.csv","", list.files(workDir, pattern = "*.DEGs.csv", full.names = FALSE), fixed = TRUE)
invisible(lapply(names(filelist), function(x) assign(x,filelist[[x]],envir=.GlobalEnv)))   # unlist from list to d/n dfs

df1 <- CM.HC.GSE117613
colnames(df1) = c("Gene", "GSE117613.CM")
df2 <- SM.HC.GSE117613
colnames(df2) = c("Gene", "GSE117613.SM")
df <- merge(df1, df2, by.y = "Gene", all= TRUE)
write_csv(df, "GSE117613.merged.csv")

##plot venn diagram
genes <- list(CM = df1$Gene, SM = df2$Gene)
GSE117613.venn <- ggvenn(genes, show_percentage = F)

## Heatmap of DEGs
df1 <- read_csv("GSE1124.merged.csv")
df2 <-  read_csv("GSE117613.merged.csv")
data <- merge(df1, df2, by.y = "Gene", all= TRUE)
df <- as.matrix(data[,2:5])
rownames(df) <- data$'Gene'
df[is.na(df)] = 0
df1 <- t(df)
heatmap <- pheatmap(df1,  cluster_rows = TRUE, cluster_cols = TRUE, 
		clustering_distance_rows = "euclidean", clustering_distance_col = "euclidean", legend = TRUE, 
		name = "logFC" , show_rownames = T,  show_colnames = F, fontsize = 10,  fontsize_row = 14,  fontsize_col = 10)


## Getting CM-genes
workDir <- getwd()
files <- list.files(path=workDir, pattern = "*.CM.csv", full.names=TRUE)
filelist <- lapply(files, read_csv)
names(filelist) <- gsub(".CM.csv","", list.files(workDir, pattern = "*.CM.csv", full.names = FALSE), fixed = TRUE)
invisible(lapply(names(filelist), function(x) assign(x,filelist[[x]],envir=.GlobalEnv)))   # unlist from list to d/n dfs

df1 <- GSE1124
df2 <- GSE117613
df <- merge(df1, df2, by.y = "Gene", all= TRUE)
#write_csv(df, "CM.merged.csv")

##plot venn diagram
genes <- list(GSE1124.CM = df1$Gene, GSE117613.CM = df2$Gene)
CM.venn2 <- ggvenn(genes, show_percentage = F)

# Plots for CM-genes
#Bar chart 
data <- read_csv("CM.genes.csv")
df <- melt(data, id = c("Gene"))
colnames(df) <- c("Gene", "Dataset", "logFC")
bar <- ggbarplot( df,  x="Gene",   y="logFC", fill = "Dataset", width = 0.2, position = position_dodge(), 
   		  xlab = "CM-gene",   ylab = "logFC",  rotate = T,  ggtheme = theme_classic()) + 
 		  theme(text = element_text(size = 14), legend.position="none") 
facet(bar, facet.by = "Dataset", scale = "free_x")

# heatmap
data <- read_csv("CM.genes.csv")
df <- as.matrix(data[,2:3])
rownames(df) <- data$'Gene'
pheatmap(df,  cluster_rows = TRUE, cluster_cols = F,  angle_col = "0",
		clustering_distance_rows = "euclidean", clustering_distance_col = "euclidean", legend = TRUE, 
		name = "logFC" , show_rownames = T,  show_colnames = T, fontsize = 13,  fontsize_row = 13,  fontsize_col = 13)

