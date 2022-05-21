library(clusterProfiler)
library(DOSE)    
library(grid)        
library(ggplot2)
library(org.Hs.eg.db)


# reading in data
df = rea.csv("red.genes.csv")

# we want the log2 fold change 
original_gene_list <- df$LogFC

# name the vector
names(original_gene_list) <- df$Gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

##GO enrichment
GO <- gseGO( geneList = gene_list ,  ont = "MF",  OrgDb = org.Hs.eg.db,
		  keyType = "SYMBOL",   exponent = 1,  minGSSize = 10,  maxGSSize = 200,
		  eps = 1e-10,   pvalueCutoff = 0.05,   pAdjustMethod = "BH",   verbose = TRUE,
		  seed = FALSE,  by = "fgsea" )

##Visualization
require(DOSE)
CC <- dotplot(GO, x = "NES", color = "p.adjust", showCategory = 20,
		font.size = 14, orderBy = "x",  split=".sign") + facet_grid(.~.sign)  
		

