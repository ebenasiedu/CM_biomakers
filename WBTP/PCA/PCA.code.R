library(ggpubr)
library(readxl)
library(PCAtools)
library(patchwork)

## PCA analysis ##
df <- read_excel("CM.HC.GSE117613.xlsx")
mat <- df
# read phenodata
df2 <- read_excel("pheno.CM.HC.GSE117613.xlsx")
metadata <- df2
row.names(metadata) <- df2$ID

# check if names match
all(colnames(mat) == rownames(metadata))

#Run PCA
p <- pca(mat, metadata = metadata, removeVar = 0.1, transposed = F)
p[["metadata"]] <- as.data.frame(p[["metadata"]])

#scree plot
screeplot(p, axisLabSize = 18, titleLabSize = 22)

#Plot PCA before outlier removal
PCA_before <- biplot(p,  colby = 'Sample',  colLegendTitle = 'Sample', legendPosition = "right",
		  legendLabSize = 12,  legendTitleSize = 14, legendIconSize = 5)
		  
# Remove outlier samples and re-plot PCA 
PCA_after <- biplot(p,  colby = 'Sample',  colLegendTitle = 'Sample',  legendPosition = "right",
		  legendLabSize = 12, legendTitleSize = 14, legendIconSize = 5)   

## Get composite plot
PCA_before | PCA_after

## END
