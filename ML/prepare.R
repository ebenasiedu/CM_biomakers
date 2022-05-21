library(GEOquery)
library(dplyr)
library(limma)

workDir = getwd()

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*1000)

#Get raw expression
#options('download.file.method.GEOquery'='auto')
#options('download.file.method.GEOquery' = 'libcurl')
gset <- getGEO("GSE72058", destdir = workDir, GSEMatrix = T, 
			AnnotGPL = T, getGPL = T, parseCharacteristics = T)

expr <- as.data.frame (gset@assayData[["exprs"]])
expr <- as.data.frame(normalizeBetweenArrays(expr)) # normalize
#Get annotation
annotation <- gset@featureData@data
anno <- select(annotation, 'ID', 'Gene Symbol')
exprs <- mutate(expr, probe = anno$ID, Gene = anno$'Gene Symbol', .before = GSM2672129)
rownames(exprs) <- NULL
exprs$probe = NULL
#deduplicate
exprs = exprs[!duplicated(exprs[c("Gene")]),]

#Get phenodata
pheno <- gset@phenoData@data
phenodata <- select(pheno, 'geo_accession', 'sample class:ch1')
rownames(phenodata) <- NULL
colnames(phenodata) <- c("Sample", "Group")

# Save raw expression
write.csv(exprs, "GSE72058.expression.csv")
write.csv(phenodata, "GSE72058.phenodata.csv")

	
## Prepare for biomarker genes
df1 <- read.csv("biomarkers.csv")
df <- merge(df1, exprs, by.y = "Gene", all= F)
write.csv(df, "GSE72058.data.csv") 
