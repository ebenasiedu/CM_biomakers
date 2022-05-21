library(GEOquery)

## Reset VROOM size
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*1000)

# For GSE117613 expression data
gset <- getGEO("GSE117613", GSEMatrix =TRUE, AnnotGPL=TRUE)
expr <- as.data.frame (gset[["GSE117613_series_matrix.txt.gz"]]@assayData[["exprs"]])
# Save raw expression
write.csv(expr, "GSE117613_expression.csv")


# For GSE1124
# GPL96
gset <- getGEO("GSE1124", GSEMatrix =TRUE, AnnotGPL=TRUE)
#Get raw expression 
# GPL96-platform
expr <- as.data.frame (gset[["GSE1124-GPL96_series_matrix.txt.gz"]]@assayData[["exprs"]])
# Save raw expression
write.csv(expr, "GSE1124_GPL96_expression.csv")

# GPL97
gset <- getGEO("GSE1124", GSEMatrix =TRUE, AnnotGPL=TRUE)
#Get raw expression 
# GPL96-platform
expr <- as.data.frame (gset[["GSE1124-GPL97_series_matrix.txt.gz"]]@assayData[["exprs"]])
# Save raw expression
write.csv(expr, "GSE1124_GPL97_expression.csv")

## END
