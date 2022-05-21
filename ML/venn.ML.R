
library(ggvenn)
library(ggpubr)

df1 <- read.csv("CM.csv")
df2 <- read.csv("benchmark.csv")

genes <- list(CM_Good_Predictors = df1$Gene,  Benchmark_Bad_Predictors = df2$Gene)
ggvenn(genes, show_percentage = F)

# Merge them
df <- merge(df1, df2, by.y = "Gene", all= TRUE)
write.csv(df, "ML.biomarker.csv")
