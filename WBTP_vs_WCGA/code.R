
library(ggvenn)
library(ggpubr)

df1 <- read.csv("WBTP.csv")
df2 <- read.csv("WCGNA.csv")

genes <- list(WBTP = df1$Gene,  WGCNA = df2$Gene)
ggvenn(genes, show_percentage = F)

# Merge them
df <- merge(df1, df2, by.y = "Gene", all= TRUE)
write.csv(df, "CM.biomarker.csv")				
