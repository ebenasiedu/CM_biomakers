library(WGCNA)

#### Data Preparations ####
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#Read in expression data set
Data = read.csv("GSE117613.expression.csv");
Data = Data[!duplicated(Data[c("Gene")]),]

# Format the data set:
datExpr = as.data.frame(t(Data[, c(2:44)]));
names(datExpr) = Data$Gene;

# Cluster samples to look for Remove outliers
sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, sub="", xlab="", cex.lab = 1.5, 
      cex.axis = 1.5, cex.main = 2) 
 abline(h = 78, col = "red");

## Removed outliers GSM3305179, GSM3305184, and GSM3305191


# Load and format clinical data
datTraits = read.csv("ClinicalTraits.csv");
rownames(datTraits) <- datTraits$Sample
datTraits$Sample = NULL

collectGarbage();

# Visualize how expression data relate to clinical data
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
#Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                   groupLabels = names(datTraits), 
                  main = "Sample dendrogram and trait heatmap")

# save prepared formates for next steps
save(datExpr, datTraits, file = "Data-Preparation.RData")



## Network construction and module detection ##
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Data-Preparation.RData");


# Soft thresholding for network topology
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#  Block-wise network construction and module detection
bwnet = blockwiseModules(datExpr, maxBlockSize = 11000,
                       power = 12, TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "DataTOM-blockwise",
                       verbose = 3)


# Observe modules and dynamic tree-cuts
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwnet$colors)

# open a graphics window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Save for further analysis
save(datExpr, datTraits, bwnet, sft, file = "Network-module-detection.RData")              



####  Module-trait correlationn hub gene identification####
##################

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
#lnames = load(file = "Network-module-detection.RData");
#The variable lnames contains the names of loaded variables.
#lnames
# Load network data saved in the second part.
lnames = load(file = "Network-module-detection.RData");

# Module-trait association
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
moduleColors <- labels2colors(bwnet$colors)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Visualize correlogram
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Gene relationship to trait
# Define variable weight containing the weight column of datTrait
Plasma.PfHRP2 = as.data.frame(datTraits$'Plasma.PfHRP2');
names(Plasma_PfHRP2) = "Plasma PfHRP2"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Plasma.PfHRP2, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(Plasma.PfHRP2), sep="");
names(GSPvalue) = paste("p.GS.", names(Plasma.PfHRP2), sep="");


# Intramodular analysis: finding hub genes
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;

#sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Plasma PfHRP2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# Save the module genes you prefer
genes <- as.data.frame(names(datExpr)[moduleColors=="red"])
colnames(genes) <- c("Gene")
write.csv(genes, "red-gene.csv")

#### END ####
