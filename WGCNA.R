library("readODS")
library("xml2")
library("BiocManager")
library("dynamicTreeCut")
library("fastcluster")
library("memoise")
library("pkgconfig")
library("backports")
library("WGCNA")

############################################# START PART 1
gff <- read.table("/WGCNABarley/data/Hv_Morex.pgsb.Jul2020.gff3",header=F,skip=9,sep="\t")
annotation <- gff[gff[,3]=="mRNA",]
gff <- gff[gff[,3]=="gene",]
dim(gff)

geneexpressiondata <- read.table("/WGCNABarley/data/leaf_reads_per_gene_total_final_deseq_normed.tsv",header=T)
rownames(geneexpressiondata) <- geneexpressiondata[,1]
geneexpressiondata <- geneexpressiondata[,-1]
geneexpressiondata <- t(geneexpressiondata)

geneexpressiondata <- geneexpressiondata[order(rownames(geneexpressiondata)),]

gsg = goodSamplesGenes(geneexpressiondata, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(geneexpressiondata)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(geneexpressiondata)[!gsg$goodSamples], collapse = ", ")));
  geneexpressiondata = geneexpressiondata[gsg$goodSamples, gsg$goodGenes]
}

dim(geneexpressiondata)

#####


# clustering samples/genotypes, see if there are outliers 
sampleTree = hclust(dist(geneexpressiondata), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)



traitData <- read.csv("/WGCNABarley/data/AEM__FT_PH.csv", sep=",",header=T)
traitDataparents <- traitData[is.element(traitData[,1],rownames(geneexpressiondata)),]
rownames(traitDataparents) <- traitDataparents[,1]
traitDataparents <- traitDataparents[,-1]
colnames(traitDataparents)[2] <- "FT"

# check order
identical(rownames(traitDataparents),rownames(geneexpressiondata))


#visualize how the phenotypic traits relate to the sample dendrogram.
# Re-cluster samples
sampleTree2 = hclust(dist(geneexpressiondata), method = "average")
# Convert traits to a color representation:                 white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitDataparents, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.

#jpeg(file="genotype_clustering_w_trait_heatmap.jpeg")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitDataparents),
                    main = "Sample dendrogram and trait heatmap")

save(geneexpressiondata, traitDataparents, file = "/WGCNABarley/WGCNA_part_1_and_2_power2_with_all_genes/Parent_Leaf_01_dataInput.RData")

#dev.off()





############Step-by-step network construction and module detection
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

#####

sft = pickSoftThreshold(geneexpressiondata, powerVector = powers, verbose = 5)
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
abline(h=0.90,col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



net = blockwiseModules(geneexpressiondata, power = 2,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "PARENT_LEAF_TOM", 
                       verbose = 3)


sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "/WGCNABarley/WGCNA_part_1_and_2_power2_with_all_genes/auto_test_power2_with_all_genes.RData")

######

library("ape")

lnames <- load("/WGCNABarley/WGCNA_part_1_and_2_power2_with_all_genes/Parent_Leaf_01_dataInput.RData")

lnames <- load("/WGCNABarley/WGCNA_part_1_and_2_power2_with_all_genes/auto_test_power2_with_all_genes.RData")

nGenes <- ncol(geneexpressiondata)
nSamples <- nrow(geneexpressiondata)

MEs0 <- moduleEigengenes(geneexpressiondata, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, traitDataparents, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix <-  paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitDataparents),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# Define variable weight containing the weight column of datTrait

weight <- as.data.frame(traitDataparents$FT)
names(weight) = "FT"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(geneexpressiondata, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("MM", modNames, sep="")
names(MMPvalue) <- paste("p.MM", modNames, sep="")

geneTraitSignificance <- as.data.frame(cor(geneexpressiondata, weight, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) <- paste("GS.", names(weight), sep="")
names(GSPvalue) <- paste("p.GS.", names(weight), sep="")



module <- "orange"
column <- match(module, modNames)
moduleGenes <- moduleColors == module


sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)




names(as.data.frame(geneexpressiondata))
length(names(as.data.frame(geneexpressiondata)))


names(as.data.frame(geneexpressiondata))[moduleColors=="orange"]
length(names(as.data.frame(geneexpressiondata))[moduleColors=="orange"])




annot <- read.gff("/WGCNABarley/data/Hv_Morex.pgsb.Jul2020.gff3",)
annot <- annot[annot[,3]=="gene",]

str(annot)

for(i in 1:nrow(annot)){
  annot$attributes[i] <- strsplit(annot$attributes[i], ";")[[1]][1]
}
annot$attributes <- gsub("ID=", "", annot$attributes, fixed = TRUE)
annot <- annot[which(annot$attributes %in% colnames(geneexpressiondata)),]

probes <- names(as.data.frame(geneexpressiondata))

probes2annot <- match(probes, annot$attributes)

sum(is.na(probes2annot))
# Should return 0.


mydata <- read.gff("/WGCNABarley/data/Hv_Morex.pgsb.Jul2020.gff3",)
mydata <- mydata[mydata[,3]=="gene",]


# Create the starting data frame
geneInfo0 <- data.frame(attributes = probes,
                        moduleColor = moduleColors,
                        geneTraitSignificance,
                        GSPvalue)


# Order modules by their significance for weight
modOrder <- order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.FT));
geneInfo = geneInfo0[geneOrder, ]



write.csv(geneInfo, file = "/WGCNABarley/WGCNA_part_1_and_2_power2_with_all_genes/geneInfo_power2.csv")



