library(WGCNA)
library(DESeq2)
options(stringsAsFactors = F)
library(preprocessCore)
#library(clusterProfiler)
#===============================================================================
#
#  Read the gene counts table and plot the sample tree
#
#===============================================================================
setwd("/projects/lu_lab/Gaoshan/MIA/RNA_seq/outputBA1")
sample_metadata = read.csv(file = "sampleBA1.csv")

# # Read traits data
traits = read.csv("traitBA1.csv", header = T)
traits_ori = read.csv("traitBA1.csv", header = T)

# Read the gene counts table 
group=read.table('total_count.txt',header=T)
genelength<-group$Length
rownames(group)<- group[,1]
group_A1<- group[,19:30]
nameA1<-colnames(group_A1)
group_B<-group[,37:48]
nameB<-colnames(group_B)

# Combine data and gene length
group_combine<-cbind(genelength,group_B,group_A1)
enhancer_def = data.frame(group_combine)

## Filter low quality data
# Discard the rows with >70% of counts is less than 5
miss <- c()
for (i in 1:nrow(enhancer_def)) {
  if(length(which(enhancer_def[i,]<5)) > 0.7*ncol(enhancer_def)) miss <- append(miss,i)
}
enhancer_def = enhancer_def[-miss,]

## Prepare the dataframe
# Obtain filtered gene length
fil_genelength<-enhancer_def[,1]
# Remove gene length from data matrix
data0<-(enhancer_def[,-1])

## Normalization with log2(FPKM+1)
# Get the coldata sheet
# sample_metadata = read.csv(file = "sampleBA.csv")
# Construct DESeq2 datamatrix
dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0,colData = sample_metadata,design = ~ Condition)
# Add filtered gene length to datamatrix (required for FPKM calculation)
mcols(dataExpr_deseq)$basepairs = fil_genelength
# Calculate FPKM
fpkm_matrix = fpkm(dataExpr_deseq)
# Normorlize FPKM by log 
dat_norm = log2(fpkm_matrix+1)

## Remove genes with MAD<0.3
#apply mad function, 1 means manipulation is performed on rows while 2 means on columns
#mad: Median Absolute Deviation
m.mad <- apply(dat_norm,1,mad)
# quantile(m.mad, probs = seq(0, 1, 0.3))
# m.mad<-na.omit(m.mad)
# max(quantile(m.mad, probs = seq(0,1,0.3))[2], 0.01)
dataExprMad <- dat_norm[which(m.mad > max(quantile(m.mad, probs = seq(0,1,0.3))[2],0.01)),]

## Construct dataframe for network analysis
# Get the transpose matrix!!!
datExpr <- as.data.frame(t(dataExprMad))
row.names(datExpr)=sample_metadata[,1]

## Check sample quality
gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}


# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average");

# plot sample tree
pdf(file = "1-n-sampleClustering.pdf", width = 12, height = 9);
par(cex = 1.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 
# Scale-free topology fit index as a function of the soft-thresholding power
pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
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
dev.off()


#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate
#power=7
TOM = TOMsimilarityFromExpr(datExpr, power = power)
dissTOM = 1-TOM 

# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()


#===============================================================================
#
#  Construct modules automatically
#
#===============================================================================
nGenes <- ncol(datExpr)
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = power,
                       corType='bicor',
                       TOMType = "signed", minModuleSize = 50,
                       maxBlockSize = nGenes,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<- stats::cor
# unsigned -> nodes with positive & negative correlation are treated equally
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)

# Show the moudules table
table(mergedColors)
pdf(file = "4-module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#===============================================================================
#
#  Construct modules Manually
#
#===============================================================================

# Module identification using dynamic tree cut
# dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
#                             pamRespectsDendro = FALSE,minClusterSize = 30);
# table(dynamicMods)
# length(table(dynamicMods)) 
# # Convert numeric labels into colors
# dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors)
# # Plot the dendrogram and colors underneath
# pdf(file = "4-module_tree.pdf", width = 8, height = 6);
# plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
#                     hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
# #dev.off()
# 
# #===============================================================================
# #
# #  Merge modules
# #
# #===============================================================================
# 
# # Merge close modules
# MEDissThres=0.25
# abline(h=MEDissThres, col = "red")
# merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
# mergedColors = merge$colors  
# mergedMEs = merge$newMEs  
# table(mergedColors)
# # Plot merged module tree
# pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
#                     c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
#                     hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
# dev.off()
# # write.table(merge$oldMEs,file="oldMEs.txt");
# # write.table(merge$newMEs,file="newMEs.txt");

#===============================================================================
#
#  Plot the heatmap of module eigen-genes and samples
#
#===============================================================================

# library("pheatmap")
# 
# # Heatmap of old module eigen-genes and samples
# pdf(file="oldMEs.pdf")
# # merge_old<-merge$oldMEs
# # row<-row.names(merge_old)
# # dim(merge_old)
# 
# rm<-colnames(mdata)
# orname<-colnames(data0)
# 
# row.names(merge$oldMEs)=colnames(data0)
# pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=T,show_colnames=T,fontsize=6)
# dev.off()
# # Heatmap of new module eigen-genes and samples
# pdf(file="newMEs.pdf",heigh=60,width=20)
# row.names(merge$newMEs)=colnames(data0)
# pheatmap(merge$newMEs,cluster_col=T,cluster_row=T,show_rownames=T,show_colnames=T,fontsize=6)
# dev.off()

#=====================================================================================
#
#  Correlation between gene modules and experimental traits
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
dim(datExpr)
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# # Read traits data
# traits = read.csv("traitBA.csv", header = T)
# traits_ori = read.csv("traitBA.csv", header = T)
# Make the first column as row names
rownames(traits) = traits[, 1]
robac<-rownames(traits)
traits = traits[, -1]
# Rownames of traits and eigengens should be the same
traits = traits[match(rownames(MEs), rownames(traits)), ]
table(rownames(MEs) == rownames(traits))
# Calculate 'pearson correlation' coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
# write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");


#=====================================================================================
#
#  Plot heatmap of module-traits relationship
#
#=====================================================================================

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("module-traits-auto_filter.pdf", width = 100, height = 30)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#  Export edge and node files for cytoscape and Gene lists for GO term
#
#=====================================================================================
table(mergedColors)
modtable<-unique(mergedColors)
probes = names(datExpr)
# Set the moudules to be exported
## Export certain modules
#intModules = c("brown", "blue", "yellow")
## Export all modules
intModules=modtable
for (modules in intModules) {
  ninModule = (mergedColors== modules);
  #ninModule = is.finite(match(moduleColors, modules));
  modProbes=probes[ninModule]

  # Export gene list for GO term
  write.csv(modProbes, paste('gene_auto_7_filter_',modules,'.csv'))

  # Get top 100 hub genes and export to cytoscape
  nTOP = 100
  IMConn = softConnectivity(datExpr[,modProbes])
  top = (rank(-IMConn) <=nTOP)
  topgenes = modProbes[top]
  # Select the corresponding Topological Overlap
  #length(TOM)
  TOM.mat = as.matrix(TOM)
  modTOM = TOM.mat[ninModule, ninModule]
  dimnames(modTOM) = list(modProbes,modProbes)

  topTOM = modTOM[top,top]
  dimnames(topTOM) = list(topgenes, topgenes)
  #  median(topTOM)
  #  quantile(topTOM)
  cyt = exportNetworkToCytoscape(
    topTOM,
    edgeFile = paste("Cytoscape-edges_auto_filter_", paste(modules, collapse="-"), ".txt", sep=""),
    nodeFile = paste("Cytoscape-nodes_auto_filter_", paste(modules, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.15,
    nodeNames = topgenes,
    nodeAttr = mergedColors[ninModule][top]
  )
  nodelist<- read.delim(paste("Cytoscape-nodes_auto_filter_", paste(modules, collapse="-"),".txt", sep=""))
  write.csv(nodelist,paste('nodelist_auto_7_filter',modules,'.csv'))
}

