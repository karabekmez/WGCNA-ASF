library(WGCNA)
options(stringsAsFactors = FALSE);

MyData <- read.csv("COVID1CSV-ASF.csv", header=FALSE, sep=",")
adjacency <- as.matrix(MyData)
 		
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM;
geneTree = hclust(as.dist(dissTOM), method = "average");
sizeGrWindow(12,9)
tiff("covid1 Gene clustering on TOM-based dissimilarity ", units="in", width=5, height=5, res=600)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
	labels = FALSE, hang = 0.04);
dev.off()

minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
	deepSplit = 2, pamRespectsDendro = FALSE,
	minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8,6)

tiff("covid1 Gene dendrogram and module colors ", units="in", width=5, height=5, res=600)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
	dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05,

	main = "Gene dendrogram and module colors")
dev.off()

# after initial modulation of WGCNA, re-allocation of genes by using raw expression data is needed
	human <- read.table("human.txt") # human.txt consists of expression data
	covidmatrix <- as.matrix(human)
	datExpr <- t(covidmatrix)	

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-adjacency;	
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

tiff("covid1 Clustering of module eigengenes ", units="in", width=5, height=5, res=600)

plot(METree, main = "Clustering of module eigengenes",
	xlab = "", sub = "")
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)

tiff("covid1 DendroAndColors ", units="in", width=5, height=5, res=600)

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)

dev.off()

# Rename moduleColors
moduleColors = mergedColors

mergedColors
write.csv(mergedColors,file = "covid1 mergedcolors.csv")

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));	# if the number of modules is greater than 50 it should be altered
moduleLabels = match(moduleColors, colorOrder)-1;

moduleLabels 
write.csv(moduleLabels,file = "covid1 modulelabels.csv")

MEs = mergedMEs;

mergedMEs
write.csv(mergedMEs,file = "covid1 mergedMEs.csv")
