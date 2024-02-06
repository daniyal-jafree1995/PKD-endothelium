## snRNA-seq analysis of Muto et al. PKD dataset
## Dataset derived from 8 PKD tissues and 5 control tissues
## Authors: Daniyal Jafree (University College London)
## Version 1: 26/09/2022

## Load packages and set working directory. Please change working directory as required.
set.seed(42)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(Matrix)
library(ggrepel)
library(patchwork)
library(tidyselect)
library(ktplots)
library(harmony)
library(SoupX)
library(EnhancedVolcano)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(singleCellNet)

#----------------------------------------------------------------------------------------------------------------#
## Load dataset and set cell types as active.ident, before creating UMAP of merged dataset
Combined <- readRDS("/Users/daniyaljafree/Desktop/scRNAseq_analyses/Muto_PKD/ENDO.rds")
metadata <- Combined@meta.data
DimPlot(Combined, pt.size = 0.5, raster= F, label = F, split.by = "disease")  # Figure 2A

## Normalisation, scaling and PCA with script to empirically decide number of PCs to carry forward (12 according to preliminary analysis)
# Normalisation of RNA counts
Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
# Scaling
all_genes <- rownames(Combined)
Combined <- ScaleData(Combined,features = all_genes)
# PCA
Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
VizDimLoadings(Combined, dims = 1:2, reduction = "pca")
DimPlot(Combined, reduction = "pca")
# Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
ElbowPlot(Combined, ndims = 50)
pct <- Combined[["pca"]]@stdev / sum(Combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs

## Unsupervised clustering using UMAP
# KNN and Louvain clustering prior to UMAP. Normally resolution should be between 0.4-1.2 for single cell datasets
Combined <- FindNeighbors(Combined, dims = 1:15)
Combined <- FindClusters(Combined, resolution = 0.4) # alternative Leiden approach from Gideon: FindClusters(wt1, resolution = 0.9, method = "igraph", algorithm = 4)
# Run and plot UMAP by cell types
Combined <- RunUMAP(Combined, dims = 1:15)
DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "disease") # can use group.by and call variable from object@meta.data to group cells in UMAP by variable

# Create raw count matrix from Seurat object
rawcountmatrix <-Combined@assays$RNA@counts
# Create metadata from Seurat object
Combined_metaData <- as.data.frame(Combined@reductions$umap@cell.embeddings)
colnames(Combined_metaData) <- c('RD1','RD2')
Combined_metaData$Cluster <- Combined@meta.data$seurat_clusters
# Start by generating SoupChannel object assuming no knowlege of empty droplets
scNoDrops <- SoupChannel(rawcountmatrix, rawcountmatrix, calcSoupProfile = FALSE)
soupProf <- data.frame(row.names = rownames(rawcountmatrix), est = rowSums(rawcountmatrix)/sum(rawcountmatrix), counts = rowSums(rawcountmatrix))
scNoDrops <- setSoupProfile(scNoDrops, soupProf)
# Then add metadata including cluster informatiton and dimension reduction
scNoDrops <- setClusters(scNoDrops, setNames(Combined_metaData$Cluster, rownames(Combined_metaData)))
# Add dimension reduction for the data
scNoDrops <- setDR(scNoDrops, Combined_metaData[colnames(rawcountmatrix), c("RD1", "RD2")])
# Estimate contamination fraction utomatically
scNoDrops <- autoEstCont(scNoDrops)
# Generated corrected count matrix
correctedcounts <- adjustCounts(scNoDrops)
plotChangeMap(scNoDrops, correctedcounts, "EMCN")  # Can optionally examine how SoupX has changed expression values for particular genes
# Make new Seurat object after SoupX correction and carry forward metadata from Part 1
Combined <- CreateSeuratObject(correctedcounts, meta.data = metadata) #N.B. metadata is re-used from earlier!

Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
all_genes <- rownames(Combined)
Combined <- ScaleData(Combined,features = all_genes)
Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
VizDimLoadings(Combined, dims = 1:2, reduction = "pca")
DimPlot(Combined, reduction = "pca")
# Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
ElbowPlot(Combined, ndims = 50)
pct <- Combined[["pca"]]@stdev / sum(Combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs

# Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
Combined <- Combined %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)

# Compute UMAP
Combined <- FindNeighbors(Combined, dims = 1:15, reduction = "harmony")
Combined <- FindClusters(Combined, resolution = 0.4)
Combined <- RunUMAP(Combined, dims = 1:15, reduction = "harmony")

# Generate UMAP and group by desired variables
DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)

Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(Combined.markers, file = "Muto_ENDO_unannotated_1stROUND.csv")

Combined <- subset(Combined, idents = c("0", "1", "3", "4", "6", "7"))
Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
# Scaling
all_genes <- rownames(Combined)
Combined <- ScaleData(Combined,features = all_genes)
# PCA
Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
VizDimLoadings(Combined, dims = 1:2, reduction = "pca")
DimPlot(Combined, reduction = "pca")
# Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
ElbowPlot(Combined, ndims = 50)

# Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
Combined <- Combined %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)
# Access and show first five Harmony embeddings for each barcode.
Combined_harmony <- Embeddings(Combined, 'harmony')
Combined_harmony[1:5, 1:5]

# Compute UMAP
Combined <- FindNeighbors(Combined, dims = 1:15, reduction = "harmony")
Combined <- FindClusters(Combined, resolution = 0.4)
Combined <- RunUMAP(Combined, dims = 1:15, reduction = "harmony")
# Generate UMAP and group by desired variables
DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(Combined, pt.size = 0.5, raster= F, label = F, split.by = "disease")

# Differential expression script by Gideon before cell type assignment
Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(Combined.markers, file = "Muto_ENDO_unannotated_2ndROUND.csv")

Combined <- subset(Combined, idents = c("0", "1", "2", "3", "4", "5", "6"))
Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
# Scaling
all_genes <- rownames(Combined)
Combined <- ScaleData(Combined,features = all_genes)
# PCA
Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
VizDimLoadings(Combined, dims = 1:2, reduction = "pca")
DimPlot(Combined, reduction = "pca")
# Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
ElbowPlot(Combined, ndims = 50)

# Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
Combined <- Combined %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 30)

# Compute UMAP
Combined <- FindNeighbors(Combined, dims = 1:13, reduction = "harmony")
Combined <- FindClusters(Combined, resolution = 0.6)
Combined <- RunUMAP(Combined, dims = 1:13, reduction = "harmony")
# Generate UMAP and group by desired variables
DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(Combined, pt.size = 0.5, raster= F, label = T)

# Differential expression script by Gideon before cell type assignment
Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(Combined.markers, file = "Muto_ENDO_unannotated_3rdROUND.csv")


Combined <- subset(Combined, idents = c("0", "1", "2", "3", "4", "5", "7", "8"))
Combined <- NormalizeData(Combined, normalization.method = "LogNormalize", scale.factor = 10000)
Combined <- FindVariableFeatures(Combined, selection.method = "vst", nfeatures = 2000)
# Scaling
all_genes <- rownames(Combined)
Combined <- ScaleData(Combined,features = all_genes)
# PCA
Combined <- RunPCA(Combined, features = VariableFeatures(object = Combined))
VizDimLoadings(Combined, dims = 1:2, reduction = "pca")
DimPlot(Combined, reduction = "pca")
# Empirically determine number of PCs to use for downstream analysis using elbow plot and script to calculate pcs
ElbowPlot(Combined, ndims = 50)

# Run Harmony and perform batch correction, returns plot with number of iterations required for convergence.
Combined <- Combined %>%
  RunHarmony("orig.ident", plot_convergence = TRUE, max.iter.harmony = 20)

# Compute UMAP
Combined <- FindNeighbors(Combined, dims = 1:13, reduction = "harmony")
Combined <- FindClusters(Combined, resolution = 0.6)
Combined <- RunUMAP(Combined, dims = 1:13, reduction = "harmony")
# Generate UMAP and group by desired variables
DimPlot(Combined, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(Combined, pt.size = 0.5, raster= F, label = T)

# Differential expression script by Gideon before cell type assignment
Combined.markers <- FindAllMarkers(Combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(Combined.markers, file = "Muto_ENDO_unannotated_4thROUND.csv")

DimPlot(Combined, pt.size = 0.5, raster= F, label = T, group.by = "subtype")
DimPlot(Combined, pt.size = 0.5, raster= F, label = T, group.by = "disease")

new.cluster.ids <- c("SPP1+ EC",  # 0
                     "GEC",  # 1
                     "Artery",  # 2
                     "AVR",  # 3
                     "DVR", # 4
                     "EC inj",  # 5
                     "PTC",  # 6
                     "Aff/Eff Art")
names(new.cluster.ids) <- levels(Combined)
Combined <- RenameIdents(Combined, new.cluster.ids)
saveRDS(Combined, file = "/Users/daniyaljafree/Desktop/Papers_and_Thesis/Papers/PKD_multimodal vascular_profiling/MUTO_ANALYSIS_FINAL/EC_PKD.rds")
#----------------------------------------------------------------------------------------------------------------#
EC_PKD <- readRDS(file = "/Users/daniyaljafree/Desktop/Papers_and_Thesis/Papers/PKD_multimodal vascular_profiling/MUTO_ANALYSIS_FINAL/EC_PKD.rds")

#Bulk_data <- subset(Combined)
#Idents(Bulk_data) <- Bulk_data@meta.data$disease
#DefaultAssay(Bulk_data) <- "RNA"
#comparebulk <- FindAllMarkers(Bulk_data, min.pct = 0.1, logfc.threshold = 0.25, only.pos = T)
#write.csv(comparebulk, file = "compare_bulk.csv")
#top20 <- comparebulk %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#DoHeatmap(Bulk_data, features = top20$gene) + NoLegend()

my_levels <- c("Artery",
               "Aff/Eff Art",
               "GEC",
               "PTC",
               "DVR",
               "AVR",
               "EC inj",
               "SPP1+ EC")
EC_PKD@active.ident <- factor(x = EC_PKD@active.ident, levels = my_levels)
DimPlot(EC_PKD, reduction = "umap", label = T, pt.size = 1, label.size = 0, cols = c("Artery" = "#F8766D",
                                                                                       "Aff/Eff Art" = "#CD9600",
                                                                                       "GEC" = "#7CAE00",
                                                                                       "PTC" = "#00BFC4",
                                                                                       "DVR" = "#00A9FF",
                                                                                       "AVR" = "#C77CFF",
                                                                                       "EC inj" = "#FF61CC",
                                                                                        "SPP1+ EC" = "grey58"))

FeaturePlot(EC_PKD, features = "FBLN2", order = F, pt.size = 1) #FIGURE 2B
FeaturePlot(EC_PKD, features = "GPM6A", order = F, pt.size = 1) #FIGURE 2B
FeaturePlot(EC_PKD, features = "SSUH2", order = F, pt.size = 1) #FIGURE 2B

DimPlot(EC_PKD,  raster= F, label = T, group.by = "disease", pt.size = 1, label.size = 0)
table(EC_PKD@meta.data$disease)
DimPlot(EC_PKD,  raster= F, label = T, group.by = "orig.ident", pt.size = 1, label.size = 0)
table(EC_PKD@active.ident, EC_PKD@meta.data$disease)

# Draw volcano plot for arterial endothelium using EnhancedVolcano
Artery_data <- subset(EC_PKD, idents = "Artery")
Idents(Artery_data) <- Artery_data@meta.data$disease
Artery.DE.lowthreshold <- FindMarkers(Artery_data,  ident.1 = "PKD", ident.2 = "control", min.pct = 0.1, logfc.threshold = 0.05)
#write.csv(Artery.DE.lowthreshold, "Artery.DE.lowthreshold.csv")
Artery.DE.lowthreshold <- slice(Artery.DE.lowthreshold, -c(1,4,11,20,22,31,32,55,83,93,99,106,122,134,144,152,161,171,180),)
EnhancedVolcano(Artery.DE.lowthreshold,
                lab = rownames(Artery.DE.lowthreshold),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 1,
                pCutoff = 0.05,
                colCustom = NULL,
                boxedLabels = F,
                drawConnectors = F,
                labSize = 4,
                col = c('black', 'black', 'blue', 'red'),
                widthConnectors = 0.5)

# Draw volcano plot for Aff/Eff Art using EnhancedVolcano
AEArt_data <- subset(EC_PKD, idents = "Aff/Eff Art")
Idents(AEArt_data) <- AEArt_data@meta.data$disease
AEArt.DE.lowthreshold <- FindMarkers(AEArt_data,  ident.1 = "PKD", ident.2 = "control", min.pct = 0.1, logfc.threshold = 0.05)
#write.csv(AEArt.DE.lowthreshold, "AEArt.DE.lowthreshold.csv")
AEArt.DE.lowthreshold <- slice(AEArt.DE.lowthreshold, -c(2,3,4,6,8,18,34,38,43,55,92,120,126,178),)
EnhancedVolcano(AEArt.DE.lowthreshold,
                lab = rownames(AEArt.DE.lowthreshold),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 1,
                pCutoff = 0.05,
                colCustom = NULL,
                boxedLabels = F,
                drawConnectors = F,
                labSize = 4,
                col = c('black', 'black', 'blue', 'red'),
                widthConnectors = 0.5)

# Draw volcano plot for GECs using EnhancedVolcano
GEC_data <- subset(EC_PKD, idents = "GEC")
Idents(GEC_data) <- GEC_data@meta.data$disease
GEC.DE.lowthreshold <- FindMarkers(GEC_data,  ident.1 = "PKD", ident.2 = "control", min.pct = 0.1, logfc.threshold = 0.05)
#write.csv(GEC.DE.lowthreshold, "GEC.DE.lowthreshold.csv")
GEC.DE.lowthreshold <- slice(GEC.DE.lowthreshold, -c(5,9,15,19,54,65,71,84,85,115,143,173),)
EnhancedVolcano(GEC.DE.lowthreshold,
                lab = rownames(GEC.DE.lowthreshold),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 1,
                pCutoff = 0.05,
                colCustom = NULL,
                boxedLabels = F,
                drawConnectors = F,
                labSize = 4,
                col = c('black', 'black', 'blue', 'red'),
                widthConnectors = 0.5)

# Draw volcano plot for AVR using EnhancedVolcano
AVR_data <- subset(EC_PKD, idents = "AVR")
Idents(AVR_data) <- AVR_data@meta.data$disease
AVR.DE.lowthreshold <- FindMarkers(AVR_data,  ident.1 = "PKD", ident.2 = "control", min.pct = 0.1, logfc.threshold = 0.05)
#write.csv(AVR.DE.lowthreshold, "AVR.DE.lowthreshold.csv")
AVR.DE.lowthreshold <- slice(AVR.DE.lowthreshold, -c(6,7,9,12,17,33,58,94,130),)
EnhancedVolcano(AVR.DE.lowthreshold,
                lab = rownames(AVR.DE.lowthreshold),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                FCcutoff = 1,
                pCutoff = 0.05,
                colCustom = NULL,
                boxedLabels = F,
                drawConnectors = F,
                labSize = 4,
                col = c('black', 'black', 'blue', 'red'),
                widthConnectors = 0.5)

# Differentially expressed markers of injured clusters for comparison to atlas
Inj.markers <- FindMarkers(EC_PKD, ident.1 = "EC inj", min.pct = 0.1, only.pos = T, logfc.threshold = 0)
write.csv(Inj.markers,"EC_PKD_inj.csv")

#Ann.markers <- FindAllMarkers(EC_PKD, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Ann.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
#write.csv(Ann.markers, file = "Ann.markers.csv")

SPP1_enriched <- FindMarkers(EC_PKD, ident.1 = "SPP1+ EC", min.pct = 0.1, only.pos = T, logfc.threshold = 0)
write.csv(SPP1_enriched, file = "SPP1_enriched.csv")
DotPlot(EC_PKD, features = c("SERPINE2",
                             "SLC6A6", "SSUH2",
                             "EHD3", "CRHBP",
                             "RGCC", "PLVAP",
                             "FBLN2", "ENPP2",
                             "GPM6A", "EDIL3"))
DotPlot(EC_PKD, features = c("KCNIP4", "SPP1", "ERBB4", "PKHD1", "KAZN", "ZFP36L1", "MACC1", "CLU", "ITGB8", "GLS"))
DotPlot(EC_PKD, features = c("ACKR1", "SELE", "VCAM1", "ICAM1", "CPE", "KCNIP4", "SPP1", "ERBB4", "PKHD1", "KAZN"))

# Differential abundance analysis of EC subsets using MiloR
  # Extracts SingleCellExperiment object from Seurat before following workflow for MiloR: https://bioconductor.org/packages/release/bioc/vignettes/miloR/inst/doc/milo_demo.html 
celltypeDJ <- Idents(EC_PKD)
EC_PKD <- AddMetaData(EC_PKD, metadata = celltypeDJ, col.name = "celltypeDJ")
EC_PKD_sce <- as.SingleCellExperiment(EC_PKD)
traj_milo <- Milo(EC_PKD_sce)
traj_milo <- buildGraph(traj_milo, k = 10, d = 30)
traj_milo <- makeNhoods(traj_milo, prop = 0.1, k = 10, d=30, refined = TRUE)
plotNhoodSizeHist(traj_milo)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="orig.ident")
head(nhoodCounts(traj_milo))
traj_design <- data.frame(colData(traj_milo))[,c("orig.ident", "disease")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]
traj_design
traj_milo <- calcNhoodDistance(traj_milo, d=30)
rownames(traj_design) <- traj_design$orig.ident
da_results <- testNhoods(traj_milo, design = ~ disease, design.df = traj_design)
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 
traj_milo <- buildNhoodGraph(traj_milo)
plotUMAP(traj_milo) + plotNhoodGraphDA(traj_milo, da_results, alpha=0.1) +
  plot_layout(guides="collect")
da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "celltypeDJ")
head(da_results)
ggplot(da_results, aes(celltypeDJ_fraction)) + geom_histogram(bins=50)
da_results$celltypeDJ <- ifelse(da_results$celltypeDJ_fraction < 0.7, "Mixed", da_results$celltypeDJ)
plotDAbeeswarm(da_results, group.by = "celltypeDJ", alpha = 0.1)


library(singleCellNet)


## SingleCellNet analysis of Muto dataset vs annotated atlas.
# Replace all instances of "EC_atlas" with the Seurat object name of the dataset you want to set as the reference (so, in this case, don't change anything)
# Replace all instances of "EC_PKD" with the Seurat object name of the dataset you want to set as the query or test (so, in this case, change it to whatever your object is called)

EC_PKD <- readRDS(file = "/Users/daniyaljafree/Desktop/Papers_and_Thesis/Papers/PKD_multimodal vascular_profiling/MUTO_ANALYSIS_FINAL/EC_PKD.rds") #data to learn from

# Set query data
expr_query <- EC_PKD@assays$RNA@counts  #raw count matrix of data to test
EC_PKD$CellID <- colnames(EC_PKD) # "CellID" is the barcode of each individual droplet
EC_PKD$newAnn <- Idents(EC_PKD) # "newAnn" is the annotation I assigned for the cell types in the query dataset
st_query <- EC_PKD@meta.data
genes_query <- rownames(EC_PKD)

# Set training data
exp_trainingraw <- EC_atlas@assays$RNA@counts
EC_atlas$CellID <- colnames(EC_atlas) # "CellID" is the barcode of each individual droplet
EC_atlas$newAnn <- Idents(EC_atlas) # "newAnn" is the annotation I assigned for the cell types in the query dataset
st_training <- EC_atlas@meta.data
st_training <-droplevels(st_training)

# Find common genes and limit analysis to these
commonGenes = intersect(rownames(exp_trainingraw), genes_query)
length(commonGenes)
exp_trainingraw = exp_trainingraw[commonGenes,]

# Split training and assessment and transform training data
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=st_training, ncells=100, dLevel="newAnn")
stTrain = stList[[1]]
expTrain = exp_trainingraw[,rownames(stTrain)]

# Train classifier and validate
system.time(class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "newAnn", colName_samp = "CellID"))
stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel="newAnn") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = exp_trainingraw[commonGenes,rownames(stTest)]
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "CellID", classTrain = "newAnn", classQuery = "newAnn", nRand = 50)
plot_PRs(tm_heldoutassessment)
nrand = 50
sla = as.vector(stTest$newAnn)
names(sla) = as.vector(stTest$CellID)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created
sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="newAnn", sid="CellID")

# Run classifier
nqRand = 50
system.time(cr_all<-scn_predict(class_info[['cnProc']], expr_query, nrand=nqRand))
sgrp = as.vector(st_query$newAnn)
names(sgrp) = as.vector(st_query$CellID)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)
sc_hmClass(cr_all, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
