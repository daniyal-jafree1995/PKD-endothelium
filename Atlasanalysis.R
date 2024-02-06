## scRNA-seq analysis of human kidney vasclar endothelium for Jafree et al. 2023
## Merged dataset derived from multiple scRNA-seq studies of human kidney (see Methods for details)
## Authors: Daniyal Jafree & Mary Ball (University College London)
## Version 1: 07/01/2023

#----------------------------------------------------------------------------------------------------------------#

## Load packages and set working directory. Please change working directory as required.
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
library(scales)

## Load dataset and set cell types as active.ident, before creating UMAP of merged dataset
load("/Users/daniyaljafree/Ben_data/Sobj/combined_data.RData")
seurat.object <- SetIdent(seurat.object, value = seurat.object@meta.data$celltype)
DimPlot(seurat.object, pt.size = 0.5, raster= F, label = F)
table(seurat.object@meta.data$donor_ID)
VlnPlot(seurat.object, features = c("FBLN2", "SSUH2", "GPM6A"), stack = T)
FeaturePlot(seurat.object, features = c("FBLN2"), raster = F, order = F, pt.size = 1)
FeaturePlot(seurat.object, features = c("SSUH2"), raster = F, order = F, pt.size = 1)
FeaturePlot(seurat.object, features = c("GPM6A"), raster = F, order = F, pt.size = 1)

## Derive blood endothelial cells from dataset, annotate accoring to disease status and tidy up using CellSelector
EC <- subset(seurat.object, idents = c("GEC", "PCE", "VRE", "arterial_endothelium", "venular_endothelium")) # Use this subset for downstream analyses

    ## Supplementary figures
    ECidents <- colnames(EC)
    DimPlot(seurat.object, cells.highlight = ECidents, pt.size = 1, raster = F) #FIGURE S1B

## Continue
EC <- SetIdent(EC, value = EC@meta.data$pathology)
pathology.grouping_EC <- c("Ctrl", "Ctrl", "CKD", "Rejection", "Rejection", "Ctrl", "Ctrl", "Ctrl", "Ctrl", "Rejection", "Rejection", "Ctrl")
names(pathology.grouping_EC) <- levels(EC)
EC <- RenameIdents(EC, pathology.grouping_EC)
EC <- SetIdent(EC, value = EC@meta.data$celltype)
DimPlot(EC, pt.size = 0.5, raster= F, label = F)
plot <- DimPlot(EC, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(EC, cells = select.cells) <- "EC_Clean"
EC <- subset(EC, idents = c("EC_Clean"))

## Initial subclustering and DEGs  of isolated EC subsets
EC <- NormalizeData(EC, normalization.method = "LogNormalize", scale.factor = 10000)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(EC), 20)
plot1 <- VariableFeaturePlot(EC)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
all.genes <- rownames(EC)
EC <- ScaleData(EC, features = all.genes)
EC <- RunPCA(EC, features = VariableFeatures(object = EC))
ElbowPlot(EC)
pct <- EC[["pca"]]@stdev / sum(EC[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs
EC <- EC %>%
  RunHarmony("donor_ID", plot_convergence = TRUE, max.iter.harmony = 20)
EC <- FindNeighbors(EC, reduction = 'harmony', dims = 1:14)
EC <- FindClusters(EC, resolution = 0.5)
EC <- RunUMAP(EC, reduction = 'harmony', dims = 1:14)
DimPlot(EC, reduction = "umap", label = TRUE, pt.size = 1)
#EC.markers <- FindAllMarkers(EC, only.pos = T, logfc.threshold = 0.5)
#EC.markers %>% group_by(cluster) %>% top_n(n = 20)
#top20 <- EC.markers %>% group_by(cluster) %>% top_n(n = 20)
#write.csv(EC.markers,"EC.markers.csv", row.names = FALSE)

## Second round of subclustering and DEGs  of isolated EC subsets
EC <- subset(EC, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "11"))
EC <- NormalizeData(EC, normalization.method = "LogNormalize", scale.factor = 10000)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(EC), 20)
plot1 <- VariableFeaturePlot(EC)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
all.genes <- rownames(EC)
EC <- ScaleData(EC, features = all.genes)
EC <- RunPCA(EC, features = VariableFeatures(object = EC))
ElbowPlot(EC)
pct <- EC[["pca"]]@stdev / sum(EC[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs
EC <- EC %>%
  RunHarmony("donor_ID", plot_convergence = TRUE, max.iter.harmony = 30)
EC <- FindNeighbors(EC, reduction = 'harmony', dims = 1:17)
EC <- FindClusters(EC, resolution = 0.5)
EC <- RunUMAP(EC, reduction = 'harmony', dims = 1:17)
DimPlot(EC, reduction = "umap", label = TRUE, pt.size = 1)
#EC2.markers <- FindAllMarkers(EC, only.pos = T, logfc.threshold = 0.5)
#EC2.markers %>% group_by(cluster) %>% top_n(n = 20)
#top20 <- EC2.markers %>% group_by(cluster) %>% top_n(n = 20)
#write.csv(EC2.markers,"EC2.markers.csv", row.names = FALSE)

## Penultimate round of subclustering and DEGs  of isolated EC subsets
EC <- subset(EC, idents = c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11"))
EC <- NormalizeData(EC, normalization.method = "LogNormalize", scale.factor = 10000)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(EC), 20)
plot1 <- VariableFeaturePlot(EC)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
all.genes <- rownames(EC)
EC <- ScaleData(EC, features = all.genes)
EC <- RunPCA(EC, features = VariableFeatures(object = EC))
ElbowPlot(EC)
pct <- EC[["pca"]]@stdev / sum(EC[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs
EC <- EC %>%
  RunHarmony("donor_ID", plot_convergence = TRUE, max.iter.harmony = 30)
EC <- FindNeighbors(EC, reduction = 'harmony', dims = 1:17)
EC <- FindClusters(EC, resolution = 0.8)
EC <- RunUMAP(EC, reduction = 'harmony', dims = 1:17)
DimPlot(EC, reduction = "umap", label = TRUE, pt.size = 1)
plot <- DimPlot(EC, reduction = "umap")
select.cells <- CellSelector(plot = plot)
Idents(EC, cells = select.cells) <- "EC_Clean"
EC <- subset(EC, idents = c("EC_Clean"))

## FINAL round of subclustering and DEGs  of isolated EC subsets, before annotation and final cluster DEGs
EC <- NormalizeData(EC, normalization.method = "LogNormalize", scale.factor = 10000)
EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(EC), 20)
plot1 <- VariableFeaturePlot(EC)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2
all.genes <- rownames(EC)
EC <- ScaleData(EC, features = all.genes)
EC <- RunPCA(EC, features = VariableFeatures(object = EC))
ElbowPlot(EC)
pct <- EC[["pca"]]@stdev / sum(EC[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2
pcs <- min(co1, co2)
pcs
EC <- EC %>%
  RunHarmony("donor_ID", plot_convergence = TRUE, max.iter.harmony = 30)
EC <- FindNeighbors(EC, reduction = 'harmony', dims = 1:18)
EC <- FindClusters(EC, resolution = 1.2)
EC <- RunUMAP(EC, reduction = 'harmony', dims = 1:18)
DimPlot(EC, reduction = "umap", label = TRUE, pt.size = 1)
EC3.markers <- FindAllMarkers(EC, only.pos = T, logfc.threshold = 0.5)
EC3.markers %>% group_by(cluster) %>% top_n(n = 20)
top20 <- EC3.markers %>% group_by(cluster) %>% top_n(n = 20)
write.csv(EC3.markers,"EC3.markers.csv", row.names = FALSE)
EC <- subset(EC, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "15"))
new.cluster.ids <- c("PTC", "PTC", "SELE+ EC", "GEC", "SELE+ EC", "SELE+ EC", "Aff Art", "GEC", "PTC", "Artery", "PTC", "Eff Art", "Eff Art", "DVR", "AVR")
names(new.cluster.ids) <- levels(EC)
EC <- RenameIdents(EC, new.cluster.ids)
saveRDS(EC, file = "/Users/daniyaljafree/Desktop/Papers_and_Thesis/Papers/PKD_multimodal vascular_profiling/EC_atlas.rds")

#----------------------------------------------------------------------------------------------------------------#

EC_atlas <- readRDS(file = "/Users/daniyaljafree/Desktop/Papers_and_Thesis/Papers/PKD_multimodal vascular_profiling/EC_atlas.rds")
DimPlot(EC_atlas, reduction = "umap", label = T, pt.size = 1, label.size = 0)
my_levels <- c("Artery",
               "Aff Art",
               "GEC",
               "Eff Art",
               "PTC",
               "DVR",
               "AVR",
               "SELE+ EC")
EC_atlas@active.ident <- factor(x = EC_atlas@active.ident, levels = my_levels)
#hue_pal()(length(levels(EC_atlas@active.ident)))
DimPlot(EC_atlas, reduction = "umap", label = T, pt.size = 1, label.size = 0, cols = c("Artery" = "#F8766D",
                                                                                 "Aff Art" = "#CD9600",
                                                                                 "GEC" = "#7CAE00",
                                                                                 "Eff Art" = "#00BE67",
                                                                                 "PTC" = "#00BFC4",
                                                                                 "DVR" = "#00A9FF",
                                                                                 "AVR" = "#C77CFF",
                                                                                 "SELE+ EC" = "#FF61CC")) #FIGURE 1B
EC.markers <- FindAllMarkers(EC_atlas, only.pos = T, logfc.threshold = 0.5)
EC.markers %>% group_by(cluster) %>% top_n(n = 20)
top20 <- EC.markers %>% group_by(cluster) %>% top_n(n = 20)
write.csv(EC.markers,"EC.markers.csv", row.names = FALSE)
EC.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(EC_atlas, features = top10$gene) #FIGURE 1C
write.csv(top10,"top10.markers.csv", row.names = FALSE)

FeaturePlot(EC_atlas, features = "FBLN2", order = F, pt.size = 1) #FIGURE 1D
FeaturePlot(EC_atlas, features = "GPM6A", order = F, pt.size = 1) #FIGURE 1E
FeaturePlot(EC_atlas, features = "SSUH2", order = F, pt.size = 1) #FIGURE 1F
DotPlot(EC_atlas, features = c("FBLN2", "SSUH2", "GPM6A"),  split.by = "dataset", cols = c("purple", "purple", "purple", "purple", "purple")) #FIGURE S1A

InjATLAS.markers <- FindMarkers(EC_atlas, ident.1 = "SELE+ EC", min.pct = 0.1, only.pos = T, logfc.threshold = 0)
write.csv(InjATLAS.markers,"InjATLAS.markers.csv")

DotPlot(EC_atlas, features = c("KCNIP4", "SPP1", "ERBB4", "PKHD1", "KAZN", "ZFP36L1", "MACC1", "CLU", "ITGB8", "GLS"))
FeaturePlot(EC_atlas, features = c("KCNIP4", "SPP1", "ERBB4", "PKHD1", "KAZN", "ZFP36L1", "MACC1", "CLU", "ITGB8", "GLS"))

