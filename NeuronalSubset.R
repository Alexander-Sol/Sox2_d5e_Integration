# Neuronal Sub Clustering

library(Seurat)
library(tidyverse)

otic.reversed <- readRDS("D5e_Sox2_Integrated.rds")

otic.reversed <- FindNeighbors(otic.reversed, 
                               k.param = 25, 
                               dims = 1:17)
otic.reversed <- FindClusters(otic.reversed,
                              resolution = 0.25)
DimPlot(otic.reversed, label = T)

neuro.sub <- subset(otic.reversed, idents = 3)
neuro.sub <- RunPCA(neuro.sub)
neuro.sub <- RunUMAP(neuro.sub, dims = 1:9)

DimPlot(neuro.sub, group.by = "cell.line")
DimPlot(neuro.sub, group.by = "previous_clusters")

neuro.sub$DayCell <- paste0(neuro.sub$Day, neuro.sub$cell.line)
DimPlot(neuro.sub, group.by = "DayCell", pt.size = 1.1)

DefaultAssay(neuro.sub) <- "RNA"
FeaturePlot(neuro.sub, features = "EPCAM", min.cutoff = 0, max.cutoff = 10, pt.size = 1.8, order = T, slot = "counts")



