# Reverse Integration Test lmao

library(Seurat)
library(tidyverse)

Sox2 <- readRDS("day20304060.rds") %>% UpdateSeuratObject() # UpdateSeuratObject fixes "no slot of name "images"" error
D5e <- readRDS("DanDay60Preprocessed.rds") %>% UpdateSeuratObject() 

downsample <- function(Seurat.object, ds.size) {
  subset( 
    Seurat.object, 
    cells = sample(
      Cells( Seurat.object ),
      ds.size
    )
  )
}

Sox2 <- downsample(Sox2, 12000)
D5e <- downsample(D5e, 3000)

# Labelling the cell line will allow us to separate them out later
Sox2$cell.line <- "Sox2"
D5e$cell.line <- "D5e"
D5e$Day <- "60"

# Both datasets have been normalized with SCTransform, which we don't want.
# So, need to work with raw data from here on out
DefaultAssay(Sox2) <- "RNA"
DefaultAssay(D5e) <- "RNA"

# Let's go one step further and erase the SCT data
Sox2@assays$SCT <- NULL
D5e@assays$SCT <- NULL

# Now we'll begin the work of integration. Following the vignette found here: https://satijalab.org/seurat/articles/integration_introduction.html#performing-integration-on-datasets-normalized-with-sctransform-1
otic.list <- list(Sox2, D5e) %>% rev()
otic.list <- lapply(otic.list, SCTransform)
otic.features <- SelectIntegrationFeatures(object.list = otic.list,
                                           nfeatures = 3000)
otic.list <- PrepSCTIntegration(object.list = otic.list,
                                anchor.features = otic.features)
otic.anchors <- FindIntegrationAnchors(object.list = otic.list,
                                       normalization.method = "SCT",
                                       anchor.features = otic.features)
otic.reversed <- IntegrateData(anchorset = otic.anchors,
                               normalization.method = "SCT")
rm(otic.anchors, otic.list)

otic.reversed <- RunPCA(otic.reversed)
otic.reversed <- RunUMAP(otic.reversed, reduction = "pca", dims = 1:17)

# Save/Load
saveRDS(otic.reversed, "D5e_Sox2_Integrated.rds")
otic.reversed <- readRDS("D5e_Sox2_Integrated.rds")

# Quick Examination Via Plotting
PlotbyDay <- DimPlot(otic.reversed, group.by = "Day")
PlotbyCellLine <- DimPlot(otic.reversed, group.by = "cell.line")

otic.reversed$previous_clusters <- otic.reversed$seurat_clusters
otic.reversed$previous_clusters[otic.reversed$cell.line == "Sox2"] <- "Sox2"
DanClusterOverlay <- DimPlot(otic.reversed, group.by = "previous_clusters", label = T, order = as.character(0:12), pt.size = 1.2,
                             cols = c("lightgrey", rainbow(13, 0.6, 0.9)))

#Save plots for sharing
plotList <- list(
  byDay = PlotbyDay,
  byCellLine = PlotbyCellLine,
  clusterOverlay = DanClusterOverlay
)

for(i in seq_len(length(plotList))) {
  ggsave(filename = paste0("RevInt_Plot", names(plotList)[i], ".png"),
         plot = plotList[[i]],
         device = "png",
         units = "in",
         height = 10,
         width = 10)
}