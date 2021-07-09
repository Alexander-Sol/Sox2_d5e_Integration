# Trial integration of 
# - Sox 2 Day 20, 30, 40, and 60 Merge
# - D5e Day 60 (EPCAM + tdT - 7AAD sorted)
#
# NOTE: Both datasets are included within the project directory as .rds files
#       However, they've been excluded from the git repo via .gitignore

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
otic.list <- list(Sox2, D5e)
otic.list <- lapply(otic.list, SCTransform)
otic.features <- SelectIntegrationFeatures(object.list = otic.list,
                                           nfeatures = 3000)
otic.list <- PrepSCTIntegration(object.list = otic.list,
                                anchor.features = otic.features)
otic.anchors <- FindIntegrationAnchors(object.list = otic.list,
                                       normalization.method = "SCT",
                                       anchor.features = otic.features)
otic.combined <- IntegrateData(anchorset = otic.anchors,
                               normalization.method = "SCT")
rm(otic.anchors, otic.list)

# Integration complete, now onto the standard workflow
otic.combined <- RunPCA(otic.combined)
otic.combined <- RunUMAP(otic.combined, reduction = "pca", dims = 1:15)
saveRDS(otic.combined, "Sox2_D5e_Integrated.rds")

# Quick Examination Via Plotting
PlotbyDay <- DimPlot(otic.combined, group.by = "Day")
PlotbyCellLine <- DimPlot(otic.combined, group.by = "cell.line")

otic.combined$previous_clusters <- otic.combined$seurat_clusters
otic.combined$previous_clusters[otic.combined$cell.line == "Sox2"] <- "Sox2"
DanClusterOverlay <- DimPlot(otic.combined, group.by = "previous_clusters", label = T, order = as.character(0:12), pt.size = 1.2,
        cols = c("lightgrey", rainbow(13, 0.6, 0.9)))

#Save plots for sharing
plotList <- list(
  byDay = PlotbyDay,
  byCellLine = PlotbyCellLine,
  clusterOverlay = DanClusterOverlay
)

for(i in seq_len(length(plotList))) {
  ggsave(filename = paste0("Plot", names(plotList)[i], ".png"),
         plot = plotList[[i]],
         device = "png",
         units = "in",
         height = 6,
         width = 6)
}

