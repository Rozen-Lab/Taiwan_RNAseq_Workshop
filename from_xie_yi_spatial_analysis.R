
# Please use following the code to test if package has been installed:
library(Seurat)

devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

options(timeout = 600000000) ### set this to avoid timeout error
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

library(spacexr)
library(Rfast)
library(ape)

# Missing package can be installed by using install.packages()​
# spacexr​ need to be installed with 

# This is the link from the tutorial to
# The code below expects it to be in "../data/mouse_hippocampus_reference.rds"
# https://www.dropbox.com/s/cs6pii5my4p3ke3/mouse_hippocampus_reference.rds?dl=0

# This is the link from the tutorial to Allen Brain atlas cortex data.
# The below expects it be in ../data/allen_cortex.rds
# 
# https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1


# Link to the tutorial:
# https://satijalab.org/seurat/articles/spatial_vignette
# setwd('C:/Users/e0205142/OneDrive - National University of Singapore/Others/tw workshop/spatial_Seurat/script/')

rm(list = ls())

SeuratData::InstallData("stxBrain")

# 10x Visium--------------
## Dataset-------------------
brain <- LoadData("stxBrain", type = "anterior1")

## Data preprocessing-------------
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

## Gene expression visualization----------------
## 1. plot in Rstudio with default parameters
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

### 2. plot and export figure as jpg
plot <- SpatialFeaturePlot(brain, features = c("Ttr")) +
  theme(legend.text = element_text(size = 0),
        legend.title = element_text(size = 20), 
        legend.key.size = unit(1, "cm"))
# Manually make sure ./output/images/" exists.
# This is changed from the original vignette, where it was ../output
jpeg(filename = "./output/images/spatial_vignette_ttr.jpg",
     height = 700, width = 1200, quality = 50)
print(plot)
dev.off()

### 3. plot in Rstudio with different parameters
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2


## Dimensionality reduction, clustering, and visualization-------------
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)
##### WATCH OUT!
### Interactive plotting---------------------
SpatialDimPlot(brain, interactive = TRUE)

SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)

LinkedDimPlot(brain)

## Identification of Spatially Variable Features----------
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

## Subset out anatomical regions-----------------
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2


## Integration with single-cell data--------------
allen_reference <- readRDS("../data/allen_cortex.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay


DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)


cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi",
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex, selection.method = "moransi"), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)


SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT",
                                        "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))


## Working with multiple slices in Seurat-------------
brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

brain.merge <- merge(brain, brain2)


DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)

DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))

SpatialDimPlot(brain.merge)

SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))



# Slide-seq----------------------
## Dataset--------------
InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")

## Data preprocessing------------
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)

plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2

SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c(1,
                                                                                             6, 13)), facet.highlight = TRUE)

## Integration with a scRNA-seq reference-----------
ref <- readRDS("../data/mouse_hippocampus_reference.rds")
ref <- UpdateSeuratObject(ref)

anchors <- FindTransferAnchors(reference = ref, query = slide.seq, normalization.method = "SCT",
                               npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE,
                                  weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay


DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Dentate Principal cells", "CA3 Principal cells", "Entorhinal cortex",
                                           "Endothelial tip", "Ependymal", "Oligodendrocyte"), alpha = c(0.1, 1))

slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c("CA3 Principal cells",
                                                                                             "Dentate Principal cells", "Endothelial tip")), facet.highlight = TRUE)

## Identification of Spatially Variable Features----------
DefaultAssay(slide.seq) <- "SCT"
slide.seq <- FindSpatiallyVariableFeatures(slide.seq, assay = "SCT", slot = "scale.data", features = VariableFeatures(slide.seq)[1:1000],
                                           selection.method = "moransi", x.cuts = 100, y.cuts = 100)

## This step will take ~ 1hr
SpatialFeaturePlot(slide.seq, features = head(SpatiallyVariableFeatures(slide.seq, selection.method = "moransi"),
                                              6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")
