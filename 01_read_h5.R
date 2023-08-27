#================================== Chapter 2.2.4 ====================================

##### Load required packages
rm(list = ls())
library(data.table)
library(Matrix)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(Seurat)
library(clustree)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
if(!dir.exists("images/")){dir.create("images/")} # Folder to store outputs

### Define color palettes and plot themes
colLib = brewer.pal(7, "Paired")                 # color for libraries
names(colLib) = c("BM0pos", "BM0neg", "BM7pos", "BM7neg", 
                  "BM9pos", "BM9neg", "BM2pos")
colDnr = colLib[c(1,3,5,7)]                      # color for donors
names(colDnr) = c("BM0", "BM7", "BM9", "BM2")
colGEX = c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
colCcy = c("black", "blue", "darkorange")        # color for cellcycle phase
plotTheme <- theme_classic(base_size = 18)
nC <- 4 # Number of threads / cores on computer


### Read in files
# Read in h5 files which are in the `0data/` directory
inpUMIs <- NULL
for(i in names(colLib)){
  tmp <- Read10X_h5(paste0("0data/", i, ".h5"))
  colnames(tmp) <- paste0(i, "_", colnames(tmp)) # Append sample prefix
  inpUMIs <- cbind(inpUMIs, tmp)
}
rm(tmp)

#==================================Chapter 2.3.4 ====================================
### A. Perform QC
# Compute cell QC metrics
oupQCcell <- data.table(
  sampleID = colnames(inpUMIs),
  library = tstrsplit(colnames(inpUMIs), "_")[[1]],
  nUMI = colSums(inpUMIs),
  nGene = colSums(inpUMIs != 0),
  pctMT = colSums(inpUMIs[grep("^MT-", rownames(inpUMIs)),]))
oupQCcell$pctMT <- 100 * oupQCcell$pctMT / oupQCcell$nUMI
oupQCcell$library <- factor(oupQCcell$library, levels = names(colLib))

# Plot cell QC metrics (all output goes into images folder)
p1 <- ggplot(oupQCcell, aes(nUMI, fill = library)) + 
  geom_histogram(binwidth = 500, color = "black") + xlim(c(0,40e3)) + 
  geom_vline(xintercept = c(1.5e3,30e3), color = "red", linetype = "dashed")+ 
  xlab("No. UMIs") + scale_fill_manual(values = colLib) + plotTheme
p2 <- ggplot(oupQCcell, aes(nGene, fill = library)) + 
  geom_histogram(binwidth = 100, color = "black") + xlim(c(0,6e3)) + 
  geom_vline(xintercept = c(0.5e3), color = "red", linetype = "dashed")+ 
  xlab("No. Detected Genes") + scale_fill_manual(values = colLib) + plotTheme
p3 <- ggplot(oupQCcell, aes(pctMT, fill = library)) + 
  geom_histogram(binwidth = 0.5, color = "black") + xlim(c(0,20)) + 
  geom_vline(xintercept = c(15), color = "red", linetype = "dashed")+ 
  xlab("Percentage MT Genes") + scale_fill_manual(values = colLib) + plotTheme
ggsave(p1 + p2 + p3 + guide_area() + plot_layout(nrow = 2, guides = "collect"), 
       width = 10, height = 8, filename = "images/basicCellQC.png")

# Remove low quality cells (39548 cells down to 33852 cells post- QC)
oupQCcell <- oupQCcell[nUMI > 1.5e3]
oupQCcell <- oupQCcell[nUMI < 30e3]
oupQCcell <- oupQCcell[nGene > 0.5e3]
oupQCcell <- oupQCcell[pctMT < 15]
inpUMIs <- inpUMIs[, as.character(oupQCcell$sampleID)]

# Compute gene QC metrics
inpUMIs <- inpUMIs[rowSums(inpUMIs) > 0,]   # Remove genes with no reads
oupQCgene <- data.table(
  gene = rownames(inpUMIs),
  meanRead = rowSums(inpUMIs) / ncol(inpUMIs),
  cellExpr = rowSums(inpUMIs > 0))
oupQCgene$log10meanRead <- log10(oupQCgene$meanRead)
oupQCgene$log10cellExpr <- log10(oupQCgene$cellExpr)

# Plot gene QC metrics
p1 <- ggplot(oupQCgene, aes(log10meanRead)) + 
  geom_histogram(binwidth = 0.05, color = "black") + 
  # geom_vline(xintercept = c(-3), color = "red", linetype = "dashed")+ 
  xlab("log10(Average UMIs)") + plotTheme
p2 <- ggplot(oupQCgene[cellExpr != 0], aes(log10cellExpr)) + 
  geom_histogram(binwidth = 0.05, color = "black") + 
  geom_vline(xintercept = log10(8), color = "red", linetype = "dashed")+ 
  xlab("log10(No. Cells Expressing Gene)") + plotTheme
ggsave(p1 + p2, width = 10, height = 4, filename = "images/basicGeneQC.png")

# Remove lowly expressed genes (24660 genes down to 19677 genes post- QC)
oupQCgene <- oupQCgene[cellExpr >= 8]
inpUMIs <- inpUMIs[as.character(oupQCgene$gene), ]

#================================== Chapter 2.4.5 ====================================
### B. Create Seurat object + Preprocessing
# Create Seurat Object
# nCount_RNA = & no. UMI counts; nFeature_RNA = no. detected genes
# We will also compute percentage MT genes here and set the library
seu <- CreateSeuratObject(inpUMIs, project = "bm")
colnames(seu@meta.data)
seu$pctMT <- 100 * colSums(inpUMIs[grep("^MT-",rownames(inpUMIs)),]) 
seu$pctMT <- seu$pctMT / seu$nCount_RNA

# Create library and donor column
seu$library = factor(seu$orig.ident, levels = names(colLib))
Idents(seu) <- seu$library
seu$donor = gsub("pos|neg", "", seu$library)
seu$donor = factor(seu$donor)

# Normalization + Feature selection
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)

# Plot gene expression mean vs variance
ggData = seu@assays$RNA@meta.features
p1 <- ggplot(ggData, aes(vst.mean, vst.variance, color = vst.variable)) + 
  geom_point() + scale_color_manual(values = c("black", "firebrick")) + 
  geom_line(aes(vst.mean, vst.variance.expected), color = "dodgerblue") + 
  xlab("Average Expression") + ylab("Gene Variance") + 
  scale_x_log10() + scale_y_log10() + plotTheme
p2 <- ggplot(ggData, aes(vst.mean, vst.variance.standardized, 
                         color = vst.variable)) + 
  geom_point() + scale_color_manual(values = c("black", "firebrick")) + 
  xlab("Average Expression") + ylab("Standardized Variance") + 
  scale_x_log10() + plotTheme + theme(legend.position = "none")
p2 <- LabelPoints(plot = p2, repel = TRUE, points = VariableFeatures(seu)[1:10])
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 4, filename = "images/basicHVG.png")

#================================== Chapter 2.5.6 ====================================
### C. PCA / tSNE / UMAP
# Run PCA
seu <- ScaleData(object = seu)    # Scale data prior to PCA
seu <- RunPCA(seu)

p1 <- DimPlot(seu, reduction = "pca", pt.size = 0.1, shuffle = TRUE, 
              cols = colLib) + plotTheme + coord_fixed()
p2 <- DimPlot(seu, reduction = "pca", pt.size = 0.1, shuffle = TRUE,
              cols = colLib, dims = c(1,3)) + plotTheme + coord_fixed()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 4, filename = "images/basicPCA.png")

# Elbow plot
p1 <- ElbowPlot(seu, ndims = 40)
ggData <- seu@reductions[["pca"]]@stdev  # Fit lines for elbow plot
lm(y~x, data = data.frame(x = 16:20, y = ggData[16:20]))
lm(y~x, data = data.frame(x = 31:35, y = ggData[31:35]))
nPC <- 23 # determined from elbow plot
p1 <- p1 + plotTheme + geom_vline(xintercept = nPC, color = "firebrick") + 
  geom_path(data = data.frame(x = 16:35, y = 16:35*-0.07535 + 3.50590),
            aes(x,y), color = "dodgerblue") +
  geom_path(data = data.frame(x = 16:35, y = 16:35*-0.0282 + 2.3996),
            aes(x,y), color = "dodgerblue")
ggsave(p1, width = 5, height = 4, filename = "images/basicPCAelbow.png")

# Run tSNE
seu <- RunTSNE(seu, dims = 1:nPC, num_threads = nC)
p1 <- DimPlot(seu, reduction = "tsne", pt.size = 0.1, shuffle = TRUE, 
              cols = colLib) + plotTheme + coord_fixed()

# Run UMAP
seu <- RunUMAP(seu, dims = 1:nPC)
p2 <- DimPlot(seu, reduction = "umap", pt.size = 0.1, shuffle = TRUE, 
              cols = colLib) + plotTheme + coord_fixed()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 4, filename = "images/basicTsUm.png")

# Plot some genes from literature
p1 <- FeaturePlot(seu, reduction = "umap", pt.size = 0.1,
                  features = c("CD34","CRHBP","GATA1",
                               "CD14","IRF8","CD19",
                               "CD4","CD8B","GZMB"), order = TRUE) &
  scale_color_gradientn(colors = colGEX) & plotTheme & coord_fixed()
ggsave(p1, width = 15, height = 12, filename = "images/basicUmapGex.png")

# Save Seurat Object at end of each section
saveRDS(seu, file = "bmSeu.rds")

#==================================Chapter 3====================================
### A. Clustering
# Perform clustering
seu <- FindNeighbors(seu, dims = 1:nPC)
seu <- FindClusters(seu, resolution = seq(0.4,1.5,0.1), verbose = FALSE)

# Find ideal resolution
p1 <- clustree(seu)
ggsave(p1 + scale_color_manual(values = brewer.pal(12, "Paired")) + 
         guides(color = guide_legend(ncol = 2)), 
       width = 10, height = 8, filename = "images/clustClustree.png")
optRes <- 1.0                # determined from clustering tree
seu$seurat_clusters <- NULL  # Remove this column to prevent confusion
seu$cluster <- seu[[paste0("RNA_snn_res.", optRes)]] 
reorderCluster = c("4","9","12","13","20","14","16",  # Prog
                   "17","10","18","23",               # Ery
                   "11","24","25","7","21","28",      # Myeloid
                   "22","15","8","27","26",           # B-Lymphoid
                   "5","1","0","19","6","3","2")      # T-Lymphoid
seu$cluster = factor(seu$cluster, levels = reorderCluster)
Idents(seu) <- seu$cluster     # Set seurat to use optRes

# Plot ideal resolution on tSNE and UMAP
nClust <- uniqueN(Idents(seu))         # Setup color palette
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)
p1 <- DimPlot(seu, reduction = "tsne", pt.size = 0.1, label = TRUE, 
              label.size = 3, cols = colCls) + plotTheme + coord_fixed()
p2 <- DimPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE, 
              label.size = 3, cols = colCls) + plotTheme + coord_fixed()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 4, filename = "images/clustDimPlot.png")

# Proportion / cell number composition per library
ggData <- as.matrix(prop.table(table(seu$cluster, seu$library), margin = 2))
pheatmap(ggData, color = colorRampPalette(colGEX)(100), 
         cluster_rows = FALSE, cutree_cols = 2, 
         display_numbers = TRUE, number_format = "%.3f", angle_col = 315,
         width = 4, height = 6, filename = "images/clustComLibH.png")

ggData = data.frame(prop.table(table(seu$cluster, seu$library), margin = 2))
colnames(ggData) = c("cluster", "library", "value")
p1 <- ggplot(ggData, aes(library, value, fill = cluster)) +
  geom_col() + xlab("Library") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colCls) + plotTheme + coord_flip()
ggData = data.frame(table(seu$cluster, seu$library))
colnames(ggData) = c("cluster", "library", "value")
p2 <- ggplot(ggData, aes(library, value, fill = cluster)) +
  geom_col() + xlab("Library") + ylab("Cell Number") +
  scale_fill_manual(values = colCls) + plotTheme + coord_flip()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 6, filename = "images/clustComLib.png")

# Proportion / cell number composition per cluster
ggData = data.frame(prop.table(table(seu$library, seu$cluster), margin = 2))
colnames(ggData) = c("library", "cluster", "value")
p1 <- ggplot(ggData, aes(cluster, value, fill = library)) +
  geom_col() + xlab("Cluster") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colLib) + plotTheme + coord_flip()
ggData = data.frame(table(seu$library, seu$cluster))
colnames(ggData) = c("library", "cluster", "value")
p2 <- ggplot(ggData, aes(cluster, value, fill = library)) +
  geom_col() + xlab("Cluster") + ylab("Cell Number") +
  scale_fill_manual(values = colLib) + plotTheme + coord_flip()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 6, filename = "images/clustComClust.png")

### B. Marker genes
# Find Markers
oupMarker <- FindAllMarkers(seu, only.pos = TRUE, 
                            logfc.threshold = 1.0, min.pct = 0.2)
oupMarker <- data.table(oupMarker)
oupMarker$pct.diff = oupMarker$pct.1 - oupMarker$pct.2
oupMarker <- oupMarker[, c("cluster","gene","avg_log2FC","pct.1","pct.2",
                           "pct.diff","p_val","p_val_adj")]
fwrite(oupMarker, sep = "\t", file = "images/clustMarkers.txt")
seu@misc$marker <- oupMarker      # Store markers into Seurat object

# Check if known genes are in the marker gene list
knownGenes <- c("CD34","CRHBP","GATA1",  "CD14","IRF8","CD19",
                "CD4","CD8B","GNLY")
oupMarker[gene %in% knownGenes]

# Get top genes for each cluster and do dotplot / violin plot
oupMarker$cluster = factor(oupMarker$cluster, levels = reorderCluster)
oupMarker = oupMarker[order(cluster, -avg_log2FC)]
genes.to.plot <- unique(oupMarker[cluster %in% reorderCluster, 
                                  head(.SD, 2), by="cluster"]$gene)
p1 <- DotPlot(seu, group.by = "cluster", features = genes.to.plot) + 
  coord_flip() + scale_color_gradientn(colors = colGEX) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))
ggsave(p1, width = 10, height = 8, filename = "images/clustMarkersDot.png")
p2 <- VlnPlot(seu, group.by = "cluster", fill.by = "ident", cols = colCls,
              features = genes.to.plot, stack = TRUE, flip = TRUE)
ggsave(p2, width = 10, height = 8, filename = "images/clustMarkersVln.png")

### C. Gene module analysis + module score
# Functional analysis using clusterProfiler and msigdb
oupMarkerFunc = data.table()
set.seed(42)
for(iDB in c("C8", "C5_GO:BP")){
  # Get reference gene sets from msigdbr
  inpGS <- tstrsplit(iDB, "_")
  msigCat <- inpGS[[1]]; msigSubCat <- NULL
  if(length(inpGS) >= 2){msigSubCat <- inpGS[[2]]}
  inpGS <- data.frame(msigdbr(species = "Homo sapiens", category = msigCat,
                              subcategory = msigSubCat))
  inpGS <- inpGS[, c("gs_name","gene_symbol")]
  
  # Start up clusterProfiler
  for(i in unique(oupMarker$cluster)){
    tmpOut <- enricher(oupMarker[cluster == i]$gene,
                       universe = rownames(seu), TERM2GENE = inpGS)
    tmpOut <- data.table(sigdb = iDB, cluster = i, data.frame(tmpOut[, -2]))
    tmpOut$mLog10Padj <- -log10(tmpOut$p.adjust)
    tmpOut <- tmpOut[order(-mLog10Padj)]
    oupMarkerFunc <- rbindlist(list(oupMarkerFunc, tmpOut))
  }
}
oupMarkerFunc$cluster <- factor(oupMarkerFunc$cluster, levels = reorderCluster)
seu@misc$markerFunc <- oupMarkerFunc   # Store func analysis into Seurat object

# Plot functional analysis results
ggData <- oupMarkerFunc[grep("BONE_MARROW", ID)]
ggData$ID <- gsub("HAY_BONE_MARROW_", "", ggData$ID)
ggData$ID <- factor(ggData$ID, levels = unique(ggData$ID))
p1 <- ggplot(ggData, aes(cluster, ID, size = Count, color = mLog10Padj)) + 
  geom_point() + ggtitle("HAY_BONE_MARROW signatures") + 
  theme_linedraw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  scale_color_gradientn(colors = colGEX, limits = c(0,10), na.value = colGEX[8])
ggsave(p1, width = 12, height = 8, filename = "images/clustMarkersFuncBM.png")

tmp <- oupMarkerFunc[sigdb == "C5_GO:BP"]
tmp <- tmp[grep("DIFF", ID)][!grep("POSITIVE|NEGATIVE", ID)]$ID
ggData <- oupMarkerFunc[ID %in% tmp]
ggData$ID <- substr(gsub("GOBP_", "", ggData$ID), 1, 30)
ggData$ID <- factor(ggData$ID, levels = unique(ggData$ID))
p1 <- ggplot(ggData, aes(cluster, ID, size = Count, color = mLog10Padj)) + 
  geom_point() + ggtitle("DIFF signatures") + 
  theme_linedraw(base_size = 18) + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
  scale_color_gradientn(colors = colGEX, limits = c(0,8), na.value = colGEX[8])
ggsave(p1, width = 12, height = 10, filename = "images/clustMarkersFuncDiff.png")

# Cell cycle + plot on UMAP
seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes, 
                        g2m.features = cc.genes$g2m.genes)
p1 <- DimPlot(seu, reduction = "umap", pt.size = 0.1, group.by = "cluster", 
              shuffle = TRUE, cols = colCls) + plotTheme + coord_fixed()
p2 <- DimPlot(seu, reduction = "umap", pt.size = 0.1, group.by = "Phase", 
              shuffle = TRUE, cols = colCcy) + plotTheme + coord_fixed()
p3 <- FeaturePlot(seu, pt.size = 0.1, feature = "S.Score") + 
  scale_color_distiller(palette = "RdYlBu") + plotTheme + coord_fixed()
p4 <- FeaturePlot(seu, pt.size = 0.1, feature = "G2M.Score") + 
  scale_color_distiller(palette = "RdYlBu") + plotTheme + coord_fixed()
ggsave(p1 + p2 + p3 + p4, 
       width = 10, height = 8, filename = "images/clustModuCellCycle.png")

# Add module score (Here, we use using HALLMARK SIGNALING gene sets)
inpGS <- data.table(msigdbr(species = "Homo sapiens", category = "H"))
inpGS <- inpGS[grep("SIGNAL", gs_name)]
inpGS$gs_name <- gsub("HALLMARK_", "", inpGS$gs_name)
inpGS$gs_name <- gsub("_SIGNALING", "", inpGS$gs_name)
inpGS <- split(inpGS$gene_symbol, inpGS$gs_name)
seu <- AddModuleScore(seu, features = inpGS, name = "HALLMARK")
colnames(seu@meta.data)[grep("HALLMARK", colnames(seu@meta.data))] <- names(inpGS)

p1 <- FeaturePlot(seu, reduction = "umap", pt.size = 0.1,
                  features = names(inpGS), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & plotTheme & coord_fixed()
ggsave(p1, width = 20, height = 12, filename = "images/clustModuScore.png")

### D. Annotating clusters
oupAnnot = data.table(cluster = reorderCluster)
# Add in top 5 marker genes
ggData <- oupMarker[, head(.SD, 8), by = "cluster"]
ggData <- ggData[, paste0(gene, collapse = ","), by = "cluster"]
colnames(ggData)[2] <- "markers"
oupAnnot <- ggData[oupAnnot, on = "cluster"]
# Add in BONE_MARROW enrichment
ggData <- oupMarkerFunc[grep("BONE_MARROW", ID)][, head(.SD, 1), by = "cluster"]
ggData <- ggData[, c("cluster", "ID")]
ggData$ID <- gsub("HAY_BONE_MARROW_", "", ggData$ID)
colnames(ggData)[2] <- "HAY_BONE_MARROW"
oupAnnot <- ggData[oupAnnot, on = "cluster"]
# Add in annotation
oupAnnot$annotation <- c("HSC","MPP","MyP1","MyP2","LyP1","LyP2","ProB",
                         "ERP","Ery","MKP","EOBM",
                         "pDC1","pDC2","cDC","Monocyte","InflamMono","Stroma",
                         "PreB","MemoryB","MatureB","BT","PlasmaCell",
                         "CD8+Naive","CD4+Naive","CD4+TCM",
                         "CD8+GZMK1","CD8+GZMK2","CD8+TEMRA","NK")
# Output table
png("images/clustReAnnotTable.png", 
    width = 12, height = 9, units = "in", res = 300)
p1 <- tableGrob(oupAnnot, rows = NULL)
grid.arrange(p1)
dev.off()

# Add new annotation into seurat and replot
tmpMap <- oupAnnot$annotation
names(tmpMap) <- oupAnnot$cluster
seu$celltype <- tmpMap[as.character(seu$cluster)]
seu$celltype <- factor(seu$celltype, levels = tmpMap)
Idents(seu) <- seu$celltype  # Set seurat to use celltype

# Plot ideal resolution on tSNE and UMAP
p1 <- DimPlot(seu, reduction = "tsne", pt.size = 0.1, label = TRUE, 
              label.size = 3, cols = colCls) + plotTheme + coord_fixed()
p2 <- DimPlot(seu, reduction = "umap", pt.size = 0.1, label = TRUE, 
              label.size = 3, cols = colCls) + plotTheme + coord_fixed()
ggsave(p1 + p2 & theme(legend.position = "none"), 
       width = 8, height = 4, filename = "images/clustReAnnotTsUm.png")



# Save Seurat Object at end of each section
saveRDS(seu, file = "bmSeu.rds")
