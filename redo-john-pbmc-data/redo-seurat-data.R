pbmc.data <- Read10X(data.dir = "seurat-tut-data/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
inpUMIs <- pbmc[["RNA"]]@counts

oupQCcell <- data.table(
  sampleID = colnames(inpUMIs),
  library = tstrsplit(colnames(inpUMIs), "_")[[1]],
  nUMI = colSums(inpUMIs),
  nGene = colSums(inpUMIs != 0),
  nMT = colSums(inpUMIs[grep("^MT-", rownames(inpUMIs)),]))
oupQCcell$pctMT <- 100 * oupQCcell$nMT / oupQCcell$nUMI
oupQCcell$library <- factor(oupQCcell$library, levels = names(colLib))


low.n.UMI <- 0.7e3
high.n.UMI <- 6e3
low.n.genes <- 0.3e3
high.MT.pct <- 4.5

p1 <- ggplot(oupQCcell, aes(nUMI, fill = library)) + 
  geom_histogram(binwidth = 500, color = "black") + xlim(c(0,40e3)) + 
  geom_vline(xintercept = c(low.n.UMI,high.n.UMI), color = "red", linetype = "dashed")+ 
  xlab("No. UMIs") + scale_fill_manual(values = colLib) + plotTheme


p2 <- ggplot(oupQCcell, aes(nGene, fill = library)) + 
  geom_histogram(binwidth = 100, color = "black") + xlim(c(0,6e3)) + 
  geom_vline(xintercept = c(low.n.genes), color = "red", linetype = "dashed")+ 
  xlab("No. Detected Genes") + scale_fill_manual(values = colLib) + plotTheme


p3 <- ggplot(oupQCcell, aes(pctMT, fill = library)) + 
  geom_histogram(binwidth = 0.5, color = "black") + xlim(c(0,20)) + 
  geom_vline(xintercept = c(high.MT.pct), color = "red", linetype = "dashed")+ 
  xlab("Percentage MT Genes") + scale_fill_manual(values = colLib) + plotTheme

ggsave(p1 + p2 + p3 + guide_area() + plot_layout(nrow = 2, guides = "collect"), 
       width = 10, height = 8, filename = "images/xxx2basicCellQC.png")


oupQCcell <- oupQCcell[nUMI > low.n.UMI]
oupQCcell <- oupQCcell[nUMI < high.n.UMI]
oupQCcell <- oupQCcell[nGene > low.n.genes]
oupQCcell <- oupQCcell[pctMT < high.MT.pct]
inpUMIs <- inpUMIs[, as.character(oupQCcell$sampleID)]

# Gene QC
inpUMIs <- inpUMIs[rowSums(inpUMIs) > 0,]   # Remove genes with no reads
oupQCgene <- data.table(
  gene = rownames(inpUMIs),
  meanRead = rowSums(inpUMIs) / ncol(inpUMIs),
  cellExpr = rowSums(inpUMIs > 0))
oupQCgene$log10meanRead <- log10(oupQCgene$meanRead)
oupQCgene$log10cellExpr <- log10(oupQCgene$cellExpr)

low.n.cells <- 13
# Plot gene QC metrics
p1 <- ggplot(oupQCgene, aes(log10meanRead)) + 
  geom_histogram(binwidth = 0.05, color = "black") + 
  # geom_vline(xintercept = c(-3), color = "red", linetype = "dashed")+ 
  xlab("log10(Average UMIs)") + plotTheme
p2 <- ggplot(oupQCgene[cellExpr != 0], aes(log10cellExpr)) + 
  geom_histogram(binwidth = 0.05, color = "black") + 
  geom_vline(xintercept = log10(low.n.cells), color = "red", linetype = "dashed")+ 
  xlab("log10(No. Cells Expressing Gene)") + plotTheme
ggsave(p1 + p2, width = 10, height = 4, filename = "images/xxx2basicGeneQC.png")

oupQCgene <- oupQCgene[cellExpr >= low.n.cells]
inpUMIs <- inpUMIs[as.character(oupQCgene$gene), ]

## Create seu object

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
p2 <- LabelPoints(plot = p2, repel = TRUE, points = VariableFeatures(seu)[1:15])
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 4, filename = "images/xxxbasicHVG.png")


## PCA

seu <- ScaleData(object = seu)    # Scale data prior to PCA
seu <- RunPCA(seu)

p1 <- DimPlot(seu, reduction = "pca", pt.size = 0.1, shuffle = TRUE, 
              cols = colLib) + plotTheme + coord_fixed()
p2 <- DimPlot(seu, reduction = "pca", pt.size = 0.1, shuffle = TRUE,
              cols = colLib, dims = c(1,3)) + plotTheme + coord_fixed()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 4, filename = "images/xxxbasicPCA.png")


# Elbow plot
p1 <- ElbowPlot(seu, ndims = 40)
ggData <- seu@reductions[["pca"]]@stdev  # Fit lines for elbow plot
lm(y~x, data = data.frame(x = 16:20, y = ggData[16:20]))
lm(y~x, data = data.frame(x = 31:35, y = ggData[31:35]))
nPC <- 9 # determined from elbow plot (below, after it was run the first time)

p1 <- p1 + plotTheme + geom_vline(xintercept = nPC, color = "firebrick") 
#+ 
#  geom_path(data = data.frame(x = 16:35, y = 16:35*-0.07535 + 3.50590),
#            aes(x,y), color = "dodgerblue") +
# geom_path(data = data.frame(x = 16:35, y = 16:35*-0.0282 + 2.3996),
#            aes(x,y), color = "dodgerblue")

ggsave(p1, width = 5, height = 4, filename = "images/xxxbasicPCAelbow.png")


# Run tSNE
seu <- RunTSNE(seu, dims = 1:nPC, num_threads = nC)
p1 <- DimPlot(seu, reduction = "tsne", pt.size = 0.1, shuffle = TRUE, 
              cols = colLib) + plotTheme + coord_fixed()

# Run UMAP
seu <- RunUMAP(seu, dims = 1:nPC)
p2 <- DimPlot(seu, reduction = "umap", pt.size = 0.1, shuffle = TRUE, 
              cols = colLib) + plotTheme + coord_fixed()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 4, filename = "images/xxxbasicTsUm.png")


p1 <- FeaturePlot(seu, reduction = "umap", pt.size = 0.1,
                  features = c("CD34","CRHBP","GATA1",
                               "CD14","IRF8","CD19",
                               "CD4","CD8B","GZMB"), order = TRUE) &
  scale_color_gradientn(colors = colGEX) & plotTheme & coord_fixed()
ggsave(p1, width = 15, height = 12, filename = "images/xxxbasicUmapGex.png")

saveRDS(seu, file = "pbmc-Seu.rds")

### A. Clustering
# Perform clustering
seu <- FindNeighbors(seu, dims = 1:nPC)
seu <- FindClusters(seu, resolution = seq(0.4,1.5,0.1), verbose = FALSE)
library(clustree)
library(RColorBrewer)
# Find ideal resolution
p1 <- clustree(seu)
ggsave(p1 + scale_color_manual(values = brewer.pal(12, "Paired")) + 
         guides(color = guide_legend(ncol = 2)), 
       width = 10, height = 8, filename = "images/xxxclustClustree.png")


optRes <- 1.1                # determined from clustering tree
# seu$seurat_clusters <- NULL  # Remove this column to prevent confusion
seu$cluster <- seu[[paste0("RNA_snn_res.", optRes)]] 


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
       width = 10, height = 4, filename = "images/xxxclustDimPlot.png")


# Proportion / cell number composition per cluster

### Marker genes

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
ggsave(p1, width = 10, height = 8, filename = "images/xxxclustMarkersDot.png")


p2 <- VlnPlot(seu, group.by = "cluster", fill.by = "ident", cols = colCls,
              features = genes.to.plot, stack = TRUE, flip = TRUE)
ggsave(p2, width = 10, height = 8, filename = "images/xxxxclustMarkersVln.png")