#####################################################
# Aim: Analyse 10x multiome (snRNA- and snATAC-seq) #
# to define tdTomato-positive cells at E9.5.        #
#####################################################

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)

e9pos <- readRDS("e9-5_P3CreTracedMultiomeSeuratObject.rds")

e9pos <- subset(
  x = e9pos,
  subset = nCount_ATAC <75000 & TSS.enrichment > 1 & nucleosome_signal <1.5 & TSS.enrichment < 10 &
    nCount_ATAC > 1e3 &
    nCount_RNA < 20000 &
    nCount_RNA > 1000 &
    percent.mt < 15
)

# e9pos RNA analysis
DefaultAssay(e9pos) <- "RNA"
e9pos <- NormalizeData(e9pos)
e9pos <- FindVariableFeatures(e9pos)
e9pos <- ScaleData(e9pos)
e9pos <- RunPCA(e9pos, npcs = 50)
ElbowPlot(e9pos, ndims = 50)

e9pos <- FindNeighbors(e9pos, dims = 1:35, reduction = "pca")
e9pos <- FindClusters(e9pos)
e9pos <- RunUMAP(e9pos, dims = 1:35, reduction = "pca")

DimPlot(e9pos, reduction = "umap")

# e9pos ATAC analysis
DefaultAssay(e9pos) <- "ATAC"
e9pos <- RunTFIDF(e9pos)
e9pos <- FindTopFeatures(e9pos, min.cutoff = 'q0')
e9pos <- RunSVD(e9pos)
e9pos <- RunUMAP(e9pos, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#Multimodal analysis
e9pos <- FindMultiModalNeighbors(e9pos, reduction.list = list("pca", "lsi"), dims.list = list(1:35, 2:50))
e9pos <- RunUMAP(e9pos, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
e9pos <- FindClusters(e9pos, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.8)

# Display results
p1 <- DimPlot(e9pos, reduction = "umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(e9pos, reduction = "umap.atac", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(e9pos, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

Idents(e9pos) <- "wsnn_res.0.8"
e9pos <- RenameIdents(e9pos,
                      '4' = "cervical somite",
                      "2" = "late PXM angioblast",
                      "0" = "late PXM angioblast",
                      "3" = "venous EC",
                      "9" = "angiogenic venous",
                      "1" = "PXM EC",
                      "7" = "trunk somite",
                      "5" = "dermomyotome",
                      "6" = "early PXM angioblast",
                      "8"= "angiogenic venous",
                      "13"= "sinus venosus",
                      "10" = "mixed mesoderm",
                      "12" = "neural crest",
                      "11" = "neural")

my_cols_preF <- c("Sprouting" = '#8FD8D3',"venous EC" = '#D88B91',"angiogenic venous"= "#65B4C4","neural tube EC" = '#FDE360',
                  "cervical somite" = '#A0A0A0',"early PXM EC 1" = '#836E96',"Pre-Arterial" = '#98C38F',"dermomyotome" = '#D8C5DE', "trunk somite" = '#CECECE',        
                  "early PXM angioblast" = "#BEA8CD","Angioblast2_PXM" = '#AB99BA',"late PXM angioblast" = '#AB99BA', "Angioblast2_G1S" = '#AB99BA',
                  "sinus venosus" = '#C97373', "Sclerotome"  = '#4527A0',"PXM EC" ="#EDA337", "neural" = '#98C38F', "neural crest" = "#74b567", "mixed mesoderm" = "#d1cdd4")

DimPlot(e9pos, reduction = "wnn.umap", label = TRUE, cols = my_cols_preF, pt.size = 1.2) & NoLegend() & NoAxes()

# Run doublet detection and thresholding ####
suppressPackageStartupMessages(library(scDblFinder))
library(SingleCellExperiment)
library(patchwork)

# Calculate RNA doublet score
DefaultAssay(e9pos) <- "RNA"
e9sce <- as.SingleCellExperiment(e9pos)
e9sce <- scDblFinder(e9sce, dbr=0.077)
table(e9sce$scDblFinder.class)
e9pos$RNAdoubletScore <- e9sce$scDblFinder.score
p1 <- FeaturePlot(e9pos, features = c('RNAdoubletScore'),reduction = 'wnn.umap',cols = c('lightgrey','darkred'))

# Calculate ATAC doublet score
DefaultAssay(e9pos) <- 'ATAC'
e9ATACsce <- as.SingleCellExperiment(e9pos)
e9ATACsce <- scDblFinder(e9ATACsce, artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
e9pos$ATACdoubletScore <- e9ATACsce$scDblFinder.score
p2 <- FeaturePlot(e9pos, features = c('ATACdoubletScore'),reduction = 'wnn.umap',cols = c('lightgrey','darkgreen'))

#Calculate multiome doublet score
e9pos$MultiomeDoubletScore <- e9pos$RNAdoubletScore + e9pos$ATACdoubletScore
p3 <- FeaturePlot(e9pos, features = c('MultiomeDoubletScore'),reduction = 'wnn.umap',cols = c('lightgrey','purple'))
wrap_plots(p1,p2,p3, ncol = 3) & NoAxes()

# Filter out doublets based on an estimated doublet rate of 7.7%
4116*0.077 #calculate estimated number of doublets
sum(e9pos$MultiomeDoubletScore > 1.1)
e9pos$MultiomeDoubletClass <- e9pos$MultiomeDoubletScore > 1.1
DimPlot(e9pos, group.by = 'MultiomeDoubletClass', pt.size = 1.3, reduction = 'wnn.umap')

e9pos_fil <- e9pos[,e9pos$MultiomeDoubletClass == FALSE]
e9pos_fil

# Display QC violin plots for doublet filtered data
VlnPlot(
  object = e9pos_fil,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1,
  cols = "#8F589F"
) & xlab(label = NULL)

VlnPlot(
  object = e9pos_fil,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 3,
  pt.size = 0.1,
  cols = "#8F589F"
) & xlab(label = NULL)


# Subset for cells in the PXM trajectory and recluster ####
e9pos_fil_select <- subset(e9pos_fil, idents = c("neural","neural crest","mixed mesoderm"), invert = T)

# e9pos_fil_select RNA analysis
DefaultAssay(e9pos_fil_select) <- "RNA"
e9pos_fil_select <- NormalizeData(e9pos_fil_select)
e9pos_fil_select <- FindVariableFeatures(e9pos_fil_select)
e9pos_fil_select <- ScaleData(e9pos_fil_select)
e9pos_fil_select <- RunPCA(e9pos_fil_select, npcs = 50)
ElbowPlot(e9pos_fil_select, ndims = 50)

e9pos_fil_select <- FindNeighbors(e9pos_fil_select, dims = 1:25, reduction = "pca")
e9pos_fil_select <- FindClusters(e9pos_fil_select)
e9pos_fil_select <- RunUMAP(e9pos_fil_select, dims = 1:25, reduction = "pca")

DimPlot(e9pos_fil_select, reduction = "umap")

# e9pos_fil_select ATAC analysis
DefaultAssay(e9pos_fil_select) <- "ATAC"
e9pos_fil_select <- RunTFIDF(e9pos_fil_select)
e9pos_fil_select <- FindTopFeatures(e9pos_fil_select, min.cutoff = 'q0')
e9pos_fil_select <- RunSVD(e9pos_fil_select)
e9pos_fil_select <- RunUMAP(e9pos_fil_select, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#Multimodal analysis
e9pos_fil_select <- FindMultiModalNeighbors(e9pos_fil_select, reduction.list = list("pca", "lsi"), dims.list = list(1:25, 2:50))
e9pos_fil_select <- RunUMAP(e9pos_fil_select, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
e9pos_fil_select <- FindClusters(e9pos_fil_select, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 1)

# Display results
p1 <- DimPlot(e9pos_fil_select, reduction = "umap", group.by = "wsnn_res.1", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(e9pos_fil_select, reduction = "umap.atac", group.by = "wsnn_res.1", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(e9pos_fil_select, reduction = "wnn.umap", group.by = "wsnn_res.1", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

Idents(e9pos_fil_select) <- "wsnn_res.1"

e9pos_fil_select <- RenameIdents(e9pos_fil_select,
                                 '0' = "cervical somite",
                                 "1" = "late angioblast",
                                 "8" = "late angioblast",
                                 "3" = "venous EC",
                                 "2" = "early EC 1",
                                 "5" = "neural tube EC",
                                 "4" = "trunk somite",
                                 "6" = "dermomyotome",
                                 "9" = "early angioblast",
                                 "11"= "angiogenic venous",
                                 "12"= "sinus venosus",
                                 "10" = "angiogenic venous",
                                 "7" = "early EC 2")

e9pos_fil_select$celltypeNew <- Idents(e9pos_fil_select)


my_cols_preF <- c("Sprouting" = '#8FD8D3',"venous EC" = '#D88B91',"angiogenic venous"= "#65B4C4","neural tube EC" = '#FDE360',
                  "cervical somite" = '#A0A0A0',"early EC 1" = '#836E96',"Pre-Arterial" = '#98C38F',"dermomyotome" = '#D8C5DE', "trunk somite" = '#CECECE',        
                  "early angioblast" = "#BEA8CD","Angioblast2_PXM" = '#AB99BA',"late angioblast" = '#AB99BA', "Angioblast2_G1S" = '#AB99BA',
                  "sinus venosus" = '#C97373', "Sclerotome"  = '#4527A0',"early EC 2" ="#EDA337", "neural" = '#98C38F', "neural crest" = "#74b567", "mixed mesoderm" = "#d1cdd4")

DimPlot(e9pos_fil_select, reduction = "wnn.umap", group.by = "celltypeNew", label = T, label.size = 2.5, cols = my_cols_preF, pt.size = 0.75) &NoLegend() & NoAxes()


# Identifying top marker genes by cell type ####
# Reorder cluster according to approximate differentiation order
my_levels <- c("cervical somite",
               "trunk somite",
               "dermomyotome",
               "early PXM angioblast",
               "late PXM angioblast",
               "early PXM EC 1",
               "early PXM EC 2",
               "venous EC",
               "angiogenic venous",
               "neural tube EC",
               "sinus venosus")

e9pos_fil_select@active.ident <- factor(x = e9pos_fil_select@active.ident, levels = my_levels)

# Find top markers for each cell type
DefaultAssay(e9pos_fil_select) <- "RNA"
markers <- FindAllMarkers(e9pos_fil_select, assay = "RNA", min.pct = 0.3)
markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) -> top2

# Order genes to match cluster order
top2 %>% arrange(factor(cluster, levels = my_levels)) -> top2

# Generate heatmap of top 2 marker genes for each celltype
DoHeatmap(e9pos_fil_select, features = top2$gene, assay = "RNA",group.colors = my_cols_preF, size = 3,lines.width = 10, raster = F) + 
  scale_fill_gradientn(colours = rev(mapal),values = scales::rescale(c(-0.3, 0, 1, 2)),aesthetics = "fill") &
  NoLegend()

# Generate figure showing motif accessibility ####
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(e9) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(e9), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
e9 <- SetAssayData(e9, assay = 'ATAC', slot = 'motifs', new.data = motif.object)
e9 <- RunChromVAR(
  object = e9,
  genome = (BSgenome.Mmusculus.UCSC.mm10)
)
motif.name <- ConvertMotifID(linPos, name = c('ETV2','NR2F2','SOX18','GATA2','FOXC2','KLF4'))
DefaultAssay(linPos) <-"chromvar"
FeaturePlot(linPos, features = motif.name[1], min.cutoff = 0,
            cols=c("#E0E0E0","darkred"), pt.size = 0.75, order=T , reduction = 'wnn.umap',
            ncol = 1) & NoLegend()

# Save RDS file ####
saveRDS(e9pos_fil_select, file = "e9pos_doubletfilteredPAPER.rds")





