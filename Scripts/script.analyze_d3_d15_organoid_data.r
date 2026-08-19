library(Seurat)
library(dplyr)
library(simspec)
source("~/Tools/scripts/cell_density_diff.r")
source("~/Tools/scripts/feature_plots.r")

seurat <- readRDS("/links/groups/treutlein/USERS/Andrea/10x/seurat_objects/integrated/seurat_organoids_integrated.rds")
seurat[['seurat_pca']] <- seurat[['pca']]
seurat[['umap_seurat']] <- seurat[['umap']]

DefaultAssay(seurat) <- "RNA"
seurat <- subset(seurat, subset=nFeature_RNA > 1500)
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
blacklist <- c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt", stringsAsFactors=F)[,1])
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
seurat <- CellCycleScoring(seurat, g2m.features=cc.genes.updated.2019$g2m.genes, s.features=cc.genes.updated.2019$s.genes, set.ident=F)
seurat <- ScaleData(seurat, vars.to.regress=c("G2M.Score","S.Score")) %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)
UMAPPlot(seurat, group.by="orig.ident")

seurat <- cluster_sim_spectrum(seurat, label_tag="orig.ident")
seurat <- RunUMAP(seurat, reduction="css", dims = 1:ncol(Embeddings(seurat,"css")), reduction.name="umap_css", reduction.key="UMAPCSS_")

par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat,"umap_css"), seurat$orig.ident, colorPal=prettyrainbow_colscheme, cex=0.6, pt_border=T)
legend("topleft", bty="n", pch=16, col=prettyrainbow_colscheme(4), legend=levels(factor(seurat$orig.ident)), cex=1.3)
FeaturePlot(seurat, features=c("SOX2","DCX","AIF1","TTR","HOXB2","OTX2","FOXG1","EMX1","MKI67"), cols=beach_colscheme(30), reduction="umap_css", order=T) & NoLegend() & NoAxes()
FeaturePlot(seurat, c("DCN","COL1A2","COL6A2","MFAP4"), reduction="umap_css", cols=beach_colscheme(30)) & NoLegend() & NoAxes()

seurat <- FindNeighbors(seurat, reduction="css", dims = 1:ncol(Embeddings(seurat,"css"))) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1)
DimPlot(seurat, group.by="RNA_snn_res.0.1", reduction="umap_css", label=T) & NoLegend() & NoAxes()
DimPlot(seurat, group.by="RNA_snn_res.0.5", reduction="umap_css", label=T) & NoLegend() & NoAxes()
DimPlot(seurat, group.by="RNA_snn_res.1", reduction="umap_css", label=T) & NoLegend() & NoAxes()

seurat$annotation <- setNames(c("Rhombencephalic NPC","Diencephalic/Mesencephalic NPC", "Rhombencephalic neuron", "Diencephalic/Mesencephalic neuron","Mesenchyme","Dorsal telencephalic NPC"), levels(seurat$RNA_snn_res.0.1))[seurat$RNA_snn_res.0.1]
seurat$annotation[seurat$RNA_snn_res.5=="57"] <- "Dorsal telencephalic neuron"

par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat,"umap_css"), seurat$annotation, colorPal=prettyrainbow_colscheme, cex=0.5, pt_border=T)
legend("topleft", pch=21, cex=1.1, pt.cex=1.5, col="#303030", pt.bg=prettyrainbow_colscheme(length(unique(seurat$annotation))), legend=sort(unique(seurat$annotation)), bty="n")
plotFeature(Embeddings(seurat,"umap_seurat"), seurat$annotation, colorPal=prettyrainbow_colscheme, cex=0.5, pt_border=T)
legend("topright", pch=21, cex=1.1, pt.cex=1.5, col="#303030", pt.bg=prettyrainbow_colscheme(length(unique(seurat$annotation))), legend=sort(unique(seurat$annotation)), bty="n")

saveRDS(seurat, file="data.organoid_integrated.seurat.rds")

knn <- RANN::nn2(Embeddings(seurat, "css"), k=50)
enrichment_time <- get_knn_fisher_enrichment(knn = knn, bool_case = seurat$orig.ident %in% c("Day 15 Organoid Co-Culture","Day 15 Organoid Control"), bool_bg = seurat$orig.ident %in% c("Day 3 Organoid Co-Culture","Day 3 Organoid Control"))
enrichment_mg <- get_knn_fisher_enrichment(knn = knn, bool_case = seurat$orig.ident %in% c("Day 15 Organoid Co-Culture","Day 3 Organoid Co-Culture"), bool_bg = seurat$orig.ident %in% c("Day 15 Organoid Control","Day 3 Organoid Control"))
enrichment_samples <- sapply(sort(unique(seurat$orig.ident)), function(x) get_knn_fisher_enrichment(knn=knn, bool_case = seurat$orig.ident == x, bool_bg=seurat$orig.ident != x) )
plotMultiFeatures(Embeddings(seurat,"umap_css"), t(enrichment_samples), colorPal=bluered_colscheme, cex=0.5, pt_border=T, cex.main=1.5)


