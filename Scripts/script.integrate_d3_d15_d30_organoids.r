library(Seurat)
library(Matrix)
library(simspec)
library(dplyr)
library(qlcMatrix)

library(SingleCellExperiment)
library(celda)
library(SoupX)

library(presto)

source("~/Tools/scripts/feature_plots.r")
source("~/Tools/scripts/py_util.r")

# load d3 and d15 data
seurat_d3_d15 <- readRDS("data.organoid_d3_d15_integrated.seurat.rds")

# load d30 data
counts_d30_cocul <- Read10X_h5("rawdata/day30_org_cocult/outs/filtered_feature_bc_matrix.h5")
counts_d30_ctrl <- Read10X_h5("rawdata/day30_org_ctrl/outs/filtered_feature_bc_matrix.h5")
seurat_d30_cocul <- CreateSeuratObject(counts = counts_d30_cocul, project="Day 30 Organoid Co-Culture")
seurat_d30_ctrl <- CreateSeuratObject(counts = counts_d30_ctrl, project="Day 30 Organoid Control")
seurat_d30 <- merge(seurat_d30_cocul, list(seurat_d30_ctrl))
seurat_d30 <- PercentageFeatureSet(seurat_d30, pattern="^MT-", col.name="percent.mt")
seurat_d30 <- NormalizeData(seurat_d30)
seurat_d30 <- subset(seurat_d30, subset = nFeature_RNA > 1500 & nFeature_RNA < 7500 & percent.mt < 12)
seurat_d30 <- CellCycleScoring(seurat_d30, g2m.features=cc.genes.updated.2019$g2m.genes, s.features=cc.genes.updated.2019$s.genes, set.ident=F)

# merge and integrate
seurat <- merge(seurat_d3_d15, list(seurat_d30))
seurat@active.ident <- setNames(factor(seurat$orig.ident), colnames(seurat))
VlnPlot(seurat, "nFeature_RNA", pt.size=0, cols = setNames(prettyrainbow_colscheme(6), sort(unique(seurat$orig.ident)))) & NoLegend()

seurat <- FindVariableFeatures(seurat, nfeatures = 5000)
blacklist <- c(unlist(cc.genes.updated.2019), read.table("~/Work/databases/GeneLists/RPgenes.txt", stringsAsFactors=F)[,1])
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
seurat <- ScaleData(seurat, vars.to.regress=c("G2M.Score","S.Score"))
seurat <- RunPCA(seurat, npcs = 50) %>% RunUMAP(dims = 1:20)
seurat <- cluster_sim_spectrum(seurat, label_tag="orig.ident", cluster_resolution=1)
seurat <- RunUMAP(seurat, reduction="css", dims = 1:ncol(Embeddings(seurat,"css")), reduction.name="umap_css", reduction.key="UMAPCSS_")
saveRDS(seurat, file="data.organoid_d3_d15_d30_integrated.seurat.rds")

(FeaturePlot(seurat, reduction="umap_css", features=c("DCN","SOX2","DCX","FOXG1","OTX2","HOXB2","LHX9","AIF1","MKI67"), cols=beach_colscheme(30), ncol=3, order=T, pt.size=0.01) & NoAxes() & NoLegend()) | (DimPlot(seurat, reduction="umap_css", group.by="orig.ident") & NoAxes())

par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat,"umap_css"), seurat$orig.ident, colorPal=prettyrainbow_colscheme, cex=0.5, pt_border=T)
legend("topleft", bty="n", pch=21, pt.cex=1.5, col="#303030", pt.bg=prettyrainbow_colscheme(length(unique(seurat$orig.ident))), legend=sort(unique(seurat$orig.ident)))
plotMultiFeatures(Embeddings(seurat,"umap_css"), seurat@assays$RNA@data[c("DCN","SOX2","DCX","FOXG1","OTX2","HOXB2","LHX9","AIF1","MKI67"),], cex=0.4, colorPal=beach_colscheme, random_order=F, sort_by_value=T, cex.main=1.5)





# ambigent RNA cleanup of d30-coculture organoid: primary cleanup - remove ambient RNA droplets
## NMF on normalized log-scaled data (single component)
nmf_d30cocul <- run_NMF(t(seurat@assays$RNA@data[VariableFeatures(seurat),which(seurat$orig.ident == "Day 30 Organoid Co-Culture")]), n_components=1, max_iter = 5000)
colnames(nmf_d30cocul$H) <- VariableFeatures(seurat)
rownames(nmf_d30cocul$W) <- colnames(seurat)[which(seurat$orig.ident == "Day 30 Organoid Co-Culture")]
saveRDS(nmf_d30cocul, file="res.nmf_d30_coculture_organoid_1comp.rds")

corr_d30cocul_nmf1 <- corSparse(seurat@assays$RNA@data[VariableFeatures(seurat),which(seurat$orig.ident == "Day 30 Organoid Co-Culture")], t(nmf_d30cocul$H))
hist(nmf_d30cocul$W[,1])
hist(corr_d30cocul_nmf1)

plot(nmf_d30cocul$W[,1], seurat$nFeature_RNA[which(seurat$orig.ident == "Day 30 Organoid Co-Culture")])
plot(nmf_d30cocul$W[,1], seurat$nCount_RNA[which(seurat$orig.ident == "Day 30 Organoid Co-Culture")])
plotFeature(Embeddings(seurat,"umap_css")[which(seurat$orig.ident == "Day 30 Organoid Co-Culture"),], nmf_d30cocul$W[,1], colorPal=beach_colscheme, cex=0.5, random_order=T, sort_by_value=T)
plotFeature(Embeddings(seurat,"umap_css")[which(seurat$orig.ident == "Day 30 Organoid Co-Culture"),], corr_d30cocul_nmf1, colorPal=beach_colscheme, cex=0.5, random_order=T, sort_by_value=T)

## NMF on normalized no-scaling data (single component)
mat <- seurat@assays$RNA@data[VariableFeatures(seurat),which(seurat$orig.ident == "Day 30 Organoid Co-Culture")]
df_mat <- summary(mat)
mat <- sparseMatrix(i = df_mat$i, j = df_mat$j, x = exp(df_mat$x)-1, dims = c(nrow(mat),ncol(mat)), dimnames = dimnames(mat))
nmf_d30cocul_exp <- run_NMF(t(mat), n_components=1, max_iter = 5000)
colnames(nmf_d30cocul_exp$H) <- VariableFeatures(seurat)
rownames(nmf_d30cocul_exp$W) <- colnames(seurat)[which(seurat$orig.ident == "Day 30 Organoid Co-Culture")]
hist(nmf_d30cocul_exp$W[,1])

plot(nmf_d30cocul_exp$W[,1], seurat$nFeature_RNA[which(seurat$orig.ident == "Day 30 Organoid Co-Culture")])
plot(nmf_d30cocul_exp$W[,1], seurat$nCount_RNA[which(seurat$orig.ident == "Day 30 Organoid Co-Culture")])
plotFeature(Embeddings(seurat,"umap_css")[which(seurat$orig.ident == "Day 30 Organoid Co-Culture"),], nmf_d30cocul_exp$W[,1], colorPal=beach_colscheme, cex=0.5, random_order=T, sort_by_value=T)

## NMF on normalized log-scaled data (multiple component)
nmf_d30cocul_10comp <- run_NMF(t(seurat@assays$RNA@data[VariableFeatures(seurat),which(seurat$orig.ident == "Day 30 Organoid Co-Culture")]), n_components=10, max_iter = 5000)
colnames(nmf_d30cocul_10comp$H) <- VariableFeatures(seurat)
rownames(nmf_d30cocul_10comp$W) <- colnames(seurat)[which(seurat$orig.ident == "Day 30 Organoid Co-Culture")]
saveRDS(nmf_d30cocul_10comp, file="res.nmf_d30_coculture_organoid_10comp.rds")

nmf_d30ctrl_10comp <- run_NMF(t(seurat@assays$RNA@data[VariableFeatures(seurat),which(seurat$orig.ident == "Day 30 Organoid Control")]), n_components=10, max_iter = 5000)
colnames(nmf_d30ctrl_10comp$H) <- VariableFeatures(seurat)
rownames(nmf_d30ctrl_10comp$W) <- colnames(seurat)[which(seurat$orig.ident == "Day 30 Organoid Control")]

logcv_nmf_d30cocul <- log(apply(nmf_d30cocul_10comp$W,1,sd)/apply(nmf_d30cocul_10comp$W,1,mean))
hist(logcv_nmf_d30cocul, breaks=50, xlab="Log(CV) 10-NMF", cex.lab=1.5, main=NA)
hist(log(apply(nmf_d30ctrl_10comp$W,1,sd)/apply(nmf_d30ctrl_10comp$W,1,mean)), breaks=50, xlab="Log(CV) 10-NMF", cex.lab=1.5, main=NA)

## nFeature / nCount: should be higher in those with ambient RNAs
plotFeature(Embeddings(seurat,"umap_css"), seurat$nFeature_RNA/seurat$nCount_RNA, colorPal=bluered_colscheme, cex=0.5)
umi_ratio_d30cocul <- (seurat$nFeature_RNA/seurat$nCount_RNA)[which(seurat$orig.ident == "Day 30 Organoid Co-Culture")]
layout(matrix(1:2,nrow=1))
hist(umi_ratio_d30cocul, breaks=50, main=NA, xlab="# detected gene per UMI", cex.lab=1.5)
hist((seurat$nFeature_RNA/seurat$nCount_RNA)[which(seurat$orig.ident == "Day 30 Organoid Control")], breaks=50, main=NA, xlab="# detected gene per UMI", cex.lab=1.5)

## use log-scaled NMF (multiple component) and UMI ratio to do filtering
layout(matrix(1)); par(mar=c(5,5,1,1))
smoothScatter(umi_ratio_d30cocul, logcv_nmf_d30cocul)

mat <- data.frame(umi_ratio=umi_ratio_d30cocul, logcv=logcv_nmf_d30cocul)
cl_kmeans <- kmeans(mat,2)
centers <- sapply(1:2, function(i) colMeans(mat[cl_kmeans$cluster==i,]))

bandwidth <- diff(apply(cbind(umi_ratio_d30cocul, logcv_nmf_d30cocul), 2, quantile, probs = c(0.05,0.95), na.rm = TRUE, names = FALSE))/25
bandwidth[bandwidth == 0] <- 1
den <- KernSmooth::bkde2D(cbind(umi_ratio_d30cocul,logcv_nmf_d30cocul), bandwidth=bandwidth)
centers_idx <- data.frame(x1 = sapply(cl_kmeans$centers[,1], function(x) max(which(den$x1<x))),
                          x2 = sapply(cl_kmeans$centers[,2], function(x) max(which(den$x2<x))))
mat_borders <- t(sapply(1:length(den$x1), function(i) sapply(1:length(den$x2), function(j){
  d1 <- sqrt(sum(c(i-centers_idx[1,1], j-centers_idx[1,2])^2))
  d2 <- sqrt(sum(c(i-centers_idx[2,1], j-centers_idx[2,2])^2))
  d12 <- as.numeric(dist(centers_idx))
  return(d12 > d1 & d12 > d2)
})))
mat_borders <- Matrix(mat_borders & (den$fhat < 1.5 & den$fhat > 0.3)+0, sparse=T)
coord_borders <- summary(mat_borders)[,1:2]
coord_borders[,1] <- den$x1[coord_borders[,1]]
coord_borders[,2] <- den$x2[coord_borders[,2]]
m <- lm(j ~ i, data=coord_borders)
cl <- setNames(sign(predict(m, newdata = data.frame(i = mat[,1])) - mat[,2]), rownames(mat))

layout(matrix(1)); par(mar=c(5,5,1,1))
smoothScatter(umi_ratio_d30cocul, logcv_nmf_d30cocul, xlab="# detected genes per UMI", ylab="Log(CV-NMF)", cex.lab=1.5)
lines(seq(0,1,length.out=50), predict(m, newdata=data.frame(i=seq(0,1,length.out=50))), lwd=2, lty=2, col="#303030")

## try soupX and decontX
### soupX (the automatic way doesn't work)
seurat_d30_cocul <- NormalizeData(seurat_d30_cocul) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 20, verbose=F) %>%
  RunUMAP(dims=1:20)
seurat_d30_cocul <- FindNeighbors(seurat_d30_cocul, dims = 1:20) %>%
  FindClusters(resolution = 0.5)
sc <- load10X("rawdata/day30_org_cocult/outs/")
sc <- setClusters(sc, seurat_d30_cocul@active.ident)
sc <- setDR(sc, Embeddings(seurat_d30_cocul, "umap"))
sc <- autoEstCont(sc)
out <- adjustCounts(sc)

### decontX
sce_d30_cocul <- SingleCellExperiment(list(counts = seurat@assays$RNA@counts[,which(seurat$orig.ident=="Day 30 Organoid Co-Culture")]))
sce_d30_cocul <- decontX(sce_d30_cocul)
saveRDS(sce_d30_cocul, file="data.organoid_d30_coculture_decontx.sce.rds")

layout(matrix(1)); par(mar=c(5,5,1,1))
smoothScatter(umi_ratio_d30cocul, colData(sce_d30_cocul)[names(umi_ratio_d30cocul),"decontX_contamination"])
smoothScatter(logcv_nmf_d30cocul, colData(sce_d30_cocul)[names(umi_ratio_d30cocul),"decontX_contamination"])

## cells filtering
cells_d30_cocul_remained <- names(cl)[cl == -1]
cells_d30_cocul_excl <- setdiff(colnames(seurat)[which(seurat$orig.ident=="Day 30 Organoid Co-Culture")], cells_d30_cocul_remained)

par(mar=c(5,5,1,1))
hist(colData(sce_d30_cocul)$decontX_contamination, xlab="DecontX estimated ambient RNA prop", cex.lab=1.1, main=NA)
hist(colData(sce_d30_cocul)[cells_d30_cocul_remained,"decontX_contamination"], xlab=NA, main=NA, col="#303030", add=T)

seurat_filtered <- subset(seurat, cells = setdiff(colnames(seurat), cells_d30_cocul_excl))
seurat_filtered <- FindVariableFeatures(seurat_filtered, nfeatures = 3000)
blacklist <- c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat_filtered), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt", stringsAsFactors=F)[,1])
VariableFeatures(seurat_filtered) <- setdiff(VariableFeatures(seurat_filtered), blacklist)
seurat_filtered <- ScaleData(seurat_filtered, vars.to.regress=c("G2M.Score","S.Score")) %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20)

seurat_filtered <- cluster_sim_spectrum(seurat_filtered, label_tag="orig.ident", cluster_resolution = 1)
seurat_filtered <- run_PCA(seurat_filtered, reduction="css", npcs=20, reduction.name="css_pca", reduction.key="CSSPCA_")
seurat_filtered <- RunUMAP(seurat_filtered, reduction="css_pca", dims=1:20, reduction.name="umap_css", reduction.key="UMAPCSS_")
saveRDS(seurat_filtered, file="data.organoid_d3_d15_d30_cleanedup_integrated.seurat.rds")

VlnPlot(seurat_filtered, "nFeature_RNA", pt.size=0, cols = setNames(prettyrainbow_colscheme(6), sort(unique(seurat_filtered$orig.ident))))

(FeaturePlot(seurat_filtered, reduction="umap_css", features=c("DCN","SOX2","DCX","FOXG1","OTX2","HOXB2","LHX9","AIF1","MKI67"), cols=beach_colscheme(30), ncol=3, order=T, pt.size=0.01) & NoAxes() & NoLegend()) | (DimPlot(seurat_filtered, reduction="umap_css", group.by="orig.ident") & NoAxes())

par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat_filtered,"umap_css"), seurat_filtered$orig.ident, colorPal=prettyrainbow_colscheme, cex=0.5, pt_border=T)
legend("bottomright", bty="n", pch=21, pt.cex=1.5, col="#303030", pt.bg=prettyrainbow_colscheme(length(unique(seurat_filtered$orig.ident))), legend=sort(unique(seurat_filtered$orig.ident)))
plotMultiFeatures(Embeddings(seurat_filtered,"umap_css"), seurat_filtered@assays$RNA@data[c("DCN","SOX2","DCX","FOXG1","OTX2","HOXB2","LHX9","AIF1","MKI67"),], cex=0.4, colorPal=beach_colscheme, random_order=F, sort_by_value=T, cex.main=1.5)

seurat_filtered <- FindNeighbors(seurat_filtered, reduction="css", dims = 1:ncol(Embeddings(seurat_filtered,"css"))) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 0.8) %>%
  FindClusters(resolution = 1) %>%
  FindClusters(resolution = 2) %>%
  FindClusters(resolution = 3) %>%
  FindClusters(resolution = 20)
DimPlot(seurat_filtered, reduction="umap_css", group.by="RNA_snn_res.0.1", label=T) & NoLegend() & NoAxes()
DimPlot(seurat_filtered, reduction="umap_css", group.by="RNA_snn_res.3", label=T) & NoLegend() & NoAxes()

seurat_filtered$annot_celltype <- setNames(c("Proliferative hindbrain NPC","Hindbrain NPC","Hindbrain NPC/neuron","Hindbrain neuron","Hindbrain neuron","Proliferative diencephalic NPC","Diencephalic NPC","AIF1+ NPC","AIF1+ NPC/neuron","Diencephalic neuron","Dorsal telencephalic NPC","Mesenchyme-1","Mesenchyme-2","Mesenchyme-3","Unknown-1","Unknown-2","Unknown-3"),c(1,0,10,5,4,8,2,3,6,7,12,9,15,11,14,16,13))[as.character(seurat_filtered$RNA_snn_res.0.5)]
seurat_filtered$annot_celltype[seurat_filtered$RNA_snn_res.3 == 49] <- "Epithelium"
seurat_filtered$annot_celltype[seurat_filtered$RNA_snn_res.3 == 47] <- "Dorsal telencephalic neuron"
ct_ordered <- sort(unique(seurat_filtered$annot_celltype))[c(15,9,10,8,14,4,3,1,2,6,5,11,12,13,7,16,17,18)]
seurat_filtered$annot_celltype <- factor(seurat_filtered$annot_celltype, levels=ct_ordered)

cols_ct <- setNames(c(prettyrainbow_colscheme(length(unique(grep("Unknown",seurat_filtered$annot_celltype,invert=T, value=T)))),"#dedede","#cdcdcd","#bdbdbd"),ct_ordered)
DimPlot(seurat_filtered, reduction="umap_css", group.by="annot_celltype", cols=cols_ct, pt.size=0.4) & NoAxes()

par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat_filtered,"umap_css"), seurat_filtered$annot_celltype, colorPal=cols_ct, cex=0.6, pt_border=T, lwd_border=0.2)
legend("bottomright", pch=21, pt.cex=1.2, col="#303030", pt.bg=cols_ct, legend=names(cols_ct), bty="n", cex=0.7)

saveRDS(seurat_filtered, file="data.organoid_d3_d15_d30_cleanedup_integrated.seurat.rds")






# secondary ambient RNA cleanup: decontaminate the remaining d30-co-culture organoid cells (naive decontX)
sce_d30_cocul_filtered <- SingleCellExperiment(list(counts = seurat_filtered@assays$RNA@counts[,which(seurat_filtered$orig.ident=="Day 30 Organoid Co-Culture")]))
sce_d30_cocul_filtered <- decontX(sce_d30_cocul_filtered)
saveRDS(sce_d30_cocul_filtered, file="data.organoid_d30_coculture_filtered_decontx.sce.rds")

d30_decont_counts <- decontXcounts(sce_d30_cocul_filtered)

## check DE changes caused by the naive decontamination
seurat_d30cocul_filtered <- subset(seurat_filtered, subset = orig.ident == "Day 30 Organoid Co-Culture") %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000)
VariableFeatures(seurat_d30cocul_filtered) <- setdiff(VariableFeatures(seurat_d30cocul_filtered), blacklist)
seurat_d30cocul_filtered <- ScaleData(seurat_d30cocul_filtered) %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5)
UMAPPlot(seurat_d30cocul_filtered, label=T) & NoLegend() & NoAxes()

seurat_d30cocul_adj <- CreateSeuratObject(counts = d30_decont_counts[,colnames(seurat_d30cocul_filtered)], meta.data = seurat_d30cocul_filtered@meta.data) %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000)
VariableFeatures(seurat_d30cocul_adj) <- setdiff(VariableFeatures(seurat_d30cocul_adj), blacklist)
seurat_d30cocul_adj <- ScaleData(seurat_d30cocul_adj) %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5)
UMAPPlot(seurat_d30cocul_adj, label=T) & NoLegend() & NoAxes()

seurat_d30cocul_adj$cl_orig <- seurat_d30cocul_filtered$RNA_snn_res.0.5
(UMAPPlot(seurat_d30cocul_filtered, label=T) | UMAPPlot(seurat_d30cocul_adj, group.by="cl_orig", label=T)) & NoLegend() & NoAxes()

DE_after_decontx <- setNames(lapply(levels(seurat_d30cocul_filtered@active.ident), function(cl){
  idx <- which(seurat_d30cocul_filtered@active.ident == cl)
  mat <- cbind(seurat_d30cocul_filtered@assays$RNA@data[,idx], seurat_d30cocul_adj@assays$RNA@data[,idx])
  lab <- rep(c("original","corrected"), each=length(idx))
  return(wilcoxauc(mat, lab))
}), levels(seurat_d30cocul_filtered@active.ident))
sapply(DE_after_decontx, function(x) sapply(c("original","corrected"), function(y) sum(x$group[which(x$auc>0.7 & x$logFC>log(1.2) & x$padj<0.01)]==y)))
lapply(DE_after_decontx, function(x) x$feature[which(x$auc>0.7 & x$logFC>log(1.2) & x$padj<0.01)])
head(sort(table(unlist(lapply(DE_after_decontx, function(x) x$feature[which(x$auc>0.7 & x$logFC>log(1.2) & x$padj<0.01)]))), decreasing=T)) 

(FeaturePlot(seurat_d30cocul_filtered,"TIMP1",reduction="umap", cols=beach_colscheme(30), order=T, pt.size=0.6) | FeaturePlot(seurat_d30cocul_adj,"TIMP1",reduction="umap",cols=beach_colscheme(30), order=T, pt.size=0.6)) & NoLegend() & NoAxes()
(FeaturePlot(seurat_d30cocul_filtered,"TTR",reduction="umap", cols=beach_colscheme(30), order=T, pt.size=0.6) | FeaturePlot(seurat_d30cocul_adj,"TTR",reduction="umap",cols=beach_colscheme(30), order=T, pt.size=0.6)) & NoLegend() & NoAxes()
plotMultiFeatures(Embeddings(seurat_filtered,"umap_css"), seurat_filtered@assays$RNA@data[c(c("TTR","TIMP1")),], colorPal=beach_colscheme, cex=0.5, cex.main=1.5, emphasize=which(seurat_filtered$orig.ident!="Day 30 Organoid Co-Culture"), ncol=2)

## integrate with other organoid data
seurat_filtered$sample <- seurat_filtered$orig.ident
counts_adj <- seurat_filtered@assays$RNA@counts
counts_adj[rownames(d30_decont_counts),colnames(d30_decont_counts)] <- d30_decont_counts
seurat_filtered_adj <- CreateSeuratObject(counts_adj, meta.data = seurat_filtered@meta.data)

seurat_filtered_adj <- NormalizeData(seurat_filtered_adj) %>%
  FindVariableFeatures(nfeatures = 3000)
blacklist <- c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat_filtered), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt", stringsAsFactors=F)[,1])
VariableFeatures(seurat_filtered_adj) <- setdiff(VariableFeatures(seurat_filtered_adj), blacklist)
seurat_filtered_adj <- CellCycleScoring(seurat_filtered_adj, g2m.features=cc.genes.updated.2019$g2m.genes, s.features=cc.genes.updated.2019$s.genes)
seurat_filtered_adj <- ScaleData(seurat_filtered_adj, vars.to.regress=c("G2M.Score","S.Score")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20) %>%
  cluster_sim_spectrum(label_tag="orig.ident", cluster_resolution = 1) %>%
  run_PCA(reduction="css", npcs=20, reduction.name="css_pca", reduction.key="CSSPCA_") %>%
  RunUMAP(reduction="css_pca", dims=1:20, reduction.name="umap_css", reduction.key="UMAPCSS_")

(FeaturePlot(seurat_filtered_adj, reduction="umap_css", features=c("DCN","SOX2","DCX","FOXG1","OTX2","HOXB2","LHX9","AIF1","MKI67"), cols=beach_colscheme(30), ncol=3, order=T, pt.size=0.01) & NoAxes() & NoLegend()) | (DimPlot(seurat_filtered_adj, reduction="umap_css", group.by="orig.ident") & NoAxes())

seurat_filtered_adj <- FindNeighbors(seurat_filtered_adj, reduction="css", dims = 1:ncol(Embeddings(seurat_filtered_adj,"css"))) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 0.8) %>%
  FindClusters(resolution = 1) %>%
  FindClusters(resolution = 2) %>%
  FindClusters(resolution = 3) %>%
  FindClusters(resolution = 20)
saveRDS(seurat_filtered_adj, file="data.organoid_d3_d15_d30_cleanedup_corrected_integrated.seurat.rds")

DimPlot(seurat_filtered_adj, reduction="umap_css", group.by="RNA_snn_res.0.1", label=T) & NoLegend() & NoAxes()
DimPlot(seurat_filtered_adj, reduction="umap_css", group.by="RNA_snn_res.0.5", label=T) & NoLegend() & NoAxes()
DimPlot(seurat_filtered_adj, reduction="umap_css", group.by="RNA_snn_res.1", label=T) & NoLegend() & NoAxes()

