library(Matrix)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(simspec)
library(presto)
source("~/Tools/scripts/feature_plots.r")
source("~/Tools/scripts/r_util.r")

setwd("GSE141862_scRNAseq")
files <- grep(".txt.gz$", list.files(), value=T)
counts_samples <- lapply(files, function(file) read.table(file, header=T, row.names=1))
names(counts_samples) <- sapply(strsplit(files, "_"),"[", 1)
genes <- sort(unique(unlist(lapply(counts_samples, rownames))))
counts_samples <- lapply(counts_samples, function(x){
  mat <- Matrix(0, nrow=length(genes), ncol=ncol(x), dimnames = list(genes,colnames(x)))
  mat[rownames(x),] <- Matrix(as.matrix(x), sparse=T)
  return(mat)
})
counts <- do.call(cbind, counts_samples)
colnames(counts) <- paste0(rep(names(counts_samples),sapply(counts_samples,ncol)),"_", unlist(sapply(counts_samples,colnames)))

meta_samples <- as.data.frame(t(read.csv("../GSE141862_family.csv", header=T, row.names=1)))
meta <- data.frame(meta_samples[rep(names(counts_samples), sapply(counts_samples, ncol)),c("geo_accession","source_name_ch1","characteristics_ch1","characteristics_ch1-1")], barcode = unlist(lapply(counts_samples, colnames)), row.names = colnames(counts))
colnames(meta)[2:4] <- c("source","fetus","timepoint")
meta$timepoint_days <- apply(sapply(strsplit(sapply(strsplit(meta$timepoint, ": "), "[", 2),"\\+"), as.numeric), 2, function(x) 7*x[1]+x[2])

save(counts, meta, file="counts_meta.RData")

genes <- read.table("/links/groups/treutlein/USERS/zhisong_he/Work/databases/Ensembl/v96/human_gencode30/human_ens_symbol_type.txt", header=T,sep="\t")
rownames(genes) <- make.unique(genes[,1])
rownames(counts) <- make.unique(genes[rownames(counts),2])
save(counts, meta, file="counts_meta_symbols.RData")

setwd("..")
seurat <- CreateSeuratObject(counts, meta.data=meta, project="human_dev_microglia")
seurat <- PercentageFeatureSet(seurat, pattern="^MT-", col.name="percent.mito")
seurat <- PercentageFeatureSet(seurat, features=read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1], col.name="percent.ribo")
seurat$log10NCount_RNA <- log10(seurat$nCount_RNA)
VlnPlot(seurat, c("nFeature_RNA","percent.mito"), group.by="orig.ident", pt.size=0)
mean_log10_count <- mean(seurat$log10NCount_RNA)
mad_log10_count <- mean(abs(seurat$log10NCount_RNA-mean(seurat$log10NCount_RNA)))
seurat <- subset(seurat, subset = log10NCount_RNA < mean_log10_count+3*mad_log10_count & log10NCount_RNA > mean_log10_count-3*mad_log10_count & percent.mito < 10) %>%
  NormalizeData() %>%
  CellCycleScoring(g2m.features=cc.genes.updated.2019$g2m.genes, s.features=cc.genes.updated.2019$s.genes) %>%
  FindVariableFeatures(nfeatures = 5000)
blacklist <- c(grep("^MT-", rownames(seurat), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1])
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
VariableFeatures(seurat) <- intersect(VariableFeatures(seurat), genes$Gene.name[genes$Gene.type=="protein_coding"])
seurat <- ScaleData(seurat, vars.to.regress = c("nCount_RNA","percent.mito","percent.ribo")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
set.seed(42)

umap_model <- uwot::umap(Embeddings(seurat,"pca")[,1:20], metric="cosine", ret_model=T, n_neighbors=30L, n_components=2L, n_epochs=NULL, learning_rate=1.0, min_dist=0.3, spread=1.0, set_op_mix_ratio=1.0, local_connectivity=1L, repulsion_strength=1, negative_sample_rate=5, a=NULL, b=NULL, fast_sgd=F)
umap_model$embedding <- Embeddings(seurat,"umap")
save_uwot(umap_model, file="uwot_model_pca.uwot")

seurat$timepoint_wks <- round(seurat$timepoint_days/7)
seurat$sample_group <- paste0(seurat$timepoint_wks, "_", seurat$source)

seurat <- FindNeighbors(seurat, dims = 1:20) %>%
  FindClusters(resolution = 0.3) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1) %>%
  FindClusters(resolution = 10)
UMAPPlot(seurat, group.by="timepoint_wks") & NoLegend() & NoAxes()
UMAPPlot(seurat, group.by="RNA_snn_res.0.5") & NoLegend() & NoAxes()
FeaturePlot(seurat, c("MAP1B","AXL","DUSP1","MKI67"), reduction="umap", cols=beach_colscheme(30), order=T) & NoLegend() & NoAxes()


# test MNN integration
seurat_mnn <- RunFastMNN(SplitObject(seurat,"timepoint_wks"))
seurat[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn,"mnn")[colnames(seurat),], key="MNN_")
seurat <- RunUMAP(seurat, reduction="mnn", dims=1:20, reduction.name="umap_mnn", reduction.key="UMAPMNN_")
DimPlot(seurat, group.by="timepoint_wks", reduction="umap_mnn") & NoLegend() & NoAxes()

# test CSS integration
seurat <- cluster_sim_spectrum(seurat, label_tag="timepoint_wks", cluster_resolution = 1)
seurat@reductions$css@cell.embeddings[is.na(seurat@reductions$css@cell.embeddings)] <- rnorm(sum(is.na(seurat@reductions$css@cell.embeddings)))
seurat <- run_PCA(seurat, reduction="css", npcs = 20, reduction.name = "css_pca", reduction.key="CSSPCA_")
seurat <- RunUMAP(seurat, reduction="css_pca", dims=1:10, reduction.name="umap_css", reduction.key="UMAPCSS_")

seurat <- FindNeighbors(seurat, reduction="css_pca", dims = 1:10)
seurat[['RNA_css_nn']] <- seurat[['RNA_nn']]
seurat[['RNA_css_snn']] <- seurat[['RNA_snn']]
seurat <- FindClusters(seurat, graph.name="RNA_css_snn", resolution = 0.3) %>%
  FindClusters(graph.name="RNA_css_snn", resolution = 0.5)

DimPlot(seurat, group.by="timepoint_wks", reduction="umap_css") & NoLegend() & NoAxes()
FeaturePlot(seurat, c("MAP1B","AXL","CX3CR1","MKI67"), reduction="umap_css", cols=beach_colscheme(30), order=T) & NoLegend() & NoAxes()
par(mar=c(1,1,1,1)); plotFeature(Embeddings(seurat,"umap_css"), factor(seurat$timepoint_wks), colorPal=bluered_colscheme, cex=0.6, pt_border=T, lwd_border=0.1, do_label=T, cex.label=1, label_round=T)
DimPlot(seurat, reduction="umap_css") & NoLegend() & NoAxes()

saveRDS(seurat, file="seurat.rds")


# use no-integration results for further annotation
seurat@active.ident <- setNames(seurat$RNA_snn_res.0.5, colnames(seurat))
DE_cl <- wilcoxauc(seurat)
markers_cl <- lapply(split(DE_cl, DE_cl$group), function(x){ res <- x[which(x$logFC > log(1.2) & x$auc > 0.7 & x$padj < 0.01 & x$pct_in > 20 & x$pct_in - x$pct_out > 20),]; res <- res[order(res$auc, decreasing=T),]; return(res) })

FeaturePlot(seurat, c("GYPA","SOX4","ITGAM","FCN1","TMEM119","PTPRC","AIF1","DUSP1","MKI67"), reduction="umap", cols=beach_colscheme(30), order=T) & NoLegend() & NoAxes()

annot <- c("9" = "Red blood cell",
           "5" = "Monocyte",
           "6" = "Perivascular macrophage",
           "7" = "SOX-hi microglia",
           "2" = "SOX-medium microglia",
           "4" = "Prolifer. microglia",
           "1" = "Homeostatic microglia (early)",
           "0" = "Homeostatic microglia (late)",
           "8" = "Activated microglia",
           "3" = "CCL-high microglia")
seurat$celltype <- factor(annot[as.character(seurat@active.ident)], levels=annot)
cols_cl <- setNames(c("#A569BD","#E59866","#E74C3C","#F4D03F","#82E0AA","#5DADE2","#52BE80","#27AE60","#48C9B0","#17A589"), annot)
par(mar=c(1,1,1,1)); plotFeature(Embeddings(seurat, "umap"), seurat$celltype, colorPal=cols_cl, cex=0.6, pt_border=T, lwd_border=0.3, do_legend = T, legend_pos = "bottomright", legend_cex = 0.7)
par(mar=c(1,1,1,1)); plotFeature(Embeddings(seurat, "umap"), factor(seurat$timepoint_wks), colorPal=prettyrainbow_colscheme, cex=0.6, pt_border=T, lwd_border=0.3, do_legend = T, legend_pos = "bottomright", legend_cex = 0.7)
FeaturePlot(seurat, c("GYPA","FCN1","CD163","SOX4","LGALS3","CCL4","MKI67","PTPRC","AIF1"), reduction="umap", cols=beach_colscheme(30), order=T) & NoLegend() & NoAxes()

saveRDS(seurat, file="seurat.rds")


# high-res clusters for projection
avg_pca_niches <- t(sapply(levels(seurat$RNA_snn_res.10), function(cl) colMeans(Embeddings(seurat,"pca")[which(seurat$RNA_snn_res.10==cl),1:20])))
avg_expr_niches <- sapply(levels(seurat$RNA_snn_res.10), function(cl) rowMeans(seurat@assays$RNA@data[,which(seurat$RNA_snn_res.10==cl),]))
umap_niche_projected <- umap_transform(X=avg_pca_niches, model=umap_model)
meta_niches <- setNames(data.frame(sapply(c("source","timepoint_days","timepoint_wks","celltype"), function(col) sapply(levels(seurat$RNA_snn_res.10), function(cl) names(which.max(table(seurat@meta.data[which(seurat$RNA_snn_res.10 == cl),col])))))), c("source","timepoint_days","timepoint_wks","celltype"))
meta_niches$timepoint_days <- as.numeric(meta_niches$timepoint_days)
meta_niches$timepoint_wks <- as.numeric(meta_niches$timepoint_days)
meta_niches$celltype <- factor(meta_niches$celltype, levels=annot)
vargenes_niches <- VariableFeatures(seurat)

## binarized normalized data
bin_counts <- seurat@assays$RNA@counts > 0
bin_counts <- t(rowNorm(t(bin_counts)))
avg_binexpr_niches <- sapply(levels(seurat$RNA_snn_res.10), function(cl) rowMeans(bin_counts[,which(seurat$RNA_snn_res.10==cl),]))

save(avg_expr_niches, avg_binexpr_niches, vargenes_niches, umap_niche_projected, meta_niches, file="niches_avgExpr_meta.RData")
par(mar=c(1,1,1,1)); plotFeature(Embeddings(seurat, "umap"), seurat$celltype, colorPal=cols_cl, cex=0.5, pt_border=F, lwd_border=0.3, do_legend = T, legend_pos = "bottomright", legend_cex = 0.7)
points(umap_niche_projected, pch=21, col="#303030", lwd=1.5, bg=cols_cl[meta_niches$celltype], cex = 2.5)

# check relationship between detection rates and gene lengths
detection_rate <- rowMeans(seurat@assays$RNA@data > 0)

gene_symbols <- read.table("~/Work/databases/Ensembl/v96/human_gencode30/human_ens_symbol_type.txt", sep="\t", header=T)
rownames(gene_symbols) <- gene_symbols[,1]
gene_len <- read.table("~/Work/databases/Ensembl/v96/human_gencode30/Homo_sapiens.GRCh38.96.chr.struc")
gene_len <- setNames(gene_len[,4], gene_len[,1])
names(gene_len) <- make.unique(gene_symbols[names(gene_len),2])
gene_len <- gene_len[rownames(seurat)]

layout(matrix(1:2,nrow=1)); par(mar=c(5,5,1,1))
boxplot(tapply(detection_rate, ceiling(rank(gene_len)/nrow(seurat)*10), list), outline=F, ylab = "Detection rate", xlab="Gene length ranks", cex.lab=1.3)
smoothScatter(log10(gene_len), log10(detection_rate), ylab = "Log10(Detection rate)", xlab="Log10(Gene length)", cex.lab=1.3)
cor(detection_rate, gene_len, use="complete.obs", method="spearman")
cor(detection_rate[gene_len>1000], gene_len[gene_len>1000], use="complete.obs", method="spearman")
cor(log(detection_rate), log(gene_len), use="complete.obs", method="pearson")