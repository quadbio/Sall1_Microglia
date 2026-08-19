library(Matrix)
library(Seurat)
library(SeuratWrappers)
library(simspec)
library(dplyr)
library(presto)
source("~/Tools/scripts/feature_plots.r")
source("~/Tools/scripts/differential_expression.r")
source("~/Tools/scripts/r_util.r")

# load d3 and d15 microglia data
seurat_d3_ctrl <- readRDS("/links/groups/treutlein/USERS/marinani/shared_data/microglia_organoid_coculture/day3/d3_mg_ctrl.rds")
seurat_d3_ctrl$sample <- "MG-d3-ctrl"
seurat_d3_ctrl$condition <- "CTRL"
seurat_d3_ctrl$timepoint <- "D3"

seurat_d3_cocul <- readRDS("/links/groups/treutlein/USERS/marinani/shared_data/microglia_organoid_coculture/day3/d3_mg_cocult.rds")
seurat_d3_cocul$sample <- "MG-d3-cocul"
seurat_d3_cocul$condition <- "COCUL"
seurat_d3_cocul$timepoint <- "D3"

seurat_d15_ctrl <- readRDS("/links/groups/treutlein/USERS/marinani/shared_data/microglia_organoid_coculture/day15/d15_mg_ctrl.rds")
seurat_d15_ctrl$sample <- "MG-d15-ctrl"
seurat_d15_ctrl$condition <- "CTRL"
seurat_d15_ctrl$timepoint <- "D15"

seurat_d15_cocul <- readRDS("/links/groups/treutlein/USERS/marinani/shared_data/microglia_organoid_coculture/day15/d15_mg_cocult.rds")
seurat_d15_cocul$sample <- "MG-d15-cocul"
seurat_d15_cocul$condition <- "COCUL"
seurat_d15_cocul$timepoint <- "D15"

# load d30 microglia data
counts_d30_ctrl <- Read10X_h5("/links/groups/treutlein/DATA/sequencing/20201207_P1380_ANDREA_10X_coculture_day_30/processed/day30_mg_ctrl/outs/filtered_feature_bc_matrix.h5")
counts_d30_cocul <- Read10X_h5("/links/groups/treutlein/DATA/sequencing/20201207_P1380_ANDREA_10X_coculture_day_30/processed/day30_mg_cocult/outs/filtered_feature_bc_matrix.h5")

seurat_d30_ctrl <- CreateSeuratObject(counts_d30_ctrl, project="coculture", meta.data = data.frame(sample = rep("MG-d30-ctrl",ncol(counts_d30_ctrl)), condition="CTRL", timepoint="D30", row.names = colnames(counts_d30_ctrl)))
seurat_d30_cocul <- CreateSeuratObject(counts_d30_cocul, project="coculture", meta.data = data.frame(sample = rep("MG-d30-cocul",ncol(counts_d30_cocul)), condition="COCUL", timepoint="D30", row.names = colnames(counts_d30_cocul)))

# merge all data
seurat <- merge(seurat_d30_ctrl, list(seurat_d30_cocul, seurat_d15_ctrl, seurat_d15_cocul, seurat_d3_ctrl, seurat_d3_cocul))
seurat <- NormalizeData(seurat) %>%
  PercentageFeatureSet("^MT-", col.name = "percent.mt")

saveRDS(seurat, file="data/data.microglia_d3_d15_d30_merged.rds")

# QC and integration
seurat <- subset(seurat, subset = nFeature_RNA > 2000 & percent.mt < 15)
seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  CellCycleScoring(g2m.features = cc.genes.updated.2019$g2m.genes, s.features = cc.genes.updated.2019$s.genes)
blacklist <- c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1])
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
seurat <- ScaleData(seurat, vars.to.regress = c("G2M.Score","S.Score")) %>%
  RunPCA(npcs = 30, verbose=F) %>%
  RunUMAP(dims = 1:20)

seurat <- cluster_sim_spectrum(seurat, label_tag="sample", cluster_resolution = 1)
seurat <- RunUMAP(seurat, reduction="css", dims=1:ncol(Embeddings(seurat,"css")), reduction.name="umap_css", reduction.key="UMAPCSS_")

seurat_mnn <- RunFastMNN(SplitObject(seurat,"sample"))
seurat[['mnn']] <- CreateDimReducObject(Embeddings(seurat_mnn,"mnn")[colnames(seurat),], key="MNN_")
seurat <- RunUMAP(seurat, reduction="mnn", dims=1:20, reduction.name="umap_mnn", reduction.key="UMAPMNN_")

(UMAPPlot(seurat, group.by="sample") & NoLegend() & NoAxes()) | (DimPlot(seurat, group.by="sample", reduction="umap_css") & NoLegend() & NoAxes()) | (DimPlot(seurat, group.by="sample", reduction="umap_mnn") & NoAxes())

seurat <- FindNeighbors(seurat, reduction="css", dims=1:ncol(Embeddings(seurat,"css"))) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 0.8) %>%
  FindClusters(resolution = 10)
seurat@active.ident <- setNames(seurat$RNA_snn_res.0.5, colnames(seurat))

saveRDS(seurat, file="data/data.microglia_d3_d15_d30_processed.rds")

par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat,"umap_css"), seurat@active.ident, colorPal=prettyrainbow_colscheme, cex=0.6, pt_border=T, lwd_border=0.1, do_label=T, cex.main = 1.5, label_round=T)

# annotate clusters by markers
DE <- wilcoxauc(seurat)
DEG_cl <- lapply(split(DE,DE$group), function(x){
  res <- x[which(x$logFC > log(1.5) & x$auc > 0.7 & x$padj < 0.01 & x$pct_in > 30 & x$pct_in - x$pct_out > 30),]
  res <- res[order(res$logFC, decreasing=T),]
})
write.table(do.call(rbind, DEG_cl), file="tab.DEG_microglia_cl.tsv", quote=F, sep="\t", row.names=F)
plotMultiFeatures(Embeddings(seurat,"umap_css"), seurat@assays$RNA@data[c("AIF1","APOE","FOSB","CCL13","NME1","MKI67","ACTA2","DCX","SOX2"),], colorPal=beach_colscheme(30), random_order=F, sort_by_value=T, cex=0.5, cex.main=1.8, mar=c(1,1,3,3))

cl_annot <- setNames(c("MG (synaptic pruning)","MG (pro-inflammatory)","MG (homeostatic/neurotrophic)","MG (pro-/anti-inflammatory)","MG (cell cycle)","MG (cell cycle)","SMC-like","NPC-like","Neuron-like","Neuron-like (phagocytosed)"), levels(seurat@active.ident)) # annotation by Marina
seurat$annot_celltype <- factor(cl_annot[seurat@active.ident], levels=unique(cl_annot)[c(1,4,2,3,5,6,7,8,9)])
par(mar=c(1,1,1,1))
plotFeature(Embeddings(seurat,"umap_css"), seurat$annot_celltype, colorPal=prettyrainbow_colscheme, cex=0.6, pt_border=T, lwd_border=0.1, do_label=T, cex.main = 1.5, label_round=F)






# compare to the human primary microglia
## logCP10k-primary vs. logUMInorm-organoid
load("~/Work/public_datasets/Kracht_Science_2020_human_fetal_microglia_scRNAseq_ATACseq/niches_avgExpr_meta.RData")
sim2_MG_niches <- ref_sim_spectrum(seurat@assays$RNA@data, avg_expr_niches[vargenes_niches,], scale=F)
dimnames(sim2_MG_niches) <- list(colnames(seurat), colnames(avg_expr_niches))
cols_MG_ct <- setNames(c("#A569BD","#E59866","#E74C3C","#F4D03F","#82E0AA","#5DADE2","#52BE80","#27AE60","#48C9B0","#17A589"), levels(meta_niches$celltype))

par(mar=c(1,1,1,1)); plotFeature(Embeddings(seurat,"umap_css"), meta_niches$celltype[apply(sim2_MG_niches,1,which.max)], colorPal=cols_MG_ct, cex=0.6, pt_border=T, lwd_border=0.1, do_legend=T, legend_pos = "bottomleft", legend_cex = 1)
par(mar=c(1,1,1,1)); plotFeature(Embeddings(seurat,"umap_css"), apply(sim2_MG_niches,1,max), colorPal=bluered_colscheme, cex=0.6, pt_border=F)
layout(matrix(1:2,ncol=1))
par(mar=c(10,5,1,1)); boxplot(tapply(apply(sim2_MG_niches,1,max), seurat$annot_celltype, list), las=2, frame=F, pch=16, cex.axis=0.8, ylab="max(R)", cex.lab=1.3)
idx_MG <- which(!seurat@active.ident %in% as.character(6:9))
par(mar=c(10,5,1,1)); boxplot(tapply(apply(sim2_MG_niches,1,max)[idx_MG], seurat$sample[idx_MG], list), las=2, frame=F, pch=16, cex.axis=0.8, ylab="max(R)", cex.lab=1.3)

avg_expr_highres <- sapply(levels(seurat$RNA_snn_res.10), function(cl) rowMeans(seurat@assays$RNA@data[,which(seurat$RNA_snn_res.10==cl)]))
sim2_highres <- ref_sim_spectrum(seurat@assays$RNA@data, avg_expr_highres[intersect(vargenes_niches,rownames(avg_expr_highres)),], scale=F)
summary(apply(sim2_highres,1,max)[idx_MG])

## binarized norm vs. binarized norm
bincounts <- seurat@assays$RNA@data > 0
sim2_MG_niches_bin <- ref_sim_spectrum(bincounts, avg_binexpr_niches[vargenes_niches,], scale=F)
dimnames(sim2_MG_niches_bin) <- list(colnames(seurat), colnames(avg_binexpr_niches))

par(mar=c(1,1,1,1)); plotFeature(Embeddings(seurat,"umap_css"), meta_niches$celltype[apply(sim2_MG_niches_bin,1,which.max)], colorPal=cols_MG_ct, cex=0.6, pt_border=T, lwd_border=0.1, do_legend=T, legend_pos = "bottomleft", legend_cex = 1)
par(mar=c(1,1,1,1)); plotFeature(Embeddings(seurat,"umap_css"), apply(sim2_MG_niches_bin,1,max), colorPal=bluered_colscheme, cex=0.6, pt_border=F)
layout(matrix(1:2,ncol=1))
par(mar=c(10,5,1,1)); boxplot(tapply(apply(sim2_MG_niches_bin,1,max), seurat$annot_celltype, list), las=2, frame=F, pch=16, cex.axis=0.8, ylab="max(R)", cex.lab=1.3)
idx_MG <- which(!seurat@active.ident %in% as.character(6:9))
par(mar=c(10,5,1,1)); boxplot(tapply(apply(sim2_MG_niches_bin,1,max)[idx_MG], seurat$sample[idx_MG], list), las=2, frame=F, pch=16, cex.axis=0.8, ylab="max(R)", cex.lab=1.3)

# check relationship between gene length and expression values
mean_expr <- rowMeans_drop0(expCounts(seurat@assays$RNA@data))

gene_symbols <- read.table("/links/groups/treutlein/DATA/sequencing/20210118_P1412_MARINA_10X_RESEQ_coculture_day_30_org_reseq/processed/day30_org_cocult/outs/filtered_feature_bc_matrix/features.tsv.gz", sep="\t")
rownames(gene_symbols) <- gene_symbols[,1]
gene_len <- read.table("~/Work/databases/cellranger_v3.struc")
gene_len <- setNames(gene_len[,4], gene_len[,1])
names(gene_len) <- make.unique(gene_symbols[names(gene_len),2])
gene_len <- gene_len[rownames(seurat)]

layout(matrix(1:2,nrow=1)); par(mar=c(5,5,1,1))
boxplot(tapply(log(mean_expr+1), ceiling(rank(gene_len)/nrow(seurat)*10), list), outline=F, ylab = "Log_p1(Average expression (nonzero))", xlab="Gene length ranks", cex.lab=1.3)
smoothScatter(log10(gene_len), log(mean_expr), ylab = "Log(Average expression (nonzero))", xlab="Log10(Gene length)", cex.lab=1.3)
cor(log(mean_expr), gene_len, use="complete.obs", method="spearman")
cor(log(mean_expr[gene_len>1000]), gene_len[gene_len>1000], use="complete.obs", method="spearman")

# check relationship between gene length and detection rates
detection_rate <- rowMeans(seurat@assays$RNA@data > 0)

layout(matrix(1:2,nrow=1)); par(mar=c(5,5,1,1))
boxplot(tapply(detection_rate, ceiling(rank(gene_len)/nrow(seurat)*10), list), outline=F, ylab = "Detection rate", xlab="Gene length ranks", cex.lab=1.3)
smoothScatter(log10(gene_len), log10(detection_rate), ylab = "Log10(Detection rate)", xlab="Log10(Gene length)", cex.lab=1.3)
cor(detection_rate, gene_len, use="complete.obs", method="spearman")
cor(detection_rate[gene_len>1000], gene_len[gene_len>1000], use="complete.obs", method="spearman")
cor(log(detection_rate), log(gene_len), use="complete.obs", method="pearson")


## DE to primary fetal microglia
seurat_MG <- readRDS("~/Work/public_datasets/Kracht_Science_2020_human_fetal_microglia_scRNAseq_ATACseq/seurat.rds")

### all co-culture microglia
idx_invitro <- which(! seurat$annot_celltype %in% c("SMC-like","NPC-like","Neuron-like","Neuron-like (phagocytosed)") & seurat$condition == "COCUL")
detection_rates_invitro <- rowMeans(seurat@assays$RNA@data[,idx_invitro] > 0)
detection_rates_invivo <- rowMeans(seurat_MG@assays$RNA@data > 0)

shared_genes <- intersect(union(names(which(detection_rates_invitro > 0.01)), names(which(detection_rates_invivo > 0.01))), intersect(names(detection_rates_invitro),names(detection_rates_invivo)))
det_mat <- cbind(seurat@assays$RNA@data[shared_genes,idx_invitro] > 0, seurat_MG@assays$RNA@data[shared_genes,] > 0)
cond <- factor(rep(c("invitro","invivo"), c(length(idx_invitro),ncol(seurat_MG))))
coverage <- c(seurat$nFeature_RNA[idx_invitro], seurat_MG$nFeature_RNA)
ct <- factor(c(as.character(meta_niches$celltype[apply(sim2_MG_niches_bin,1,which.max)[idx_invitro]]),
               as.character(seurat_MG$celltype)),
             levels = levels(seurat_MG$celltype))

DE_invitro_invivo <- binom_group_test(det_mat,
                                      group = cond,
                                      coverage = coverage,
                                      covar = data.frame(ct = ct),
                                      num_threads = 50,
                                      return_coef_group = "invivo")
saveRDS(DE_invitro_invivo, file="res.DE-binom_cocultured_invivo_microglia.rds")
blacklist <- c(grep("^MT-", rownames(det_mat), value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1])
DEG_invitro_invivo <- setdiff(rownames(DE_invitro_invivo),blacklist)[which(p.adjust(DE_invitro_invivo$p_Resi[!rownames(DE_invitro_invivo) %in% blacklist]*2, method="BH")<0.05)]

layout(matrix(1)); par(mar=c(5,5,1,1))
plot(detection_rates_invitro[shared_genes],
     detection_rates_invivo[shared_genes],
     xlab = "Detection rate (in vitro, co-cultured)", ylab = "Detection rate (in vivo)", cex.lab = 1.5, bty = "n", col = "#69696960", pch=16, cex=1)
points(detection_rates_invitro[DEG_invitro_invivo], detection_rates_invivo[DEG_invitro_invivo], col="#303030", bg="red", pch=21, cex=1.5)

plot(log10(detection_rates_invitro[shared_genes_d30]),
     log10(detection_rates_invivo[shared_genes_d30]),
     xlab = "Log10-Detection rate (in vitro, co-cultured)", ylab = "Log10-Detection rate (in vivo)", cex.lab = 1.5, bty = "n", col = "#69696960", pch=16, cex=1)
points(log10(detection_rates_invitro[DEG_invitro_invivo]), log10(detection_rates_invivo[DEG_invitro_invivo]), col="#303030", bg="red", pch=21, cex=1.5)

human_TFs <- read.table("~/Work/databases/animalTFDB/animalTFDB_v3/Homo_sapiens_TF.tsv", sep="\t", header=T)
DEG_TF_invitro_invivo <- sort(intersect(DEG_invitro_invivo, human_TFs$Symbol))
sort(setNames(DE_invitro_invivo[DEG_TF_invitro_invivo,"coef"], DEG_TF_invitro_invivo))

### only d30 co-culture microglia
idx_invitro_d30 <- which(! seurat$annot_celltype %in% c("SMC-like","NPC-like","Neuron-like","Neuron-like (phagocytosed)") & seurat$condition == "COCUL" & seurat$timepoint == "D30")
detection_rates_invitro_d30 <- rowMeans(seurat@assays$RNA@data[,idx_invitro_d30] > 0)

shared_genes_d30 <- intersect(union(names(which(detection_rates_invitro_d30 > 0.01)), names(which(detection_rates_invivo > 0.01))), intersect(names(detection_rates_invitro_d30),names(detection_rates_invivo)))
det_mat <- cbind(seurat@assays$RNA@data[shared_genes_d30,idx_invitro_d30] > 0, seurat_MG@assays$RNA@data[shared_genes_d30,] > 0)
cond <- factor(rep(c("invitro","invivo"), c(length(idx_invitro_d30),ncol(seurat_MG))))
coverage <- c(seurat$nFeature_RNA[idx_invitro_d30], seurat_MG$nFeature_RNA)
ct <- factor(c(as.character(meta_niches$celltype[apply(sim2_MG_niches_bin,1,which.max)[idx_invitro_d30]]),
               as.character(seurat_MG$celltype)),
             levels = levels(seurat_MG$celltype))

DE_invitro_d30_invivo <- binom_group_test(det_mat,
                                          group = cond,
                                          coverage = coverage,
                                          covar = data.frame(ct = ct),
                                          num_threads = 50,
                                          return_coef_group = "invivo")
saveRDS(DE_invitro_d30_invivo, file="res.DE-binom_cocultured-30d_invivo_microglia.rds")
DEG_invitro_d30_invivo <- setdiff(rownames(DE_invitro_d30_invivo),blacklist)[which(p.adjust(DE_invitro_d30_invivo$p_Resi[!rownames(DE_invitro_d30_invivo) %in% blacklist]*2, method="BH")<0.05)]

layout(matrix(1)); par(mar=c(5,5,1,1))
plot(detection_rates_invitro_d30[shared_genes_d30],
     detection_rates_invivo[shared_genes_d30],
     xlab = "Detection rate (in vitro, co-cultured 30d)", ylab = "Detection rate (in vivo)", cex.lab = 1.5, bty = "n", col = "#69696960", pch=16, cex=1)
points(detection_rates_invitro_d30[DEG_invitro_d30_invivo], detection_rates_invivo[DEG_invitro_d30_invivo], col="#303030", bg="red", pch=21, cex=1.5)

plot(log10(detection_rates_invitro_d30[shared_genes_d30]),
     log10(detection_rates_invivo[shared_genes_d30]),
     xlab = "Log10-Detection rate (in vitro, co-cultured 30d)", ylab = "Log10-Detection rate (in vivo)", cex.lab = 1.5, bty = "n", col = "#69696960", pch=16, cex=1)
points(log10(detection_rates_invitro_d30[DEG_invitro_d30_invivo]), log10(detection_rates_invivo[DEG_invitro_d30_invivo]), col="#303030", bg="red", pch=21, cex=1.5)

DEG_TF_invitro_d30_invivo <- sort(intersect(DEG_invitro_d30_invivo, human_TFs$Symbol))
sort(setNames(DE_invitro_d30_invivo[DEG_TF_invitro_d30_invivo,"coef"], DEG_TF_invitro_d30_invivo))

plot(log10(detection_rates_invitro_d30[shared_genes_d30]),
     log10(detection_rates_invivo[shared_genes_d30]),
     xlab = "Log10-Detection rate (in vitro, co-cultured 30d)", ylab = "Log10-Detection rate (in vivo)", cex.lab = 1.5, bty = "n", col = "#90909060", pch=16, cex=1)
points(log10(detection_rates_invitro_d30[DEG_invitro_d30_invivo]), log10(detection_rates_invivo[DEG_invitro_d30_invivo]), bg=ifelse(DE_invitro_d30_invivo[DEG_invitro_d30_invivo,"coef"]>0,"#EC7063","#85C1E9"), col=ifelse(DE_invitro_d30_invivo[DEG_invitro_d30_invivo,"coef"]>0, "#EC7063","#3498DB"), pch=21, cex=1.2)
points(log10(detection_rates_invitro_d30[DEG_TF_invitro_d30_invivo]), log10(detection_rates_invivo[DEG_TF_invitro_d30_invivo]), col="#D4AC0D", bg="#F4D03F", pch=21, cex=1.5)
text(log10(detection_rates_invitro_d30[DEG_TF_invitro_d30_invivo]), log10(detection_rates_invivo[DEG_TF_invitro_d30_invivo])+0.1, cex=0.6, labels=DEG_TF_invitro_d30_invivo)

#### check specifically the microglia core TFs and regulators
genes <- c("RUNX1","SPI1","SALL1","MAFB","MEF2A","MEF2C","IRF8","TGFB1","CSF1","TNF","CSF2")
plot(log10(detection_rates_invitro_d30[shared_genes_d30]),
     log10(detection_rates_invivo[shared_genes_d30]),
     xlab = "Log10-Detection rate (in vitro, co-cultured 30d)", ylab = "Log10-Detection rate (in vivo)", cex.lab = 1.5, bty = "n", col = "#69696960", pch=16, cex=0.6)
points(log10(detection_rates_invitro_d30[DEG_invitro_d30_invivo]), log10(detection_rates_invivo[DEG_invitro_d30_invivo]), bg=ifelse(DE_invitro_d30_invivo[DEG_invitro_d30_invivo,"coef"]>0,"#EC7063","#85C1E9"), col=ifelse(DE_invitro_d30_invivo[DEG_invitro_d30_invivo,"coef"]>0, "#EC7063","#3498DB"), pch=21, cex=0.8)
points(log10(detection_rates_invitro_d30[genes]), log10(detection_rates_invivo[genes]), col="#303030", bg=ifelse(genes %in% DEG_invitro_d30_invivo, ifelse(DE_invitro_d30_invivo[genes,"coef"]>0,"#EC7063","#85C1E9"), "#F9E79F"), pch=21, cex=2)
text(log10(detection_rates_invitro_d30[genes]), log10(detection_rates_invivo[genes])+0.12, labels=genes, cex=1, font=2)


### no-coculture microglia
idx_nocult <- which(! seurat$annot_celltype %in% c("SMC-like","NPC-like","Neuron-like","Neuron-like (phagocytosed)") & seurat$condition == "CTRL")
detection_rates_nocult <- rowMeans(seurat@assays$RNA@data[,idx_nocult] > 0)

shared_genes_nocult <- intersect(union(names(which(detection_rates_nocult > 0.01)), names(which(detection_rates_invivo > 0.01))), intersect(names(detection_rates_nocult),names(detection_rates_invivo)))
det_mat <- cbind(seurat@assays$RNA@data[shared_genes_nocult,idx_nocult] > 0, seurat_MG@assays$RNA@data[shared_genes_nocult,] > 0)
cond <- factor(rep(c("invitro","invivo"), c(length(idx_nocult),ncol(seurat_MG))))
coverage <- c(seurat$nFeature_RNA[idx_nocult], seurat_MG$nFeature_RNA)
ct <- factor(c(as.character(meta_niches$celltype[apply(sim2_MG_niches_bin,1,which.max)[idx_nocult]]),
               as.character(seurat_MG$celltype)),
             levels = levels(seurat_MG$celltype))

DE_nocult_invivo <- binom_group_test(det_mat,
                                     group = cond,
                                     coverage = coverage,
                                     covar = data.frame(ct = ct),
                                     num_threads = 50,
                                     return_coef_group = "invivo")
saveRDS(DE_nocult_invivo, file="res.DE-binom_nococulture_invivo_microglia.rds")
DEG_nocult_invivo <- setdiff(rownames(DE_nocult_invivo),blacklist)[which(p.adjust(DE_nocult_invivo$p_Resi[!rownames(DE_nocult_invivo) %in% blacklist]*2, method="BH")<0.05)]

layout(matrix(1)); par(mar=c(5,5,1,1))
plot(log10(detection_rates_nocult[shared_genes_nocult]),
     log10(detection_rates_invivo[shared_genes_nocult]),
     xlab = "Log10-Detection rate (in vitro, no co-culture)", ylab = "Log10-Detection rate (in vivo)", cex.lab = 1.5, bty = "n", col = "#69696960", pch=16, cex=1)
points(log10(detection_rates_nocult[DEG_nocult_invivo]), log10(detection_rates_invivo[DEG_nocult_invivo]), col="#303030", bg="red", pch=21, cex=1.5)

human_TFs <- read.table("~/Work/databases/animalTFDB/animalTFDB_v3/Homo_sapiens_TF.tsv", sep="\t", header=T)
DEG_TF_nocult_invivo <- sort(intersect(DEG_nocult_invivo, human_TFs$Symbol))
sort(setNames(DE_invitro_invivo[DEG_TF_nocult_invivo,"coef"], DEG_TF_nocult_invivo))

#### check specifically the microglia core TFs and regulators
genes <- c("RUNX1","SPI1","SALL1","MAFB","MEF2A","MEF2C","IRF8","TGFB1","CSF1","TNF","CSF2")
plot(log10(detection_rates_nocult[shared_genes_nocult]),
     log10(detection_rates_invivo[shared_genes_nocult]),
     xlab = "Log10-Detection rate (in vitro, no co-culture)", ylab = "Log10-Detection rate (in vivo)", cex.lab = 1.5, bty = "n", col = "#69696960", pch=16, cex=0.6)
points(log10(detection_rates_nocult[DEG_nocult_invivo]), log10(detection_rates_invivo[DEG_nocult_invivo]), bg=ifelse(DE_nocult_invivo[DEG_nocult_invivo,"coef"]>0,"#EC7063","#85C1E9"), col=ifelse(DE_nocult_invivo[DEG_nocult_invivo,"coef"]>0, "#EC7063","#3498DB"), pch=21, cex=0.8)
points(log10(detection_rates_nocult[genes]), log10(detection_rates_invivo[genes]), col="#303030", bg=ifelse(genes %in% DEG_nocult_invivo, ifelse(DE_nocult_invivo[genes,"coef"]>0,"#EC7063","#85C1E9"), "#F9E79F"), pch=21, cex=2)
text(log10(detection_rates_nocult[genes]), log10(detection_rates_invivo[genes])+0.12, labels=genes, cex=1, font=2)

#### compare with co-cultured microglia
table(sign(DE_nocult_invivo[intersect(rownames(DE_nocult_invivo),rownames(DE_invitro_d30_invivo)),"coef"]),
      sign(DE_invitro_d30_invivo[intersect(rownames(DE_nocult_invivo),rownames(DE_invitro_d30_invivo)),"coef"]))
table(sign(DE_nocult_invivo[intersect(DEG_nocult_invivo,rownames(DE_invitro_d30_invivo)),"coef"]),
      sign(DE_invitro_d30_invivo[intersect(DEG_nocult_invivo,rownames(DE_invitro_d30_invivo)),"coef"]))

all_shared_genes <- intersect(rownames(DE_nocult_invivo), intersect(rownames(DE_invitro_d30_invivo),rownames(DE_invitro_invivo)))
par(mar=c(5,5,1,1))
plot(DE_nocult_invivo[all_shared_genes,"coef"],
     DE_invitro_d30_invivo[all_shared_genes,"coef"],
     cex = 0.8, pch=16, col="#90909060", bty="n",
     xlab="Effect size (No coculture vs. in vivo)", ylab="Effect size (co-cultured 30d vs. in vivo)", cex.lab=1.2)
abline(h=0,v=0,lty=2, col="#636363")
points(DE_nocult_invivo[intersect(DEG_nocult_invivo,all_shared_genes),"coef"],
       DE_invitro_d30_invivo[intersect(DEG_nocult_invivo,all_shared_genes),"coef"],
       cex=1.2, pch=21, col="#303030", bg="red", lwd=0.3)
points(DE_nocult_invivo[intersect(setdiff(DEG_nocult_invivo,DEG_invitro_d30_invivo),all_shared_genes),"coef"],
       DE_invitro_d30_invivo[intersect(setdiff(DEG_nocult_invivo,DEG_invitro_d30_invivo),all_shared_genes),"coef"],
       cex=1.2, pch=21, col="#303030", bg="#F1948A", lwd=0.3)
lines(c(-30,30),c(-30,30), lty=3, col="#636363", lwd=2)

par(mar=c(8,5,1,1))
boxplot(list(abs(DE_nocult_invivo[intersect(DEG_nocult_invivo,all_shared_genes),"coef"])-abs(DE_invitro_invivo[intersect(DEG_nocult_invivo,all_shared_genes),"coef"]),
             abs(DE_nocult_invivo[intersect(DEG_nocult_invivo,all_shared_genes),"coef"])-abs(DE_invitro_d30_invivo[intersect(DEG_nocult_invivo,all_shared_genes),"coef"]),
	     abs(DE_invitro_invivo[intersect(DEG_nocult_invivo,all_shared_genes),"coef"])-abs(DE_invitro_d30_invivo[intersect(DEG_nocult_invivo,all_shared_genes),"coef"])),
	names = c("No vs. Co-all","No vs. Co-30d","Co-all vs. Co-30d"), cex.axis = 1,
        ylab="DEG effect size changes", cex.lab = 1.5, xlab=NA, frame=F, pch=16, las=2)
abline(h=0, lty=2, col="#909090")


# coculture-duration-related expression changes
idx_invitro <- which(! seurat$annot_celltype %in% c("SMC-like","NPC-like","Neuron-like","Neuron-like (phagocytosed)") & seurat$condition == "COCUL")
detection_rates_invitro <- rowMeans(seurat@assays$RNA@data[,idx_invitro] > 0)
gene_cands <- names(which(detection_rates_invitro>0))

DE_cocul_duration <- ancova_group_test(expr = seurat@assays$RNA@data[gene_cands,idx_invitro],
                                       group = as.numeric(gsub("D","",seurat$timepoint[idx_invitro])),
                                       covar = data.frame(coverage = log(seurat$nCount_RNA[idx_invitro]), ct=seurat$annot_celltype[idx_invitro],row.names=colnames(seurat)[idx_invitro]),
                                       num_threads = 50)
DE_cocul_duration_noct <- ancova_group_test(expr = seurat@assays$RNA@data[gene_cands,idx_invitro],
                                            group = as.numeric(gsub("D","",seurat$timepoint[idx_invitro])),
                                            covar = data.frame(coverage = log(seurat$nCount_RNA[idx_invitro]), row.names=colnames(seurat)[idx_invitro]),
                                            num_threads = 50)
DEG_cocult <- setdiff(rownames(DE_cocul_duration_noct)[which(p.adjust(DE_cocul_duration_noct$p_ANOVA, method="bonferroni")<0.01)], blacklist)
