pbmc.data <- Read10X(data.dir = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/Ctrl_FC/new_reference/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "Ctrl_FC", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt:")
pbmc[["percent.rp"]] <- PercentageFeatureSet(pbmc, pattern = "^Rp")
pbmc <- subset(pbmc, subset = nCount_RNA < 12000 & nFeature_RNA < 1800 & nFeature_RNA > 650 & percent.mt < 10)
pbmc <- NormalizeData(pbmc)
cc.genes <- readLines(con = "/home/dchatterjee/cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68]
s.genes <- cc.genes[69:124]
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc$CC.Difference <- pbmc$S.Score - pbmc$G2M.Score
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pbmc <- FindNeighbors(pbmc, dims = 1:60)
pbmc <- FindClusters(pbmc, resolution = 1)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:60, n.neighbors = 40, min.dist = 0.25)
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = TRUE)
pbmc <- subset(pbmc, idents = c(19,17,8,13,15,9), invert = T)

pbmc -> ctrl

ldat <- ReadVelocity(file = "/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglIR/new_reference/LglIR.loom")
pbmc <- as.Seurat(x = ldat)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt:")
pbmc[["percent.rp"]] <- PercentageFeatureSet(pbmc, pattern = "^Rp")
pbmc <- subset(pbmc, subset = nCount_spliced < 15000 & nFeature_spliced < 1800 & nFeature_spliced > 600 & percent.mt < 4.5)
pbmc <- NormalizeData(pbmc)
cc.genes <- readLines(con = "/home/dchatterjee/cell_cycle_genes.txt")
g2m.genes <- cc.genes[1:68]
s.genes <- cc.genes[69:124]
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc$CC.Difference <- pbmc$S.Score - pbmc$G2M.Score
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_spliced"), verbose = TRUE)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pbmc <- FindNeighbors(pbmc, dims = 1:60)
pbmc <- FindClusters(pbmc, resolution = 1)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:60, n.neighbors = 40, min.dist = 0.20)
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = TRUE)

pbmc <- subset(pbmc, idents = c(11,8,16,12), invert = T)
pbmc@assays$spliced -> pbmc@assays$RNA
pbmc@meta.data$nFeature_spliced -> pbmc@meta.data$nFeature_RNA
pbmc@meta.data$nCount_spliced -> pbmc@meta.data$nCount_RNA
pbmc@meta.data$orig.ident <- "LglIR_FC"

pbmc -> stim

all.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), anchor.features = 2000, dims = 1:50)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:50)
all.combined <- ScaleData(all.combined, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
all.combined -> pbmc

pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pbmc <- FindNeighbors(pbmc, dims = 1:60)
pbmc <- FindClusters(pbmc, resolution = 1.5)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:60, n.neighbors = 40, min.dist = 0.20)
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = F, group.by = "orig.ident")
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = T)

pbmc <- subset(pbmc, idents = c(24), invert = T)

split.pbmc <- SplitObject(pbmc, split.by = "orig.ident")
split.pbmc
ctrl <- split.pbmc$Ctrl_FC
stim <- split.pbmc$LglIR_FC
DefaultAssay(ctrl) <- "RNA"
DefaultAssay(stim) <- "RNA"
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)
ctrl <- NormalizeData(ctrl)
stim <- NormalizeData(stim)
all.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), anchor.features = 2000, dims = 1:50)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:50)
DefaultAssay(all.combined) <- "integrated"
all.combined <- ScaleData(all.combined, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), verbose = TRUE)
all.combined -> pbmc
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), npcs = 100, nfeature.print = 10, ndims.print = 1:5, verbose = T)
pbmc <- FindNeighbors(pbmc, dims = 1:60)
pbmc <- FindClusters(pbmc, resolution = 2.5)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:60, n.neighbors = 25, min.dist = 0.20)
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = TRUE)

pbmc <- FindNeighbors(pbmc, dims = 1:60)
pbmc <- FindClusters(pbmc, resolution = 1.5)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:60, n.neighbors = 25, min.dist = 0.20)
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = TRUE)
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = F, group.by = "orig.ident")

pbmc <- FindNeighbors(pbmc, dims = 1:60)
pbmc <- FindClusters(pbmc, resolution = 2)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:60, n.neighbors = 25, min.dist = 0.20)
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = TRUE)

pbmc <- FindNeighbors(pbmc, dims = 1:60)
pbmc <- FindClusters(pbmc, resolution = 1.75)
table(pbmc@active.ident,pbmc@meta.data$orig.ident)
pbmc <- RunUMAP(pbmc, dims = 1:60, n.neighbors = 25, min.dist = 0.20)
DimPlot(pbmc, reduction = "umap", pt.size = 1.2, label = TRUE)


#LOOM and SCVELO PREPARATION
library(SeuratDisk)
pbmc <- readRDS("/home/dchatterjee/Documents/Lgl_Upd_NICD_FC/LglKD+CtrlWT_FC/new_reference/Reproducible/CtrlFC_veloLglIRFC_ONLYFC.rds")

DefaultAssay(pbmc) <- "RNA"
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(pbmc)
DefaultAssay(pbmc) <- "spliced"
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(pbmc)
DefaultAssay(pbmc) <- "unspliced"
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(pbmc)
DefaultAssay(pbmc) <- "ambiguous"
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(pbmc)
DefaultAssay(pbmc) <- "integrated"
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(pbmc)
DefaultAssay(pbmc) <- "RNA"
pbmc.loom <- as.loom(pbmc, filename = "/home/cluster7ONLY_vel.loom", verbose = TRUE, overwrite = TRUE)
pbmc.loom

SaveH5Seurat(pbmc, filename = "lgl.h5Seurat", overwrite = TRUE)
Convert("lgl.h5Seurat", dest = "h5ad", overwrite = TRUE)

pbmc.loom$close_all()

#IN PYTHON
import scvelo as scv
adata = scv.read("lgl.h5ad")
adata

#scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
#scv.pp.moments(adata, n_pcs=60, n_neighbors=25)
#scv.tl.velocity(adata)
scv.tl.velocity(adata, mode='steady_state', min_r2=1)
print(adata.var['velocity_genes'].sum(), adata.n_vars, adata.n_obs)

scv.pl.proportions(adata, groupby='seurat_clusters')

scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters")
scv.pl.velocity_embedding(adata, basis="umap", color="seurat_clusters", arrow_length=3, arrow_size=2, dpi=120)
scv.pl.velocity_graph(adata, threshold=5, color='seurat_clusters')
scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])
scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(50)


scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters')
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', cmap='coolwarm', size=80)
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100, cmap='coolwarm')
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, color='seurat_clusters', frameon=False)
scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(20)

var_names = ['Mmp1', 'Keap1', 'cnc', 'sn']
scv.pl.scatter(adata, var_names, color='seurat_clusters', frameon=True)
#scv.pl.scatter(adata, x='latent_time', y=var_names, color='seurat_clusters', frameon=True)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
scv.pl.velocity_embedding_stream(adata, basis='umap', color='seurat_clusters')
fig.savefig('/home/dchatterjee/myimage.svg', format='svg', dpi=1200)

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='seurat_clusters', n_convolve=100)

top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False, color='seurat_clusters')

scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
df.head(20)

scv.tl.differential_kinetic_test(adata, groupby='seurat_clusters')
df = scv.get_df(adata, ['fit_diff_kinetics', 'fit_pval_kinetics'], precision=2)
df.head(20)
kwargs = dict(linewidth=2, add_linfit=True, frameon=False)
scv.pl.scatter(adata, basis=['Ets21C'], add_outline='fit_diff_kinetics', **kwargs)

#
diff_clusters=list(adata.var['fit_diff_kinetics'])
scv.pl.scatter(adata, legend_loc='right', size=60, title='diff kinetics', add_outline=diff_clusters, outline_width=(.8, .2))




