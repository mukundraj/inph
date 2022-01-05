# Script for converting processing integrated data through Seurat pipeline
# Load initial function definitions from code_mukund.R

data=.myReadDataFn_cellClass(exNonMicCells=NULL,argList=.ArgList,convert_to_gene_level=F,returnDataM=F)

# misc commands start
class(data[[1]][[1]])
rownames(rowData(data[[1]][[1]]))
rownames(rowData(data[[1]][[1]]))
colnames(colData(data[[1]][[1]]))
# misc commands ends

tmp = unlist(lapply(data$data, function(x) {x$ds_batch[1]}))
inph_data = data$data[tmp=="human_Tushar_iNPH"]
tst=.mycBindFn(inputList = inph_data,verbose = T)

# misc commands part2

rowData(tst)
colData(tst)
colData(tst)$clusters

colnames(colData(tst))
table(colData(tst)$ds_batch)
table(colData(tst)$clusters)
table(colData(tst)$anno_sex)
table(colData(tst)$anno_tissue)
table(colData(tst)$status)
table(colData(tst)$dataset)

# misc commands part2 ends



tst_seurat=.extraExport2SeuratFn(tst)

# misc commands part3
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final
SaveH5Seurat(pbmc3k.final, filename = "pbmc3k.h5Seurat")
Convert("pbmc3k.h5Seurat", dest = "h5ad")

## python
# >>> import zarr
# >>> z2 = zarr.open('my_store.zarr', mode='r')
# >>> z2.tree()
# >>> z2["/obs/seurat_clusters"]
# >>> z2.info
# >>> z2["/obs/seurat_annotations"].info
# >>> z2.obs.seurat_clusters
# >>> z2.obs.seurat_clusters[:5]


## SaveH5Seurat(tst_seurat, filename = "tst_seurat.h5Seurat") # currently throwing error

# misc commands part3 ends



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
tst_seurat[["percent.mt"]] <- PercentageFeatureSet(tst_seurat, pattern = "^MT-")


# Visualize QC metrics as a violin plot
VlnPlot(tst_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(tst_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tst_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


tst_seurat <- subset(tst_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

tst_seurat <- NormalizeData(tst_seurat, normalization.method = "LogNormalize", scale.factor = 10000)



tst_seurat <- FindVariableFeatures(tst_seurat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(tst_seurat), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(tst_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


all.genes <- rownames(tst_seurat)
tst_seurat <- ScaleData(tst_seurat, features = all.genes)

tst_seurat <- RunPCA(tst_seurat, features = VariableFeatures(object = tst_seurat))

print(tst_seurat[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(tst_seurat, dims = 1:2, reduction = "pca")

DimPlot(tst_seurat, reduction = "pca")

DimHeatmap(tst_seurat, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(tst_seurat, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
# tst_seurat <- JackStraw(tst_seurat, num.replicate = 100)
# tst_seurat <- ScoreJackStraw(tst_seurat, dims = 1:20)

# JackStrawPlot(tst_seurat, dims = 1:15)

ElbowPlot(tst_seurat)


# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
tst_seurat <- RunUMAP(tst_seurat, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(tst_seurat, reduction = "umap")

saveRDS(tst_seurat, file = "./inph_tutorial.rds")

### clustering

tst_seurat <- FindNeighbors(tst_seurat, dims = 1:10)
tst_seurat <- FindClusters(tst_seurat, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(tst_seurat), 5)


###


# find all markers of cluster 2
cluster2.markers <- FindMarkers(tst_seurat, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(tst_seurat, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

## Useful

# For SingleCellExperiment class
# https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html