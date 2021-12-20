# Script for converting Seurat object to Vitessce json 

# Adapted from https://github.com/vitessce/vitessce-r/blob/HEAD/R/wrappers.R

# clear env
rm(list=ls())

library(Seurat)
library(SeuratObject)
install.packages("devtools")
devtools::install_github("vitessce/vitessce-r")
library(vitessce)


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
tst_seurat[["percent.mt"]] <- PercentageFeatureSet(tst_seurat, pattern = "^MT-")

tst_seurat <- subset(tst_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

tst_seurat <- NormalizeData(tst_seurat)
tst_seurat <- FindVariableFeatures(tst_seurat, selection.method = "vst", nfeatures = 2000)



all.genes <- rownames(tst_seurat)
tst_seurat <- ScaleData(tst_seurat, features = all.genes)


tst_seurat <- RunPCA(tst_seurat, features = VariableFeatures(object = tst_seurat))
ElbowPlot(tst_seurat)

tst_seurat <- FindNeighbors(tst_seurat, dims = 1:10)
tst_seurat <- FindClusters(tst_seurat, resolution = 0.5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
tst_seurat <- RunUMAP(tst_seurat, dims = 1:10)
DimPlot(tst_seurat, reduction = "umap")

wrapped_seurat <- SeuratWrapper$new(tst_seurat, assay = "RNA", cell_set_meta_names="seurat_clusters")
wrapped_seurat$convert_and_save(0,0)
print (wrapped_seurat$get_out_dir(0,0))

cells_list <- wrapped_seurat$create_cells_list()
cells_json <- jsonlite::toJSON(cells_list)
write(cells_json, file = "inph_data.json")

# scp mraj@34.66.221.119:/tmp/RtmpPC8JjO/0/0/cells.json /Users/mraj/Desktop/work/data/testportal_data
# scp mraj@34.66.221.119:/tmp/RtmpPC8JjO/0/0/cell-sets.json /Users/mraj/Desktop/work/data/testportal_data
# scp mraj@34.66.221.119:/tmp/RtmpPC8JjO/0/0/expression-matrix.json /Users/mraj/Desktop/work/data/testportal_data
## --- converting seurat to vitessce json ends ---

## misc commands

colnames(tst_seurat@meta.data)