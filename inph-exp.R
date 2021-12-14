

library(qs)
data<-qread("/home/mraj/Desktop/work/data/Tushar_iNPH/Astro_data_arranged_updatedId_final_batches.qs")

dim(assay(data))

metaData <- colData(data)
dim(metaData)
head(metaData) ## tushar meet
table(metaData$clusters)

pData <- rowData(data)
dim(pData)


tst_seurat=.extraExport2SeuratFn(data)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
tst_seurat[["percent.mt"]] <- PercentageFeatureSet(tst_seurat, pattern = "^MT-")

tst_seurat <- subset(tst_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Visualize QC metrics as a violin plot
VlnPlot(tst_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(tst_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(tst_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
  

## --- converting seurat to vitessce json starts ---
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

wrapped_seurat <- SeuratWrapper$new(tst_seurat, assay = "RNA")
wrapped_seurat$convert_and_save(0,0)
print (wrapped_seurat$get_out_dir(0,0))

cells_list <- wrapped_seurat$create_cells_list()
cells_json <- jsonlite::toJSON(cells_list)
write(cells_json, file = "inph_data.json")

#scp mraj@34.136.224.248:/tmp/RtmpFFPrAk/0/0/cells.json /Users/mraj/Desktop/work/data/testportal_data
#scp mraj@34.136.224.248:/tmp/RtmpFFPrAk/0/0/cell-sets.json /Users/mraj/Desktop/work/data/testportal_data
#scp mraj@34.136.224.248:/tmp/RtmpFFPrAk/0/0/expression-matrix.json /Users/mraj/Desktop/work/data/testportal_data

## --- converting seurat to vitessce json ends ---