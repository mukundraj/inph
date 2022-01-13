# Script for converting loading inph data and converting to anndata format 
# (for input to python script to eventually convert to zarr format)
# Loads initial function definitions from code_mukund.R
# Instruction: run scripts from git repo's root directory
# Started on 13 Jan 2022 by Mukund

# load functions from code_mukund.R
source('base_code_vahid.R')

# clear env
rm(list=ls())

# loading the data .. here Astro.. (faster) / another e.g. Oligo
.ArgList=.myArgListCreatorFn(cellClass = "Astro")

data=.myReadDataFn_cellClass(exNonMicCells=NULL,argList=.ArgList,convert_to_gene_level=F,returnDataM=F)
tmp = unlist(lapply(data$data, function(x) {x$ds_batch[1]}))
inph_data = data$data[tmp=="human_Tushar_iNPH"]
tst=.mycBindFn(inputList = inph_data,verbose = T) # takes time to complete
tst_seurat=.extraExport2SeuratFn(tst)

# Edit Vahid's seurat object to conform to requirement of sceasy library
tst_seurat@assays$RNA@meta.features <- subset (tst_seurat@assays$RNA@meta.features, select = -entrezid)

# Compute the coordinates and clustering info using Seurat pipeline

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
tst_seurat <- RunUMAP(tst_seurat, dims = 1:10)

# Convert seurat object to anndata

library(sceasy)
library(reticulate)
sceasy::convertFormat(tst_seurat, from="seurat", to="anndata",
                      outFile='output/tmp-anndata.h5ad')

# Call python script to convert anndata to zarr
system('python scripts/anndata_to_zarr.py ./output/tmp-anndata.h5ad ./output/my-zarr-store.zarr')

# Call gsutil function to copy to google bucket