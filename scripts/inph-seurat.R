# Script for processing iNPH data through the Seurat pipeline

library(qs)

# Data downloaded from gs://macosko_data/vgazesta/SC_data/Human/brain/snRNA/Tushar_iNPH
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