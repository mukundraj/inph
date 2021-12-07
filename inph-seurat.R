
# Load initial function definitions from code_mukund.R

data=.myReadDataFn_cellClass(exNonMicCells=NULL,argList=.ArgList,convert_to_gene_level=F,returnDataM=F)

tmp = unlist(lapply(data$data, function(x) {x$ds_batch[1]}))
inph_data = data$data[tmp=="human_Tushar_iNPH"]
tst=.mycBindFn(inputList = inph_data,verbose = T)
tst_seurat=.extraExport2SeuratFn(tst)



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
