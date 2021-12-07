rm(list=ls())
.ArgList=.myArgListCreatorFn(cellClass = "Astro")
#exNonMicCells=NULL;argList=.ArgList;convert_to_gene_level=F;returnDataM=F
data=.myReadDataFn_cellClass(exNonMicCells=NULL,argList=.ArgList,convert_to_gene_level=F,returnDataM=F)
data_grouped=.myReadDataFn_cellClass(exNonMicCells=NULL,argList=.ArgList,convert_to_gene_level=F,returnDataM=F,breakSubjects=F)
names(data)
fd=.myGeneAnnoFn(inputData = data)
cell_phenoData=.myCellPhenoDataFn(inputArg = .ArgList,inputData = data)
data=cell_phenoData$data
pd=cell_phenoData$pd
rm(cell_phenoData)
cell_phenoData=.myCellPhenoDataFn(inputArg = .ArgList,inputData = data_grouped)
data_grouped=cell_phenoData$data
pd_grouped=cell_phenoData$pd
rm(cell_phenoData)
rm(list=ls())
.ArgList=.myArgListCreatorFn(cellClass = "Astro")
#exNonMicCells=NULL;argList=.ArgList;convert_to_gene_level=F;returnDataM=F
data=.myReadDataFn_cellClass(exNonMicCells=NULL,argList=.ArgList,convert_to_gene_level=F,returnDataM=F)
data_grouped=.myReadDataFn_cellClass(exNonMicCells=NULL,argList=.ArgList,convert_to_gene_level=F,returnDataM=F,breakSubjects=F)
names(data)
fd=.myGeneAnnoFn(inputData = data)
cell_phenoData=.myCellPhenoDataFn(inputArg = .ArgList,inputData = data)
data=cell_phenoData$data
pd=cell_phenoData$pd
rm(cell_phenoData)
cell_phenoData=.myCellPhenoDataFn(inputArg = .ArgList,inputData = data_grouped)
data_grouped=cell_phenoData$data
pd_grouped=cell_phenoData$pd
rm(cell_phenoData)
data_seurat=.myHighVarGeneSlFn(data,dataorganism="Human",argList = .ArgList)
data_grouped=.myHighVarGeneSlFn(data_grouped,dataorganism="Human",argList = .ArgList)
p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = "anno_cluster_res")
p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = "final_anno")
ggsave(plot=p,file="~/myBucket/torm.pdf",width=12,height=12)
res_delist=.myDEanalysisFn(inputPd=refined_anno,inputFd=fd,batch_names=batch_names,inputData=data,cluster_col="anno_cluster_res",do.pval=F)
names(res_delist)
p=res_delist$umap_plot
p
#######################
#Marker analysis module
#seurat_clustering: anno_cluster_res
#inph clusters: inferred_prop_iNPH_lbl
res_delist=.myDEanalysisFn(inputPd=refined_anno,inputFd=fd,batch_names=batch_names,inputData=data,cluster_col="anno_cluster_res",do.pval=F)
names(res_delist)
p=.my2dPlotArrangerFn(inputPd=pd,inputGene="Gfap",inputExpData=data_seurat$data,argList=.ArgList,ncores=10,FeaturePlotOnly=F)
p
ggsave(plot=p,file="~/myBucket/torm2.pdf",height=20,width=20)
res_delist=.myDEanalysisFn(inputPd=refined_anno,inputFd=fd,batch_names=batch_names,inputData=data,cluster_col="anno_cluster_res",do.pval=F)
names(res_delist)
p=res_delist$umap_plot
p
res_delist=.myDEanalysisFn(inputPd=refined_anno,inputFd=fd,inputData=data,cluster_col="anno_cluster_res",do.pval=F)
names(res_delist)
p=res_delist$umap_plot
p
#######################
#Marker analysis module
#seurat_clustering: anno_cluster_res
#inph clusters: inferred_prop_iNPH_lbl
res_delist=.myDEanalysisFn(inputPd=refined_anno,inputFd=fd,inputData=data,cluster_col="anno_cluster_res",do.pval=F)
#######################
#Marker analysis module
#seurat_clustering: anno_cluster_res
#inph clusters: inferred_prop_iNPH_lbl
res_delist=.myDEanalysisFn(inputPd=pd,inputFd=fd,inputData=data,cluster_col="anno_cluster_res",do.pval=F)
pd[1:4,1:4]
pd[1:4, 1:4]
fd[1:4, 1:4]
dim(fd)
fd[1:4, 1:3]
class(data)
class(data[[1]])
length(data)
names(data)
data$summary_statistic
length(data$data)
tmp = unlist(lapply(data$data, function(x) {x$ds_batch[1]}))
table(tmp)
inph_data = data$data[tmp="human_Tushar_iNPH"]
length(inph_data)
inph_data = data$data[tmp=="human_Tushar_iNPH"]
length(inph_data)
inph_data=.mycBindFn(inph_data)
dim(inph_data)
dims(inph_data)
class(inph_data)
class(inph_data[[1]])
names(inph_data)
tst=.mycBindFn(inputList = inph_data,verbose = T)
dim(tst)
class(tst)
pd2=colData(tst)
dim(pd2)
pd2[1:4,1:4]
pd2[1:4,1:10]
fd=rowData(tst)
dim(fd)
head(fd)
tail(fd)
expData=counts(tst)
expData[1:4,1:4]
head(row.names(fd))
head(colnames(expData))
head(row.names(pd))
head(pd$sample)
head(pd2)
head(row.names(pd2))
head(colnames(expData))
all(row.names(pd2)==colnames(expData))
table(row.names(pd2)==colnames(expData))
table(row.names(pd2)!=colnames(expData))
head(row.names(pd2)==colnames(expData))
class(expData)
class(tst)
tst_seurat=.extraExport2SeuratFn(tst)
class(tst_seurat)
.myDimPlotFn(pd,attCol = "final_anno")
.myDimPlotFn(pd,dimCols = c("UMAP_2","UMAP_1"),attCol = "final_anno")
Seurat::DimPlot(tst_seurat)
slotNames(tst_seurat)
names(tst_seurat@reductions)
pd[1:4,1:4]
head(pd$final_anno)
head(pd$UMAP_1)
head(pd$UMAP_2)
.myDimPlotFn(pd,dimCols = c("UMAP_2","UMAP_1"),attCol = "ds_batch")
history(100)
