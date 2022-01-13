# Seprated base code that was initially in code_mukund.R from Vahid to be able 
# call from other scripts without having to also run data loading and analysis.
# Started on 13 Jan 2022 by Mukund

library(googleCloudStorageR)
library(qs)
library(future)
library(googleAuthR)


# gcs_source("vgazesta/code/sconline_code.R", bucket="macosko_data")
# gcs_source("vgazesta/code/Human_mouse_alignmnet_code.R")

source("sconline_code.R")
source("Human_mouse_alignmnet_code.R")


plan("multicore", workers = 8)
options(future.globals.maxSize = 1000 * 1024^4)

#inputPd=pd;inputFd=fd;inputData=data;cluster_col="anno_cluster_res";do.pval=F;selected_clusters=NULL;contrast_clusters=NULL
.myDEanalysisFn=function(inputPd,inputFd=NULL,inputData,cluster_col,do.pval=F,selected_clusters=NULL,contrast_clusters=NULL,ncores=10){
  
  myLocalDEfn=function(iclust,tmp,contrast_clusters){
    if(!is.null(contrast_clusters)){
      tmp2=.myEvalMarkers(object=tmp, cells.1=colnames(tmp)[Idents(tmp)==iclust], cells.2=colnames(tmp)[!tolower(as.character(Idents(tmp))) %in% tolower(c(contrast_clusters,iclust))], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=iclust)
    } else {
      tmp2=.myEvalMarkers(object=tmp, cells.1=colnames(tmp)[Idents(tmp)==iclust], cells.2=colnames(tmp)[Idents(tmp)!=iclust], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=iclust)
    }
    
    tmp2$gene=row.names(tmp2)
    return(tmp2)
  }
  
  if(!is.null(contrast_clusters)&do.pval){
    warning("do.pval is set to FALSE")
    do.pval=F
  }
  
  batch_names=unlist(lapply(inputData$data,function(x) x$ds_batch[1]))
  
  umap_plot=.myDimPlotFn(object=inputPd, dimCols = c("UMAP_1", "UMAP_2"), attCol = cluster_col)
  
  res_DE_analysis=list()
  for(ids in unique(batch_names)){
    tmp=inputData$data[which(batch_names==ids)]
    tmp_cols=unlist(lapply(tmp,function(x) colnames(x)))
    if(sum(inputPd$sample %in% tmp_cols)>0){
      
      if(length(tmp)>1){
        tmp=.mycBindFn(inputList=tmp,verbose=T)
      } else {
        tmp=tmp[[1]]
      }
      
      tmp_pd=inputPd[match(colnames(tmp),inputPd$sample),]
      tmp$cluster=tmp_pd[,cluster_col]
      tmp$inferred_lbl=tmp_pd$inferred_prop_iNPH_lbl
      tmp=tmp[,!is.na(tmp$cluster)]
      if(class(tmp$cluster)==class("")){
        tmp$cluster=as.factor(tmp$cluster)
      }
      tmp$cluster=droplevels(tmp$cluster)
      tmp=.extraExport2SeuratFn(inputData=tmp,project="scRNA")
      tmp=Seurat::NormalizeData(tmp)
      Idents(tmp)=tmp$cluster
      all_genes_res=list()
      tmp_sl_clusters=unique(tmp$cluster)
      if(!is.null(selected_clusters)){
        tmp_sl_clusters=tmp_sl_clusters[which(tolower(tmp_sl_clusters) %in% tolower(selected_clusters))]
      }
      if(F){
        for(iclust in unique(tmp$cluster)){
          if(!is.null(contrast_clusters)){
            tmp2=.myEvalMarkers(object=tmp, cells.1=colnames(tmp)[Idents(tmp)==iclust], cells.2=colnames(tmp)[!tolower(as.character(Idents(tmp))) %in% tolower(c(contrast_clusters,iclust))], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=iclust)
          } else {
            tmp2=.myEvalMarkers(object=tmp, cells.1=colnames(tmp)[Idents(tmp)==iclust], cells.2=colnames(tmp)[Idents(tmp)!=iclust], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=iclust)
          }
          
          tmp2$gene=row.names(tmp2)
          all_genes_res=c(all_genes_res,list(tmp2))
        }
      } else {
        all_genes_res=parallel::mclapply(unique(tmp$cluster),myLocalDEfn,tmp=tmp,contrast_clusters=contrast_clusters,mc.cores=ncores)
      }
      
      all_genes_res=do.call("rbind",all_genes_res)
      
      if(do.pval){
        tmp=Seurat::FindAllMarkers(tmp,logfc.threshold = 0.1,min.diff.pct =0.3)
        tmp=merge(all_genes_res,tmp[,c("cluster","gene","p_val_adj")],by.x=c("cluster_name","gene"),by.y=c("cluster","gene"),all.x=T)
      } else {
        tmp=all_genes_res
      }
      
      
      tmp$dataset=ids
      if(!is.null(inputFd)){
        tmp_fd=inputFd[match(tmp$gene,inputFd$rowname),]
        #tst=dplyr::left_join(tmp,inputFd[,c("rowname","gene_short_name")],by.x="gene",by.y="rowname")
        tmp$gene_short_name=tmp_fd$gene_short_name#=merge(tmp,inputFd[,c("rowname","gene_short_name")],by.x="gene",by.y="rowname",all.x=T)
      }
      res_DE_analysis=c(res_DE_analysis,list(new=tmp))
      names(res_DE_analysis)[names(res_DE_analysis)=="new"]=ids
      
    }
    
  }
  res_delist=res_DE_analysis
  res_delist_m=do.call("rbind",res_delist)
  res_delist_m$expressed=0
  res_delist_m$expressed[res_delist_m$cell.count.1>10]=1
  res_delist_m$is.de=0
  res_delist_m$is.de[res_delist_m$avg_logFC>0.1&(res_delist_m$pct.1 - res_delist_m$pct.2)>0.1&res_delist_m$expressed==1]=1
  res_delist_m$organism=unlist(lapply(strsplit(res_delist_m$dataset,"_"),function(x) x[1]))
  
  res_delist_expressed_human=aggregate(expressed~cluster_name+gene_short_name+gene,data=res_delist_m[which(res_delist_m$organism=="human"),],sum)
  res_delist_expressed_mouse=aggregate(expressed~cluster_name+gene_short_name+gene,data=res_delist_m[which(res_delist_m$organism=="mouse"),],sum)
  res_delist_expressed_human$organism="human"
  res_delist_expressed_mouse$organism="mouse"
  res_delist_expressed=rbind(res_delist_expressed_human,res_delist_expressed_mouse)
  
  res_delist_de_human=aggregate(is.de~cluster_name+gene_short_name+gene,data=res_delist_m[which(res_delist_m$organism=="human"),],sum)
  res_delist_de_mouse=aggregate(is.de~cluster_name+gene_short_name+gene,data=res_delist_m[which(res_delist_m$organism=="mouse"),],sum)
  res_delist_de_human$organism="human"
  res_delist_de_mouse$organism="mouse"
  res_delist_de=rbind(res_delist_de_human,res_delist_de_mouse)
  
  res_delist_arranged=merge(res_delist_de,res_delist_expressed,by=c("cluster_name","gene_short_name","gene","organism"),all=T)
  if(length(unique(res_delist_arranged$organism))>1){
    res_delist_arranged=split(res_delist_arranged,res_delist_arranged$organism)
    colnames(res_delist_arranged[[1]])[colnames(res_delist_arranged[[1]]) %in% c("expressed","is.de")]=paste0(names(res_delist_arranged)[1],"_",colnames(res_delist_arranged[[1]])[colnames(res_delist_arranged[[1]]) %in% c("expressed","is.de")])
    colnames(res_delist_arranged[[2]])[colnames(res_delist_arranged[[2]]) %in% c("expressed","is.de")]=paste0(names(res_delist_arranged)[2],"_",colnames(res_delist_arranged[[2]])[colnames(res_delist_arranged[[2]]) %in% c("expressed","is.de")])
    res_delist_arranged=merge(res_delist_arranged[[1]],res_delist_arranged[[2]],by=c("cluster_name","gene_short_name","gene"),all=T)
  }
  res_delist_arranged=res_delist_arranged[,!grepl("organism",colnames(res_delist_arranged))]
  
  inph_data=res_delist_m[res_delist_m$dataset=="human_Tushar_iNPH",]
  colnames(inph_data)[!colnames(inph_data) %in% c("gene","gene_short_name","cluster_name")]=paste0("iNPH_",colnames(inph_data)[!colnames(inph_data) %in% c("gene","gene_short_name","cluster_name")])
  
  res_delist_arranged=merge(inph_data[,!colnames(inph_data) %in% c("iNPH_organism","iNPH_dataset")],res_delist_arranged,by=c("gene","gene_short_name","cluster_name"),all=T)
  if(sum(colnames(res_delist_arranged)=="human_expressed")>0){
    res_delist_arranged=res_delist_arranged[order(res_delist_arranged$human_is.de/res_delist_arranged$human_expressed,res_delist_arranged$human_is.de,decreasing = T),]
  } else if(sum(colnames(res_delist_arranged)=="mouse_expressed")>0){
    res_delist_arranged=res_delist_arranged[order(res_delist_arranged$mouse_is.de/res_delist_arranged$mouse_expressed,res_delist_arranged$mouse_is.de,decreasing = T),]
  } else if(sum(colnames(res_delist_arranged)=="expressed")>0){
    res_delist_arranged=res_delist_arranged[order(res_delist_arranged$is.de/res_delist_arranged$expressed,res_delist_arranged$is.de,decreasing = T),]
  }
  
  
  return(list(res_delist_m=res_delist_m,res_delist_arranged=res_delist_arranged,umap_plot=umap_plot))
}

.myArgListCreatorFn=function(cellClass){
  if(cellClass %in% c("MG","Astro","InN","ExN","OPC","Endo")){
    .ArgList=qread(paste0("~/myBucket/integrative_Tushar/SC_data/meta_analysis/",cellClass,"_ArgList_HVGlist_IndScaling_Feb16_MGIsymbol_Final.qs"))
  } else if(cellClass=="Oligo"){
    .ArgList=qread("~/myBucket/integrative_Tushar/SC_data/meta_analysis/Oligo_ArgList_HVGlist_IndScaling_Feb16_MGIsymbol_correct_nocovariate.qs")
  } else if(cellClass=="Astro_Human"){
    .ArgList=qread("~/myBucket/integrative_Tushar/SC_data/meta_analysis/Astro_Human_ArgList_HVGlist_IndScaling_Feb16_MGIsymbol_correct_nocovariate.qs")
    
  } else if(cellClass=="Astro_Projection"){
    .ArgList=qread(paste0("~/myBucket/SC_data/meta_analysis/Astro_HumanMouse_ArgList_HVGlist","_IndScaling_Feb16_MGIsymbol_projection.qs"))
    
  } else {
    stop("Unrecoginzed input cellClass")
  }
  
  if(sum(names(.ArgList)=="do.split.prop")==0){
    .ArgList$do.split.prop=T
  }
  
  .ArgList$prop.n.neighbors=4
  
  
  if(sum(names(.ArgList)=="Leng200")==0){
    .ArgList$Leng200=F
  }
  
  if(sum(names(.ArgList)=="pseudocell_selection_method")==0){
    .ArgList$pseudocell_selection_method="kmeans"
  }
  
  return(.ArgList)
}
.myGeneAnnoFn=function(inputData){
  fd=as.data.frame(rowData(inputData$data[[1]]))
  fd=fd[,c("ID","gene_short_name")]
  fd$rowname=row.names(fd)
  for(i in 2:length(inputData$data)){
    tmp=as.data.frame(rowData(inputData$data[[i]]))
    tmp=tmp[,c("ID","gene_short_name")]
    tmp$rowname=row.names(tmp)
    fd=rbind(fd,tmp)
    fd=fd[!duplicated(fd$rowname),]
  }
  
  return(fd)
}
.myCellPhenoDataFn=function(inputArg,inputData){
  
  if(sum(names(inputArg)=="projection")==0){
    inputArg$projection=F
  }
  
  {
    tmp=qread(paste0("~/myBucket/integrative_Tushar/SC_data/meta_analysis/refined_anno/",.ArgList$cellClass,"_Final_anno.qs"))
    pd_umap=tmp$pd
  }
  
  
  pd=list()
  for(i in 1:length(inputData$data)){
    tmp=inputData$data[[i]]
    if(sum(colnames(tmp) %in% pd_umap$sample)>0){
      tmp=tmp[,colnames(tmp) %in% pd_umap$sample]
      tmp_pd=pd_umap[match(colnames(tmp),pd_umap$sample),]
      colData(tmp)=colData(tmp)[,which(!colnames(colData(tmp)) %in% colnames(tmp_pd))]
      colData(tmp)=cbind(colData(tmp),tmp_pd)
      inputData$data[[i]]=tmp
      pd=c(pd,list(as.data.frame(colData(tmp))))
    } else {
      stop("Here")
      inputData$data[[i]]=NULL
    }
  }
  
  getfun<-function(x) {
    if(length(grep("::", x))>0) {
      parts<-strsplit(x, "::")[[1]]
      getExportedValue(parts[1], parts[2])
    } else {
      x
    }
  }
  
  pd2=do.call(getfun("plyr::rbind.fill"), pd)
  
  return(list(data=inputData,pd=pd2))
}

.myDotPlotArrangerFn=function(inputGene,inputPd,inputExpData,cluster_col="anno_cluster_res"){
  if(length(inputGene)>1){
    stop("Only one gene can be visualized at a time!")
  }
  tmp2=inputPd
  tmp2$cluster_name=tmp2[,cluster_col]
  tmp2=tmp2[!is.na(tmp2$cluster_name),]
  tmp2$cluster_name=gsub("inferred_","",as.character(tmp2$cluster_name))
  
  row.names(tmp2)=tmp2$sample
  #inputDataList=inputExpData;single_gene=inputGene;inputPd=tmp2;cluster_col="cluster_name";normalizeExp=F;gene_name_col="gene_short_name";dataset_name_col="ds_batch";return_plot=F;ignore_missing_anno=T;return_plot=F;ncores=10
  tmp=.my2d_dotPlot(inputDataList=inputExpData,single_gene=inputGene,inputPd=tmp2,cluster_col="cluster_name",normalizeExp=F,gene_name_col="gene_short_name",dataset_name_col="ds_batch",ignore_missing_anno=T,return_plot=F,ncores=10)#+ggtitle(igene)
  
  inputPd=inputPd[which(inputPd$ds_batch=="human_Tushar_iNPH"),]
  tmp2=.my2d_dotPlot(inputDataList=inputExpData,single_gene=inputGene,inputPd=inputPd,cluster_col="anno_cellState",normalizeExp=F,gene_name_col="gene_short_name",dataset_name_col="ds_batch",ignore_missing_anno=T,return_plot=F,ncores=10)
  
  if(!is.null(tmp2$p_data)){
    if(nrow(tmp2$p_data)>0){
      tmp2$p_data$dataset=paste0("org_",tmp2$p_data$dataset)
      p_data=rbind(tmp$p_data,tmp2$p_data)
      p=tmp2$plot
      
      p_data$features.plot=inputGene
      p_data$avg.exp.scaled[which(is.nan(p_data$avg.exp.scaled))]=0
      tmp_avg.exp=aggregate(avg.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      tmp_pct.exp=aggregate(pct.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      tmp_avg.exp.scaled=aggregate(avg.exp.scaled~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      
      tmp_data=merge(tmp_avg.exp,tmp_avg.exp.scaled,by=c("id","features.plot","dataset"),all=T)
      tmp_data=merge(tmp_data,tmp_pct.exp,by=c("id","features.plot","dataset"),all=T)
      
      
      tmp_feature.plot=tmp_data$features.plot
      tmp_data$features.plot=tmp_data$dataset
      tmp_data$gene=tmp_feature.plot
      p$data=tmp_data
      p=p+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5))+xlab("Dataset")+ggtitle(inputGene)
      
      
    }
  } else {
    if(!is.null(tmp$p_data)){
      p_data=tmp$p_data
      p=tmp$plot
      
      p_data$features.plot=inputGene
      p_data$avg.exp.scaled[which(is.nan(p_data$avg.exp.scaled))]=0
      tmp_avg.exp=aggregate(avg.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      tmp_pct.exp=aggregate(pct.exp~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      tmp_avg.exp.scaled=aggregate(avg.exp.scaled~id+features.plot+dataset,data=p_data,function(x) mean(x,na.rm = T))
      
      tmp_data=merge(tmp_avg.exp,tmp_avg.exp.scaled,by=c("id","features.plot","dataset"),all=T)
      tmp_data=merge(tmp_data,tmp_pct.exp,by=c("id","features.plot","dataset"),all=T)
      
      
      tmp_feature.plot=tmp_data$features.plot
      tmp_data$features.plot=tmp_data$dataset
      tmp_data$gene=tmp_feature.plot
      p$data=tmp_data
      p=p+theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust=0.5))+xlab("Dataset")+ggtitle(inputGene)
    }
    
  }
  return(p)
}

.my2dPlotArrangerFn=function(inputPd,inputGene,inputExpData,argList,ncores,FeaturePlotOnly=F){
  
  
  pd2=inputPd#[match(colnames(tmp),row.names(pd)),]
  
  ds_list=unlist(lapply(inputExpData,function(x) x$ds_batch[1]))
  
  res=list()
  for(ids in unique(ds_list)){
    tmp=inputExpData[which(ds_list==ids)]
    if(sum(inputPd$ds_batch==ids)>0){
      if(length(tmp)==1|FeaturePlotOnly){
        tmp=.mycBindFn(tmp)
        tmp=tmp[,tmp$sample %in% (inputPd$sample)]
        pd2=inputPd[match(tmp$sample,inputPd$sample),]
        
        if(class(tmp)=="SingleCellExperiment"){
          tmp=.extraExport2SeuratFn(tmp)
        }
        row.names(pd2)=colnames(tmp)
        p=.myFeaturePlot(inputSeurat = tmp,inputDimData = pd2,inputGenes = inputGene,combine_figs = F,order=F)
      } else {
        pd2=inputPd[inputPd$ds_batch==ids,]
        row.names(pd2)=pd2$sample
        p=.my2dPlot_counts(inputPCA=pd2,batch_values="batch_merging",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=list(inputGene),geneNameCol="gene_short_name",expData=tmp,ncolumns=6,ncores=ncores,combine_figs = F)
      }
      
      res=c(res,list(p))
      names(res)[length(res)]=ids
    }
    
  }
  
  geneNameList=unlist(lapply(res,names))
  geneNameList=unique(geneNameList)
  p_combined=NULL
  if(length(geneNameList)==1){
    for(geneName in geneNameList){
      
      for(ik in 1:length(res)){
        p=ggplot(data=NULL)
        if(sum(names(res[[ik]])==geneName)>0){
          p=res[[ik]][[geneName]]
        }
        p=p+ggtitle(names(res)[ik])
        if(ik==1){
          p_combined=p
        } else {
          p_combined=p_combined+p
        }
      }
      
      #ggsave(plot=p_combined,file=paste0("~/myBucket/torm.png"),width = 29,height = 29)
      
    }
  }
  
  
  
  return(p_combined)
}

.myDEloadFn=function(cellClass){
  seurat_marker=qread(paste0("~/myBucket/integrative_Tushar/",cellClass,"_arranged_markers.qs"))
  inph_marker=qread(paste0("~/myBucket/integrative_Tushar/",cellClass,"_iNPHanno_arranged_markers.qs"))
  all_markers=qread(paste0("~/myBucket/integrative_Tushar/",cellClass,"_all_genes.qs"))
  return(list(inph_marker=inph_marker,seurat_marker=seurat_marker,all_markers=all_markers))
}

