# Saved files -> variable name

## | file name prefix     | variable name                           |
## |----------------------+-----------------------------------------|
## | "UMAP_anno"          | pd                                      |
## | "pca_centroids"      | pca_centroid                            |
## | "harmony-embeddings" | harmony_embeddings                      |
## | "exp_merged"         | expData                                 |
## | "UMAP_centroid"      | UMAP_centroid                           |
## | "pca_anno"           | pca_res                                 |
## | "UMAP_anno"          | pd=as.matrix(pd[,c("UMAP_1","UMAP_2")]) |
## | "res_DE_wZscore"     | res_arranged                            |
## | "res_dataset_array"  | data                                    |
## | "res_meta"           | data                                    |
## |----------------------+-----------------------------------------|
## | global/varGenes      | tmp                                     |
## | pca/pca_Anno         | pca_res, pd                             |
## | pca_centroids        | pca_centroid                            |
## | kmeans_res_clustes   | res_clusters                            |
## | UMAP_res             | resUMAP                                 |
## | centroid_scores      | resDistances                            |
## | res_DE               | res_arranged                            |
## |----------------------+-----------------------------------------|

#set the main_source and utilities source code addresses
#change the harmony path in .sconline.embeddingFn(), if desired

library(googleCloudStorageR)
#gcs_source("vgazesta/code/mySC.R")
#gcs_source("vgazesta/code/my_SC_metaAnalysis.R")

source("mySC.R")
source("my_SC_metaAnalysis.R")

library(rliger)
library(qs)

#Accessory functions
.extra_sconline.scatterPlot_summary2d=function(object,reductionCols,n=300){
  
  if(sum(colnames(object) %in% reductionCols)!=2){
    stop("Specified reduction columns couldn't be found ")
  }
  
  
  yval=object[,reductionCols]
  xval=yval[,1]
  yval=yval[,2]
  
  rescoord=expand.grid(seq(min(xval),max(xval),length.out = n),seq(min(yval),max(yval),length.out = n))
  colnames(rescoord)=reductionCols
  
  knet=RANN::nn2(data=rescoord,query = object[,reductionCols],k=1,eps=0)
  counts=as.data.frame(table(knet$nn.idx[,1]))
  counts$Var1=as.character(counts$Var1)
  rescoord$id=1:nrow(rescoord)
  rescoord=merge(rescoord,counts,by.x="id",by.y="Var1")
  rescoord=rescoord[,-which(colnames(rescoord)=="id")]
  
  return(rescoord)
}

#assay = NULL; rev.pca = FALSE; weight.by.var = TRUE;
#verbose = TRUE; ndims.print = 1:5; nfeatures.print = 30; 
#reduction.key = "PC_"; seed.use = 42; approx = TRUE
.mySeuratRunPCA=function (object,projection_data=NULL, assay = NULL, npcs = 50, rev.pca = FALSE, weight.by.var = TRUE, 
          verbose = TRUE, ndims.print = 1:5, nfeatures.print = 30, 
          reduction.key = "PC_", seed.use = 42, approx = TRUE, ...) {
  require(Seurat)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  projection_pca=NULL
  if (rev.pca) {
    npcs <- min(npcs, ncol(x = object) - 1)
    pca.results <- irlba(A = object, nv = npcs, ...)
    total.variance <- sum(Seurat:::RowVar(x = t(x = object)))
    sdev <- pca.results$d/sqrt(max(1, nrow(x = object) - 
                                     1))
    if (weight.by.var) {
      feature.loadings <- pca.results$u %*% diag(pca.results$d)
    }
    else {
      feature.loadings <- pca.results$u
    }
    cell.embeddings <- pca.results$v
  } else {
    total.variance <- sum(Seurat:::RowVar(x = object))
    if (approx) {
      npcs <- min(npcs, nrow(x = object) - 1)
      pca.results <- irlba(A = t(x = object), nv = npcs, ...)
      feature.loadings <- pca.results$v
      
      if(!is.null(projection_data)){
        projection_pca=t(projection_data) %*% feature.loadings
        projection_pca=sweep(projection_pca,2,pca.results$d,"/")
      }
      
      
      sdev <- pca.results$d/sqrt(max(1, ncol(object) -1))
      if (weight.by.var) {
        cell.embeddings <- pca.results$u %*% diag(pca.results$d)
        if(!is.null(projection_data)){
          projection_pca = projection_pca %*% diag(pca.results$d)
        }
        
      }
      else {
        cell.embeddings <- pca.results$u
      }
    } else {
      npcs <- min(npcs, nrow(x = object))
      pca.results <- prcomp(x = t(object), rank. = npcs, ...)
      feature.loadings <- pca.results$rotation
      sdev <- pca.results$sdev
      if (weight.by.var) {
        cell.embeddings <- pca.results$x
      }
      else {
        cell.embeddings <- pca.results$x/(pca.results$sdev[1:npcs] * 
                                            sqrt(x = ncol(x = object) - 1))
      }
    }
  }
  rownames(x = feature.loadings) <- rownames(x = object)
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:npcs)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- Seurat::CreateDimReducObject(embeddings = cell.embeddings, 
                                         loadings = feature.loadings, assay = assay, stdev = sdev, 
                                         key = reduction.key, misc = list(total.variance = total.variance))
  if (verbose) {
    msg <- capture.output(print(x = reduction.data, dims = ndims.print, 
                                nfeatures = nfeatures.print))
    message(paste(msg, collapse = "\n"))
  }
  
  return(list(reduction.data=reduction.data,projection_pca=projection_pca))
}

.myPCAfn=function(data, argList,projection_data=NULL,saveFiles=T,...){
  
  library(future)
  require(purrr)
  plan("multicore", workers = 5)
  plan()
  options(future.globals.maxSize = 1000 * 1024^4)
  
  UMI_cor_thr=argList$UMI_cor_thr
  
  myPrepDR=function (scaledData, features, verbose = TRUE,projection_state=F) {
    
    data.use <- scaledData
    if (nrow(x = data.use) == 0) {
      stop("Data has not been scaled. Please run ScaleData and retry")
    }
    features.keep <- unique(x = features[features %in% rownames(x = data.use)])
    if (length(x = features.keep) < length(x = features)) {
      features.exclude <- setdiff(x = features, y = features.keep)
      if (verbose) {
        warning(paste0("The following ", length(x = features.exclude),
                       " features requested have not been scaled (running reduction without them): ",
                       paste0(features.exclude, collapse = ", ")))
      }
    }
    features <- features.keep
    # TODO jonah parallize buyt make sure chunked
    features.var <- apply(X = data.use[features, ], MARGIN = 1,
                          FUN = var)
    if(!projection_state){
      features.keep <- features[features.var > 0]
    } else {
      features.keep <- features
    }
    
    if (length(x = features.keep) < length(x = features)) {
      features.exclude <- setdiff(x = features, y = features.keep)
      if (verbose) {
        warning(paste0("The following ", length(x = features.exclude),
                       " features requested have zero variance (running reduction without them): ",
                       paste0(features.exclude, collapse = ", ")))
      }
    }
    features <- features.keep
    features <- features[!is.na(x = features)]
    data.use <- data.use[features, ]
    return(data.use)
  }
  
  mySeuratFn2_org=function(inputData,varFeatures){
    plan("sequential")
    inputData@assays$RNA@var.features=varFeatures
    inputData = ScaleData(inputData, verbose = FALSE,features=varFeatures)
    return(inputData)
  }
  
  mySeuratFn2=function(inputData,varFeatures){
    plan("sequential")
    inputData@assays$RNA@var.features=varFeatures
    inputData = Seurat:::ScaleData.default(inputData@assays$RNA@data, verbose = FALSE,features=varFeatures)
    return(inputData)
  }
  if(saveFiles){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_anno",argList = argList,pseudoImportant = F));F}, error=function(e) {return(T)})
  } else {
    reRunCheck=T
  }
  
  pca_final_res="Done"
  
  if(reRunCheck|!saveFiles){
    tmpInd=c()
    if(sum(!is.null(argList$covariates))>0&!is.null(data$data_m)){
      for(i in argList$covariates[!is.null(argList$covariates)]){
        tmpInd=c(tmpInd,which(is.na(data$data_m@meta.data[,i])))
      }
      if(length(tmpInd)>0){
        data$data_m=data$data_m[,-unique(tmpInd)]
      }
    }
    
    if(argList$indScaling){
      dataList=data$data
      
      varFeatures=c()
      if(!is.null(argList$HVG_list)){
        for(iHVG in argList$HVG_list){
          varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))])
        }
      } else {
        varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(argList$HVG_count-1))])
      }
      
      varFeatures=unique(varFeatures)
      
      # browser()
      dataList=parallel::mclapply(dataList,mySeuratFn2,varFeatures=varFeatures,mc.cores = argList$ncores)
      if(!is.null(projection_data)){
        projectionDataList=parallel::mclapply(projection_data$data,mySeuratFn2,varFeatures=varFeatures,mc.cores = argList$ncores)
      }
      gc()
      .mycBindFillFn=function(mat1,mat2){
        
        mat1c=setdiff(row.names(mat2),row.names(mat1))
        mat2c=setdiff(row.names(mat1),row.names(mat2))
        if(length(mat1c)>0){
          mat1cc=matrix(0,nrow=length(mat1c),ncol=ncol(mat1))
          row.names(mat1cc)=mat1c
          mat1=rbind(mat1,mat1cc)
        }
        if(length(mat2c)>0){
          mat2cc=matrix(0,nrow=length(mat2c),ncol=ncol(mat2))
          row.names(mat2cc)=mat2c
          mat2=rbind(mat2,mat2cc)
        }
        mat2=mat2[match(row.names(mat1),row.names(mat2)),]
        mat=cbind(mat1,mat2)
        return(mat)
      }
      
      if(F){
        #Jonah's
        rowNamesTable = dataList %>%
          map(~row.names(.@assays$RNA@scale.data)) %>%
          unlist %>% as.character() %>% table
        
        if(max(rowNamesTable) != min(rowNamesTable)){
          stop("Weird, asymetric row names. Make sure that works with cbind below")
        }
        
        resScaled = dataList %>%
          map(~.@assays$RNA@scale.data) %>%
          do.call(cbind, .)
      }
      
      resScaled=list()
      all_genes=unique(unlist(lapply(dataList,function(x) row.names(x))))
      for(i in 1:length(dataList)){
        tmp=dataList[[i]]
        tmp=tmp[match(all_genes,row.names(tmp)),]
        tmp[is.na(tmp)]=0
        row.names(tmp)=all_genes
        resScaled=c(resScaled,list(tmp))
        #dataList[[i]]=tmp
      }
      resScaled=do.call("cbind",resScaled)
      
      
      if(!is.null(projection_data)){
        projectionScaled=list()
        for(i in 1:length(projectionDataList)){
          tmp=projectionDataList[[i]]
          tmp=tmp[match(all_genes,row.names(tmp)),]
          tmp[is.na(tmp)]=0
          row.names(tmp)=all_genes
          projectionScaled=c(projectionScaled,list(tmp))
          
        }
        projectionScaled=do.call("cbind",projectionScaled)
        row.names(projectionScaled)=all_genes
      }
      
      pd=lapply(1:length(dataList),function(i){
        tmp=data$data[[i]]@meta.data
        tmp$sample=colnames(dataList[[i]])
        tmp
      })
      pd=do.call(eval(parse(text='plyr::rbind.fill')), pd)
      
      if(!is.null(projection_data)){
        pd_projection=lapply(1:length(projection_data$data),function(i){
          tmp=projection_data$data[[i]]@meta.data
          tmp$sample=colnames(projection_data$data[[i]])
          tmp
        })
        pd_projection=do.call(eval(parse(text='plyr::rbind.fill')), pd_projection)
        pd_all_projection=pd_projection
      }
      
      pd_all=pd
      gc()
      if(is.null(argList$HVG_list)){
        argList$HVG_list=argList$HVG_count
      }
      
      if(!saveFiles){
        argList$HVG_list=argList$HVG_count
      }
      
      if(is.null(argList$input_highly_var_genes)){
        for(iHVG in argList$HVG_list){
          
          argList$HVG_count=iHVG
          tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]
          
          #pca_res=.extraPCAfn(object=resScaled[which(row.names(resScaled) %in% tmp_varFeatures),], npcs = 30,findElbowPoint=F,seed.use = 42,reduction.key = "PC_",weight.by.var=T, approx = TRUE)
          #pca_res=pca_res$embeddings
          
          if(F){
            #sklearn: too slow
            library(reticulate)
            pd <- import("pandas", delay_load = TRUE)
            np <- import("numpy", delay_load = TRUE)
            skl_lr <- import("sklearn.decomposition", delay_load = TRUE)
            
            
            
            ipc=skl_lr$IncrementalPCA(n_components=50L, batch_size=10L)
            for(i in 1:length(dataList)){
              ipc = ipc$partial_fit(t(dataList[[i]]))
            }
          }
          
          if(F){
            #rsvd
            tst=rsvd::rsvd(A=t(myPrepDR(scaledData=resScaled,
                                        features=tmp_varFeatures, verbose = TRUE)),k=50)
          }
          
          object=myPrepDR(scaledData=resScaled,
                          features=tmp_varFeatures, verbose = TRUE)
          
          if(!is.null(projection_data)){
            projection_data=myPrepDR(scaledData=projectionScaled[row.names(object),],
                                     features=tmp_varFeatures, verbose = TRUE,projection_state = T)
          }
          
          projection_pca=NULL
          if(F){
            tst=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled,
                                                            features=tmp_varFeatures, verbose = TRUE),
                                            npcs = max(50,argList$nPCs+10))
            tst=tst@cell.embeddings
          } else {
            pca_res=.mySeuratRunPCA(object=object, npcs = max(50,argList$nPCs+10),projection_data=projection_data)
            projection_pca=pca_res$projection_pca
            pca_res=pca_res$reduction.data
          }
          
          pca_res=pca_res@cell.embeddings
          
          if(!is.null(projection_pca)){
            pca_res=rbind(pca_res,projection_pca)
            pd=plyr::rbind.fill(pd_all,pd_all_projection)
          } else {
            pd=pd_all
          }
          
          pd=pd[pd$sample %in% row.names(pca_res),]
          row.names(pd)=pd$sample
          pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
          if(saveFiles){
            save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
          } else {
            pca_final_res=list(pd=pd,pca_res=pca_res)
          }
          
          gc()
        }
      } else {
        print("Highly variable genes are provided in argList, using it!")
        {
          
          tmp_varFeatures=as.character(argList$input_highly_var_genes)
          
          #pca_res=.extraPCAfn(object=resScaled[which(row.names(resScaled) %in% tmp_varFeatures),], npcs = 30,findElbowPoint=F,seed.use = 42,reduction.key = "PC_",weight.by.var=T, approx = TRUE)
          #pca_res=pca_res$embeddings
          pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled,
                                                          features=tmp_varFeatures, verbose = TRUE,
                                                          ...),
                                          npcs = max(50,argList$nPCs+10))
          pca_res=pca_res@cell.embeddings
          
          pd=pd_all[pd_all$sample %in% row.names(pca_res),]
          row.names(pd)=pd$sample
          pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
          if(saveFiles){
            save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
          } else {
            pca_final_res=list(pd=pd,pca_res=pca_res)
          }
          
          gc()
        }
      }
      
      
    }
    
    if(!argList$indScaling){
      
      varFeatures=c()
      if(is.null(argList$input_highly_var_genes)){
        if(!is.null(argList$HVG_list)){
          for(iHVG in argList$HVG_list){
            varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))])
          }
        } else {
          varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(argList$HVG_count-1))])
        }
        varFeatures=unique(varFeatures)
        
        data$data_m <- ScaleData(data$data_m, verbose = FALSE,features=varFeatures, vars.to.regress =argList$covariates)
        
        if(is.null(argList$HVG_list)){
          argList$HVG_list=HVG_count
        }
        {
          for(iHVG in argList$HVG_list){
            argList$HVG_count=iHVG
            tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]
            data$data_m = RunPCA(data$data_m,features=tmp_varFeatures,npcs = max(50,argList$nPCs+10))
            pca_res=data$data_m@reductions$pca@cell.embeddings
            pd=data$data_m@meta.data
            
            cordata=c()
            for(ik in 1:ncol(pca_res)){
              cordata=c(cordata,cor(pca_res[,ik],pd$QC_Gene_unique_count))
            }
            if(length(which(abs(cordata)>UMI_cor_thr))>0){
              pca_res=pca_res[,-which(abs(cordata)>UMI_cor_thr)]
            }
            
            if(saveFiles){
              save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
            }  else {
              pca_final_res=list(pd=pd,pca_res=pca_res)
            }
            
          }
        }
      } else {
        
        tmp_varFeatures=as.character(argList$input_highly_var_genes)
        data$data_m <- ScaleData(data$data_m, verbose = FALSE,features=tmp_varFeatures, vars.to.regress =argList$covariates)
        
        data$data_m = RunPCA(data$data_m,features=tmp_varFeatures,npcs = max(50,argList$nPCs+10))
        pca_res=data$data_m@reductions$pca@cell.embeddings
        pd=data$data_m@meta.data
        
        cordata=c()
        for(ik in 1:ncol(pca_res)){
          cordata=c(cordata,cor(pca_res[,ik],pd$QC_Gene_unique_count))
        }
        if(length(which(abs(cordata)>UMI_cor_thr))>0){
          pca_res=pca_res[,-which(abs(cordata)>UMI_cor_thr)]
        }
        
        if(saveFiles){
          save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        }  else {
          pca_final_res=list(pd=pd,pca_res=pca_res)
        }
        
      }
      
      
      
    }
  }
  
  
  return(pca_final_res)
}

.extra_sconline_LigerToExpSet=function(inputData,organism,limit=NULL){
  require(rliger)
  require(scater)
  require(scran)
  
  expList=list()
  for(i in 1:length(inputData@raw.data)){
    tmp=inputData@raw.data[[i]]
    if(!is.null(limit)){
      if(ncol(tmp)>limit){
        tmp=tmp[,sample(ncol(tmp),min(limit,ncol(tmp)))]
      }
    }
    exp=SingleCellExperiment(assays = list(counts = tmp),colData = data.frame(sampleName=colnames(tmp),stringsAsFactors = F),rowData=data.frame(gene=row.names(tmp),stringsAsFactors = F))
    expList=c(expList,list(exp))
    names(expList)[length(expList)]=names(inputData@raw.data)[i]
  }
  
  expList=.extra_sconline_cBindFn(expList,batchNames = names(expList))
  
  pd=as.data.frame(inputData@cell.data)
  
  if(length(inputData@clusters)>0){
    if(!all(row.names(pd)==names(inputData@clusters))){
      stop("Error!")
    }
    pd$clusters=as.character(inputData@clusters)
  }
  
  pd=pd[match(expList$sampleName,row.names(pd)),]
  if(!all(row.names(pd)==expList$sample)){
    stop("Error!")
  }
  
  
  if(sum(colnames(pd)=="batch_merging")>0){
    pd$batch_merging_prv=pd$batch_merging
  }
  
  pd$batch_merging=expList$batch_merging
  
  
  row.names(pd)=as.character(expList$sampleName)
  
  data=.myExpSetCreatorFn(inputExpData=counts(expList),organism=organism,minExpCells=0,inputPdata=pd,inputFdata=as.data.frame(rowData(expList)),addExtraAnno=F,server=T)
  return(data)
}

.extra_sconline_SplitObject=function(object,colName){
  if(class(object)=="Seurat"){
    res=Seurat::SplitObject(object,split.by = colName)
  } else{
    pd=colData(object)[,colName]
    res=list()
    for(i in unique(pd)){
      res=c(res,list(object[,which(pd==i)]))
      names(res)[length(res)]=i
    }
  }
  return(res)
}

.extra_sconline_cBindFn=function(inputList,batchNames=NULL,verbose=F){
  #cbinds multiple singleCellExpression datasets with differring number of rows.
  #inputList: the list of datasets to be merged
  #batchNames: the batch name to be assigned to each dataset. length(batchNames)==length(inputList)
  res_m=""
  if(!is.null(batchNames)){
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
    } else if(length(inputList)>2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
      for(i in 3:length(inputList)){
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",batchNames[i]))
      }
    }
  } else {
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
    } else if(length(inputList)>2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
      for(i in 3:length(inputList)){
        
        #x1=res_m;x2=inputList[[i]];batchNames=c("","")
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",""))
        if(verbose){
          print(paste("dataset:",i,"; nrow:",nrow(res_m),"; ncol",ncol(res_m)))
        }
      }
    }
  }
  
  print("batch information is in the batch_merging variable")
  return(res_m)
  
}

.extra_sconline.exp_creatorFn=function(argList,inputExpData,batch_variable,organism="Mouse",addAnno=F){
  
  if(class(inputExpData)!=class(list())){
    if(class(inputExpData)=="Seurat"){
      inputExpData=.myExpSetCreatorFn(inputExpData=inputExpData@assays$RNA@counts,
                                      organism=organism,
                                      minExpCells=0,
                                      inputPdata=as.data.frame(inputExpData@meta.data),
                                      inputFdata=as.data.frame(inputExpData@assays$RNA@meta.features),
                                      addExtraAnno=addAnno,server=T)
    } else if(class(inputExpData)=="liger"){
      inputExpData=.extra_sconline_LigerToExpSet(inputData=data,organism=organism,limit=NULL)
    } else if(class(inputExpData)!="SingleCellExperiment"){
      stop("Unrecognized input expression data format!")
    }
    
    if(sum(duplicated(row.names(inputExpData)))>0){
      print(paste(sum(duplicated(row.names(inputExpData))),"Duplicate gene ids were found! duplicates were randomly removed from the data"))
      inputExpData=inputExpData[!duplicated(row.names(inputExpData)),]
    }
    
    inputExpData$batch_merging=as.character(colData(inputExpData)[,batch_variable])
    if(class(inputExpData)!=class(list())){
      inputExpData=.myReadData_spliterFn(inputData=inputExpData,removeHighExp=argList$excludeHighExp)
    }
  } else {
    if(length(inputExpData)>2&sum(names(inputExpData)=="data_m")==0){
      inputExpData=list(data=inputExpData,data_m=NULL)
    }
  }
  
  
  
  return(inputExpData)
}

.extra_sconline.AffinityFn=function(sim_mat){
  
  thr=apply(sim_mat,1,function(x){
    x=x[!is.na(x)]
    y=NA
    if(length(x)>2){
      x=x[order(x,decreasing = T)]
      y=x[3]
    }
    y
  })
  thr=median(thr,na.rm=T)
  affinity=exp(-3*(pmax(thr-sim_mat,0)))
  #affinity=res_cells/thr
  #thr=apply(affinity,1,function(x) {
  #  x=x[order(x,decreasing = T)]
  #  x[min(20,length(x))]
  #})
  #thr=matrix(thr,nrow=nrow(affinity),ncol=ncol(affinity),byrow = F)
  #affinity[affinity<thr]=0
  #affinity=affinity+t(affinity)
  
  #affinity=t(apply(affinity,1,function(x) x/sum(x)))
  
  return(affinity)
  
}

.extra_sconline.NetVisFn=function(net,input_pd,input_umap_centroid){
  net$weight=net$score
  #net=net[net$score>0.8,]
  net$Fnode=net$source
  net$Snode=net$target
  res_net=net
  net=res_net[order(res_net$weight,decreasing = T),]
  #net=net[!duplicated(net$Snode),]
  net$Fnode=gsub("C","",net$Fnode)
  net$Snode=gsub("C","",net$Snode)
  net=net[!duplicated(paste0(net$Fnode,"_",net$Snode)),]
  net = network::network(net[,c("Fnode","Snode")], directed = T,matrix.type="edgelist")
  
  netVerNames=network::network.vertex.names(net)
  network::set.edge.attribute(net, "weight", ((res_net$weight-min(res_net$weight)+0.05)/(1-min(res_net$weight)+0.05))^3)
  
  centroid_layout=input_umap_centroid[match(gsub("C","",as.character(netVerNames)),as.character(input_umap_centroid$centroid)),-1]
  centroid_layout=as.matrix(centroid_layout)
  colnames(centroid_layout)=c("x","y")
  
  net=ggnetwork:::fortify.network(net,layout = centroid_layout)
  
  pd_summary=input_pd
  
  
  scale_factor=net[!duplicated(net$vertex.names),]
  scale_factor=merge(scale_factor,input_umap_centroid,by.x="vertex.names",by.y="centroid")
  scale_factor1=lm(x~UMAP_1,data=scale_factor)
  pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
  scale_factor2=lm(y~UMAP_2,data=scale_factor)
  pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
  
  pd_summary$xend=pd_summary$UMAP_1
  pd_summary$yend=pd_summary$UMAP_2
  pd_summary$vertex.names=""
  pd_summary$color="gray"
  
  #predicting the background color
  library(ggnetwork)
  #centroids_ideal=c("173","192","191","187","127","194")
  p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray")+
    geom_edges( color = "black",aes(size=weight),arrow = arrow(length = unit(6, "pt"), type = "closed")) +
    geom_label(aes(label=vertex.names))+
    theme_blank()+scale_size_continuous(range = c(0.06,1))+scale_color_identity()+scale_fill_identity()+theme(legend.position = "none")
  return(p)
}

varibow <- function(n_colors) {
  sats <- rep_len(c(0.55,0.7,0.85,1),length.out = n_colors)
  vals <- rep_len(c(1,0.8,0.6),length.out = n_colors)
  sub("FF$","",grDevices::rainbow(n_colors, s = sats, v = vals))
}


V<-View

library(stringr)
.mcsaveRDS <- function(object,file,mc.cores=min(parallel::detectCores(),10, na.rm=T)) {
  file = str_replace_all(file, " ", "\\\\ ")
  con <- pipe(paste0("pigz -p",mc.cores," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}
.mcreadRDS <- function(file,mc.cores=min(parallel::detectCores(),10, na.rm=T)) {
  file = str_replace_all(file, " ", "\\\\ ")
  con <- pipe(paste0("pigz -d -c -p",mc.cores," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}
.printDot <- function(){
  cat(".")
  # TODO check if exists first, for jupyter/console
  flush.console()
}
varSizes <- function(){
  sizes <- sapply(ls(envir=.GlobalEnv), function(n) object.size(get(n)), simplify = FALSE);
  print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))
}

.extra_sconline.NetIdFn=function(inputNet,col1="source1",col2="source2"){
  id=paste0(inputNet[,col1],"_",inputNet[,col2])
  id[inputNet[,col2]>inputNet[,col1]]=paste0(inputNet[,col2],"_",inputNet[,col1])[inputNet[,col2]>inputNet[,col1]]
  return(id)
}


.extra_sconline.purityAnalysisFn=function(argList,inputEmbeddings=NULL,inputPhenoData=NULL,expData=NULL,run_harmony=F,organism,batch_variable="batch_merging",addAnno=F,collapse_datasets=T,minCellCountThr=4,analysis_seed=1,extendedMode=F,umap.method='umap-learn',L2Norm=T,mergePseudocells=T,hierarchical_refinement=T,colNormalize=F,merging_strength=0.15){
  #collapse_datasets=F;minCellCountThr=4;analysis_seed=1;extendedMode=F;umap.method='umap-learn'
  require(purrr)
  require(furrr)
  
  options(future.globals.maxSize= 750*1024^4)
  
  res_embeddings=.sconline.embeddingFn(argList,inputEmbeddings=inputEmbeddings,run_harmony=run_harmony,pd=inputPhenoData,inputBatchCol=batch_variable)
  
  res_umap=.sconline.umapFn(argList,umap.method=umap.method)
  
  
  set.seed(analysis_seed)
  supportingFractionThr=argList$DE_supportingFractionThr
  n.adaptiveKernel=argList$DE_n.adaptiveKernel
  nPropIter=argList$DE_nPropIter
  
  #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("pca_centroids",argList=argList))
  
  
  harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),]
  if(sum(is.na(harmony_embeddings))>0){
    stop("Error in matching Names!")
  }
  
  
  
  if(is.null(expData)){
    tmp=qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
    expData=SplitObject(tmp, split.by = "batch_merging")
  } else {
    data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData = expData,batch_variable = batch_variable,organism = organism,addAnno=addAnno)
    expData = data$data
    
    tmp_pd=NULL
    for(i in 1:length(expData)){
      tmp_pd=rbind(tmp_pd,data.frame(sample=colnames(expData[[i]]),batch_merging=expData[[i]]$batch_merging,stringsAsFactors = F))
    }
    #pd=pd[row.names(pd) %in% tmp_pd$sample,]
    tmp_pd=tmp_pd[match(row.names(pd),tmp_pd$sample),]
    if(sum(is.na(tmp_pd$batch_merging))>0){
      stop("Error! phenoData doesn't match with the expression data")
    }
    pd$batch_merging=tmp_pd$batch_merging
  }
  
  pcaList=split(as.data.frame(harmony_embeddings),pd$batch_merging)
  
  
  dataArranged = parallel::mclapply(expData, function(thisExp){
    thisExp=thisExp[,colnames(thisExp) %in% row.names(pd)]
    
    return(list(
      dsName=as.character(thisExp$batch_merging[1]),
      pcaData=pcaList[[as.character(thisExp$batch_merging[1])]]))
  })
  
  if(sum(unlist(lapply(dataArranged,function(x) nrow(x$pcaData))))<sum(unlist(lapply(expData,ncol)))){
    warning("Annotation was found for only a subset of exp data, that subset was only used in the analysis!")
  }
  
  
  #if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))){
  cat("           Constructing the propagation matrices ...\n")
  dataArranged_old = dataArranged
  dataArranged = dataArranged %>% keep(~dim(.$pcaData)[[1]] >= argList$min_ds_size)
  print(paste0("Going from data size ", length(dataArranged_old), " to ",
               length(dataArranged)))
  rm(dataArranged_old)
  
  if(argList$do.split.prop){
    res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final_v2_split_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
  } else {
    #centroidPCAdata=pca_centroid;exCentroids=NULL;runIndx=1;n.neighbors=argList$prop.n.neighbors;batchPCAdata=harmony_embeddings;n.trees=50;NNmethod="annoy";L2Norm=T;mergePseudocells=T;hierarchical_refinement=T;batch_variable="batch_merging"
    res=.myConcensusDEFn_step2_detail_newprop3_final_v14(dataArranged=dataArranged,centroidPCAdata=pca_centroid,argList=argList,exCentroids=NULL,runIndx=1,batchPCAdata=harmony_embeddings,n.neighbors=argList$prop.n.neighbors,L2Norm=L2Norm,mergePseudocells=T,batch_variable=batch_variable,colNormalize=F,hierarchical_refinement=T,merging_strength=merging_strength)
  }
  
  cat("           Performing purity analysis ...\n")
  #res_prop=res;annoCol="anno_cellState";return_plot=T;min_effective_size=5
  p=.sconline.anno2pseudocell_tmp(res_prop=res$dataArranged,argList=argList,annoCol="anno_cellState",collapse_datasets=collapse_datasets,return_plot=T,min_effective_size=5)
  p$prop_mat=res$prop_mat
  return(p)
}

.myConcensusDEFn_step2=function(argList,expData=NULL,minCellCountThr=4,meta_method="Stouffer",addClusteringModule=F,analysis_seed=1,extendedMode=F,L2Norm=T,mergePseudocells=T,merging_strength=0.3,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2){
  
  #argList=.ArgList;inputEmbeddings=NULL;inputPhenoData=NULL;inputExpData=NULL;organism="Mouse";batch_variable="batch_merging";run_harmony=F;addAnno=F;addClusteringModule=F;L2Norm=T;mergePseudocells=T;generateUMAP=F;extendedMode=F;merging_strength=0.3
  
  #minCellCountThr=4;meta_method="Stouffer";addClusteringModule=F;analysis_seed=1;extendedMode=F;L2Norm=T;mergePseudocells=T;merging_strength=0.3;analysis_seed=1;analysis_seed=1;include.singletons=T;colNormalize=T
  #expData=NULL
  
  
  myMetaAnalysisFn=function(res,argList,prop_mat,secondRun=F,extendedMode=F){
    
    my_pb <- function(total){
      progress::progress_bar$new(
        format = " [:bar] :current/:total (:percent) eta: :eta (already :elapsed)",
        total = total, clear = FALSE, width= 80)
    }
    
    myStoufferFn_overheadissue=function(zscore_array,weight_array,argList){
      weight_array=sqrt(weight_array)
      z_w=apply(zscore_array*weight_array, c(2,3), sum)
      
      w2=matrix(0,nrow=nrow(z_w),ncol=ncol(z_w))
      
      myStoufferWeightFn=function(i,zscore_array,weight_array){
        
        res_stouffer=rep(1,dim(zscore_array)[2])
        tmp_zscore=t(zscore_array[,i,])
        
        empty_ind=which(colSums(tmp_zscore)==0)
        if(length(empty_ind)<(ncol(tmp_zscore))){
          tmp_zscore[tmp_zscore==0]=NA
          tmp_zscore_cor=cor(tmp_zscore,method = "spearman",use="pairwise.complete.obs")
          tmp_zscore_cor[is.na(tmp_zscore_cor)]=0
          tmp_zscore_cor[which(tmp_zscore_cor<=0)]=0
          tmp_w2=unlist(lapply(1:dim(weight_array)[3],function(j) {
            sum(rowSums(Matrix::crossprod(t(weight_array[,i,j]),weight_array[,i,j])*tmp_zscore_cor))
          }))
          res_stouffer=tmp_w2
          
        }
        return(res_stouffer)
      }
      
      #w2=parallel::mclapply(1:dim(zscore_array)[2],myStoufferWeightFn,zscore_array=zscore_array,weight_array=weight_array,mc.cores = argList$ncores)
      w2=lapply(1:dim(zscore_array)[2],myStoufferWeightFn,zscore_array=zscore_array,weight_array=weight_array)
      w2=do.call("rbind",w2)
      w2=sqrt(w2)
      if(sum(w2>0&w2<0.99)>0){
        #warning("Possible issue in the weights")
        w2[which(w2>0&w2<1)]=1
      }
      z_w=z_w/w2
      z_w[which(w2==0)]=0
      return(z_w)
    }
    
    myStoufferFn_org=function(zscore_array,weight_array,argList){
      weight_array=sqrt(weight_array)
      z_w=apply(zscore_array*weight_array, c(2,3), sum)
      
      w2=matrix(0,nrow=nrow(z_w),ncol=ncol(z_w))
      
      myStoufferWeightFn=function(i,zscore_array,weight_array){
        
        res_stouffer=rep(1,dim(zscore_array)[2])
        tmp_zscore=t(zscore_array[,i,])
        
        empty_ind=which(colSums(tmp_zscore)==0)
        if(length(empty_ind)<(ncol(tmp_zscore))){
          tmp_zscore[tmp_zscore==0]=NA
          tmp_zscore_cor=cor(tmp_zscore,method = "spearman",use="pairwise.complete.obs")
          tmp_zscore_cor[is.na(tmp_zscore_cor)]=0
          tmp_zscore_cor[which(tmp_zscore_cor<=0)]=0
          tmp_w2=unlist(lapply(1:dim(weight_array)[3],function(j) {
            sum(rowSums(Matrix::crossprod(t(weight_array[,i,j]),weight_array[,i,j])*tmp_zscore_cor))
          }))
          res_stouffer=tmp_w2
          
        }
        return(res_stouffer)
      }
      
      #w2t=parallel::mclapply(1:dim(zscore_array)[2],myStoufferWeightFn,zscore_array=zscore_array,weight_array=weight_array,mc.cores = argList$ncores)
      w2=list()
      windowsize=100
      if(dim(zscore_array)[2]-max(seq(1,dim(zscore_array)[2],100))<2){
        windowsize=101
      }
      for(i in seq(1,dim(zscore_array)[2],windowsize)){
        tmp=parallel::mclapply(1:length(i:min(i+windowsize-1,dim(zscore_array)[2])),myStoufferWeightFn,zscore_array=zscore_array[,i:min(i+windowsize-1,dim(zscore_array)[2]),],weight_array=weight_array[,i:min(i+windowsize-1,dim(zscore_array)[2]),],mc.cores = argList$ncores)
        w2=c(w2,tmp)
      }
      w2=do.call("rbind",w2)
      w2=sqrt(w2)
      if(sum(w2>0&w2<0.99)>0){
        #warning("Possible issue in the weights")
        w2[which(w2>0&w2<1)]=1
      }
      z_w=z_w/w2
      z_w[which(w2==0)]=0
      return(z_w)
    }
    
    myStoufferFn=function(zscore_array,weight_array,argList){
      
      
      for(i in 1:length(weight_array)){
        weight_array[[i]]=sqrt(weight_array[[i]])
      }
      
      z_w=zscore_array[[1]]*weight_array[[1]]
      for(i in 2:length(zscore_array)){
        z_w=z_w+zscore_array[[i]]*weight_array[[i]]
      }
      
      
      
      myStoufferWeightFn=function(i,zscore_array,weight_array){
        
        res_stouffer=rep(1,dim(zscore_array)[2])
        tmp_zscore=t(zscore_array[,i,])
        
        empty_ind=which(colSums(tmp_zscore)==0)
        if(length(empty_ind)<(ncol(tmp_zscore))){
          tmp_zscore[tmp_zscore==0]=NA
          tmp_zscore_cor=cor(tmp_zscore,method = "spearman",use="pairwise.complete.obs")
          tmp_zscore_cor[is.na(tmp_zscore_cor)]=0
          tmp_zscore_cor[which(tmp_zscore_cor<=0)]=0
          tmp_w2=unlist(lapply(1:dim(weight_array)[3],function(j) {
            sum(rowSums(Matrix::crossprod(t(weight_array[,i,j]),weight_array[,i,j])*tmp_zscore_cor))
          }))
          res_stouffer=tmp_w2
          
        }
        return(res_stouffer)
      }
      
      zscore_array2=array(0,dim = c(length(zscore_array),nrow(zscore_array[[1]]),ncol(zscore_array[[2]])))
      weight_array2=array(0,dim = c(length(weight_array),nrow(weight_array[[1]]),ncol(weight_array[[2]])))
      
      for(itr in 1:length(zscore_array)){
        zscore_array2[itr,,]=as.matrix(zscore_array[[itr]])
      }
      for(itr in 1:length(weight_array)){
        weight_array2[itr,,]=as.matrix(weight_array[[itr]])
      }
      
      rm(zscore_array,weight_array)
      gc()
      
      #w2t=parallel::mclapply(1:dim(zscore_array)[2],myStoufferWeightFn,zscore_array=zscore_array,weight_array=weight_array,mc.cores = argList$ncores)
      w2=list()
      windowsize=100
      if(dim(zscore_array2)[2]-max(seq(1,dim(zscore_array2)[2],100))<2){
        windowsize=101
      }
      for(i in seq(1,dim(zscore_array2)[2],windowsize)){
        
        tmp=parallel::mclapply(1:length(i:min(i+windowsize-1,dim(zscore_array2)[2])),myStoufferWeightFn,zscore_array=zscore_array2[,i:min(i+windowsize-1,dim(zscore_array2)[2]),],weight_array=weight_array2[,i:min(i+windowsize-1,dim(zscore_array2)[2]),],mc.cores = argList$ncores)
        w2=c(w2,tmp)
      }
      w2=do.call("rbind",w2)
      w2=sqrt(w2)
      if(sum(w2>0&w2<0.99)>0){
        #warning("Possible issue in the weights")
        w2[which(w2>0&w2<1)]=1
      }
      z_w=z_w/w2
      z_w[which(w2==0)]=0
      return(z_w)
    }
    
    
    myWeightedMeanFn=function(data_array,weight_array){
      z_w=apply(data_array*weight_array, c(2,3), sum)
      w=apply(weight_array, c(2,3), sum)
      if(sum(w>0&w<0.99)>0){
        #warning("Possible issue in the weights")
        w[which(w>0&w<1)]=1
      }
      z_w=z_w/w
      z_w[which(w==0)]=0
      return(z_w)
    }
    
    myPctAvgfn_org=function(pct_array,centroid_weight_array,effective_size_array){
      myWeightedMeanFn=function(data_array,weight_array){
        z_w=apply(data_array*weight_array, c(2,3), function(x) sum(x,na.rm = T))
        w=apply(weight_array, c(2,3), function(x) sum(x,na.rm = T))
        if(sum(w>0&w<0.99)>0){
          #warning("Possible issue in the weights")
          w[which(w>0&w<1)]=1
        }
        z_w=z_w/w
        z_w[which(w==0)]=0
        return(z_w)
      }
      
      effective_size_array=centroid_weight_array*effective_size_array
      #pct_array[centroid_weight_array<0.9]=NA
      med_pct.1=myWeightedMeanFn(data_array=pct_array,weight_array=effective_size_array)
      row.names(med_pct.1)=dimnames(pct_array)[[2]]
      return(med_pct.1)
      
    }
    
    myPctAvgfn=function(pct_array,effective_size_array){
      
      all(names(pct_array)==names(effective_size_array))
      
      z_w=pct_array[[1]]*effective_size_array[[1]]
      for(itr in 2:length(pct_array)){
        z_w=z_w+pct_array[[itr]]*effective_size_array[[itr]]
      }
      
      w=effective_size_array[[1]]
      for(witr in 2:length(effective_size_array)){
        w=w+effective_size_array[[witr]]
      }
      if(sum(w@x>0&w@x<0.99)>0){
        #warning("Possible issue in the weights")
        w@x[which(w@x>0&w@x<1)]=1
      }
      w_ident=w
      w_ident@x=rep(1,length(w_ident@x))
      z_w=Matrix::drop0(z_w*w_ident)
      z_w@x=z_w@x/w@x
      
      return(z_w)
      
    }
    
    
    myHarmonicPvalFn=function(zscore_array,weight_array,offset=10^(-25)){
      pval=pnorm(zscore_array,lower.tail = F)
      weight_normFactor=apply(weight_array,c(2,3),sum)
      weight_array_n=sweep(weight_array,c(2,3),weight_normFactor,"/")
      
      pval[which(pval<0.1)]=pval[which(pval<0.1)]+offset
      pval[which(pval>0.9)]=pval[which(pval>0.9)]-offset
      resHarmony=apply(weight_array_n/pval,c(2,3),sum)
      resHarmony=1/resHarmony
      resHarmony[is.na(resHarmony)]=0.5
      resHarmony=pmin(resHarmony,1-offset)
      resHarmony=qnorm(resHarmony,lower.tail = F)
      
      return(resHarmony)
    }
    
    library(magrittr)
    library(purrr)
    library(stringr)
    library(future)
    library(furrr)
    
    if(F){
      #Jonah's
      exp_final_fp = file.path(argList$saveDir, "detail_exp_final")
      all_exp_files = list.files(exp_final_fp)
      
      all_exp_files %<>% set_names(., str_remove(., "[.]RDS"))
      
      my_pb <- function(total){
        progress::progress_bar$new(
          format = " [:bar] :current/:total (:percent) eta: :eta (already :elapsed)",
          total = total, clear = FALSE, width= 80)
      }
      this_pb = my_pb(length(all_exp_files))
      res_exp_final = map(all_exp_files, ~{
        this_pb$tick()
        # don't do mc when testing bc doesn't like swap space
        return(.mcreadRDS(file.path(exp_final_fp, .)))
      })
      this_pb$terminate()
      res = res_exp_final
      rm(res_exp_final)
    }
    
    gc()
    plan(sequential)
    
    
    cat("           Arranging z-score results ...\n")
    geneList=list()
    for(gitr in 1:length(res)){
      geneList=c(geneList,list(data.frame(gene=colnames(res[[gitr]]$matWeights),weight=colMaxs(res[[gitr]]$matWeights),stringsAsFactors = F)))
    }
    
    geneList=as.data.frame(data.table::rbindlist(geneList))
    geneList=aggregate(weight~gene,data=geneList,function(x) sum(x>0.9))
    geneList=geneList[which(geneList$weight>(length(res)*argList$DE_supportingFractionThr)),]
    geneList=as.character(geneList$gene)
    
    divby = ceiling(length(res)/50)
    this_pb = my_pb(ceiling(length(res)/divby))
    for(gitr in 1:length(res)){
      if(gitr %% divby == 0){
        this_pb$tick()
      }
      
      set_name_list=setdiff(names(res[[gitr]]),c("dsName","prop_mat","pseudocell_sim_mat","gene_name_list","pseudocell_name_list"))
      if(!secondRun){
        set_name_list=c("pct.1","zscore","matWeights","matEffectiveSize")
      }
      
      for(dsitr in set_name_list){
        tmp=res[[gitr]][[dsitr]]
        tmp_c=setdiff(geneList,colnames(tmp))
        if(length(tmp_c)>0){
          tmp_c=as(matrix(0,nrow=nrow(tmp),ncol=length(tmp_c),dimnames = list(row.names(tmp),tmp_c)),"dgCMatrix")
          tmp=cbind(tmp,tmp_c)
        }
        tmp=tmp[,match(geneList,colnames(tmp))]
        colnames(tmp)=geneList
        
        res[[gitr]][[dsitr]]=as(tmp,"dgCMatrix")
      }
      
      if(gitr %% 10 == 0){
        gc()
      }
    }
    this_pb$terminate()
    
    if(secondRun){
      matWeights=data.frame(centroid=row.names(prop_mat),weight=apply(res[[1]]$matWeights,1,max),effectiveSize=apply(res[[1]]$matEffectiveSize,1,max),dataset=res[[1]]$dsName,stringsAsFactors = F)
      
      if(length(res)>1){
        for(idsScore in 2:length(res)){
          tmp=data.frame(centroid=row.names(prop_mat),weight=apply(res[[idsScore]]$matWeights,1,max),effectiveSize=apply(res[[idsScore]]$matEffectiveSize,1,max),dataset=res[[idsScore]]$dsName,stringsAsFactors = F)
          matWeights=rbind(matWeights,tmp)
        }
      }
      
      matEffectiveSize=reshape2::dcast(centroid~dataset,data = matWeights,value.var = "effectiveSize")
      matWeights=reshape2::dcast(centroid~dataset,data = matWeights,value.var = "weight")
      resDistances=list(matWeights=matWeights,matEffectiveSize=matEffectiveSize)
      save(resDistances,file=do.call('.myFilePathMakerFn',args=c(list(filename="centroid_scores",uniformImportant=T,propImportant=T,qsFormat=T),argList)))
      rm(matEffectiveSize,matWeights,resDistances)
    }
    
    if(extendedMode){
      z_q.8=apply(res_array$zscore,c(2,3),function(x) {
        if(sum(x>3)>0){
          sl_thr=quantile(x[x!=0],0.8)
          sl_ind=which(x>sl_thr)
          x[sl_ind]=sl_thr
        } 
        if(sum(x<(-3))>0){
          sl_thr=quantile(x[x!=0],0.2)
          sl_ind=which(x<sl_thr)
          x[sl_ind]=sl_thr
        }
        x
      })
    }
    
    
    #rearranging the res object
    res_rearranged=NULL
    for(i in setdiff(names(res[[1]]),"dsName")){
      tmp=list()
      for(j in 1:length(res)){
        tmp=c(tmp,list(res[[j]][[i]]))
        names(tmp)[length(tmp)]=res[[j]]$dsName
      }
      res_rearranged=c(res_rearranged,list(tmp))
      names(res_rearranged)[length(res_rearranged)]=i
    }
    
    if(secondRun){
      for(i in names(res_rearranged)){
        tmp=res_rearranged[[i]]
        qsave(tmp,file=.myFilePathMakerFn(paste0("res_array_",i),argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
      }
    }
    
    res_rearranged=res_rearranged[c("zscore","pct.1","matWeights","matEffectiveSize")]
    rm(res)
    
    res_eff_size=list()
    for(i in 1:length(res_rearranged[["matWeights"]])){
      res_eff_size=c(res_eff_size,list(as(res_rearranged$matWeights[[i]]*res_rearranged$matEffectiveSize[[i]],"dgCMatrix")))
      names(res_eff_size)[length(res_eff_size)]=names(res_rearranged$matEffectiveSize)[i]
    }
    
    res_rearranged$matEffectiveSize=res_eff_size
    rm(res_eff_size)
    res_rearranged=res_rearranged[c("zscore","pct.1","matEffectiveSize")]
    
    gc()
    cat("           Calculating meta pct.1 & pct.2 values ...\n")
    med_pct.1=myPctAvgfn(pct_array=res_rearranged$pct.1,effective_size_array=res_rearranged$matEffectiveSize)
    med_pct.2=NULL
    med_logFC=NULL
    res_rearranged$pct.1=NULL
    
    cat("           Calculating meta z-scores ...\n")
    if(meta_method=="Stouffer"){
      meta_z=myStoufferFn(zscore_array = res_rearranged$zscore,weight_array = res_rearranged$matEffectiveSize,argList = argList)
      
      if(extendedMode){
        meta_z.8=myStoufferFn(zscore_array = z_q.8,weight_array = res_array$matEffectiveSize,argList = argList)
      }
      
    } else if(meta_method=="Harmony"){
      meta_z=myHarmonicPvalFn(zscore_array = res_array$zscore,weight_array = res_array$matWeights)
      meta_z[which(meta_z<(-10))]=(-10)
      if(extendedMode){
        meta_z.8=myHarmonicPvalFn(zscore_array = z_q.8,weight_array = res_array$matWeights)
        meta_z.8[which(meta_z.8<(-10))]=(-10)
      }
    } else if(meta_method=="average"){
      meta_z=myPctAvgfn(pct_array=res_array$zscore,centroid_weight_array=res_array$matWeights,effective_size_array=res_array$matEffectiveSize)
    }
    
    return(list(med_pct.1=med_pct.1,med_pct.2=med_pct.2,med_logFC=med_logFC,meta_z=meta_z))
    
  }
  
  require(qs)
  require(purrr)
  require(furrr)
  
  reRunCheck=T
  if(file.exists(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))&file.exists(.myFilePathMakerFn("res_DE_wZscore_pathwayAnalysis",argList=argList,uniformImportant=T,propImportant = T))&!argList$newRun){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T));F}, error=function(e) {return(T)})
    if(!reRunCheck){
      reRunCheck=tryCatch({load(.myFilePathMakerFn("res_DE_wZscore_pathwayAnalysis",argList=argList,uniformImportant=T,propImportant = T));F}, error=function(e) {return(T)})
      
    }
    if(!reRunCheck){
      reRunCheck=tryCatch({tst=qs::qread(.myFilePathMakerFn("res_meta",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T));F}, error=function(e) {return(T)})
    }
    
  }
  
  res=""
  if(reRunCheck){
    
    set.seed(analysis_seed)
    supportingFractionThr=argList$DE_supportingFractionThr
    n.adaptiveKernel=argList$DE_n.adaptiveKernel
    nPropIter=argList$DE_nPropIter
    
    #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    
    harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),]
    if(sum(is.na(harmony_embeddings))>0){
      stop("Error in matching Names!")
    }
    
    if(is.null(expData)){
      if(!file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))){
        stop("Expression data is missing!")
      }
      tmp=qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
      expData=SplitObject(tmp, split.by = "batch_merging")
    } else {
      tmp_pd=NULL
      for(i in 1:length(expData)){
        tmp_pd=rbind(tmp_pd,data.frame(sample=colnames(expData[[i]]),batch_merging=expData[[i]]$batch_merging,stringsAsFactors = F))
      }
      #pd=pd[row.names(pd) %in% tmp_pd$sample,]
      tmp_pd=tmp_pd[match(row.names(pd),tmp_pd$sample),]
      if(sum(is.na(tmp_pd$batch_merging))>0){
        stop("Error! phenoData doesn't match with the expression data")
      }
      pd$batch_merging=tmp_pd$batch_merging
    }
    
    pcaList=split(as.data.frame(harmony_embeddings),pd$batch_merging)
    
    if(argList$do.split.prop){
      dataArranged = parallel::mclapply(expData, function(thisExp){
        thisExp=thisExp[,colnames(thisExp) %in% row.names(pd)]
        if(class(thisExp)=="SingleCellExperiment"){
          tmp=counts(thisExp)
          #tmp=.extraExport2SeuratFn(thisExp)
          #plan("sequential")
          # plan("multicore", workers = min(parallel::detectCores(), 8, na.rm=T))
          #library(Seurat)
          #normalExp=Seurat::NormalizeData(tmp,verbose=F) # TODO JONAH Maybe not sequential
          # expData[[i]]=tmp
        }else{
          #normalExp = thisExp
          tmp=thisExp@assays$RNA@counts
        }
        
        #tmp=Seurat::NormalizeData(normalExp,normalization.method ="RC",verbose=F)
        
        return(list(
          countData=tmp,
          #logNormData=normalExp@assays$RNA@data,
          #countData=normalExp@assays$RNA@counts,
          #scaleData=tmp@assays$RNA@data,
          dsName=as.character(thisExp$batch_merging[1]),
          pcaData=pcaList[[as.character(thisExp$batch_merging[1])]]))
      })
    } else {
      dataArranged = parallel::mclapply(expData, function(thisExp){
        thisExp=thisExp[,colnames(thisExp) %in% row.names(pd)]
        if(class(thisExp)=="SingleCellExperiment"){
          tmp=counts(thisExp)
          #tmp=.extraExport2SeuratFn(thisExp)
          #plan("sequential")
          # plan("multicore", workers = min(parallel::detectCores(), 8, na.rm=T))
          #library(Seurat)
          #normalExp=Seurat::NormalizeData(tmp,verbose=F) # TODO JONAH Maybe not sequential
          # expData[[i]]=tmp
        }else{
          #normalExp = thisExp
          tmp=thisExp@assays$RNA@counts
        }
        
        #tmp=Seurat::NormalizeData(normalExp,normalization.method ="RC",verbose=F)
        
        return(list(
          countData=tmp,
          dsName=as.character(thisExp$batch_merging[1])))
      })
    }
    
    
    if(sum(unlist(lapply(dataArranged,function(x) ncol(x$countData))))<sum(unlist(lapply(expData,ncol)))){
      warning("Annotation was found for only a subset of exp data, that subset was only used in the analysis!")
    }
    
    res_fd=NULL
    for(i in 1:length(expData)){
      .printDot()
      if(class(expData[[i]])[1]=="SingleCellExperiment"){
        fd=as.data.frame(rowData(expData[[i]]))
      } else {
        fd=as.data.frame(expData[[i]]@assays$RNA@meta.features)
      }
      
      fd$ensembl_gene_id=gsub("_","-",fd$ensembl_gene_id)
      slCols=intersect(colnames(fd),c("gene_name","gene_biotype","symbol","gene_short_name","ensembl_gene_id"))
      fd=fd[,colnames(fd) %in% slCols]
      if(sum(is.na(fd$ensembl_gene_id))>0){
        fd$ensembl_gene_id[is.na(fd$ensembl_gene_id)]=row.names(fd)[is.na(fd$ensembl_gene_id)]
      }
      
      if(!is.null(res_fd)){
        fd=fd[!fd$ensembl_gene_id %in% res_fd$ensembl_gene_id,]
        if(nrow(fd)>0){
          res_fd=rbind(res_fd,fd)
        }
      } else {
        res_fd=fd
      }
    }
    fd=res_fd
    rm(res_fd,expData,pcaList,harmony_embeddings)
    gc()
    
    #if(!file.exists(do.call('.myFilePathMakerFn',args=c(list(filename="prop_mat_list_step2",uniformImportant=T),argList)))){
    cat("           Constructing the propagation matrices ...\n")
    dataArranged_old = dataArranged
    dataArranged = dataArranged %>% keep(~dim(.$countData)[[2]] >= argList$min_ds_size)
    if(length(dataArranged_old)!=length(dataArranged)){
      cat(paste0("Excluding datasets with less than ", argList$min_ds_size," cells => Retaining ", length(dataArranged), " out of ",
                   length(dataArranged_old)," datasets\n"))
    }
    rm(dataArranged_old)
    
    #centroidPCAdata=pca_centroid;nPropIter=1;n.neighbors=argList$prop.n.neighbors
    
    if(sum(is.na(pca_centroid))>0){
      stop("Error: NA values were found in the pca of pseudocells!")
    }
    
    if(argList$do.split.prop&F){
      res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final_v2_split_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
    } else {
      #batchPCAdata=NULL;dataArranged=dataArranged;centroidPCAdata=pca_centroid;argList=argList;exCentroids=NULL;runIndx=1;n.neighbors=argList$prop.n.neighbors;n.trees=50;NNmethod="annoy"
      res=.myConcensusDEFn_step2_detail_newprop3_final_v14(dataArranged=dataArranged,centroidPCAdata=pca_centroid,argList=argList,exCentroids=NULL,runIndx=1,batchPCAdata=NULL,n.neighbors=argList$prop.n.neighbors,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength)
    }
    
    #res_limited=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final_v2_limited,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
    
    #res1=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
    
    #Evaluating the performance of the propagation step
    #res_prop=res;argList=argList;annoCol="anno_cellState";collapse_datasets=F;return_plot=T
    #p=.sconline.anno2pseudocell_tmp(res_prop=res,argList=argList,annoCol="anno_cellState",collapse_datasets=T,return_plot=T,min_effective_size=5)
    #ggsave(plot=p$plot,file="~/myBucket/torm.pdf",width=20,height=20)
    #summary(.myEffSizePropMat(res[[1]]$prop_mat)$effective_sample_size/.myEffSizePropMat(res1[[1]]$prop_mat)$effective_sample_size)
    #tst_res=.myEffSizePropMat(res[[1]]$prop_mat)$effective_sample_size
    #tst_res_limited=.myEffSizePropMat(res_limited[[1]]$prop_mat)$effective_sample_size
    
    if(F){
      p=.sconline.anno2pseudocell_tmp(res_prop=res,argList=argList,annoCol="anno_cellState",collapse_datasets=T,return_plot=T,min_effective_size=5)
      
      tst=p$results
      tst=apply(tst[,-c(1,ncol(tst))],1,function(x){
        y=colnames(tst)[-c(1,ncol(tst))]
        y=y[which(x==max(x))[1]]
      })
      tst=which(grepl("purkinje",tolower(tst)))
      
      load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
      
      UMAP_centroid=.sconline.fetch_data("umap_centroids",argList=argList)
      
      pd_summary=.sconline.fetch_data("annotation",argList=argList)
      
      pd_summary=pd_summary[grepl("purkinje",tolower(pd_summary$anno_cellState)),]
      
      inputData=p$results[tst,]
      piechart_data=merge(inputData,UMAP_centroid,by.x="pseudocell",by.y="centroid",all=T)
      
      anno_cols=setdiff(colnames(inputData),c("pseudocell","effective_size"))
      pd_summary$anno_cellState=factor(as.character(pd_summary$anno_cellState),levels=anno_cols)
      p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=anno_cellState),size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                cols=anno_cols,r=2) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(anno_cols))))+scale_color_manual(values=c("gray",hues::iwanthue(length(anno_cols))))
      
      ggsave(plot=p,file="~/myBucket/torm.pdf",width=12,height = 12)
       
    }
    
    prop_mat=res$prop_mat
    res=res[["dataArranged"]]
    
    tmpValCheck=(unlist(lapply(res,length)))
    if(sum(tmpValCheck==1)>0){
      stop(res[[which(tmpValCheck==1)[1]]])
    }
    rm(tmpValCheck)
    
    ###Checking for outlier centroids
    
    
    
    
    cat("           Calculating dataset specific z-scores ...\n")
    #tst=.myConcensusDEFn_step2_detail_exp_final(res[[4]],argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode)
    #dim(tst$prop_mat)
    
    res2=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.cores = argList$ncores)
    while(sum(unlist(lapply(res2,class))=="NULL")>0){
      gc()
      null_ind=which(unlist(lapply(res2,class))=="NULL")
      cat(paste0("           Redoing the z-score calculations for ",length(null_ind)," datasets ...\n"))
      res3=parallel::mclapply(res[null_ind],.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.cores = 3)
      res2[null_ind]=res3
      rm(res3)
    }
    res=res2
    rm(res2)
    gc()
    
    #merging pseudocells in each batch
    
    
    
    
    tmpValCheck=(unlist(lapply(res,length)))
    if(sum(tmpValCheck==1)>0){
      stop(res[[which(tmpValCheck==1)[1]]])
    }
    rm(tmpValCheck)
    
    #qsave(dataArranged,file="~/torm_dataArranged.qs")
    res_meta=myMetaAnalysisFn(res=res,prop_mat=prop_mat,argList = argList,secondRun=F,extendedMode=F)
    
    med_pct.1=res_meta$med_pct.1
    #med_pct.2=res_meta$med_pct.2
    #med_logFC=res_meta$med_logFC
    meta_z=res_meta$meta_z
    #res_array=res_meta$res_array
    rm(res_meta,res)
    gc()
    
    
    #three_d_image = tf.concat(0, [[r], [g], [b]])
    #three_d_image = tf.stack([r,g,b], axis=2)
    
    
    ##################
    
    #MetaQC: initial implementation
    
    if(F){
      require(GeneMeta)
      meta_z=matrix(0,nrow=dim(res_array$t_adj)[2],ncol=dim(res_array$t_adj)[3])
      for(i in 1:dim(res_array$t_adj)[2]){
        tmp_t=t(res_array$t_adj[,i,])
        tmp_sigmad=t(res_array$sigmad_adj[,i,])
        
        empty_ind=which(colSums(tmp_sigmad)==0|colSums(tmp_t)==0)
        if(length(empty_ind)<(ncol(tmp_t)-1)){
          if(length(empty_ind)>0){
            tmp_t=tmp_t[,-empty_ind]
            tmp_sigmad=tmp_sigmad[,-empty_ind]
          }
          
          my.Q   <- f.Q(tmp_t, tmp_sigmad)
          
          my.tau2.DL<-tau2.DL(my.Q, ncol(tmp_t), my.weights=1/tmp_sigmad)
          #obtain new variances s^2+tau^2
          myvarsDL <- tmp_sigmad + my.tau2.DL
          #compute
          muREM <- mu.tau2(tmp_t, myvarsDL)
          #cumpute mu(tau)
          varREM <- var.tau2(myvarsDL)
          ZREM <- muREM/sqrt(varREM)
          meta_z[i,]=ZREM
        } else if(length(empty_ind)==(ncol(tmp_t)-1)) {
          meta_z[i,]=res_array$zscore[-empty_ind,i,]
        }
        
        
      }
    }
    
    if(F){
      .mycheck.R=function (R, checksym = TRUE, checkna = TRUE, checkpd = FALSE, 
                           nearpd = FALSE, checkcor = FALSE, checkdiag = TRUE, isbase = TRUE, 
                           k, adjust, fun) {
        if (nearpd) 
          checkpd <- TRUE
        if (checkpd) {
          checkna <- TRUE
          checksym <- TRUE
        }
        if (inherits(R, "dpoMatrix")) 
          R <- as.matrix(R)
        if (inherits(R, "data.frame")) 
          R <- as.matrix(R)
        if (checksym && !(is.matrix(R) && isSymmetric(unname(R)))) 
          stop("Argument 'R' must be a symmetric matrix.", call. = FALSE)
        if (checkna && any(is.na(R))) 
          stop("Values in 'R' must not contain NAs.", call. = FALSE)
        if (checkpd && any(eigen(R)$values <= 0)) {
          if (nearpd) {
            warning("Matrix 'R' is not positive definite. Used Matrix::nearPD() to make 'R' positive definite.", 
                    call. = FALSE)
            R <- as.matrix(.find.nonegmat(R))
          }
          else {
            stop("Matrix 'R' is not positive definite.", call. = FALSE)
          }
        }
        if (checkcor && any(abs(R) > 1, na.rm = TRUE)) 
          stop("Argument 'R' must be a correlation matrix, but contains values outside [-1,1].", 
               call. = FALSE)
        if (checkdiag && (any(is.na(diag(R))) || any(diag(R) != 1))) 
          stop("Diagonal values in 'R' must all be equal to 1.", 
               call. = FALSE)
        if (isbase) {
          if (k != nrow(R)) 
            stop("Length of 'p' vector (", k, ") does not match the dimensions of the 'R' matrix (", 
                 nrow(R), ",", ncol(R), ").", call. = FALSE)
          if (adjust == "none") 
            warning("Argument 'R' was specified, but no adjustment method was chosen via the 'adjust' argument.\nTo account for dependence, specify an adjustment method. See help(", 
                    fun, ") for details.", call. = FALSE)
          if (adjust == "user") 
            warning("When 'm' is specified, argument 'R' is irrelevant and ignored.")
        }
        return(R)
      }
      
      .myMefffn=function (R, eigen, method="nyholt", ...) {
        .is.numeric.vector=function (x){
          is.atomic(x) && is.numeric(x) && !is.matrix(x) && !is.null(x)
        }
        
        
        method <- match.arg(method, c("nyholt", "liji", "gao", "galwey"))
        if (missing(eigen)) {
          if (missing(R)) 
            stop("Argument 'R' must be specified.", call. = FALSE)
          R <- .mycheck.R(R, checksym = TRUE, checkna = TRUE, checkpd = FALSE, 
                          nearpd = FALSE, checkcor = TRUE, checkdiag = TRUE, 
                          isbase = FALSE)
          evs <- base::eigen(R)$values
        } else {
          if (!.is.numeric.vector(eigen)) 
            stop("Argument 'eigen' must be a numeric vector.", 
                 call. = FALSE)
          evs <- eigen
        }
        if (any(evs < 0)) 
          warning(paste0("One or more eigenvalues ", ifelse(missing(eigen), 
                                                            "derived from the 'R' matrix ", ""), "are negative."), 
                  call. = FALSE)
        if (method == "nyholt") {
          k <- length(evs)
          m <- 1 + (k - 1) * (1 - var(evs)/k)
        }
        if (method == "liji") {
          abs.evs <- abs(evs) + sqrt(.Machine$double.eps)
          m <- sum(ifelse(abs.evs >= 1, 1, 0) + (abs.evs - floor(abs.evs)))
        }
        if (method == "gao") {
          ddd <- list(...)
          if (!is.null(ddd$C)) {
            C <- ddd$C
          }
          else {
            C <- 0.995
          }
          if (C < 0 || C >= 1) 
            warning("Value of 'C' should be >= 0 and < 1.", call. = FALSE)
          m <- which(cumsum(sort(evs, decreasing = TRUE))/sum(evs) > 
                       C)[1]
        }
        if (method == "galwey") {
          if (any(evs < 0)) {
            warning(paste0("Negative eigenvalues ", ifelse(missing(eigen), 
                                                           "derived from the 'R' matrix ", ""), "were set to 0."), 
                    call. = FALSE)
            evs[evs < 0] <- 0
          }
          m <- sum(sqrt(evs))^2/sum(evs)
        }
        m <- floor(m)
        return(m)
      }
      
      myStoufferFn_org=function(zscore_array,weight_array){
        z_w=apply(zscore_array*weight_array, c(2,3), sum)
        w2=sqrt(apply(weight_array*weight_array, c(2,3), sum))
        if(sum(w2>0&w2<0.99)>0){
          #warning("Possible issue in the weights")
          w2[which(w2>0&w2<1)]=1
        }
        z_w=z_w/w2
        z_w[which(w2==0)]=0
        return(z_w)
      }
    }
    
    
    if(F){
      #alternate slower solution
      myPctAvgfn=function(pct_array,centroid_weight_array,effective_size_array){
        myWeightedMeanFn=function(data_array,weight_array){
          z_w=multiApply::Apply(data_array*weight_array, 1, function(x) sum(x,na.rm = T))[[1]]
          w=multiApply::Apply(weight_array, 1, function(x) sum(x,na.rm = T))[[1]]
          if(sum(w>0&w<0.99)>0){
            #warning("Possible issue in the weights")
            w[which(w>0&w<1)]=1
          }
          z_w=z_w/w
          z_w[which(w==0)]=0
          return(z_w)
        }
        
        pct_array[centroid_weight_array<0.9]=NA
        med_pct.1=myWeightedMeanFn(data_array=pct_array,weight_array=effective_size_array)
        row.names(med_pct.1)=dimnames(pct_array)[[2]]
        return(med_pct.1)
        
      }
      
      tst=myPctAvgfn(pct_array=res_array$pct.1,centroid_weight_array=res_array$matWeights,effective_size_array=res_array$matEffectiveSize)
      
    }
    
    
    
    #tst=list(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F)
    #pct_mat=med_pct.1;meta_z_mat=meta_z;pct_mat_ref=NULL;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
    
    pct_diff_count=.extra_sconline.PctScoreFn(pct_mat=med_pct.1,argList = argList,meta_z_mat=meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    
    
    
    diff=prop_mat
    diff@x=rep(1,length(diff@x))
    diff=diff %*% t(diff)
    diff=sweep(diff,1,diag(diff),"/")
    diff=as.matrix(diff)
    diff[diff<0.100001]=0.100001
    diff=abs(log10(diff))
    diff=diff+pct_diff_count
    
    diff=hclust(as.dist(diff),method = "complete")
    diff=cutree(diff,h=0.999)
    
    if(F){
      diff=data.frame(pseudocell=names(diff),cluster=diff,stringsAsFactors = F)
      adj2=prop_mat
      load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
      pd=pd[colnames(adj2),]
      
      table(is.na(pd$anno_cellState))
      tmp= adj2 %*% as.matrix(.myOneHotFn(pd$anno_cellState))
      
      tmp=apply(tmp,1,function(x){
        y=which(x==max(x))[1]
        tmp2=data.frame(cluster=colnames(tmp)[y],purity=x[y])
        tmp2
      })
      
      tmp=do.call("rbind",tmp)
      tmp2=aggregate(purity~cluster,data=tmp,median)
      
      table(tmp2[,2]>0.8)
      
      
      
      diff=merge(diff,data.frame(pseudocell=row.names(tmp),tmp,stringsAsFactors = F),by="pseudocell")
      diff2=aggregate(cluster.y~cluster.x,data=diff,function(x){
        y=as.numeric(table(x))
        max(y)/sum(y)
      })
    }
    
    
    if(sum(duplicated(diff))>0){
      cat(paste0("           Merging pseudocells based on their DE count; retaining ",length(unique(diff))," out of ",length(diff),"\n"))
      diff=data.frame(pseudocell=names(diff),cluster=diff,stringsAsFactors = F)
      diff_id=diff[!duplicated(diff$cluster),]
      diff=merge(diff,diff_id,by="cluster")
      diff2=reshape2::dcast(pseudocell.x~pseudocell.y,data=diff,fun.aggregate=length)
      diff2=diff2[match(row.names(prop_mat),diff2[,1]),]
      diff2=t(as.matrix(diff2[,-1]))
      
      diff2=as(diff2,"dgCMatrix")
      diff2=diff2 %*% prop_mat
      diff2 <- Matrix::Diagonal(x = 1 / (rowSums(diff2)+0.000000000001)) %*% diff2
      
      .diff2=diff2
      
      diff2=.diff2
      
      rowMeans_drop0 <- function (dgCMat) {
        RowInd <- dgCMat@i + 1
        sapply(split(dgCMat@x, RowInd), function(x)quantile(x,0.95))
      }
      
      diff2 =  Matrix::Diagonal(x = 1 / (rowMeans_drop0(diff2))) %*% diff2
      #diff2@x=pmin(diff2@x,1)
      colMax_vals=c()
      for(i in seq(1,ncol(diff2),5000)){
        tmp_max=as.numeric(qlcMatrix::colMax(diff2[,i:min(i+4999,ncol(diff2))]))
        colMax_vals=c(colMax_vals,as.numeric(tmp_max))
      }
      diff2 =  diff2 %*% Matrix::Diagonal(x = 1 / colMax_vals)
      diff2 <- Matrix::Diagonal(x = 1 / rowSums(diff2)) %*% diff2
      
      prop_mat=diff2
      
      matWeights=.myEffSizePropMat(diff2)
      
      matEffectiveSize=matWeights$effective_sample_size
      matWeights=matWeights$centroid_weights
      matWeights=matWeights[match(row.names(diff2),names(matWeights))]
      matEffectiveSize=matEffectiveSize[match(row.names(diff2),names(matEffectiveSize))]
      
      res=dataArranged
      for(i in 1:length(res)){
        tmp_prop_mat=diff2[,match(colnames(res[[i]]$countData),colnames(diff2))]
        #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
        tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
        tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
        tmp_weights[tmp_effsize<4]=0
        tmp=list(prop_mat=tmp_prop_mat,data=res[[i]],matWeights=tmp_weights,matEffectiveSize=tmp_effsize)
        res[[i]]=tmp
      }
      
      gc()
      plan("sequential")
      cat(paste0("           Redoing the propagation after the merging\n"))
      res=parallel::mclapply(res,.myConcensusDEFn_step2_detail_exp_final,argList=argList,sd_offset=0.005,reg_sd=T,extendedMode=extendedMode,mc.cores = 4)
      
      cat(paste0("           Redoing the meta analysis of zscores\n"))
      res_meta=myMetaAnalysisFn(res=res,prop_mat=prop_mat,argList = argList,secondRun=T,extendedMode=extendedMode)
      
      med_pct.1=res_meta$med_pct.1
      #med_pct.2=res_meta$med_pct.2
      #med_logFC=res_meta$med_logFC
      meta_z=res_meta$meta_z
      #res_array=res_meta$res_array
      rm(res_meta)
      gc()
    }
    
    
    de_pct_res=NULL
    if(addClusteringModule&F){
      
      cosine_dist=.extra_sconline.CosineDistFn(inputMatrix = meta_z,sig_thr = 2)
      
      de_pct_res=.extra_sconline.PctScoreFn(argList = argList,pct_mat=med_pct.1,meta_z_mat=meta_z,pct_mat_ref=NULL,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F)
      
    }
    #meta_z_mat=meta_z$meta_z;cosine_dist=cosine_dist;de_dist=de_pct_res;pct_mat=meta_z$med_pct.1;min_marker_thr=20;sig1_thr=3;sig2_thr=1;pct_de_count_thr=1;pct_diff_thr=0.2;pct2_thr=0.3
    #####################
    
    
      
    #####################
    
    
    if(extendedMode){
      res_meta=list(meta_z=meta_z,meta_z.8=meta_z.8,med_pct.1=med_pct.1,med_pct.2=med_pct.2,med_logFC=med_logFC,med_n=med_n,fd=fd)
    } else {
      res_meta=list(meta_z=meta_z,med_pct.1=med_pct.1,fd=fd)
    }
    
    qsave(res_meta,file=.myFilePathMakerFn("res_meta",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
    rm(res_meta)
    
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    
    
    x_z=apply(meta_z,2,max)
    x_z=meta_z[,which(x_z>2)]
    x_z=Seurat::RunPCA(as.matrix(t(x_z)),npcs = min(max(50,argList$nPCs+10),nrow(x_z)-2),verbose = FALSE)
    network=.myConcensusDEFn_step2_FindNeighbors(inputCentroids = x_z@cell.embeddings,argList=argList,verbose = F)
    
    snnNet=as.matrix(network$snn)
    ind=which(snnNet>0,arr.ind = T)
    
    snnNet=data.frame(Fnode=row.names(snnNet)[ind[,1]],Snode=row.names(snnNet)[ind[,2]],weight=snnNet[ind],stringsAsFactors = F)
    snnNet=snnNet[snnNet$Fnode!=snnNet$Snode,]
    snnNet$weight=1-snnNet$weight
    snnNet2=igraph::graph_from_data_frame(snnNet, directed = F, vertices = NULL)
    
    snnNet2=igraph::distances(snnNet2, mode = c("all"), weights = snnNet$weight, algorithm = "dijkstra")
    if(length(which(snnNet2==Inf))>0){
      snnNet2[which(snnNet2==Inf)]=(max(snnNet2[snnNet2!=Inf],na.rm = T)+1)
    }
    
    geneList=colnames(meta_z)
    
    
    cat("           Calculating cor and MI ...\n")
    resCorMI=NULL
    slGenes=apply(meta_z,2,function(x) max(abs(x)))
    slGenes=names(slGenes)[which(slGenes>2)]
    
    myCorMIfn=function(inputGenes,z_score_mat,snnNet2){
      resCorMI=NULL
      for(i in inputGenes){
        tmp=z_score_mat[,which(colnames(z_score_mat)==i)]
        tmpNet=snnNet2[row.names(snnNet2) %in% names(tmp),colnames(snnNet2) %in% names(tmp)]
        tmp=tmp[match(colnames(tmpNet),names(tmp))]
        tmp=data.frame(zscore=tmp,centroid=names(tmp),gene=i,stringsAsFactors = F)
        if(max(abs(tmp$zscore))>2){
          zOrdered=tmp$zscore[order(tmp$zscore,decreasing = T)]
          
          tmpInd=which(row.names(tmpNet) %in% tmp$centroid[which(tmp$zscore>=max(zOrdered[3],1))])
          if(length(tmpInd)>0){
            if(length(tmpInd)>1){
              distMat=tmpNet[tmpInd,]
              distMat=apply(distMat,2,mean)
            } else {
              distMat=tmpNet[tmpInd,]
            }
            
            tmp=data.frame(distance=distMat,heat=tmp$zscore)
            tmpCor=cor(tmp$distance,tmp$heat)
            
            resCorMI=rbind(resCorMI,data.frame(gene=i,cor=tmpCor,stringsAsFactors = F))
            
          } else{
            resCorMI=rbind(resCorMI,data.frame(gene=i,cor=0,stringsAsFactors = F))
          }
        }
      }
      return(resCorMI)
    }
    
    if(argList$ncores>1){
      geneList=split(as.character(slGenes), cut(1:length(as.character(slGenes)), argList$ncores, labels = FALSE)) 
    } else {
      geneList=list(as.character(slGenes))
    }
    
    resCor=parallel::mclapply(geneList,myCorMIfn,z_score_mat=meta_z,snnNet2=snnNet2,mc.cores = argList$ncores)
    resCor=do.call("rbind",resCor)
    
    slInd=which(abs(meta_z)>2,arr.ind = T)
    if(extendedMode){
      if((!all(row.names(meta_z)==row.names(meta_z.8)))|(!all(colnames(meta_z)==colnames(meta_z.8)))){
        stop("Error in matching!")
      }
      if((!all(row.names(meta_z)==row.names(med_n)))|(!all(colnames(meta_z)==colnames(med_n)))){
        stop("Error in matching!")
      }
    }
    
    
    if((!all(row.names(meta_z)==row.names(med_pct.1)))|(!all(colnames(meta_z)==colnames(med_pct.1)))){
      stop("Error in matching!")
    }
    
    
    if(extendedMode){
      res_arranged=data.frame(gene=colnames(meta_z)[slInd[,2]],centroid=row.names(meta_z)[slInd[,1]],zscore=meta_z[slInd],zscore.8=meta_z.8[slInd],pct.1_median=med_pct.1[slInd],pct.2_median=med_pct.2[slInd],logFC_median=med_logFC[slInd],n_median=med_n[slInd],stringsAsFactors = F)
    } else {
      res_arranged=data.frame(gene=colnames(meta_z)[slInd[,2]],centroid=row.names(meta_z)[slInd[,1]],zscore=meta_z[slInd],pct.1_median=med_pct.1[slInd],stringsAsFactors = F)
    }
    
    res_arranged=merge(res_arranged,resCor,by="gene",all.x=T)
    
    #zscore_count measures in how many centroids each gene as a nominal significance
    zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(x>2))
    colnames(zscore_counts)=c("gene","zscore_count")
    res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
    zscore_counts=aggregate(zscore~gene,data=res_arranged,function(x) sum(abs(x)>2))
    colnames(zscore_counts)=c("gene","zscore_count_abs")
    res_arranged=merge(res_arranged,zscore_counts,by="gene",all.x=T)
    if(sum(res_arranged$zscore_count>0)>0){
      res_arranged=res_arranged[res_arranged$zscore_count>0,]
    }
    
    res_arranged=res_arranged[order(res_arranged$zscore,decreasing = T),]
    
    save(res_arranged,file=.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))
    
    qsave(prop_mat,file=.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    
    
    if(F){
      prop_mat_list=list()
      for(i in 1:length(res)){
        prop_mat_list=c(prop_mat_list,res[[i]]$prop_mat)
        names(prop_mat_list)[length(prop_mat_list)]=res[[i]]$dsName
      }
      
      pseudocell_sim_list=list()
      for(i in 1:length(res)){
        pseudocell_sim_list=c(pseudocell_sim_list,res[[i]]$pseudocell_sim_mat)
        names(pseudocell_sim_list)[length(pseudocell_sim_list)]=res[[i]]$dsName
      }
      
      #res_array$prop_mat=prop_mat_list
      #res_array$pseudocell_sim_list=pseudocell_sim_list
      
    }
    
    
    
    
    
    if(F){
      qsave(res_array,file=.myFilePathMakerFn("res_dataset_array",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
    }
    
    
    
    gc()
    
  }
  
  return("Done")
}

.myConcensusDEFn_step2_detail_exp_final_org=function(inputData,argList,sd_offset=0.001,reg_sd=T,zscore_cap=10,minCellCountThr=4,extendedMode=F){
  
  #inputData=res[[1]];sd_offset=0.005;reg_sd=T;zscore_cap=10;minCellCountThr=4;extendedMode=F
  if(F){
    #Jonah's
    write.table(data.frame(), file=paste0("/tmp/touch_", as.character(runif(1)*10^20)),
                col.names=FALSE)
    
    library("RhpcBLASctl")
    omp_set_num_threads(4)
    blas_set_num_threads(4)
    
  }
  
  
  #inputData=res[[1]];sd_offset=0.001;reg_sd=T;zscore_cap=10;minCellCountThr=4
  
  #tmp_Propweights=rowSums(as.matrix(inputData$prop_mat))
  
  if(sum(names(argList)=="do.split.prop")==0){
    argList$do.split.prop=T
  }
  
  if(!argList$do.split.prop){
    inputData$prop_mat=Matrix::Diagonal(x = 1 / (rowSums(inputData$prop_mat)+0.000000000001)) %*% inputData$prop_mat
  }
  
  
  require(Matrix,quietly = T)
  
  zscore_cap=abs(zscore_cap)
  if(zscore_cap<2){
    warning("zscore cap was set to below 2. It was changed to 5.")
    zscore_cap=5
  }
  sd_offset=max(sd_offset,0)
  
  resGeneMeanSd=.extra_gmean(inputData$data$countData)
  resGeneMeanSd=data.frame(gene=names(resGeneMeanSd),geo_mean_count=resGeneMeanSd,stringsAsFactors = F)
  
  
  myWeightedVar=function(inputWeight,inputExp,min_cell_count){
    require(Matrix,quietly = T)
    inputExp=t(inputExp)
    res_binary=(sign(inputWeight)%*% sign(inputExp))
    exp_sq=inputExp^2
    res=(inputWeight%*% exp_sq)
    res=sweep(res,1,rowSums(inputWeight),"*")
    res=res- (inputWeight%*% inputExp)^2
    res=sweep(res,1,rowSums(inputWeight)^2,"/")
    
    singletons=apply(inputWeight,1,function(x) sum(x>0))
    singletons[singletons<=min_cell_count]=0
    singletons[singletons>0]=1
    if(sum(singletons<1)>0){
      res=sweep(res,1,singletons,"*")
    }
    res[res_binary<=min_cell_count]=0
    return(res)
  }
  
  prop_mat=inputData$prop_mat
  
  
  prop_mat_c=prop_mat
  if(quantile(rowSums(prop_mat_c > 0,na.rm = T), 0.25) < (0.85*ncol(prop_mat))){
    prop_mat_c[prop_mat_c>0]=1
  }
  
  prop_mat_c=1-prop_mat_c
  prop_mat_c <- Matrix::Diagonal(x = 1 / (rowSums(prop_mat_c))) %*% prop_mat_c
  
  logNormData=Seurat:::NormalizeData.default(inputData$data$countData,normalization.method = "LogNormalize",verbose = F)
  
  
  x_exp=prop_mat %*% t(logNormData)
  if(extendedMode){
    x2_exp=Seurat:::NormalizeData.default(inputData$data$countData,normalization.method = "RC",verbose = F)
    x2_exp=prop_mat %*% t(x2_exp)
  }
  
  y_exp=prop_mat_c %*% t(logNormData)
  
  
  #inputWeight = inputData$prop_mat;inputExp = inputData$data$logNormData
  var_x=myWeightedVar(inputWeight = inputData$prop_mat,inputExp = logNormData,min_cell_count=minCellCountThr)
  var_y=myWeightedVar(inputWeight = prop_mat_c,inputExp = logNormData,min_cell_count=minCellCountThr)
  
  n_x=.myEffSizePropMat(prop_mat = prop_mat)
  n_x=n_x$effective_sample_size
  
  n_y=.myEffSizePropMat(prop_mat = prop_mat_c)
  n_y=n_y$effective_sample_size
  
  vxn=sweep(var_x,1,n_x,"/")
  vyn=sweep(var_y,1,n_y,"/")
  
  vxn2=sweep(var_x,1,n_x-1,"*")
  vyn2=sweep(var_y,1,n_y-1,"*")
  
  
  df=(vxn + vyn)^2
  df=df/(sweep(vxn^2,1,(n_x - 1),"/") + 
           sweep(vyn^2,1,n_y - 1,"/"))
  if(sum(vxn==0,na.rm = T)>0){
    df[Matrix::which(vxn==0)]=matrix(n_y,nrow=nrow(vxn),ncol=ncol(vxn),byrow = F)[Matrix::which(vxn==0)]
  }
  
  vxn[is.na(as.matrix(vxn))]=1
  sxy=sqrt(vxn+vyn)
  
  d_sxy=sqrt(sweep(vxn2+vyn2,1,n_x+n_y-2,"/"))
  rm(vxn2,vyn2)
  
  sd_reg=sxy
  d_sd_reg=d_sxy
  if(reg_sd){
    
    if(F){
      #Vahid's version
      sd_reg=lapply(1:nrow(sxy),function(x){
        y=sxy[x,]
        tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=sxy[x,])
        slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
        if(length(slInd)>0){
          tmp_data=tmp_data[slInd,]
          meanExp=tmp_data$exp
          design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
          fit2=lm.fit(design,log10(tmp_data$sd))
          tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
          
          y[slInd]=tmp_data$fit2
        } else {
          y=pmax(y,1)
        }
        
        return(as.numeric(y))
      })
      sd_reg=do.call("rbind",sd_reg)
    }
    
    
    #Jonah's version
    sd_reg=apply(sxy,1, function(y){
      # browser()
      # y=sxy[x,]
      tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=y)
      slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
      if(length(slInd)>0){
        tmp_data=tmp_data[slInd,]
        meanExp=tmp_data$exp
        design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
        fit2=lm.fit(design,log10(tmp_data$sd))
        tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
        
        y[slInd]=tmp_data$fit2
      } else {
        y=pmax(y,1)
      }
      
      return(as.numeric(y))
    })
    ## browser()
    sd_reg=t(sd_reg)
    
    
    if(F){
      #Vahid's version
      d_sd_reg=lapply(1:nrow(d_sxy),function(x){
        y=d_sxy[x,]
        tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=d_sxy[x,])
        slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
        if(length(slInd)>0){
          tmp_data=tmp_data[slInd,]
          meanExp=tmp_data$exp
          design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
          fit2=lm.fit(design,log10(tmp_data$sd))
          tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
          
          y[slInd]=tmp_data$fit2
        } else {
          y=pmax(y,1)
        }
        
        return(as.numeric(y))
      })
      d_sd_reg=do.call("rbind",d_sd_reg)
    }
    
    
    d_sd_reg=apply(d_sxy, 1, function(y){
      # y=d_sxy[x,]
      tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=y)
      slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
      if(length(slInd)>0){
        tmp_data=tmp_data[slInd,]
        meanExp=tmp_data$exp
        design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
        fit2=lm.fit(design,log10(tmp_data$sd))
        tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
        
        y[slInd]=tmp_data$fit2
      } else {
        y=pmax(y,1)
      }
      
      return(as.numeric(y))
    })
    d_sd_reg = t(d_sd_reg)
    # d_sd_reg=do.call("rbind",d_sd_reg)
    
  }
  rm(d_sxy)
  
  sd_reg=sxy
  sd_offset=0
  
  
  t <- (x_exp - y_exp)/(sd_reg+sd_offset)
  t[which(vxn==0)]=0
  if(F){
    t_adj=sweep(t,1,sqrt((n_x + n_y)/(n_x * n_y)),"*")
    t_adj=t_adj- sweep((3 * t_adj),1,(4 * ((n_x+n_y) - 2) - 1),"/")
    sigmad_adj=1/n_x+1/n_y+sweep((t_adj^2),1,(2*(n_x+n_y)),"/")
  }
  
  if(extendedMode){
    cohen_d=(x_exp - y_exp)/d_sd_reg
    hodge_g=cohen_d*(1-3/(4*(n_x+n_y)-9))
    if(sum(is.na(hodge_g))>0){
      hodge_g[is.na(hodge_g)]=0
    }
    
    se.g=0.5*sweep(hodge_g^2,1,n_x+n_y-3.94,"/")
    se.g <- sqrt(sweep(se.g,1, (n_x+n_y)/(n_x*n_y),"+"))
    if(sum(is.na(se.g))>0){
      se.g[is.na(se.g)]=0
    }
    
  }
  
  
  
  #t2 <- (x_exp - y_exp)/(sd_reg)
  
  t_vec=as.numeric(t)
  zscore=qnorm(pt(as.numeric(abs(t_vec)), as.numeric(df),lower.tail = F,log.p = T),lower.tail = F,log.p = T)
  zscore[is.na(zscore)]=0
  zscore=zscore*sign(t_vec)
  zscore=matrix(zscore,nrow=nrow(t),ncol=ncol(t))
  
  colnames(zscore)=row.names(logNormData)
  row.names(zscore)=row.names(inputData$prop_mat)
  
  exp_binary=t(logNormData)
  exp_binary@x=rep(1,length(exp_binary@x))
  
  pct.1=as.matrix(prop_mat %*% exp_binary)
  pct.2=as.matrix(prop_mat_c %*% exp_binary)
  
  exp_norm=t(expm1(logNormData))
  
  exp.1=log2(prop_mat %*% exp_norm+1)
  exp.2=log2(prop_mat_c %*% exp_norm+1)
  
  logFC=as.matrix(exp.1 - exp.2)
  
  if(extendedMode){
    n=sqrt(sweep(pct.1,1,n_x,"*")*sweep(pct.2,1,n_y,"*"))
  }
  
  
  if(sum(zscore>zscore_cap,na.rm = T)>0){
    zscore[which(zscore>zscore_cap)]=zscore_cap
  }
  
  if(sum(zscore<(-1*zscore_cap),na.rm = T)>0){
    zscore[which(zscore<(-1*zscore_cap))]=(-1*zscore_cap)
  }
  
  
  if(F){
    #moved to the propagation step
    matWeights=.myEffSizePropMat(prop_mat)
    
    matEffectiveSize=matWeights$effective_sample_size
    matWeights=matWeights$centroid_weights
    matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
    matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
    
  }
  
  
  matWeights=matrix(inputData$matWeights,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
  matEffectiveSize=matrix(inputData$matEffectiveSize,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
  row.names(matWeights)=row.names(matEffectiveSize)=row.names(prop_mat)
  colnames(matWeights)=colnames(matEffectiveSize)=colnames(zscore)
  matWeights[which(pct.1*matEffectiveSize<minCellCountThr&pct.2<0.001)]=0
  matWeights[which(zscore==0|var_x==0)]=0
  
  matEffectiveSize[matWeights<0.0001]=0
  if(extendedMode){
    res=list(zscore=zscore,pct.1=pct.1,pct.2=pct.2,logFC=logFC,se.g=se.g,hodge_g=hodge_g,n=n,matWeights=matWeights,matEffectiveSize=matEffectiveSize,exp1=x_exp,exp2=x2_exp)
  } else {
    res=list(zscore=zscore,pct.1=pct.1,pct.2=pct.2,logFC=logFC,matWeights=matWeights,matEffectiveSize=matEffectiveSize,exp1=x_exp)
  }
  if(F){
    #Jonah's
    output = c(res,list(dsName=inputData$data$dsName,prop_mat=inputData$prop_mat,pseudocell_sim_mat=inputData$pseudocell_sim_mat))
    
    
    res_non_sparse = res
    for(dsitr in c("zscore", "pct.1", "pct.2", "logFC", "n", "matWeights", "matEffectiveSize")){
      res[[dsitr]] = as.sparse(res[[dsitr]])
    }
    gc()
    
    dir.create(file.path(argList$saveDir, "detail_exp_final"), showWarnings = FALSE)
    FN = file.path(argList$saveDir, "detail_exp_final", paste0(inputData$data$dsName, ".RDS"))
    
    .mcsaveRDS <- function(object,file,mc.cores=min(parallel::detectCores(),3, na.rm=T)) {
      con <- pipe(paste0("pigz -p",mc.cores," > ",file),"wb")
      saveRDS(object, file = con)
      close(con)
    }
    print(FN)
    
    .mcsaveRDS(output, FN)
    
    
  }
  
  setdiff(names(inputData) ,c("matWeights","matEffectiveSize"))
  return(c(res,list(dsName=inputData$data$dsName)))
  #return(T)
}

.myConcensusDEFn_step2_detail_exp_final=function(inputData,argList,sd_offset=0.001,reg_sd=T,zscore_cap=10,minCellCountThr=4,extendedMode=F){
  
  #inputData=res[[1]];sd_offset=0.005;reg_sd=T;zscore_cap=10;minCellCountThr=4;extendedMode=F
  if(F){
    #Jonah's
    write.table(data.frame(), file=paste0("/tmp/touch_", as.character(runif(1)*10^20)),
                col.names=FALSE)
    
    library("RhpcBLASctl")
    omp_set_num_threads(4)
    blas_set_num_threads(4)
    
  }
  #require(bigmemory)
  
  #inputData=res[[1]];sd_offset=0.001;reg_sd=T;zscore_cap=10;minCellCountThr=4
  
  #tmp_Propweights=rowSums(as.matrix(inputData$prop_mat))
  
  if(sum(names(argList)=="do.split.prop")==0){
    argList$do.split.prop=T
  }
  
  if(!argList$do.split.prop){
    inputData$prop_mat=Matrix::Diagonal(x = 1 / (rowSums(inputData$prop_mat)+0.000000000001)) %*% inputData$prop_mat
  }
  
  
  require(Matrix,quietly = T)
  
  zscore_cap=abs(zscore_cap)
  if(zscore_cap<2){
    warning("zscore cap was set to below 2. It was changed to 5.")
    zscore_cap=5
  }
  sd_offset=max(sd_offset,0)
  
  resGeneMeanSd=.extra_gmean(inputData$data$countData)
  resGeneMeanSd=data.frame(gene=names(resGeneMeanSd),geo_mean_count=resGeneMeanSd,stringsAsFactors = F)
  
  
  myWeightedVar=function(inputWeight,inputExp,min_cell_count){
    require(Matrix,quietly = T)
    #inputExp=t(inputExp)
    res_binary=(sign(inputWeight)%*% sign(inputExp))
    exp_sq=inputExp^2
    res=(inputWeight%*% exp_sq)
    res=sweep(res,1,rowSums(inputWeight),"*")
    res=res- (inputWeight%*% inputExp)^2
    res=sweep(res,1,rowSums(inputWeight)^2,"/")
    
    singletons=apply(inputWeight,1,function(x) sum(x>0))
    singletons[singletons<=min_cell_count]=0
    singletons[singletons>0]=1
    if(sum(singletons<1)>0){
      res=sweep(res,1,singletons,"*")
    }
    res[res_binary<=min_cell_count]=0
    return(res)
  }
  
  prop_mat=inputData$prop_mat
  
  
  prop_mat_c=prop_mat
  if(quantile(rowSums(prop_mat_c > 0,na.rm = T), 0.25) < (0.85*ncol(prop_mat))){
    prop_mat_c=Matrix::drop0(prop_mat_c)
    prop_mat_c@x=rep(1,length(prop_mat_c@x))
  }
  
  prop_mat_c=Matrix::drop0(1-prop_mat_c)
  prop_mat_c <- Matrix::Diagonal(x = 1 / (rowSums(prop_mat_c))) %*% prop_mat_c
  
  logNormData=t(Seurat:::NormalizeData.default(inputData$data$countData,normalization.method = "LogNormalize",verbose = F))
  logNormData=logNormData[colnames(prop_mat),]
  
  x_exp=prop_mat %*% logNormData
  if(extendedMode){
    x2_exp=Seurat:::NormalizeData.default(inputData$data$countData,normalization.method = "RC",verbose = F)
    x2_exp=prop_mat %*% t(x2_exp)
  }
  
  gene_name_list=colnames(logNormData)
  pseudocell_name_list=row.names(prop_mat)
  
  #inputWeight = inputData$prop_mat;inputExp = inputData$data$logNormData
  var_x=myWeightedVar(inputWeight = inputData$prop_mat,inputExp = logNormData,min_cell_count=minCellCountThr)
  var_y=myWeightedVar(inputWeight = prop_mat_c,inputExp = logNormData,min_cell_count=minCellCountThr)
  
  n_x=.myEffSizePropMat(prop_mat = prop_mat)
  n_x=n_x$effective_sample_size
  
  n_y=.myEffSizePropMat(prop_mat = prop_mat_c)
  n_y=n_y$effective_sample_size
  
  vxn=sweep(var_x,1,n_x,"/")
  vyn=sweep(var_y,1,n_y,"/")
  
  vxn2=sweep(var_x,1,n_x-1,"*")
  vyn2=sweep(var_y,1,n_y-1,"*")
  
  
  df=(vxn + vyn)^2
  df=df/(sweep(vxn^2,1,(n_x - 1),"/") + 
           sweep(vyn^2,1,n_y - 1,"/"))
  if(sum(vxn==0,na.rm = T)>0){
    df[Matrix::which(vxn==0)]=matrix(n_y,nrow=nrow(vxn),ncol=ncol(vxn),byrow = F)[Matrix::which(vxn==0)]
  }
  
  vxn[is.na(as.matrix(vxn))]=1
  sxy=sqrt(vxn+vyn)
  rm(vyn)
  
  d_sxy=sqrt(sweep(vxn2+vyn2,1,n_x+n_y-2,"/"))
  rm(vxn2,vyn2)
  
  sd_reg=sxy
  d_sd_reg=d_sxy
  if(F){
    
    if(reg_sd& extendedMode){
      
      if(F){
        #Vahid's version
        sd_reg=lapply(1:nrow(sxy),function(x){
          y=sxy[x,]
          tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=sxy[x,])
          slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
          if(length(slInd)>0){
            tmp_data=tmp_data[slInd,]
            meanExp=tmp_data$exp
            design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
            fit2=lm.fit(design,log10(tmp_data$sd))
            tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
            
            y[slInd]=tmp_data$fit2
          } else {
            y=pmax(y,1)
          }
          
          return(as.numeric(y))
        })
        sd_reg=do.call("rbind",sd_reg)
      }
      
      
      #Jonah's version
      sd_reg=apply(sxy,1, function(y){
        # browser()
        # y=sxy[x,]
        tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=y)
        slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
        if(length(slInd)>0){
          tmp_data=tmp_data[slInd,]
          meanExp=tmp_data$exp
          design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
          fit2=lm.fit(design,log10(tmp_data$sd))
          tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
          
          y[slInd]=tmp_data$fit2
        } else {
          y=pmax(y,1)
        }
        
        return(as.numeric(y))
      })
      ## browser()
      sd_reg=t(sd_reg)
      
      
      if(F){
        #Vahid's version
        d_sd_reg=lapply(1:nrow(d_sxy),function(x){
          y=d_sxy[x,]
          tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=d_sxy[x,])
          slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
          if(length(slInd)>0){
            tmp_data=tmp_data[slInd,]
            meanExp=tmp_data$exp
            design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
            fit2=lm.fit(design,log10(tmp_data$sd))
            tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
            
            y[slInd]=tmp_data$fit2
          } else {
            y=pmax(y,1)
          }
          
          return(as.numeric(y))
        })
        d_sd_reg=do.call("rbind",d_sd_reg)
      }
      
      
      if(extendedMode){
        d_sd_reg=apply(d_sxy, 1, function(y){
          # y=d_sxy[x,]
          tmp_data=data.frame(exp=resGeneMeanSd$geo_mean_count,sd=y)
          slInd=which(tmp_data$exp>0.001&tmp_data$sd>0)
          if(length(slInd)>0){
            tmp_data=tmp_data[slInd,]
            meanExp=tmp_data$exp
            design=splines::ns(log10(meanExp), df = 3, intercept = TRUE)
            fit2=lm.fit(design,log10(tmp_data$sd))
            tmp_data$fit2=pmax(tmp_data$sd,10^fit2$fitted.values)
            
            y[slInd]=tmp_data$fit2
          } else {
            y=pmax(y,1)
          }
          
          return(as.numeric(y))
        })
        d_sd_reg = t(d_sd_reg)
      }
      
      # d_sd_reg=do.call("rbind",d_sd_reg)
      
    }
    
  }
  
  
  rm(sxy,d_sxy)
  sd_offset=0
  
  y_exp=prop_mat_c %*% logNormData#t(logNormData)
  t <- (x_exp - y_exp)/(sd_reg+sd_offset)
  t[which(vxn==0)]=0
  if(F){
    t_adj=sweep(t,1,sqrt((n_x + n_y)/(n_x * n_y)),"*")
    t_adj=t_adj- sweep((3 * t_adj),1,(4 * ((n_x+n_y) - 2) - 1),"/")
    sigmad_adj=1/n_x+1/n_y+sweep((t_adj^2),1,(2*(n_x+n_y)),"/")
  }
  
  
  if(extendedMode){
    cohen_d=(x_exp - y_exp)/d_sd_reg
    hodge_g=cohen_d*(1-3/(4*(n_x+n_y)-9))
    if(sum(is.na(hodge_g))>0){
      hodge_g[is.na(hodge_g)]=0
    }
    
    se.g=0.5*sweep(hodge_g^2,1,n_x+n_y-3.94,"/")
    se.g <- sqrt(sweep(se.g,1, (n_x+n_y)/(n_x*n_y),"+"))
    if(sum(is.na(se.g))>0){
      se.g[is.na(se.g)]=0
    }
    
  }
  rm(vxn,sd_reg,y_exp)
  
  x_exp=as(x_exp, "dgCMatrix")
  gc()
  #t2 <- (x_exp - y_exp)/(sd_reg)
  
  t_vec=as.numeric(t)
  zscore=qnorm(pt(as.numeric(abs(t_vec)), as.numeric(df),lower.tail = F,log.p = T),lower.tail = F,log.p = T)
  zscore[is.na(zscore)]=0
  zscore=zscore*sign(t_vec)
  zscore=as(matrix(zscore,nrow=nrow(t),ncol=ncol(t)), "dgCMatrix")
  zscore_identity=zscore
  zscore_identity@x=rep(1,length(zscore_identity@x))
  rm(t_vec,df)
  
  colnames(zscore)=colnames(logNormData)
  row.names(zscore)=row.names(inputData$prop_mat)
  
  exp_binary=logNormData#t(logNormData)
  exp_binary@x=rep(1,length(exp_binary@x))
  
  pct.1=prop_mat %*% exp_binary
  pct.2=(prop_mat_c %*% exp_binary)*zscore_identity
  pct.2=Matrix::drop0(pct.2)
  
  rm(exp_binary)
  
  exp_norm=expm1(logNormData)#t(expm1(logNormData))
  
  logFC=log2(prop_mat %*% exp_norm+1) - log2(prop_mat_c %*% exp_norm+1)
  rm(exp_norm,logNormData)
  
  if(extendedMode){
    n=sqrt(sweep(pct.1,1,n_x,"*")*sweep(pct.2,1,n_y,"*"))
  }
  rm(n_y)
  
  if(sum(zscore>zscore_cap,na.rm = T)>0){
    zscore[which(zscore>zscore_cap)]=zscore_cap
  }
  
  if(sum(zscore<(-1*zscore_cap),na.rm = T)>0){
    zscore[which(zscore<(-1*zscore_cap))]=(-1*zscore_cap)
  }
  
  
  if(F){
    #moved to the propagation step
    matWeights=.myEffSizePropMat(prop_mat)
    
    matEffectiveSize=matWeights$effective_sample_size
    matWeights=matWeights$centroid_weights
    matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
    matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
    
  }
  
  
  matWeights=matrix(inputData$matWeights,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
  matEffectiveSize=matrix(inputData$matEffectiveSize,byrow = F,nrow=nrow(zscore),ncol=ncol(zscore))
  row.names(matWeights)=row.names(matEffectiveSize)=row.names(prop_mat)
  colnames(matWeights)=colnames(matEffectiveSize)=colnames(zscore)
  matWeights[which(pct.1*matEffectiveSize<minCellCountThr&pct.2<0.001)]=0
  matWeights[which(zscore==0|var_x==0)]=0
  
  
  matEffectiveSize[matWeights<0.0001]=0
  if(extendedMode){
    res=list(zscore=as(zscore, "dgCMatrix"),pct.1=as.big.matrix(pct.1),pct.2=as.big.matrix(pct.2),logFC=logFC,se.g=as.big.matrix(se.g),hodge_g=hodge_g,n=n,matWeights=matWeights,matEffectiveSize=matEffectiveSize,exp1=x_exp,exp2=x2_exp)
  } else {
    res=list(zscore=zscore,pct.1=pct.1,pct.2=pct.2,logFC=logFC,matWeights=matWeights,matEffectiveSize=matEffectiveSize,exp1=x_exp,gene_name_list=gene_name_list,pseudocell_name_list=pseudocell_name_list)
  }
  if(F){
    #Jonah's
    output = c(res,list(dsName=inputData$data$dsName,prop_mat=inputData$prop_mat,pseudocell_sim_mat=inputData$pseudocell_sim_mat))
    
    
    res_non_sparse = res
    for(dsitr in c("zscore", "pct.1", "pct.2", "logFC", "n", "matWeights", "matEffectiveSize")){
      res[[dsitr]] = as.sparse(res[[dsitr]])
    }
    gc()
    
    dir.create(file.path(argList$saveDir, "detail_exp_final"), showWarnings = FALSE)
    FN = file.path(argList$saveDir, "detail_exp_final", paste0(inputData$data$dsName, ".RDS"))
    
    .mcsaveRDS <- function(object,file,mc.cores=min(parallel::detectCores(),3, na.rm=T)) {
      con <- pipe(paste0("pigz -p",mc.cores," > ",file),"wb")
      saveRDS(object, file = con)
      close(con)
    }
    print(FN)
    
    .mcsaveRDS(output, FN)
    
    
  }
  
  gc()
  setdiff(names(inputData) ,c("matWeights","matEffectiveSize"))
  return(c(res,list(dsName=inputData$data$dsName)))
  #return(T)
}


.extra_sconline.argFn=function(runIndx,do.split.prop=T,min_cluster_size,pseudocell_selection_method,saveDir,exNonMicCells=F,prop.n.neighbors = 4,ncores=7,conservativeMapping=F,oldMapping=F,MGIsymbol=F,includeHuman,includeMouse,FinnishSbjBased,DE_supportingFractionThr,DE_n.adaptiveKernel,DE_nPropIter,Leng200=F,uniformZscore=F,dist_zscore_gamma,dist_zscore_norm,dist_zscore_nbinom,regularize,geoMean=F,prefix=NULL,newRun=F,inputDf=NULL,internal_pseudocell_count,pseudocell_size,sensitiveSearch,external_DE_path=NULL,external_DE_name=NULL,do.liger=NULL,singleton.method=singleton.method,include.singletons=T,colNormalize=T){
  
  if(is.null(inputDf)){
    if(includeMouse&(!includeHuman)){
      .dfSetting=.myRunSettingsComplete_mouse()
    } else if((!includeMouse)&(includeHuman)){
      .dfSetting=.myRunSettingsComplete_human()
    } else {
      .dfSetting=.myRunSettingsComplete_mouse()
    }
  } else {
    .dfSetting=inputDf
  }
  
  
  .commonExpressed=.dfSetting$commonExpressed[runIndx]
  .includeHuman=includeHuman
  .includeMouse=includeMouse
  .nPCs=.dfSetting$nPCs[runIndx]
  .HVG_count=.dfSetting$HVG_count[runIndx] #4
  .HVG_method="vst" #"vst" and "sctransform"
  .exNonMicCells=exNonMicCells
  .covariates=.dfSetting$covariates[runIndx]
  .removeHighExp=.dfSetting$removeHighExp[runIndx]
  .UMI_cor_thr=.dfSetting$UMI_cor_thr[runIndx]
  .slMicrogliaClusters=.dfSetting$slMicrogliaClusters[runIndx]
  .breakHammond=T
  if(sum(colnames(.dfSetting)=="breakHammond")>0){
    .dfSetting$breakHammond[runIndx]
  }
  if(.slMicrogliaClusters==""){
    .slMicrogliaClusters=NULL
  }
  
  if(!is.null(.slMicrogliaClusters)){
    .exNonMicCells=T
    .slMicrogliaClusters=unlist(strsplit(.slMicrogliaClusters,","))
  }
  
  if(!is.null(.covariates)){
    if(all(.covariates=='3')){
      .covariates=c('QC_Gene_total_count','QC_top50_pct','QC_Gene_unique_count')
    } else if(all(.covariates=='2')) {
      .covariates=c('QC_top50_pct','QC_Gene_unique_count')
    } else if(all(.covariates=='1')){
      .covariates=c('QC_Gene_unique_count')
    } else if(all(.covariates=="")){
      .covariates=NULL
    }
  }
  
  
  .exCellStates=c("CAM","Mnc","BAM","T cells","Neutrophils","cDC","B cells","NK cells","Proliferation","T/NKT cells","cDC2","Monocyte/Mdc","ILC","ydT cells","pDC","cDC1","migDC","MCs","Non-cl. monocytes")
  
  #"depth_per_gene"#'QC_Gene_unique_count'#c('QC_Gene_total_count','QC_top50_pct',)#, "organism","depth_per_gene" 'QC_Gene_total_count','QC_top50_pct','QC_Gene_unique_count'
  .indScaling=.dfSetting$indScaling[runIndx]
  
  if(is.null(prefix)){
    prefix=""
  }
  
  if(sum(c(.includeHuman,.includeMouse))==2){
    prefix=paste0(prefix,"human-mouse")
  } else if(.includeMouse){
    prefix=paste0(prefix,"mouse")
  } else if(.includeHuman){
    prefix=paste0(prefix,"human")
  }
  
  if(is.null(do.liger)){
    do.liger=.dfSetting$do.liger[runIndx]
  }
  
  if(do.liger){
    prefix=paste0(prefix,"_liger")
  }
  
  if(!is.null(external_DE_name)){
    prefix=paste0(prefix,"_",external_DE_name)
  }
  
  if(Leng200){
    prefix=paste0(prefix,"_Leng200")
  }
  
  .saveDir=.mySaveDirMaker(paste0(saveDir,"/",prefix),nPCs = .nPCs,cov=.covariates,exNonMicCells=.exNonMicCells,FinnishSbjBased)
  if(.exNonMicCells){
    .saveDirGlobal=paste0(saveDir,"/",prefix,"-Global-exNonMicCells")
  } else {
    .saveDirGlobal=paste0(saveDir,"/",prefix,"-Global")
  }
  
  if(FinnishSbjBased&includeHuman){
    .saveDirGlobal=paste0(.saveDirGlobal,"-FinnishSbjBased")
  }
  
  argList=list(commonExpressed=.dfSetting$commonExpressed[runIndx],min_cluster_size=min_cluster_size,pseudocell_selection_method=pseudocell_selection_method,do.split.prop=do.split.prop,prop.n.neighbors = prop.n.neighbors,includeHuman=.includeHuman,includeMouse=.includeMouse,conservativeMapping=conservativeMapping,oldMapping=oldMapping,MGIsymbol=MGIsymbol,nPCs=.dfSetting$nPCs[runIndx],HVG_count=.dfSetting$HVG_count[runIndx],HVG_method="vst",include.singletons=include.singletons,colNormalize=colNormalize,exNonMicCells=exNonMicCells,covariates=.covariates,indScaling=.dfSetting$indScaling[runIndx],saveDir=.saveDir,saveDirGlobal=.saveDirGlobal,ncores=ncores,excludeHighExp=.removeHighExp,slMicrogliaClusters=.slMicrogliaClusters,breakHammond=.breakHammond,UMI_cor_thr=.UMI_cor_thr,onlyRankBased=.dfSetting$onlyRankBased[runIndx],varScore.thr=.dfSetting$varScore.thr[runIndx],do.liger=.dfSetting$do.liger[runIndx],allGenesFraction=.dfSetting$allGenesFraction[runIndx],FinnishSbjBased=FinnishSbjBased,uniformZscore=uniformZscore,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=DE_nPropIter,dist_zscore_gamma=dist_zscore_gamma,dist_zscore_norm=dist_zscore_norm,dist_zscore_nbinom=dist_zscore_nbinom,regularize=regularize,geoMean=geoMean,prefix=prefix,internal_pseudocell_count=internal_pseudocell_count,pseudocell_size=pseudocell_size,sensitiveSearch=sensitiveSearch,external_DE_path=external_DE_path,external_DE_name=external_DE_name,singleton.method=singleton.method)
  
  .extraCreateDirFn(argList)
  
  .extraNewRun(argList=argList,newRun=newRun)
  
  return(argList)
  
}

#inputData=res;argList=.ArgList;min_effective_size=5;pd=pd;cell_annoCol=annoCol
.extra_sconline.visPseudocellAnno=function(inputData,argList,min_effective_size=5,pd=NULL,cell_annoCol=NULL,pie_scale=1){
  require(ggplot2)
  require(scatterpie)
  require(hues)
  
  if(F){
    if(sum(colnames(inputData)=="dataset")>0){
      inputData=inputData[,-which(colnames(inputData)=="dataset")]
      absent_pseudocells=apply(inputData[,-which(colnames(inputData) %in% c("pseudocell","effective_size"))],1,sum)
      absent_pseudocells=which(absent_pseudocells<0.2)
      if(length(absent_pseudocells)>0){
        inputData[absent_pseudocells,-which(colnames(inputData)=="pseudocell")]=NA
      }
      
      inputData=aggregate(.~pseudocell,data=inputData,function(x) mean(x,na.rm=T))
    }
  }
  
  
  UMAP_centroid=.sconline.fetch_data("umap_centroids",argList=argList)
  
  if(is.null(pd)){
    pd=.sconline.fetch_data("annotation",argList=argList)
  }
  
  if(is.null(cell_annoCol)){
    pd_summary=.extra_sconline.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"),n=300)
  } else {
    pd_summary=pd
  }
  
  
  if(sum(colnames(inputData)=="dataset")>0){
    inputData=inputData[which(inputData$effective_size>min_effective_size),]
    
  }
  piechart_data=merge(inputData,UMAP_centroid,by.x="pseudocell",by.y="centroid",all.x=T)
  
  anno_cols=setdiff(colnames(inputData),c("pseudocell","effective_size","dataset"))
  do.scatterPlot=F
  if(length(anno_cols)>1){
    do.scatterPlot=T
  }
  
  if(is.null(cell_annoCol)){
    if(do.scatterPlot){
      pd$cluster_anno_res=factor(as.character(pd$cluster_anno_res),levels = anno_cols)
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(anno_cols))))+facet_wrap(~dataset)
        
        
      } else {
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(anno_cols))))
        
      }
    } else {
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_point(aes(x=UMAP_1, y=UMAP_2, color=cluster), data=piechart_data) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(unique(pd$cluster_anno_res)))))+facet_wrap(~dataset)+scale_color_manual(values=c("gray",hues::iwanthue(length(unique(pd$cluster_anno_res)))))
        
        
      } else {
        
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2),color="lightgray",size=0.3) + geom_point(aes(x=UMAP_1, y=UMAP_2, fill=cluster),color="black", data=piechart_data) + coord_equal()+theme_classic()+scale_fill_manual(values=c("gray",hues::iwanthue(length(unique(pd$cluster_anno_res)))))+scale_color_manual(values=c("gray",hues::iwanthue(length(unique(pd$cluster_anno_res)))))
        
      }
    }
    
  } else {
    pd$cluster_anno_res=pd[,cell_annoCol]
    if(do.scatterPlot){
      pd$cluster_anno_res=factor(as.character(pd$cluster_anno_res),levels = anno_cols)
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(length(anno_cols))))+scale_color_manual(values=c(hues::iwanthue(length(anno_cols))))+facet_wrap(~dataset)
        
        
      } else {
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_scatterpie(aes(x=UMAP_1, y=UMAP_2, group=pseudocell), data=piechart_data,
                                                                                                                  cols=anno_cols,pie_scale = pie_scale) + coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(length(anno_cols))))+scale_color_manual(values=c(hues::iwanthue(length(anno_cols))))
        
      }
    } else {
      if(sum(colnames(inputData)=="dataset")>0){
        p <- ggplot()+geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=cluster_anno_res),size=0.3) + geom_point(aes(x=UMAP_1, y=UMAP_2, fill=cluster),color="black",shape=21)+
                                                                                                             coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))+facet_wrap(~dataset)+scale_color_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))
        
        
      } else {
        p=.myDimPlotFn(object=pd_summary, dimCols = c("UMAP_1", "UMAP_2"), attCol = "cluster_anno_res",set_col=F) + geom_point(data=piechart_data,aes(x=UMAP_1, y=UMAP_2, fill=cluster),color="black",shape=21)+ 
                                                                                                             coord_equal()+theme_classic()+scale_fill_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))+scale_color_manual(values=c(hues::iwanthue(length(unique(pd$cluster_anno_res)))))
        
      }
    }
  }
  
  
  return(p)
  
}

.extra_sconline.CosineDistFn=function(inputMatrix,sig_thr=3){
  
  mymakeCosine=function(inputMatrix,wJaccardDistance=F,sig_only=T,sig_thr=3){
    
    res_mat=NULL
    if(wJaccardDistance){
      ##########################
      res_mat1=res_mat2=inputMatrix
      res_mat1[which(inputMatrix<3)]=0
      res_mat1[which(res_mat1>0)]=1
      
      res_mat2[which(inputMatrix<2)]=0
      res_mat2[which(res_mat2>0)]=1
      
      res_mat=res_mat1 %*% t(res_mat2)
      res_mat=sweep(res_mat,1,pmax(rowSums(res_mat1),1),"/")
      row.names(res_mat)=colnames(res_mat)=row.names(inputMatrix)
      
      res_mat2=t(res_mat)
      res_mat[res_mat<res_mat2]=res_mat2[res_mat<res_mat2]
      rm(res_mat1,res_mat2)
      
      ##########################
    }
    
    if(sig_only){
      inputMatrix=mymakeCosine2(inputMatrix = inputMatrix,sig_thr = sig_thr)
    } else {
      inputMatrix_cosine=inputMatrix*inputMatrix
      inputMatrix_cosine=sqrt(rowSums(inputMatrix_cosine))
      inputMatrix=sweep(inputMatrix,1,inputMatrix_cosine,"/")
      inputMatrix=inputMatrix
      inputMatrix=inputMatrix %*% t(inputMatrix)
    }
    
    
    if(wJaccardDistance){
      inputMatrix=res_mat
    }
    
    return(inputMatrix)
  }
  
  
  res_mat=NULL
  
  inputMatrix_org=inputMatrix
  inputMatrix_cosine=inputMatrix*inputMatrix
  inputMatrix_cosine=sqrt(rowSums(inputMatrix_cosine))
  inputMatrix=sweep(inputMatrix,1,inputMatrix_cosine,"/")
  
  
  inputMatrix2=matrix(0,nrow=nrow(inputMatrix),ncol=nrow(inputMatrix))
  
  for(i in 1:nrow(inputMatrix)){
    #print(i)
    for(j in i:nrow(inputMatrix)){
      sl_ind=apply(inputMatrix_org[c(i,j),],2,max)
      sl_ind=which(sl_ind>=sig_thr)
      if(length(sl_ind)>10){
        sl_ind=mymakeCosine(inputMatrix_org[c(i,j),sl_ind],wJaccardDistance = F,sig_only = F,sig_thr = sig_thr)
        inputMatrix2[i,j]=sl_ind[1,2]
      } else {
        inputMatrix2[i,j]=NA
      }
      
    }
  }
  inputMatrix2=inputMatrix2+t(inputMatrix2)
  diag(inputMatrix2)=1
  
  row.names(inputMatrix2)=colnames(inputMatrix2)=row.names(inputMatrix)
  
  return(inputMatrix2)
}

.extra_sconline.PctScoreFn_slow=function(argList,pct_mat,meta_z_mat,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F){
  
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  meta_z_mat=apply(meta_z_mat,2,function(x) as.numeric(x>sig1_thr))
  
  pctDiffCountFn_org=function(i, pct_mat,pct_mat_ref,centers,pct_diff_thr,pct2_thr,meta_z_mat,sig1_thr){
    res_counts=unlist(lapply(1:length(centers),function(j){
      sum((.extra_sconline.PctDiffscoreFn(pct_mat[i,]-pct_mat_ref[jList[j],],pct_diff_thr = pct_diff_thr)*.extra_sconline.Pct2scoreFn(pct_mat_ref[jList[j],],pct2_thr = pct2_thr)>0.99)&meta_z_mat[i,]>sig1_thr)
    }))
    return(res_counts)
  }
  
  
  pctDiffCountFn=function(ilist, pct_mat2,pct_diff_thr,pct2_thr,meta_z_mat2){
    res=matrix(0,nrow=length(ilist),ncol=nrow(pct_mat2))
    counter=0
    for(i in ilist){
      counter=counter+1
      p1=sweep(pct_mat2,2,pct_mat2[i,],"-")*(-1)
      #p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=unlist(lapply(1:nrow(p1),function(x) {
        x1=.extra_sconline.PctDiffscoreFn(p1[x,],pct_diff_thr = pct_diff_thr)
        x2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2[x,],pct2_thr = pct2_thr)
        sum(as.numeric(x1*x2>0.99)*meta_z_mat2[i,])
      }))
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  tmp_groups=split(1:nrow(pct_mat),ceiling(1:nrow(pct_mat)/(nrow(pct_mat)/(argList$ncores*3))))
  
  tstFn=function(tmp_groups){
    res_mat=parallel::mclapply(tmp_groups,pctDiffCountFn,pct_mat2=pct_mat2,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat2,mc.cores = argList$ncores)
    res_mat=do.call("rbind",res_mat)
    return(res_mat)
  }
  
  #pct_mat2=pct_mat;meta_z_mat2=meta_z_mat
  my.env <- new.env()
  my.env$meta_z_mat2 <- meta_z_mat
  my.env$pct_mat2 <- pct_mat
  my.env$pct_diff_thr <- pct_diff_thr
  my.env$pct2_thr <- pct2_thr
  my.env$argList <- argList
  
  with_env <- function(f, e=parent.frame()) {
    stopifnot(is.function(f))
    environment(f) <- e
    f
  }
  
  res_mat=with_env(tstFn,my.env)(tmp_groups)
  
  #res_mat=lapply(1:nrow(prop_mat),pctDiffCountFn_org,centers=centers,pct_mat_ref=pct_mat,pct_mat=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat=meta_z,sig1_thr=sig1_thr)
  #tst=pctDiffCountFn(ilist=tmp_groups[[1]],pct_mat2=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat)
  #res_mat=do.call("rbind",res_mat)
  gc()
  if(symmetric){
    res_mat2=res_mat+t(res_mat)
  } else {
    res_mat2=res_mat
  }
  
  row.names(res_mat2)=row.names(pct_mat)
  colnames(res_mat2)=row.names(pct_mat)
  
  return(res_mat2)
}

.extra_sconline.PctScoreFn=function(argList,pct_mat,meta_z_mat,sig1_thr=3,centers=NULL,pct2_thr=0.3,pct_diff_thr=0.2,symmetric=F){
  
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  meta_z_mat=apply(meta_z_mat,2,function(x) as.numeric(x>sig1_thr))
  
  p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat,pct2_thr = pct2_thr)
  
  res_mat=matrix(0,nrow=nrow(pct_mat),ncol=length(centers))
  pctDiffCountFn_org=function(i, pct_mat,pct_mat_ref,centers,pct_diff_thr,pct2_thr,meta_z_mat,sig1_thr){
    res_counts=unlist(lapply(1:length(centers),function(j){
      sum((.extra_sconline.PctDiffscoreFn(pct_mat[i,]-pct_mat_ref[jList[j],],pct_diff_thr = pct_diff_thr)*.extra_sconline.Pct2scoreFn(pct_mat_ref[jList[j],],pct2_thr = pct2_thr)>0.99)&meta_z_mat[i,]>sig1_thr)
    }))
    return(res_counts)
  }
  
  
  pctDiffCountFn=function(ilist, pct_mat, pct_mat2,pct_diff_thr,meta_z_mat){
    res=matrix(0,nrow=length(ilist),ncol=ncol(pct_mat))
    counter=0
    
    for(i in ilist){
      counter=counter+1
      #print(counter)
      p1=Matrix::drop0(pct_mat[,i] - pct_mat)
      #p1=sweep(pct_mat,2,pct_mat[i,],"-")*(-1)
      p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=colSums((p1*pct_mat2>0.99)*meta_z_mat[,i])
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  tmp_groups=split(1:nrow(pct_mat),ceiling(1:nrow(pct_mat)/(nrow(pct_mat)/(argList$ncores*3))))
  
  
  tstFn=function(tmp_groups){
    res_mat=parallel::mclapply(tmp_groups,pctDiffCountFn,pct_mat=pct_mat,pct_mat2=pct_mat2,pct_diff_thr=pct_diff_thr,meta_z_mat=meta_z_mat,mc.cores = argList$ncores)
    res_mat=do.call("rbind",res_mat)
    return(res_mat)
  }
  
  #pct_mat2=pct_mat;meta_z_mat2=meta_z_mat
  my.env <- new.env()
  my.env$meta_z_mat <- t(meta_z_mat)
  my.env$pct_mat2 <- t(p2)
  my.env$pct_mat <- t(pct_mat)
  my.env$pct_diff_thr <- pct_diff_thr
  my.env$argList <- argList
  
  with_env <- function(f, e=parent.frame()) {
    stopifnot(is.function(f))
    environment(f) <- e
    f
  }
  
  res_mat=with_env(tstFn,my.env)(tmp_groups)
  
  #res_mat=lapply(1:nrow(prop_mat),pctDiffCountFn_org,centers=centers,pct_mat_ref=pct_mat,pct_mat=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat=meta_z,sig1_thr=sig1_thr)
  #tst=pctDiffCountFn(ilist=tmp_groups[[1]],pct_mat2=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat)
  #res_mat=do.call("rbind",res_mat)
  gc()
  if(symmetric){
    res_mat2=res_mat+t(res_mat)
  } else {
    res_mat2=res_mat
  }
  
  row.names(res_mat2)=row.names(pct_mat)
  colnames(res_mat2)=row.names(pct_mat)
  
  return(res_mat2)
}

.extra_sconline.PctScoreFn2=function(){
  
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  meta_z_mat=apply(meta_z_mat,2,function(x) as.numeric(x>sig1_thr))
  
  res_mat=matrix(0,nrow=nrow(pct_mat),ncol=length(centers))
  pctDiffCountFn_org=function(i, pct_mat,pct_mat_ref,centers,pct_diff_thr,pct2_thr,meta_z_mat,sig1_thr){
    res_counts=unlist(lapply(1:length(centers),function(j){
      sum((.extra_sconline.PctDiffscoreFn(pct_mat[i,]-pct_mat_ref[jList[j],],pct_diff_thr = pct_diff_thr)*.extra_sconline.Pct2scoreFn(pct_mat_ref[jList[j],],pct2_thr = pct2_thr)>0.99)&meta_z_mat[i,]>sig1_thr)
    }))
    return(res_counts)
  }
  
  
  pctDiffCountFn=function(ilist, pct_mat2,pct_diff_thr,pct2_thr,meta_z_mat2){
    res=matrix(0,nrow=length(ilist),ncol=nrow(pct_mat2))
    counter=0
    for(i in ilist){
      counter=counter+1
      p1=sweep(pct_mat2,2,pct_mat2[i,],"-")*(-1)
      #p1=.extra_sconline.PctDiffscoreFn(p1,pct_diff_thr = pct_diff_thr)
      #p2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2,pct2_thr = pct2_thr)
      p1=unlist(lapply(1:nrow(p1),function(x) {
        x1=.extra_sconline.PctDiffscoreFn(p1[x,],pct_diff_thr = pct_diff_thr)
        x2=.extra_sconline.Pct2scoreFn(pct2_value = pct_mat2[x,],pct2_thr = pct2_thr)
        sum(as.numeric(x1*x2>0.99)*meta_z_mat2[i,])
      }))
      #p1=do.call("rbind",p1)
      #p1=Matrix::rowSums(p1)
      res[counter,]=p1
    }
    row.names(res)=ilist
    
    return(res)
  }
  
  tmp_groups=split(1:nrow(pct_mat),ceiling(1:nrow(pct_mat)/(nrow(pct_mat)/(argList$ncores*3))))
  
  tstFn=function(tmp_groups){
    res_mat=parallel::mclapply(tmp_groups,pctDiffCountFn,pct_mat2=pct_mat2,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat2,mc.cores = argList$ncores)
    res_mat=do.call("rbind",res_mat)
    return(res_mat)
  }
  
  #pct_mat2=pct_mat;meta_z_mat2=meta_z_mat
  my.env <- new.env()
  my.env$meta_z_mat2 <- meta_z_mat
  my.env$pct_mat2 <- pct_mat
  my.env$pct_diff_thr <- pct_diff_thr
  my.env$pct2_thr <- pct2_thr
  my.env$argList <- argList
  
  with_env <- function(f, e=parent.frame()) {
    stopifnot(is.function(f))
    environment(f) <- e
    f
  }
  
  res_mat=with_env(tstFn,my.env)(tmp_groups)
  
  #res_mat=lapply(1:nrow(prop_mat),pctDiffCountFn_org,centers=centers,pct_mat_ref=pct_mat,pct_mat=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat=meta_z,sig1_thr=sig1_thr)
  #tst=pctDiffCountFn(ilist=tmp_groups[[1]],pct_mat2=pct_mat,pct_diff_thr=pct_diff_thr,pct2_thr=pct2_thr,meta_z_mat2=meta_z_mat)
  #res_mat=do.call("rbind",res_mat)
  gc()
  if(symmetric){
    res_mat2=res_mat+t(res_mat)
  } else {
    res_mat2=res_mat
  }
  
  row.names(res_mat2)=row.names(pct_mat)
  colnames(res_mat2)=row.names(pct_mat)
  
  return(res_mat2)
}


.extra_sconline.PctDiffscoreFn=function(pct_diff_value,pct_diff_thr){
  #pct_diff_value=1-pct_diff_value
  #pct_diff_thr=1-pct_diff_thr
  pct_diff_value=pct_diff_thr - pct_diff_value# - pct_diff_thr
  pct_diff_value@x=(-1)*pmax(pct_diff_value@x,0)*20 - pmin(pct_diff_value@x,0)*2
  pct_diff_value@x=1-exp(pct_diff_value@x)
  pct_diff_value=1-as.matrix(pct_diff_value)
  #pct_diff=exp(-pmax(pct_diff_value,0)*20-pmin(pct_diff_value,0)*2)
  return(pct_diff_value)
}

.extra_sconline.Pct2scoreFn=function(pct2_value,pct2_thr){
  pct2=round(exp(-(pmax((pct2_value-pct2_thr),0)/(0.1))^2),3)
  return(pct2)
}


######################
#main sconline functions

#Performs UMAP on the PC space and maps pseudocells to it
.sconline.umapFn_org=function(argList,umap.method='umap-learn'){
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if(!reRunCheck){
    reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_centroid",argList=argList));F}, error=function(e) {return(T)})
  }
  
  if(reRunCheck|argList$newRun){
    cat("Running UMAP\n")
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    
    load(.myFilePathMakerFn("pca_centroids",argList=argList))
    
    pca_centroid=pca_centroid[,1:argList$nPCs]
    
    
    
    tst=.reductionUMAPFn(harmony_embeddings,umap.method=umap.method,testPCAembeddings=pca_centroid)
    resUMAP=tst$embedding
    .mySaveFn(resUMAP,file=.myFilePathMakerFn("UMAP_res",argList=argList))
    
    x=as.data.frame(resUMAP)
    x=cbind(x,pd)
    row.names(x)=row.names(pd)
    pd=x
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    tst=tst[["test"]]
    UMAP_centroid=data.frame(centroid=row.names(pca_centroid),UMAP_1=tst[,1],UMAP_2=tst[,2],stringsAsFactors = F)
    .mySaveFn(UMAP_centroid,file=.myFilePathMakerFn("UMAP_centroid",argList=argList))
    
    
  }
  return("Done!")
}

#argList=.ArgList;inputExpData=data;inputGenes="DDR3"
.sconline.markerPlot=function(argList,inputExpData,inputGenes){
  library(dplyr)
  library(purrr)
  library(cowplot)
  library(patchwork)
  
  load(.myFilePathMakerFn("consensusDE_markers",argList=argList,uniformImportant=T))
  load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
  
  cat("           Generating marker figs (cell based) ...\n")
  load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  load(do.call('.myFilePathMakerFn',args=c(list(filename="exp_merged",expData=T),argList)))
  
  tmp=tmp[,colnames(tmp) %in% row.names(pd)]
  pd2=pd[match(colnames(tmp),row.names(pd)),]
  
  
  plotFn_integrated=function(index,geneNameList,expData,annotationFile,input_cell_bkg_umap,prefix){
    
    genes=unlist(geneNameList[[index]]$gene_short_name)
    
    figWidth=49
    if(length(genes)<21){
      figWidth=49/3*ceiling(length(genes)/7)
    }
    
    figHieght=35
    if(length(genes)<7){
      figHieght=35/7*length(genes)
    }
    
    expData=expData[match(c(geneNameList[[index]]$gene,setdiff(row.names(expData),geneNameList[[index]]$gene)),row.names(expData)),]
    p=list()
    for(i in seq(1,length(genes),7)){
      p1=.myFeaturePlot(inputSeurat = expData,inputDimData = annotationFile,inputGenes = genes[i:min((i+6),length(genes))],combine_figs=F)
      p2=.my2dPlot_counts(inputPCA=annotationFile,batch_values="batch_merging",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=genes[i:min((i+6),length(genes))],geneNameCol="gene_short_name",expData=expData,combine_figs=F)
      
      tmp_UMAP=list()
      for(ij in i:min((i+6),length(genes))){
        tmp=UMAP_centroid
        tmp$gene_id=geneNameList[[index]]$gene[which(geneNameList[[index]]$gene_short_name==genes[ij])]
        tmp$gene=genes[ij]
        tmp_res=res_arranged[which(res_arranged$gene==unique(tmp$gene_id)),]
        tmp_res=tmp_res[match(as.character(tmp$centroid),as.character(tmp_res$centroid)),]
        tmp$zscore=tmp_res$zscore
        tmp$effectiveSize=tmp_res$effective_size
        tmp$size=tmp_res$count
        if(sum(is.na(tmp$zscore))>0){
          tmp$zscore[is.na(tmp$zscore)]=0
          tmp$effectiveSize[is.na(tmp$effectiveSize)]=0
          tmp$size[is.na(tmp$size)]=0
        }
        tmp$effectiveSize=tmp$effectiveSize+1
        tmp_UMAP=c(tmp_UMAP,list(tmp))
      }
      tmp_UMAP=do.call("rbind",tmp_UMAP)
      
      tmp_UMAP$gene=factor(as.character(tmp_UMAP$gene),levels=genes[i:min((i+6),length(genes))])
      
      tmp_UMAP$effectiveSize=as.numeric(as.character(tmp_UMAP$effectiveSize))
      tmp_UMAP$effectiveSize[is.na(tmp_UMAP$effectiveSize)]=min(min(tmp_UMAP$effectiveSize,na.rm = T),1)
      if(sum(tmp_UMAP$zscore>5)>0){
        tmp_UMAP$zscore[which(tmp_UMAP$zscore>5)]=5
      }
      if(sum(tmp_UMAP$zscore<(-5))>0){
        tmp_UMAP$zscore[which(tmp_UMAP$zscore<(-5))]=(-5)
      }
      p3=tmp_UMAP %>% 
        group_split(gene) %>% 
        map(
          ~ggplot(.,aes(UMAP_1,UMAP_2,color=zscore,size=effectiveSize))+geom_point(data=input_cell_bkg_umap,aes(UMAP_1,UMAP_2),color="lightgrey",size=0.1)+geom_point()+theme_cowplot()+ theme(plot.title = element_text(hjust = 0.5))+ggtitle(toupper(unique(.$gene)))+scale_color_gradientn(colors = c("lightblue","lightblue","black","yellow","orange","red"),breaks=c(-5,-2,0,2,3,5),limits=c(-5,5))+scale_size_continuous(limits=c(0,max(tmp_UMAP$effectiveSize)),breaks=seq(0,max(tmp_UMAP$effectiveSize),5))
        ) #%>% 
      #plot_grid(plotlist = ., align = 'hv',ncol=1)
      
      if(length(genes)<(i+6)&length(genes)>7){
        dfNull=data.frame()
        pNull=ggplot(dfNull)+geom_point()
        for(inull in 1:(i+6-length(genes))){
          p1=c(p1,list(pNull))
          p2=c(p2,list(pNull))
          p3=c(p3,list(pNull))
        }
      }
      
      p=c(p,p1,p3,p2)
      
    }
    p=wrap_plots(p,nrow = min(length(genes),7),byrow = F)
    if(is.null(prefix)){
      ggsave(plot=p,file=.myFilePathMakerFn(filename = paste0("marker_",names(geneNameList)[index]),argList = argList,pdf = T,uniformImportant=T,makePlotDir=T),width = figWidth,height = figHieght)
    } else {
      ggsave(plot=p,file=.myFilePathMakerFn(filename = paste0(prefix,"_marker_",names(geneNameList)[index]),argList = argList,pdf = T,uniformImportant=T,makePlotDir=T),width = figWidth,height = figHieght)
    }
    
    return("Done")
  }
  
  if(!is.null(inputGenes)){
    geneEffectSizes=geneEffectSizes[which(toupper(geneEffectSizes$gene) %in% toupper(inputGenes)),]
  }
  slGenes=geneEffectSizes[order(geneEffectSizes$effectSize,decreasing = T),]
  id=rep(0,nrow(slGenes))
  for(ic in unique(slGenes$cluster)){
    tmpId=which(slGenes$cluster==ic)
    tmpId=split(tmpId,ceiling(seq_along(tmpId)/21))
    counter=0
    for(ij in 1:length(tmpId)){
      counter=counter+1
      slGenes$cluster[tmpId[[ij]]]=paste0("cluster",slGenes$cluster[tmpId[[ij]]],"_",counter)
    }
  }
  slGenes=split(slGenes, slGenes$cluster)
  
  input_cell_bkg_umap=.scatterPlot_summary2d(object=pd2,reductionCols=c("UMAP_1","UMAP_2"))
  
  plotMI=parallel::mclapply(1:length(slGenes),plotFn_integrated,geneNameList=slGenes,expData=tmp[unique(c(1:100,which(row.names(tmp) %in% geneEffectSizes$gene))),],annotationFile=pd2,input_cell_bkg_umap=input_cell_bkg_umap,prefix=prefix,mc.cores = argList$ncores)
  save(plotMI,file="final_errors.rda")
  
  
  return("Done!")
}



.sconline.umapFn=function(argList,umap.method='umap-learn',generateUMAP=T,saveFiles=T,input_UMAP_embedding=NULL){
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F));F}, error=function(e) {return(T)})
  
  if((reRunCheck|argList$newRun)&generateUMAP&is.null(input_UMAP_embedding)){
    cat("Running UMAP\n")
    load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    
    #load(.myFilePathMakerFn("pca_centroids",argList=argList))
    #pca_centroid=pca_centroid[,1:argList$nPCs]
    
    
    
    tst=.reductionUMAPFn(harmony_embeddings,umap.method=umap.method,testPCAembeddings=NULL,n.neighbors =40)
    resUMAP=tst$embedding
    .mySaveFn(resUMAP,file=.myFilePathMakerFn("UMAP_res",argList=argList,pseudoImportant = F))
    
    x=as.data.frame(resUMAP)
    x=cbind(x,pd)
    row.names(x)=row.names(pd)
    pd=x
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
  } else if(!is.null(input_UMAP_embedding)){
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    if(sum(row.names(pd) %in% row.names(input_UMAP_embedding))<nrow(pd)){
      stop("some of the cells don't have UMAP embedding. row.names of the input_UMAP_embedding should contain all cell names")
    }
    input_UMAP_embedding=input_UMAP_embedding[match(row.names(pd),row.names(input_UMAP_embedding)),]
    pd$UMAP_1=input_UMAP_embedding[,"UMAP_1"]
    pd$UMAP_2=input_UMAP_embedding[,"UMAP_2"]
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else if(!generateUMAP){
    load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    pd$UMAP_1=1
    pd$UMAP_2=2
    .mySaveFn(pd,file=.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  }
  
  reRunCheck=tryCatch({load(.myFilePathMakerFn("UMAP_centroid",argList=argList));F}, error=function(e) {return(T)})
  if((reRunCheck|argList$newRun)|(!is.null(input_UMAP_embedding))){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    UMAP_centroid=data.frame(centroid=sl_pseudo$cluster,pd[sl_pseudo$pseudocell,c("UMAP_1","UMAP_2")],stringsAsFactors = F)
    row.names(UMAP_centroid)=UMAP_centroid$centroid
    .mySaveFn(UMAP_centroid,file=.myFilePathMakerFn("UMAP_centroid",argList=argList))
  }
  
  
  return("Done!")
}

.extra_sconline.JaccardDistanceFn=function(meta_z_mat,meta_z_mat_bkg=NULL,sig1_thr=3,sig2_thr=1,mode="ratio"){
  res_mat1=meta_z_mat
  if(is.null(meta_z_mat_bkg)){
    meta_z_mat_bkg=meta_z_mat
  }
  res_mat2=meta_z_mat_bkg
  res_mat1[which(meta_z_mat<sig1_thr)]=0
  res_mat1[which(res_mat1>0)]=1
  
  res_mat2[which(meta_z_mat_bkg<sig2_thr)]=0
  res_mat2[which(res_mat2>0)]=1
  
  res_mat=res_mat1 %*% t(res_mat2)
  if(mode=="ratio"){
    res_mat=sweep(res_mat,1,pmax(rowSums(res_mat1),1),"/")
  } else if(mode=="difference"){
    res_mat=sweep(res_mat,1,pmax(rowSums(res_mat1),1),"-")*(-1)
  }
  
  row.names(res_mat)=colnames(res_mat)=row.names(meta_z_mat)
  return(res_mat)
}

.extra_sconline.densityPeakClustering=function(meta_z_mat,cosine_dist,de_dist,pct_de_count_thr=0,sig1_thr,min_marker_thr,pseudocells=NULL){
  
  if(!is.null(pseudocells)){
    cosine_dist=cosine_dist[pseudocells,pseudocells]
    de_dist=de_dist[pseudocells,pseudocells]
    if(length(pseudocells)<3){
      return(list(D=pseudocells[1],cluster_assignments=data.frame(local_density=0,centroid=pseudocells,sigma=0,source_name=pseudocells[1])))
    }
  }
  
  local_density=cosine_dist#de_dist
  local_density[is.na(local_density)]=0
  local_density[t(de_dist)>pct_de_count_thr]=0
  local_density=.extra_sconline.AffinityFn(sim_mat = local_density)
  local_density_thr=apply(local_density,1,function(x) quantile(x,0.95))
  local_density2=apply(local_density,1,function(x) sum(x[x>=quantile(local_density_thr,0.5)],na.rm = T))
  local_density=apply(local_density,1,function(x) sum(x>=quantile(local_density_thr,0.75)))
  
  local_density=local_density[order(local_density,local_density2,decreasing = T)]
  #ggplot(data=data$input_umap_centroid,aes(UMAP_1,UMAP_2,label=centroid))+geom_point()+geom_label(data=data$input_umap_centroid[data$input_umap_centroid$centroid %in% gsub("C","",names(local_density[1:40])),])
  
  sig_thr=apply(meta_z_mat,1,function(x) sum(x>=sig1_thr))
  local_density=local_density[names(local_density) %in% row.names(meta_z_mat)[sig_thr>=min_marker_thr]]
  
  de_dist_sym=de_dist+t(de_dist)
  
  res_min_dist=0
  res_source_name=names(local_density)[1]
  for(i in 2:length(local_density)){
    source_name_list=names(local_density)[1:(i-1)]
    tmp_score=de_dist_sym[names(local_density[i]),source_name_list]
    sl_id=which(tmp_score==min(tmp_score))[1]
    if(length(sl_id)>1){
      stop("Here")
    }
    res_source_name=c(res_source_name,source_name_list[sl_id])
    res_min_dist=c(res_min_dist,tmp_score[sl_id])
  }
  res_min_dist[1]=max(res_min_dist)
  
  
  res=data.frame(local_density=local_density,centroid=names(local_density),sigma=res_min_dist,source_name=res_source_name,stringsAsFactors = F)
  res=res[order(res$sigma,decreasing = T),]
  
  #ggplot(data=data$input_umap_centroid,aes(UMAP_1,UMAP_2,label=centroid))+geom_point()+geom_label(data=data$input_umap_centroid[data$input_umap_centroid$centroid %in% gsub("C","",res$centroid[1:3]),])
  
  #.de_dist=de_dist
  #de_dist=de_dist+t(de_dist)
  checkFlag=T
  counter=1
  D=res$centroid[1]
  while(checkFlag&counter<nrow(de_dist_sym)){
    counter=counter+1
    #print(counter)
    
    tmp=de_dist_sym[res$centroid[counter],D]
    
    if(min(tmp)>=pct_de_count_thr&res$sigma[counter]>=pct_de_count_thr){
      D=c(D,res$centroid[counter])
      
    } else {
      checkFlag=F
    }
  }
  
  for(center in D){
    res$source_name[res$centroid==center]=center
  }
  
  for(center in D){
    res$source_name[res$centroid==center]=center
    centerList=res$centroid[res$source_name==center]
    preCenters=""
    while(!identical(centerList,preCenters)){
      preCenters=centerList
      res$source_name[res$source_name %in% centerList]=center
      centerList=res$centroid[res$source_name==center]
    }
  }
  
  return(list(D=D,cluster_assignments=res))
}

.extra_sconline.AdjancyMatConstructor=function(meta_z_mat,pct_mat,de_dist,cosine_dist,pct_diff_thr,centers=NULL,symmetric=F){
  
  if(is.null(centers)){
    centers=row.names(pct_mat)
  }
  
  if(length(centers)==1){
    res=matrix(1,nrow=nrow(pct_mat),ncol=1)
    row.names(res)=row.names(pct_mat)
    colnames(res)=centers
    return(res)
  }
  
  diag(de_dist)=max(as.numeric(de_dist),na.rm = T)
  res_mat2=matrix(de_dist[,centers],ncol=length(centers))
  row.names(res_mat2)=row.names(de_dist)
  colnames(res_mat2)=centers
  
  res_mat2=t(apply(res_mat2,1,function(x) {
    x[x>min(x)]=NA
    x
  }))
  res_mat2[!is.na(res_mat2)]=1
  res_mat2[is.na(res_mat2)]=0
  if(sum(rowSums(res_mat2,na.rm = T)>1)>0){
    for(i in which(rowSums(res_mat2,na.rm = T)>1)){
      sl_ind=which(res_mat2[i,]>0)
      sl_scores=cosine_dist[i,sl_ind]
      sl_ind2=sl_ind[which(sl_scores==max(sl_scores))]
      res_mat2[i,sl_ind]=0
      res_mat2[i,sl_ind2]=1
    }
  }
  
  res_mat2_binary=matrix(pct_mat[centers,],nrow=length(centers))
  row.names(res_mat2_binary)=centers
  colnames(res_mat2_binary)=colnames(meta_z_mat)
  res_mat_residual=res_mat2 %*% res_mat2_binary
  
  res_mat1=pct_mat
  res_mat1[which(.extra_sconline.PctDiffscoreFn(res_mat1-res_mat_residual,pct_diff_thr = pct_diff_thr)*.sconline_extra_Pct2scoreFn(res_mat_residual,pct2_thr = pct2_thr)>0.99)]=0
  
  de_dist2=.sconline_extra_PctScoreFn(meta_z_mat = meta_z_mat,pct_mat=res_mat1,pct_mat_ref = pct_mat,pct2_thr = pct2_thr,pct_diff_thr = pct_diff_thr,centers=centers,symmetric = symmetric)
  de_dist2=de_dist2+max(as.numeric(de_dist2))*res_mat2
  de_dist2[cbind(match(centers,row.names(de_dist2)),match(centers,colnames(de_dist2)))]=max(as.numeric(de_dist2))
  thr=apply(de_dist2,1,min)
  res_mat3=t(apply(de_dist2,1,function(x) {
    x[x>min(x)]=NA
    x
  }))
  res_mat3[!is.na(res_mat3)]=1
  res_mat3[is.na(res_mat3)]=0
  if(sum(rowSums(res_mat3,na.rm = T)>1)>0){
    for(i in which(rowSums(res_mat3,na.rm = T)>1)){
      sl_ind=which(res_mat3[i,]>0)
      sl_scores=cosine_dist[i,sl_ind]
      sl_ind2=sl_ind[which(sl_scores==max(sl_scores))]
      res_mat3[i,sl_ind]=0
      res_mat3[i,sl_ind2]=1
    }
  }
  
  res_mat_net=res_mat2+res_mat3
  res_mat_net[cbind(match(centers,row.names(res_mat_net)),match(centers,colnames(res_mat_net)))]=0
  final_residual=apply(de_dist[,centers],1,min)
  thr=thr[match(names(final_residual),names(thr))]
  final_residual=data.frame(fdeg=final_residual,sdeg=thr,stringsAsFactors = F)
  return(list(res_mat_net=res_mat_net,residual=final_residual))
}

.sconline_ClusteringFn_archive=function(argList,min_marker_thr=20,sig1_thr=3,sig2_thr=1,pct_de_count_thr=1,pct_diff_thr=0.2,pct2_thr=0.3){
  
  #min_marker_thr=50;sig1_thr=3;sig2_thr=1;pct_de_count_thr=2;pct_diff_thr=0.2;pct2_thr=0.3
  
  meta_z_list=.sconline.fetch_data("meta_z",argList = argList)
  
  meta_z_mat=meta_z_list$meta_z
  cosine_dist=meta_z_list$cosine_dist
  de_dist=meta_z_list$de_pct_dist
  pct_mat=meta_z_list$med_pct.1
  
  j_dist=.extra_sconline.JaccardDistanceFn(meta_z_mat=meta_z_mat,meta_z_mat_bkg=NULL,sig1_thr=sig1_thr,sig2_thr=sig2_thr,mode="difference")
  j_dist=j_dist[row.names(de_dist),colnames(de_dist)]
  
  de_dist[j_dist<min_marker_thr]=0
  
  
  net_components=which(de_dist<=pct_de_count_thr,arr.ind = T)
  net_components=data.frame(from=row.names(de_dist)[net_components[,1]],to=colnames(de_dist)[net_components[,2]],stringsAsFactors = F)
  net_components=igraph::graph_from_data_frame(net_components,directed = F)
  net_components=igraph::components(net_components)
  net_components=data.frame(pseudocell=names(net_components$membership),membership=net_components$membership,stringsAsFactors = F)
  D=c()
  cluster_assignments=NULL
  for(icomponent in unique(net_components$membership)){
    sl_pseudocell=net_components$pseudocell[net_components$membership==icomponent]
    res_d=.extra_sconline.densityPeakClustering(meta_z_mat=meta_z_mat,cosine_dist=cosine_dist,de_dist=de_dist,pseudocells=sl_pseudocell,pct_de_count_thr = pct_de_count_thr,sig1_thr=sig1_thr,min_marker_thr=min_marker_thr)
    D=c(D,res_d$D)
    cluster_assignments=rbind(cluster_assignments,res_d$cluster_assignments[,c("centroid","source_name")])
  }
  
  
  #tst=.extra_sconline.AdjancyMatConstructor(meta_z_mat = meta_z_mat,pct_mat=pct_mat,de_dist=de_dist+t(de_dist),cosine_dist=cosine_dist,pct_diff_thr=pct_diff_thr,centers=D,symmetric=F)
  
  res_intermediates=NULL
  for(i in setdiff(row.names(pct_mat),D)){
    m2 =glmnet::glmnet(t(pct_mat[setdiff(D,i),]), y=pct_mat[i,], lower.limits = 0,lambda = 0,alpha=1,family="gaussian", intercept = FALSE)
    m2=coef(m2)
    m2=m2[row.names(m2)!="(Intercept)",]
    m2=m2[order(m2,decreasing = T)]
    
    m2 =glmnet::glmnet(t(meta_z_mat[names(m2)[1:2],]), y=meta_z_mat[i,], lower.limits = 0,lambda = 0,alpha=1,family="gaussian", intercept = FALSE)
    m2=coef(m2)
    m2=m2[row.names(m2)!="(Intercept)",]
    m2=m2[order(m2,decreasing = T)]
    
    
    res_intermediates=rbind(res_intermediates,data.frame(source1=names(m2)[1],source2=names(m2)[2],target=i,value1=m2[1],value2=m2[2],ratio=m2[1]/m2[2],cosine_count=sum(cosine_dist[i,names(m2)[1:2]]>0),stringsAsFactors = F))
  }
  
  res_intermediates=res_intermediates[which(res_intermediates$ratio<7&res_intermediates$cosine_count>1),]
  res_intermediates=res_intermediates[order(res_intermediates$ratio,decreasing = F),]
  
  res_intermediates=res_intermediates[!duplicated(.extra_sconline.NetIdFn(res_intermediates)),]
  tormList=c()
  for(i in 1:nrow(res_intermediates)){
    if(sum(de_dist[res_intermediates$target[i],cluster_assignments$centroid[cluster_assignments$source_name==res_intermediates$source1[i]]]>0)==0){
      tormList=c(tormList,i)
    }
    if(sum(de_dist[res_intermediates$target[i],cluster_assignments$centroid[cluster_assignments$source_name==res_intermediates$source2[i]]]>0)==0){
      tormList=c(tormList,i)
    }
  }
  if(length(tormList)>0){
    tormList=unique(tormList)
    res_intermediates=res_intermediates[-tormList,]
  }
  
  net_directional=data.frame(source=c(res_intermediates$source1,res_intermediates$source2),target=c(res_intermediates$target,res_intermediates$target),score=1,stringsAsFactors = F)
  
  
  #exNodes=row.names(tst$residual)[which(tst$residual[,1]==0)]
  
  #ggplot(data=data$input_umap_centroid,aes(UMAP_1,UMAP_2,label=centroid))+geom_point()+geom_label(data=data$input_umap_centroid[data$input_umap_centroid$centroid %in% gsub("C","",D),])
  
  gc()
  
  
  
  #cosine_dist=mymakeCosine(meta_z$meta_z,wJaccardDistance = F)
  #cosine_dist[D,D]
  
  #min_marker_thr=10
  
  
  
  #res_mat_net=.sconline_extra_AdjancyMatConstructor(pct_mat = pct_mat[D,],de_dist = de_dist[D,D],cosine_dist = cosine_dist[D,D])
  #res_mat_net=res_mat_net$res_mat_net
  #res_mat_net_indx=which(res_mat_net>0,arr.ind = T)
  #net=data.frame(target=row.names(res_mat_net)[res_mat_net_indx[,1]],source=row.names(res_mat_net)[res_mat_net_indx[,2]],score=1)
  #.extra_sconline.NetVisFn(net[!duplicated(myNetIdFn(net)),])
  
  #D=c("C16","C13","C126")
  
  colnames(cluster_assignments)=c("pseudocell","cluster")
  
  
  
  #cluster_assignments=.sconline.extra_clustering_assignments(res_mat=cosine_dist,centroids=D,not_sig_centroids_org=NULL,binary_thr=0.8,return_affinities=F)
  #cluster_assignments$cluster[cluster_assignments$pseudocell=="C160"]="C_14"
  #cluster_assignments$cluster[cluster_assignments$cluster=="C_3"]="C_1"
  #cluster_assignments$cluster[cluster_assignments$pseudocell %in% c("C26","C157")]="C_4"
  
  if(length(unique(cluster_assignments$cluster))<12&length(unique(cluster_assignments$cluster))>3){
    cluster_assignments$color=rcartocolor::carto_pal(length(unique(cluster_assignments$cluster)), "Vivid")[as.numeric(factor(cluster_assignments$cluster))]
  } else {
    cluster_assignments$color=hues::iwanthue(length(unique(cluster_assignments$cluster)))[as.numeric(factor(cluster_assignments$cluster))]
  }
  
  cluster_assignments$color[which(cluster_assignments$cluster=="C_0")]="gray"
  
  if(sum(is.na(cluster_assignments$color))>0){
    colPallette=cluster_assignments[cluster_assignments$pseudocell %in% D,]
    clust_affinities=.sconline.extra_clustering_assignments(res_mat=cosine_dist,centroids=D,not_sig_centroids_org=NULL,binary_thr=0.6,return_affinities=T)
    clust_affinities=clust_affinities[which(row.names(clust_affinities) %in% cluster_assignments$pseudocell[is.na(cluster_assignments$cluster)]),]
    
    col_set=.sconline_extra_clustering_color_gradient(affinity_mat = clust_affinities,col_rgb_pallette = colPallette)
    col_set=data.frame(pseudocell=row.names(clust_affinities),color2=col_set,stringsAsFactors = F)
    col_set=col_set[match(cluster_assignments$pseudocell,col_set$pseudocell),]
    cluster_assignments$color[is.na(cluster_assignments$color)]=col_set$color2[is.na(cluster_assignments$color)]
  }
  
  
  p2=""
  p1=""
  {
    
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
    input_umap_centroid=UMAP_centroid
    
    require(ggnetwork)
    res_net=which(upper.tri(de_dist)&de_dist<1&cosine_dist>0.8,arr.ind = T)
    res_net=data.frame(Fnode=row.names(de_dist)[res_net[,1]],Snode=colnames(de_dist)[res_net[,2]],score=de_dist[res_net],stringsAsFactors = F)
    res_net$score=(2.5*(res_net$score-0.5))^3
    
    
    #res_net=which(res_mat>=0.5,arr.ind = T)
    #res_net=data.frame(Fnode=row.names(res_mat)[res_net[,1]],Snode=colnames(res_mat)[res_net[,2]],score=res_mat[res_net],stringsAsFactors = F)
    #res_net$score=(2.5*(res_net$score-0.5))^3
    
    #id=paste0(res_net$Fnode,"_",res_net$Snode)
    #id[res_net$Fnode>res_net$Snode]=paste0(res_net$Snode,"_",res_net$Fnode)[res_net$Fnode>res_net$Snode]
    #res_net=res_net[!duplicated(id),]
    
    #res_net=sg2
    #res_net$Fnode=res_net$from
    #res_net$Snode=res_net$to
    #res_net$score=res_net$weight
    #res_net=res_net[,c("Fnode","Snode","score")]
    
    if(length(setdiff(colnames(de_dist),c(res_net$Fnode,res_net$Snode)))>0){
      res_net=rbind(res_net,data.frame(Fnode=setdiff(colnames(de_dist),c(res_net$Fnode,res_net$Snode)),Snode=setdiff(colnames(de_dist),c(res_net$Fnode,res_net$Snode)),score=0,stringsAsFactors = F))
    }
    net = network::network(res_net, directed = FALSE,matrix.type="edgelist")
    
    
    library(gplots)
    library(devtools)
    library(ggnetwork)
    require(network)
    require(sna)
    require(ggplot2)
    require(Matrix)
    
    netVerNames=network::network.vertex.names(net)
    resClustering=cluster_assignments
    #resClustering$cluster[which(resClustering$cluster=="C_0")]=NA
    net %v% "cluster" = resClustering$cluster[match(netVerNames,as.character(resClustering$pseudocell))]
    network::set.edge.attribute(net, "weight", res_net$score)
    
    
    centroid_layout=input_umap_centroid[match(gsub("C","",as.character(netVerNames)),as.character(input_umap_centroid$centroid)),-1]
    centroid_layout=as.matrix(centroid_layout)
    colnames(centroid_layout)=c("x","y")
    
    net=ggnetwork:::fortify.network(net,layout = centroid_layout)
    
    clusterCol=resClustering[!duplicated(resClustering$cluster),]
    
    pd_summary=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"))
    
    p1=.extra_sconline.NetVisFn(net=net_directional,input_pd=pd,input_umap_centroid=input_umap_centroid)
    
    resClustering$pseudocell=gsub("C","",resClustering$pseudocell)
    lda.training.data=merge(resClustering,input_umap_centroid,by.x="pseudocell",by.y="centroid")
    #lda.fit=MASS::lda(cluster~UMAP_1+UMAP_2,data=lda.training.data)
    #lda.pred=predict(lda.fit,pd_summary)
    
    knet=RANN::nn2(query=pd_summary[,c("UMAP_1","UMAP_2")],data = lda.training.data[,c("UMAP_1","UMAP_2")],k=12,eps=0)
    affinity_mat=matrix(0,nrow=nrow(knet$nn.idx),ncol=length(unique(lda.training.data$pseudocell)))
    colnames(affinity_mat)=unique(lda.training.data$pseudocell)
    
    tmp_affinities=t(apply(knet$nn.dists,1,function(x) exp((-1)*(x/x[2])^2)))
    tmp_affinities=t(apply(tmp_affinities,1,function(x) x/sum(x)))
    for(itr in 1:ncol(knet$nn.idx)){
      tst=.myOneHotFn(inputVector=factor(lda.training.data$pseudocell[knet$nn.idx[,itr]],levels=unique(lda.training.data$pseudocell)))
      tst=tst[,match(colnames(affinity_mat),colnames(tst))]
      tst=as.matrix(tst)
      tst=tmp_affinities[,itr]*tst
      affinity_mat=affinity_mat+tst
      #print(paste(itr,":",affinity_mat[6151,2]))
      rm(tst)
    }
    
    
    col_rgb_pallette=resClustering
    col_rgb_pallette=col_rgb_pallette[match(colnames(affinity_mat),col_rgb_pallette$pseudocell),]
    col_rgb1=matrix(grDevices::col2rgb(col_rgb_pallette$color)[1,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb2=matrix(grDevices::col2rgb(col_rgb_pallette$color)[2,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    col_rgb3=matrix(grDevices::col2rgb(col_rgb_pallette$color)[3,],nrow=nrow(affinity_mat),ncol=nrow(col_rgb_pallette),byrow = T)
    
    colnames(col_rgb1)=col_rgb_pallette$cluster
    colnames(col_rgb2)=col_rgb_pallette$cluster
    colnames(col_rgb3)=col_rgb_pallette$cluster
    
    col_rgb1=affinity_mat*col_rgb1
    col_rgb2=affinity_mat*col_rgb2
    col_rgb3=affinity_mat*col_rgb3
    
    col_rgb1=rowSums(col_rgb1)
    col_rgb2=rowSums(col_rgb2)
    col_rgb3=rowSums(col_rgb3)
    
    pd_summary$color=rgb(red=col_rgb1,green=col_rgb2,blue=col_rgb3,maxColorValue = 255)
    pd_summary_main=pd_summary
    
    scale_factor=net[!duplicated(net$vertex.names),]
    scale_factor$vertex.names=gsub("C","",scale_factor$vertex.names)
    scale_factor=merge(scale_factor,input_umap_centroid,by.x="vertex.names",by.y="centroid")
    scale_factor1=lm(x~UMAP_1,data=scale_factor)
    pd_summary$UMAP_1=predict(scale_factor1,newdata=pd_summary)
    scale_factor2=lm(y~UMAP_2,data=scale_factor)
    pd_summary$UMAP_2=predict(scale_factor2,newdata=pd_summary)
    
    pd_summary$xend=pd_summary$UMAP_1
    pd_summary$yend=pd_summary$UMAP_2
    
    #predicting the background color
    net=merge(net,cluster_assignments[,c("pseudocell","color")],by.x="vertex.names",by.y="pseudocell")
    #net$pct.sig="#FFFFFF"
    #net$pct.sig[net$vertex.names %in% cluster_assignments$pseudocell[cluster_assignments$cluster %in% cluster_assignments$cluster[cluster_assignments$pseudocell %in% res_clustering$centroids2]]]="#FCFAD2"
    
    
    if(F){
      p=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) + geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=color))+
        geom_edges( color = "grey50",aes(size=weight)) +
        geom_point(data=net[!duplicated(net$vertex.names),],size=2,aes(fill=color),shape=21)+
        geom_nodelabel(data=net,aes(color = color, label = as.character(cluster),fill=pct.sig),fontface = "bold")+
        theme_blank()+scale_size_continuous(range = c(0.04,1.3))+scale_color_identity()+scale_fill_identity()+theme(legend.position = "none")
      
    }
    
    p2=ggplot(net, aes(x = x, y = y, xend = xend, yend = yend))+ geom_point(data=pd_summary,aes(UMAP_1,UMAP_2,color=color))+
      geom_edges( color = "black") +
      #geom_point(data=net[!duplicated(net$vertex.names)&is.na(net$cluster),],aes(fill=color,size=2),shape=21)+
      geom_point(data=net[!duplicated(net$vertex.names)&!is.na(net$cluster),],aes(fill=color,size=5),shape=21,stroke = 1)+
      #geom_label(aes(label=cluster))+
      theme_blank()+scale_size_continuous(range = c(0.3,2))+scale_fill_identity()+scale_color_identity()+theme(legend.position = "none")
    
    
  }
  
  #p
  #ggsave("~/Desktop/FullCB_clusters.png",width = 8,height=8)
  
  #net=data.frame(target=row.names(res_mat_net)[res_mat_net_indx[,1]],source=row.names(res_mat_net)[res_mat_net_indx[,2]],score=1)
  #.extra_sconline.NetVisFn(net[!duplicated(myNetIdFn(net)),])
  
  #net=data.frame(source=c("C16","C9"),target=c("C9","C126"),score=1)
  
  #myNetVisFn2(net,data_pd=pd_summary_main,convertUMAP = T)
  #ggsave("~/Desktop/UBC_pattern.png",width = 8,height=8)
  
  
  
  
  return(list(net=net,plot1=p1,plot2=p2,cluster_assignments=cluster_assignments,centers=D))
}

.sconline.clusteringFn=function(argList,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,affinity_param=5,scaling_factor=10,createPlot=T,anno_col="anno_cellState",clustering_coefficient=NULL,n_clusters=NULL,cluster_selection_regularExp=NULL,inputPd=NULL,clustering_method="average"){

  #argList=.ArgList;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;affinity_param=5;scaling_factor=10;createPlot=T;anno_col='genotype_dissection';clustering_coefficient=1.05;cluster_selection_regularExp=NULL;n_clusters=NULL
  
  if(createPlot){
    require(ggraph)
    library(igraph)
    require(tidyverse)
    require(ape)
    require(seriation)
  } 
  
  if(!is.null(n_clusters)&!is.null(clustering_coefficient)){
    warning("Exact number of clusters is provided; clustering_coefficient is being ignored")
  }
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  if(sum(colnames(pd)==anno_col)==0){
    warning("Provided anno_col was not identified!")
    anno_col=NULL
  }
  
  if(is.null(n_clusters)&!is.null(anno_col)){
    
    n_clusters=clustering_coefficient*length(unique(pd[,anno_col]))
  }
  
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  
  
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(F){
    diff=prop_mat
    diff@x=rep(1,length(diff@x))
    diff=diff %*% t(diff)
    diff=sweep(diff,1,diag(diff),"/")
    diff=as.matrix(diff)
    diff[diff<1/scaling_factor]=1/scaling_factor
    diff=abs(log10(diff))
    diff=diff+scaling_factor*(1-exp(-1*de_pct_res/affinity_param))
  } else {
    myL2normFn=function(inputMat){
      prop_mat2=rowSums(inputMat^2)
      prop_mat2=sqrt(prop_mat2)
      res=Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
      return(res)
    }
    tst=prop_mat
    #tst@x=rep(1,length(tst@x))
    tst=myL2normFn(inputMat = tst)
    diff=t(tst)
    diff=tst %*% diff
    #diff@x=2*(exp(diff@x/max(quantile(diff@x,0.95),0.1))/(1+exp(diff@x/max(quantile(diff@x,0.95),0.1)))-0.5)
    diff=scaling_factor*(1-diff)
    diff=diff+scaling_factor*(1-exp(-1*(de_pct_res)/affinity_param))
    
  }
  
  if(F){
    matWeights=.myEffSizePropMat(prop_mat)
    matEffectiveSize=matWeights$effective_sample_size
  }
  
  
  diff=diff + t(diff)
  
  diff_clust=hclust(as.dist(diff),method = clustering_method)
  
  
  p_dendogram=NULL
  p_confusionMatrix=NULL
  p_anno_celllevel=NULL
  p_coverage=NULL
  cdf_analysis=NULL
  if(createPlot&!is.null(anno_col)){
    
    require(ggraph)
    library(igraph)
    require(tidyverse) 
    require(ape)
    require(seriation)
    
    #diff_clust=seriation:::reorder.hclust(x=diff_clust, dist=as.dist(diff), method = "OLO")
    adj2=prop_mat
    
    #adj2=diff2
    
    
    pd=pd[colnames(adj2),]
    pd[,anno_col]=gsub(" ",".",as.character(pd[,anno_col]))
    
    
    
    tmp= adj2 %*% as.matrix(.myOneHotFn(pd[,anno_col]))
    
    tmp=apply(tmp,1,function(x){
      y=which(x==max(x))[1]
      tmp2=data.frame(cluster=colnames(tmp)[y],purity=x[y])
      tmp2
    })
    
    tmp=do.call("rbind",tmp)
    
    if(!is.null(cluster_selection_regularExp)){
      sl_ind=which(grepl(tolower(cluster_selection_regularExp),tolower(tmp$cluster)))
      diff=diff[sl_ind,sl_ind]
      tmp=tmp[sl_ind,]
      
      diff_clust=hclust(as.dist(diff),method = "average")
      diff_clust=seriation:::reorder.hclust(x=diff_clust, dist=as.dist(diff), method = "OLO")
      pd=pd[which(grepl(cluster_selection_regularExp,tolower(pd[,anno_col]))),]
      
    }
    #tmp2=aggregate(purity~cluster,data=tmp,median)
    
    
    
    graph_tree = as.phylo(diff_clust)
    edges = graph_tree$edge
    node_names=graph_tree$tip.label[edges[,2]]
    node_names[is.na(node_names)]=edges[is.na(node_names),2]
    edges[,2]=node_names
    edges=as.data.frame(edges)
    colnames(edges)=c("from","to")
    # Let's add a column with the group of each name. It will be useful later to color points
    vertices = data.frame(
      name = unique(c(as.character(edges$from), as.character(edges$to)))
    ) 
    vertices$group=tmp$cluster[match(vertices$name,row.names(tmp))]
    #Let's add information concerning the label we are going to add: angle, horizontal adjustement and potential flip
    #calculate the ANGLE of the labels
    vertices$id=NA
    myleaves=which(!is.na(vertices$group))
    nleaves=length(myleaves)
    vertices$id[ myleaves ] = 1:nleaves
    vertices$angle= 90 - 360 * vertices$id / nleaves
    
    # calculate the alignment of labels: right or left
    # If I am on the left part of the plot, my labels have currently an angle < -90
    vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)
    
    # flip angle BY to make them readable
    vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)
    
    # Create a graph object
    vertices$group=as.factor(vertices$group)
    mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )
    
    # Make the plot
    p_dendogram=ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
      geom_edge_elbow(color="gray") +
      scale_edge_colour_distiller(palette = "RdPu") +
      geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=group, angle = angle, hjust=hjust, colour=group), size=2.7, alpha=1) +
      geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, alpha=0.2)) +
      scale_colour_manual(values= hues::iwanthue(length(unique(vertices$group)))) +
      scale_size_continuous( range = c(0.1,10) ) +
      theme_void() +
      theme(
        legend.position="none",
        plot.margin=unit(c(0,0,0,0),"cm"),
      ) +
      expand_limits(x = c(-1.5, 1.5), y = c(-1.5, 1.5))
    #ggsave(plot=p_dendogram,file="~/myBucket/torm.pdf")
    
    #Create confusion matrix
    
    
    d_conMat=cutree(diff_clust,k=n_clusters)
    
    
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
    
    colMax_vals=c()
    for(i in seq(1,ncol(prop_merged),5000)){
      tmp_max=as.numeric(qlcMatrix::colMax(prop_merged[,i:min(i+4999,ncol(prop_merged))]))
      colMax_vals=c(colMax_vals,as.numeric(tmp_max))
    }
    prop_merged =  prop_merged %*% Matrix::Diagonal(x = 1 / (colMax_vals+0.00000000001))
    #prop_mat=Matrix::drop0(prop_mat,0.1)
    prop_merged <- Matrix::Diagonal(x = 1 / rowSums(prop_merged)) %*% prop_merged
    
    matWeights=.myEffSizePropMat(prop_merged)
    matEffectiveSize=matWeights$effective_sample_size
    
    anno=prop_merged %*% as.matrix(.myOneHotFn(inputVector=pd[,anno_col]))
    matWeights=.myEffSizePropMat(prop_merged)
    matEffectiveSize=matWeights$effective_sample_size
    tmp_coverage=apply(anno,1,function(x) colnames(anno)[which(x==max(x))])
    tmp_coverage=data.frame(sconline=names(tmp_coverage),eff_size=matEffectiveSize,anno=tmp_coverage,stringsAsFactors = F)
    tmp_coverage=merge(tmp_coverage,as.data.frame(table(pd[,anno_col])),by.x="anno",by.y="Var1")
    p_coverage=ggplot(tmp_coverage,aes(anno,eff_size/Freq))+geom_boxplot()+theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.8,hjust=1,color="black"))
    
    anno=apply(anno,2,function(x) x/sum(x))
    
    anno=as.data.frame(anno)
    anno$ps=row.names(prop_merged)
    anno=reshape::melt(anno,id.vars="ps")
    anno=anno[order(anno$value,decreasing = T),]
    anno$ps=factor(as.character(anno$ps),levels=unique(anno$ps))
    anno=anno[order(anno$ps),]
    tst_lbl=as.character(anno$variable)[anno$value!=0]
    tst_lbl=tst_lbl[!duplicated(tst_lbl)]
    tst_lbl=c(tst_lbl,setdiff(as.character(anno$variable),tst_lbl))
    anno$variable=factor(as.character(anno$variable),levels=tst_lbl)
    recovered_clusters=anno[order(anno$value,decreasing = T),]
    recovered_clusters=recovered_clusters[!duplicated(recovered_clusters$ps),]
    
    tmp_anno=anno[order(anno$value,decreasing = T),]
    tmp_anno=tmp_anno[!duplicated(tmp_anno$ps),]
    
    
    
    tmp_recall=unlist(lapply(1:nrow(tmp_anno),function(x){
      length(unique(tmp_anno$variable[1:x]))/length(unique(anno$variable))
    }))
    tmp_recall=data.frame(purity=tmp_anno$value,Recall=tmp_recall,stringsAsFactors = F)
    
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=0.99)
    anno_m=colMax_vals_m %*% as.matrix(.myOneHotFn(inputVector=pd[,anno_col]))
    anno_m=as.matrix(Matrix::Diagonal(x=1/rowSums(anno_m)) %*% anno_m)
    anno_m=as.data.frame(anno_m)
    anno_m$cluster=row.names(prop_merged)
    anno_m=reshape::melt(anno_m,id.vars="cluster")
    anno_m=anno_m[order(anno_m$value,decreasing = T),]
    anno_m=anno_m[!duplicated(anno_m$cluster),]
    
    cdf_binary=unlist(lapply(1:nrow(anno_m), function(x) length(unique(anno_m$variable[1:x]))))
    cdf_binary=data.frame(purity=anno_m$value,Recall=cdf_binary/length(unique(anno$variable)),stringsAsFactors = F)
    cdf_binary=cdf_binary[order(cdf_binary$Recall,decreasing = T),]
    cdf_binary=cdf_binary[!duplicated(cdf_binary$purity),]
    cdf_binary=cdf_binary[order(cdf_binary$purity,decreasing = T),]
    cdf_binary$status="Binary"
    tmp_recall$status="Continuous"
    cdf_analysis=rbind(tmp_recall,cdf_binary)
    
    p_anno_celllevel=ggplot(data=anno,aes(ps,variable,fill=value))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+ggtitle(paste("Total # recoved clusters:",nrow(recovered_clusters)))
    #ggsave(plot=p_anno_celllevel,file="~/myBucket/torm.pdf")
    
    d_conMat=data.frame(tmp,cluster_hclust=d_conMat,stringsAsFactors = F)
    d_conMat=rbind(d_conMat)
    d_conMat=reshape2::dcast(cluster~cluster_hclust,data=d_conMat,value.var = "purity",fun.aggregate =length)
    d_conMat=d_conMat[match(unique(pd[,anno_col]),d_conMat[,1]),]
    row.names(d_conMat)=unique(pd[,anno_col])
    d_conMat=d_conMat[,-1]
    d_conMat[is.na(d_conMat)]=(0)
    d_conMat=sweep(d_conMat,2,colSums(d_conMat),"/")
    d_conMat$cluster=row.names(d_conMat)
    d_conMat=reshape::melt(d_conMat,id.vars="cluster")
    colnames(d_conMat)[colnames(d_conMat)=="variable"]="sconline_cluster"
    d_conMat=d_conMat[order(d_conMat$value,as.character(d_conMat$cluster),decreasing = T),]
    d_conMat$sconline_cluster=factor(as.character(d_conMat$sconline_cluster),levels=unique(d_conMat$sconline_cluster))
    d_conMat=d_conMat[order(d_conMat$sconline_cluster),]
    tst_lbl=d_conMat$cluster[d_conMat$value!=0]
    tst_lbl=tst_lbl[!duplicated(tst_lbl)]
    tst_lbl=c(tst_lbl,setdiff(as.character(d_conMat$cluster),tst_lbl))
    d_conMat$cluster=factor(as.character(d_conMat$cluster),levels=tst_lbl[!duplicated(tst_lbl)])
    p_confusionMatrix=ggplot(data=d_conMat,aes(sconline_cluster,cluster,fill=value))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+ylab("")
    #ggsave(plot=p_confusionMatrix,file="~/myBucket/torm.pdf")
  }
  
  return(list(cluster_object=diff_clust,distance_matrix=diff,plot_dendrogram=p_dendogram,plot_confusionMatrix=p_confusionMatrix,p_coverage=p_coverage,p_anno_celllevel=p_anno_celllevel,cdf_analysis=cdf_analysis,n_clusters=n_clusters))
}

.sconline.mastDEanalysis=function(argList,inputExpData,clust_obj,maineffect="genotype",covList=c("dissection"),randEffect_var="batch_merging",contrast_df,n_clusters=NULL,inputPd=NULL,verbose = TRUE,ncores=10,tol_level=0.95){
  #argList=.ArgList;inputExpData=xpo_data;clust_obj=clust_obj;maineffect="genotype";covList=c("dissection");randEffect_var="batch_merging";n_clusters=NULL;inputPd=NULL;verbose = TRUE;ncores=10;tol_level=0.95
  #contrast_df=data.frame(g1=c("WT","WT"),g2=c("HET","KO"),stringsAsFactors = F)
  #formula_str="~genotype+dissection+(1|batch_merging)"
  
  require(scater)
  require(MAST)
  
  #contrast_df=data.frame(g1=c("WT","WT"),g2=c("HET","KO"),stringsAsFactors = F)
  #formula_str="~genotype+dissection+(1|batch_merging)"
  
  
  
  if(ncores>1){
    options("mc.cores"=ncores)
  }
  
  if(is.null(n_clusters)){
    n_clusters=clust_obj$n_clusters
  }
  
  if(is.null(inputPd)){
    inputPd=.sconline.fetch_data("annotation",argList)
  }
  
  min_obs_size=length(unique(inputPd[,maineffect]))
  if(!is.null(covList)){
    for(icov in covList){
      min_obs_size=max(min_obs_size,length(unique(inputPd[,icov])))
    }
  }
  
  if(!is.null(randEffect_var)){
    min_obs_size=max(min_obs_size,length(unique(inputPd[,randEffect_var])))
  }
  
  min_obs_size=max(min_obs_size,10)
  
  d_conMat=cutree(clust_obj$cluster_object,k=n_clusters)
  
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(inputPd)]
  inputPd=inputPd[match(colnames(inputExpData),row.names(inputPd)),]
  
  {
    for(icov in 1:ncol(colData(inputExpData))){
      colData(inputExpData)[,icov]=gsub(" ",".",as.character(colData(inputExpData)[,icov]))
    }
    
  }
  
  
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=prop_mat[,row.names(inputPd)]
  prop_mat=Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),row.names(inputPd)]
  prop_merged=Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(inputPd),prop_m_hardCluster$j),]
  inputExpData$cluster_anno_res=inputPd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
  
  inputExpData=.mySplitObject(inputExpData,"cluster_anno_res")
  
  #formula_str=paste0(varList,collapse = "+")
  formula_str=paste0("~",maineffect)
  randomEffectModel=F
  if(!is.null(randEffect_var)){
    randomEffectModel=T
    randEffect_var=paste0("(1|",randEffect_var,")")
  }
  
  
  
  de_res=list()
  
  for(iclust in names(inputExpData)){
    for(icontrast in 1:nrow(contrast_df)){
      
      tmpexp=inputExpData[[iclust]]
      tmpexp=tmpexp[,colData(tmpexp)[,maineffect] %in% c(contrast_df[icontrast,1],contrast_df[icontrast,2])]
      if(length(unique(colData(tmpexp)[,maineffect]))>1){
        colData(tmpexp)[,maineffect]=factor(as.character(colData(tmpexp)[,maineffect]),levels=c(contrast_df[icontrast,1],contrast_df[icontrast,2]))
        tmpgroup=inputPd[match(colnames(tmpexp),row.names(inputPd)),]
        tmp_formula_str=formula_str
        for(icov in covList){
          if(length(unique(tmpgroup[,icov]))>1){
            tmp_formula_str=paste0(tmp_formula_str,"+",icov)
          }
        }
        
        if(randomEffectModel){
          tmp_formula_str=paste0(tmp_formula_str,"+",randEffect_var)
        }
        
        tmpexp_binary=counts(inputExpData[[iclust]])
        
        tmpexp_binary@x=rep(1,length(tmpexp_binary@x))
        tmpCount=rowSums(tmpexp_binary)
        tmpexp=tmpexp[tmpCount>min_obs_size,]
        
        
        tmpexp = scater::logNormCounts(tmpexp)
        zlm.res=list()
        for(igene in seq(1,nrow(tmpexp),5000)){
          sca = MAST::SceToSingleCellAssay(tmpexp[igene:min(igene+4999,nrow(tmpexp)),])
          
          if(randomEffectModel){
            #,fitArgsD=list(nAGQ=0)
            zlmCond <- MAST::zlm(as.formula(tmp_formula_str), sca = sca,method='glmer', ebayes=FALSE,parallel = T)
          } else {
            zlmCond <- MAST::zlm(as.formula(tmp_formula_str), sca = sca,parallel = T)
          }
          
          zlm.lr <- MAST::lrTest(zlmCond, maineffect)
          zlm.lr=zlm.lr[,,"Pr(>Chisq)"]
          zlm.lr=as.data.frame(zlm.lr)
          zlm.res=c(zlm.res,list(zlm.lr))
        }
        zlm.res=do.call("rbind",zlm.res)
        
        tmp_seurat=.extraExport2SeuratFn(tmpexp)
        tmp_seurat=Seurat::NormalizeData(tmp_seurat)
        tmp_fc_data=.myEvalMarkers(object=tmp_seurat, cells.1=colnames(tmp_seurat)[tmp_seurat@meta.data[,maineffect]==contrast_df[icontrast,2]], cells.2=colnames(tmp_seurat)[tmp_seurat@meta.data[,maineffect]==contrast_df[icontrast,1]], slot = "data", features = NULL,thresh.min = 0, pseudocount.use = 1,cluster_name=iclust)
        tmp_arranged=merge(data.frame(gene=row.names(tmp_fc_data),tmp_fc_data,stringsAsFactors = F),data.frame(gene=row.names(zlm.res),zlm.res,stringsAsFactors = F),by="gene",all=T)
        tmp_arranged$cluster=iclust
        tmp_arranged$contrast=paste0(contrast_df[icontrast,2],"_vs_",contrast_df[icontrast,1])
        de_res=c(de_res,list(tmp_arranged))
      }
      
    }
  }
  de_res=do.call("rbind",de_res)
  
  return(de_res)
}

.sconline.cluster_markers=function(argList,cluster_obj,inputExpData,n_clusters=NULL,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,only.pos.logFC=T){
  
  #argList=.ArgList;cluster_obj=.tmp;n_clusters=NULL;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;inputExpData=.xpo_data
  
  #cluster_obj: output of .sconline.clusteringFn()
  #n_clusters: number of desired clusters
  #sig1_thr: z-score threshold
  #pct2_thr: marker pct.2 threshold (% expressed in other clusters)
  #pct_diff_thr: pct.1 - pct.2 threshold
  #only.pos.logFC: report only upregulated genes in the cluster
  
  {
    require(ggraph)
    library(igraph)
    require(tidyverse)
    require(ape)
    require(seriation)
  } 
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  diff_clust=cluster_obj$cluster_object
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(!file.exists(.myFilePathMakerFn("res_marker_data",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))){
    prop_mat_c=prop_mat
    if(quantile(rowSums(prop_mat_c > 0,na.rm = T), 0.25) < (0.85*ncol(prop_mat))){
      prop_mat_c=Matrix::drop0(prop_mat_c)
      prop_mat_c@x=rep(1,length(prop_mat_c@x))
    }
    
    prop_mat_c=Matrix::drop0(1-prop_mat_c)
    prop_mat_c <- Matrix::Diagonal(x = 1 / (rowSums(prop_mat_c))) %*% prop_mat_c
    
    
    logNormData=t(Seurat:::NormalizeData.default(counts(inputExpData)[,colnames(prop_mat)],normalization.method = "RC",verbose = F))
    exp_binary=logNormData
    exp_binary@x=rep(1,length(exp_binary@x))
    
    fc_1=prop_mat %*% logNormData
    batch_size=100
    if(nrow(prop_mat_c) %% 100 ==0){
      batch_size=99
    }
    fc_2=parallel::mclapply(seq(1,nrow(prop_mat_c),batch_size),function(x) {prop_mat_c[x:min(x+batch_size-1,nrow(prop_mat_c)),] %*% logNormData},mc.cores = argList$ncores)
    fc_2=do.call("rbind",fc_2)
    
    fc_1@x=log2(fc_1@x+1)
    fc_2=log2(fc_2@x+1)
    
    logFC_res=fc_1 - fc_2#log2(prop_mat %*% logNormData+1) - log2(prop_mat_c %*% logNormData+1)
    
    
    pct.1=prop_mat %*% exp_binary
    pct.2=parallel::mclapply(seq(1,nrow(prop_mat_c),batch_size),function(x) {prop_mat_c[x:min(x+batch_size-1,nrow(prop_mat_c)),] %*% exp_binary},mc.cores = argList$ncores)#(prop_mat_c %*% exp_binary)
    pct.2=do.call("rbind",pct.2)
    pct.2=Matrix::drop0(pct.2)
    qsave(list(pct.1=pct.1,pct.2=pct.2,logFC=logFC_res),file=.myFilePathMakerFn("res_marker_data",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    marker_data=qread(.myFilePathMakerFn("res_marker_data",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    pct.1=marker_data$pct.1
    pct.2=marker_data$pct.2
    logFC_res=marker_data$logFC
  }
  
  
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  prop_merged=Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  
  
  z_mat=meta_data$meta_z
  matWeights=.myEffSizePropMat(prop_mat)
  matEffectiveSize=matWeights$effective_sample_size
  z_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  z_merged=sweep(z_merged,2,matEffectiveSize,"*")
  z_merged=Matrix::Diagonal(x=1/rowSums(z_merged)) %*% z_merged
  pct.1=z_merged %*% pct.1[,colnames(meta_data$meta_z)]
  pct.2=z_merged %*% pct.2[,colnames(meta_data$meta_z)]
  logFC_res=z_merged %*% logFC_res[,colnames(meta_data$meta_z)]
  z_merged= z_merged %*% meta_data$meta_z
  
  
  pct_diff=pct.1 - pct.2
  pct_diff=Matrix::drop0(pct_diff,tol=pct_diff_thr)
  
  
  
  pct_diff=pct_diff[,colnames(z_merged)]
  pct.1=pct.1[,colnames(z_merged)]
  pct.2=pct.2[,colnames(z_merged)]
  logFC_res=logFC_res[,colnames(z_merged)]
  z_merged=Matrix::drop0(z_merged,tol=sig1_thr)
  
  if(only.pos.logFC){
    z_merged@x[z_merged@x<0]=0
    z_merged=Matrix::drop0(z_merged,tol=sig1_thr)
  }
  
  pct_diff@x=rep(1,length(pct_diff@x))
  pct.2_binary=pct.2
  pct.2_binary@x[pct.2_binary@x>pct2_thr]=0
  pct.2_binary=Matrix::drop0(pct.2_binary)
  pct.2_binary@x=rep(1,length(pct.2_binary))
  z_merged=z_merged*pct_diff*pct.2_binary
  
  z_binary=z_merged=Matrix::drop0(z_merged)
  z_binary@x=rep(1,length(z_binary@x))
  pct.1=Matrix::drop0(pct.1*z_binary)
  pct.2=Matrix::drop0(pct.2*z_binary)
  logFC_res=Matrix::drop0(logFC_res*z_binary)
  pct_diff=Matrix::drop0(pct.1 - pct.2)
  
  my_sparseMat_summary=function(inputMat,val_col_name=NULL){
    res=as.data.frame(summary(inputMat))
    res$i=row.names(inputMat)[res$i]
    res$j=colnames(inputMat)[res$j]
    if(!is.null(val_col_name)){
      colnames(res)[colnames(res)=="x"]=val_col_name
    }
    return(res)
  }
  
  z_merged=my_sparseMat_summary(z_merged,val_col_name = "zscore")
  pct.1=my_sparseMat_summary(pct.1,val_col_name = "pct.1")
  pct.2=my_sparseMat_summary(pct.2,val_col_name = "pct.2")
  logFC_res=my_sparseMat_summary(logFC_res,val_col_name = "avg_log2FC")
  
  res=merge(z_merged,pct.1,by=c("i","j"))
  res=merge(res,pct.2,by=c("i","j"))
  res=merge(res,logFC_res,by=c("i","j"))
  colnames(res)[colnames(res)=="i"]="cluster"
  colnames(res)[colnames(res)=="j"]="gene"
  res=merge(data.frame(rowname=row.names(meta_data$fd),meta_data$fd,stringsAsFactors = F),res,by.y="gene",by.x="rowname",all.y=T)
  #res=res[,colnames(res)!="rowname"]
  
  return(res)
}

.sconline.clusteringVis=function(argList,cluster_obj,annoCol=NULL,n_clusters=NULL,pie_scale=1,tol_level=0.9){
  
  #argList=.ArgList;cluster_obj=.tmp;n_clusters=NULL
  
  {
    require(ggraph)
    library(igraph)
    require(tidyverse)
    require(ape)
    require(seriation)
  } 
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  
  diff_clust=cluster_obj$cluster_object
  
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  if(!is.null(annoCol)){
    
    
    
    if(is.numeric(pd[,annoCol])){
      x=matrix(pd[,annoCol],ncol=1)
    } else {
      x=as.matrix(.myOneHotFn(inputVector=pd[,annoCol]))
    }
    
    if(sum(is.na(pd$UMAP_1))>0|sum(colnames(pd)=="UMAP_1")==0){
      warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
    }
    res=prop_mat %*% x
    
    res=as.data.frame(res)
    
    res$pseudocell=row.names(prop_mat)
    pd$cluster_anno_res=pd[,annoCol]
    
  } else {
    prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
    prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
    prop_merged=Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
    #prop anno
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    colMax_vals_m=prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
    pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
    res=data.frame(pseudocell=names(d_conMat),cluster=paste0("C",d_conMat),stringsAsFactors=F)
  }
  
  #inputData=res;argList=argList;min_effective_size=5;pd=pd;cell_annoCol="cluster_anno_res";pie_scale=1
  p=.extra_sconline.visPseudocellAnno(inputData=res,argList=argList,min_effective_size=5,pd=pd,cell_annoCol="cluster_anno_res",pie_scale=pie_scale)
  
  
  
  return(p)
}

.sconline.MeanExp=function(argList,inputExpdata,normalization.method = "RC",return.only.unique.cellTypes=F,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2,affinity_param=5,scaling_factor=10,createPlot=T,anno_col="anno_cellState",clustering_coefficient=1.05,n_clusters=NULL,cluster_selection_regularExp=NULL){
  #normalization.method: NULL, "LogNormalize", "RC"
  #return.only.unique.cellTypes: filter the clustering results to include only one representative from each cell type.
  #argList=.ArgList;inputExpdata=data;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2;affinity_param=5;scaling_factor=10;createPlot=T;anno_col="anno_cellState";clustering_coefficient=1.05;n_clusters=NULL;normalization.method = "RC"
  
  if(!is.null(n_clusters)){
    warning("Exact number of clusters is provided; clustering_coefficient is being ignored")
  }
  #res_array=.sconline.fetch_data("sconline_arrays",argList = argList)
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  de_prefix=paste0("de_sig",sig1_thr,"_pctdiff",pct_diff_thr,"_pct2",pct2_thr)
  #pct_mat=meta_data$med_pct.1;argList = argList;meta_z_mat=meta_data$meta_z;sig1_thr=3;centers=NULL;pct2_thr=0.3;pct_diff_thr=0.2;symmetric=F
  if(!file.exists(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))){
    de_pct_res=.extra_sconline.PctScoreFn(pct_mat=meta_data$med_pct.1,argList = argList,meta_z_mat=meta_data$meta_z,sig1_thr=sig1_thr,centers=NULL,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr,symmetric=F)
    qsave(de_pct_res,file=.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  } else {
    de_pct_res=qread(.myFilePathMakerFn(de_prefix,argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  }
  
  
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(F){
    diff=prop_mat
    diff@x=rep(1,length(diff@x))
    diff=diff %*% t(diff)
    diff=sweep(diff,1,diag(diff),"/")
    diff=as.matrix(diff)
    diff[diff<1/scaling_factor]=1/scaling_factor
    diff=abs(log10(diff))
    diff=diff+scaling_factor*(1-exp(-1*de_pct_res/affinity_param))
  } else {
    myL2normFn=function(inputMat){
      prop_mat2=rowSums(inputMat^2)
      prop_mat2=sqrt(prop_mat2)
      res=Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
      return(res)
    }
    tst=prop_mat
    #tst@x=rep(1,length(tst@x))
    tst=myL2normFn(inputMat = tst)
    diff=t(tst)
    diff=tst %*% diff
    #diff@x=2*(exp(diff@x/max(quantile(diff@x,0.95),0.1))/(1+exp(diff@x/max(quantile(diff@x,0.95),0.1)))-0.5)
    diff=scaling_factor*(1-diff)
    diff=diff+scaling_factor*(1-exp(-1*(de_pct_res)/affinity_param))
    
  }
  
  
  diff=diff + t(diff)
  
  diff_clust=hclust(as.dist(diff),method = "average")
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  if(is.null(n_clusters)){
    d_conMat=cutree(diff_clust,k=clustering_coefficient*length(unique(pd[,anno_col])))
  } else {
    d_conMat=cutree(diff_clust,k=n_clusters)
  }
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  prop_merged=Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  anno=prop_merged %*% as.matrix(.myOneHotFn(inputVector=pd$anno_cellState))
  
  if(return.only.unique.cellTypes){
    tmp=apply(anno,1,function(x){
      y=which(x==max(x))[1]
      tmp2=data.frame(cluster=colnames(anno)[y],purity=x[y])
      tmp2
    })
    
    tmp=do.call("rbind",tmp)
    tmp$pseudocell=row.names(tmp)
    tmp=tmp[order(tmp$purity,decreasing = T),]
    tmp=tmp[!duplicated(tmp$cluster),]
    prop_merged=prop_merged[row.names(prop_merged) %in% tmp$pseudocell,]
    anno=prop_merged %*% as.matrix(.myOneHotFn(inputVector=pd$anno_cellState))
  }
  
  
  anno=as.data.frame(anno)
  anno$dominant_type=apply(anno,1,function(x) colnames(anno)[which(x==max(x))])
  anno$ps=row.names(prop_merged)
  
  
  colMax_vals_m=as.numeric(qlcMatrix::colMax(prop_merged))
  colMax_vals_m[which(colMax_vals_m==0)]=1
  colMax_vals_m=prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=0.95)
  anno_m=colMax_vals_m %*% as.matrix(.myOneHotFn(inputVector=pd$anno_cellState))
  anno_m=as.matrix(Matrix::Diagonal(x=1/rowSums(anno_m)) %*% anno_m)
  anno_m=as.data.frame(anno_m)
  anno_m$dominant_type=apply(anno_m,1,function(x) colnames(anno_m)[which(x==max(x))])
  anno_m$ps=row.names(prop_merged)
  
  colMax_vals_m=colMax_vals_m[,colSums(colMax_vals_m)==1]
  colMax_vals_m=as.data.frame(as.matrix(t(colMax_vals_m)))
  colMax_vals_m=apply(colMax_vals_m,1,function(x) colnames(colMax_vals_m)[x==1])
  colMax_vals_m=data.frame(cell=names(colMax_vals_m),cluster=colMax_vals_m,stringsAsFactors = F)
  tmp_pd=pd[match(colMax_vals_m$cell,pd$sample),]
  tmp_pd$anno_cluster_res=colMax_vals_m$cluster
  p=.myDimPlotFn(object=tmp_pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = 'anno_cluster_res')
  
  if(!is.null(normalization.method)){
    logNormData=t(Seurat:::NormalizeData.default(counts(inputExpdata),normalization.method = normalization.method,verbose = F))
  } else {
    logNormData=t(counts(inputExpdata))
  }
  
  x_exp=t(prop_merged %*% logNormData[colnames(prop_merged),])
  x_exp_binary=t(prop_m_hardCluster %*% logNormData[colnames(prop_merged),])
  
  return(list(cluster_plot=p,anno_hardCluster=anno_m,anno_softCluster=anno,meanExp_hardCluster=x_exp_binary,meanExp_softCluster=x_exp))
}

.sconline.markerPlot=function(argList,inputExpData,geneName,pseudocell_size_offset=0.2,pseudocell_size_pwr=1){
  #metric: zscore, pct.1
  library(dplyr)
  library(purrr)
  library(cowplot)
  library(patchwork)
  
  if(length(which(toupper(geneName)==toupper(row.names(inputExpData))))==0){
    stop("input gene was not found!")
  }
  
  geneName=row.names(inputExpData)[which(toupper(geneName)==toupper(row.names(inputExpData)))]
  
  meta_data=.sconline.fetch_data("meta_z",argList = argList)
  load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(pd)]
  pd=pd[match(colnames(inputExpData),row.names(pd)),]
  
  if(class(inputExpData)!="Seurat"){
    warning("converting the inputExpData to seurat object speeds up the process!")
    inputExpData=.extraExport2SeuratFn(inputData = inputExpData)
    inputExpData=Seurat::NormalizeData(inputExpData,verbose=F)
    
  }
  
  input_cell_bkg_umap=.scatterPlot_summary2d(object=pd,reductionCols=c("UMAP_1","UMAP_2"))
  
  umap.reduction <- Seurat:::CreateDimReducObject(embeddings = as.matrix(pd[,c("UMAP_1","UMAP_2")]), 
                                         key = "UMAP_", assay = "RNA", global = TRUE)
  
  inputExpData[["umap"]]=umap.reduction
  
  p1=Seurat::FeaturePlot(inputExpData, features=geneName,cols=c("lightgrey","#0D0887FF","#9C179EFF","#ED7953FF","#F0F921FF"))
  #p2=.my2dPlot_counts(inputPCA=pd,batch_values="batch_merging",reductionCols=c("UMAP_1","UMAP_2"),geneNameList=geneName,geneNameCol=NULL,expData=inputExpData,combine_figs=F)
  
  zscore_df=data.frame(pseudocell=row.names(meta_data$meta_z),zscore=meta_data$meta_z[,geneName],stringsAsFactors = F)
  pct_df=data.frame(pseudocell=row.names(meta_data$med_pct.1),pct.1=meta_data$med_pct.1[,geneName],stringsAsFactors = F)
  sc_df=merge(zscore_df,pct_df,by="pseudocell")
  sc_df=merge(sc_df,UMAP_centroid,by.x="pseudocell",by.y="centroid",all.x=T)
  sc_df=sc_df[sc_df$zscore!=0,]
  sc_df$scaled_pct.1=(sc_df$pct.1+pseudocell_size_offset)^pseudocell_size_pwr
  p3=ggplot(sc_df,aes(UMAP_1,UMAP_2,color=zscore,size=scaled_pct.1))+geom_point(data=input_cell_bkg_umap,aes(UMAP_1,UMAP_2),color="lightgrey",size=0.1)+geom_point()+theme_cowplot()+ theme(plot.title = element_text(hjust = 0.5))+ggtitle(geneName)+scale_color_gradientn(colors = c("lightblue","darkblue","black","yellow","orange","red"),breaks=c(min(sc_df$zscore),-2,0,3,5,max(sc_df$zscore)))+scale_size_identity()
  
  p=p1+p3+ plot_layout(nrow=1,ncol=2)
  
  
  return(p)
}


.sconline.clusteringARI=function(n_clusters,argList,cluster_obj,cell_anno_col="anno_cellState", inputPd=NULL,tol_level=0.9) {
  
  if(is.null(n_clusters)){
    n_clusters=cluster_obj$n_clusters
  }
  if(is.null(n_clusters)){
    stop("Number of desired clusters should be provided!")
  }
  
  
  diff_clust=cluster_obj$cluster_object
  
  d_conMat=cutree(diff_clust,k=n_clusters)
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  if(is.null(inputPd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPd
  }
  
  prop_mat=prop_mat[,colnames(prop_mat) %in% row.names(pd)]
  pd=pd[match(row.names(pd),colnames(prop_mat)),]
  
  
  prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
  prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
  prop_merged=Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
  #prop anno
  
  colMax_vals_m=qlcMatrix::colMax(prop_merged)
  colMax_vals_m=prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
  prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
  prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
  #prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
  prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
  prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
  pd$anno_cluster_res=paste0("C",as.character(prop_m_hardCluster$i))
  
  
  
  ARI = mclust::adjustedRandIndex(pd$anno_cluster_res, pd[,cell_anno_col])
  
  ari = mclust::adjustedRandIndex(pd$anno_cluster_res, pd[,cell_anno_col])
  nmi=aricode::NMI(pd$anno_cluster_res, pd[,cell_anno_col])
  ami=aricode::AMI(pd$anno_cluster_res, pd[,cell_anno_col])
  ariScores=data.frame(cluster_count=n_clusters,ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  
  
  tmp_pd=pd
  
  tmp_clust_size=as.data.frame(table(tmp_pd$anno_cellState))
  tmp_clust_size=tmp_clust_size[scale(tmp_clust_size[,2])<2,]
  sl_ind=tmp_pd$anno_cellState %in% as.character(tmp_clust_size[,1])
  tmp_pd=tmp_pd[sl_ind,]
  
  
  ari = mclust::adjustedRandIndex(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col])
  nmi=aricode::NMI(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col])
  ami=aricode::AMI(tmp_pd$anno_cluster_res, tmp_pd[,cell_anno_col])
  ariScores_balanced=data.frame(cluster_count=n_clusters,ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  
  return(list(all_cellTypes=ariScores,rm_LargeCellTypes=ariScores_balanced))
}

.sconline.seuratClusteringBN=function(argList,ncores=NULL){
  load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd=pd[match(row.names(pd),row.names(harmony_embeddings)),]
  library(future)
  options(future.globals.maxSize = 1000 * 1024^4)
  
  if(is.null(ncores)){
    plan("multicore", workers = argList$ncores)
  } else {
    plan("multicore", workers = ncores)
  }
  
  
  seurat_knn=Seurat:::FindNeighbors(object = harmony_embeddings[,1:argList$nPCs])
  seurat_clusters=Seurat:::FindClusters.default(object = seurat_knn[["snn"]],resolution = c(seq(0.01,0.1,0.01),seq(0.15,4,0.05)))
  
  clust_count=apply(seurat_clusters,2,function(x) length(unique(x)))
  seurat_clusters=seurat_clusters[,!duplicated(clust_count)]
  
  ariScores=lapply(1:ncol(seurat_clusters),function(x) {
    ari = mclust::adjustedRandIndex(seurat_clusters[,x], pd$anno_cellState)
    nmi=aricode::NMI(seurat_clusters[,x], pd$anno_cellState)
    ami=aricode::AMI(seurat_clusters[,x], pd$anno_cellState)
    data.frame(resolution=colnames(seurat_clusters)[x],cluster_count=length(unique(seurat_clusters[,x])),ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  })
  
  ariScores=do.call("rbind",ariScores)
  
  tmp_pd=pd
  
  tmp_clust_size=as.data.frame(table(tmp_pd$anno_cellState))
  tmp_clust_size=tmp_clust_size[scale(tmp_clust_size[,2])<2,]
  sl_ind=tmp_pd$anno_cellState %in% as.character(tmp_clust_size[,1])
  tmp_pd=tmp_pd[sl_ind,]
  tmp_seurat=seurat_clusters[sl_ind,]
  
  ariScores_balanced=lapply(1:ncol(seurat_clusters),function(x) {
    ari = mclust::adjustedRandIndex(tmp_seurat[,x], tmp_pd$anno_cellState)
    nmi=aricode::NMI(tmp_seurat[,x], tmp_pd$anno_cellState)
    ami=aricode::AMI(tmp_seurat[,x], tmp_pd$anno_cellState)
    data.frame(resolution=colnames(tmp_seurat)[x],cluster_count=length(unique(tmp_seurat[,x])),ari=ari,nmi=nmi,ami=ami,stringsAsFactors = F)
  })
  
  ariScores_balanced=do.call("rbind",ariScores_balanced)
  #.res_seurat=list(all_cellTypes=ariScores,rm_LargeCellTypes=ariScores_balanced)
  return(list(all_cellTypes=ariScores,rm_LargeCellTypes=ariScores_balanced))
  
}


#query = NULL; k.param = 20; prune.SNN = 1/15; nn.method = "annoy"; n.trees = 50; annoy.metric = "euclidean"; 
#nn.eps = 0; verbose = TRUE; force.recalc = FALSE; 
#cache.index = FALSE; index = NULL


.extra_sconline.sl_pseudocell.densityPeakFn_org=function (object,argList, query = NULL, k.param = 20, 
                                                      prune.SNN = 1/15, nn.method = "annoy", n.trees = 50, annoy.metric = "euclidean", 
                                                      nn.eps = 0, verbose = TRUE, force.recalc = FALSE, 
                                                      cache.index = FALSE, index = NULL) {
  plan("multicore", workers = argList$ncores)
  options(future.globals.maxSize = 1000 * 1024^4)
  cat("Performing densityPeak clustering\n")
  if (is.null(x = dim(x = object))) {
    warning("Object should have two dimensions, attempting to coerce to matrix", 
            call. = FALSE)
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.", 
            call. = FALSE)
    k.param <- n.cells - 1
  }
  
  if(is.null(query)){
    query=object
  }
  
  {
    nn.ranked <- Seurat:::NNHelper(data = object, query = query, k = k.param, 
                                   method = nn.method, n.trees = n.trees, searchtype = "standard", 
                                   eps = nn.eps, metric = annoy.metric, cache.index = cache.index, 
                                   index = index)
    nn.ranked <- Indices(object = nn.ranked)
  }
  
  snn.matrix <- Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
  rownames(x = snn.matrix) <- rownames(x = object)
  colnames(x = snn.matrix) <- rownames(x = object)
  
  tmp=(snn.matrix+t(snn.matrix))
  tmp=Matrix::Diagonal(x = 1 / (rowSums(tmp)+0.000000000001)) %*% tmp
  tmp=tmp %*% t(tmp)
  tmp=Matrix::Diagonal(x = 1 / (Matrix::diag(tmp)+0.000000000001)) %*% tmp
  tmp@x[which(tmp@x>=1)]=1
  
  tmp_density=Matrix::rowSums(snn.matrix)
  tmp_density=tmp_density[order(tmp_density,decreasing = T)]
  
  x=tmp[names(tmp_density),names(tmp_density)]
  
  x_dim=subset(summary(x), j < i)
  x_dim=sparseMatrix(x_dim[,"i"],x_dim[,"j"],x=1,dims = c(nrow(x), ncol(x)))
  x_dim=x*x_dim
  x_dim2=qlcMatrix::rowMax(x_dim)
  x_dim2=data.frame(density=tmp_density,similarity=as.matrix(x_dim2),cell=names(tmp_density),stringsAsFactors = F)
  x_dim2=x_dim2[which(x_dim2$density>1),]
  x_dim2=x_dim2[order(x_dim2$similarity*(-1),x_dim2$density,decreasing = T),]
  
  sim_thr=x_dim2$similarity[argList$internal_pseudocell_count]
  tst=x_dim2[x_dim2$similarity<=sim_thr,]
  sl_pseudocells=row.names(tst)
  
  
  
  return(sl_pseudocells)
}

.extra_sconline.sl_pseudocell.densityPeakFn_new=function (object,argList, query = NULL, k.param = 20, 
                                                      prune.SNN = 1/15, nn.method = "annoy", n.trees = 50, annoy.metric = "euclidean", 
                                                      nn.eps = 0, verbose = TRUE, force.recalc = FALSE, 
                                                      cache.index = FALSE, index = NULL,priority_list=NULL) {
  require(qs)
  plan("multicore", workers = argList$ncores)
  options(future.globals.maxSize = 1000 * 1024^4)
  if(is.null(priority_list)){
    cat("Performing densityPeak clustering\n")
  }
  
  if (is.null(x = dim(x = object))) {
    warning("Object should have two dimensions, attempting to coerce to matrix", 
            call. = FALSE)
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.", 
            call. = FALSE)
    k.param <- n.cells - 1
  }
  
  if(is.null(query)){
    query=object
  }
  
  
  if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))){
    idx=Seurat:::AnnoyBuildIndex(data = object, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = query,k=k.param,include.distance = T,search.k = -1)
    qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
  } else {
    nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
  }
  
  if(F){
    tst=Seurat::FindNeighbors(object,prune.SNN =0)
    tst_nn=as(tst[["nn"]],"dgCMatrix")
    tst_snn=as(tst[["snn"]],"dgCMatrix")
    x=tst_nn*tst_snn
    length(x@x)
    
    x2=adj*tst_nn
    length(x2@x)
    length(adj@x)
  }
  
  
  
  
  
  affinities=Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  affinities[,1]=affinities[,2]
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%k.param + 1
  x=as.numeric(t(affinities))
  adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(object),nrow(object)))
  rownames(adj) <- row.names(object)
  colnames(adj)=c(row.names(object))
  adj=Matrix::drop0(adj,tol=0.01)
  adj= Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  #adj=Matrix::drop0(adj,tol=quantile(adj@x,0.8))
  adj_t=t(adj)
  adj_t= Matrix::Diagonal(x=1/rowSums(adj_t)) %*% adj_t
  adj_t=adj %*% adj_t
  
  if(F){
    for(iitr in 1:4){
      adj=adj %*% adj_t*0.7+adj_t *0.3
      adj=Matrix::drop0(adj,tol=quantile(adj@x,0.1))
    }
  }
  
  #adj=adj %*% adj_t
  
  tst=Matrix::drop0(adj,tol=quantile(adj_t@x,0.8))
  tst@x=rep(1,length(tst@x))
  
  tmp_density=Matrix::rowSums(tst)
  if(F){
    if(!is.null(priority_list)){
      tmp_density[names(tmp_density) %in% priority_list]=max(tmp_density)+tmp_density[names(tmp_density) %in% priority_list]
    }  
  }
  
  tmp_density=tmp_density[order(tmp_density,decreasing = T)]
  if(F){
    pd=pd[match(names(tmp_density),row.names(pd)),]
    harmony_embeddings=harmony_embeddings[match(row.names(pd),row.names(harmony_embeddings)),]
    knet=RANN::nn2(data=harmony_embeddings,query = harmony_embeddings,k=100,eps=0)
    res_arranged=apply(as.data.frame((knet$nn.idx)),2,function(x) pd$anno_cellState[x])
    res_arranged=res_arranged[,-1]
    res_counts=apply(res_arranged,1,function(x) {x=as.numeric(table(x)); max(x)/sum(x)})
    summary(res_counts)
    summary(res_counts>0.9)
    summary(res_counts>0.5)
    pd$purity=res_counts
    
    
    pd2=pd[match(names(tmp_density),row.names(pd)),]
    pd2$anno_cellState=as.character(pd2$anno_cellState)
    head(pd2$anno_cellState,20)
    head(pd2$purity,20)
    pd2$purity[which(pd2$anno_cellState=="Group1")[1:3]]
    pd2$purity[which(pd2$anno_cellState=="Group3")[1:3]]
    which(pd2$anno_cellState=="Group1")[1:3]
  }
  
  
  res_dist=1
  adj2=adj[match(names(tmp_density),row.names(adj)),match(names(tmp_density),colnames(adj))]
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  while(nrow(adj2)>max(5*argList$internal_pseudocell_count,11000)){
    sl_ps_list=c()
    batch_size=max(10000,5*argList$internal_pseudocell_count)
    if(nrow(adj2) %% batch_size<3){
      batch_size=batch_size+3
    }
    for(i in seq(1,nrow(adj2),batch_size)){
      adj3=adj2[i:min(nrow(adj2),i+batch_size-1),]
      prop_mat2=myL2normFn(inputMat=adj3)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      
      
      
      tst=summary(c_c_aff)
      tst=tst[tst[,1]>tst[,2],]
      tst=sparseMatrix(i = tst[,1], j = tst[,2], x =tst[,3], dims = c(nrow(adj3),nrow(adj3)))
      row.names(tst)=row.names(adj3)
      tst=as.numeric(qlcMatrix::rowMax(tst))
      res_dist=tst
      names(res_dist)=row.names(adj3)
      res_dist=res_dist[order(res_dist,decreasing = F)]
      sl_ps_list=c(sl_ps_list,names(res_dist)[1:min(length(res_dist),argList$internal_pseudocell_count)])
    }
    
    sl_ps_list=sl_ps_list[!is.na(sl_ps_list)]
    adj2=adj2[sl_ps_list,sl_ps_list]
  }
  
  
  prop_mat2=myL2normFn(inputMat=adj2)
  c_c_aff=t(prop_mat2)
  c_c_aff=prop_mat2 %*% c_c_aff
  
  tst=summary(c_c_aff)
  tst=tst[tst[,1]>tst[,2],]
  tst=sparseMatrix(i = tst[,1], j = tst[,2], x =tst[,3], dims = c(nrow(adj2),nrow(adj2)))
  row.names(tst)=row.names(adj2)
  tst=as.numeric(qlcMatrix::rowMax(tst))
  res_dist=tst
  names(res_dist)=row.names(adj2)
  res_dist=res_dist[order(res_dist,decreasing = F)]
  
  if(F){
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    length(unique(pd$anno_cellState[pd$sample %in% names(res_dist)[1:200]]))
    setdiff(pd$anno_cellState,unique(pd$anno_cellState[pd$sample %in% names(res_dist)[1:200]]))
    
    table(pd$anno_cellState[pd$sample %in% names(res_dist)[1:400]])
    
  }
  
  sl_pseudocells=names(res_dist)[1:min(argList$internal_pseudocell_count,sum(res_dist<0.2))]
  
  return(sl_pseudocells)
}



.extra_sconline.sl_pseudocell.densityPeakFn=function (object,argList, query = NULL, k.param = 20, 
                                                      prune.SNN = 1/15, nn.method = "annoy", n.trees = 50, annoy.metric = "euclidean", 
                                                      nn.eps = 0, verbose = TRUE, force.recalc = FALSE, 
                                                      cache.index = FALSE, index = NULL,priority_list=NULL) {
  require(qs)
  plan("multicore", workers = argList$ncores)
  options(future.globals.maxSize = 1000 * 1024^4)
  if(is.null(priority_list)){
    cat("Performing densityPeak clustering\n")
  }
  
  if (is.null(x = dim(x = object))) {
    warning("Object should have two dimensions, attempting to coerce to matrix", 
            call. = FALSE)
    object <- as.matrix(x = object)
  }
  if (is.null(rownames(x = object))) {
    stop("Please provide rownames (cell names) with the input object")
  }
  n.cells <- nrow(x = object)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.", 
            call. = FALSE)
    k.param <- n.cells - 1
  }
  
  if(is.null(query)){
    query=object
  }
  
  
  if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))){
    idx=Seurat:::AnnoyBuildIndex(data = object, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = query,k=k.param,include.distance = T,search.k = -1)
    qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
  } else {
    nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
  }
  
  affinities=Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
  affinities@x=-1*affinities@x^2
  affinities@x=exp(affinities@x)
  affinities[,1]=affinities[,2]
  j <- as.numeric(t(nn.ranked.1$nn.idx))
  i <- ((1:length(j)) - 1)%/%k.param + 1
  x=as.numeric(t(affinities))
  adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(object),nrow(object)))
  rownames(adj) <- row.names(object)
  colnames(adj)=c(row.names(object))
  adj= Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
  
  adj_t=t(adj)
  adj_t= Matrix::Diagonal(x=1/rowSums(adj_t)) %*% adj_t
  adj=adj %*% adj_t
  
  tmp_density=Matrix::rowSums(adj)
  if(!is.null(priority_list)){
    tmp_density[names(tmp_density) %in% priority_list]=max(tmp_density)+tmp_density[names(tmp_density) %in% priority_list]
  }
  tmp_density=tmp_density[order(tmp_density,decreasing = T)]
  
  res_dist=1
  adj2=adj[match(names(tmp_density),row.names(adj)),match(names(tmp_density),colnames(adj))]
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  
  while(nrow(adj2)>max(5*argList$internal_pseudocell_count,11000)){
    sl_ps_list=c()
    batch_size=max(10000,5*argList$internal_pseudocell_count)
    if(nrow(adj2) %% batch_size<3){
      batch_size=batch_size+3
    }
    for(i in seq(1,nrow(adj2),batch_size)){
      adj3=adj2[i:min(nrow(adj2),i+batch_size-1),]
      prop_mat2=myL2normFn(inputMat=adj3)
      c_c_aff=t(prop_mat2)
      c_c_aff=prop_mat2 %*% c_c_aff
      
      
      
      tst=summary(c_c_aff)
      tst=tst[tst[,1]>tst[,2],]
      tst=sparseMatrix(i = tst[,1], j = tst[,2], x =tst[,3], dims = c(nrow(adj3),nrow(adj3)))
      row.names(tst)=row.names(adj3)
      tst=as.numeric(qlcMatrix::rowMax(tst))
      res_dist=tst
      names(res_dist)=row.names(adj3)
      res_dist=res_dist[order(res_dist,decreasing = F)]
      sl_ps_list=c(sl_ps_list,names(res_dist)[1:min(length(res_dist),argList$internal_pseudocell_count)])
    }
    
    sl_ps_list=sl_ps_list[!is.na(sl_ps_list)]
    adj2=adj2[sl_ps_list,sl_ps_list]
  }
  
  
  myL2normFn=function(inputMat){
    prop_mat2=rowSums(inputMat^2)
    prop_mat2=sqrt(prop_mat2)
    res=Matrix::Diagonal(x=1/(prop_mat2+0.00000001)) %*% inputMat
    return(res)
  }
  prop_mat2=myL2normFn(inputMat=adj2)
  c_c_aff=t(prop_mat2)
  c_c_aff=prop_mat2 %*% c_c_aff
  
  tst=summary(c_c_aff)
  tst=tst[tst[,1]>tst[,2],]
  tst=sparseMatrix(i = tst[,1], j = tst[,2], x =tst[,3], dims = c(nrow(adj2),nrow(adj2)))
  row.names(tst)=row.names(adj2)
  tst=as.numeric(qlcMatrix::rowMax(tst))
  res_dist=tst
  names(res_dist)=row.names(adj2)
  res_dist=res_dist[order(res_dist,decreasing = F)]
  
  if(F){
    load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
    
    length(unique(pd$anno_cellState[pd$sample %in% names(res_dist)[1:200]]))
    setdiff(pd$anno_cellState,unique(pd$anno_cellState[pd$sample %in% names(res_dist)[1:200]]))
    
    table(pd$anno_cellState[pd$sample %in% names(res_dist)[1:400]])
    
  }
  
  sl_pseudocells=names(res_dist)[1:min(argList$internal_pseudocell_count,sum(res_dist<0.2))]
  
  return(sl_pseudocells)
}




.extra_sconline.sl_pseudocell.kmeansFn=function(harmony_embeddings,argList){
  
  myPseudoAffinityMakerFn=function(harmony_embeddings,k.param=20,prune.SNN=1/15,n.trees = 50){
    #nn.ranked.1 <- RANN::nn2(harmony_embeddings, k = 10, eps = 0)
    
    if(!file.exists(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))){
      idx=Seurat:::AnnoyBuildIndex(data = harmony_embeddings, metric = "euclidean", 
                                   n.trees = n.trees)
      nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = harmony_embeddings,k=k.param,include.distance = T,search.k = -1)
      qsave(nn.ranked.1,file=.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    } else {
      nn.ranked.1=qread(.myFilePathMakerFn("knn_net",argList=argList,qsFormat = T))
    }
    
    nn.ranked=nn.ranked.1$nn.idx
    graph= Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(graph) <- rownames(x = harmony_embeddings)
    colnames(graph) <- rownames(x = harmony_embeddings)
    
    
    #graph=Seurat::FindNeighbors(harmony_embeddings,compute.SNN=T)
    #graph=as(graph[["snn"]], "dgCMatrix")
    
    #j <- as.numeric(t(nn.ranked.1$nn.idx))
    #i <- ((1:length(j)) - 1)%/%ncol(nn.ranked.1$nn.idx) + 1
    #k=1#as.numeric(t(nn.ranked.1$nn.dists))
    #graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(harmony_embeddings), nrow(harmony_embeddings)))
    
    #graph=graph+t(graph)
    #graph@x=rep(1,length(graph@x))
    
    if(F){
      graph=graph %*% t(graph)
      diag(graph)=0
      
      
      listCols_sparse<-function(X){
        #converts a sparse Matrix into a list of its columns
        #each list item contains only the nonzero elements of the column
        X<-as(X,"CsparseMatrix")
        res<-split(X@x, findInterval(seq_len(nnzero(X)), X@p, left.open=TRUE))
        names(res)<-colnames(X)
        res
      }
      
      colapply_sparse_nonzero<-function(X,FUN,...,mc.cores=1){
        #apply a function FUN to NONZERO elements of each column of sparse Matrix X
        #for an alternative that operates on all values, see colapply_full
        #mc: should parallel processing be used? Only recommended if FUN is slow
        #... additional args passed to mclapply or to FUN
        #this function always returns a list of length ncol(X)
        if(mc.cores>1){
          res=mclapply(listCols_sparse(X),FUN,...,mc.cores=mc.cores)
        } else {
          res=lapply(listCols_sparse(X),FUN,...)
        }
        res=unlist(res)
        X@x=res
        return(X)
      }
      
      affinities=colapply_sparse_nonzero(X=t(graph),FUN=function(x) exp((-3)*((max(x)-x)/(max(x)+1))^2),mc.cores=.ArgList$ncores)
      affinities=sqrt(t(affinities)*affinities)
      diag(affinities)=0
    } else {
      affinities=graph
    }
    
    return(affinities)
  }
  
  set.seed(1)
  doClustering=T
  itrClustering=0
  while(doClustering&itrClustering<10){
    itrClustering=itrClustering+1
    res_clust=kmeans(harmony_embeddings[,1:argList$nPCs],argList$internal_pseudocell_count,iter.max = 1000) #,algorithm = "Lloyd")
    
    if(sum(is.na(res_clust$centers))==0){
      doClustering=F
    }
  }
  
  if(sum(is.na(res_clust$centers))>0){
    stop("Error in identification of the pseudocells")
  }
  
  res_clusters=data.frame(sample=names(res_clust$cluster),cluster_id=res_clust$cluster,stringsAsFactors = F)
  save(res_clusters,file=.myFilePathMakerFn("kmeans_res_clusters",argList=argList))
  
  pseudo_affinity=myPseudoAffinityMakerFn(harmony_embeddings = harmony_embeddings)
  row.names(pseudo_affinity)=colnames(pseudo_affinity)=row.names(harmony_embeddings)
  
  sl_pseudo=NULL
  for(x in unique(res_clusters$cluster_id)){
    if(sum(res_clusters$cluster_id==x)>5){
      tmp=res_clusters$sample[res_clusters$cluster_id==x]
      tmp=rowSums(as.matrix(pseudo_affinity[tmp,tmp]))
      tmp=tmp[order(tmp,decreasing = T)]
      tmp=names(tmp)[1]
      x=data.frame(cluster=x,pseudocell=tmp,stringsAsFactors = F)
      sl_pseudo=rbind(sl_pseudo,x)
    }
  }
  
  #pca_centroid=res_clust$centers
  pca_centroid=harmony_embeddings[sl_pseudo$pseudocell,]
  row.names(pca_centroid)=sl_pseudo$cluster
  pca_centroid=pca_centroid[row.names(res_clust$centers)[row.names(res_clust$centers) %in% row.names(pca_centroid)],]
  
  output=sl_pseudo$pseudocell
  return(output)
}

#prepare the embedding space and defines the pseudocells in it using k-means
#embedding can be provided as input or it can be calculated by the method using pca and optionally harmony
#inputEmbeddings: a matrix of cell x embedding
#run_harmony: logical; if true, approach performs harmony to do batch correction
.sconline.embeddingFn=function(argList,inputEmbeddings,run_harmony,inputBatchCol='batch_merging',pd=NULL,saveFiles=T){
  #argList=.ArgList;inputEmbeddings=NULL;run_harmony=F;inputBatchCol='batch_merging';pd=NULL;saveFiles=F
  
  reRunCheck=tryCatch({load(.myFilePathMakerFn("pca_centroids",argList=argList));F}, error=function(e) {return(T)})
  
  if(!(reRunCheck|argList$newRun)){
    return("Done")
  }
  
  res="Done"
  if(is.null(inputEmbeddings)){
    cat('Creating Embeddings\n')
    
    library(future)
    plan("multicore", workers = min(parallel::detectCores(), 8, na.rm=T))
    # plan()
    options(future.globals.maxSize = 1000 * 1024^4)
    
    if(is.null(pd)){
      load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    }
    
    if(run_harmony){
      
      
      
      # CHANGED iter.max per https://github.com/immunogenomics/harmony/issues/25 and per line below
      # also changed lloyd
      #devtools::load_all("~/anaconda3/envs/psuedocells/harmony/")
      
      harmony_embeddings <- harmony::HarmonyMatrix(pca_res[,1:argList$nPCs], pd, 'batch_merging', do_pca = FALSE, verbose=FALSE)
      
    } else {
      harmony_embeddings=pca_res[,1:argList$nPCs]
    }
    
    if(saveFiles){
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    } else {
      res=harmony_embeddings
    }
    
    
    
  } else {
    inputEmbeddings=as.matrix(inputEmbeddings)
    
    if(sum(is.na(as.matrix(inputEmbeddings)))>0){
      stop("Error in the inputEmbeddings")
    }
    
    if(!is.null(pd)&!is.null(inputEmbeddings)){
      pca_res=inputEmbeddings
      save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
    } else {
      stop("Run .sconline.pca function first or provide the inputEmbeddings")
    }
    
    if(run_harmony){
      
      # CHANGED iter.max per https://github.com/immunogenomics/harmony/issues/25 and per line below
      # also changed lloyd
      #devtools::load_all("~/anaconda3/envs/psuedocells/harmony/")
      if(is.null(pd)){
        stop("phenoData should be provided")
      }
      inputEmbeddings <- harmony::HarmonyMatrix(inputEmbeddings[,1:argList$nPCs], pd, inputBatchCol, do_pca = FALSE, verbose=FALSE)
      
    }
    
    harmony_embeddings <- inputEmbeddings[,1:argList$nPCs]
    if(saveFiles){
      .mySaveFn(harmony_embeddings,file=.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
    }
    
    
  }
  
  argList$internal_pseudocell_count=round(nrow(harmony_embeddings)/argList$pseudocell_size)
  
  
  if(argList$pseudocell_selection_method=="kmeans"){
    pseudocell_names=.extra_sconline.sl_pseudocell.kmeansFn(argList=argList,harmony_embeddings=harmony_embeddings)#,kmeansMethod=kmeans_method)
    #.pseudocell_names=pseudocell_names
    if(F){
      knet=RANN::nn2(data=harmony_embeddings,query = harmony_embeddings[pseudocell_names,],k=100,eps=0)
      res_arranged=apply(as.data.frame((knet$nn.idx)),2,function(x) pd$anno_cellState[x])
      res_arranged=res_arranged[,-1]
      res_counts=apply(res_arranged,1,function(x) {x=as.numeric(table(x)); max(x)/sum(x)})
      summary(res_counts)
      summary(res_counts>0.9)
      summary(res_counts>0.5)
      pd_sl=pd[pseudocell_names,]
      table(pd_sl$anno_cellState[res_counts<.8])
    }
    
    #object=harmony_embeddings;argList= argList;priority_list=pseudocell_names
    pseudocell_names=.extra_sconline.sl_pseudocell.densityPeakFn(object=harmony_embeddings,argList= argList,priority_list=pseudocell_names)
  } else {
    pseudocell_names=.extra_sconline.sl_pseudocell.densityPeakFn(object=harmony_embeddings,argList= argList)
  }
  
  
  #pca_centroid=res_clust$centers
  
  pca_centroid=harmony_embeddings[pseudocell_names,]
  row.names(pca_centroid)=paste0("ps_",1:nrow(pca_centroid))
  sl_pseudo=data.frame(cluster=paste0("ps_",1:nrow(pca_centroid)),pseudocell=pseudocell_names,stringsAsFactors = F)
  save(pca_centroid,file=.myFilePathMakerFn("pca_centroids",argList=argList))
  save(sl_pseudo,file=.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  
  
  
  
  return("Done!")
}

#Constructs the argument list specifying the algorithm parameters
#Path to which save the results or load the results from
#n.prop: number of propagations to do
#newRun: to redo all analyses
#min_ds_size: datasets smaller than this size (based on the defined batch structure) are excluded from the analysis
.sconline.arg_creator=function(prefix,do.split.prop=T,min_cluster_size=20,saveDir,n.prop,nPCs,HVG_count=3,HVG_list=NULL,indScaling=T,ncores=min(parallel::detectCores(), 8, na.rm=T),pseudocell_selection_method="kmeans",pseudocell_size=200,input_highly_var_genes=NULL,newRun=F,min_ds_size=50,Leng200=F,include.singletons=T,colNormalize=T,singleton.method="fast"){
  
  dfSettings=NULL
  #HVG_count: for a gene to be selected as highly variable, how many datasets should support its variability. can't be larger than the number of datasets 
  #indScaling: performe independent scaling before pca analysis
  #nn.method: "fast", "snn"
  dfSettings=rbind(dfSettings,data.frame(commonExpressed=T,
                                         nPCs=nPCs,
                                         HVG_count=HVG_count,covariates="",indScaling=indScaling,removeHighExp=F,slMicrogliaClusters="",breakHammond=T,UMI_cor_thr=0.3,onlyRankBased=F,varScore.thr=5000,allGenesFraction=0.5,do.liger=F,stringsAsFactors = F))
  inputExpData=data
  runSetting=dfSettings
  runIndx=1
  DE_supportingFractionThr=0.1
  DE_n.adaptiveKernel=20;pval_thr=0.001
  
  
  if(sum(duplicated(row.names(inputExpData)))>0){
    print(paste(sum(duplicated(row.names(inputExpData))),"Duplicate gene ids were found! duplicates were randomly removed from the data"))
    inputExpData=inputExpData[!duplicated(row.names(inputExpData)),]
  }
  
  .ArgList=.extra_sconline.argFn(runIndx=runIndx,min_cluster_size=min_cluster_size,saveDir=saveDir,do.split.prop=do.split.prop,exNonMicCells=F,ncores=ncores,sensitiveSearch=1,include.singletons=include.singletons,colNormalize=colNormalize,includeHuman=F,includeMouse=F,FinnishSbjBased=F,uniformZscore=F,Leng200=Leng200,DE_supportingFractionThr=DE_supportingFractionThr,DE_n.adaptiveKernel=DE_n.adaptiveKernel,DE_nPropIter=n.prop,pseudocell_selection_method=pseudocell_selection_method,prop.n.neighbors=n.prop,dist_zscore_gamma=F,dist_zscore_norm=T,dist_zscore_nbinom=F,regularize=T,geoMean=F,prefix=prefix,newRun = newRun,inputDf=runSetting,internal_pseudocell_count=NULL,pseudocell_size=pseudocell_size,external_DE_path=NULL,external_DE_name=NULL,singleton.method=singleton.method)
  .ArgList$exclude_non_freq_pseudocells=F
  .ArgList$input_highly_var_genes=input_highly_var_genes
  .ArgList$min_ds_size=min_ds_size
  .ArgList$HVG_list=.ArgList$HVG_count
  
  # JONAH ADDED 
  # .ArgList$nHighlyVar = 3500
  
  #saveDir is the address to the location that main results are saved
  .ArgList$saveDir
  #saveDirGlobal includes address to the expression dataset
  .ArgList$saveDirGlobal
  .ArgList$newRun=newRun
  return(.ArgList)
}

#Performs HVG (Highly variable gene selection) and PCA analsyis
#if HVG_list is provided, it performs pca on different sets of HVGs
#batch_variable: the column name in the meta data that specifies the batch structure of the data
.sconline.pca=function(inputExpData,argList,batch_variable="batch_merging",organism="unknown",addAnno=F){
  library(rliger)
  library(Seurat)
  library(scater)
  library(scran)
  
  data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData,batch_variable=batch_variable,organism = organism,addAnno=addAnno)
  
  
  
  cat('Loading and preprocessing the datasets\n')
  
  if(!file.exists(.myFilePathMakerFn("pca_anno",argList = argList,pseudoImportant = F))|argList$newRun){
    
    cat('Selecting the highly variable genes\n')
    {
      data=.myHighVarGeneSlFn(data,dataorganism=organism,argList = argList)
      
    }
    
    tmp=data[c("varFeatures","allGenes" )]
    if(!file.exists(.myFilePathMakerFn("varGenes",argList=argList,varGenes = T))|argList$newRun){
      save(tmp,file=.myFilePathMakerFn("varGenes",argList=argList,varGenes =T))
    }
    
    tmp=data$data_m
    if(!file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))|argList$newRun){
      if(!is.null(data$data_m)){
        qs::qsave(tmp,file=.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
      }
    }
    
    cat('Scaling the data and PCA\n')
    .myPCAfn(data,argList = argList)
    
  } else {
    if(!argList$newRun){
      stop("pca file already exists. if need reRun, set the reRun variable in the argList to TRUE")
    }
  }
  
  return("Done!")
}

#running the sconline propagation method
#inputEmbeddings: optional; a matrix of [cells] x [PC] dimension
#inputExpData: optional; gene expression data; if specified, batch_variable should be defined too (specifies the column in the meta data that indicate the batch structure of the dataset)
#inputPhenoData should be provided if we want to do harmony on the inputEmbeddings
#batchVariable: the column name of the batchVariable in the phenoData
.sconline.runPropagation=function(argList,inputEmbeddings=NULL,inputPhenoData=NULL,inputExpData=NULL,input_UMAP_embedding=NULL,organism,batch_variable="batch_merging",run_harmony=F,addAnno=F,addClusteringModule=F,umap.method='umap-learn',extendedMode=F,L2Norm=T,mergePseudocells=T,generateUMAP=T,merging_strength=0.3,sig1_thr=3,pct2_thr=0.3,pct_diff_thr=0.2){
  
  
  #argList=.ArgList;inputEmbeddings=NULL;inputPhenoData=NULL;inputExpData=NULL;organism="Mouse";batch_variable="batch_merging";run_harmony=F;addAnno=F;addClusteringModule=F;L2Norm=T;mergePseudocells=T;generateUMAP=F;extendedMode=F;merging_strength=0.3;sig1_thr=3;pct2_thr=0.3;pct_diff_thr=0.2
  #run_harmony=T
  #umap.method="uwot";generateUMAP=T
  #pd=inputPhenoData;inputBatchCol=batch_variable
  
  if(!is.null(input_UMAP_embedding)){
    if(sum(colnames(input_UMAP_embedding) %in% c("UMAP_1","UMAP_2"))!=2){
      stop("input_UMAP_embedding is missing the UMAP_1 and UMAP_2 columns!")
    }
  }
  
  res_embeddings=.sconline.embeddingFn(argList,inputEmbeddings=inputEmbeddings,run_harmony=run_harmony,pd=inputPhenoData,inputBatchCol=batch_variable)
  
  res_umap=.sconline.umapFn(argList,umap.method=umap.method,generateUMAP = generateUMAP,input_UMAP_embedding=input_UMAP_embedding)
  #pd=.sconline.fetch_data("annotation",argList)
  #p=.myDimPlotFn(object=pd, dimCols = c("UMAP_1", "UMAP_2"), attCol = "anno_cellState")
  #ggsave(plot=p,file="~/myBucket/torm.pdf")
  
  if(!is.null(inputExpData)){
    data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData = inputExpData,batch_variable = batch_variable,organism = organism,addAnno=addAnno)
    data=data$data
    #expData = data
    res_prop=.myConcensusDEFn_step2(argList,expData = data,addClusteringModule=addClusteringModule,extendedMode = extendedMode,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr)
  } else {
    res_prop=.myConcensusDEFn_step2(argList,extendedMode = extendedMode,L2Norm=L2Norm,mergePseudocells=mergePseudocells,merging_strength=merging_strength,sig1_thr=sig1_thr,pct2_thr=pct2_thr,pct_diff_thr=pct_diff_thr)
  }
  
  return(res_prop)
}


#Creates a Seurat object from the cells
.sconline.create_seurat=function(argList,inputExpData=NULL,inputEmbeddings=NULL){
  if(is.null(inputExpData)){
    inputExpData=qs::qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
  } else if(class(inputExpData)=="SingleCellExperiment"){
    inputExpData=.extraExport2SeuratFn(inputExpData)
  } else if(class(inputExpData)=="liger"){
    inputExpData=.extra_sconline_LigerToExpSet(inputExpData)
    inputExpData=.extraExport2SeuratFn(inputExpData)
  } else {
    stop("unrecognized inputExpData format")
  }
  
  if(is.null(inputEmbeddings)){
    inputEmbeddings=.sconline.fetch_data(dataType = "embeddings",argList = argList)
  }
  
  pd=.sconline.fetch_data("annotation",argList = argList)
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(pd)]
  if(length(setdiff(colnames(inputExpData),row.names(pd)))>0){
    stop("Error! phenotype info were not found for some of the cells")
  } else {
    pd=pd[match(colnames(inputExpData),row.names(pd)),]
  }
  
  sl_meta_data=setdiff(colnames(pd),colnames(inputExpData@meta.data))
  if(length(sl_meta_data)>0){
    if(length(sl_meta_data)==1){
      tmp_pd=data.frame(new=pd[,sl_meta_data])
      colnames(tmp_pd)=sl_meta_data
      AddMetaData(inputExpData,tmp_pd) 
    } else {
      inputExpData=AddMetaData(inputExpData,pd[,sl_meta_data]) 
    }
  }
  
  umap_data=.sconline.fetch_data(dataType = "umap",argList = argList)
  umap_data=umap_data[match(colnames(inputExpData),row.names(umap_data)),]
  umap.reduction <- Seurat::CreateDimReducObject(embeddings = umap_data, 
                                         key = "UMAP_", assay = "RNA", global = TRUE)
  
  inputExpData[["umap"]]=umap.reduction
  
  #Adding in the embedding data
  embeddings=inputEmbeddings[match(colnames(inputExpData),row.names(inputEmbeddings)),]
  colnames(embeddings)=gsub("PC_","",colnames(embeddings))
  reduction.data <- Seurat::CreateDimReducObject(embeddings = embeddings, 
                                         assay = "RNA",  
                                         key = "PC_")
  inputExpData[["pca"]] <- reduction.data
  
  return(inputExpData)
}

#Propagates the cell annotations to the pseudocells
#res_prop: output of .sconline.runPropagation(); .sconline.fetch_data("sconline_arrays")
#collapse_datasets: logical; combine the results across datasets or provide the propagation results at the dataset level
#return_plot: return a ggplot object of the results
.sconline.anno2pseudocell=function(argList,annoCol="anno_cellState",collapse_datasets=T,return_plot=T,min_effective_size=5){
  
  #argList=.ArgList;annoCol="anno_cellState";collapse_datasets=T;return_plot=T;min_effective_size=5
  
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  
  prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
  prop_mat=Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
  
  matWeights=.myEffSizePropMat(prop_mat)
  
  matEffectiveSize=matWeights$effective_sample_size
  matWeights=matWeights$centroid_weights
  matWeights=matWeights[match(row.names(prop_mat),names(matWeights))]
  matEffectiveSize=matEffectiveSize[match(row.names(prop_mat),names(matEffectiveSize))]
  
  if(collapse_datasets&!argList$do.split.prop){
    
    if(is.numeric(pd[,annoCol])){
      x=matrix(pd[,annoCol],ncol=1)
    } else {
      x=as.matrix(.myOneHotFn(inputVector=pd[,annoCol]))
    }
    
    if(sum(is.na(pd$UMAP_1))>0|sum(colnames(pd)=="UMAP_1")==0){
      warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
    }
    res=prop_mat %*% x
    
    res=as.data.frame(res)
    
    res$pseudocell=row.names(prop_mat)
    res$effective_size=matEffectiveSize
    
  } else {
    res=list()
    for(ids in unique(pd$batch_merging)){
      tmp_pd=pd[pd$batch_merging==ids,]
      tmp_prop_mat=prop_mat[,match(tmp_pd$sample,colnames(prop_mat))]
      #tmp_prop_mat=Matrix::Diagonal(x = 1 / (rowSums(tmp_prop_mat)+0.000000000001)) %*% tmp_prop_mat
      tmp_effsize=matEffectiveSize*rowSums(tmp_prop_mat)
      
      tmp_weights=1/(1+exp((6-1/(1/(tmp_effsize+0.000000001)))))
      tmp_weights[tmp_effsize<4]=0
      
      if(is.numeric(tmp_pd[,annoCol])){
        x=matrix(tmp_pd[,annoCol],ncol=1)
      } else {
        x=as.matrix(.myOneHotFn(inputVector=tmp_pd[,annoCol]))
      }
      
      if(sum(is.na(tmp_pd$UMAP_1))>0|sum(colnames(tmp_pd)=="UMAP_1")==0){
        warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
      }
      tmp_res=tmp_prop_mat %*% x
      
      if(median(rowSums(tmp_prop_mat))>0.9){
        tmp_res=sweep(tmp_res,1,tmp_weight,"*")
      }
      
      tmp_res=as.data.frame(tmp_res)
      tmp_res$dataset=ids
      
      tmp_res$pseudocell=row.names(tmp_prop_mat)
      tmp_res$effective_size=tmp_effsize
      res=c(res,list(tmp_res))
      
    }
    
    
    res=do.call(eval(parse(text='plyr::rbind.fill')), res)
    res[is.na(res)]=0
    res=res[,c(setdiff(colnames(res),c("dataset","pseudocell","effective_size")),c("dataset","pseudocell","effective_size"))]
    
    if(collapse_datasets){
      res=res[,-which(colnames(res)=="dataset")]
      absent_pseudocells=apply(res[,-which(colnames(res) %in% c("pseudocell","effective_size"))],1,sum)
      absent_pseudocells=which(absent_pseudocells<0.2)
      if(length(absent_pseudocells)>0){
        res[absent_pseudocells,-which(colnames(res)=="pseudocell")]=NA
      }
      
      res=aggregate(.~pseudocell,data=res,function(x) sum(x,na.rm=T))
    }
  }
  
  
  
  p=""
  if(return_plot){
    p=.extra_sconline.visPseudocellAnno(inputData=res,argList = argList,min_effective_size=min_effective_size)
  }
  return(list(results=res,plot=p))
  
}


.sconline.anno2pseudocell_tmp=function(res_prop,argList,annoCol="anno_cellState",collapse_datasets=T,return_plot=T,min_effective_size=5,pd=NULL,subset_pd=F){
  if(is.null(pd)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  }
  
  
  res=list()
  for(i in 1:length(res_prop)){
    
    effective_size=res_prop[[i]]$matEffectiveSize
    
    tmp_pd=pd
    
    if(is.numeric(tmp_pd[,annoCol])){
      x=matrix(tmp_pd[,annoCol],ncol=1)
    } else {
      x=as.matrix(.myOneHotFn(inputVector=tmp_pd[,annoCol]))
    }
    tmp_pd=tmp_pd[match(colnames(res_prop[[i]]$prop_mat),row.names(pd)),]
    x=x[match(colnames(res_prop[[i]]$prop_mat),row.names(pd)),]
    if(sum(is.na(tmp_pd$UMAP_1))>0){
      warning("UMAP coordinates are missing for some annotations. there might an issue in the data!")
    }
    tmp_res=res_prop[[i]]$prop_mat %*% x
    tmp_effective_size=res_prop[[i]]$matEffectiveSize
    tmp_effective_size=tmp_effective_size[match(row.names(tmp_res),names(tmp_effective_size))]
    if(median(rowSums(res_prop[[i]]$prop_mat))>0.9){
      tmp_weight=res_prop[[i]]$matWeights
      tmp_weight=tmp_weight[match(row.names(tmp_res),names(tmp_weight))]
      
      tmp_res=sweep(tmp_res,1,tmp_weight,"*")
    }
    
    
    tmp_res=as.data.frame(tmp_res)
    if(sum(names(res_prop)=="zscore")==0){
      tmp_res$dataset=res_prop[[i]]$data$dsName
    } else {
      tmp_res$dataset=names(res_prop)[i]
    }
    
    tmp_res$pseudocell=row.names(res_prop[[i]]$prop_mat)
    tmp_res$effective_size=tmp_effective_size
    res=c(res,list(tmp_res))
  }
  
  res=as.data.frame(data.table::rbindlist(res))
  
  if(collapse_datasets){
    res=res[,-which(colnames(res)=="dataset")]
    absent_pseudocells=apply(res[,-which(colnames(res) %in% c("pseudocell","effective_size"))],1,sum)
    absent_pseudocells=which(absent_pseudocells<0.2)
    if(length(absent_pseudocells)>0){
      res[absent_pseudocells,-which(colnames(res)=="pseudocell")]=NA
    }
    
    res=aggregate(.~pseudocell,data=res,function(x) sum(x,na.rm=T))
  }
  
  p=""
  if(return_plot){
    if(subset_pd){
      p=.extra_sconline.visPseudocellAnno(inputData=res,argList = argList,min_effective_size=min_effective_size,pd=pd)
    } else {
      p=.extra_sconline.visPseudocellAnno(inputData=res,argList = argList,min_effective_size=min_effective_size,pd=NULL)
    }
    
  }
  return(list(results=res,plot=p))
  
}


#Fetch the data based on the argList
#if dataType=NULL, prints the datatypes that can be return by the object
.sconline.fetch_data=function(dataType=NULL,argList=NULL){
  result=""
  if(is.null(dataType)){
    cat("*********\nannotation\nembeddings\npca\nexpression\npca_centroids\numap_centroids\nmarkers\nsconline_arrays\nmeta_z\n*********\n")
  } else {
    if(is.null(argList)){
      stop("argList should be provided!")
    }
    result = switch(  
      dataType,  
      "annotation"= {
        if(file.exists(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))){
          load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
        } else {
          load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        }
        
        pd
      },  
      "pca_centroids"= {
        load(.myFilePathMakerFn("pca_centroids",argList=argList))
        pca_centroid
      },  
      "embeddings"= {
        load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
        harmony_embeddings
      },  
      "expression"= {
        expData=qs::qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
        expData
      },
      "umap_centroids"= {
        load(.myFilePathMakerFn("UMAP_centroid",argList=argList))
        UMAP_centroid
      },
      "pca"= {
        load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
        pca_res
      },
      "umap"= {
        load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
        pd=as.matrix(pd[,c("UMAP_1","UMAP_2")])
        pd
      },
      "markers"= {
        load(.myFilePathMakerFn("res_DE_wZscore",argList=argList,uniformImportant=T,propImportant = T))
        res_arranged
      },
      "sconline_arrays"= {
        data=qread(.myFilePathMakerFn("res_dataset_array",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
        data
      },
      "meta_z"= {
        data=qread(.myFilePathMakerFn("res_meta",argList=argList,uniformImportant=T,propImportant = T,qsFormat = T))
        data
      }
    )  
  }
  return(result)
}

#Correlation of pseudocells at the cell level
#res_prop: output of .sconline.runPropagation(); .sconline.fetch_data("sconline_arrays")
#plot_save_path: if provided, saves a heatmap in the path
.sconline.pseudocellCor=function(res_prop,argList,annoCol=NULL,plot_save_path=NULL){
  
  res_cor=lapply(1:length(res_prop$prop_mat),function(x){
    tmp_weights=.myEffSizePropMat(res_prop$prop_mat[[x]])$centroid_weights
    if(sum(tmp_weights<0.5)>0){
      tmp_weights[tmp_weights<0.5]=NA
    }
    as.data.frame(as.matrix(t(sweep(res_prop$prop_mat[[x]],1,tmp_weights,"*"))))
  }) 
  res_cor=as.matrix(data.table::rbindlist(res_cor))
  
  res_cor=cor(cor(res_cor,use="pairwise.complete.obs",method="spearman"))
  axis_col=rep("white",nrow(res_cor))
  if(!is.null(annoCol)){
    anno_data=.sconline.anno2pseudocell(res_prop = res_prop,argList=argList,annoCol = annoCol,collapse_datasets = T,return_plot = F)
    anno_data=anno_data$results[match(row.names(res_cor),anno_data$results$pseudocell),]
    anno_data=anno_data[,-which(colnames(anno_data) %in% c("pseudocell","effective_size"))]
    anno_data=colnames(anno_data)[apply(anno_data,1,function(x) {
      if(max(x)>0.5){
        which(x==max(x))[1]
      } else {
        #NA
        which(x==max(x))[1]
      }
    })]
    anno_data[is.na(anno_data)]="white"
    anno_data2=as.numeric(factor(anno_data))
    axis_col=c("gray",hues::iwanthue(length(unique(anno_data2))))[anno_data2]
    axis_col[anno_data=="white"]="white"
  }
  
  res_cor[which(is.na(res_cor))]=0
  res_cor[which(abs(res_cor)<0.05)]=0
  colramp = colorRampPalette(c("dodgerblue","black","yellow"))(9)
  
  singleton_pseudocells=which(rowSums(res_cor)==1)
  if(length(singleton_pseudocells)>0&length(singleton_pseudocells)<nrow(res_cor)){
    res_cor=res_cor[-singleton_pseudocells,-singleton_pseudocells]
    axis_col=axis_col[-singleton_pseudocells]
    
  }
  
  if(!is.null(plot_save_path)){
    pdf(file = "~/myBucket/torm.pdf")
    heatmap(res_cor,col=colramp, RowSideColors = axis_col,distfun = function(x) as.dist(1-cor(t(x))),hclustfun = function(x) hclust(x, method="ward.D2"), ColSideColors = axis_col, margins = c(5,10),scale="none")
    dev.off()
  }
  
  return(res_cor)
}

#Plot the pseudocells in the umap space
.sconline.plot_pseudocell_umap=function(argList,selectedPseudocells=NULL){
  require(ggplot2)
  umap_data=.sconline.fetch_data(dataType="umap_centroids",argList=argList)
  p=ggplot(umap_data,aes(UMAP_1,UMAP_2,label=centroid))
  
  if(!is.null(selectedPseudocells)){
    p=p+geom_point()+geom_label(data=umap_data[which(umap_data$centroid %in% as.character(selectedPseudocells)),],color="red")
  } else {
    p+geom_label()
  }
  
  p=p+theme_classic()
  
  return(p)
  
}


.sconline.HVGselection=function(argList,inputExpData,batch_variable="batch_merging",run_harmony=F,L2Norm=T){
  
  #argList=.ArgList;inputExpData=data;batch_variable="batch_merging";run_harmony=F;L2Norm=T
  #argList$input_highly_var_genes=NULL
  
  if(!is.null(argList$input_highly_var_genes)){
    stop("Highly variable genes are already specified in the argList")
  }
  
  argList$prefix="tmp_HVG_analysis"
  data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData=inputExpData,batch_variable=batch_variable,organism = "unknown",addAnno=F)
  cat('running pca step\n')
  
  cat('-- Selecting the highly variable genes\n')
  data=.myHighVarGeneSlFn(data,dataorganism="unknwon",argList = argList)
  
  cat('-- Scaling the data and PCA\n')
  pca_data=.myPCAfn(data,argList = argList,saveFiles=F)
  
  argList$pseudocell_size=round((sum(unlist(lapply(data$data,ncol))))/300)
  #inputExpData = data;inputEmbeddings=pca_data$pca_res;inputPhenoData=as.data.frame(pca_data$pd);run_harmony=run_harmony;batch_variable=batch_variable;generateUMAP = F;saveFiles = F;mergePseudocells=T;hierarchical_refinement=F;colNormalize=F
  prop_mat=.sconline.runPropagation(argList = argList,inputExpData = data,inputEmbeddings=pca_data$pca_res,inputPhenoData=as.data.frame(pca_data$pd),run_harmony=run_harmony,batch_variable=batch_variable,generateUMAP = F,saveFiles = F,mergePseudocells=T,hierarchical_refinement=F,colNormalize=F)
  
    
    return(res_prop)
}



.sconline.subset_purityAnalysis=function(argList,sl_cells=NULL,inputPCAdata=NULL,batch_variable="batch_merging",collapse_datasets=T,minCellCountThr=4,analysis_seed=1,extendedMode=F,cluster_count=NULL){
  #Running pca
  library(rliger)
  library(Seurat)
  library(scater)
  library(scran)
  
  myPrepDR=function (scaledData, features, verbose = TRUE) {
    
    data.use <- scaledData
    if (nrow(x = data.use) == 0) {
      stop("Data has not been scaled. Please run ScaleData and retry")
    }
    features.keep <- unique(x = features[features %in% rownames(x = data.use)])
    if (length(x = features.keep) < length(x = features)) {
      features.exclude <- setdiff(x = features, y = features.keep)
      if (verbose) {
        warning(paste0("The following ", length(x = features.exclude),
                       " features requested have not been scaled (running reduction without them): ",
                       paste0(features.exclude, collapse = ", ")))
      }
    }
    features <- features.keep
    # TODO jonah parallize buyt make sure chunked
    features.var <- apply(X = data.use[features, ], MARGIN = 1,
                          FUN = var)
    features.keep <- features[features.var > 0]
    if (length(x = features.keep) < length(x = features)) {
      features.exclude <- setdiff(x = features, y = features.keep)
      if (verbose) {
        warning(paste0("The following ", length(x = features.exclude),
                       " features requested have zero variance (running reduction without them): ",
                       paste0(features.exclude, collapse = ", ")))
      }
    }
    features <- features.keep
    features <- features[!is.na(x = features)]
    data.use <- data.use[features, ]
    return(data.use)
  }
  
  #organism="unknown";addAnno=F
  
  if(is.null(inputPCAdata)){
    inputExpData=.sconline.fetch_data("expression",argList = argList)
    if(!is.null(sl_cells)){
      print(paste0("Selecting ",sum(colnames(inputExpData) %in% sl_cells)," out of ",ncol(inputExpData)," cells."))
      if(sum(colnames(inputExpData) %in% sl_cells)==0){
        stop("None of the selected cells was identified in the dataset!")
      }
      inputExpData=inputExpData[,colnames(inputExpData) %in% sl_cells]
    } else {
      stop("sl_cells argument was not provided!")
    }
    
    data=.extra_sconline.exp_creatorFn(argList=argList,inputExpData=inputExpData,batch_variable=batch_variable,organism = "unknown",addAnno=F)
    argList2=argList
    argList2$HVG_list=argList2$HVG_count=min(max(round(length(data$data)/2),1),argList$HVG_count)
    argList2$input_highly_var_genes=NULL
    cat('Re-running pca step\n')
    
    cat('-- Selecting the highly variable genes\n')
    data=.myHighVarGeneSlFn(data,dataorganism="unknwon",argList = argList2)
    
    cat('-- Scaling the data and PCA\n')
    pca_data=.myPCAfn(data,argList = argList2,saveFiles=F)
    
  } else {
    argList$nPCs=ncol(inputPCAdata)-5
    pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=t(inputPCAdata),
                                                    features=colnames(inputPCAdata), verbose = FALSE),
                                    npcs = max(50,argList$nPCs+10))
    pca_data=list(pca_res=pca_res@cell.embeddings)
  }
  
  
  
  cat("Running Propagation step on the subset\n")
  
  set.seed(analysis_seed)
  supportingFractionThr=argList$DE_supportingFractionThr
  n.adaptiveKernel=argList$DE_n.adaptiveKernel
  nPropIter=argList$DE_nPropIter
  
  #load(.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
  load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  pd=pd[which(row.names(pd) %in% row.names(pca_data$pca_res)),]
  
  load(.myFilePathMakerFn("pca_centroids_assignments",argList=argList))
  sl_pseudo=sl_pseudo[which(sl_pseudo$pseudocell %in% row.names(pca_data$pca_res)),]
  pca_centroid=pca_data$pca_res[sl_pseudo$pseudocell,]
  row.names(pca_centroid)=sl_pseudo$cluster
  
  
  harmony_embeddings=pca_data$pca_res[row.names(pd),]
  if(sum(is.na(harmony_embeddings))>0){
    stop("Error in matching Names!")
  }
  
  pcaList=split(as.data.frame(harmony_embeddings),pd$batch_merging)
  
  
  dataArranged = parallel::mclapply(names(pcaList), function(thisExp){
    
    return(list(
      dsName=thisExp,
      pcaData=pcaList[[thisExp]]))
  })
  
  if(argList$do.split.prop&F){
    res=parallel::mclapply(dataArranged,.myConcensusDEFn_step2_detail_newprop3_final_v2_split_prop,centroidPCAdata=pca_centroid,argList=argList,n.adaptiveKernel=n.adaptiveKernel,nPropIter=1,n.neighbors=argList$prop.n.neighbors,mc.cores = argList$ncores)
  } else {
    #centroidPCAdata=pca_centroid;argList=argList;exCentroids=NULL;runIndx=1;n.neighbors=argList$prop.n.neighbors;batchPCAdata=harmony_embeddings
    res=.myConcensusDEFn_step2_detail_newprop3_final_v11subset(dataArranged=dataArranged,centroidPCAdata=pca_centroid,argList=argList,exCentroids=NULL,n.neighbors=argList$prop.n.neighbors,batchPCAdata=harmony_embeddings,returnPropMat=T,cluster_count=cluster_count)
  }
  
  return(res)
  
}


