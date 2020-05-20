#' @import Seurat
#' @import umap
#' @import SingleCellExperiment
#' @import dplyr
#' 



#' Dimensional reduction, UMAP calculation.
#'
#' @param object Seurat object
#' @param resolutions list of resolutions to test
#' @return Seurat object with PCA, UMAP and clustering calculations
#' @export
umapClustering<- function(object,resolutions=seq(0.1,0.5,0.1)){
  object <-RunPCA(object, npcs = 30, verbose = F)
  object <- RunUMAP(object,reduction = 'pca', dims = 1:25) 
  object<- FindNeighbors(object = object, reduction = "pca", dims = 1:25)
  object <- FindClusters(object = object, resolution = resolutions, print.output = FALSE, save.SNN = TRUE, force.recalc = TRUE)
  
  return(object)
}







#' Automatic assigner of cluster name.
#'
#'
#' @param object Seurat object
#' @param clusters.names list of names for each cluster in order
#' @param clusters.column metadata column containing the clusters to be named. Must have numbers. The names will be given in the order of the number/factors
#' @param output if an output with the cell information, including cluster name and number is desire please indicate output file
#' @return Seurat object but with the cluster names assigned.
#' @export
clusterIDassign<- function(object,clusters.names,clusters.column='SCT_snn_res.0.1',output=NULL,new_colName='cell_type'){
  temp<- object@meta.data[,clusters.column] #obtaining the cluster numbers
  temp<- as.numeric(as.character(temp))
  
  for(cluster in unique(sort(temp))){
    temp[temp==cluster]<- clusters.names[cluster+1] #assigning cluster name to each cluster number
    #e.g. cluster 0 = Monocytes
  }
  
  object@meta.data[,new_colName]<-factor(temp, levels=unique(clusters.names)) #setting levels
  
  if(!is.null(output)){
    #if an output file is specified, writing a csv with the metadata
    write.csv(object@meta.data,file = output)
  }
  
  return(object)
}










#' Re-calculated the cell cycle scoring of a seurat data set.
#'
#'
#' @param object Seurat object
#' @param SCTransform Logical. If it is necessary to regress the cell cycle and other variables in vars_regress. Default TRUE.
#' @param vars_regress variables to regress. Includes cell cycle and mitochondrial DNA, but more can be added or it can be modified.
#' @return Seurat object with calculated cell cycle scoring and regressed data.
#' @export
cellcycle_reScoring<- function(object, SCTransform=T, vars_regress=c("percent.mito", "S.Score", "G2M.Score")){
  DefaultAssay(object) = 'RNA'
  object <- NormalizeData(object, verbose = FALSE)
  s_genes <- cc.genes$s.genes[cc.genes$s.genes %in% rownames(object@assays$RNA@data)]
  g2m_genes <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% rownames(object@assays$RNA@data)]
  object = CellCycleScoring(object, s.features = s_genes, g2m.features = g2m_genes, set.ident=TRUE)
  
  if(SCTransform){
    object <- SCTransform(object, vars.to.regress = vars_regress, verbose = T, variable.features.n = nrow(object@assays$RNA@counts))
  }
  
  return(object)
}







#' Pre_process a seurat object by applying different cut-offs defined by the user.
#'
#' @param data_dir directory with the matrix, barcodes and features. String.
#' @param nGenes_min minimun number of genes a cell must have in order to keep it. integer.
#' @param nCounts_min minimun number of counts a cell must have in order to keep it. integer.
#' @param PercentMitoDNA maximun percentage of mitoDNA a cell must have in order to keep it. integer. e.g. 10
#' @return Filtered seurat object.
#' @export

pre_process<- function(data_dir,Project_name='Seurat_object'){
  
  temp<- Read10X(data_dir)
  df_t<- CreateSeuratObject(temp, min.cells = 3, min.features = 200, project =Project_name)
  
  #calculating mitochondrial DNA percentage
  df_t[["percent.mito"]] <- PercentageFeatureSet(object = df_t, pattern = "^MT-")
  #Filtering cells with low number of genes, counts and high mito DNA
  df_t <- subset(x = df_t, subset = nFeature_RNA > 300 & percent.mito < 10 & nCount_RNA > 200)
  
  
  return(df_t)
}


