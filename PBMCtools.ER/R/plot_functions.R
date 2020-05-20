
#' @import Seurat
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import cowplot
#' @import umap
#' @import gridExtra
#' @import grid
#' @import RColorBrewer
#' @import corrplot
#' @import purrr
#' @import comprehenr
#' @import SingleCellExperiment
#' @import ggrepel
#' @import cowplot
#' @import RColorBrewer
#' @import scales





#' Dimensional reduction plot for UMAP
#'
#' Graphs the output of the dimensional reduction technique UMAP on a 2D scatter plot where each point is a
#' cell and it's positioned based on the cell embeddings determined by the reduction technique.
#'
#' @param object Seurat object
#' @param group.by Name of one or more metadata columns to group (color) cells by. REQUIRED.
#' @param cols Colors for the different levels in the group.by parameter. The default is a color pallette.
#' @param na.value color for NA values
#' @param alpha.val Value for semitransparencies in the dots.
#' @param pt.size Dot size
#' @param label.DotSize Dot size in the legend.
#' @param legend.name Name to use as legend title.
#' @return ggplot 2D scatter plot where each point is a cell and it's positioned based on the cell embeddings determined by the reduction technique.
#' @export
umapPlot<- function(object,
                    group.by = NULL,
                    cols = NULL,

                    label = FALSE,
                    repel = FALSE,
                    na.value = 'grey50',
                    alpha.val=1,
                    pt.size = 0.5,
                    label.DotSize = 2,
                    palette=NULL,
                    legend.name=NULL
){


  #Checking if the user provided the groupby parameter
  if(is.null(group.by)){
    print('group required')
    return(NULL)
  }

  if(is.null(legend.name)){
    #If legend name is not provided it will use the group.by as parameter
    legend.name<- group.by
  }


  umap_coordinates <- data.frame(cbind(object@reductions$umap@cell.embeddings, object@meta.data))


  ##Creating a data frame with the median UMAP1 and UMAP2 coordinates for cells in each cluster
  ## These centroids will be used to add "ggrepel" labels for each cluster if it's requested and to add legends
  umap_centroids <- data.frame(cbind(t(sapply(levels(as.factor(umap_coordinates[,group.by])), FUN=function(n){apply(umap_coordinates[umap_coordinates[,group.by] == n,c("UMAP_1","UMAP_2")],MARGIN=2,median)}))),
                               id=levels(as.factor(umap_coordinates[,group.by])))

  #The column include in it's levels the unsassigned and doublet categories. We will only consider the categories
  # that are really in this columns
  umap_centroids<- umap_centroids[umap_centroids$id%in%as.vector(unique(umap_coordinates[,group.by])),]


  if(is.null(x = cols)){
    ## Defining a colour palette
    #cols <- RColorBrewer::brewer.pal(n=length(umap_centroids$id),name="Paired")
    if(is.null(palette)){ cols<- hue_pal()(length(umap_centroids$id)) } #not paletter specifies? using ggplot default colors
    else{cols <- RColorBrewer::brewer.pal(n=length(umap_centroids$id),name=palette)}#paletted specified. Using it.
  }


  ## Plotting cells, colour coding by group.by. Transparency and dot size as specified
  plot<- ggplot(umap_coordinates, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color=get(group.by)),alpha=alpha.val,size=pt.size) +
    #OPTIONAL: Adds ggrepel labels to each cluster
    scale_color_manual(values=cols,labels=as.character(umap_centroids$id),na.value = na.value, name = legend.name) +
    theme_bw() +
    theme(panel.grid = element_blank())+
    guides(colour = guide_legend(override.aes = list(size=label.DotSize, alpha = 1)))

  #if(labes_plot){
  #  plot+ geom_label_repel(data=umap_centroids, aes(x=UMAP_1, y=UMAP_2, label=id, fill=id))+
  #  scale_alpha_manual(values=cols) + #OPTIONAL: Only if adding ggrepel labels

  #}

  return(plot)
}










#' Dimensional reduction plot for UMAP coloring the expression of marker genes. Creates PNG plots.
#'
#' The gene panels are the following:
#' Granulocytes CEACAM8, LTF,MME,ITGAM,CD63,ENPP3
#' T cells: CD3D, CD4, CD8A, CD8B
#' NKT: CD1d, CD2, CD3, CD94 (KLRD1), CD160, KLRG1, TCRa/b, ZAP70
#' B cell: CD19, CD24, CD72, CD40
#'
#' Monocyte: CD14, CD16 (FCGR3A), CD2, CD11b (ITGAM), CD31 (PECAM1), CD56 (NCAM1), CD62L (SELL), CD115 (CSF1R), CD192 (CCR2), CX3CR1, CXCR3, CXCR4,CD45 (PTPRC)
#'   Classical: CD14++ CD16-, LYZ, CD45 (PTPRC), CD36(neg), CD64(FCGR1A, neg)
#'   Non classical: CD14(dim) CD16++, MS4A7,CD11c (ITGAX)
#'   Intermediate: CD14+, CD16+, LYZ, IFI30, CD74, CD86, HLA-DRA, S100A8, S100A10
#'
#' Macrophague: CD14, CD16 (FCGR3A), CD64 (FCGR1A), CD68, CD71 (TFRC) and CCR5
#'   M1 : CD86, CD80, CD68, MHCII, IL1R, TLR2, TLR4, iNOS (NOS2), SOCS3
#'   M2a: CD163, MHCII, SR, MMR/CD206 (MRC1), CD200R (CD200R1), TGM2, IL1R II
#'
#' Dendritic cells: absence of CD3 (T cell), CD14 (Monocyte), CD19 (B cell), CD56 (NCAM1, NK cells) and CD66b (CEACAM8, granulocytes)
#'   cDC: BDCA-1 (CD1C), CD8, CD8alpha, CD11b (ITGAM), CD11c (ITGAX), CD103 (ITGAE), CD205 (LY75), MHCII (HLA-DRA)
#'   pDC: BDCA-2 (CLEC4C), CD45RA (PTPRC), CD123 (IL3RA), ILT-7 (LILRA4),TLR7, TLR9, CD11clow (ITGAX), MHCII (HLA-DRA) low
#'
#' NK cells: KLRB1, CD56 (NCAM1), CD62L (SELL), CD57 (B3GAT1), CD94/NKG2A (KLRC1), CD16 (FCGR3A), NKR-P1 (CD161, KLRB1), NKG2D (KLRK1, CD314), NCR, FCGR3A TYROBP CLIC3 NKY7 PRF1 TRBC2(neg), TRAC(neg)
#' Granulocytes 'CEACAM8', 'LTF','MME','ITGAM','CD63','ENPP3'
#'
#'
#'
#' @param object Seurat object
#' @param filePreffix Preffix to use for the file name, including the path to the saving folder
#' @param pt.size Dot size
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if cells expressing given feature are getting buried.
#' @export


markerGenesPlot<- function(object, filePreffix, pt.size=NULL,order=TRUE){
  ####T cells: CD3D, CD4, CD8A, CD8B
  plot<- FeaturePlot(object = object, features = c('CD3D', 'CD4','CD8A', 'CD8B'), min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_Tcell_markerGenes.png", sep=""),res = 300, width=3000,height=3000)
  print(plot)
  dev.off()

  ##NKT: CD1d, CD2, CD3, CD94 (KLRD1), CD160, KLRG1, TCRa/b, ZAP70
  plot<- FeaturePlot(object = object, features = c('CD3D', 'CD8A', 'CD8B', 'CD1D', 'CD2', 'KLRD1', 'CD160', 'KLRG1', 'ZAP70'), min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_NKT_markerGenes.png", sep=""),res = 300, width=4000,height=3000)
  print(plot)
  dev.off()



  ####B cell: CD19, CD24, CD72, CD40
  plot<- FeaturePlot(object = object, features = c('CD19', 'CD24', 'CD72', 'CD40'), min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_Bcells_markerGenes.png", sep=""),res = 300, width=3000,height=3000)
  print(plot)
  dev.off()




  ####Monocyte: CD14, CD16 (FCGR3A), CD2, CD11b (ITGAM), CD31 (PECAM1), CD56 (NCAM1), CD62L (SELL), CD115 (CSF1R), CD192 (CCR2), CX3CR1, CXCR3, CXCR4,CD45 (PTPRC)

  # classical: CD14++ CD16-, LYZ, CD45 (PTPRC), CD36(neg), CD64(FCGR1A, neg)
  # non classical: CD14(dim) CD16++, MS4A7,CD11c (ITGAX)
  # intermediate: CD14+, CD16+, LYZ, IFI30, CD74, CD86, HLA-DRA, S100A8, S100A10


  plot<- FeaturePlot(object = object, features = c('CD14','FCGR3A','LYZ','ITGAM','ITGAX','PECAM1','PTPRC','CD2', 'NCAM1', 'SELL', 'CSF1R', 'CCR2', 'CX3CR1', 'CXCR3', 'CXCR4'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_Monocytes_markerGenes.png", sep=""),res = 300, width=6000,height=4000)
  print(plot)
  dev.off()



  plot<- FeaturePlot(object = object, features = c('CD14','FCGR3A','LYZ','PTPRC','CD36'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_MonocytesClassic_markerGenes.png", sep=""),res = 300, width=3500,height=3000)
  print(plot)
  dev.off()



  plot<- FeaturePlot(object = object, features = c('CD14','FCGR3A','MS4A7','ITGAX'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_MonocytesNonClassic_markerGenes.png", sep=""),res = 300, width=3000,height=3000)
  print(plot)
  dev.off()


  plot<- FeaturePlot(object = object, features = c('CD14','FCGR3A','LYZ','IFI30','CD74','CD86','HLA-DRA','S100A8','S100A10'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_MonocytesInter_markerGenes.png", sep=""),res = 300, width=5000,height=3000)
  print(plot)
  dev.off()


  ####Macrophague: CD14, CD16 (FCGR3A), CD64 (FCGR1A), CD68, CD71 (TFRC) and CCR5
  plot<- FeaturePlot(object = object, features = c('CD14', 'FCGR3A', 'FCGR1A', 'CD68', 'TFRC', 'CCR5'),   min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_Macrophague_markerGenes.png", sep=""),res = 300, width=3000,height=3000)
  print(plot)
  dev.off()


  #M1 : CD86, CD80, CD68, MHCII, IL1R, TLR2, TLR4, iNOS (NOS2), SOCS3
  plot<- FeaturePlot(object = object, features = c('CD14', 'FCGR3A', 'CD86', 'CD80', 'TLR2', 'TLR4', 'NOS2','SOCS3', 'HLA-DRA', 'HLA-DRB1'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_MacrophagueM1_markerGenes.png", sep=""),res = 300, width=5200,height=3000)
  print(plot)
  dev.off()



  #M2a: CD163, MHCII, SR, MMR/CD206 (MRC1), CD200R (CD200R1), TGM2, IL1R II
  plot<- FeaturePlot(object = object, features = c('CD14', 'FCGR3A', 'CD163', 'HLA-DRA', 'HLA-DRB1',  'MRC1', 'CD200R1', 'TGM2', 'IL1R2'), min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_MacrophagueM2a_markerGenes.png", sep=""),res = 300, width=5500,height=3000)
  print(plot)
  dev.off()



  ####Dendritic cells: absence of CD3 (T cell), CD14 (Monocyte), CD19 (B cell), CD56 (NCAM1, NK cells) and CD66b (CEACAM8, granulocytes)
  plot<- FeaturePlot(object = object, features = c('CD3D', 'CD14', 'CD19', 'NCAM1'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_DCnegativeMarkers_markerGenes.png", sep=""),res = 300, width=3000,height=3000)
  print(plot)
  dev.off()

  #cDC: BDCA-1 (CD1C), CD8, CD8alpha, CD11b (ITGAM), CD11c (ITGAX), CD103 (ITGAE), CD205 (LY75), MHCII (HLA-DRA)
  plot<- FeaturePlot(object = object, features = c('ID2', 'IRF8','BATF3','ZEB2','IRF4','CD1C', 'CD8A', 'ITGAM', 'ITGAX', 'ITGAE', 'LY75', 'HLA-DRA'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_DCc_markerGenes.png", sep=""),res = 300, width=5500,height=3000)
  print(plot)
  dev.off()


  #pDC: BDCA-2 (CLEC4C), CD45RA (PTPRC), CD123 (IL3RA), ILT-7 (LILRA4),TLR7, TLR9, CD11clow (ITGAX), MHCII (HLA-DRA) low
  plot<- FeaturePlot(object = object, features = c( 'TCF4','CLEC4C', 'PTPRC', 'NRP1','IL3RA','IRF8','IRF4','ZEB2', 'LILRA4', 'TLR7', 'TLR9', 'ITGAX','HLA-DRA'), min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_DCp_markerGenes.png", sep=""),res = 300, width=5500,height=3000)
  print(plot)
  dev.off()



  ####NK cells: KLRB1, CD56 (NCAM1), CD62L (SELL), CD57 (B3GAT1), CD94/NKG2A (KLRC1), CD16 (FCGR3A), NKR-P1 (CD161, KLRB1), NKG2D (KLRK1, CD314), NCR
  plot<-   FeaturePlot(object = object, features = c('KLRB1', 'NCAM1', 'SELL', 'B3GAT1', 'FCGR3A',  'KLRK1', 'NCR1', 'NCR2', 'NCR3'), min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_NK_markerGenes.png", sep=""),res = 300, width=5000,height=3000)
  print(plot)
  dev.off()


  ####NK cells: eddie markers FCGR3A TYROBP CLIC3 NKY7 PRF1 TRBC2(neg) TRAC(neg)
  plot<- FeaturePlot(object = object, features = c('FCGR3A', 'TYROBP', 'CLIC3', 'NKY7', 'PRF1','TRBC2','TRAC'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_NK_markerGenes_2.png", sep=""),res = 300, width=4500,height=3000)
  print(plot)
  dev.off()



  ####Granulocytes 'CEACAM8', 'LTF','MME','ITGAM','CD63','ENPP3'
  plot<- FeaturePlot(object = object, features = c('CEACAM8', 'LTF','MME','ITGAM','CD63','ENPP3'), min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_Granulocytes_markerGenes.png", sep=""),res = 300, width=3000,height=3000)
  print(plot)
  dev.off()

}














#' Dimensional reduction plot for UMAP coloring the expression of pre-defined marker genes for the monocyte and NK subsets. Creates PNG plots
#'
#' The gene panels are the following:
#' Monocyte: CD14, CD16 (FCGR3A), CD2, CD11b (ITGAM), CD31 (PECAM1), CD56 (NCAM1), CD62L (SELL), CD115 (CSF1R), CD192 (CCR2), CX3CR1, CXCR3, CXCR4,CD45 (PTPRC)
#'   Classical: CD14++ CD16-, LYZ, CD45 (PTPRC), CD36(neg), CD64(FCGR1A, neg)
#'   Non classical: CD14(dim) CD16++, MS4A7,CD11c (ITGAX)
#'   Intermediate: CD14+, CD16+, LYZ, IFI30, CD74, CD86, HLA-DRA, S100A8, S100A10
#'
#' NK cells: KLRB1, CD56 (NCAM1), CD62L (SELL), CD57 (B3GAT1), CD94/NKG2A (KLRC1), CD16 (FCGR3A), NKR-P1 (CD161, KLRB1), NKG2D (KLRK1, CD314), NCR, FCGR3A TYROBP CLIC3 NKY7 PRF1 TRBC2(neg), TRAC(neg)
#' Dendritic cells: absence of CD3 (T cell), CD14 (Monocyte), CD19 (B cell), CD56 (NCAM1, NK cells) and CD66b (CEACAM8, granulocytes)
#'   cDC: BDCA-1 (CD1C), CD8, CD8alpha, CD11b (ITGAM), CD11c (ITGAX), CD103 (ITGAE), CD205 (LY75), MHCII (HLA-DRA)
#'   pDC: BDCA-2 (CLEC4C), CD45RA (PTPRC), CD123 (IL3RA), ILT-7 (LILRA4),TLR7, TLR9, CD11clow (ITGAX), MHCII (HLA-DRA) low
#'

#'
#' @param object Seurat object
#' @param filePreffix Preffix to use for the file name, including the path to the saving folder
#' @param pt.size Dot size
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if cells expressing given feature are getting buried.
#' @param cell_type Monocyte or NK
#' @return PNG files with the plots of the markers genes of Monocytes or NK cells. DCs included.
#' @export
markerGenesPlot_MonoNK<- function(object, filePreffix, pt.size=NULL,order=TRUE,cell_type='Monocyte'){

  if(cell_type== 'Monocyte'){
    ####Monocyte: CD14, CD16 (FCGR3A), CD2, CD11b (ITGAM), CD31 (PECAM1), CD56 (NCAM1), CD62L (SELL), CD115 (CSF1R), CD192 (CCR2), CX3CR1, CXCR3, CXCR4,CD45 (PTPRC)

    # classical: CD14++ CD16-, LYZ, CD45 (PTPRC), CD36(neg), CD64(FCGR1A, neg)
    # non classical: CD14(dim) CD16++, MS4A7,CD11c (ITGAX)
    # intermediate: CD14+, CD16+, LYZ, IFI30, CD74, CD86, HLA-DRA, S100A8, S100A10


    plot<- FeaturePlot(object = object, features = c('CD14','FCGR3A','LYZ','ITGAM','ITGAX','PECAM1','PTPRC','CD2', 'NCAM1', 'SELL', 'CSF1R', 'CCR2', 'CX3CR1', 'CXCR3', 'CXCR4'),  min.cutoff = "q9",pt.size=pt.size,order=order)

    png(file = paste(filePreffix,"_Monocytes_markerGenes.png", sep=""),res = 300, width=6000,height=4000)
    print(plot)
    dev.off()



    plot<- FeaturePlot(object = object, features = c('CD14','FCGR3A','LYZ','PTPRC','CD36'),  min.cutoff = "q9",pt.size=pt.size,order=order)

    png(file = paste(filePreffix,"_MonocytesClassic_markerGenes.png", sep=""),res = 300, width=3500,height=3000)
    print(plot)
    dev.off()



    plot<- FeaturePlot(object = object, features = c('CD14','FCGR3A','MS4A7','ITGAX'),  min.cutoff = "q9",pt.size=pt.size,order=order)

    png(file = paste(filePreffix,"_MonocytesNonClassic_markerGenes.png", sep=""),res = 300, width=3000,height=3000)
    print(plot)
    dev.off()


    plot<- FeaturePlot(object = object, features = c('CD14','FCGR3A','LYZ','IFI30','CD74','CD86','HLA-DRA','S100A8','S100A10'),  min.cutoff = "q9",pt.size=pt.size,order=order)

    png(file = paste(filePreffix,"_MonocytesInter_markerGenes.png", sep=""),res = 300, width=5000,height=3000)
    print(plot)
    dev.off()
  }else if (cell_type== 'NK') {

    ####NK cells: KLRB1, CD56 (NCAM1), CD62L (SELL), CD57 (B3GAT1), CD94/NKG2A (KLRC1), CD16 (FCGR3A), NKR-P1 (CD161, KLRB1), NKG2D (KLRK1, CD314), NCR
    plot<-   FeaturePlot(object = object, features = c('KLRB1', 'NCAM1', 'SELL', 'B3GAT1', 'FCGR3A',  'KLRK1', 'NCR1', 'NCR2', 'NCR3'), min.cutoff = "q9",pt.size=pt.size,order=order)

    png(file = paste(filePreffix,"_NK_markerGenes.png", sep=""),res = 300, width=5000,height=3000)
    print(plot)
    dev.off()


    ####NK cells: eddie markers FCGR3A TYROBP CLIC3 NKY7 PRF1 TRBC2(neg) TRAC(neg)
    plot<- FeaturePlot(object = object, features = c('FCGR3A', 'TYROBP', 'CLIC3', 'NKY7', 'PRF1','TRBC2','TRAC'),  min.cutoff = "q9",pt.size=pt.size,order=order)

    png(file = paste(filePreffix,"_NK_markerGenes_2.png", sep=""),res = 300, width=4500,height=3000)
    print(plot)
    dev.off()
  }


  #cDC: BDCA-1 (CD1C), CD8, CD8alpha, CD11b (ITGAM), CD11c (ITGAX), CD103 (ITGAE), CD205 (LY75), MHCII (HLA-DRA)
  plot<- FeaturePlot(object = object, features = c('ID2', 'IRF8','BATF3','ZEB2','IRF4','CD1C', 'CD8A', 'ITGAM', 'ITGAX', 'ITGAE', 'LY75', 'HLA-DRA'),  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_DCc_markerGenes.png", sep=""),res = 300, width=5500,height=3000)
  print(plot)
  dev.off()


  #pDC: BDCA-2 (CLEC4C), CD45RA (PTPRC), CD123 (IL3RA), ILT-7 (LILRA4),TLR7, TLR9, CD11clow (ITGAX), MHCII (HLA-DRA) low
  plot<- FeaturePlot(object = object, features = c( 'TCF4','CLEC4C', 'PTPRC', 'NRP1','IL3RA','IRF8','IRF4','ZEB2', 'LILRA4', 'TLR7', 'TLR9', 'ITGAX','HLA-DRA'), min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_DCp_markerGenes.png", sep=""),res = 300, width=5500,height=3000)
  print(plot)
  dev.off()


}










#' Dotplot of marker genes. Creates PNG plots
#'
#'
#' @param object Seurat object
#' @param filePreffix Preffix to use for the file name, including the path to the saving folder
#' @param group.by Name of the metadata column to group for.
#' @param cols Colors to scale the plot. Vector of two colors.
#' @return PNG file. Dot plot of the marker genes of each cell type.
#' @export
markerGenesDotplot<- function(object, filePreffix, group.by=NULL, cols=c("lightgrey", "#cc0000")){
  ####T cells: CD3D, CD4, CD8A, CD8B
  plot<- DotPlot(object, features = c('CD3D', 'CD4','CD8A', 'CD8B'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_Tcell_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()

  ##NKT: CD1d, CD2, CD3, CD94 (KLRD1), CD160, KLRG1, TCRa/b, ZAP70
  plot<- DotPlot(object, features = c('CD3D', 'CD8A', 'CD8B', 'CD1D', 'CD2', 'KLRD1', 'CD160', 'KLRG1', 'ZAP70'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_NKT_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()



  ####B cell: CD19, CD24, CD72, CD40
  plot<- DotPlot(object, features = c('CD19', 'CD24', 'CD72', 'CD40'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_Bcells_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()




  ####Monocyte: CD14, CD16 (FCGR3A), CD2, CD11b (ITGAM), CD31 (PECAM1), CD56 (NCAM1), CD62L (SELL), CD115 (CSF1R), CD192 (CCR2), CX3CR1, CXCR3, CXCR4,CD45 (PTPRC)

  # classical: CD14++ CD16-, LYZ, CD45 (PTPRC), CD36(neg), CD64(FCGR1A, neg)
  # non classical: CD14(dim) CD16++, MS4A7,CD11c (ITGAX)
  # intermediate: CD14+, CD16+, LYZ, IFI30, CD74, CD86, HLA-DRA, S100A8, S100A10

  plot<- DotPlot(object, features = c('CD14','FCGR3A','LYZ','ITGAM','ITGAX','PECAM1','PTPRC','CD2', 'NCAM1', 'SELL', 'CSF1R', 'CCR2', 'CX3CR1', 'CXCR3', 'CXCR4'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_Monocytes_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()


  #Classical
  plot<- DotPlot(object, features = c('CD14','FCGR3A','LYZ','PTPRC','CD36'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_MonocytesClassic_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()


  #NonClassic
  plot<- DotPlot(object, features = c('CD14','FCGR3A','MS4A7','ITGAX'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_MonocytesNonClassic_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()

  #Intermediate
  plot<- DotPlot(object, features = c('CD14','FCGR3A','LYZ','IFI30','CD74','CD86','HLA-DRA','S100A8','S100A10'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_MonocytesInter_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()


  ####Macrophague: CD14, CD16 (FCGR3A), CD64 (FCGR1A), CD68, CD71 (TFRC) and CCR5
  plot<- DotPlot(object, features = c('CD14', 'FCGR3A', 'FCGR1A', 'CD68', 'TFRC', 'CCR5'),cols = cols ,group.by = group.by) + RotatedAxis()


  png(file = paste(filePreffix,"_Macrophague_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()


  #M1 : CD86, CD80, CD68, MHCII, IL1R, TLR2, TLR4, iNOS (NOS2), SOCS3
  plot<- DotPlot(object, features = c('CD14', 'FCGR3A', 'CD86', 'CD80', 'TLR2', 'TLR4', 'NOS2','SOCS3', 'HLA-DRA', 'HLA-DRB1'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_MacrophagueM1_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()



  #M2a: CD163, MHCII, SR, MMR/CD206 (MRC1), CD200R (CD200R1), TGM2, IL1R II
  plot<- DotPlot(object, features = c('CD14', 'FCGR3A', 'CD163', 'HLA-DRA', 'HLA-DRB1',  'MRC1', 'CD200R1', 'TGM2', 'IL1R2'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_MacrophagueM2a_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()



  ####Dendritic cells: absence of CD3 (T cell), CD14 (Monocyte), CD19 (B cell), CD56 (NCAM1, NK cells) and CD66b (CEACAM8, granulocytes)
  plot<- DotPlot(object, features = c('CD3D', 'CD14', 'CD19', 'NCAM1'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_DCnegativeMarkers_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()

  #cDC: BDCA-1 (CD1C), CD8, CD8alpha, CD11b (ITGAM), CD11c (ITGAX), CD103 (ITGAE), CD205 (LY75), MHCII (HLA-DRA), IRF8 high, IRF4 low,
  plot<- DotPlot(object, features = c('ID2', 'IRF8','BATF3','ZEB2','IRF4','CD1C', 'CD8A', 'ITGAM', 'ITGAX', 'ITGAE', 'LY75', 'HLA-DRA'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_DCc_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()


  #pDC: BDCA-2 (CLEC4C), CD11clow (ITGAX), CD45RA (PTPRC), CD123 (IL3RA), ILT-7 (LILRA4), MHCIIlow, TLR7, TLR9, IRF8 high, IRF4 high
  plot<- DotPlot(object, features = c( 'TCF4','CLEC4C', 'PTPRC', 'NRP1','IL3RA','IRF8','IRF4','ZEB2', 'LILRA4', 'TLR7', 'TLR9', 'ITGAX','HLA-DRA'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_DCp_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()




  ####NK cells: KLRB1, CD56 (NCAM1), CD62L (SELL), CD57 (B3GAT1), CD94/NKG2A (KLRC1), CD16 (FCGR3A), NKR-P1 (CD161, KLRB1), NKG2D (KLRK1, CD314), NCR
  plot<- DotPlot(object, features = c('KLRB1', 'NCAM1', 'SELL', 'B3GAT1', 'FCGR3A',  'KLRK1', 'NCR1', 'NCR2', 'NCR3'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_NK_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()


  ####NK cells: eddie markers FCGR3A TYROBP CLIC3 NKY7 PRF1 TRBC2(neg) TRAC(neg)
  plot<- DotPlot(object, features = c('FCGR3A', 'TYROBP', 'CLIC3', 'NKY7', 'PRF1','TRBC2','TRAC'), cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_NK_markerGenes_2.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()



  ####Granulocytes
  plot<- DotPlot(object, features = c('CEACAM8', 'LTF','MME','ITGAM','CD63','ENPP3'),cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_Granulocytes_markerGenes.png", sep=""),res = 300, width=3000,height=2000)
  print(plot)
  dev.off()


}










#' Dotplot of marker genes. Creates PNG plots
#'
#'
#' @param object Seurat object
#' @param filePreffix Preffix to use for the file name, including the path to the saving folder
#' @param group.by Name of the metadata column to group for.
#' @param cols Colors to scale the plot. Vector of two colors.
#' @export
GenesDotplot<- function(object, filePreffix, genes, group.by=NULL, cols=c("lightgrey", "#cc0000"),width=3000,height=2000 ){
  ####T cells: CD3D, CD4, CD8A, CD8B
  plot<- DotPlot(object, features = genes, cols = cols ,group.by = group.by) + RotatedAxis()

  png(file = paste(filePreffix,"_markerGenes.png", sep=""),res = 300, width=width, height=height)
  print(plot)
  dev.off()
}








#' Dimensional reduction plot for UMAP coloring the expression of marker genes given by the user. Creates PNG plots only for the indicated gene.
#'
#'
#' @param object Seurat object
#' @param genes Genes to plot
#' @param pt.size Dot size
#' @param order Boolean determining whether to plot cells in order of expression. Can be useful if cells expressing given feature are getting buried.
#' @param width width
#' @param height height
#' @return PNG file. UMAP plot showing the indicated. marker genes
#' @export
FeaturePlot_saver<- function(object,genes=NULL,filePreffix,pt.size=NULL,order=TRUE,width=5000,height=3000){
  plot<- FeaturePlot(object = object, features = genes,  min.cutoff = "q9",pt.size=pt.size,order=order)

  png(file = paste(filePreffix,"_markerGenes.png", sep=""),res = 300, width=width, height=height)
  print(plot)
  dev.off()
}













#' Plot the resolution umaps. Resolutions required 0.1, 0.2, 0.3, 0.4, 0.5
#'
#'
#' @param object Seurat object
#' @param clusters.names list of names for each cluster in order
#' @param clusters.column metadata clumn containing the clusters to considere. Must have numbers.
#' @param output if an output with the cell information, including cluster name and number is desire please indicate output file
#' @return PNG file. UMAP plots showing the clustering of the different resolutions.
#' @export
resolutionPlots<- function(object,reduction='umap',output.image,integrated=F){

  if(integrated){
    plot1 <- DimPlot(object = object, reduction = "umap", group.by = "integrated_snn_res.0.1",)
    plot2 <- DimPlot(object = object, reduction = "umap", group.by = "integrated_snn_res.0.2")
    plot3 <- DimPlot(object = object, reduction = "umap", group.by = "integrated_snn_res.0.3")
    plot4 <- DimPlot(object = object, reduction = "umap", group.by = "integrated_snn_res.0.4")
    plot5 <- DimPlot(object = object, reduction = "umap", group.by = "integrated_snn_res.0.5")

  } else{
    plot1 <- DimPlot(object = object, reduction = "umap", group.by = "SCT_snn_res.0.1",)
    plot2 <- DimPlot(object = object, reduction = "umap", group.by = "SCT_snn_res.0.2")
    plot3 <- DimPlot(object = object, reduction = "umap", group.by = "SCT_snn_res.0.3")
    plot4 <- DimPlot(object = object, reduction = "umap", group.by = "SCT_snn_res.0.4")
    plot5 <- DimPlot(object = object, reduction = "umap", group.by = "SCT_snn_res.0.5")
  }

  ggsave(CombinePlots(plots = list(plot1, plot2, plot3, plot4,plot5)),file = output.image, units = "cm",dpi = 300, width=30,height=20)
}








#' Plot the umaps coloring the mitochondrial DNA %.
#'
#'
#' @param object Seurat object
#' @param output_image file name of the image
#' @param width image width
#' @param height image height
#' @return PNG file. UMAP plot.
#' @export
percent_mitoPlot<- function(object, output_image,width=22, height=20){
  # Determine metrics to plot present in seurat_control@meta.data
  metrics <-  c("percent.mito")

  # Extract the UMAP coordinates for each cell and include information about the metrics to plot
  qc_data <- FetchData(object, vars = c(metrics, "UMAP_1", "UMAP_2"))


  # Plot a UMAP plot for each metric
  p<- ggplot(qc_data, aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color="percent.mito"), alpha = 0.7,size=0.15) +
    scale_color_gradient(low = "grey90", high = "blue")  +
    ggtitle("percent.mito")+ theme_bw()+ theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5))

  ggsave(p,file = output_image,units = "cm",dpi = 300, width=width,height=height)

}




#' Plot the umaps coloring the desired metric as a gradient
#'
#'
#' @param object Seurat object
#' @param metrics metric in metadata to plot with a color gradient
#' @param output_image file name of the image
#' @param width image width
#' @param height image height
#' @return PNG file. UMAP plot.
#' @export
metadataMetrics_gradientPlot<- function(object,metrics='nFeature_RNA', output_image,width=22, height=20, pt.size=0.15, alpha=0.7){
  # Determine metrics to plot present in seurat_control@meta.data

  # Extract the UMAP coordinates for each cell and include information about the metrics to plot
  qc_data <- FetchData(object, vars = c(metrics, "UMAP_1", "UMAP_2"))


  # Plot a UMAP plot for each metric
  p<- ggplot(qc_data, aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=metrics), alpha = alpha,size=pt.size) +
    scale_color_gradient(low = "grey90", high = "blue")  +
    ggtitle(metrics)+ theme_bw()+ theme(panel.grid = element_blank(),plot.title = element_text(hjust=0.5))

  ggsave(p,file = output_image,units = "cm",dpi = 300, width=width,height=height)

}










#' Plot the violin plot of the distribution per group (y) of a desired metric in the metadata
#'
#'
#' @param object Seurat object
#' @param x Groups to divide the violin plot. Violing use this color.
#' @param y metric in metadata to plot distribution
#' @param output_image file name of the image
#' @param width image width
#' @param height image height
#' @return PNG file. violin plots.
#' @export
metadataMetrics_violin<- function(object,x,y,output_image,width=17, height=15){
  plot<- ggplot(object@meta.data, aes(x=get(x), y=get(y),fill=get(x))) +
    geom_violin()+ geom_jitter(position=position_jitter(0.2),size=0.1)+ theme_minimal()+
    labs(fill=x,x =x, y = y)

  ggsave(plot,file = output_image, units = "cm",dpi = 700, width=width,height=height)
}




