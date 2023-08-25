#' Getting expression data grouped by clusters
#' @description This function divides the expression matrix into different matrixes according to clusters
#' @param SeuratObject a Seurat object of spatial transcriptomics
#' @param ExpMatrix a matrix containing the expression values of all genes in all spots
#' @return A list containing expression matrixes divided by clusters
#' @export
#' @examples exp <- getExpMatrix(mouse_brain)
#' clustered_exp <- getClusteredExp(mouse_brain,exp)

getClusteredExp <- function(SeuratObject,ExpMatrix){
  cluster_cell <- list()
  for (i in 1:length(levels(SeuratObject@active.ident))) {
    n = i-1
    cell_name = names(SeuratObject@active.ident[SeuratObject@active.ident==n])
    cluster_cell[[i]] <- ExpMatrix[,colnames(ExpMatrix) %in% cell_name]
  }
 return(cluster_cell)
}
