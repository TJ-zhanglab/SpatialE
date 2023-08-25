#' Getting expression matrix
#' @description Get the expression matrix from Seurat object of spatial transcriptomics
#' @param SeuratObject a Seurat object of spatial transcriptomics
#' @return A matrix containing the expression values of all genes in all spots
#' @export
#' @examples
#' data("mouse_brain")
#' mouse_brain <- SCTransform(mouse_brain, assay = "Spatial") %>%
#'                RunPCA(assay = "SCT", verbose = FALSE) %>%
#'                FindNeighbors(reduction = "pca", dims = 1:30) %>%
#'                FindClusters(resolution = 0.8) %>%
#'                RunUMAP(reduction = "pca", dims = 1:30)
#' exp <- getExpMatrix(mouse_brain)

getExpMatrix <- function(SeuratObject){
  exp_matrix <- as.matrix(SeuratObject@assays$SCT@counts)
  return(exp_matrix)
}
