#' Calculating the mean expression of all genes in each clusters
#' @description Calculate the mean expression of each genes in each clusers separately
#' @param clusteredExp a list containing expression matrixes divided by clusters
#' @return A dataframe containing mean expression of all genes in each clusters
#' @export
#' @examples clustered_exp <- getClusteredExp(mouse_brain,exp)
#' clustered_mean_exp <- Clustered_mean_Exp(clustered_exp)

Clustered_mean_Exp <- function(clusteredExp){
  gene_name = rownames(clusteredExp[[1]])
  gene_cluster <- data.frame()
  for (i in 1:length(gene_name)) {
    for (j in 1:length(clusteredExp)) {
     gene_cluster[i,j] <- sum(clusteredExp[[j]][gene_name[i],])/ncol(clusteredExp[[j]])
    }
    rownames(gene_cluster)[i] <- gene_name[i]
  }
  colnames(gene_cluster) <- paste0('cluster', c(1:length(clusteredExp)))
  return(gene_cluster)
}
