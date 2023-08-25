#' filtering out genes which expression levels are too low
#' @description Filter out genes which mean expression levels are lower than threshhold
#' @param clusteredMeanExp a dataframe containing mean expression of all genes in each clusters, rows are genes, columns are clusters
#' @param threshhold Limits for filtering genes.Default is 0.2
#' @return A dataframe containing mean expression of genes which expression levels are hingher than threshhold
#' @export
#' @examples clustered_mean_exp_filtered <- genefilter(clustered_mean_exp, threshhold = 0.2)

genefilter <- function(clusteredMeanExp,threshhold = 0.2){
  keep.list <- rownames(dplyr::filter(data.frame(max = apply(clusteredMeanExp, 1, max)), max >= threshhold))
  clusteredMeanExp <- clusteredMeanExp[rownames(clusteredMeanExp) %in% keep.list,]
  return(clusteredMeanExp)
}
