#' Calculate the weight of all genes
#' @description calculate the weight of all genes according to the degree of dispersion
#' @param entropy the entropy dataframe which to use
#' @return A dataframe containing weight of all genes
#' @export
#' @examples weight_gene <- weight(entopy_gene)

weight <- function(entropy){
  d = 1-entropy
  w = d/sum(d)
  colnames(w) = 'weight'
  w$gene <- rownames(w)
  return(w)
}
