#' Percentage calcultion
#' @description calculate the expression percentages of all genes in each clusters
#' @param rescale a rescale dataframe which to use
#' @return A dataframe containing expression percentages of all genes in each clusters
#' @export
#' @examples percentage_gene <- percentage(rescale_gene)

percentage <- function(rescale){
  y = data.frame()
  rescale = t(rescale)
  for (i in 1:nrow(rescale) ){
    for (j in 1:ncol(rescale)) {
      y[i,j] = rescale[i,j]/sum(rescale[i,])
    }
  }
  colnames(y) = colnames(rescale)
  rownames(y) = rownames(rescale)
  return(y)
}
