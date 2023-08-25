#' Calculate information entropies of all genes
#' @description calculate information entropies of all genes according to the degree of dispersion
#' @param percentage the percentage dataframe which to use
#' @return A dataframe containing information entropies of all genes
#' @export
#' @examples entropy_gene <- entropy(percentage_gene)

entropy <- function(percentage){
  ent <- function(y){
    a = 0
    for (i in 1:length(y)) {
      if (y[i]==0) {
        z = 0
      }
      else{
        z = y[i]*log(y[i])
      }
      a = a + z
    }
    E = -a/log(length(y))
  }
  entropy_gene <- data.frame(apply(percentage, 1, ent))
  colnames(entropy_gene) <- 'entropy'
  return(entropy_gene)
}
