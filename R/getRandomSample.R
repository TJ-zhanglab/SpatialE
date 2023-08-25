#' Getting the randomly sampling gene set and corresponding delta matrix
#' @description get 10,000 randomly sampling gene set with same length of target gene set and corresponding delta matrix
#' @param delta the delta dataframe of all genes
#' @param target_geneset the target gene set which to test
#' @return A list containing 10,000 randomly sampling gene set and corresponding delta matrix
#' @export
#' @examples random_delta <- getRandomSample(delta, target_geneset)

getRandomSample <- function(delta,target_geneset){
  genesetlength <- length(target_geneset[target_geneset %in% rownames(delta)])
  sample_set <- list()
  for (i in 1:10000) {
    sample_set[[i]] = sample(rownames(delta), genesetlength, replace = F)
    names(sample_set)[i] = paste0('sample',i)
  }
  sample_geneset_exp <- list()
  for (i in 1:10000) {
    sample_geneset_exp[[i]] = delta[rownames(delta) %in% sample_set[[i]],]
    names(sample_geneset_exp)[i] = paste0('sample',i)
  }
  return(sample_geneset_exp)
}
