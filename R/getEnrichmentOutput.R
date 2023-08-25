#' Getting the enrichment output of target gene set
#' @description get the enrichment output of target gene set by comparing the weighted deviation accumulation of target gene set with randomly sampling gene set
#' @param target_gama the weighted deviation accumulation dataframe of target gene set
#' @param random_gama the list containing 10,000 weighted deviation accumulation dataframe of randomly sampling gene set
#' @return A dataframe containing the enrichment output of target gene set and associated statistics (pval, fdr, bonferroni, etc.)
#' @export
#' @examples enrichment_output <- getEnrichmentOutput(target_gama, random_gama)

getEnrichmentOutput <- function(target_gama,random_gama){
  probability <- data.frame()
  for (i in 1:ncol(target_gama)) {
    geneset_clusteri = target_gama[1,i]
    a=0
    for (j in 1:length(random_gama)) {
      sample_clusteri = random_gama[[j]][1,i]
      if (geneset_clusteri > sample_clusteri) {
        a=a+1
      }
    }
    if (a == length(random_gama)) {
      a = a-0.1
    }
    probability[i,1] <- a
  }
  probability <- probability %>% dplyr::mutate(P = V1/10000)
  probability <- probability %>% dplyr::mutate(pval = 1-P)
  probability <- probability %>% dplyr::mutate(cluster = c(0:(ncol(target_gama)-1)))
  data.table::setorder(probability, cols = 'pval')
  probability <- probability %>% dplyr::mutate(qval = p.adjust(pval,'fdr'))
  probability <- probability %>% dplyr::mutate(bonferroni = p.adjust(pval,'bonferroni'))
  #probability <- probability %>% dplyr::mutate(graphdata = -log10(pval))
  probability <- probability %>% dplyr::mutate(graphdata = -log10(bonferroni))
  return(probability)
}
