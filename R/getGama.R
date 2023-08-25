#' Calculate the weighted deviation accumulation of gene set
#' @description calculate the weighted deviation accumulation of provided gene set in each clusters
#' @param geneset_delta a dataframe containing delta values of genes in gene set, the last column is the weight of gene
#' @param type provided gene set type. Available options are "target" for taget gene set; "random" for randomly sampling gene set
#' @param cores number of cores used for parallel computation
#' @return A dataframe or a list containing the weighted deviation accumulation of gene set in each clusters
#' @export
#' @examples target_gama <- getGama(target_delta, type = "target")
#' random_gama <- getGama(random_delta, type = "random", cores = 10)

getGama <- function(geneset_delta, type = c('target','random'), cores = n){
  gama <- function(y){
    target_geneset_exp <- data.frame()
    for (i in 1:(ncol(y)-2)) {
      sum = 0
      for(j in 1:nrow(y)){
        g = y[j,i]*y[j,ncol(y)]
        sum = sum + g
      }
      target_geneset_exp[1,i] = sum
      colnames(target_geneset_exp)[i] = colnames(y)[i]
    }
    return(target_geneset_exp)
  }
  if(type == 'target')
    {target_geneset_exp <- gama(geneset_delta)
    return(target_geneset_exp)}
  else if(is.na(cores)==FALSE){
      no_cores = cores
      cl <- parallel::makeCluster(no_cores)
      sampleset_exp <- parallel::parLapply(cl, geneset_delta, gama)
      parallel::stopCluster(cl)
    }
       else{sampleset_exp <- lapply(geneset_delta, gama)}
    return(sampleset_exp)
}
