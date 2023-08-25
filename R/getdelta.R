#' Calculating the delta values of all genes
#' @description Calculate the delta values (differences of one gene from the average expression in cluster to the overall average expression) of all genes, two-sides delta, high expression delta, low expression delta can be calculated
#' @param cluster_mean_exp_filtered a filtered clustered gene mean expression dataframe
#' @param ExpMatrix a matrix containing the expression values of all genes in all spots
#' @param type calculation type which to use. Available options are "two sides":calculating delta values in all clusters; "high expression":only calculating delta values of in clusters where mean expressions are above overall averages; "low expression":only calculating delta values of in clusters where mean expressions are below overall averages
#' @return A dataframe containing delta values of all genes
#' @export
#' @examples delta <- getdelta(cluster_mean_exp_filtered, exp, type = "two sides")

getdelta <- function(cluster_mean_exp_filtered,ExpMatrix,type = c('two sides','high expression','low expression')){
  ExpMatrix = ExpMatrix[rownames(ExpMatrix) %in% rownames(cluster_mean_exp_filtered),]
  mean_gene <- data.frame(mean_expression = apply(ExpMatrix, 1, mean))
  delta <- data.frame()
  if(type == 'two sides'){
   for (i in 1:nrow(cluster_mean_exp_filtered)) {
     delta <- rbind(delta, abs(cluster_mean_exp_filtered[i,]-mean_gene[i,1]))
   }
  }
  if(type == 'high expression'){
    for (i in 1:nrow(cluster_mean_exp_filtered)) {
      for (j in 1:ncol(cluster_mean_exp_filtered)) {
        if(cluster_mean_exp_filtered[i,j]>mean_gene[i,1]){
          delta[i,j] = cluster_mean_exp_filtered[i,j]-mean_gene[i,1]
        }
        else{
          delta[i,j] = 0
        }
      }
    }
  }
  if(type == 'low expression'){
    for (i in 1:nrow(cluster_mean_exp_filtered)) {
      for (j in 1:ncol(cluster_mean_exp_filtered)) {
        if(cluster_mean_exp_filtered[i,j]<mean_gene[i,1]){
          delta[i,j] = mean_gene[i,1]-cluster_mean_exp_filtered[i,j]
        }
        else{
          delta[i,j] = 0
        }
      }
    }
  }
  rownames(delta) <- rownames(cluster_mean_exp_filtered)
  colnames(delta) <- colnames(cluster_mean_exp_filtered)
  return(delta)
}
