#' Rescale calculation
#' @description rescale the delta values
#' @param delta the delta dataframe which to use
#' @return A dataframe containing rescale values
#' @export
#' @examples rescale_gene <- rescale(delta)

rescale = function(delta){
  res = function(y){
    rng = range(y, na.rm = T)
    z = (y-rng[1])/(rng[2]-rng[1])
    return(z)
  }
  rescale_gene = apply(delta, 1, res)
  return(rescale_gene)
}
