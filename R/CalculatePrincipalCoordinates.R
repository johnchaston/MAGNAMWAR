#' Show Principal Components Breakdown
#' 
#' Function to show Principal Components statistics based on the OrthoMCL presence absence groupings.
#' @param mcl_data output of FormatAfterOrtho --list of 2 things-- 1: binary matrix indicating the presence / absence of genes in each OG and 2: vector of names of OGs
#' @return returns a named list of principal components and accompanying proportion of variance for each
#' @examples
#' CalculatePrincipalCoordinates(after_ortho_format)
#' @export

CalculatePrincipalCoordinates <- function(mcl_data){
  pa <- t(mcl_data$pa_matrix)
  
  # turn to numeric
  x <- mapply(pa, FUN=as.numeric)
  
  rows <- dimnames(pa)[[1]]
  cols <- dimnames(pa)[[2]]
  
  m <- matrix(data=x, ncol=length(cols), nrow=length(rows))
  
  rownames(m) <- rows
  colnames(m) <- cols
  
  nm <- m[ , apply(m, 2, var) != 0]
  
  prm <- prcomp(nm, scale=T)
  pr_coor_mtx = prm$x
  return(summary(prm)$importance[2,])
}
