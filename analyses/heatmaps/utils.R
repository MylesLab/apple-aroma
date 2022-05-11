# Title     : Utility Functions
# Objective : This script contains all the functions used by this module


#' Takes an nxm matrix with n features and m samples and performs pair-wise correlation
#' 
#' @param cormat the data matrix to run correlation on
#' 
#' @return a list with two variables r and pval which are nxn matrices of r-values
#'         and p-values respectively
perform_matrix_correlation <- function(cormat) {
  cormat.cor.r <- matrix(,ncol(cormat),ncol(cormat))
  cormat.cor.pval <- matrix(,ncol(cormat),ncol(cormat))
  for(col in 1:ncol(cormat)) {
    for(row in 1:ncol(cormat)){
      test <- cor.test(
        as.numeric(unlist(cormat[,col])),
        as.numeric(unlist(cormat[,row])),
        method = "pearson"
      )
      cormat.cor.r[row,col] <- test$estimate
      cormat.cor.pval[row,col] <- test$p.value
    }
  }
  colnames(cormat.cor.r) <- colnames(cormat)
  rownames(cormat.cor.r) <- colnames(cormat)
  colnames(cormat.cor.pval) <- colnames(cormat)
  rownames(cormat.cor.pval) <- colnames(cormat)
  
  return(
    list(
      r=cormat.cor.r,
      pval=cormat.cor.pval
    )
  )
}
