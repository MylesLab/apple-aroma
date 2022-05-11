#######################
## TABLE OF CONTENTS ##
#######################

## 0. IMPORTS & DATA LOADING
## 1. PCA HELPER FUNCTIONS
##    1.1. MAIN PCA ANALYSIS

###############################
## 0. IMPORTS & DATA LOADING ##
###############################

#############################
## 1. PCA HELPER FUNCTIONS ##
#############################

## 1.1 MAIN PCA ANALYSIS

#' Given a n x m matrix where there are n samples and m variables, this function
#' runs a principal functions anlaysis on the matrix and returns eigenvectors
#' representing the variance as well as other statistics of the PCA.
#' 
#' @param mat - The n x m matrix containg the data on which to run the PCA.
#'
#' @return results -- the results object contining various statistics
#'
#' @examples
#' # Basic example
#' run_pca(mat)
#'
run_pca <- function(mat) {

  utopian_mat <- mat[, c(1:100)]

  # scale the matrix
  pca_matrix <- scale(utopian_mat)

  # perform PCA on the matrix 
  eig <- stats::prcomp(pca_matrix, center = TRUE)

  #bi-plot


  min(eig$rotation[,1])

  scores <- eig$x
  
  ncol(eig$rotation)
  
  plot(scores[,1],scores[,2])
  
  
  # eigen_values = (eig$sdev)^2
  # 
  # percent_variance = ( eigen_values / sum(eigen_values) ) * 100
  # 
  # # the eigen_values represent the proportion of variance explained by the 
  # # principal functions.
  # 
  # eigen_vectors = eig$rotation
  # # these eigen_vectors (a.k.a loadings) determine the direction of 
  # # the new feature space which accounts for variance in the original matrix.
  # 
  # scores = eig$x
  # # the coordinates of individuals on principal functions
  # 
  # pc1_scores = scores[,1]
  # pc2_scores = scores[,1]
  # pc3_scores = scores[,1]
  # 
  # plot(pc1_scores,pc2_scores, xlab="PC1", ylab="PC2")
  # plot(pc1_scores,pc3_scores, xlab="PC1", ylab="PC3")
  # plot(pc2_scores,pc3_scores, xlab="PC1", ylab="PC3")
  # 
  return(list(
    eig=eig, 
    eigen_values=eigen_values, 
    percent_variance=percent_variance, 
    eigen_vectors=eigen_vectors,
    scores=scores
  ))
}
