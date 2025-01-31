#' Title
#'
#' @param semi_def_matrix A list containing the three matrices, householder matrix, kinship matrix, identity matrix.
#' @param sigma2_vec A vector containing the variances of the random effects in the same order as the three matrices.
#'
#' @returns The variance/covariance matrix Sigma
#' @export
#'
#' @examples
#' mat1 <- matrix(1:9, ncol = 3)
#' mat2 <- matrix(10:18, ncol = 3)
#' mat3 <- matrix(19:27, ncol = 3)
#' omega_func(semi_def_matrix = list(mat1, mat2, mat3), sigma2_vec = c(1,2,3))
omega_func <- function(semi_def_matrix, sigma2_vec){

  omega <- 0

  for (i in 1:length(semi_def_matrix)){
    omega <- omega + semi_def_matrix[[i]] * sigma2_vec[i]
  }

  return(omega)
}
