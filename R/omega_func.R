#' @title omega_func
#'
#' @description
#' Calculates the variance/covariance matrix Omega or V.
#'
#'
#' @param semi_def_matrix A list containing the three matrices, householder matrix, kinship matrix, identity matrix.
#' @param sigma2_vec A vector containing the variances of the random effects in the same order as the three matrices.
#'
#' @returns The variance/covariance matrix Omega
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





#' Residual Sum of Squares
#'
#' @description
#' Calculates the residual sum of squares.
#'
#' @param X design matrix.
#' @param semi_def_matrix variance/covariance matrix.
#' @param y outcome vector.
#' @param params parameter vector.
#'
#' @returns Residual sum of squares.
#' @export
RSS_func <-  function(X, semi_def_matrix, y, params){

  #Calculating inverse covariance matrix
  V <- omega_func(semi_def_matrix, params)
  Vinv <- chol2inv(chol(V))
  XtVinv <- t(X) %*% Vinv
  XtVinvX <- XtVinv %*% X

  #Calculating betahat
  betahat <- chol2inv(chol(XtVinvX)) %*% (XtVinv %*% y)

  res <- y - X %*% betahat

  #RSS
  RSS <- t(res) %*% Vinv %*% res

  return(RSS)
}



#' XtVinvX
#'
#' @param X design matrix.
#' @param semi_def_matrix variance/covariance matrix.
#' @param y outcome vector.
#' @param params parameter vector.
#'
#' @returns XtVinvX
#' @export
XtVinvX_func <- function(X, semi_def_matrix, y, params){

  #Calculating inverse covariance matrix
  V <- omega_func(semi_def_matrix, params)
  Vinv <- chol2inv(chol(V))
  XtVinv <- t(X) %*% Vinv
  XtVinvX <- XtVinv %*% X

  return(XtVinvX)
}
