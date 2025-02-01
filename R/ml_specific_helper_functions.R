# ---------------Calculate residual-------------------
#' Residual
#'
#' @description
#' Calculates the residuals given observed values, beta coefficient and a design matrix.
#'
#'
#' @param outcomes vector of outcomes.
#' @param design_matrix design matrix.
#' @param beta vector of beta coefficients.
#'
#' @returns Residuals.
#' @export
residual_function <- function(outcomes, design_matrix, beta){
  return(outcomes - design_matrix %*% beta)
}




# -------------Calculate ML S matrix----------------- (27.22)

#' S matrix
#'
#' @description
#' Calculates the S matrix for the ML fisher scoring algorithm.
#'
#'
#' @param semi_def_matrix list of covariance component matrices
#' @param omega_inv the inverse of omega as calculated by omega_func()
#'
#' @returns The S matrix.
#' @export
S_matrix_function <- function(semi_def_matrix, omega_inv){

  # Multiply each covariance component matrix by the inverse of the covariance matrix
  A <- multiply_list_by_matrix(omega_inv, semi_def_matrix)

  # Construct empty matrix to fill with values
  S <- matrix(data = NA, nrow = length(semi_def_matrix), ncol = length(semi_def_matrix))

  for (i in 1:length(semi_def_matrix)){
    for (j in i:length(semi_def_matrix)){
      S[i,j] <- 0.5 * sum(A[[i]] * A[[j]])
      S[j,i] <- S[i,j]
    }
  }

  return(S)

}




# --------------Py matrix----------------- (27.17c)

#' Py matrix
#'
#' @param omega_inv the inverse of omega as calculated by omega_func().
#' @param residual_vec vector of residuals as calculated by residual_function().
#'
#' @returns Py matrix
#' @export
Py_func <- function(omega_inv, residual_vec){
  return(omega_inv %*% residual_vec)
}


#-------------Calculate ML scores-------------------
#' ML Score function
#'
#' @description
#' Function to calculate the maximum likelihood score.
#'
#'
#' @param XtVinv the XtVinv matrix.
#' @param semi_def_matrix list of covariance component matrices.
#' @param omega_inv the inverse of omega as calculated by omega_func().
#' @param Py Py matrix as calculated by Py_func().
#' @param residual_vec vector of residuals as calculated by residual_function().
#'
#' @returns Maximum likelihood score.
#' @export
parameter_score <- function(XtVinv, semi_def_matrix, omega_inv, Py, residual_vec){

  #Fixed effects score (27.10)
  score_beta <- XtVinv %*% residual_vec


  #Random effects score (27.14b, 27.17c)
  score_sigma <- rep(NA, length(semi_def_matrix))

  for (i in 1:length(semi_def_matrix)){
    trVinvVi <- sum(omega_inv * semi_def_matrix[[i]])
    PytViPy <- crossprod(Py, semi_def_matrix[[i]]) %*% Py

    score_sigma[i] <- 0.5 * (PytViPy-trVinvVi)
  }

  score <- c(score_beta, score_sigma)

  return(score)
}
