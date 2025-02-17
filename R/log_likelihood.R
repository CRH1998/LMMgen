
#' Block-wise log-likelihood
#'
#' @param design_matrix fixed effects design matrix for a single block.
#' @param semi_def_matrix list of covariance structure matrices for each block.
#' @param outcomes vector of outcomes for the block.
#' @param fixed_effect_par vector of fixed effect parameters.
#' @param random_effect_par vector of random effect parameters.
#'
#' @returns The log-likelihood for a single block.
#' @export
log_likelihood_block <- function(design_matrix, semi_def_matrix, outcomes, fixed_effect_par, random_effect_par){

  #Determining n_i
  n_i <- nrow(design_matrix)

  #Calculating inverse covariance matrix
  omega <- omega_func(semi_def_matrix, random_effect_par)

  #Inverse omega
  omega_inv <- chol2inv(chol(omega))


  #Calculating log-likelihood
  res <- -n_i/2 * log(2 * pi) - 1/2 * log(det(omega)) - 1/2 * t(outcomes - design_matrix %*% fixed_effect_par) %*% omega_inv %*% (outcomes - design_matrix %*% fixed_effect_par)

  #Returning log-likelihood
  return(res)
}




#--------------------------------------------------------------------------------------------
#       Calculate full log-likelihood using block-wise log-likelihood
#--------------------------------------------------------------------------------------------

#' Log-likelihood
#'
#' @param design_matrices list of fixed effects design matrices - one for each block.
#' @param semi_def_matrices list of lists of covariance structure matrices for each block.
#' @param outcome_list list of vector of outcomes for each block.
#' @param fixed_effect_par vector of fixed effect parameters.
#' @param random_effect_par vector of random effect parameters.
#'
#' @returns The log-likelihood.
#' @export
log_likelihood <- function(design_matrices, semi_def_matrices, outcome_list, fixed_effect_par, random_effect_par){

  #Applying log-likehood function to each element of lists, parameter vector and sigma vector is the same for each individual
  res <- Map(log_likelihood_block, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(fixed_effect_par, random_effect_par))

  return(Reduce('+', res))
}






