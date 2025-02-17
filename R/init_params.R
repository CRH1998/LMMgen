
#' Initialize parameters
#'
#' @param design_matrices list of design matrices
#' @param outcome_list list of outcomes
#'
#' @returns initial parameters for the fisher scoring algorithm
#' @export
init_params <- function(design_matrices, outcome_list){

  beta_estimates <- mapply(function(X, y) chol2inv(chol(t(X) %*% X)) %*% t(X) %*% y,
                           design_matrices, outcome_list, SIMPLIFY = FALSE)
  beta_matrix <- do.call(cbind, beta_estimates)
  beta_init <- rowMeans(beta_matrix)

  resid <- mapply(residual_function, outcome_list, design_matrices, beta_estimates, SIMPLIFY = F)
  sigma_R_init <- mean(unlist(resid)^2)
  sigma_H_init <- 1
  sigma_K_init <- 1

  return(list(init_fixed_effect_par = beta_init, init_random_effect_par = c(sigma_R_init, sigma_H_init, sigma_K_init)))
}
