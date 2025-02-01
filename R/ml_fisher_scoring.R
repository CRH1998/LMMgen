
#####################################################################
#                                                                   #
#     This script contains the algorithm for running the            #
#     ML fisher-scoring algorithm in a block-multivariate           #
#     normal gaussian where the covariance matrix is block-         #
#     diagonal and a linear combination of known semi-              #
#     definite matrices                                             #
#                                                                   #
#####################################################################









#-------------------------------------------
#       ML fisher scoring function
#-------------------------------------------
#' @title Fisher scoring algorithm
#'
#' @description
#' Computes the relevant parts of the fisher scoring algorithm for a single block.
#'
#' @param design_matrix fixed effects design matrix for a single block.
#' @param semi_def_matrix list of covariance structure matrices for each block.
#' @param outcomes vector of outcomes for the block.
#' @param params vector of initial parameters.
#' @param small_value_threshold values in matrices below this threshold is set to 0.
#' @param add_small_constant adds this constant to diagonal of matrices with small values for stability.
#'
#' @returns List containing scores and fisher information.
#' @export
score_fisher_function <- function(design_matrix, semi_def_matrix, outcomes, params, small_value_threshold = 1e-12, add_small_constant = 1e-9){

  #-------------Parameters-----------------
  beta <- params[1:ncol(design_matrix)]                             #extracting mean-value parameters
  sigma2_vec <- params[(ncol(design_matrix) + 1):length(params)]    #extracting variance parameters


  # Calculating omega inverse
  omega <- omega_func(semi_def_matrix, sigma2_vec)

  # Setting very small values to 0
  omega[omega < small_value_threshold] <- 0

  # Adding small value to diagonal if diagonal values are very small
  omega <- omega + (diag(omega) < small_value_threshold) * add_small_constant * diag(length(diag(omega)))

  # Inverting omega
  omega_inv <- chol2inv(chol(omega))

  #Calculating residual vector
  residual_vec <- residual_function(outcomes = outcomes, design_matrix = design_matrix, beta = beta)


  #--------------Calculating Py matrix-----------------
  Py <- Py_func(omega_inv = omega_inv, residual_vec = residual_vec)



  #-------------Calculating mean value parameters----- (27.21)
  XtVinv <- crossprod(design_matrix, omega_inv)

  M <- XtVinv %*% design_matrix



  #-------------Calculating S matrix-----------------
  S <- S_matrix_function(semi_def_matrix = semi_def_matrix, omega_inv = omega_inv)


  #-------------Calculating scores-------------------


  score <- parameter_score(XtVinv = XtVinv,
                           omega_inv = omega_inv,
                           semi_def_matrix = semi_def_matrix,
                           Py = Py,
                           residual_vec = residual_vec)


  return(list('M' = M, 'S' = S, 'score' = score))
}






#-------------------------------------------
#       ML fisher scoring algorithm
#-------------------------------------------
#' Find MLE parameters
#'
#' @description
#' Finds the MLE parameters using the fisher scoring algorithm.
#'
#'
#' @param init_params vector of initial parameters, first the fixed effects and then the random effects.
#' @param design_matrices list of design matrices for fixed effects for each block.
#' @param semi_def_matrices list of lists of variance/covariance matrices.
#' @param outcome_list list of vectors of outcomes for each block.
#' @param update_step_size update steps are multiplied with this value.
#' @param max_iter an integer determining the maximum number of iterations.
#' @param tolerance numeric value to determine convergence tolerance.
#' @param small_value_threshold values in matrices below this threshold is set to 0.
#' @param add_small_constant adds this constant to diagonal of matrices with small values for stability.
#'
#' @returns The parameter estimates, the inverse fisher information, the score and a summary.
#' @export
find_mle_parameters <- function(init_params, design_matrices, semi_def_matrices, outcome_list, update_step_size = 1, max_iter = 10000, tolerance = 1e-12,
                                small_value_threshold = 1e-9, add_small_constant = 1e-12){

  for (iter in 1:max_iter) {

    out <- Map(score_fisher_function, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(init_params))


    # Sum blocks
    M_sum <- 0
    for (i in 1:length(out)){
      M_sum <- M_sum + out[[i]]$M
    }

    S_sum <- 0
    for (i in 1:length(out)){
      S_sum <- S_sum + out[[i]]$S
    }

    #M_sum <- Reduce('+',lapply(out, function(x) x$M))
    #S_sum <- Reduce('+',lapply(out, function(x) x$S))


    # Setting very small values to 0
    M_sum[M_sum < small_value_threshold] <- 0
    S_sum[S_sum < small_value_threshold] <- 0

    # Adding small value to diagonal if diagonal values are very small
    M_sum <- M_sum + (diag(M_sum) < small_value_threshold) * add_small_constant * diag(length(diag(M_sum)))
    S_sum <- S_sum + (diag(S_sum) < small_value_threshold) * add_small_constant * diag(length(diag(S_sum)))


    # Define inverse fisher information
    fisher_inv <- Matrix::bdiag(chol2inv(chol(M_sum)), chol2inv(chol(S_sum)))

    # Sum scores
    score <- 0
    for (i in 1:length(out)){
      score <- score + out[[i]]$score
    }

    #score <- rowSums(sapply(out, function(x) x$score))

    #Calculate update step
    update_step <- fisher_inv %*% score


    # Check for convergence
    if (sum((update_step)^2) < tolerance) {
      break
    }

    # Update parameters for the next iteration
    init_params <- init_params + update_step_size * update_step
    #init_params[init_params < 0] <- 0.001
  }


  #--------------------Summary--------------------
  n_fixed_effects <- ncol(design_matrices[[1]])
  n_random_effects <- length(semi_def_matrices[[1]])
  n_effects <- list('Number of fixed effects' = n_fixed_effects, 'Number of random effects' = n_random_effects)


  #Summary dataframe
  summary <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(summary) <- c("Parameter", "Estimate", "SE", "Lower_CI_bound", "Upper_CI_bound")

  for (i in 1:length(init_params)){

    if (i <= n_fixed_effects){
      Parameter <- paste0('beta_',i-1)
    } else {
      Parameter <- paste0('sigma_',i-n_fixed_effects-1)
    }
    Estimate <- init_params[i]
    SE <- sqrt(fisher_inv[i,i])
    Lower_CI_bound <- Estimate - 1.96*SE
    Upper_CI_bound <- Estimate + 1.96*SE

    summary[nrow(summary) + 1,] <- c(Parameter, Estimate, SE, Lower_CI_bound, Upper_CI_bound)
  }


  return(list('Estimates' = init_params, 'Inverse fisher' = fisher_inv, 'Score' = score, 'Summary' = summary, 'Number of effects' = n_effects))
}




