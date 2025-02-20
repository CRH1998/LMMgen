
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
#' @param fixed_effect_par vector of initial fixed effect parameters.
#' @param random_effect_par vector of initial random effect parameters.
#' @param small_value_threshold values in matrices below this threshold is set to 0.
#' @param add_small_constant adds this constant to diagonal of matrices with small values for stability.
#'
#' @returns List containing scores and fisher information.
#' @export
score_fisher_function <- function(design_matrix, semi_def_matrix, outcomes, fixed_effect_par, random_effect_par, small_value_threshold = 1e-12, add_small_constant = 1e-9){

  if (ncol(design_matrix) != length(fixed_effect_par)){
    stop("Number of columns in design matrix does not match the prespecified number of fixed effects.")
  }

  if (length(semi_def_matrix) != length(random_effect_par)){
    stop("Number of variance/covariance matrices does not match the number of random effects. Please specify the residual matrix, the householder matrix and the kinship matrix. Check semi_def_matrix and random_effect_par.")
  }


  #-------------Parameters-----------------
  beta <- fixed_effect_par                              # mean-value parameters
  sigma2_vec <- random_effect_par                       # variance parameters


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
#' @param init_fixed_effect_par vector of initial fixed effects parameters. Make sure the order of the initial fixed effects match the order of the columns in the design matrices.
#' @param init_random_effect_par vector of initial random effects parameters. Make sure the order of the initial random effects match the order of the variance/covariance effects in the semi_def_matrices.
#' @param design_matrices list of design matrices for fixed effects for each block.
#' @param semi_def_matrices list of lists of variance/covariance matrices.
#' @param outcome_list list of vectors of outcomes for each block.
#' @param update_step_size update steps are multiplied with this value.
#' @param max_iter an integer determining the maximum number of iterations.
#' @param tolerance numeric value to determine convergence tolerance.
#' @param small_value_threshold values in matrices below this threshold is set to 0.
#' @param add_small_constant adds this constant to diagonal of matrices with small values for stability.
#'
#' @useDynLib LMMgen
#' @importFrom Rcpp sourceCpp
#'
#' @returns The parameter estimates, the inverse fisher information, the score and a summary.
#' @export
find_mle_parameters <- function(init_fixed_effect_par = NULL, init_random_effect_par = NULL, design_matrices, semi_def_matrices, outcome_list, update_step_size = 1, max_iter = 10000, tolerance = 1e-12,
                                small_value_threshold = 1e-9, add_small_constant = 1e-12){

  if (is.null(init_fixed_effect_par)){
    init_fixed_effect_par <- init_params(design_matrices = design_matrices, outcome_list = outcome_list)$init_fixed_effect_par
  }

  if (is.null(init_random_effect_par)){
    init_random_effect_par <- init_params(design_matrices = design_matrices, outcome_list = outcome_list)$init_random_effect_par
  }

  # Count number of iterations
  fisher_scoring_iterations <- 0

  for (iter in 1:max_iter) {

    fisher_scoring_iterations <- fisher_scoring_iterations + 1

    out <- Map(score_fisher_function, design_matrices, semi_def_matrices, outcome_list, MoreArgs = list(init_fixed_effect_par, init_random_effect_par))

    #Sum blocks
    M_sum <- 0
    for (i in 1:length(out)){
      M_sum <- M_sum + out[[i]]$M
    }

    S_sum <- 0
    for (i in 1:length(out)){
      S_sum <- S_sum + out[[i]]$S
    }


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

    #Calculate update step
    update_step <- fisher_inv %*% score

    # Check for convergence
    if (sum((update_step)^2) < tolerance) {
      converged <- TRUE
      break
    }

    # Update parameters for the next iteration
    params <- c(init_fixed_effect_par, init_random_effect_par) + update_step_size * update_step
    init_fixed_effect_par <- params[1:length(init_fixed_effect_par)]
    init_random_effect_par <- params[(length(init_fixed_effect_par)+1):length(params),]
  }

  #--------------------Warnings--------------------
  if (any(init_random_effect_par <= 0)){
    init_random_effect_par[init_random_effect_par < 0] <- 0
    warning("Random effect variance was estimated to be negative and set to 0.")
  }

  if (fisher_scoring_iterations == max_iter){
    converged <- FALSE
    warning("Algorithm did not converge. Coefficients may be misleading!")
  }


  #--------------------Summary--------------------
  n_fixed_effects <- length(init_fixed_effect_par)
  n_random_effects <- length(init_random_effect_par)
  n_effects <- list('Number of fixed effects' = n_fixed_effects, 'Number of random effects' = n_random_effects)


  #Summary dataframe
  summary <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(summary) <- c("Parameter", "Estimate", "SE", "Lower_CI_bound", "Upper_CI_bound")

  for (i in 1:length(init_fixed_effect_par)){

    Parameter <- paste0('beta_',i-1)
    Estimate <- init_fixed_effect_par[i]
    SE <- as.numeric(sqrt(fisher_inv[i,i]))
    Lower_CI_bound <- Estimate - 1.96*SE
    Upper_CI_bound <- Estimate + 1.96*SE

    summary[nrow(summary) + 1,] <- c(Parameter, Estimate, SE, Lower_CI_bound, Upper_CI_bound)
  }

  for (i in 1:length(init_random_effect_par)){

    Parameter <- paste0('sigma_',i-1)
    Estimate <- init_random_effect_par[i]
    SE <- sqrt(fisher_inv[i+n_fixed_effects,i+n_fixed_effects])
    Lower_CI_bound <- Estimate - 1.96*SE
    Upper_CI_bound <- Estimate + 1.96*SE

    summary[nrow(summary) + 1,] <- c(Parameter, Estimate, SE, Lower_CI_bound, Upper_CI_bound)
  }


  return(list('n_fixed_effects' = n_fixed_effects,
              'fixef' = init_fixed_effect_par,
              'n_random_effects' = n_random_effects,
              'ranef' = init_random_effect_par,
              'fisher_inv' = fisher_inv,
              'score' = score,
              'fisher_scoring_iterations' = fisher_scoring_iterations,
              'converged' = converged,
              'summary' = summary))
}





#' Find MLE parameters
#'
#' @description
#' Finds the MLE parameters using the fisher scoring algorithm.
#'
#'
#' @param init_fixed_effect_par vector of initial fixed effects parameters. Make sure the order of the initial fixed effects match the order of the columns in the design matrices.
#' @param init_random_effect_par vector of initial random effects parameters. Make sure the order of the initial random effects match the order of the variance/covariance effects in the semi_def_matrices.
#' @param design_matrices list of design matrices for fixed effects for each block.
#' @param semi_def_matrices list of lists of variance/covariance matrices.
#' @param outcome_list list of vectors of outcomes for each block.
#' @param update_step_size update steps are multiplied with this value.
#' @param max_iter an integer determining the maximum number of iterations.
#' @param tolerance numeric value to determine convergence tolerance.
#' @param small_value_threshold values in matrices below this threshold is set to 0.
#' @param add_small_constant adds this constant to diagonal of matrices with small values for stability.
#'
#' @useDynLib LMMgen
#' @importFrom Rcpp sourceCpp
#'
#' @returns The parameter estimates, the inverse fisher information, the score and a summary.
#' @export
lmmgen <- function(init_fixed_effect_par = NULL, init_random_effect_par = NULL,
                   design_matrices, semi_def_matrices, outcome_list,
                   update_step_size = 1, max_iter = 100, tolerance = 1e-12,
                   small_value_threshold = 1e-9, add_small_constant = 1e-12){



  ml_output <- find_mle_parameters(init_fixed_effect_par = init_fixed_effect_par,
                                   init_random_effect_par = init_random_effect_par,
                                   design_matrices = design_matrices,
                                   semi_def_matrices = semi_def_matrices,
                                   outcome_list = outcome_list,
                                   update_step_size = update_step_size,
                                   max_iter = max_iter,
                                   tolerance = tolerance,
                                   small_value_threshold = small_value_threshold,
                                   add_small_constant = add_small_constant)

  log_lik <- log_likelihood(design_matrices = design_matrices,
                            semi_def_matrices = semi_def_matrices,
                            outcome_list = outcome_list,
                            fixed_effect_par = ml_output$fixef,
                            random_effect_par = ml_output$ranef)

  output <- structure(
    list('n_fixed_effects' = ml_output$n_fixed_effects,
         'fixef' = ml_output$fixef,
         'n_random_effects' = ml_output$n_random_effects,
         'ranef' = ml_output$ranef,
         'fisher_inv' = ml_output$fisher_inv,
         'score' = ml_output$score,
         'fisher_scoring_iterations' = ml_output$fisher_scoring_iterations,
         'converged' = ml_output$converged,
         'summary' = ml_output$summary,
         log_lik = log_lik),
    class = "lmmgen"
  )
  return(output)
}




#' Print lmmgen model likelihood estimates
#'
#' @param x an object of class lmmgen
#' @param ... additional arguments
#'
#' @returns Print for object of class lmmgen
#' @export
#'
print.lmmgen <- function(x, ...){

  cat(paste0("Object of class ", class(x)), "\n\n")

  #browser()
  # Converged
  if(x$converged){
    cat(paste0("Log-likelihood: ", x$log_lik, "\n"))
    cat(paste0("Converged after: ", x$fisher_scoring_iterations, " fisher scoring iterations.\n\n"))
  } else {
    cat(paste0("Log-likelihood: ", x$log_lik, "\n"))
    cat(paste0("Did not converge after: ", x$fisher_scoring_iterations, " fisher scoring iterations.\n\n"))
  }

  cat(paste0("Number of fixed effects parameters: ", x$n_fixed_effects, "\n"))
  cat(paste0("Number of random effects parameters: ", x$n_random_effects, "\n\n"))



  # Print Coefficients
  cat("Fixed effects coefficients:\n")
  fixef <- x$fixef
  fixef_se <- sqrt(abs(x$fisher_inv[1:length(fixef), 1:length(fixef)]))

  # Column headers
  cat(paste0("               Estimate   Std. Error    Lower Wald CI Upper Wald CI\n"))

  # Formatting coefficients
  for (i in 1:length(fixef)) {
    cat(sprintf("%-12s %8.4f %10.4f %13.4f %14.4f\n",
                paste0('beta_',i-1), fixef[i], as.matrix(fixef_se)[i,i], fixef[i] - 1.96 * as.matrix(fixef_se)[i,i], fixef[i] + 1.96 * as.matrix(fixef_se)[i,i]))
  }

  cat("\n")

  cat("Random effects coefficients:\n")
  ranef <- x$ranef
  ranef_se <- sqrt(abs(x$fisher_inv[(length(fixef)+1):ncol(x$fisher_inv), (length(fixef)+1):ncol(x$fisher_inv)]))

  # Column headers
  cat(paste0("               Estimate   Std. Error    Lower Wald CI Upper Wald CI\n"))

  # Formatting coefficients
  for (i in 1:length(ranef)) {
    cat(sprintf("%-12s %8.4f %10.4f %13.4f %14.4f\n",
                paste0('sigma_',i-1), ranef[i], ranef_se[i,i], ranef[i] - 1.96 * ranef_se[i,i], ranef[i] + 1.96 * ranef_se[i,i]))
  }

}







