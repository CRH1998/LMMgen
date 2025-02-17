
#' @title pair_parents_with_children
#' @description pairs parents with unassigned individuals
#'
#' @param parents Data frame. Data frame with four columns containing dadid, sex, momid, sex
#' @param unassigned_individuals Data frame. Data frame containing id and sex
#'
#' @returns list containing the current family dataframe (sub_family_df) and the remaining unassigned individuals (unassigned_individuals)
#' @export

pair_parents_with_children <- function(parents, unassigned_individuals){

  sub_family_df <- data.frame()

  for (i in 1:nrow(parents)){

    # Assign potential children to parents
    no_children <- rbinom(1,nrow(unassigned_individuals),0.25) #Number of children to be sampled from unassigned individuals
    children <- unassigned_individuals[sample(nrow(unassigned_individuals), no_children),] #Unassigned individuals sampled as children


    # Construct sub_family dataframe
    ifelse(no_children == 0, #If the number of sampled children is 0, the sub_family dataframe is empty
           sub_family_df <- rbind(sub_family_df, data.frame(matrix(ncol = 4, nrow = 0))),
           sub_family_df <- rbind(sub_family_df, data.frame(children, parents[i,1], parents[i,3]))
    )

    #Remove children who were assigned in previous statement
    unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% children$id),]
  }

  names(sub_family_df) <- c("id", "sex","dadid", "momid")
  return(list('sub_family_df' = sub_family_df, 'unassigned_individuals' = unassigned_individuals))
}





# Construct function that pairs children with spouses

#' @title pair_children_with_spouses
#' @description pairs children with unassigned individuals (spouses)
#'
#' @param family_df Data frame. Data frame containing four columns: id, sex, momid, dadid
#' @param unassigned_individuals Data frame. Data frame containing id and sex
#'
#' @returns List. List containing parents data frame, the remaining unassigned individuals data frame and assigned individuals data frame
#' @export

pair_children_with_spouses <- function(family_df, unassigned_individuals){

  # Get sex of children
  male_children <- family_df[family_df$sex == 1,]
  female_children <- family_df[family_df$sex == 2,]


  # Get ids of men and women in the id column
  unassigned_men <- unassigned_individuals[unassigned_individuals$sex == 1,]
  unassigned_women <- unassigned_individuals[unassigned_individuals$sex == 2,]


  # Sample spouses for children among the remaining
  no_male_spouses <- sample(0:min(nrow(female_children), nrow(unassigned_men)), 1)
  no_female_spouses <- sample(0:min(nrow(male_children), nrow(unassigned_women)), 1)

  married_male_unassigned <- unassigned_men[sample(nrow(unassigned_men), no_male_spouses),]
  married_female_unassigned <- unassigned_women[sample(nrow(unassigned_women), no_female_spouses),]


  # Sample children
  married_male_children <- male_children[sample(nrow(male_children), no_female_spouses),][,c("id", "sex")]
  married_female_children <- female_children[sample(nrow(female_children), no_male_spouses),][,c("id", "sex")]


  # Construct parent dataframe
  parents <- data.frame(rbind(cbind(married_male_children, married_female_unassigned),
                              cbind(married_male_unassigned, married_female_children)))


  # Update assigned individuals
  ifelse(no_male_spouses == 0 & no_female_spouses == 0,
         assigned_individuals <- data.frame(matrix(ncol = 4, nrow = 0)),
         assigned_individuals <- data.frame(rbind(married_male_unassigned, married_female_unassigned), 0, 0))


  names(assigned_individuals) <- c("id", "sex","dadid", "momid")


  # Update unassigned individuals
  unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% assigned_individuals$id),]

  return(list('parents' = parents, 'unassigned_individuals' = unassigned_individuals, 'assigned_individuals' = assigned_individuals))
}






# Construct function that generates kinship matrix

#' @title generate_kinship_matrix
#' @description function that generate kinship matrix
#'
#' @param n_individuals Integer. Number of individuals in given cluster
#' @param seed NA or integer. If left unspecified no seed is set. Otherwise set a seed.
#'
#' @returns Matrix. Returns a random kinship matrix of dimensions n_individuals x n_individuals
#' @export

generate_kinship_matrix <- function(n_individuals, seed = NA){

  if (!is.na(seed)){
    set.seed(seed)
  }


  # Specify id of each individual in the family
  id  <- 1:n_individuals

  # Randomly assign sex to each individual in the family
  sex <- sample(c(1,2), n_individuals, prob = c(1/2,1/2), replace = T)

  # Unassigned individuals
  unassigned_individuals <- data.frame(id,sex)

  # Get ids of men and women in the id column
  men <- unassigned_individuals[unassigned_individuals$sex == 1,]
  women <- unassigned_individuals[unassigned_individuals$sex == 2,]

  # Pair random man with random woman to create first parents
  first_dad <- men[sample(nrow(men), 1),]
  first_mom <- women[sample(nrow(women), 1),]

  parents <- data.frame(first_dad, first_mom)

  # Remove first parents from unassigned individuals
  unassigned_individuals <- unassigned_individuals[!(unassigned_individuals$id %in% c(first_dad$id, first_mom$id)),]

  full_family_df <- data.frame(rbind(first_dad, first_mom), 0, 0)
  names(full_family_df) <- c("id", "sex","dadid", "momid")

  while(nrow(unassigned_individuals) > 0){

    sub_family <- pair_parents_with_children(parents, unassigned_individuals)

    # Define sub_family_df
    sub_family_df <- sub_family$sub_family_df

    # Update unassigned individuals
    unassigned_individuals <- sub_family$unassigned_individuals

    # Assign spouses to children
    parents_list <- pair_children_with_spouses(sub_family_df, unassigned_individuals)

    # Update potential parents
    if (nrow(parents_list$parents) > 0){
      parents <- parents_list$parents
    }

    # Update assigned individuals, that is unassigned individuals who were married to one of the children
    assigned_individuals <- parents_list$assigned_individuals
    unassigned_individuals <- parents_list$unassigned_individuals
    full_family_df <- rbind(full_family_df, sub_family_df, assigned_individuals)

  }

  full_family_df <- full_family_df %>% dplyr::arrange(id)

  kinship_matrix <- kinship2::kinship(id = full_family_df$id,
                            dadid = full_family_df$dadid,
                            momid = full_family_df$momid,
                            sex = full_family_df$sex)

  return(kinship_matrix*2) #Multiply by 2 to get diagonal of ones
}








#-----------------------------------------------
#               Generate family data
#-----------------------------------------------


# Generate data set for testing

#' @title Family data set generator
#' @description Generate genetics data.
#'
#' @param n_clusters an integer denoting the number of clusters.
#' @param n_individuals an integer denoting the number of individuals in given cluster.
#' @param variance_param a vector of characters taking values 'R', 'H' and/or 'K' specifying covariance structure.
#' @param n_mean_param an integer specifying number of fixed effects. Taking values 1,2,3.
#' @param beta_0 a real denoting the first fixed effect.
#' @param beta_1 a real denoting the second fixed effect.
#' @param beta_2 a real denoting the third fixed effect.
#' @param sigma_R a real denoting the residual variance.
#' @param sigma_H a real denoting the household variance.
#' @param sigma_K a real denoting the kinship variance.
#' @param seed NA or integer. If left unspecified no seed is set. Otherwise set a seed.
#'
#' @importFrom stats rbinom rnorm
#'
#' @returns List. Returns a list consisting of a list of outcomes, list of gamma-matrices and list of design matrices
#' @export

family_dataset_generator <- function(n_clusters = 100, n_individuals = 10, variance_param = c('R', 'H', 'K'), n_mean_param = 1,
                                     beta_0 = 1, beta_1 = 3, beta_2 = 3, sigma_R = 1, sigma_H = 3, sigma_K = 5, seed = NA){

  # Set seed if specified for sampling later on
  if (!is.na(seed)){
    set.seed(seed)
  }



  #------------------------Generate covariance structure------------------------

  sigma2_vec <- c(sigma_R, sigma_H, sigma_K)


  if(("R" %in% variance_param) & ("H" %in% variance_param) & ("K" %in% variance_param)){

    # R matrix
    residual_matrix <- as.matrix(diag(1, nrow = n_individuals))

    # H matrix
    household_matrix <- matrix(1, nrow = n_individuals, ncol = n_individuals)

    # K matrix and repeat R and H matrix
    semi_def_matrices <- replicate(n_clusters,
                                   list(residual_matrix, household_matrix, generate_kinship_matrix(n_individuals = n_individuals)),
                                   simplify = F)

  } else if (("R" %in% variance_param) & ("H" %in% variance_param)){

    # R matrix
    residual_matrix <- as.matrix(diag(1, nrow = n_individuals))

    # H matrix
    household_matrix <- matrix(1, nrow = n_individuals, ncol = n_individuals)

    # K matrix and repeat R and H matrix
    semi_def_matrices <- replicate(n_clusters,
                                   list(residual_matrix, household_matrix),
                                   simplify = F)

  } else if (("R" %in% variance_param) & ("K" %in% variance_param)){

    # R matrix
    residual_matrix <- as.matrix(diag(1, nrow = n_individuals))

    # K matrix and repeat R matrix
    semi_def_matrices <- replicate(n_clusters,
                                   list(residual_matrix, generate_kinship_matrix(n_individuals = n_individuals)),
                                   simplify = F)

  } else if (("H" %in% variance_param) & ("K" %in% variance_param)){

    # H matrix
    household_matrix <- matrix(1, nrow = n_individuals, ncol = n_individuals)

    # K matrix and repeat H matrix
    semi_def_matrices <- replicate(n_clusters,
                                   list(household_matrix, generate_kinship_matrix(n_individuals = n_individuals)),
                                   simplify = F)

  } else if ("R" %in% variance_param){

    # R matrix
    residual_matrix <- as.matrix(diag(1, nrow = n_individuals))

    # Repeat R matrix
    semi_def_matrices <- replicate(n_clusters,
                                   list(residual_matrix),
                                   simplify = F)

  } else if ("H" %in% variance_param){

    # H matrix
    household_matrix <- matrix(1, nrow = n_individuals, ncol = n_individuals)

    # Repeat H matrix
    semi_def_matrices <- replicate(n_clusters,
                                   list(household_matrix),
                                   simplify = F)

  } else {

    # K matrix and repeat
    semi_def_matrices <- replicate(n_clusters,
                                   list(generate_kinship_matrix(n_individuals = n_individuals)),
                                   simplify = F)
  }


  #-----------------------------------------------------------------------------




  #-----------------------------Mean value structure----------------------------
  # Defining the design_matrices
  if (n_mean_param == 1){
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals))), simplify = F)
  } else if (n_mean_param == 2){
    #In this case the design matrices are an intercept and one covariate x1 which is drawn from a normal distribution with mean 0 and sd = 1.
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals),
                                                             rnorm(n_individuals))), simplify = F)
  } else if (n_mean_param == 3) {
    #In this case the design matrices are an intercept and one covariate x1 which is drawn from a normal distribution with mean 0 and sd = 1.
    design_matrices <- replicate(n_clusters, as.matrix(cbind(rep(1, n_individuals),
                                                             rnorm(n_individuals),
                                                             rnorm(n_individuals))), simplify = F)
  }
  #-----------------------------------------------------------------------------




  #-------------------Generate data based specifications above------------------


  # Sampling n_clusters outcomes from a multivariate normal distribution with
  # specified mean and covariance structure with n_individuals in each cluster

  # Adding mean-value parameters if specified

  outcome_list <- list()

  if (n_mean_param == 1){
    for (i in 1:n_clusters){
      outcome_list[[i]] <- MASS::mvrnorm(n = 1, mu = rep(beta_0, n_individuals), Sigma = omega_func(semi_def_matrices[[i]], sigma2_vec = sigma2_vec))
    }
  } else if (n_mean_param == 2){
    for (i in 1:n_clusters){
      outcome_list[[i]] <- MASS::mvrnorm(n = 1, mu = rep(beta_0, n_individuals) + beta_1 * design_matrices[[i]][,2],
                                   Sigma = omega_func(semi_def_matrices[[i]], sigma2_vec = sigma2_vec))
    }
  } else if (n_mean_param == 3){
    for (i in 1:n_clusters){
      outcome_list[[i]] <- MASS::mvrnorm(n = 1, mu = rep(beta_0, n_individuals) + beta_1 * design_matrices[[i]][,2] + beta_2 * design_matrices[[i]][,3],
                                   Sigma = omega_func(semi_def_matrices[[i]], sigma2_vec = sigma2_vec))
    }
  }

  return(list('design_matrices' = design_matrices, 'semi_def_matrices' = semi_def_matrices, 'outcome_list' = outcome_list))
}








#-------------------------------FDG object--------------------------------------



#' @title Family data set generator object
#' @description Object for generating genetics data.
#'
#' @param n_clusters an integer denoting the number of clusters.
#' @param n_individuals an integer denoting the number of individuals in given cluster.
#' @param variance_param a vector of characters taking values 'R', 'H' and/or 'K' specifying covariance structure.
#' @param n_mean_param an integer specifying number of fixed effects. Taking values 1,2,3.
#' @param beta_0 a real denoting the first fixed effect.
#' @param beta_1 a real denoting the second fixed effect.
#' @param beta_2 a real denoting the third fixed effect.
#' @param sigma_R a real denoting the residual variance.
#' @param sigma_H a real denoting the household variance.
#' @param sigma_K a real denoting the kinship variance.
#' @param seed NA or integer. If left unspecified no seed is set. Otherwise set a seed.
#'
#' @returns Family data set
#' @export
FDG <- function(n_clusters = 100, n_individuals = 10,
               variance_param = c('R', 'H', 'K'),
               n_mean_param = 1, beta_0 = 1, beta_1 = 3, beta_2 = 3,
               sigma_R = 1, sigma_H = 3, sigma_K = 5, seed = NA) {

  family_dataset <- family_dataset_generator(n_clusters = n_clusters,
                                             n_individuals = n_individuals,
                                             variance_param = variance_param,
                                             n_mean_param = n_mean_param,
                                             beta_0 = beta_0, beta_1 = beta_1, beta_2 = beta_2,
                                             sigma_R = sigma_R, sigma_H = sigma_H, sigma_K = sigma_K,
                                             seed = seed)

  outcome_mean <- get_outcome_mean(family_dataset$outcome_list)
  outcome_var <- get_outcome_variance(family_dataset$outcome_list)



  output <- structure(
    list(
      design_matrices = family_dataset$design_matrices,
      semi_def_matrices = family_dataset$semi_def_matrices,
      outcome_list = family_dataset$outcome_list,
      n_clusters = n_clusters,
      n_individuals = n_individuals,
      variance_param = variance_param,
      sigma_R = sigma_R, sigma_H = sigma_H, sigma_K = sigma_K,
      outcome_mean = outcome_mean,
      outcome_var = outcome_var),
    class = "FDG"
  )
  return(output)
}


#' Family data set generator
#'
#' @param x an object of class FDG.
#' @param ... additional arguments
#'
#' @returns Print for object of class FDG.
#' @export
#'
print.FDG <- function(x, ...){

  cat(paste0("Object of class ", class(x)), "\n\n")
  cat(paste0("Number of clusters: ", x$n_clusters, " with ", x$n_individuals, " individuals in each cluster\n\n"))


  cat("Mean and variance of simulated outcomes:\n")
  # Column headers
  cat(paste0("               Mean      Variance\n"))
  # Formatting coefficients
  cat(sprintf("%-12s %8.4f %10.4f\n", "y", round(x$outcome_mean, 3), round(x$outcome_var, 3)),"\n")


  if(("R" %in% x$variance_param) & ("H" %in% x$variance_param) & ("K" %in% x$variance_param)){

    cat("Variance parameters used: \nResidual (R), Household (H) and Kinship (K)\n\n")
    cat(paste0("Values used: ", "\nsigmaR: ", x$sigma_R, "  sigmaH: ", x$sigma_H, "  sigmaK: ", x$sigma_K))

  } else if (("R" %in% x$variance_param) & ("H" %in% x$variance_param)){

    cat("Variance parameters used: \nResidual (R) and Household (H)\n\n")
    cat(paste0("Values used: ", "\nsigmaR: ", x$sigma_R, "  sigmaH: ", x$sigma_H))

  } else if (("R" %in% x$variance_param) & ("K" %in% x$variance_param)){

    cat("Variance parameters used: \nResidual (R) and Kinship (K)\n\n")
    cat(paste0("Values used: ", "\nsigmaR: ", x$sigma_R, "  sigmaK: ", x$sigma_K))

  } else if (("H" %in% x$variance_param) & ("K" %in% x$variance_param)){

    cat("Variance parameters used: \nHousehold (H) and Kinship (K)\n\n")
    cat(paste0("Values used: ", "\nsigmaH: ", x$sigma_H, "  sigmaK: ", x$sigma_K))

  } else if ("R" %in% x$variance_param){

    cat("Variance parameters used: \nResidual (R)\n\n")
    cat(paste0("Values used: ", "\nsigmaR: ", x$sigma_R))

  } else if ("H" %in% x$variance_param){

    cat("Variance parameters used:\nHousehold (H)\n\n")
    cat(paste0("Values used: ", "\nsigmaH: ", x$sigma_H))

  } else {

    cat("Variance parameters used: \nKinship (K)\n\n")
    cat(paste0("Values used: ", "\nsigmaK: ", x$sigma_K))

  }
}




