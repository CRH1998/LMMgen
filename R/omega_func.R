omega_func <- function(semi_def_matrix, sigma2_vec){

  omega <- 0

  for (i in 1:length(semi_def_matrix)){
    omega <- omega + semi_def_matrix[[i]] * sigma2_vec[i]
  }

  return(omega)
}
