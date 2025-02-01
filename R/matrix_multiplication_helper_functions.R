#-------------------------------------------
#   Matrix multiplication helper functions
#-------------------------------------------


# Multiply list by matrix
multiply_list_by_matrix <- function(matrix, list){

  matrix_mult <- function(matrix1){
    return(function(matrix2){matrix1 %*% matrix2})
  }

  matrix_mult_matrix <- matrix_mult(matrix)

  return(lapply(list, matrix_mult_matrix))
}



# Multiply list by matrix and take trace of each matrix product
tr_multiply_list_by_matrix <- function(matrix, list){

  tr_matrix_mult <- function(matrix1){
    return(function(matrix2){sum(matrix1 * matrix2)})
  }

  tr_matrix_mult_matrix <- tr_matrix_mult(matrix)

  return(lapply(list, tr_matrix_mult_matrix))
}



# Multiply list by matrix and specify if multiplication is on the right or the left side
matrix_mult_specify <- function(matrix1, mult_by_right = F){

  if(mult_by_right == F){
    return(function(matrix2){matrix1 %*% matrix2})
  } else {
    return(function(matrix2){matrix2 %*% matrix1})
  }
}

matrix_mult_list_by_matrix <- function(matrix, list, mult_by_right = F){
  matrix_mult_matrix <- matrix_mult_specify(matrix, mult_by_right)
  return(lapply(list, matrix_mult_matrix))
}
