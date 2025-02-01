#-------------------------------------------
#         List helper functions
#-------------------------------------------


#' Title
#'
#' @param outcome list of list of outcomes.
#'
#' @returns the mean of the list of lists of outcomes.
#' @export
#'
get_outcome_mean <- function(outcome){
  return(mean(as.vector(sapply(outcome, function(x){return(x)}))))
}

#' Title
#'
#' @param outcome list of list of outcomes.
#'
#' @importFrom stats var
#'
#' @returns the variance of the list of lists of outcomes.
#' @export
#'
get_outcome_variance <- function(outcome){
  return(var(as.vector(sapply(outcome, function(x){return(x)}))))
}
