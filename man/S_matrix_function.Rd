% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_specific_helper_functions.R
\name{S_matrix_function}
\alias{S_matrix_function}
\title{S matrix}
\usage{
S_matrix_function(semi_def_matrix, omega_inv)
}
\arguments{
\item{semi_def_matrix}{list of covariance component matrices}

\item{omega_inv}{the inverse of omega as calculated by omega_func()}
}
\value{
The S matrix.
}
\description{
Calculates the S matrix for the ML fisher scoring algorithm.
}
