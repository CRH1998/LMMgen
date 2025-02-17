#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List sum_blocks(List out) {
  // Extract the first element to determine matrix dimensions
  List first_element = out[0];
  NumericMatrix M_sum = clone(as<NumericMatrix>(first_element["M"])); // Clone for correct dimensions
  NumericMatrix S_sum = clone(as<NumericMatrix>(first_element["S"]));

  // Ensure sum starts at 0
  std::fill(M_sum.begin(), M_sum.end(), 0);
  std::fill(S_sum.begin(), S_sum.end(), 0);

  for (int i = 0; i < out.size(); i++) {
    List element = out[i];

    NumericMatrix M_mat = element["M"];
    NumericMatrix S_mat = element["S"];

    // Ensure matrices are the same size
    if (M_mat.nrow() != M_sum.nrow() || M_mat.ncol() != M_sum.ncol() ||
        S_mat.nrow() != S_sum.nrow() || S_mat.ncol() != S_sum.ncol()) {
      stop("Matrix dimensions are inconsistent across list elements");
    }

    // Element-wise summation
    M_sum += M_mat;
    S_sum += S_mat;
  }

  return List::create(Named("M_sum") = M_sum,
                      Named("S_sum") = S_sum);
}
