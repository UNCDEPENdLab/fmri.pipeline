#include "fmri.pipeline.h"

//' Sigmoid transform
//'
//' @name sigmoid
//' @param x value to be transformed by sigmoid
//' @param beta beta slope (steepness) of sigmoid transform
//'
//' @keywords internal

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec sigmoid(const arma::vec& x, double beta=1) {
  arma::vec y = 1/(1+exp(-beta*x));
  return(y);
}

double sigmoid(double x, double beta=1) {
  return 1.0 / (1.0 + std::exp(-beta*x));
}

/*** R
sigmoid(seq(-10, 10), .5)
*/

