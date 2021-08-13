#include "fmri.pipeline.h"

//' Dsigmoid transform
//'
//' @name dsigmoid
//' @param x value to be transformed
//' @param beta slope (steepness) of sigmoid transform
//' @keywords internal

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec dsigmoid(const arma::vec& x, double beta=1) {
  arma::vec y = (1 - sigmoid(x, beta)) * sigmoid(x, beta);
  return(y);
}

double dsigmoid(double x, double beta=1) {
  return (1.0 - sigmoid(x, beta)) * sigmoid(x, beta);
}

/*** R
dsigmoid(seq(-10, 10), .5)
*/

