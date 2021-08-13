//#include <Rcpp.h>
//#include <RcppArmadillo.h>
#include "fmri.pipeline.h"

//' Internal function to convolve two vectors using nested for loop
//'
//' @name convolve_cpp
//' @param a A vector
//' @param b A vector
//' @return A vector containing the convolution of a and b
//'
//' @details This is an internal function
//'       Same result as, in R:
//'       convolve(a, b, conj=TRUE, type="open")

//' @author Dirk Eddelbuettel

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec convolve_cpp(const arma::vec &a, const arma::vec &b) {
  /* Define the numeric vectors and put in the data from R */
  //Rcpp::NumericVector xa(a);
  //Rcpp::NumericVector xb(b);

  /* Get the sizes of the vectors */
  int n_xa = a.n_elem, n_xb = b.n_elem;

  /* The length of the overlap */
  int nab = n_xa + n_xb - 1;

  arma::vec xab(nab);
  for (int i = 0; i < n_xa; i++)
    for (int j = 0; j < n_xb; j++)
      xab[i + j] += a[i] * b[j];
  return xab;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
convolve_cpp(rnorm(100), rnorm(100))
*/

//Rcpp::NumericVector convolve_cpp(Rcpp::NumericVector a, Rcpp::NumericVector b) {
