#include "fmri.pipeline.h"

//' This function convolves a stimulus vector with the double-gamma hrf
//'
//' @name convolve_double_gamma
//' @param x A vector of volume numbers used to evaluate the function at each value
//' @param a1 The a1 parameter of the double gamma
//' @param a2 The a2 parameter of the double gamma
//' @param b1 The b1 parameter of the double gamma
//' @param b2 The b2 parameter of the double gamma
//' @param cc The cc parameter of the double gamma
//' @return A vector of the double-gamma HRF at each value of \code{x}
//'
//' @details This is an internal function that is used by convolve_hrf
//' @author Michael Hallquist

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec convolve_double_gamma(const arma::vec& stimulus, double a1=6.0, double a2=12.0,
                           double b1=0.9, double b2=0.9, double cc=0.35) {
  double d1 = a1 * b1; //time to peak
  double d2 = a2 * b2; //time to undershoot

  arma::vec x = regspace(1, stimulus.n_elem); //1-based seq_along stimulus
  //arma::vec c1 = pow(x/d1, a1);
  //arma::vec c2 = cc * pow(x/d2, a2);
  arma::vec res = pow(x/d1, a1) % exp(-(x-d1)/b1) - (cc * pow(x/d2, a2)) % exp(-(x-d2)/b2);

  //Rcpp::Rcout << "c1 is: " << c1 << endl;
  //Rcpp::Rcout << "d1 is: " << d1 << endl;
  //Rcpp::Rcout << "d2 is: " << d2 << endl;

  //this is the equivalent of 'open' in R convolve()
  //arma::vec cvector = conv(stimulus, res);
  arma::vec cvector = convolve_cpp(stimulus, res);

  //need to subset down to just the first n_elem values post-convolution
  arma::vec central = join_cols(zeros(1), cvector.elem(regspace<uvec>(0, stimulus.n_elem - 2)));

  return(central);
}

/*** R
convolve_double_gamma(seq(50, 100, 10))
*/

//const arma::vec& x,
