//#include <Rcpp.h>
//#include <RcppArmadillo.h>
#include "fmri.pipeline.h"

//' Internal port of fsl do_convolve
//'
//' @name do_convolve
//' @param input the input
//' @param kernel the kernel to convolve
//' @param phase an integer for phase-shifting the HRF
//' @param renorm boolean indicating whether to renormalize the output by the sum
//' @return A vector containing the convolution of input and kernel
//'
//' @details This is an internal function used for testing alignment with FSL HRF convolution

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec do_convolve(const arma::vec &input, const arma::vec &kernel, int phase = 0, const int renorm = true) {
  int ts_len = input.n_elem;
  int kern_len = kernel.n_elem;
  arma::vec output(ts_len, fill::zeros);
  //phase=5;
  double kernel_norm;
  
  for (int t = 0; t < ts_len; t++) {
    kernel_norm = 0.0;
    //int b_i = std::max(0, t + phase - ts_len + 1);
    //int e_i = std::min(kern_len, t + phase + 1);
    //Rcout << "start: " << b_i << ", end: " << e_i << endl;
    for (int i = std::max(0, t + phase - ts_len + 1); i < std::min(kern_len, t + phase + 1); i++) {
      output(t) += input(t+phase-i) * kernel(i); // also remove +1?
      kernel_norm += kernel(i);
    }
    
    if (renorm)
      output(t) /= kernel_norm;
  }
  
  return output;
}

/*** R
v1 <- rnorm(100)
v2 <- rnorm(30)
v3 <- rnorm(100)
r1_fsl <- do_convolve(v1, v2)
r1_r <- convolve(v1, v2, type = "open")
r2_fsl <- do_convolve(v1, v3)
r2_r <- convolve(v1, v3)
*/

//Rcpp::NumericVector convolve_cpp(Rcpp::NumericVector a, Rcpp::NumericVector b) {
