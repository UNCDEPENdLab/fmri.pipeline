//#include <Rcpp.h>
//#include <RcppArmadillo.h>
#include "fmri.pipeline.h"

//' Internal port of fsl gammapdf
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec gammapdf(const arma::vec &vals, const float mu, const float var) {
  arma::vec res(vals.n_elem, fill::zeros);
  
  if ((mu > 0) && (var > 0.00001)) {
    float a = std::pow(mu,2)/var;
    float b = mu/var;
    float c = std::lgamma(a);
    //Rcout << "c: " << c << endl;
    if(std::abs(c) < 150){
      for (int mc = 0; mc < res.n_elem; mc++){
        if (vals(mc) > 0.000001) {
          res(mc) = std::exp(a*std::log(b) +
            (a-1) * std::log(vals(mc)) - b*vals(mc) - c);
        }
      }
    }
  }
  
  return res;

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
gammapdf(1:10, mu=5, var=2)
*/

//Rcpp::NumericVector convolve_cpp(Rcpp::NumericVector a, Rcpp::NumericVector b) {
