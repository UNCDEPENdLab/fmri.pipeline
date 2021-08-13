#include "fmri.pipeline.h"

//#include <Rcpp.h>
//using namespace Rcpp;

//' This function creates K shifts of a neural events vector according to the kernel length, K.
//'
//' @name generate_feature
//' @param encoding The neural events vector (same length as BOLD time series)
//' @param K The length of the kernel
//' @return A matrix of length(encoding) rows and K columns, where each column contains a successively
//'    lagged copy of the encoding vector
//'
//' @details This is an internal function that is used inside a while loop by deconvolve_nlreg.
//'   Profiling of the algorithm revealed that this is the primary bottleneck, so I ported it to
//'   an Rcpp function
//' @author Michael Hallquist

// [[Rcpp::export]]
NumericMatrix generate_feature(NumericVector encoding, int K) {
  int n = encoding.size();
  NumericMatrix fmatrix(n, K); //preallocate matrix

  fmatrix(_,0) = encoding; //unshifted variant placed in first column

  for (int i=1; i < K; i++) {
    NumericVector shift_col(n); //populates with zeros anyhow, so we just need to overwrite relevant elements

    //Fill in shifted encoding for this column in a nested loop.
    //This is about 30% slower than the vectorized approach below
    //But on the ICS clusters, the vectorized approach segfaults regularly (not totally clear why)
    for (int j=0; j < n - i; j++) {
      shift_col(j+i) = encoding(j);
    }

    //Rcpp quietly refuses to handle subvectors on LHS and RHS
    //shift_col[Rcpp::seq(i, n)] = encoding[Rcpp::seq(0, n - i)];

    //Use of a temporary variable works
    //NumericVector shift_part = encoding[Rcpp::seq(0, n - i)];
    //shift_col[Rcpp::seq(i, n)] = shift_part;

    fmatrix(_,i) = shift_col; //assign shifted column back to matrix

  }

  return fmatrix;
}

//some misguided attempts to use submatrices (to no avail)
//NumericVector zpad(i-1, 0); // padding
//NumericMatrix::Sub zpad_ref = fmatrix(Rcpp::Range(0, i-1), Rcpp::Range(i, i)); //grab a reference to the first i-1 row of the ith column
//NumericMatrix zpad(i-1,1); //grab a reference to the first i-1 row of the ith column
//zpad_ref = zpad;
//NumericMatrix::Sub zpad = fmatrix(Rcpp::Range(0, 1), Rcpp::Range(0, 2)); //grab a reference to the first i-1 row of the ith column

//zpad = 0;

//NumericMatrix::Column shift = fmatrix(Rcpp::Range(0, shift.size() + zpad.size()), Rcpp::Range(i,i));
//shift = encoding[Range(0, n - i - 1)];

//fmatrix(Rcpp::Range(0, zpad.size()),i) = zpad;
//fmatrix(Rcpp::Range(0, shift.size() + zpad.size()), i) = shift;

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
generate_feature_cpp(1:100, 10)
*/
