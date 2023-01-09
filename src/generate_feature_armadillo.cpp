#include "fmri.pipeline.h"

//' This function creates K shifts of a neural events vector according to the kernel length, K.
//'
//' @name generate_feature_armadillo
//' @param encoding The neural events vector (same length as BOLD time series)
//' @param K The length of the kernel
//' @return A matrix of length(encoding) rows and K columns, where each column contains a successively
//'    lagged copy of the encoding vector
//'
//' @details This is an internal function that is used inside a while loop by deconvolve_nlreg.
//'   Profiling of the algorithm revealed that this is the primary bottleneck, so I ported it to
//'   an Rcpp function
//' @author Michael Hallquist

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat generate_feature_armadillo(const arma::vec& encoding, int K) {

  //Rcpp::NumericMatrix generate_feature_armadillo(const arma::vec& encoding, int K) {
  int n = encoding.n_elem;
  arma::mat fmatrix(n, K, fill::zeros); //preallocate matrix

  //fmatrix.col(0) = encoding; //unshifted variant placed in first column
  arma::vec first_col(fmatrix.colptr(0), n, false); //in-memory variant of preceding line
  first_col = encoding;

  for (int i=1; i < K; i++) {
    //these don't pan out well because they don't return a vector type
    //double* shift_col = fmatrix.colptr(i); //pointer to column
    //const uword* shift_col = fmatrix.colptr(i); //pointer to column

    //use a pointer to the memory of the column to avoid allocating and copying in memory
    //the first argument is the column pointer plus the positional offset
    //the second argument is the length of the vector (how many elements of the column)
    arma::vec shift_col(fmatrix.colptr(i) + i, n-i, false);
    shift_col = encoding.elem(linspace<uvec>(0, n-i-1, n-i)); //by setting this, these elements of fmatrix will be changed directly

    //The alternative is to allocate a vector, then copy the relevant elements of encoding, respecting the offset
    //We then copy this back to fmatrix below
    //arma::vec shift_col = zeros<vec>(n);
    //shift_col.elem(linspace<uvec>(i,n - 1, n - i)) = encoding.elem(linspace<uvec>(0, n-i-1, n-i));

    //debugging that the prototype works
    //shift_col.elem(linspace<uvec>(1,4,4)) = encoding.elem(linspace<uvec>(0,3,4));

    //need to copy the new column back into the matrix
    //fmatrix.col(i) = shift_col;

  }

  //return wrap(fmatrix); //if numeric matrix
  return fmatrix;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
generate_feature_armadillo(1:100, 10)
*/
