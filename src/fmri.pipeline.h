#ifndef _fmripipeline_FMRIPIPELINE_h
#define _fmripipeline_FMRIPIPELINE_h

#define ARMA_NO_DEBUG
#include <RcppArmadillo.h>
//#include <omp.h>
using namespace Rcpp;
using namespace arma;

//function definitions

arma::vec convolve_cpp(const arma::vec &a, const arma::vec &b);
arma::vec convolve_double_gamma(const arma::vec &stimulus, double a1, double a2, double b1, double b2, double cc);
//Rcpp::NumericMatrix generate_feature_armadillo(const arma::vec& encoding, int K);
arma::mat generate_feature_armadillo(const arma::vec& encoding, int K);
Rcpp::NumericMatrix generate_feature(Rcpp::NumericVector encoding, int K);

arma::mat deconvolve_nlreg(arma::mat BOLDobs, const arma::vec& kernel,
                             double nev_lr, double epsilon, double beta,
                             bool normalize, bool trim_kernel);

arma::vec sigmoid(const arma::vec& x, double beta);
arma::vec dsigmoid(const arma::vec& x, double beta);
double sigmoid(double x, double beta);
double dsigmoid(double x, double beta);


//RcppExport arma::vec convolve_cpp(arma::vec a, arma::vec b);
//extern "C" arma::vec convolve_cpp(arma::vec a, arma::vec b);

#endif
