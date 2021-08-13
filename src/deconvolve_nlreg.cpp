#include "fmri.pipeline.h"

//' C++ port of Bush and Cisler 2013, Magnetic Resonance Imaging
//' Adapted from the original provided by Keith Bush
//' as well as C++ code from Jiang Bian
//'
//' @name deconvolve_nlreg
//' @param BOLDobs matrix of observed BOLD timeseries (n_timepoints x n_signals)
//' @param kernel  assumed kernel of the BOLD signal (e.g., from spm_hrf)
//' @param nev_lr  learning rate for the assignment of neural events. Default: .01
//' @param epsilon relative error change (termination condition). Default: .005
//' @param beta slope of the sigmoid transfer function (higher = more nonlinear)
//' @param normalize whether to unit-normalize (z-score) \code{BOLDobs} before deconvolution. Default: TRUE
//' @param trim_kernel whether to remove the first K time points from the deconvolved vector, corresponding to
//'            kernel leftovers from convolution. Default: TRUE
//'
//' @details
//' This function deconvolves the BOLD signal using Bush 2011 method
//'
//' Author:      Keith Bush, PhD
//' Institution: University of Arkansas at Little Rock
//' Date:        Aug. 9, 2013
//'
//' The original code did not unit normalize the BOLD signal in advance, but in my testing, this
//' proves useful in many cases (unless you want to mess with the learning rate a lot), especially
//' when the time series has a non-zero mean (e.g., mean 100).
//'
//' @return A time series of the same length containing reconstructed neural events
//' @author Michael Hallquist
//' @export

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat deconvolve_nlreg(arma::mat BOLDobs, const arma::vec& kernel,
                               double nev_lr=.01, double epsilon=.005, double beta=40,
                               bool normalize=true, bool trim_kernel=true) {

  //since we may need to transform BOLDobs (z-score columns), don't use const + & in argument
  //this is a little slower since we need a deep copy of BOLDobs, but that matrix shouldn't be too
  //big and the copy time should pale in comparison to deconvolution itself

  //mat must be timepoints x signals (e.g., regions)
  int N = BOLDobs.n_rows;
  int P = BOLDobs.n_cols;

  //Calc simulation steps related to simulation time
  int K = kernel.n_elem;
  int A = K - 1 + N;

  //double sd_i;
  //double m_i;

  arma::vec col_means = arma::mean(BOLDobs);
  //cout << "means" << col_means << endl;

  if (normalize == true) {
    //cout << "normalizing data" << endl;
    // for (int i=0; i < P; i++) {
    //   arma::vec col_i(BOLDobs.colptr(i), N, false); //get memory reference to ith column
    //   sd_i = arma::stddev(col_i);
    //   m_i = arma::mean(col_i);
    // }

    arma::vec col_sds = arma::stddev(BOLDobs);
    //cout << "sds" << col_sds << endl;

    BOLDobs.each_row() -= col_means; //subtract means
    BOLDobs.each_row() /= col_sds; //divide by sds
  } else {
    //always demean columns
    BOLDobs.each_row() -= col_means; //subtract means
  }

  //arma::vec col_i = Y.row(i);

  //Construct initial random activation vector (fluctuate slightly around zero between -2e-9 and 2e-9)
  //This can be done once and re-used for each time series
  //arma::vec uu(A, fill::ones); //fill with ones
  //arma::vec rr(A, fill::randu); //random uniform initialization
  //arma::vec activation = (2e-9*uu)*rr - 1e-9; //uu contains fixed 2e-9, rr is rand unif number
  arma::vec init_activation = arma::randu<arma::vec>(A) * (2e-9) - (1e-9);

  int max_hrf_id_adjust = kernel.index_max() - 1; // element of kernel 1 before max

  //arma::mat dEde(K, K, fill::zeros);
  //dEde.diag().ones(); //make identity matrix
  // shouldn't just use the kernel? why need to get identity matrix?
  arma::mat dEde = arma::eye<arma::mat>(K, K) * kernel;

  arma::mat decon_results(A, P, fill::zeros); //result

  // loop over time series (signals)
  //#pragma omp parallel for
  for (int i=0; i < P; i++) {
    arma::vec thisBOLD(BOLDobs.colptr(i), N, false); //get memory reference to ith column
    arma::vec activation = init_activation; //start with initial (random) activation vector

    //Termination Params (reset for each time series)
    double preverror = 1e9; //previous error
    double currerror = 0;   //current error

    //Presolve activations to fit target_adjust as encoding
    //arma::vec BOLDobs_adjust = BOLDobs.submat(max_hrf_id_adjust, i, N, i); //BOLDobs_adjust = BOLDobs[max_hrf_id_adjust:N];
    arma::vec BOLDobs_adjust = thisBOLD.elem(regspace<arma::uvec>(max_hrf_id_adjust, N-1)); //BOLDobs_adjust = BOLDobs[max_hrf_id_adjust:N];
    arma::vec pre_encoding = BOLDobs_adjust - min(BOLDobs_adjust);
    pre_encoding = pre_encoding/arma::max(pre_encoding); // unit normalize
    arma::vec encoding = pre_encoding;

    //update elements of activation
    activation.elem(regspace<arma::uvec>(K, K-1+BOLDobs_adjust.n_elem) - 1) = (1/beta)*log(pre_encoding/(1-pre_encoding));

    //preallocate vectors
    arma::vec ytilde(N, fill::zeros);
    arma::vec brf(N, fill::zeros);
    arma::vec dEdbrf(N, fill::zeros);
    arma::vec dEdy(N, fill::zeros);

    while (std::abs(preverror-currerror) > epsilon) {

      // Compute encoding vector
      encoding = sigmoid(activation, beta);

      // Construct feature space
      //cout << "encoding is vec: " << encoding << ", len:" << encoding.n_elem << ", K is: " << K << endl;
      arma::mat feature = generate_feature_armadillo(encoding, K);

      //Generate virtual bold response by multiplying feature (N x K) by kernel (K x 1) to get N x 1 estimated response
      ytilde = feature.rows(K - 1, feature.n_rows - 1) * kernel;

      //Convert to percent signal change
      double meanCurrent = arma::mean(ytilde);
      brf = (ytilde - meanCurrent)/meanCurrent; //unit normalize

      // Compute dEdbrf
      dEdbrf = brf - thisBOLD;

      // Assume normalization does not impact deriv much.
      dEdy = dEdbrf;

      // Precompute derivative components
      arma::vec back_error = arma::zeros(dEdy.n_elem + 2 * (K-1));
      back_error(regspace<arma::uvec>(K-1, K-2+dEdy.n_elem)) = dEdy; //0 padding on left and right tails

      // Backpropagate Errors
      arma::vec delta(A, fill::zeros);

      for (int j = 0; j < A; j++) {
        double active = activation(j);
        double deda = dsigmoid(active, beta);
        arma::vec dEda = dEde * deda;
        arma::vec this_error = back_error(regspace<arma::uvec>(j, (j-1+K)));
        delta(j) = arma::sum(dEda % this_error);
      }

      // Update estimate
      activation = activation - nev_lr * delta;

      // Iterate Learning
      preverror = currerror;
      currerror = arma::sum(dEdbrf%dEdbrf);
    }

    decon_results.col(i) = encoding;
  }

  if (trim_kernel == true) {
    // remove the initial timepoints corresponding to HRF (so that returned signal matches in time and length)
    decon_results.shed_rows(0, K-2);
  }

  return(decon_results);
}
