#define RNIFTI_NIFTILIB_VERSION 2

#include "Rcpp.h"
#include "RNifti.h"
#include "RNiftiAPI.h"

using namespace Rcpp;
using namespace RNifti;

// needed for NiftiImageData dim 
typedef int64_t dim_t;

//' Get Dimensions of a NIfTI Image
//'
//' Reads the header of a NIfTI file and returns its image dimensions.
//'
//' This function uses the RNifti C++ API to efficiently extract the dimensions
//' (e.g., x, y, z, time) of a NIfTI image without loading the entire image into memory.
//'
//' @name get_nifti_dim
//' @param infile Character string. Path to a valid NIfTI file (e.g., `.nii` or `.nii.gz`).
//'
//' @return A numeric vector containing the dimensions of the image. For a 4D image,
//'         the result will be a vector of length 4: \code{c(x, y, z, t)}.
//'
//' @examples
//' \dontrun{
//'   dims <- get_nifti_dim("sub-001_task-rest_bold.nii.gz")
//' }
//' @keywords internal

// [[Rcpp::export]]
NumericVector get_nifti_dim(std::string infile) {
  // Load header and extract dimensions
  NiftiImage image(infile, false); // read nifti header only
  NiftiImageData data(image);      // get pointer to image data
  std::vector<dim_t> dims = image.dim();

  return NumericVector(dims.begin(), dims.end());
}

