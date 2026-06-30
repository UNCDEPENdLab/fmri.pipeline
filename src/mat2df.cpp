#include "fmri.pipeline.h"
#include <algorithm>
#include <numeric>

//' Fast conversion of a 2D numeric matrix to a 3-column data.frame
//'
//' @name mat2df
//' @description Converts a 2D numeric matrix into a 3-column data.frame
//' @details This function is a fast matrix-to-long-data-frame converter for the simple 2D case.
//' @param mat A \code{matrix} to convert to data.frame
//' @param na_zeros Whether near-zero values should be converted to \code{NA}
//' @param varnames Names for the two dimension key columns
//' @param value_name Name for the value column
//' @return A 3-column data.frame with two dimension key columns and one value column
//' @keywords internal
//' @author Michael Hallquist
// [[Rcpp::export]]
DataFrame mat2df(
    NumericMatrix mat,
    bool na_zeros = false,
    CharacterVector varnames = CharacterVector::create("dim1", "dim2"),
    std::string value_name = "value"
) {
  if (varnames.size() != 2) {
    stop("varnames must contain exactly two column names");
  }

  const int nrow = mat.nrow();
  const int ncol = mat.ncol();
  const R_xlen_t ncell = static_cast<R_xlen_t>(nrow) * static_cast<R_xlen_t>(ncol);
  const std::string dim1_name = as<std::string>(varnames[0]);
  const std::string dim2_name = as<std::string>(varnames[1]);
  SEXP rownames_attr = R_NilValue;
  SEXP colnames_attr = R_NilValue;

  SEXP dimnames_attr = mat.attr("dimnames");
  if (!Rf_isNull(dimnames_attr)) {
    List dimnames(dimnames_attr);
    if (dimnames.size() >= 2) {
      rownames_attr = dimnames[0];
      colnames_attr = dimnames[1];
      if (Rf_isNull(rownames_attr)) {
        rownames_attr = R_NilValue;
      }
      if (Rf_isNull(colnames_attr)) {
        colnames_attr = R_NilValue;
      }
    }
  }

  IntegerVector dim1 = no_init(ncell);
  IntegerVector dim2 = no_init(ncell);
  NumericVector value = no_init(ncell);
  std::copy(mat.begin(), mat.end(), value.begin());

  IntegerVector::iterator dim1_begin = dim1.begin();
  IntegerVector::iterator dim2_begin = dim2.begin();
  for (int col = 1; col <= ncol; col++) {
    const R_xlen_t offset = static_cast<R_xlen_t>(col - 1) * nrow;
    std::iota(dim1_begin + offset, dim1_begin + offset + nrow, 1);
    std::fill(dim2_begin + offset, dim2_begin + offset + nrow, col);
  }

  if (na_zeros) {
    for (R_xlen_t i = 0; i < ncell; i++) {
      if (std::abs(value[i]) < 1e-4) {
        value[i] = NA_REAL;
      }
    }
  }

  if (!Rf_isNull(rownames_attr)) {
    CharacterVector rownames(rownames_attr);
    dim1.attr("class") = "factor";
    dim1.attr("levels") = rownames;
  }

  if (!Rf_isNull(colnames_attr)) {
    CharacterVector colnames(colnames_attr);
    dim2.attr("class") = "factor";
    dim2.attr("levels") = colnames;
  }

  return DataFrame::create(
    Named(dim1_name.c_str()) = dim1,
    Named(dim2_name.c_str()) = dim2,
    Named(value_name.c_str()) = value
  );
}
