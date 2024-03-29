% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_helper_functions.R
\name{get_mr_abspath}
\alias{get_mr_abspath}
\title{small helper function to pull absolute paths to a given column in run_data or trial_data}
\usage{
get_mr_abspath(mr_df, col = "run_nifti")
}
\arguments{
\item{mr_df}{a data.frame from a gpa object, following standard variable mapping nomenclature}

\item{col}{a character string denoting the column in \code{mr_df} to be used for looking up
absolute paths}
}
\value{
Absolute paths to all files in the specified \code{col}
}
\description{
small helper function to pull absolute paths to a given column in run_data or trial_data
}
\details{
Note that if a given value of the requested column is an absolute path, it will
  not be combined with $mr_dir to generate the combined path. Thus, the \code{col} in
  \code{mr_df} can contain a mixture of relative an absolute paths. The relative paths will
  be combined with
}
\keyword{internal}
