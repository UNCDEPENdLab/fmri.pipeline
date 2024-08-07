% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_helper_functions.R
\name{names_to_internal}
\alias{names_to_internal}
\title{helper function to rename columns of input data.frame to internal nomenclature
  based on the variable mapping (vm) vector}
\usage{
names_to_internal(df, vm)
}
\arguments{
\item{df}{a data.frame containing columns to be renamed to internal standards}

\item{vm}{a named vector of columns in \code{df} that identify internal constructs
such as id, session, and run.}
}
\value{
a modified version of \code{df} with column names modified to use internal names
}
\description{
helper function to rename columns of input data.frame to internal nomenclature
  based on the variable mapping (vm) vector
}
\keyword{internal}
