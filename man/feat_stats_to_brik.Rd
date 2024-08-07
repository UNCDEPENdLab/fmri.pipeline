% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsl_helper_functions.R
\name{feat_stats_to_brik}
\alias{feat_stats_to_brik}
\title{internal function to convert a feat level 1 directory to an AFNI brik/head file}
\usage{
feat_stats_to_brik(
  feat_dir,
  out_filename = "feat_stats",
  what = c("cope_files", "z_files", "varcope_files"),
  label_prefix = NULL
)
}
\arguments{
\item{feat_dir}{a .feat directory}

\item{out_filename}{the directory and filename prefix for the output brik/head file}

\item{what}{which elements of the feat structure to add to the brik.}

\item{label_prefix}{an optional character string to add as a prefix to the labels}
}
\value{
a list containing the feat_info for the directory, the full path of the output image,
  a vector of the brik names, and a vector of the contrast names
}
\description{
internal function to convert a feat level 1 directory to an AFNI brik/head file
}
\keyword{internal}
