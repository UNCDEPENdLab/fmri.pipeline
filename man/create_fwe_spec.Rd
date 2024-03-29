% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_fwe_correction.R
\name{create_fwe_spec}
\alias{create_fwe_spec}
\title{helper function to create an FWE object that specifies the type of FWE correction and the 
  model outputs to which it is applied}
\usage{
create_fwe_spec(gpa, level = level, lg = NULL)
}
\arguments{
\item{gpa}{a \code{glm_pipeline_arguments} object that already contains model output information in $l3_model_setup$fsl}

\item{level}{an integer indicating what modeling level this FWE correction applies to (1, 2, or 3)}

\item{lg}{a logger object for logging messages}
}
\description{
helper function to create an FWE object that specifies the type of FWE correction and the 
  model outputs to which it is applied
}
\details{
note that this only works for level = 3 for now
}
