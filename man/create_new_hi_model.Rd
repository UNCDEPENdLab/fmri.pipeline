% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_l2_models.R
\name{create_new_hi_model}
\alias{create_new_hi_model}
\title{Internal function to walk through l2/l3 model setup, including populated from a specification (YAML) file}
\usage{
create_new_hi_model(
  data,
  to_modify = NULL,
  level = NULL,
  cur_model_names = NULL,
  spec_list = NULL,
  lg = NULL
)
}
\arguments{
\item{data}{a data.frame containing data that can be included in the model}

\item{to_modify}{an existing `hi_model_spec` object to be modified}

\item{level}{the model level: 2 or 3}

\item{cur_model_names}{a character vector of the names of all current models in the model list}

\item{spec_list}{a list containing model specifications read in from a YAML/JSON file}

\item{lg}{the current logger}
}
\value{
a `hi_model_spec` object containing the L2/L3 model to be added to gpa$l2_models or gpa$l3_models
}
\description{
Internal function to walk through l2/l3 model setup, including populated from a specification (YAML) file
}
\details{
if `spec_list` is passed in, the resulting model object will be populated from the fields
  of the specification list, rather than prompting the user for input
}
\keyword{internal}
