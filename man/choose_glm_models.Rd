% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm_helper_functions.R
\name{choose_glm_models}
\alias{choose_glm_models}
\title{helper function to ask user to choose models at a given level for further processing}
\usage{
choose_glm_models(gpa, model_names, level, lg = NULL)
}
\arguments{
\item{gpa}{a \code{glm_pipeline_arguments} with models specified at a given level}

\item{model_names}{a user-specified string of models to process/use}

\item{level}{the level of GLM analysis to be specified (1, 2, or 3)}
}
\value{
a character vector of user-specified model names at this level
}
\description{
helper function to ask user to choose models at a given level for further processing
}
\keyword{internal}
