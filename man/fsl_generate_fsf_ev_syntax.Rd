% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsl_generate_fsf_ev_syntax.R
\name{fsl_generate_fsf_ev_syntax}
\alias{fsl_generate_fsf_ev_syntax}
\title{This function generates syntax for FSL Feat .fsf files for the EVs tab of a higher-level fMRI analysis.
It accepts a numeric design matrix whose colum names correspond to individual EVs in the model}
\usage{
fsl_generate_fsf_ev_syntax(inputs, dmat, group_membership = NULL)
}
\arguments{
\item{inputs}{A vector of inputs (usually lower-level .feat directories) corresponding to rows of \code{dmat}.}

\item{dmat}{A numeric matrix whose rows specify individual inputs and columns specify EVs. Columns should be named by EV!}

\item{group_membership}{An optional vector specifying group membership of each row in \code{dmat} in order to
estimate separate variances for multiple groups. If not specified, variances will be modeled as one group.}
}
\value{
A character vector containing .fsf syntax for the contrasts portion of a Feat analysis based on \code{cmat}.
}
\description{
This function generates syntax for FSL Feat .fsf files for the EVs tab of a higher-level fMRI analysis.
It accepts a numeric design matrix whose colum names correspond to individual EVs in the model
}
\details{
This allows you to generate a design matrix dynamically based on a numeric matrix that has been
setup in R using a function such as lm() to handle more complex designs appropriately. The syntax
generated by this function can be combined with a .fsf template file to implement the entire LVL2 or LVL3
analysis setup.
}
\examples{
  dmat <- matrix(rnorm(1000), ncol=10)
  colnames(dmat) <- letters[1:10]

  inputs <- replicate(100, tempfile())

  result <- fsl_generate_fsf_ev_syntax(inputs, dmat)

}
\author{
Michael Hallquist
}
