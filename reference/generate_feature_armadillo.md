# This function creates K shifts of a neural events vector according to the kernel length, K.

This function creates K shifts of a neural events vector according to
the kernel length, K.

## Arguments

- encoding:

  The neural events vector (same length as BOLD time series)

- K:

  The length of the kernel

## Value

A matrix of length(encoding) rows and K columns, where each column
contains a successively lagged copy of the encoding vector

## Details

This is an internal function that is used inside a while loop by
deconvolve_nlreg. Profiling of the algorithm revealed that this is the
primary bottleneck, so I ported it to an Rcpp function

## Author

Michael Hallquist
