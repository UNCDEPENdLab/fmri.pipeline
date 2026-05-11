# Core function for computing multivariate time series compression scores by principal components

Core function for computing multivariate time series compression scores
by principal components

## Usage

``` r
compress_mts_pca(mts, pexp_target = 0.9, scale_columns = TRUE)
```

## Arguments

- mts:

  multivariate time series structured time x signals

- pexp_target:

  Proportion of variance explained by principal components. This can be
  a vector, in which case compression is calculated at different
  thresholds.

- scale_columns:

  whether to z-score the time series prior to eigendecomposition
  (recommended)

## Value

a list containing compression estimates of the matrix. For each
pexp_target value, two values are included, one representing the
compression calculated using integer

## Details

This function accepts a time x signals (e.g., voxels) time series
matrix. It computes the singular value decomposition (SVD) and then
examines how many eigenvectors are needed to explain at least
`pexp_target` proportion of variance.

Compression scores are normalized 0 – 1.0 by the equation: 1 -
(n_components / n_timeseries). Thus, if 6 components explain 91

Given that the number of eigenvectors is an integer, a linear
approximation to the exact proportion of variance explained is also
calculated. For example, if 3 components explain 84 92 get us to
precisely 90

## Author

Michael Hallquist
