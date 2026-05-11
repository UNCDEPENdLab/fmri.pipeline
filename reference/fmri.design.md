# fmri.design function from `fmri` R package. Brought into this package to avoid dependency on `fmri` package, especially the `gsl` package dependency in `aws`, which is finicky on clusters if the GSL module isn't available.

fmri.design function from `fmri` R package. Brought into this package to
avoid dependency on `fmri` package, especially the `gsl` package
dependency in `aws`, which is finicky on clusters if the GSL module
isn't available.

## Usage

``` r
fmri.design(stimulus, order = 2, cef = NULL, verbose = FALSE)
```

## Arguments

- stimulus:

  The regressor of interest

- order:

  The polynomial order of baseline/drift terms to incorporate in the
  design matrix

- cef:

  confound regressors to add

- verbose:

  if TRUE, output detailed progress about each step in the design setup
