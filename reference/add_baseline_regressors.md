# Helper function to append baseline regressors to the convolved design matrix. These are especially useful for GLMs that use concatenated data for multiple runs, where run-specific intercepts and drift terms are needed.

Helper function to append baseline regressors to the convolved design
matrix. These are especially useful for GLMs that use concatenated data
for multiple runs, where run-specific intercepts and drift terms are
needed.

## Usage

``` r
add_baseline_regressors(
  dmat_convolved,
  baseline_coef_order = -1L,
  baseline_parameterization = "Legendre"
)
```

## Arguments

- dmat_convolved:

  The design matrix containing convolved regressors for all events,
  created by `build_design_matrix`

- baseline_coef_order:

  The polynomial order of baseline regressors to be added

- baseline_parameterization:

  The parameterization of baseline regressors. Default is "Legendre".
  Alternative is "orthogonal.polynomials".

## Value

A modified version of `dmat_convolved` that contains the baseline
regressors requested by the user.
