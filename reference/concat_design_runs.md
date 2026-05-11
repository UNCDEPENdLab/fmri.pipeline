# Concatenate design matrices for each run to form a single design with unique baselines per run (ala AFNI)

Concatenate design matrices for each run to form a single design with
unique baselines per run (ala AFNI)

## Usage

``` r
concat_design_runs(d, convolved = TRUE)
```

## Arguments

- d:

  A design matrix object created by `build_design_matrix`. The
  \$design_convolved element will be used for concatenation.

- convolved:

  If `TRUE` (the default), concatenate the convolved design matrix. If
  `FALSE`, use the unconvolved.
