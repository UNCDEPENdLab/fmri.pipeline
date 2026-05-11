# compute spike regressors from a data.frame or matrix of motion parameters and a set of expressions that are evaluated against this parameter matrix.

compute spike regressors from a data.frame or matrix of motion
parameters and a set of expressions that are evaluated against this
parameter matrix.

## Usage

``` r
compute_spike_regressors(mot = NULL, spike_volume = NULL, lg = NULL)
```

## Arguments

- mot:

  a data.frame or matrix with volumes on rows and named motion
  parameters on columns

- spike_volume:

  a character vector of expressions to evaluate against `mot`. Resulting
  columns will be prefixed with the names of each expression. If the
  expressions are unnamed, prefixes will be expr1\_, expr2\_, etc.

## Value

a matrix of spike regressors (volumes on rows, spike regressors on
columns)
