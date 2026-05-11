# internal helper function to setup a linear model for a given l1, l2, or l3 model

internal helper function to setup a linear model for a given l1, l2, or
l3 model

## Usage

``` r
mobj_fit_lm(mobj = NULL, model_formula = NULL, data, id_cols = NULL, lg = NULL)
```

## Arguments

- mobj:

  a model object to be populated or modified

- model_formula:

  a character string specifying the formula of the model to be fit

- data:

  a data.frame containing all columns used in model fitting

- id_cols:

  a character vector of column names in `data` that identify the
  observations and can be used for merging the model against related
  datasets

## Value

a model object containing the fitted model
