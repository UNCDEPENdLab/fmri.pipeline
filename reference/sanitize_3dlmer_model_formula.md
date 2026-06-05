# Internal function to construct the 3dLMEr command string

Internal function to construct the 3dLMEr command string

## Usage

``` r
sanitize_3dlmer_model_formula(model_formula)
```

## Arguments

- model_formula:

  fixed effects formula string

- prefix:

  output prefix for 3dLMEr

- qVars:

  character vector of quantitative variables

- glt_codes:

  list of named GLT code strings

- data_table_file:

  path to the data table text file

- mask:

  path to the brain mask file

- njobs:

  number of parallel jobs for 3dLMEr

- ss_type:

  sum of squares type (default 3)

## Value

a character string containing the 3dLMEr command
