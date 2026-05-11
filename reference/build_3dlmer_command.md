# Internal function to construct the 3dLMEr command string

Internal function to construct the 3dLMEr command string

## Usage

``` r
build_3dlmer_command(
  prefix,
  model_formula,
  qVars = NULL,
  glt_codes = NULL,
  data_table_file,
  mask = NULL,
  njobs = 1,
  ss_type = 3
)
```

## Arguments

- prefix:

  output prefix for 3dLMEr

- model_formula:

  fixed effects formula string

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
