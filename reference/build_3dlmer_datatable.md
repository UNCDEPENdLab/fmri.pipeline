# Internal function to build the data table for AFNI 3dLMEr

Internal function to build the data table for AFNI 3dLMEr

## Usage

``` r
build_3dlmer_datatable(subject_data, input_files, model_variables)
```

## Arguments

- subject_data:

  a data.frame containing subject/session-level covariates with unique
  id/session rows

- input_files:

  a data.frame with columns 'id', 'session', and 'InputFile'

- model_variables:

  a character vector of variables that must be included in the table

## Value

a data.frame formatted for 3dLMEr -dataTable
