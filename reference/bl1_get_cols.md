# Onset, duration, value column selection helper function

Onset, duration, value column selection helper function

## Usage

``` r
bl1_get_cols(
  l1_model_set,
  trial_data,
  field_name = NULL,
  field_desc = NULL,
  select_cols = NULL,
  select_regex = NULL,
  limit_cols = NULL,
  force_selection = TRUE,
  alpha_sort = TRUE,
  prompt_input = TRUE
)
```

## Arguments

- l1_model_set:

  an l1_model_set object containing onsets etc.

- trial_data:

  a trial + subjects x events data.frame that contains potential onset
  columns

- field_name:

  the element of `l1_model_set` containing columns of a certain purpose
  (onset, duration, value)

- field_desc:

  the text description of the field being modified (e.g., 'parametric
  value')

- select_cols:

  a character vector of current columns specified by the user to be
  added/included

- select_regex:

  a PCRE-compatible regular expression for identifying columns

- limit_cols:

  a vector of variable names in `trial_data` that constrains what can be
  chosen

- force_selection:

  do not allow an empty return for this field

- alpha_sort:

  whether to display eligible columns in alphabetical order

- prompt_input:

  whether to ask user to confirm selections

## Value

a modified version of l1_model_set that has `field_name` updated
according to user specification
