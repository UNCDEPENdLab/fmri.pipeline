# Interactive function to build an l1 model specification for setup_glm_pipeline

Interactive function to build an l1 model specification for
setup_glm_pipeline

## Usage

``` r
build_l1_models(
  gpa = NULL,
  trial_data = NULL,
  ppi_data = NULL,
  l1_model_set = NULL,
  from_spec_file = NULL,
  onset_cols = NULL,
  onset_regex = ".*(onset|time).*",
  duration_cols = NULL,
  duration_regex = ".*duration.*",
  value_cols = NULL,
  value_regex = NULL,
  isi_cols = NULL,
  isi_regex = "^(iti|isi).*"
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object containing an analysis pipeline to
  which \$l1_models should be added. If \$l1_models is already present,
  these will be amended.

- trial_data:

  a data.frame containing trial-level data for one or more subjects

- ppi_data:

  an optional data.frame containing physiological time-series inputs for
  PPI signal construction

- l1_model_set:

  optional existing l1_model_set to be modified

- from_spec_file:

  optional YAML or JSON file containing settings to populated into l1
  models

- onset_cols:

  an optional character vector of columns in `trial_data` that should be
  in the set of event onsets

- onset_regex:

  an optional PCRE-compatible regular expression for identifying
  potential event onset columns in `trial_data`

- duration_cols:

  an optional character vector of columns in `trial_data` that should be
  in the set of event durations

- duration_regex:

  an optional PCRE-compatible regular expression for identifying
  potential event duration columns in `trial_data`

- value_cols:

  an optional character vector of columns in `trial_data` that should be
  in the set of signal values

- value_regex:

  an optional PCRE-compatible regular expression for identifying
  potential event value columns in `trial_data`

- isi_cols:

  an optional character vector of columns in `trial_data` that should be
  in the set of signal isi/iti

- isi_regex:

  an optional PCRE-compatible regular expression for identifying
  potential isi/iti columns in `trial_data`

## Value

a `l1_model_set` object containing events, signals, and models,
compatible with build_design_matrix

## Details

if `gpa` is not passed in, then we will work from trial_data and
l1_model_set.

## Author

Michael Hallquist
