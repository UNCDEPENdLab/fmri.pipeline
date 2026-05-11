# Generate SPM contrasts from an fmri.pipeline L1 model contrast matrix

Generate SPM contrasts from an fmri.pipeline L1 model contrast matrix

## Usage

``` r
generate_spm_contrasts_from_model(
  output_dir,
  mobj,
  spm_path = "/gpfs/group/mnh5174/default/lab_resources/spm12",
  execute = FALSE,
  matlab_cmd = "matlab",
  matlab_args = "-batch",
  average_across_runs = TRUE,
  projection_interaction_terms = character(0),
  projection_interaction_contrast_modes = character(0),
  projection_interaction_run_labels = NULL,
  projection_main_effect_terms = character(0),
  projection_main_effect_weights = NULL
)
```

## Arguments

- output_dir:

  location for SPM outputs and scripts for estimating contrasts. Must
  contain SPM.mat already

- mobj:

  l1_model_spec object containing \$contrasts and regressor names

- spm_path:

  see generate_spm_mat

- execute:

  whether to run contrast setup. This depends on SPM.mat having been
  created already. Default: FALSE

- matlab_cmd:

  see generate_spm_mat

- matlab_args:

  see generate_spm_mat

- average_across_runs:

  whether to average regressor weights across run-specific columns in
  SPM.xX

- projection_interaction_terms:

  optional character vector of projected L2-by-L1 interaction term names

- projection_interaction_contrast_modes:

  optional character vector. Allowed values are 'pooled',
  'session_specific', and 'session_differences'

- projection_interaction_run_labels:

  optional run labels aligned to SPM sessions for naming

- projection_main_effect_terms:

  optional character vector of projected L2 term names to estimate in
  contrast-space against each L1 contrast row

- projection_main_effect_weights:

  optional data.frame or matrix containing projected L2 term values
  aligned to SPM session/run order
