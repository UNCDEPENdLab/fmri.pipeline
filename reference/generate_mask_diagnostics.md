# Generate a set of diagnostic images and plots for brain masks relative to one another and a template

Generate a set of diagnostic images and plots for brain masks relative
to one another and a template

## Usage

``` r
generate_mask_diagnostics(
  input = NULL,
  reference_mask = NULL,
  reference_t1 = NULL,
  output_directory = NULL,
  generate_automask = NULL,
  generate_run_plots = TRUE,
  ncores = 8L
)
```

## Arguments

- input:

  Can be a character vector of mask files, a `gpa` object, or a
  `data.frame` containing a `brain_mask` field, or at least a
  `$run_nifti` field. If a `gpa` object is provided, we will use the
  `$run_data` `data.frame`.

- reference_mask:

  a string pointing to the reference mask used for run brain mask
  comparisons. Must be in the same stereotaxic space and have the same
  dimensions as the run brain masks.

- reference_t1:

  a string pointing to the T1w image of the reference brain. Used to aid
  interpretability of brain mask visualizations.

- output_directory:

  a string pointing to a directory for mask diagnostics. Defaults to
  "mask_diagnostics" in the current directory. This directory will be
  created if it does not exist.

- generate_automask:

  a logical indicating whether to use AFNI 3dAutomask to generate brain
  masks if they are not provided in the \$brain_mask field. If NULL
  (default), the user will be prompted for whether they wish to generate
  missing brain masks. If `FALSE`, masks will not be generate, and if
  `TRUE`, they will be generated.
