# This function reads an SPM.mat file and generate contrasts based on the design matrix specification

This function reads an SPM.mat file and generate contrasts based on the
design matrix specification

## Usage

``` r
generate_spm_contrasts(
  output_dir,
  condition_contrasts = TRUE,
  unit_contrasts = TRUE,
  effects_of_interest_F = TRUE,
  spm_path = "/gpfs/group/mnh5174/default/lab_resources/spm12",
  execute = FALSE,
  matlab_cmd = "matlab",
  matlab_args = "-batch"
)
```

## Arguments

- output_dir:

  location for SPM outputs and scripts for estimating contrasts. Must
  contain SPM.mat already

- condition_contrasts:

  see generate_spm_mat

- unit_contrasts:

  see generate_spm_mat

- effects_of_interest_F:

  see generate_spm_mat

- spm_path:

  see generate_spm_mat

- execute:

  whether to run contrast setup. This depends on SPM.mat having been
  created already. Default: FALSE

- matlab_cmd:

  see generate_spm_mat

- matlab_args:

  see generate_spm_mat

## Author

Michael Hallquist
