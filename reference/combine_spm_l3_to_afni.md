# Combine SPM L3 outputs into AFNI BRIK+HEAD files for visualization

Combine SPM L3 outputs into AFNI BRIK+HEAD files for visualization

## Usage

``` r
combine_spm_l3_to_afni(
  gpa,
  spm_l3_combined_filename = NULL,
  spm_l3_combined_briknames = NULL,
  template_brain = NULL
)
```

## Arguments

- gpa:

  a glm_pipeline_arguments object having a populated l3_model_setup
  field.

- spm_l3_combined_filename:

  a glue expression for the path and filename prefix.

- spm_l3_combined_briknames:

  a glue expression for naming the subbriks in the AFNI output.

- template_brain:

  an optional filename for the MNI template that should be used as an
  underlay in AFNI.

## Value

A data.frame describing the images that were combined.

## Details

This function mirrors the logic of combine_feat_l3_to_afni, but for SPM
L3 outputs. It detects contrast maps (con\_\*.nii) and statistic maps
(spmT\_\*.nii / spmF\_\*.nii) and concatenates them into AFNI BRIK/HEAD
files with informative labels.
