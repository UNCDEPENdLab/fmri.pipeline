# Function to combine Feat L3 analyses into AFNI BRIK+HEAD files according to user-specified combinations

Function to combine Feat L3 analyses into AFNI BRIK+HEAD files according
to user-specified combinations

## Usage

``` r
combine_feat_l3_to_afni(
  gpa,
  feat_l3_combined_filename = NULL,
  feat_l3_combined_briknames = NULL,
  template_brain = NULL
)
```

## Arguments

- gpa:

  a glm_pipeline_arguments object having a populated l3_model_setup
  field. The L3 analyses should also be complete (i.e., after
  run_glm_pipeline).

- feat_l3_combined_filename:

  a glue expression for the path and filename prefix

- feat_l3_combined_briknames:

  a glue expression for naming the subbriks in the AFNI output

- template_brain:

  an optional filename for the MNI template that should be used as an
  underlay in AFNI. This image will be symbolically linked into each
  directory created by `combine_feat_l3_to_afni`.

## Details

To specify the folder and filename structure for the combined feat
analyses, use a `glue` expression that indicates how outputs should be
structured. In particular, variables in gpa\$l3_model_setup\$fsl can be
used for dynamically naming of afni outputs. Rows that evaluate to the
same `feat_l3_combined_filename` are stitched into the same AFNI file,
so filename variables control grouping. The multi-run default
intentionally excludes `l2_cope_name` from the filename and includes it
in `feat_l3_combined_briknames`; this keeps per-cope L2 inputs grouped
into a compact AFNI file per L1 contrast, with L2/L3 contrast names
preserved as sub-brik labels.

- l1_model:

  the name of the level 1 model

- l2_model:

  the name of the level 2 model

- l3_model:

  the name of the level 3 model

- l1_cope_name:

  the level 1 contrast name

- l2_cope_name:

  the level 2 contrast name

- l3_cope_name:

  the level 3 contrast name

## Examples

``` r
if (FALSE) { # \dontrun{
  feat_l3_combined_filename <- "{gpa$output_directory}/afni_combined/L1m-{l1_model}/l1c-{l1_cope_name}/L3m-{l3_model}_stats"
  feat_l3_combined_briknames = "l2c-{l2_cope_name}_l3c-{l3_cope_name}"
  template_brain <- "/proj/mnhallqlab/lab_resources/standard/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_brain.nii"

  combine_feat_l3_to_afni(gpa, feat_l3_combined_filename, feat_l3_combined_briknames, template_brain)
} # }
```
