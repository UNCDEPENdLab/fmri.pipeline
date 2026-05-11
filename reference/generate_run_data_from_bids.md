# Function to generate a run_data object from a BIDS-compliant folder

Function to generate a run_data object from a BIDS-compliant folder

## Usage

``` r
generate_run_data_from_bids(
  bids_dir,
  modality = "func",
  task_name = "ridl",
  desc = "postproc",
  suffix = "bold",
  space = NULL,
  anat_root = NULL,
  fmap_root = NULL
)
```

## Arguments

- bids_dir:

  a directory containing BIDS-compliant processed data for analysis

- modality:

  the subfolder within `bids_dir` that contains data of a certain
  modality. Almost always 'func', which is the default.

- task_name:

  the name of the task, which is appended with "task-"

- suffix:

  an optional suffix in the expected filename (just before the file
  extension)

## Value

a data.frame containing all run_nifti and confound_input_file results
for subjects in the folder

## Details

The files should generally have a name like
sub-220256_task-ridl3_space-MNI152NLin2009cAsym_desc-preproc_bold_postprocessed.nii.gz
and be located in a folder like:
/proj/mnhallqlab/proc_data/sub-220256/func/ where 'func' is the
`modality`, 'ridl' is the `task_name`, and '\_postprocessed' is the
`suffix`. If the expected `modality` folder is missing, the function
will also search session-level folders (e.g., `sub-<id>/ses-<id>/`) or
the subject root for NIfTI files, and fall back to any `.nii/.nii.gz`
files found directly in those directories.

## Examples

``` r
if (FALSE) { # \dontrun{
  run_df <- generate_run_data_from_bids("/proj/mnhallqlab/no_backup/flanker-fmriprep", task_name = "flanker", desc = "postproc")
} # }
```
