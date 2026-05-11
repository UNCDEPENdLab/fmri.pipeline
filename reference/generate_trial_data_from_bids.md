# Function to generate a trial_data object from a BIDS-compliant folder

Function to generate a trial_data object from a BIDS-compliant folder

## Usage

``` r
generate_trial_data_from_bids(bids_dir, modality = "func", task_name = "ridl")
```

## Arguments

- bids_dir:

  a directory containing BIDS-compliant processed data for analysis

- modality:

  the subfolder within `bids_dir` that contains data of a certain
  modality. Almost always 'func', which is the default.

- task_name:

  the name of the task, which is appended with `type`

- type:

  at present, always 'task' to denote this part of the BIDS filename...
  Not totally sure what else it could be

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
`modality`, 'task' is the `type`, 'ridl' is the `task_name`, and
'\_postprocessed' is the `suffix`.

## Examples

``` r
if (FALSE) { # \dontrun{
  df <- generate_trial_data_from_bids("/proj/mnhallqlab/no_backup/flanker-fmriprep", task_name = "flanker")
} # }
```
