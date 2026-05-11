# R6 class for running 3dFWHMx on a group of input files using a scheduler/cluster

R6 class for running 3dFWHMx on a group of input files using a
scheduler/cluster

R6 class for running 3dFWHMx on a group of input files using a
scheduler/cluster

## Methods

### Public methods

- [`afni_3dfwhmx_list$new()`](#method-afni_3dfwhmx_list-new)

- [`afni_3dfwhmx_list$refresh()`](#method-afni_3dfwhmx_list-refresh)

- [`afni_3dfwhmx_list$submit()`](#method-afni_3dfwhmx_list-submit)

- [`afni_3dfwhmx_list$get_batch()`](#method-afni_3dfwhmx_list-get_batch)

- [`afni_3dfwhmx_list$get_acf_average()`](#method-afni_3dfwhmx_list-get_acf_average)

- [`afni_3dfwhmx_list$get_acf_df()`](#method-afni_3dfwhmx_list-get_acf_df)

- [`afni_3dfwhmx_list$get_effective_fwhm()`](#method-afni_3dfwhmx_list-get_effective_fwhm)

- [`afni_3dfwhmx_list$is_complete()`](#method-afni_3dfwhmx_list-is_complete)

- [`afni_3dfwhmx_list$clone()`](#method-afni_3dfwhmx_list-clone)

------------------------------------------------------------------------

### Method `new()`

create a new afni_3dfwhmx_list object

#### Usage

    afni_3dfwhmx_list$new(
      input_files = NULL,
      mask_files = NULL,
      scheduler = NULL,
      wall_time = "10:00:00"
    )

#### Arguments

- `input_files`:

  A vector of input files that should each be passed to 3dFWHMx (usually
  first-level residuals)

- `mask_files`:

  A vector of mask files corresponding to `input_files` (i.e., these
  should align), use to set the volume over which smoothness is
  estimated by 3dFWHMx

- `scheduler`:

  The HPC scheduler to be used: 'local', 'slurm', or 'torque'. Default:
  'slurm'

- `wall_time`:

  The total time required to run all 3dFWHMx objects on the job
  scheduler/HPC

------------------------------------------------------------------------

### Method `refresh()`

recreate afni_3dfwhmx objects and completion status

#### Usage

    afni_3dfwhmx_list$refresh(...)

#### Arguments

- `...`:

  Passes through any arguments to afni_3dfwhmx\$new()

#### Details

useful for updating job status and ACF params after a run completes

------------------------------------------------------------------------

### Method `submit()`

Submit the 3dFWHMx batch for all inputs to the scheduler

#### Usage

    afni_3dfwhmx_list$submit(force = FALSE)

#### Arguments

- `force`:

  If `TRUE`, completed 3dFWHMx runs will be included in the batch.
  Default: FALSE

------------------------------------------------------------------------

### Method `get_batch()`

Return R_batch_job objects used to submit 3dFWHMx for all inputs

#### Usage

    afni_3dfwhmx_list$get_batch(include_complete = FALSE)

#### Arguments

- `include_complete`:

  If TRUE, complete 3dFWHMx jobs will be included in the batch, leading
  these to be re-run

------------------------------------------------------------------------

### Method `get_acf_average()`

method to calculate the overall ACF average across all input datasets

#### Usage

    afni_3dfwhmx_list$get_acf_average(allow_incomplete_fwhmx = FALSE)

#### Arguments

- `allow_incomplete_fwhmx`:

  If `TRUE`, the average will be calculated even if some datasets do not
  have 3dFWHMx parameter estimates (e.g., if they crashed or were never
  run). Default. `FALSE`.

#### Returns

a three-element vector containing the ACF estimates averaged across all
datasets

------------------------------------------------------------------------

### Method `get_acf_df()`

return the 3dFWHMx results as a data.frame

#### Usage

    afni_3dfwhmx_list$get_acf_df()

------------------------------------------------------------------------

### Method `get_effective_fwhm()`

return the effective FWHM estimated by 3dFWHMx -ACF

#### Usage

    afni_3dfwhmx_list$get_effective_fwhm()

------------------------------------------------------------------------

### Method `is_complete()`

Simple method to return whether 3dFWHMx is complete for all input files

#### Usage

    afni_3dfwhmx_list$is_complete()

#### Returns

returns `TRUE` if 3dFWHMx has completed for all runs

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    afni_3dfwhmx_list$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
