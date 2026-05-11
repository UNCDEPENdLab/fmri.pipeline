# R6 class for running 3dFWHMx on a single input file based on user specification

R6 class for running 3dFWHMx on a single input file based on user
specification

R6 class for running 3dFWHMx on a single input file based on user
specification

## Details

N.B. This class doesn't even expose the Gaussian ACF options given false
positive problems

## Methods

### Public methods

- [`afni_3dfwhmx$new()`](#method-afni_3dfwhmx-new)

- [`afni_3dfwhmx$run()`](#method-afni_3dfwhmx-run)

- [`afni_3dfwhmx$get_call()`](#method-afni_3dfwhmx-get_call)

- [`afni_3dfwhmx$get_acf_params()`](#method-afni_3dfwhmx-get_acf_params)

- [`afni_3dfwhmx$get_acf_by_radius()`](#method-afni_3dfwhmx-get_acf_by_radius)

- [`afni_3dfwhmx$get_fwhm_by_volume()`](#method-afni_3dfwhmx-get_fwhm_by_volume)

- [`afni_3dfwhmx$get_input_file()`](#method-afni_3dfwhmx-get_input_file)

- [`afni_3dfwhmx$get_mask_file()`](#method-afni_3dfwhmx-get_mask_file)

- [`afni_3dfwhmx$get_outputs()`](#method-afni_3dfwhmx-get_outputs)

- [`afni_3dfwhmx$is_complete()`](#method-afni_3dfwhmx-is_complete)

- [`afni_3dfwhmx$delete_outputs()`](#method-afni_3dfwhmx-delete_outputs)

- [`afni_3dfwhmx$clone()`](#method-afni_3dfwhmx-clone)

------------------------------------------------------------------------

### Method `new()`

#### Usage

    afni_3dfwhmx$new(
      input_file = NULL,
      mask_file = NULL,
      out_dir = NULL,
      demed = NULL,
      unif = NULL,
      average = "geometric",
      ncpus = 1L
    )

#### Arguments

- `input_file`:

  The input dataset whose smoothness should be calculated (often
  first-level GLM residuals)

- `mask_file`:

  Only compute smoothness within this mask (if not provided, -automask
  will be used)

- `out_dir`:

  The output directory for fwhmx files

- `demed`:

  If `TRUE`, subtract the median of each voxels time series before
  calculating FWHM (wraps -demed). Default: FALSE.

- `unif`:

  If `TRUE`, normalize each voxel's time series to have the same MAD
  before calculating FWHM (wraps -unif). Default: FALSE

- `average`:

  For multi-volume data, compute the averaged ACF estimates by either
  the geometric or arithmetic mean. Default: "geometric".

- `ncpus`:

  The number of threads/cores to use for running 3dFWHMx. Default is 1.
  Controls OMP_NUM_THREADS.

------------------------------------------------------------------------

### Method `run()`

runs 3dFWHMx on this input dataset

#### Usage

    afni_3dfwhmx$run(force = FALSE)

#### Arguments

- `force`:

  If `TRUE`, 3dFWHMx will be run even if the expected output files
  already exist.

------------------------------------------------------------------------

### Method `get_call()`

return the 3dFWHMx call used for the specified input

#### Usage

    afni_3dfwhmx$get_call()

#### Details

this is useful if you want to call 3dFWHMx yourself directly or if you
want to debug the 3dFWHMx call specification.

#### Returns

a character string with the 3dFWHMx call

------------------------------------------------------------------------

### Method `get_acf_params()`

return the estimated ACF parameters for this run of data, averaged over
volumes

#### Usage

    afni_3dfwhmx$get_acf_params()

#### Details

Will issue a warning if 3dFWHMx has not run successfully on this dataset
already

#### Returns

a three-element vector containing the ACF estimates for the dataset,
averaging over volumes.

------------------------------------------------------------------------

### Method `get_acf_by_radius()`

return the estimated ACF parameters for this run of data

#### Usage

    afni_3dfwhmx$get_acf_by_radius()

#### Details

Will issue a warning if 3dFWHMx has not run successfully on this dataset
already

#### Returns

a three-element vector containing the ACF estimates for the dataset,
averaging over volumes.

------------------------------------------------------------------------

### Method `get_fwhm_by_volume()`

return the estimated ACF parameters for this run of data

#### Usage

    afni_3dfwhmx$get_fwhm_by_volume()

#### Details

Will issue a warning if 3dFWHMx has not run successfully on this dataset
already

#### Returns

a three-element vector containing the ACF estimates for the dataset,
averaging over volumes.

------------------------------------------------------------------------

### Method `get_input_file()`

return the input file (NIfTI) used for 3dFWHMx

#### Usage

    afni_3dfwhmx$get_input_file()

#### Returns

a character string of the input file location

------------------------------------------------------------------------

### Method `get_mask_file()`

return the mask file used for 3dFWHMx estimation

#### Usage

    afni_3dfwhmx$get_mask_file()

#### Returns

a character string of the mask file location

------------------------------------------------------------------------

### Method `get_outputs()`

return the expected output files related to this 3dFWHMx object

#### Usage

    afni_3dfwhmx$get_outputs()

#### Returns

a character vector containing expected output files

------------------------------------------------------------------------

### Method `is_complete()`

method to indicate whether 3dFWHMx has already run and completed for
this input

#### Usage

    afni_3dfwhmx$is_complete()

#### Returns

TRUE if expected 3dFWHMx output file exists, FALSE if it does not

------------------------------------------------------------------------

### Method `delete_outputs()`

method to delete any/all files generated by this object

#### Usage

    afni_3dfwhmx$delete_outputs(prompt = FALSE)

#### Arguments

- `prompt`:

  if TRUE, user will have to confirm deletion of each file. If FALSE,
  files are deleted without prompting.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    afni_3dfwhmx$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
