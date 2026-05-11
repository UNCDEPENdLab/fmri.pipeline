# helper class to generate 3dttest++ null datasets using -randomsign

helper class to generate 3dttest++ null datasets using -randomsign

helper class to generate 3dttest++ null datasets using -randomsign

## Methods

### Public methods

- [`simulate_null_3dttest$new()`](#method-simulate_null_3dttest-new)

- [`simulate_null_3dttest$is_complete()`](#method-simulate_null_3dttest-is_complete)

- [`simulate_null_3dttest$submit()`](#method-simulate_null_3dttest-submit)

- [`simulate_null_3dttest$get_3dttest_calls()`](#method-simulate_null_3dttest-get_3dttest_calls)

- [`simulate_null_3dttest$get_permutation_files()`](#method-simulate_null_3dttest-get_permutation_files)

- [`simulate_null_3dttest$get_batch()`](#method-simulate_null_3dttest-get_batch)

- [`simulate_null_3dttest$get_use_sdat()`](#method-simulate_null_3dttest-get_use_sdat)

- [`simulate_null_3dttest$clone()`](#method-simulate_null_3dttest-clone)

------------------------------------------------------------------------

### Method `new()`

initialize a simulate_null_3dttest object to support 3dttest++
-randomsign

#### Usage

    simulate_null_3dttest$new(
      residuals_file = NULL,
      mask_file = NULL,
      njobs = NULL,
      n_permutations = NULL,
      use_sdat = NULL,
      wall_time = NULL,
      memgb_per_3dttest = NULL,
      memgb_combine = NULL
    )

#### Arguments

- `residuals_file`:

  the 4D file containing voxelwise residuals for all subjects (e.g.,
  res4d.nii.gz in FEAT)

- `mask_file`:

  A mask file used to specify which voxels should be analyzed/permuted

- `njobs`:

  The number of independent jobs across which permutations are
  distributed

- `n_permutations`:

  The total number of null datasets to be computed by sign-flipping

- `use_sdat`:

  a logical indicating whether to output null datasets in sdat format
  (single-precision, serialized, I think)

- `wall_time`:

  The amount of compute time needed for this job as a dd-hh:mm:ss string

- `memgb_per_3dttest`:

  The number of gigabytes of memory requested for each 3dttest++
  permutation process. If not specified, 8 will be requested

- `memgb_combine`:

  The number of gigabytes of memory requested for the parent 3dttest++
  job that combines permutation outputs. If not specified, 32 will be
  requested

------------------------------------------------------------------------

### Method `is_complete()`

test whether the expected outputs of the permutation are complete

#### Usage

    simulate_null_3dttest$is_complete()

#### Returns

a boolean indicating whether the expected output file exists

------------------------------------------------------------------------

### Method `submit()`

submit the permutation job to the cluster

#### Usage

    simulate_null_3dttest$submit(force = FALSE)

#### Arguments

- `force`:

  whether to run the permutation job even if the outputs are already
  complete

------------------------------------------------------------------------

### Method `get_3dttest_calls()`

build the 3dttest++ command strings for each permutation job

#### Usage

    simulate_null_3dttest$get_3dttest_calls(include_complete = FALSE)

#### Arguments

- `include_complete`:

  If TRUE, include calls whose output files already exist. If FALSE,
  return only calls with missing outputs.

#### Returns

A character vector of 3dttest++ calls.

------------------------------------------------------------------------

### Method `get_permutation_files()`

return the permutation output files once available

#### Usage

    simulate_null_3dttest$get_permutation_files()

#### Returns

A named character vector with entries for the mask file and permutation
file.

------------------------------------------------------------------------

### Method `get_batch()`

return the R_batch_job used to run 3dttest permutations

#### Usage

    simulate_null_3dttest$get_batch()

#### Returns

An R_batch_job object.

------------------------------------------------------------------------

### Method `get_use_sdat()`

indicate whether sdat output is enabled

#### Usage

    simulate_null_3dttest$get_use_sdat()

#### Returns

A logical scalar.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    simulate_null_3dttest$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
