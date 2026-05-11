# R6 class for 3dClustSim automation

R6 class for 3dClustSim automation

R6 class for 3dClustSim automation

## Public fields

- `fwhmx_set`:

  A afni_3dfwhmx_list object containing 3dFWHMx information for all
  fwhmx_input_files

- `inset_files`:

  A character vector of files to use directly as volumes to threshold
  and clusterize

- `insdat_file`:

  An sdat file containing permutations to be passed to 3dClustSim
  through -insdat

- `insdat_mask_file`:

  A mask file corresponding to the insdat_file data that indicates where
  each value is in space

- `null_3dttest_obj`:

  only used if residuals_file is passed in, this contains the object for
  running the permutations

## Methods

### Public methods

- [`afni_3dclustsim$new()`](#method-afni_3dclustsim-new)

- [`afni_3dclustsim$submit()`](#method-afni_3dclustsim-submit)

- [`afni_3dclustsim$get_clustsim_df()`](#method-afni_3dclustsim-get_clustsim_df)

- [`afni_3dclustsim$get_clustsim_output_files()`](#method-afni_3dclustsim-get_clustsim_output_files)

- [`afni_3dclustsim$get_call()`](#method-afni_3dclustsim-get_call)

- [`afni_3dclustsim$get_out_dir()`](#method-afni_3dclustsim-get_out_dir)

- [`afni_3dclustsim$get_ncpus()`](#method-afni_3dclustsim-get_ncpus)

- [`afni_3dclustsim$use_fwhmx_acf()`](#method-afni_3dclustsim-use_fwhmx_acf)

- [`afni_3dclustsim$is_complete()`](#method-afni_3dclustsim-is_complete)

- [`afni_3dclustsim$refresh()`](#method-afni_3dclustsim-refresh)

- [`afni_3dclustsim$apply_clustsim()`](#method-afni_3dclustsim-apply_clustsim)

- [`afni_3dclustsim$clone()`](#method-afni_3dclustsim-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new afni_3dclustsim object

#### Usage

    afni_3dclustsim$new(
      out_dir = NULL,
      prefix = NULL,
      fwhmx_input_files = NULL,
      fwhmx_mask_files = NULL,
      residuals_file = NULL,
      residuals_mask_file = NULL,
      residuals_njobs = NULL,
      inset_files = NULL,
      insdat_file = NULL,
      insdat_mask_file = NULL,
      dxyz = NULL,
      nxyz = NULL,
      clustsim_mask = NULL,
      acf_params = NULL,
      nopad = NULL,
      pthr = NULL,
      athr = NULL,
      iter = NULL,
      nodec = NULL,
      seed = NULL,
      scheduler = NULL,
      ncpus = NULL
    )

#### Arguments

- `out_dir`:

  the intended output directory for 3dClustSim files

- `prefix`:

  the prefix to be included in the names of 3dClustSim output files

- `fwhmx_input_files`:

  A character vector of input files to be passed through 3dFWHMx (-ACF
  method)

- `fwhmx_mask_files`:

  A character vector of masks containing the volume over which to
  estimate the ACF in 3dFWHMx. Must match 1:1 with `fwhmx_input_files`.

- `residuals_file`:

  The filename of the group residuals file to be used in null dataset
  generation (permutation approach).

- `residuals_mask_file`:

  The volume over which null datasets should be generated from the
  residuals

- `residuals_njobs`:

  The number of independent jobs for splitting up the residuals
  permutations. If NULL, 32 jobs will be used.

- `inset_files`:

  A character vector of dataset files to pass directly to 3dClustSim via
  -inset.

- `insdat_file`:

  An sdat file containing null datasets to be passed to 3dClustSim via
  -insdat.

- `insdat_mask_file`:

  A mask file corresponding to the insdat_file that maps values to
  space.

- `dxyz`:

  the size of voxels in x y z (vector of 3 numbers)

- `nxyz`:

  the number of voxels along x y z (vector of 3 positive integers)

- `clustsim_mask`:

  This controls the volume over which to correct for FWE using
  3dClustSim. If you give a whole-brain mask, then your cluster
  thresholds reflect whole-brain FWE correction. If you give a smaller
  mask (e.g., a single region or network), you are correcting only over
  that volume (i.e., a small-volume correction).

- `acf_params`:

  a vector of 3 autocorrelation parameters (a, b, c) to be used to
  simulate smoothness. Usually produced by 3dFWHMx. The 'a' parameter
  must be between 0 and 1. The 'b' and 'c' parameters (scale radii) must
  be positive. The spatial autocorrelation function is given by:
  \`ACF(r) = a \* exp(-r\*r/(2\*b\*b)) + (1-a)\*exp(-r/c)\`

- `nopad`:

  If TRUE, disable 3dClustSim's default 'padding' slices along each face
  to allow for edge effects of the smoothing process. Default: FALSE

- `pthr`:

  A vector of voxelwise p-values to be tested in 3dClustSim (-pthr). Can
  be a space-separated string or a numeric vector. Default: ".01 .005
  .002 .001 .0005 .0002 .0001".

- `athr`:

  A vector of cluster p-values to be tested in 3dClustSim (-athr). Can
  be a space-separated string or a numeric vector. Default: ".05 .02 .01
  .005 .002 .001 .0005 .0002 .0001".

- `iter`:

  The number of iterations to use in simulating null datasets in
  3dClustSim. Default: 30000.

- `nodec`:

  If TRUE, clusters will be printed without decimal places, rounding up
  (e.g., 27.2 becomes 28). Default: FALSE.

- `seed`:

  The seed to use when starting 3dClustSim random number generation.
  Default: 0 (sets a random seed)

- `scheduler`:

  The HPC scheduler to use. Can be 'local', 'slurm', or 'torque'.

- `ncpus`:

  The number of cores to use in the 3dClustSim job. This sets
  OMP_NUM_THREADS in the 3dClustSim job to speed up computation.

------------------------------------------------------------------------

### Method `submit()`

submit the 3dClustSim compute job to the cluster

#### Usage

    afni_3dclustsim$submit(force = FALSE)

#### Arguments

- `force`:

  if TRUE, re-estimate 3dClustSim even though its output files already
  exist

#### Details

Note that if fwhmx_input_files are provided at the corresponding 3dFWHMx
has not been run yet, this job will be submitted and serve as a
dependency for the 3dClustsim job. Likewise, if a residuals_file is
provided, then the 3dttest++ permutation step will be submitted as a job
and 3dClustSim will be dependent on this finishing.

------------------------------------------------------------------------

### Method `get_clustsim_df()`

return a data.frame containing the results of the 3dClustSim for all NN
levels and sidedness.

#### Usage

    afni_3dclustsim$get_clustsim_df()

------------------------------------------------------------------------

### Method `get_clustsim_output_files()`

return a character vector of output text files generated by 3dClustSim

#### Usage

    afni_3dclustsim$get_clustsim_output_files()

------------------------------------------------------------------------

### Method `get_call()`

return the 3dClustSim call related to these inputs. This is what will be
passed to the scheduler.

#### Usage

    afni_3dclustsim$get_call()

------------------------------------------------------------------------

### Method `get_out_dir()`

return the output directory for 3dClustSim files

#### Usage

    afni_3dclustsim$get_out_dir()

------------------------------------------------------------------------

### Method `get_ncpus()`

return the number of cpus (cores) to be used by 3dClustsim (via
OMP_NUM_THREADS)

#### Usage

    afni_3dclustsim$get_ncpus()

------------------------------------------------------------------------

### Method `use_fwhmx_acf()`

return TRUE/FALSE for whether this clustsim relies on the ACF estimates
from 3dFWHMx

#### Usage

    afni_3dclustsim$use_fwhmx_acf()

------------------------------------------------------------------------

### Method `is_complete()`

return whether 3dClustSim has completed for this input

#### Usage

    afni_3dclustsim$is_complete()

------------------------------------------------------------------------

### Method `refresh()`

Simple method to refresh the clustsim_df, fwhmx files, and 3dttest++
permutation files

#### Usage

    afni_3dclustsim$refresh()

#### Details

The refresh is useful if the object needs to be updated just in time to
determine whether permutations or ACF params are available.

------------------------------------------------------------------------

### Method `apply_clustsim()`

method to apply 3dClustSim results to a statistic image to create an
integer-valued cluster mask and/or a thresholded statistic image

#### Usage

    afni_3dclustsim$apply_clustsim(
      statistic_nifti = NULL,
      NN = 1,
      sided = "bi",
      athr = 0.05,
      pthr = 0.001,
      voxelwise_stat = list(stat_type = "z"),
      output_cluster_mask = TRUE,
      output_thresholded_image = FALSE,
      add_whereami = TRUE,
      whereami_atlases = NULL
    )

#### Arguments

- `statistic_nifti`:

  A filename to a NIfTI image containing voxelwise statistics that
  should be used to calculate the threshold for pthr

- `NN`:

  The cluster definition from 3dClustSim to be applied when thresholding
  the image (1, 2, or 3)

- `sided`:

  Whether to apply the cluster threshold for one-sided, two-sided, or
  bi-sided tests ('1', '2', or 'bi')

- `athr`:

  The clusterwise p-value to be applied. Default: .05

- `pthr`:

  The voxelwise threshold to be applied to the statistic_file. Default:
  .001.

- `voxelwise_stat`:

  A list object specifying the statistic contained in the
  `statistic_nifti`. At present, this consists of the `stat_type`: 'z',
  't', or 'p'. If a 't' statistic is passed, also include \$df in the
  list specifying the degrees of freedom.

- `output_cluster_mask`:

  If `TRUE`, output an integer-valued mask containing any whole
  brain-significant clusters (calculated by 3dClusterize).

- `output_thresholded_image`:

  If `TRUE`, output a copy of `statistic_nifti` that has been
  thresholded according to the settings specified here.

- `add_whereami`:

  Whether to also call AFNI whereami to get anatomical landmarks of
  interest. Default: TRUE

- `whereami_atlases`:

  An optional character vector of atlases to be requested in whereami.

#### Returns

an afni_3dclusterize object containing clusters in `statistic_nifti`
that survive the 3dClustSim correction requested.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    afni_3dclustsim$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
