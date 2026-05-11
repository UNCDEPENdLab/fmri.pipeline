# wrapper class for 3dClusterize

wrapper class for 3dClusterize

wrapper class for 3dClusterize

## Active bindings

- `sided`:

  Whether to clusterize 1-sided ('one'), two-sided ('two'), or bi-sided
  ('bi')

- `NN`:

  the cluster definition basis: 1 = faces touch; 2 = edges touch; 3 =
  corners touch

- `clust_nvox`:

  The minimum number of voxels required for each cluster

- `clust_vol`:

  The minimum volume (in microliters) required for each cluster

- `lower_thresh`:

  For two/bi-sided clusterizing, the lower threshold for the left tail
  of the distribution

- `upper_thresh`:

  For two/bi-sided clusterizing, the upper threshold for the right tail
  of the distribution

- `one_thresh`:

  For one-sided clusterizing, the threshold for the test statistic
  distribution

- `mask`:

  the mask within which 3dClusterize searches for clusters

- `pref_map`:

  The name/location of the pref_map (aka cluster_map) file containing an
  integer-valued mask of identified clusters

- `pref_dat`:

  The name/location of the pref_data (aka cluster_masked_data) file
  containing the input data masked by the clusters

- `clusterize_output_file`:

  The name of the text file containing the output of 3dClusterize (i.e.,
  the table of clusters)

- `whereami`:

  passthrough access to whereami object if that has been setup

## Methods

### Public methods

- [`afni_3dclusterize$new()`](#method-afni_3dclusterize-new)

- [`afni_3dclusterize$run()`](#method-afni_3dclusterize-run)

- [`afni_3dclusterize$get_clust_df()`](#method-afni_3dclusterize-get_clust_df)

- [`afni_3dclusterize$get_cluster_map_nifti()`](#method-afni_3dclusterize-get_cluster_map_nifti)

- [`afni_3dclusterize$get_call()`](#method-afni_3dclusterize-get_call)

- [`afni_3dclusterize$get_outputs()`](#method-afni_3dclusterize-get_outputs)

- [`afni_3dclusterize$get_orient()`](#method-afni_3dclusterize-get_orient)

- [`afni_3dclusterize$add_whereami()`](#method-afni_3dclusterize-add_whereami)

- [`afni_3dclusterize$is_complete()`](#method-afni_3dclusterize-is_complete)

- [`afni_3dclusterize$generate_subclusters()`](#method-afni_3dclusterize-generate_subclusters)

- [`afni_3dclusterize$run_subclustering()`](#method-afni_3dclusterize-run_subclustering)

- [`afni_3dclusterize$has_clusters()`](#method-afni_3dclusterize-has_clusters)

- [`afni_3dclusterize$reset_cache()`](#method-afni_3dclusterize-reset_cache)

- [`afni_3dclusterize$get_subclust_list()`](#method-afni_3dclusterize-get_subclust_list)

- [`afni_3dclusterize$subset_atlas_against_clusters()`](#method-afni_3dclusterize-subset_atlas_against_clusters)

- [`afni_3dclusterize$delete_outputs()`](#method-afni_3dclusterize-delete_outputs)

- [`afni_3dclusterize$clone()`](#method-afni_3dclusterize-clone)

------------------------------------------------------------------------

### Method `new()`

initialization function for a new afni_3dclusterize object. Arguments
largely mirror the 3dClusterize parameters.

#### Usage

    afni_3dclusterize$new(
      inset = NULL,
      mask = NULL,
      threshold_file = NULL,
      data_file = NULL,
      mask_from_hdr = NULL,
      out_mask = NULL,
      ithr = NULL,
      idat = NULL,
      onesided = NULL,
      twosided = NULL,
      bisided = NULL,
      lower_thresh = NULL,
      upper_thresh = NULL,
      one_thresh = NULL,
      one_tail = NULL,
      NN = NULL,
      clust_nvox = NULL,
      clust_vol = NULL,
      pref_map = "default",
      pref_dat = "default",
      quiet = NULL,
      orient = NULL,
      binary = NULL,
      clusterize_output_file = NULL
    )

#### Arguments

- `inset`:

  A 4D dataset containing the statistic to use for thresholding (ithr)
  and, optionally, the data value to output/retain

- `mask`:

  If specified, the volume will be masked by `mask` prior to
  clusterizing

- `threshold_file`:

  A 3D dataset containing the statistic to use for thresholding.
  Mutually exclusive with `inset` If passed, `ithr` and `idat` are
  ignored because the `inset` file is generated internally.

- `data_file`:

  A 3D dataset containing the data value to be retained in clusters
  post-thresholding. Must be passed with `threshold_file` and will be
  stitched together with it internally. Mutually exclusive with `inset`.

- `mask_from_hdr`:

  passes through as -mask_from_hdr

- `out_mask`:

  passes through as -out_mask

- `ithr`:

  sub-brik number for the voxelwise threshold. Passes through as -ithr

- `idat`:

  sub-brik number for the voxelwise data to be output in cluster table.
  Passes through as -idat

- `onesided`:

  if TRUE, clusterizing will be conducted on one tail of the statistic
  distribution (-ithr)

- `twosided`:

  if TRUE, clusterizing will be conducted on both tails of the statistic
  distribution (-ithr)

- `bisided`:

  if TRUE, clusterizing will be conducted on each tail of the
  distribution individually

- `lower_thresh`:

  the lower tail cutoff for two/bi-sided testing

- `upper_thresh`:

  the upper tail cutoff for two/bi-sided testing

- `one_thresh`:

  The threshold value for one-sided testing

- `one_tail`:

  For one-sided clusterizing, whether to threshold the LEFT or RIGHT
  tail of the distribution

- `NN`:

  1, 2, 3. Default: 1. Passes through as -NN.

- `clust_nvox`:

  The minimum number of voxels allowed in a cluster. Passes through as
  -clust_nvol

- `clust_vol`:

  The minimum volume in (microliters) allowed in a cluster (mutually
  exclusive with clust_nvox). Passes through as -clust_vol

- `pref_map`:

  File name for the integer-valued mask containing each cluster, ordered
  by descending voxel size. Passes through as -pref_map.

- `pref_dat`:

  File name for the clusterized and thresholded data. Passes through as
  -pref_dat.

- `quiet`:

  passes through as -quiet.

- `orient`:

  passes through as -orient. 'RAI' or 'LPI'. Default is LPI.

- `binary`:

  if TRUE, the pref_map (cluster mask) will be output as a 1/0 binary
  image instead of integer-valued. Passes through as -binary.

- `clusterize_output_file`:

  The name/location of the 3dClusterize output file containing a table
  of identified clusters. Defaults to adding the suffix '\_clusters.1D'
  to the input image and placing the file in the same folder as the
  input.

------------------------------------------------------------------------

### Method `run()`

run the 3dClusterize command relevant to this object

#### Usage

    afni_3dclusterize$run(force = FALSE, quiet = FALSE)

#### Arguments

- `force`:

  if TRUE, 3dClusterize will be re-run

- `quiet`:

  if TRUE, don't output messages as object is run or checked

------------------------------------------------------------------------

### Method `get_clust_df()`

return the 3dClusterize table of clusters as a data.frame

#### Usage

    afni_3dclusterize$get_clust_df(
      include_whereami = TRUE,
      include_subclusters = TRUE,
      include_overlap = TRUE
    )

#### Arguments

- `include_whereami`:

  If TRUE and if \$add_whereami() is already complete, merge the
  whereami data into the cluster data.frame that is returned by this
  function.

- `include_subclusters`:

  If TRUE and if \$generate_subclusters() is already complete, merge the
  subcluster data into the cluster data.frame that is returned by this
  function.

- `include_overlap`:

  If TRUE and if \$add_whereami() is already complete, merge the mask
  overlap data into the cluster data.frame as a nested list-column
  called overlap (since there are many rows for each ROI)

#### Details

This function will return an empty data.frame if the 3dClusterize output
file cannot be found.

------------------------------------------------------------------------

### Method `get_cluster_map_nifti()`

method to read and return the integer-valued clusterized mask (aka
-pref_map) as an oro.nifti object

#### Usage

    afni_3dclusterize$get_cluster_map_nifti()

#### Returns

an oro.nifti object containing the clusterized mask from 3dClusterize

------------------------------------------------------------------------

### Method `get_call()`

returns the 3dClusterize call for this specification

#### Usage

    afni_3dclusterize$get_call()

------------------------------------------------------------------------

### Method `get_outputs()`

Provides a vector of expected output files that correspond to this
3dClusterize setup

#### Usage

    afni_3dclusterize$get_outputs(exclude_missing = TRUE)

#### Arguments

- `exclude_missing`:

  if TRUE (default), any output file that cannot be found will be
  returned as NA.

#### Returns

a named vector of output files related to this 3dClusterize setup

------------------------------------------------------------------------

### Method `get_orient()`

returns the orientation code for this 3dClusterize call (LPI or RAI)

#### Usage

    afni_3dclusterize$get_orient()

------------------------------------------------------------------------

### Method `add_whereami()`

Add's an afni_whereami object to this class in the \$whereami slot. The
corresponding whereami command is also run when this is added so that
coordinates and labels can be obtained immediately. To access the
whereami object and its methods, use \$whereami()

#### Usage

    afni_3dclusterize$add_whereami(atlases = NULL)

#### Arguments

- `atlases`:

  An optional character vector of atlases to be requested in whereami.

------------------------------------------------------------------------

### Method `is_complete()`

returns TRUE if all expected output files exist for this 3dClusterize
call

#### Usage

    afni_3dclusterize$is_complete()

------------------------------------------------------------------------

### Method `generate_subclusters()`

break up large clusters into subclusters

#### Usage

    afni_3dclusterize$generate_subclusters(
      break_nvox = 400,
      min_subclust_nvox = 25,
      max_subclust_nvox = NULL,
      min_n_subclust = 2,
      max_n_subclust = NULL,
      step_size = 0.1,
      max_iter = 50,
      add_whereami = TRUE,
      whereami_atlases = NULL,
      print_progress = FALSE
    )

#### Arguments

- `break_nvox`:

  Break up any clusters larger than this value into subclusters.
  Default: 400

- `min_subclust_nvox`:

  The smallest number of voxels allowed in a subcluster. Default: 25.

- `max_subclust_nvox`:

  The largest numver of voxels allowed in a subcluster. If NULL, no
  upper limit is set.

- `min_n_subclust`:

  The smallest number of subclusters that will be allowed. Must be 2 or
  greater. Default: 2

- `max_n_subclust`:

  The maximum number of subclusters that will be allowed. If NULL, no
  upper limit is set.

- `step_size`:

  The step size used to change the threshold values in the test
  statistic map being clusterized. Default: 0.1.

- `max_iter`:

  The maximum number of steps to be taken for subcluster search.
  Default: 50.

- `add_whereami`:

  If TRUE, whereami will be run for each subcluster. Default: TRUE

- `whereami_atlases`:

  Passes through to afni_whereami for specifying which atlases to use in
  lookup

- `print_progress`:

  If TRUE, the user will see the thresholds being used to subcluster
  each region.

------------------------------------------------------------------------

### Method `run_subclustering()`

runs a subclustering algorithm on this object, increasing the thresholds
until the desired constraints are satisfied

#### Usage

    afni_3dclusterize$run_subclustering(
      min_clust = NULL,
      max_clust = NULL,
      min_nvox = NULL,
      max_nvox = NULL,
      step_size = NULL,
      refine_steps = 5,
      max_iter = 50,
      print_progress = TRUE
    )

#### Arguments

- `min_clust`:

  The minimum number of clusters that will be accepted

- `max_clust`:

  The maximum number of clusters that will be accepted

- `min_nvox`:

  The minimum number of voxels in a subcluster that will be accepted

- `max_nvox`:

  The maximum number of voxels in a subcluster that will be accepted

- `step_size`:

  The increments in the threshold values from one step to the next.

- `refine_steps`:

  The number of steps backward from a winning solution. This maximizes
  the subcluster sizes. Default: 5.

- `max_iter`:

  The maximum number of increment steps that will be taken before giving
  up. Default: 50.

- `print_progress`:

  If TRUE, the user will see the thresholds being used to subcluster
  each region.

#### Details

this is intended to be used internally

------------------------------------------------------------------------

### Method `has_clusters()`

check whether clusteres were found

#### Usage

    afni_3dclusterize$has_clusters()

#### Details

returns TRUE if clusters were found, FALSE if they were not found, and
NULL if the expected cluster output file does not exist (e.g., if
3dClusterize has not been run yet)

------------------------------------------------------------------------

### Method `reset_cache()`

not intended to be called by user, this resets the cluster data.frame
and whereami objects to NULL

#### Usage

    afni_3dclusterize$reset_cache()

#### Details

this is used internally when cloning the parent clusterize object for
subclustering

------------------------------------------------------------------------

### Method `get_subclust_list()`

return list of subcluster details for each large ROI broken up by
generate_subclusters()

#### Usage

    afni_3dclusterize$get_subclust_list()

------------------------------------------------------------------------

### Method `subset_atlas_against_clusters()`

Compares the clusters generated by 3dClusterize to an atlas of interest,
then returns the subset of the atlas that overlaps sufficiently with the
map from 3dClusterize

#### Usage

    afni_3dclusterize$subset_atlas_against_clusters(
      atlas_file = NULL,
      atlas_lower_threshold = 0,
      atlas_upper_threshold = Inf,
      minimum_overlap = 0.8,
      mask_by_overlap = FALSE,
      output_atlas = "default",
      roi_stats = c("mean", "max", "min")
    )

#### Arguments

- `atlas_file`:

  A NIfTI file containing parcels or perhaps meta-analytic statistics

- `atlas_lower_threshold`:

  Only retain values greater than this threshold in the comparison
  against the clusters. Default: 0

- `atlas_upper_threshold`:

  Only retain values less than this threshold in the comparison against
  the clusters. Default: Inf (retain all ROIs above the lower threshold)

- `minimum_overlap`:

  The proportion overlap of an atlas parcel with a cluster required for
  the parcel to be retained. Default: 0.8. If an integer \> 1 is
  provided, then the function treats this as the minimum number of
  \*voxels\* that must overlap between the parcel and the data-derived
  cluster.

- `mask_by_overlap`:

  If TRUE, only voxels in the atlas that overlapped with a cluster are
  retained. In essence, this erodes the retained atlas parcels to only
  include voxels that were in a cluster. Default: FALSE

- `output_atlas`:

  the name/location of the file to output containing the subset of
  parcels in `atlas_file` that are retained by this function. If
  `'default'` or `TRUE`, the subset atlas will be placed in the same
  folder as the 3dClusterize input image, with a filename that combines
  the atlas file name with the input/threshold image name. To disable
  creation of this file, set `output_atlas = FALSE`.

- `roi_stats`:

  A character vector of summary stats to extract from the cluster-masked
  data file. Supported values: `"mean"`, `"min"`, `"max"`. Set to NULL
  to skip extraction.

------------------------------------------------------------------------

### Method `delete_outputs()`

method to delete any/all files generated by this object

#### Usage

    afni_3dclusterize$delete_outputs(prompt = FALSE)

#### Arguments

- `prompt`:

  if TRUE, user will have to confirm deletion of each file. If FALSE,
  files are deleted without prompting.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    afni_3dclusterize$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
