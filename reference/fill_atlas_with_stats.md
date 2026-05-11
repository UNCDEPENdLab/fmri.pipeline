# Populates parcel-wise statistics from an atlas or set of clusters back into NIfTIs in the original space

Populates parcel-wise statistics from an atlas or set of clusters back
into NIfTIs in the original space

## Usage

``` r
fill_atlas_with_stats(
  atlas_nifti,
  stat_dt,
  stat_cols = c("t", "p"),
  stat_labels = NULL,
  atlas_col = "atlas_value",
  split_on = NULL,
  stack_along = NULL,
  img_prefix = "atlasfill",
  out_dir = getwd(),
  overwrite = FALSE,
  afni_dir = "~/abin"
)
```

## Arguments

- atlas_nifti:

  The filename of the NIfTI image containing atlas values to match
  against `stat_dt`

- stat_dt:

  a `data.table` object containing the statistics to write to parcels in
  the atlas

- stat_cols:

  The column names in `stat_dt` that should be written to parcels in the
  NIfTI. By default, each column will be written as a separate sub-brik
  (i.e., along the 4th dimension of the NIfTI), but this can be modified
  by the `stack_along` argument.

- stat_labels:

  A character vector of the same length as `stat_cols` giving the names
  for each statistic that will be used in the output NIfTI files. By
  default the labels are simply the `stat_cols`.

- atlas_col:

  The column name in `stat_dt` containing the integer values that
  correspond to parcel numbers in `atlas_nifti`. This is used to map the
  statistics back into the `atlas_nifti` space.

- split_on:

  An optional character vector that yields an output image for each
  unique combination of splits. For example, if `stat_dt` contains
  multiple GLM contrasts in `"contrast"` and multiple subjects in
  `"subject"`, then `split_on = c("contrast", "subject")` would yield a
  separate NIfTI with the statistics for each contrast and subject. This
  allows the fill operation to be repeated over an arbitrarily stacked
  dataset.

- stack_along:

  By default, this is `NULL`, in which case each output NIfTI will have
  a sub-briks for each statistic in `stat_cols`. Alternatively,
  `stack_along` can be a character vector of other columns in `stat_dt`.
  In this case,

- img_prefix:

  A string prefixing the filenames for all NIfTIs

- out_dir:

  The directory in which to place the stat NIfTIs. If not present, it
  will be created

- overwrite:

  Whether to allow existing files to be overwritten. Default: FALSE

- afni_dir:

  The location of your AFNI installation, which is used to cobble
  together NIfTIs
