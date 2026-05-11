# internal function to convert a feat level 1 directory to an AFNI brik/head file

internal function to convert a feat level 1 directory to an AFNI
brik/head file

## Usage

``` r
feat_stats_to_brik(
  feat_dir,
  out_filename = "feat_stats",
  what = c("cope_files", "z_files", "varcope_files"),
  label_prefix = NULL
)
```

## Arguments

- feat_dir:

  a .feat directory

- out_filename:

  the directory and filename prefix for the output brik/head file

- what:

  which elements of the feat structure to add to the brik.

- label_prefix:

  an optional character string to add as a prefix to the labels

## Value

a list containing the feat_info for the directory, the full path of the
output image, a vector of the brik names, and a vector of the contrast
names
