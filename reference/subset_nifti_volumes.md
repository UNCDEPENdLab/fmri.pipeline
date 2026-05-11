# Subset Timepoints from a 4D NIfTI Image

This function keeps or removes specified timepoints (volumes) from a 4D
NIfTI image and saves the resulting image to a new file. Timepoints are
specified using 1-based indexing, consistent with R conventions.

## Arguments

- infile:

  Character string. Path to the input 4D NIfTI file.

- tpts:

  Integer vector. Timepoints (1-based) to remove or keep.

- mode:

  Character string. Either "remove" or "keep".

- outfile:

  Character string. Path to save the output NIfTI file with selected
  volumes.

## Value

None. The function writes a new NIfTI file to `outfile`.

## Details

This function uses the `volumes` argument in RNifti to efficiently read
only the retained timepoints from disk. If all volumes are removed, an
error is thrown. The input image must be 4-dimensional (i.e., include a
time dimension).

## Examples

``` r
if (FALSE) { # \dontrun{
subset_nifti_volumes("input_bold.nii.gz", tpts = c(1, 2, 100),
  mode = "remove", outfile = "trimmed_bold.nii.gz")
} # }
```
