# Get Dimensions of a NIfTI Image

Reads the header of a NIfTI file and returns its image dimensions.

## Arguments

- infile:

  Character string. Path to a valid NIfTI file (e.g., \`.nii\` or
  \`.nii.gz\`).

## Value

A numeric vector containing the dimensions of the image. For a 4D image,
the result will be a vector of length 4: `c(x, y, z, t)`.

## Details

This function uses the RNifti C++ API to efficiently extract the
dimensions (e.g., x, y, z, time) of a NIfTI image without loading the
entire image into memory.

## Examples

``` r
if (FALSE) { # \dontrun{
  dims <- get_nifti_dim("sub-001_task-rest_bold.nii.gz")
} # }
```
