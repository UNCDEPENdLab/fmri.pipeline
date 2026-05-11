# Helper function to obtain the number of voxels in a 4D file for populating FEAT FSF

Helper function to obtain the number of voxels in a 4D file for
populating FEAT FSF

## Usage

``` r
lookup_nvoxels(nifti)
```

## Arguments

- nifti:

  a 4D nifti file

## Value

the number of voxels in `nifti` as calculated by x \* y \* z \* t
