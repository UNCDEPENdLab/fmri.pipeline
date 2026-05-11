# function to extract voxelwise betas from a given level of FSL analysis

function to extract voxelwise betas from a given level of FSL analysis

## Usage

``` r
extract_fsl_betas(
  gpa,
  extract = NULL,
  level = NULL,
  what = c("cope", "zstat"),
  aggregate = TRUE,
  aggFUN = mean,
  remove_zeros = TRUE,
  mask_file = NULL,
  ncores = 1L,
  scheduler = "local",
  lg = NULL
)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object

- level:

  which level to extract: 1, 2, 3

- what:

  which elements of the FEAT output should be extracted. Default is cope
  and zstat

- aggregate:

  whether to aggregate voxels within each parcel

- aggFUN:

  the function to use for aggregating voxels within a parcel

- remove_zeros:

  whether to remove statistics that are very close to zero (and may thus
  represent masked-out voxels)

- mask_file:

  a nifti image containing integer values for each parcel from which we
  should extract statistics

- ncores:

  the number of cores to use for image extraction. Default: 1

- lg:

  a Logger object for logging beta extraction
