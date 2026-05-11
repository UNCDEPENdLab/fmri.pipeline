# This is a wrapper around the spm_extract_anatomical_rois.m script in the inst directory

This is a wrapper around the spm_extract_anatomical_rois.m script in the
inst directory

## Usage

``` r
spm_extract_anatomical_rois(
  l1spmdirs,
  masks,
  threshold = 0.2,
  threshdesc = "none",
  session = 1,
  extent = 0,
  adjust_F_index = 1,
  contrast_index = NULL,
  ncores = 1,
  spm_path = "/gpfs/group/mnh5174/default/lab_resources/spm12",
  matlab_path = "/opt/aci/sw/matlab/R2017b/bin"
)
```

## Arguments

- l1spmdirs:

  character vector of level 1 SPM directories containing SPM.mat files

- masks:

  character vector of NIfTI mask images for each anatomical ROI of
  interest

- threshold:

  p-value threshold applied to contrast within mask before extraction

- threshdesc:

  multiple comparisons correction on p-value. 'none' or 'FWE'

- session:

  which session (run) to use for extracting time series

- extent:

  exclude clusters having fewer than voxels than extent

- adjust_F_index:

  index of F-test in SPM.mat to adjust for all effects of interest

- contrast_index:

  index of t-test contrast in SPM.mat that is of interest

- ncores:

  number of cores to use in a parallel approach; function parallelizes
  over `l1spmdirs`

- spm_path:

  path to spm12 installation; added to MATLAB path at runtime

- matlab_path:

  location of MATLAB binary; used with matlabr for run_matlab_code()

## Author

Michael Hallquist
