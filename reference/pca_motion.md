# small helper function for compressing motion parameters using PCA

small helper function for compressing motion parameters using PCA

## Usage

``` r
pca_motion(motion_df, num_pcs = 3L, zscore = TRUE, verbose = FALSE)
```

## Arguments

- motion_df:

  volumes x motion parameters data frame for PCA compression

- num_pcs:

  number of principal components to extract

- zscore:

  whether to standardize motion parameters prior to PCA (this is a good
  idea)

- verbose:

  whether to print the variance explained by PCs and the multiple
  correlations of among motion parameters

## Value

A matrix of PCA-compressed motion regressors
