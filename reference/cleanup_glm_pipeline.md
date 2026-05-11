# helper function to refresh l3 model status and save gpa object from batch pipeline back to its cache

helper function to refresh l3 model status and save gpa object from
batch pipeline back to its cache

## Usage

``` r
cleanup_glm_pipeline(gpa, backend_cache_paths = NULL)
```

## Arguments

- gpa:

  a glm_pipeline_arguments object

- backend_cache_paths:

  optional named character vector of backend cache files to merge

## Value

a refreshed version of the gpa object
