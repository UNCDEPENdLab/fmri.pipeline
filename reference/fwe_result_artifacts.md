# Select artifacts from a collected FWE result

Select artifacts from a collected FWE result

## Usage

``` r
fwe_result_artifacts(
  x,
  artifact_role = NULL,
  fwe_alpha = NULL,
  target_id = NULL,
  l3_cope_name = NULL,
  require_existing = TRUE
)
```

## Arguments

- x:

  an FWE plan, run, result, artifact manifest, or persisted path.

- artifact_role:

  optional artifact roles to retain.

- fwe_alpha:

  optional FWE alpha values to retain.

- target_id, l3_cope_name:

  optional semantic target selectors.

- require_existing:

  if `TRUE`, retain only existing artifacts.

## Value

a filtered artifact-manifest data frame.
