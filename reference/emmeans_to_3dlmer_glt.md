# Translator from emmeans-style contrast spec to 3dLMEr -gltCode

Translator from emmeans-style contrast spec to 3dLMEr -gltCode

## Usage

``` r
emmeans_to_3dlmer_glt(
  mobj,
  data,
  qVars = NULL,
  raw_glt_codes = NULL,
  context = NULL
)
```

## Arguments

- mobj:

  high-level model specification object (hi_model_spec)

- data:

  the data.frame used for the model (to check factor levels)

## Value

a list of 3dLMEr gltCode strings
