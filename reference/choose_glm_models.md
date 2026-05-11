# helper function to ask user to choose models at a given level for further processing

helper function to ask user to choose models at a given level for
further processing

## Usage

``` r
choose_glm_models(gpa, model_names, level, lg = NULL)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` with models specified at a given level

- model_names:

  a user-specified string of models to process/use

- level:

  the level of GLM analysis to be specified (1, 2, or 3)

## Value

a character vector of user-specified model names at this level
