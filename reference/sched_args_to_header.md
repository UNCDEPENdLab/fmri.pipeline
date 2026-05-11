# helper function to convert the \$sched_args field to a vector of directives that can be included dynamically in the header of PBS or SBATCH scripts.

helper function to convert the \$sched_args field to a vector of
directives that can be included dynamically in the header of PBS or
SBATCH scripts.

## Usage

``` r
sched_args_to_header(gpa)
```

## Arguments

- gpa:

  a `glm_pipeline_arguments` object containing the
  \$parallel\$sched_args field

## Value

a character vector where the \$sched_args elements are converted to
corresponding \# SBATCH or \# PBS directives for inclusion in scripts
