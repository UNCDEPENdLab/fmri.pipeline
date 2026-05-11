# Wrapper for running an AFNI command safely within R

Wrapper for running an AFNI command safely within R

## Usage

``` r
run_afni_command(
  args,
  afnidir = NULL,
  stdout = NULL,
  stderr = NULL,
  echo = TRUE,
  omp_num_threads = 1L,
  ...
)
```

## Arguments

- args:

  AFNI command string to be run

- afnidir:

  Location of AFNI installation. If NULL, the function will search the
  environment for AFNIDIR or afni on PATH

- stdout:

  File target for redirecting stdout. If NULL, stdout will not be
  written to file.

- stderr:

  File target for redirecting stderr. If NULL, stderr will not be
  written to file.

- echo:

  Whether to print AFNI command to the screen. Default: TRUE

- omp_num_threads:

  sets the number of OpenMP threads used for this AFNI command (if
  supported)

- ...:

  Arguments passed through to the \`system\` command.

## Value

The exit status of the executed AFNI command. 0 for success, non-zero
for failure

## Details

This command ensures that AFNI commands are run in an environment that
is setup correctly

## Author

Michael Hallquist

## Examples

``` r

if (FALSE) { # \dontrun{
run_afni_command("3dcopy test_data copy_data")
} # }
```
