# Helper to write 3dLMEr script and data table

Helper to write 3dLMEr script and data table

## Usage

``` r
write_3dlmer_files(output_dir, dt, cmd, script_name = "run_3dlmer.sh")
```

## Arguments

- output_dir:

  directory where files should be written

- dt:

  the data table data.frame

- cmd:

  the 3dLMEr command string

- script_name:

  name of the shell script to write

## Value

list with paths to script and data table
