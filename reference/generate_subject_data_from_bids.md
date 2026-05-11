# Function to generate a subject_data object from a BIDS-compliant folder

Function to generate a subject_data object from a BIDS-compliant folder

## Usage

``` r
generate_subject_data_from_bids(bids_dir)
```

## Arguments

- bids_dir:

  a directory containing BIDS-compliant processed data for analysis

## Value

a data.frame containing subject-level data contained in participants.tsv

## Details

This function basically just reads participants.tsv from the root of the
BIDS folder and renames \`participant_id“ to \`id\` to match this
pipeline's conventions
