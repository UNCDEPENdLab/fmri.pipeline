# Add base logger appenders for top-level setup logs This is for backwards compatibility now that logs are moved into a new folder for new gpa objects

Add base logger appenders for top-level setup logs This is for backwards
compatibility now that logs are moved into a new folder for new gpa
objects

## Usage

``` r
add_base_logger_appenders(
  lg,
  gpa,
  log_txt_path,
  log_json_path,
  txt_appender_name,
  json_appender_name
)
```

## Arguments

- lg:

  A logger instance

- gpa:

  glm_pipeline_arguments object containing logging settings

- log_txt_path:

  Path to the text log file

- log_json_path:

  Path to the json log file

- txt_appender_name:

  Appender name for text logs

- json_appender_name:

  Appender name for json logs
