scheduler_safe_token <- function(x, fallback = "job", max_chars = 120L) {
  if (is.null(x) || length(x) == 0L || is.na(x[1L]) || !nzchar(x[1L])) {
    x <- fallback
  }
  x <- as.character(x[1L])
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) x <- fallback
  substr(x, 1L, max_chars)
}

scheduler_log_file <- function(log_directory, scheduler, job_name = NULL, job_id = "%j", extra = NULL) {
  checkmate::assert_string(log_directory)
  scheduler <- tolower(scheduler %||% "local")
  prefix <- if (scheduler %in% c("slurm", "sbatch")) {
    "slurm"
  } else if (scheduler %in% c("torque", "qsub", "pbs")) {
    "pbs"
  } else {
    "local"
  }

  parts <- c(prefix, job_id, scheduler_safe_token(job_name), scheduler_safe_token(extra, fallback = "", max_chars = 80L))
  parts <- parts[nzchar(parts)]
  file.path(log_directory, paste0(paste(parts, collapse = "-"), ".out"))
}

scheduler_output_directives <- function(scheduler, log_directory, job_name = NULL, extra = NULL) {
  scheduler <- tolower(scheduler %||% "local")
  if (scheduler %in% c("slurm", "sbatch")) {
    log_file <- scheduler_log_file(log_directory, scheduler = "slurm", job_name = job_name, job_id = "%j", extra = extra)
    return(c(
      paste0("#SBATCH --output=", log_file),
      paste0("#SBATCH --error=", log_file)
    ))
  }

  if (scheduler %in% c("torque", "qsub", "pbs")) {
    log_file <- scheduler_log_file(log_directory, scheduler = "pbs", job_name = job_name, job_id = "pending", extra = extra)
    return(c(
      paste("#PBS -o", log_file),
      "#PBS -j oe"
    ))
  }

  character(0)
}

scheduler_child_log_directory <- function(gpa = NULL) {
  candidates <- character(0)

  if (!is.null(gpa)) {
    candidates <- c(
      candidates,
      gpa$batch_run$batch_directory,
      gpa$output_locations$active_batch_directory
    )
  }

  candidates <- c(
    candidates,
    Sys.getenv("FMRI_PIPELINE_BATCH_DIRECTORY", unset = ""),
    Sys.getenv("SLURM_SUBMIT_DIR", unset = ""),
    Sys.getenv("PBS_O_WORKDIR", unset = ""),
    getwd()
  )

  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  candidates[1L]
}

scheduler_runtime_log_assignment <- function(scheduler, log_directory, job_name = NULL, extra = NULL, job_var = "job_id") {
  scheduler <- tolower(scheduler %||% "local")
  scheduler_name <- if (scheduler %in% c("slurm", "sbatch")) {
    "slurm"
  } else if (scheduler %in% c("torque", "qsub", "pbs")) {
    "pbs"
  } else {
    "local"
  }

  pattern <- scheduler_log_file(
    log_directory,
    scheduler = scheduler_name,
    job_name = job_name,
    job_id = paste0("${", job_var, "}"),
    extra = extra
  )
  pattern <- gsub("\\\\", "\\\\\\\\", pattern)
  pattern <- gsub('"', '\\"', pattern, fixed = TRUE)
  c(
    paste0("scheduler_log=\"", pattern, "\""),
    if (scheduler_name == "pbs") {
      "exec > \"$scheduler_log\" 2>&1"
    } else {
      NULL
    }
  )
}

job_manifest_shell_function <- function() {
  c(
    "function write_job_manifest() {",
    "  artifact_dir=\"$1\"",
    "  artifact_type=\"$2\"",
    "  artifact_status=\"$3\"",
    "  artifact_source=\"$4\"",
    "  manifest=\"${artifact_dir}/job_manifest.tsv\"",
    "  mkdir -p \"${artifact_dir}\"",
    "  if [ ! -f \"${manifest}\" ]; then",
    "    printf 'artifact_dir\\tartifact_type\\tstatus\\tjob_id\\tjob_name\\tscheduler\\tscheduler_log\\tbatch_directory\\tbatch_file\\tsource_file\\tsubmitted_at\\n' > \"${manifest}\"",
    "  fi",
    "  printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \"${artifact_dir}\" \"${artifact_type}\" \"${artifact_status}\" \"${job_id}\" \"${job_name}\" \"${scheduler_name}\" \"${scheduler_log}\" \"${batch_directory}\" \"${batch_file}\" \"${artifact_source}\" \"$(date -Is)\" >> \"${manifest}\"",
    "}"
  )
}
