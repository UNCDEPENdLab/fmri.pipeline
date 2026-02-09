# GLM backend specification and resolution helpers

glm_backend_not_implemented <- function(backend, feature) {
  force(backend)
  force(feature)
  fn <- function(...) {
    stop(sprintf("GLM backend '%s' does not implement '%s' yet.", backend, feature), call. = FALSE)
  }
  attr(fn, "glm_backend_not_implemented") <- TRUE
  fn
}

resolve_backend_fn <- function(fn_name, pkg = "fmri.pipeline") {
  if (is.null(fn_name)) return(NULL)
  checkmate::assert_string(fn_name)
  checkmate::assert_string(pkg)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required namespace '%s' not available for backend resolution.", pkg), call. = FALSE)
  }
  get(fn_name, envir = asNamespace(pkg))
}

default_glm_backend_specs <- function() {
  list(
    fsl = list(
      name = "fsl",
      l1_setup = "fsl_l1_model",
      l2_setup = "fsl_l2_model",
      l3_setup = "fsl_l3_model",
      l1_status = "get_feat_status",
      l2_status = "get_feat_status",
      l3_status = "get_feat_status",
      l1_status_inputs = c("feat_dir", "feat_fsf"),
      l2_status_inputs = c("feat_dir", "feat_fsf"),
      l3_status_inputs = c("feat_dir", "feat_fsf"),
      output_dir = "get_output_directory",
      l1_run = "run_feat_sepjobs",
      l2_run = "run_feat_sepjobs",
      l3_run = "run_feat_sepjobs"
    ),
    spm = list(
      name = "spm",
      l1_setup = "spm_l1_model",
      l2_setup = NULL,
      l3_setup = "spm_l3_model",
      l1_status = "get_spm_status",
      l2_status = NULL,
      l3_status = "get_spm_status",
      l1_status_inputs = c("spm_dir"),
      l2_status_inputs = character(0),
      l3_status_inputs = c("spm_dir"),
      output_dir = "get_output_directory",
      l1_run = "run_spm_sepjobs",
      l2_run = NULL,
      l3_run = "run_spm_sepjobs"
    ),
    afni = list(
      name = "afni",
      l1_setup = "__not_implemented__",
      l2_setup = "__not_implemented__",
      l3_setup = "__not_implemented__",
      l1_status = "__not_implemented__",
      l2_status = "__not_implemented__",
      l3_status = "__not_implemented__",
      l1_status_inputs = character(0),
      l2_status_inputs = character(0),
      l3_status_inputs = character(0),
      output_dir = "get_output_directory",
      l1_run = "__not_implemented__",
      l2_run = "__not_implemented__",
      l3_run = "__not_implemented__"
    )
  )
}

resolve_glm_backends <- function(specs, pkg = "fmri.pipeline") {
  checkmate::assert_list(specs)
  out <- list()
  for (backend_name in names(specs)) {
    spec <- specs[[backend_name]]
    if (is.null(spec) || !is.list(spec)) next
    resolved <- spec
    for (field in c("l1_setup", "l2_setup", "l3_setup", "l1_status", "l2_status", "l3_status", "output_dir", "l1_run", "l2_run", "l3_run")) {
      val <- resolved[[field]]
      if (is.null(val)) next
      if (identical(val, "__not_implemented__")) {
        resolved[[field]] <- glm_backend_not_implemented(backend_name, field)
      } else if (is.character(val)) {
        resolved[[field]] <- resolve_backend_fn(val, pkg = pkg)
      }
    }
    out[[backend_name]] <- resolved
  }
  out
}

get_glm_backends <- function(gpa, must_exist = TRUE) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_flag(must_exist)

  glm_software <- unique(gpa$glm_software)
  specs <- gpa$glm_backend_specs
  if (is.null(specs)) specs <- default_glm_backend_specs()

  # Always resolve lazily from specs to avoid persisting stale function handles
  resolved <- resolve_glm_backends(specs)

  out <- list()
  for (name in glm_software) {
    if (!name %in% names(resolved)) {
      if (isTRUE(must_exist)) {
        stop(sprintf("Unknown GLM backend '%s'. Available: %s", name, paste(names(resolved), collapse = ", ")), call. = FALSE)
      }
      next
    }
    out[[name]] <- resolved[[name]]
  }

  out
}

initialize_glm_backends <- function(gpa) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  if (is.null(gpa$glm_backend_specs)) {
    gpa$glm_backend_specs <- default_glm_backend_specs()
  }
  # Do not persist function handles in the gpa object
  gpa$glm_backends <- NULL
  gpa
}
