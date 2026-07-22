#' @keywords internal 
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @useDynLib fmri.pipeline, .registration = TRUE
## usethis namespace: end
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ":=", ".", "..split_on", ".filename", ".I", ".SD",
    "afni_complete", "atlas_value", "boundary", "child_level",
    "cope_df",
    "dd", "dt_split", "emm", "emt", "event_list",
    "evt_onset", "evt_time", "expect_complete",
    "expect_func_dir", "expect_func_file", "expect_mprage_bet",
    "expect_mprage_dir", "expect_mprage_file", "expect_mprage_warpcoef",
    "feat_complete", "feat_dir", "feat_dir_exists", "feat_fsf",
    "fname", "Freq", "grouping_scope",
    "i", "id", "ii", "img", "img_chunk", "img_exists", "img_stats", "is_complete",
    "isi", "iti_post", "l1_model", "l2_model", "l3_model", "label",
    "lg", "mask", "mask_value",
    "model", "model_matrix",
    "model_name", "motion_params_file", "name", "outcome", "p.value",
    "pct_true", "rhs",
    "run", "run_nifti", "run_nifti_present", "run_number", "sequence_id",
    "session", "si", "spm_complete", "STAT", "term", "time_started", "to_run", "trial", "val_n",
    "value", "vnum", "volume", "vox_ts", "x", "y", "z"
  ))
}
