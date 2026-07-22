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
    "abs_correlation", "afni_complete", "atlas_value", "boundary", "child_level",
    "cmd", "commands_per_cpus", "contrast_df", "cope", "cope_df", "correlation",
    "dd", "dt_split", "duration", "emm", "emt", "event_list", "events_file",
    "evt_onset", "evt_time", "expect_complete",
    "expect_func_dir", "expect_func_file", "expect_mprage_bet",
    "expect_mprage_dir", "expect_mprage_file", "expect_mprage_warpcoef",
    "feat_complete", "feat_dir", "feat_dir_exists", "feat_fsf",
    "fname", "Freq", "grouping_scope", "high_correlation", "high_vif",
    "i", "id", "ii", "img", "img_chunk", "img_exists", "img_stats", "is_complete",
    "isi", "iti_post", "j", "k", "l1_model", "l2_model", "l3_model", "label",
    "lg", "mask", "mask_value",
    "mean_abs_cor", "mean_cor", "mean_vif", "model", "model_matrix",
    "model_name", "motion_params_file", "n_flagged", "n_runs", "name", "onset", "outcome", "p.value",
    "pct_true", "regressor", "regressor1", "regressor2", "rhs", "roi",
    "run", "run_nifti", "run_nifti_present", "run_number", "sequence_id",
    "session", "si", "spm_complete", "STAT", "term", "time_started", "to_run", "trial", "val_n",
    "value", "vnum", "volume", "vox_ts", "x", "y", "z"
  ))
}
