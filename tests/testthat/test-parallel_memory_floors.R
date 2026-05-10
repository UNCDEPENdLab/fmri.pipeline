test_that("setup_parallel_settings enforces memory floors for estimation jobs", {
  gpa <- list(
    parallel = list(
      l1_setup_cores = 1L,
      l1_setup_memgb = "8G",
      fsl = list(
        l1_feat_memgb = "8",
        l2_feat_memgb = "12",
        l3_feat_memgb = "4"
      ),
      spm = list(
        l1_spm_memgb = "16",
        l3_spm_memgb = "24"
      ),
      compute_environment = list()
    ),
    l1_models = list(models = NULL),
    scheduler = "slurm",
    nodename = "unit-test-node"
  )

  out <- setup_parallel_settings(gpa, lg = lgr::get_logger("test/mem_floors"))

  expect_identical(out$parallel$l1_setup_memgb, "16G")
  expect_identical(out$parallel$fsl$l1_feat_memgb, "16")
  expect_identical(out$parallel$fsl$l2_feat_memgb, "16")
  expect_identical(out$parallel$fsl$l3_feat_memgb, "16")
  expect_identical(out$parallel$spm$l1_spm_memgb, "32")
  expect_identical(out$parallel$spm$l3_spm_memgb, "32")
})
