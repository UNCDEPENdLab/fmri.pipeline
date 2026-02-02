# tests for merging backend caches into shared gpa

# Helper to create setup objects with proper class
make_l_setup <- function(level, ...) {
  setup <- list(...)
  class(setup) <- c(paste0("l", level, "_setup"), "list")
  setup
}

make_gpa <- function(...) {
 gpa <- list(...)
  class(gpa) <- c("glm_pipeline_arguments", "list")
  gpa
}

test_that("merge_backend_caches merges backend-specific l1 setups", {
  shared_fsl <- data.frame(
    id = "sub1", session = 1L, l1_model = "m1",
    feat_dir = "feat_old", feat_fsf = "feat_old.fsf", feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  spm_df <- data.frame(
    id = "sub1", session = 1L, l1_model = "m1",
    spm_dir = "spm_dir", spm_complete = TRUE,
    stringsAsFactors = FALSE
  )

  shared_gpa <- make_gpa(l1_model_setup = make_l_setup(1, fsl = shared_fsl))
  spm_gpa <- make_gpa(l1_model_setup = make_l_setup(1, spm = spm_df))

  spm_cache <- tempfile(pattern = "spm_cache_", fileext = ".RData")
  gpa <- spm_gpa
  save(gpa, file = spm_cache)

  merged <- fmri.pipeline:::merge_backend_caches(shared_gpa, c(spm = spm_cache))

  expect_true(inherits(merged$l1_model_setup, "l1_setup"))
  expect_equal(merged$l1_model_setup$fsl, shared_fsl)
  expect_equal(merged$l1_model_setup$spm, spm_df)
})

test_that("merge_backend_caches replaces backend data with cache contents", {
  shared_fsl <- data.frame(
    id = "sub1", session = 1L, l1_model = "m1",
    feat_dir = "feat_old", feat_fsf = "feat_old.fsf", feat_complete = TRUE,
    stringsAsFactors = FALSE
  )
  cache_fsl <- data.frame(
    id = "sub1", session = 1L, l1_model = "m1",
    feat_dir = "feat_new", feat_fsf = "feat_new.fsf", feat_complete = TRUE,
    stringsAsFactors = FALSE
  )

  shared_gpa <- make_gpa(l1_model_setup = make_l_setup(1, fsl = shared_fsl))
  fsl_gpa <- make_gpa(l1_model_setup = make_l_setup(1, fsl = cache_fsl))

  fsl_cache <- tempfile(pattern = "fsl_cache_", fileext = ".RData")
  gpa <- fsl_gpa
  save(gpa, file = fsl_cache)

  merged <- fmri.pipeline:::merge_backend_caches(shared_gpa, c(fsl = fsl_cache))
  expect_equal(merged$l1_model_setup$fsl, cache_fsl)
})

test_that("merge_backend_caches handles missing cache file gracefully", {
  shared_fsl <- data.frame(
    id = "sub1", session = 1L, l1_model = "m1",
    feat_dir = "feat_dir", stringsAsFactors = FALSE
  )
  shared_gpa <- make_gpa(l1_model_setup = make_l_setup(1, fsl = shared_fsl))

  # Non-existent cache file
  missing_cache <- tempfile(pattern = "missing_", fileext = ".RData")

  # Should not error, just warn and skip the missing cache
  merged <- fmri.pipeline:::merge_backend_caches(shared_gpa, c(spm = missing_cache))

  # Original data should be preserved
 expect_equal(merged$l1_model_setup$fsl, shared_fsl)
  expect_null(merged$l1_model_setup$spm)
})

test_that("merge_backend_caches merges multiple levels (l1, l2, l3)", {
  # Shared GPA has FSL L1 and L2
  shared_fsl_l1 <- data.frame(id = "sub1", feat_dir = "l1_feat", stringsAsFactors = FALSE)
  shared_fsl_l2 <- data.frame(id = "sub1", feat_dir = "l2_feat", stringsAsFactors = FALSE)

  shared_gpa <- make_gpa(
    l1_model_setup = make_l_setup(1, fsl = shared_fsl_l1),
    l2_model_setup = make_l_setup(2, fsl = shared_fsl_l2)
  )

  # SPM cache has L1 and L3
  spm_l1 <- data.frame(id = "sub1", spm_dir = "l1_spm", stringsAsFactors = FALSE)
  spm_l3 <- data.frame(id = "sub1", spm_dir = "l3_spm", stringsAsFactors = FALSE)

  spm_gpa <- make_gpa(
    l1_model_setup = make_l_setup(1, spm = spm_l1),
    l3_model_setup = make_l_setup(3, spm = spm_l3)
  )

  spm_cache <- tempfile(pattern = "spm_multi_", fileext = ".RData")
  gpa <- spm_gpa
  save(gpa, file = spm_cache)

  merged <- fmri.pipeline:::merge_backend_caches(shared_gpa, c(spm = spm_cache))

  # L1: FSL preserved, SPM merged
  expect_equal(merged$l1_model_setup$fsl, shared_fsl_l1)
  expect_equal(merged$l1_model_setup$spm, spm_l1)

  # L2: FSL preserved, no SPM (SPM doesn't use L2)
  expect_equal(merged$l2_model_setup$fsl, shared_fsl_l2)

  # L3: SPM merged from cache
  expect_equal(merged$l3_model_setup$spm, spm_l3)
})

test_that("merge_backend_caches handles cache without gpa object", {
  shared_fsl <- data.frame(id = "sub1", feat_dir = "feat_dir", stringsAsFactors = FALSE)
  shared_gpa <- make_gpa(l1_model_setup = make_l_setup(1, fsl = shared_fsl))

  # Create cache with wrong object name
  bad_cache <- tempfile(pattern = "bad_cache_", fileext = ".RData")
  wrong_obj <- list(foo = "bar")
  save(wrong_obj, file = bad_cache)

  # Should not error, just skip invalid cache
  merged <- fmri.pipeline:::merge_backend_caches(shared_gpa, c(spm = bad_cache))

  # Original data preserved
  expect_equal(merged$l1_model_setup$fsl, shared_fsl)
})

test_that("merge_backend_caches handles cache with non-gpa object", {
  shared_fsl <- data.frame(id = "sub1", feat_dir = "feat_dir", stringsAsFactors = FALSE)
  shared_gpa <- make_gpa(l1_model_setup = make_l_setup(1, fsl = shared_fsl))

  # Create cache with gpa that's not glm_pipeline_arguments
  bad_cache <- tempfile(pattern = "bad_class_", fileext = ".RData")
  gpa <- list(foo = "bar")  # No class
  save(gpa, file = bad_cache)

  # Should not error, just skip invalid cache
  merged <- fmri.pipeline:::merge_backend_caches(shared_gpa, c(spm = bad_cache))

  # Original data preserved
  expect_equal(merged$l1_model_setup$fsl, shared_fsl)
})

test_that("merge_backend_caches merges both FSL and SPM caches", {
  # Empty shared GPA
  shared_gpa <- make_gpa()

  # FSL cache
  fsl_l1 <- data.frame(id = "sub1", feat_dir = "fsl_feat", stringsAsFactors = FALSE)
  fsl_gpa <- make_gpa(l1_model_setup = make_l_setup(1, fsl = fsl_l1))
  fsl_cache <- tempfile(pattern = "fsl_", fileext = ".RData")
  gpa <- fsl_gpa
  save(gpa, file = fsl_cache)

  # SPM cache
  spm_l1 <- data.frame(id = "sub1", spm_dir = "spm_dir", stringsAsFactors = FALSE)
  spm_gpa <- make_gpa(l1_model_setup = make_l_setup(1, spm = spm_l1))
  spm_cache <- tempfile(pattern = "spm_", fileext = ".RData")
  gpa <- spm_gpa
  save(gpa, file = spm_cache)

  # Merge both
  merged <- fmri.pipeline:::merge_backend_caches(
    shared_gpa,
    c(fsl = fsl_cache, spm = spm_cache)
  )

  expect_equal(merged$l1_model_setup$fsl, fsl_l1)
  expect_equal(merged$l1_model_setup$spm, spm_l1)
})

test_that("merge_backend_caches preserves class of setup objects", {
  shared_gpa <- make_gpa()

  fsl_l1 <- data.frame(id = "sub1", feat_dir = "feat", stringsAsFactors = FALSE)
  fsl_l2 <- data.frame(id = "sub1", feat_dir = "feat_l2", stringsAsFactors = FALSE)
  fsl_l3 <- data.frame(id = "sub1", feat_dir = "feat_l3", stringsAsFactors = FALSE)

  fsl_gpa <- make_gpa(
    l1_model_setup = make_l_setup(1, fsl = fsl_l1),
    l2_model_setup = make_l_setup(2, fsl = fsl_l2),
    l3_model_setup = make_l_setup(3, fsl = fsl_l3)
  )

  fsl_cache <- tempfile(pattern = "fsl_class_", fileext = ".RData")
  gpa <- fsl_gpa
  save(gpa, file = fsl_cache)

  merged <- fmri.pipeline:::merge_backend_caches(shared_gpa, c(fsl = fsl_cache))

  expect_true(inherits(merged$l1_model_setup, "l1_setup"))
  expect_true(inherits(merged$l2_model_setup, "l2_setup"))
  expect_true(inherits(merged$l3_model_setup, "l3_setup"))
})

