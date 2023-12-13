test_that("Basic GLM pipeline setup works.", {
	gpa <- get_gpa("testdata")

	finalize_pipeline_configuration(gpa)

	expect_equal(1, 1)
})

# entropy
# intercept