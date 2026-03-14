library(testthat)
library(dplyr)

test_that("build_3dlmer_command constructs valid strings", {
  cmd <- fmri.pipeline:::build_3dlmer_command(
    prefix = "myset",
    model_formula = "Group*Session + (1|Subj)",
    qVars = "Age",
    glt_codes = list(con1 = "Group : 1*patient -1*control"),
    data_table_file = "dt.txt",
    mask = "mask.nii.gz",
    njobs = 4
  )
  
  expect_match(cmd, "3dLMEr")
  expect_match(cmd, "-prefix myset")
  expect_match(cmd, "-model 'Group\\*Session \\+ \\(1\\|Subj\\)'")
  expect_match(cmd, "-qVars 'Age'")
  expect_match(cmd, "con1 'Group : 1\\*patient -1\\*control'")
  expect_match(cmd, "-dataTable @dt.txt")
  expect_match(cmd, "-mask mask.nii.gz")
  expect_match(cmd, "-jobs 4")
})

test_that("build_3dlmer_datatable handles factor merging", {
  dt <- data.frame(id = c("s1", "s2"), session = 1, InputFile = c("f1.nii.gz", "f2.nii.gz"), stringsAsFactors = FALSE)
  subject_data <- data.frame(id = c("s1", "s2"), session = 1, Group = c("A", "B"), Age = c(20, 30), stringsAsFactors = FALSE)
  
  res <- fmri.pipeline:::build_3dlmer_datatable(subject_data, dt, model_variables = c("Group", "Age"))
  
  expect_equal(nrow(res), 2)
  expect_true(all(c("Group", "Age") %in% names(res)))
  # In modern R, these are character vectors by default
  expect_equal(as.character(res$Group), c("A", "B"))
})

test_that("emmeans_to_3dlmer_glt translates correctly", {
  # Mock a model object
  data <- data.frame(
    Group = factor(rep(c("A", "B"), each = 5)),
    Y = rnorm(10)
  )
  m <- lm(Y ~ Group, data = data)
  
  mobj <- list(
    lmfit = m,
    contrast_spec = list(
      cond_means = "Group",
      weights = "equal"
    )
  )
  class(mobj) <- "hi_model_spec"
  
  glts <- fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, data)
  
  expect_true("Group.A" %in% names(glts))
  expect_true("Group.B" %in% names(glts))
  expect_match(glts$Group.A, "Group : 1\\*A")
  expect_match(glts$Group.B, "Group : 1\\*B")
})
