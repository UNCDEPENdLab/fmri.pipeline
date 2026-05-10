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
  expect_true(grepl(paste("-prefix", shQuote("myset")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-model", shQuote("Group*Session+(1|Subj)")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-qVars", shQuote("Age")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-gltCode", "con1", shQuote("Group : 1*patient -1*control")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-dataTable", shQuote("@dt.txt")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-mask", shQuote("mask.nii.gz")), cmd, fixed = TRUE))
  expect_true(grepl("-jobs 4", cmd, fixed = TRUE))
})

test_that("build_3dlmer_command removes whitespace from model formulas for AFNI parsing", {
  cmd <- fmri.pipeline:::build_3dlmer_command(
    prefix = "myset",
    model_formula = " session_label + ( 1 | Subj ) ",
    data_table_file = "dt.txt"
  )

  model_arg <- sub(".*-model '([^']+)'.*", "\\1", cmd)
  expect_equal(model_arg, "session_label+(1|Subj)")
  expect_no_match(model_arg, " ")
})

test_that("validate_3dlmer_formula_datatable checks formula variables against AFNI dataTable", {
  dt <- data.frame(
    Subj = c("s1", "s2"),
    session = c(1L, 2L),
    age = c(20, 30),
    group = c("A", "B"),
    InputFile = c("cope1.nii.gz", "cope2.nii.gz"),
    stringsAsFactors = FALSE
  )

  vars <- fmri.pipeline:::validate_3dlmer_formula_datatable(
    "age + group + session + (1 | Subj)",
    dt
  )
  expect_equal(vars, c("age", "group", "session", "Subj"))
  expect_error(
    fmri.pipeline:::validate_3dlmer_formula_datatable(
      "age + treatment + (1 | Subj)",
      dt,
      context = "test model"
    ),
    "test model.*treatment.*Available dataTable columns"
  )
  expect_error(
    fmri.pipeline:::validate_3dlmer_formula_datatable(
      "age + group + (1 | id)",
      dt
    ),
    "Use 'Subj' instead of 'id'"
  )
})

test_that("validate_3dlmer_glt_table validates variables, levels, qVars, and labels", {
  dt <- data.frame(
    Subj = c("s1", "s2"),
    session_label = c("pre", "post"),
    age = c(20, 30),
    InputFile = c("cope1.nii.gz", "cope2.nii.gz"),
    stringsAsFactors = FALSE
  )

  glts <- fmri.pipeline:::validate_3dlmer_glt_table(
    list(
      "session post vs pre" = "session_label : 1*post -1*pre",
      "age slope" = "age : 1"
    ),
    datatable = dt,
    qVars = "age"
  )

  expect_true(is.data.frame(glts))
  expect_equal(glts$label_afni, c("session_post_vs_pre", "age_slope"))
  expect_true(all(glts$valid))
  expect_equal(glts$source, c("user", "user"))

  expect_error(
    fmri.pipeline:::validate_3dlmer_glt_table(
      list(bad = "session_label : 1*followup -1*pre"),
      datatable = dt
    ),
    "followup"
  )
  expect_error(
    fmri.pipeline:::validate_3dlmer_glt_table(
      list(bad = "age : 1*post"),
      datatable = dt,
      qVars = "age"
    ),
    "Quantitative variable 'age'"
  )
  expect_error(
    fmri.pipeline:::validate_3dlmer_glt_table(
      list(bad = "treatment : 1*A"),
      datatable = dt
    ),
    "No 'variable :' block found"
  )
})

test_that("build_3dlmer_command quotes shell paths and sanitizes GLT names", {
  cmd <- fmri.pipeline:::build_3dlmer_command(
    prefix = "/tmp/lmer outputs/my set",
    model_formula = "Group + (1 | Subj)",
    qVars = c("Age", "Days"),
    glt_codes = list(
      "Group A/B" = "Group : 1*A -1*B",
      "1 follow-up contrast" = "Days : 1",
      "Group A/B" = "Group : -1*A 1*B"
    ),
    data_table_file = "data Table.txt",
    mask = "/tmp/masks/main mask.nii.gz"
  )

  expect_true(grepl(paste("-prefix", shQuote("/tmp/lmer outputs/my set")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-mask", shQuote("/tmp/masks/main mask.nii.gz")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-dataTable", shQuote("@data Table.txt")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-qVars", shQuote("Age,Days")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-gltCode", "Group_A_B", shQuote("Group : 1*A -1*B")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-gltCode", "glt_1_follow_up_contrast", shQuote("Days : 1")), cmd, fixed = TRUE))
  expect_true(grepl(paste("-gltCode", "Group_A_B_1", shQuote("Group : -1*A 1*B")), cmd, fixed = TRUE))
})

test_that("build_3dlmer_datatable handles factor merging", {
  dt <- data.frame(
    id = c("s1", "s1", "s2"),
    session = 1,
    InputFile = c("f1.nii.gz", "f1.nii.gz", "f2.nii.gz"),
    stringsAsFactors = FALSE
  )
  subject_data <- data.frame(
    id = c("s1", "s1", "s2"),
    session = 1,
    Group = c("A", "A", "B"),
    Age = c(20, 20, 30),
    stringsAsFactors = FALSE
  )
  
  res <- fmri.pipeline:::build_3dlmer_datatable(subject_data, dt, model_variables = c("Group", "Age"))
  
  expect_equal(nrow(res), 2)
  expect_true(all(c("Group", "Age") %in% names(res)))
  # In modern R, these are character vectors by default
  expect_equal(as.character(res$Group), c("A", "B"))
})

test_that("emmeans_to_3dlmer_glt translates correctly", {
  data <- data.frame(
    Group = factor(rep(c("A", "B"), each = 5)),
    InputFile = paste0("cope", seq_len(10), ".nii.gz"),
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
  
  expect_true(is.data.frame(glts))
  expect_true("Group_A" %in% glts$label_afni)
  expect_true("Group_B" %in% glts$label_afni)
  expect_match(glts$code[glts$label_afni == "Group_A"], "Group : 1\\*A")
  expect_match(glts$code[glts$label_afni == "Group_B"], "Group : 1\\*B")
  expect_true(all(glts$valid))
})

test_that("emmeans_to_3dlmer_glt supports simple pairwise and raw GLTs but rejects generated interactions", {
  data <- data.frame(
    Group = factor(rep(c("A", "B"), each = 5)),
    Session = factor(rep(c("pre", "post"), 5)),
    Age = rep(c(20, 30), each = 5),
    InputFile = paste0("cope", seq_len(10), ".nii.gz"),
    Y = rnorm(10)
  )
  m <- lm(Y ~ Group * Session + Age, data = data)

  mobj <- list(
    lmfit = m,
    contrast_spec = list(
      pairwise_diffs = "Group",
      weights = "equal"
    ),
    lmer_glt_codes = list("raw age" = "Age : 1")
  )
  class(mobj) <- "hi_model_spec"

  glts <- fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, data, qVars = "Age")
  expect_true("Group_pw_A_vs_B" %in% glts$label_afni)
  expect_true("raw_age" %in% glts$label_afni)
  expect_match(glts$code[glts$label_afni == "Group_pw_A_vs_B"], "Group : 1\\*A -1\\*B")
  expect_identical(glts$source[glts$label_afni == "raw_age"], "user")

  mobj$contrast_spec <- list(
    pairwise_diffs = "Group:Session",
    weights = "equal"
  )
  expect_error(
    fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, data, qVars = "Age"),
    "provide it explicitly with lmer_glt_codes"
  )
})

test_that("emmeans_to_3dlmer_glt does not emit treatment-coded one-sided pairwise GLTs", {
  data <- data.frame(
    session_label = factor(
      rep(c("analgesic", "placebo", "opioid"), each = 4),
      levels = c("analgesic", "placebo", "opioid")
    ),
    InputFile = paste0("cope", seq_len(12), ".nii.gz"),
    Y = rnorm(12)
  )
  m <- lm(Y ~ session_label, data = data)

  mobj <- list(
    lmfit = m,
    contrast_spec = list(
      pairwise_diffs = "session_label",
      weights = "equal"
    )
  )
  class(mobj) <- "hi_model_spec"

  glts <- fmri.pipeline:::emmeans_to_3dlmer_glt(mobj, data)
  analgesic_vs_placebo <- glts$code[glts$label_afni == "session_label_pw_analgesic_vs_placebo"]
  analgesic_vs_opioid <- glts$code[glts$label_afni == "session_label_pw_analgesic_vs_opioid"]

  expect_equal(analgesic_vs_placebo, "session_label : 1*analgesic -1*placebo")
  expect_equal(analgesic_vs_opioid, "session_label : 1*analgesic -1*opioid")
  expect_false(any(grepl("^session_label : -1\\*", glts$code)))
})

test_that("validate_3dlmer_glt_table rejects generated one-sided pairwise GLTs", {
  dt <- data.frame(
    Subj = paste0("s", 1:6),
    session_label = rep(c("analgesic", "placebo", "opioid"), 2),
    InputFile = paste0("cope", seq_len(6), ".nii.gz"),
    stringsAsFactors = FALSE
  )
  bad <- data.frame(
    label_raw = "session_label.pw.analgesic.vs.opioid",
    code = "session_label : -1*opioid",
    source = "generated",
    contrast_type = "pairwise_diff",
    stringsAsFactors = FALSE
  )

  expect_error(
    fmri.pipeline:::validate_3dlmer_glt_table(bad, datatable = dt),
    "must use explicit cell-level weights"
  )
})
