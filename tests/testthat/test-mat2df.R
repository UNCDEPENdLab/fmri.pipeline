test_that("mat2df converts numeric matrices like matrix melt", {
  mat <- matrix(as.numeric(seq_len(6)), nrow = 2)
  row_key <- c(1L, 2L, 1L, 2L, 1L, 2L)
  col_key <- c(1L, 1L, 2L, 2L, 3L, 3L)
  signal <- as.numeric(seq_len(6))

  expect_identical(
    fmri.pipeline:::mat2df(mat, value_name = "signal", varnames = c("vnum", "volume")),
    data.frame(
      vnum = row_key,
      volume = col_key,
      signal = signal
    )
  )

  rownamed <- mat
  rownames(rownamed) <- c("row_a", "row_b")
  expect_identical(
    fmri.pipeline:::mat2df(rownamed, value_name = "signal", varnames = c("vnum", "volume")),
    data.frame(
      vnum = factor(c("row_a", "row_b", "row_a", "row_b", "row_a", "row_b"),
        levels = c("row_a", "row_b")),
      volume = col_key,
      signal = signal
    )
  )

  colnamed <- mat
  colnames(colnamed) <- c("col_a", "col_b", "col_c")
  expect_identical(
    fmri.pipeline:::mat2df(colnamed, value_name = "signal", varnames = c("vnum", "volume")),
    data.frame(
      vnum = row_key,
      volume = factor(c("col_a", "col_a", "col_b", "col_b", "col_c", "col_c"),
        levels = c("col_a", "col_b", "col_c")),
      signal = signal
    )
  )

  fully_named <- mat
  dimnames(fully_named) <- list(c("row_a", "row_b"), c("col_a", "col_b", "col_c"))
  expect_identical(
    fmri.pipeline:::mat2df(fully_named, value_name = "signal", varnames = c("vnum", "volume")),
    data.frame(
      vnum = factor(c("row_a", "row_b", "row_a", "row_b", "row_a", "row_b"),
        levels = c("row_a", "row_b")),
      volume = factor(c("col_a", "col_a", "col_b", "col_b", "col_c", "col_c"),
        levels = c("col_a", "col_b", "col_c")),
      signal = signal
    )
  )
})
