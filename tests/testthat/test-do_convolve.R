fsl_do_convolve_reference <- function(input, kernel, phase = 0L, renorm = 1L) {
  n_input <- length(input)
  n_kernel <- length(kernel)
  output <- rep(0, n_input)

  for (t in seq.int(0L, n_input - 1L)) {
    kernel_norm <- 0
    i_start <- max(0L, 1L + t + phase - n_input)
    i_end <- min(n_kernel, t + phase + 1L)

    if (i_start < i_end) {
      for (i in seq.int(i_start, i_end - 1L)) {
        output[t + 1L] <- output[t + 1L] + input[t + phase - i + 1L] * kernel[i + 1L]
        kernel_norm <- kernel_norm + kernel[i + 1L]
      }
    }

    if (renorm) {
      output[t + 1L] <- output[t + 1L] / kernel_norm
    }
  }

  output
}

fsl_feat_model_path <- function() {
  candidates <- unique(c(
    unname(Sys.which("feat_model")),
    file.path(Sys.getenv("FSLDIR", unset = ""), "bin", "feat_model")
  ))
  candidates <- candidates[nzchar(candidates)]
  candidates[file.exists(candidates)][1L]
}

read_fsl_design_matrix <- function(path) {
  lines <- readLines(path, warn = FALSE)
  matrix_line <- which(lines == "/Matrix")
  stopifnot(length(matrix_line) == 1L)
  utils::read.table(text = paste(lines[(matrix_line + 1L):length(lines)], collapse = "\n"))
}

replace_fsf_setting <- function(fsf, key, value) {
  pattern <- paste0("^set ", gsub("([()\\.])", "\\\\\\1", key), " ")
  replacement <- paste("set", key, value)
  fsf[grepl(pattern, fsf)] <- replacement
  fsf
}

write_do_convolve_feat_fsf <- function(root, timing_file, tr, n_vols, tempfilt = 0L) {
  template <- readLines(system.file("feat_lvl1_nparam_template.fsf", package = "fmri.pipeline"), warn = FALSE)
  template <- gsub(".OUTPUTDIR.", paste0(root, ".feat"), template, fixed = TRUE)
  template <- gsub(".TR.", format(tr, scientific = FALSE), template, fixed = TRUE)
  template <- gsub(".NVOL.", as.character(n_vols), template, fixed = TRUE)
  template <- gsub(".NVOXELS.", "1", template, fixed = TRUE)
  template <- gsub(".FUNCTIONAL.", "", template, fixed = TRUE)
  template <- gsub(".CONFOUNDS.", "", template, fixed = TRUE)
  template <- replace_fsf_setting(template, "fmri(temphp_yn)", "0")
  template <- replace_fsf_setting(template, "fmri(templp_yn)", "0")
  template <- replace_fsf_setting(template, "fmri(confoundevs)", "0")

  ev_syntax <- fsl_generate_fsf_lvl1_ev_syntax(list(
    list(
      name = "event",
      waveform = "custom_3",
      timing_file = timing_file,
      convolution = "double_gamma",
      convolve_phase = 0,
      tempfilt = tempfilt,
      deriv = 0
    )
  ))

  contrast <- matrix(1, nrow = 1L, dimnames = list("event", "event"))
  contrast_syntax <- fsl_generate_fsf_contrast_syntax(contrast)
  writeLines(c(template, "", ev_syntax, "", contrast_syntax), paste0(root, ".fsf"))
}

fsl_double_gamma_kernel <- function() {
  mult <- 1 / 0.05
  sigma1 <- 2.449
  delay1 <- 6
  sigma2 <- 4
  delay2 <- 16
  ratio <- 6
  fw <- as.integer((delay2 + sigma2 * 5) * mult)
  x <- (0:(fw - 1L)) / mult
  fmri.pipeline:::gammapdf(x, delay1, sigma1 * sigma1) -
    fmri.pipeline:::gammapdf(x, delay2, sigma2 * sigma2) / ratio
}

expected_unfiltered_feat_design <- function(events, tr, n_vols) {
  mult <- 1 / 0.05
  trmult <- tr * mult
  negpts <- as.integer(30 * mult)
  mnpts <- as.integer((n_vols + 10) * trmult + negpts)
  high_res <- rep(0, mnpts)

  for (event_i in seq_len(nrow(events))) {
    start <- events$onset[event_i] * mult + negpts
    stop <- start + events$duration[event_i] * mult

    if (stop > 0 && stop - start < 1) {
      stop <- start + 1.1
    }
    if (start < 0) {
      start <- 0
    }
    if (start > mnpts) {
      stop <- 0
    }
    if (stop > mnpts) {
      stop <- mnpts
    }

    if (start < stop) {
      high_res[seq.int(as.integer(start), as.integer(stop) - 1L) + 1L] <- events$value[event_i]
    }
  }

  high_res <- high_res - mean(high_res)
  convolved <- fmri.pipeline:::do_convolve(high_res, fsl_double_gamma_kernel(), phase = 0L, renorm = 1L)
  downsample_index <- as.integer(((0:(n_vols - 1L)) + 0.5) * trmult) + negpts + 1L
  design <- convolved[downsample_index]
  design - mean(design)
}

test_that("do_convolve matches FSL causal convolution without renormalization", {
  input <- 1:5
  kernel <- c(2, 1)

  out <- fmri.pipeline:::do_convolve(input, kernel, phase = 0L, renorm = FALSE)

  expect_equal(as.numeric(out), c(2, 5, 8, 11, 14), tolerance = 1e-12)
})

test_that("do_convolve matches FSL renormalization for the available kernel support", {
  input <- 1:5
  kernel <- c(2, 1)

  out <- fmri.pipeline:::do_convolve(input, kernel, phase = 0L, renorm = TRUE)

  expect_equal(
    as.numeric(out),
    c(2 / 2, 5 / 3, 8 / 3, 11 / 3, 14 / 3),
    tolerance = 1e-12
  )
  expect_equal(fmri.pipeline:::do_convolve(input, kernel), out, tolerance = 1e-12)
})

test_that("do_convolve applies positive FSL phase with end-of-series truncation", {
  input <- 1:5
  kernel <- c(2, 1)

  raw <- fmri.pipeline:::do_convolve(input, kernel, phase = 1L, renorm = FALSE)
  renorm <- fmri.pipeline:::do_convolve(input, kernel, phase = 1L, renorm = TRUE)

  expect_equal(as.numeric(raw), c(5, 8, 11, 14, 5), tolerance = 1e-12)
  expect_equal(
    as.numeric(renorm),
    c(5 / 3, 8 / 3, 11 / 3, 14 / 3, 5 / 1),
    tolerance = 1e-12
  )
})

test_that("do_convolve preserves FSL behavior for negative phase before support", {
  input <- 1:5
  kernel <- c(2, 1)

  raw <- fmri.pipeline:::do_convolve(input, kernel, phase = -2L, renorm = FALSE)
  renorm <- fmri.pipeline:::do_convolve(input, kernel, phase = -2L, renorm = TRUE)

  expect_equal(as.numeric(raw), c(0, 0, 2, 5, 8), tolerance = 1e-12)
  expect_true(all(is.nan(renorm[1:2])))
  expect_equal(as.numeric(renorm[3:5]), c(2 / 2, 5 / 3, 8 / 3), tolerance = 1e-12)
})

test_that("do_convolve matches a direct FSL source translation across randomized inputs", {
  set.seed(2407)

  cases <- expand.grid(
    n_input = c(3L, 5L, 12L, 31L),
    n_kernel = c(1L, 2L, 4L, 9L),
    phase = -5L:5L,
    renorm = c(0L, 1L)
  )

  for (case_i in seq_len(nrow(cases))) {
    n_input <- cases$n_input[case_i]
    n_kernel <- cases$n_kernel[case_i]
    phase <- cases$phase[case_i]
    renorm <- cases$renorm[case_i]
    input <- stats::rnorm(n_input)
    kernel <- stats::runif(n_kernel, min = 0.2, max = 2)

    expected <- fsl_do_convolve_reference(input, kernel, phase = phase, renorm = renorm)
    observed <- fmri.pipeline:::do_convolve(input, kernel, phase = phase, renorm = renorm)

    expect_equal(
      as.numeric(observed),
      expected,
      tolerance = 1e-6,
      info = sprintf(
        "n_input=%d, n_kernel=%d, phase=%d, renorm=%d",
        n_input, n_kernel, phase, renorm
      )
    )
  }
})

test_that("do_convolve reproduces FSL feat_model design.mat before temporal filtering", {
  feat_model <- fsl_feat_model_path()
  skip_if(is.na(feat_model), "FSL feat_model is not available")

  old_fsldir <- Sys.getenv("FSLDIR", unset = NA_character_)
  if (identical(old_fsldir, NA_character_) || !nzchar(old_fsldir)) {
    Sys.setenv(FSLDIR = normalizePath(dirname(dirname(feat_model)), mustWork = FALSE))
  }
  on.exit({
    if (identical(old_fsldir, NA_character_)) {
      Sys.unsetenv("FSLDIR")
    } else {
      Sys.setenv(FSLDIR = old_fsldir)
    }
  }, add = TRUE)

  tmp <- tempfile("do-convolve-feat-")
  dir.create(tmp)
  root <- file.path(tmp, "design")
  timing_file <- file.path(tmp, "events.txt")
  tr <- 2
  n_vols <- 40L
  events <- data.frame(
    onset = c(2, 11, 23.5, 48),
    duration = c(1.2, 2.5, 0.3, 3),
    value = c(1, 0.5, 1.5, -0.75)
  )
  utils::write.table(events, timing_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
  write_do_convolve_feat_fsf(root, timing_file, tr = tr, n_vols = n_vols, tempfilt = 0L)

  status <- system2(feat_model, root, stdout = TRUE, stderr = TRUE)
  expect_true(is.null(attr(status, "status")), info = paste(status, collapse = "\n"))
  expect_true(file.exists(paste0(root, ".mat")), info = paste(status, collapse = "\n"))

  fsl_design <- read_fsl_design_matrix(paste0(root, ".mat"))[[1L]]
  expected <- expected_unfiltered_feat_design(events, tr = tr, n_vols = n_vols)

  expect_equal(fsl_design, as.numeric(expected), tolerance = 1e-5)
})
