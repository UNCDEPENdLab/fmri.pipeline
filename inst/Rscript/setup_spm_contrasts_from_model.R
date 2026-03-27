#!/usr/bin/env Rscript

# This script generates SPM contrasts from an fmri.pipeline L1 model contrast matrix.

args <- commandArgs(trailingOnly = FALSE)

scriptpath <- dirname(sub("--file=", "", grep("--file=", args, fixed = TRUE, value = TRUE), fixed = TRUE))
argpos <- grep("--args", args, fixed = TRUE)
if (length(argpos) > 0L) {
  args <- args[(argpos + 1):length(args)]
} else {
  args <- c()
}

matfile <- NULL
contrast_rds <- NULL
spm_path <- NULL
average_across_runs <- TRUE

argpos <- 1
while (argpos <= length(args)) {
  if (args[argpos] == "-mat_file") {
    matfile <- args[argpos + 1]
    stopifnot(file.exists(matfile))
    argpos <- argpos + 2
  } else if (args[argpos] == "-contrast_rds") {
    contrast_rds <- args[argpos + 1]
    stopifnot(file.exists(contrast_rds))
    argpos <- argpos + 2
  } else if (args[argpos] == "-spm_path") {
    spm_path <- args[argpos + 1]
    stopifnot(dir.exists(spm_path))
    argpos <- argpos + 2
  } else if (args[argpos] == "-average_across_runs") {
    average_across_runs <- as.logical(args[argpos + 1])
    argpos <- argpos + 2
  } else {
    stop("Not sure what to do with argument: ", args[argpos])
  }
}

if (is.null(matfile) || is.null(contrast_rds)) {
  stop("Both -mat_file and -contrast_rds are required.")
}

library(R.matlab)

spec <- readRDS(contrast_rds)
cmat <- spec$contrast_matrix
if (is.null(cmat) || !inherits(cmat, "matrix") || nrow(cmat) == 0L) {
  stop("contrast_matrix must be a non-empty matrix in contrast_rds.")
}

regressors <- colnames(cmat)
if (is.null(regressors) || length(regressors) != ncol(cmat)) {
  stop("contrast_matrix must have column names matching regressor names.")
}

spm_info <- readMat(matfile)
mnames <- unlist(spm_info$mnames)
cpos <- as.vector(spm_info$cpos)

output_dir <- normalizePath(dirname(matfile))

spm_preamble <- c(
  ifelse(is.null(spm_path), "", paste0("addpath('", spm_path, "');")),
  "spm('defaults', 'fmri');",
  "spm_jobman('initcfg');",
  ""
)

# Only consider columns of interest
col_names <- mnames[cpos]
col_clean <- sub("Sn\\(\\d+\\)\\s+", "", col_names, perl = TRUE)
col_has_bf <- grepl("\\*bf\\(\\d+\\)$", col_clean)
col_base <- sub("\\*bf\\(\\d+\\)$", "", col_clean, perl = TRUE)

# Parse condition and pmod names from SPM labels.
# SPM can encode pmods as either:
# 1) "cond x pmod-pmod^1" (with spaces), or
# 2) "condxpmod-pmod^1" (no spaces).
parse_cond_pmod <- function(x) {
  out <- list(cond = x, pmod = NA_character_)

  if (grepl("\\sx\\s", x, perl = TRUE)) {
    out$cond <- trimws(sub("\\sx\\s.*", "", x, perl = TRUE))
    out$pmod <- trimws(sub(".*\\sx\\s", "", x, perl = TRUE))
  } else {
    parsed <- regexec("^(.*)x(.*-pmod(?:\\^\\d+)?)$", x, perl = TRUE)
    parts <- regmatches(x, parsed)[[1]]
    if (length(parts) == 3L) {
      out$cond <- trimws(parts[2])
      out$pmod <- trimws(parts[3])
    }
  }

  out$pmod <- sub("\\^\\d+$", "", out$pmod, perl = TRUE)
  out
}

parsed_labels <- lapply(col_base, parse_cond_pmod)
cond_name <- vapply(parsed_labels, function(x) x$cond, character(1))
pmod_name <- vapply(parsed_labels, function(x) x$pmod, character(1))

match_columns <- function(reg) {
  idx <- which(
    col_base == reg |
      cond_name == reg |
      pmod_name == reg |
      pmod_name == paste0(reg, "-pmod") |
      col_base == paste0(reg, "-pmod")
  )
  if (length(idx) == 0L) return(integer(0))

  # Prefer canonical basis function if present
  idx_bf1 <- idx[grepl("\\*bf\\(1\\)$", col_clean[idx])]
  if (length(idx_bf1) > 0L) idx <- idx_bf1
  idx
}

full_cmat <- matrix(0, nrow = nrow(cmat), ncol = length(mnames))
rownames(full_cmat) <- rownames(cmat)
colnames(full_cmat) <- mnames

missing_regs <- character(0)
for (ii in seq_len(nrow(cmat))) {
  for (jj in seq_len(ncol(cmat))) {
    w <- cmat[ii, jj]
    if (is.na(w) || w == 0) next
    reg <- regressors[jj]
    idx <- match_columns(reg)
    if (length(idx) == 0L) {
      missing_regs <- union(missing_regs, reg)
      next
    }
    if (isTRUE(average_across_runs)) {
      w <- w / length(idx)
    }
    full_cmat[ii, cpos[idx]] <- full_cmat[ii, cpos[idx]] + w
  }
}

if (length(missing_regs) > 0L) {
  warning("No SPM columns matched these regressors: ", paste(missing_regs, collapse = ", "))
}

# Drop contrasts with all zero weights
nonzero_rows <- apply(full_cmat, 1, function(x) any(abs(x) > 1e-12))
if (!any(nonzero_rows)) {
  stop("All contrasts are empty after mapping to SPM columns.")
}
full_cmat <- full_cmat[nonzero_rows, , drop = FALSE]

cn <- 1
m_string <- c(
  spm_preamble,
  "% SETUP BATCH JOB STRUCTURE",
  "contrast = struct;",
  "% spmmat",
  paste0("contrast.matlabbatch{1}.spm.stats.con.spmmat = { ['", output_dir, "' filesep 'SPM.mat']};"),
  ""
)

for (ii in seq_len(nrow(full_cmat))) {
  cname <- rownames(full_cmat)[ii]
  cname <- gsub("'", "", cname, fixed = TRUE)
  m_string <- c(
    m_string,
    paste0("% consess ", cn),
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.tcon.name = '", cname, "';"),
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.tcon.convec = [ ", paste(full_cmat[ii, ], collapse = ", "), " ];"),
    paste0("contrast.matlabbatch{1}.spm.stats.con.consess{", cn, "}.tcon.sessrep = 'none';")
  )
  cn <- cn + 1
}

m_string <- c(
  m_string,
  "",
  "% delete?",
  "contrast.matlabbatch{1}.spm.stats.con.delete = 0;",
  "% RUN BATCH JOB",
  "spm_jobman('run',contrast.matlabbatch);"
)

cat(m_string, file = file.path(output_dir, "estimate_glm_contrasts.m"), sep = "\n")
