#' Function to extract subject-level coefficients from an FSL group analysis
#' 
#' @param gpa a \code{glm_pipeline_arguments} object for which results of GLM analyses are already available
#' @param mask_files vector of filenames for NIfTI-based atlases. The function will loop over each and extract coefficients from
#'   each unique value in each mask.
#' @param what which statistics to extract from each parcel. Default is 'cope' (aka 'beta') and 'zstat'.
#' @param out_dir the directory to which statistics are written as .csv.gz files.
#' @param extract_l1 a character vector of l1 models from which to extract run-level (level 1) statistics. If "none",
#'   l1 statistic extraction will be skipped. If "all", then all l1 models will be extracted.
#' @param extract_l2 a character vector of l2 models from which to extract subject-level (level 2) statistics. If "none",
#'   l2 statistic extraction will be skipped. If "all", then all l2 models will be extracted (you will still be prompted
#'   for which l1 models you want to extract from each l2 model).
#' @param extract_l3 a character vector of l3 models from which to extract group-level (level 3) statistics. If "none",
#'   l3 statistic extraction will be skipped. If "all", then all l3 models will be extracted (you will still be prompted
#'   for which l1 and l2 models you want to extract from each l3 model).
#' @param ncores The number of cores to use for extracting statistics from lower-level imgages. If NULL, lookup value in
#'   gpa$parallel$extract_glm_betas_ncores. If this value is > 1 and scheduler = "local", then mclapply will be used locally to extract.
#' @param scheduler The scheduler to use for extracting statistics. If "local", use lapply/mclapply within the current compute session.
#'   If 'slurm', use doFuture with multiple slurm jobs to extract.
#' @param aggregate whether to take the average (or other central tendency measure) of voxels within a given mask value. This
#'   only pertains to integer-valued masks, not continuous ones.
#' @param aggFUN the function used to aggregate statistics for voxels within a given mask value. Default is mean.
#' @param return_data if TRUE, then extracted statistics will be returned as a list of data.frames with elements l1, l2, and l3.
#' @param write_data if TRUE, then extracted statistics will be written as .csv.gz files to \code{out_dir}.
#' 
#' @importFrom checkmate test_integerish assert_class assert_data_frame assert_character assert_logical
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ foreach
#' @importFrom iterators iter
#' @importFrom future plan
#' @importFrom future.batchtools batchtools_slurm
#' @importFrom doFuture registerDoFuture
#' @importFrom data.table fwrite rbindlist
#' @importFrom glue glue_data
#' @importFrom RNifti readNifti
#' @export
extract_glm_betas_in_mask <- function(gpa, mask_files, what=c("cope", "zstat"), out_dir=getwd(),
  extract_l1="all", extract_l2="all", extract_l3="all",
  ncores = NULL, scheduler = "local", aggregate=TRUE, aggFUN=mean, return_data=TRUE, write_data=TRUE) {

  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(gpa$subject_data)
  checkmate::assert_character(mask_files)
  checkmate::assert_file_exists(mask_files)
  checkmate::assert_string(extract_l1)
  checkmate::assert_string(extract_l2)
  checkmate::assert_string(extract_l3)
  checkmate::assert_logical(return_data, len = 1L)
  checkmate::assert_logical(write_data, len = 1L)
  if (isFALSE(return_data) && isFALSE(write_data)) {
    stop("Nothing will happen if both write_data and return_data are FALSE! Set at least one to TRUE to proceed.")
  }

  l3_expression <- "{out_dir}/L1m-{l1_model}/{mask_name}_{statistic}_l3.csv.gz"
  l2_expression <- "{out_dir}/L1m-{l1_model}/{mask_name}_{statistic}_l2.csv.gz"
  l1_expression <- "{out_dir}/L1m-{l1_model}/{mask_name}_{statistic}_l1.csv.gz"

  if (is.null(ncores)) {
    if (!is.null(gpa$parallel$extract_glm_betas_ncores)) {
      ncores <- gpa$parallel$extract_glm_betas_ncores
    } else {
      ncores <- 1L # default to serial execution
    }
  }

  checkmate::assert_integerish(ncores, lower=1L)

  lg <- lgr::get_logger("glm_pipeline/extract_glm_betas")
  lg$set_threshold(gpa$lgr_threshold)

  # enforce available outputs across GLM levels
  enforce_glms_complete(gpa, level = 1L, lg)
  enforce_glms_complete(gpa, level = 3L, lg)

  # in multi-run setup, an l2_model_setup must be present
  if (isTRUE(gpa$multi_run)) enforce_glms_complete(gpa, level = 2L, lg)

  # ask user to specify the models from which betas should be extracted
  extract_model_spec <- build_beta_extraction(gpa, extract_l1, extract_l2, extract_l3, lg)

  # helper function to write out data.frames to csv files for each element in a list (split by key factors like contrast)
  write_list <- function(data_split) {
    for (ff in seq_along(data_split)) {
      out_file <- names(data_split)[ff]
      if (!dir.exists(dirname(out_file))) dir.create(dirname(out_file), recursive=TRUE)
      data.table::fwrite(data_split[[ff]], file=names(data_split)[ff])
    }
    return(names(data_split))
  }

  # subfunction to process and write betas from a given level
  extract_and_write <- function(level=1L) {
    to_extract <- switch(level,
      `1` = extract_l1,
      `2` = extract_l2,
      `3` = extract_l3
    )

    lev_name <- switch(level,
      `1` = "run-level",
      `2` = "subject-level",
      `3` = "group-level"
    )

    espec <- switch(level,
      `1` = extract_model_spec$l1,
      `2` = extract_model_spec$l2,
      `3` = extract_model_spec$l3
    )

    ospec <- switch(level,
      `1` = l1_expression,
      `2` = l2_expression,
      `3` = l3_expression
    )

    res_list <- list()

    # if nothing to extract, return empty list
    # a bit too accommodating on mixed types, but don't want to fix up build_beta_extraction
    if (is.null(espec) || (is.data.frame(espec) && nrow(espec) == 0L) || (is.list(espec) && length(espec) == 0L))  {
      return(res_list)
    }

    for (aa in mask_files) {
      lg$info("Extracting statistics from mask: %s", aa)
      lg$info("Extracting %s statistics for model combinations:", lev_name)
      lg$info("%s", capture.output(print(espec)))
      res_list[[aa]] <- extract_fsl_betas(gpa,
        extract = espec, level = level, what = what,
        aggregate = aggregate, aggFUN = aggFUN, mask_file = aa, ncores = ncores, scheduler = scheduler, lg = lg
      )
    }

    res_list <- rbindlist(res_list)
    if (isTRUE(write_data)) {
      res_split <- res_list %>%
        mutate(out_file = glue_data(., !!ospec)) %>%
        split(by = "out_file", keep.by = FALSE) %>%
        write_list()
    }

    return(res_list)
  }

  l1_df <- extract_and_write(1L)
  l2_df <- extract_and_write(2L)
  l3_df <- extract_and_write(3L)

  if (isTRUE(return_data)) {
    return(list(l1 = l1_df, l2 = l2_df, l3 = l3_df))
  } else {
    return(invisible(NULL))
  }

  # if (ncpus > 1L) {
  #   cl <- makeCluster(ncpus, type="FORK")
  #   registerDoParallel(cl)
  #   on.exit(try(stopCluster(cl))) #cleanup pool upon exit of this function
  # } else {
  #   registerDoSEQ()
  # }

}

#' function to extract voxelwise betas from a given level of FSL analysis
#' @param gpa a \code{glm_pipeline_arguments} object
#' @param level which level to extract: 1, 2, 3
#' @param what which elements of the FEAT output should be extracted. Default is cope and zstat
#' @param aggregate whether to aggregate voxels within each parcel
#' @param aggFUN the function to use for aggregating voxels within a parcel
#' @param remove_zeros whether to remove statistics that are very close to zero (and may thus represent masked-out voxels)
#' @param mask_file a nifti image containing integer values for each parcel from which we should extract statistics
#' @param ncores the number of cores to use for image extraction. Default: 1
#' @param lg a Logger object for logging beta extraction
#'
#' @keywords internal
#' @importFrom oro.nifti readNIfTI translateCoordinate
#' @importFrom checkmate assert_class assert_integerish assert_subset assert_logical assert_file_exists
#' @importFrom dplyr group_by summarize ungroup left_join inner_join bind_cols filter n
#' @importFrom tidyr pivot_longer unnest
#' @importFrom glue glue_data
#' @importFrom tidyselect all_of
#' @importFrom data.table copy setDT set
#' @importFrom RNifti readNifti
#' @importFrom parallel mclapply
extract_fsl_betas <- function(gpa, extract=NULL, level=NULL, what = c("cope", "zstat"),
  aggregate = TRUE, aggFUN = mean, remove_zeros = TRUE, mask_file = NULL, ncores=1L, scheduler = "local", lg=NULL) {

  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1, upper = 3, len = 1)
  checkmate::assert_subset(what, c("cope", "varcope", "zstat", "tstat"))
  checkmate::assert_logical(aggregate, len = 1)
  checkmate::assert_file_exists(mask_file)
  checkmate::assert_integerish(ncores, lower=1)
  checkmate::assert_class(lg, "Logger")

  if (is.null(extract)) {
    lg$debug("Nothing to extract in extract_fsl_betas.")
    return(NULL)
  } else {
    checkmate::assert_data_frame(extract)
  }

  # helper subfunction to extract voxel statistics from a 3D image within a mask
  get_img_stats <- function(img_name, mask, is_int_mask = TRUE, aggregate = TRUE, aggFUN=mean, remove_zeros = TRUE) {
    # use readNifti from the RNifti package for speed (about 8x faster than oro.nifti)
    img <- readNifti(img_name, internal = TRUE)

    stopifnot(identical(mask$dim, dim(img))) # enforce identical dimensions for mask and target image

    # use data.table for faster aggregate + join operations
    coef_df <- copy(mask$coordinates)
    setDT(coef_df) # convert to data.table without resorting. Danger: setting keys will sort data.table before we use mask indices!
    coef_df[, value := img[mask$indices]]
    if (isTRUE(remove_zeros)) {
      zvals <- which(abs(coef_df$value) < 2 * .Machine$double.eps)
      if (length(zvals) > 0L) {
        data.table::set(coef_df, i = zvals, j = "value", value = NA_real_) # use fast set methods in data.table to set NAs
      }
    }

    if (isTRUE(is_int_mask) && isTRUE(aggregate)) {
      # use na.omit to remove any NAs generated from the remove_zeros step
      coef_df <- merge(
        mask$agg_coordinates,
        coef_df[, .(value = aggFUN(na.omit(value))), by = .(mask_value)],
        by = "mask_value"
      )
    }

    return(coef_df)
  }

  if (level == 1L) {
    to_extract <- extract %>% dplyr::left_join(gpa$l1_model_setup$fsl, by = c("l1_model"))
  } else if (level == 2L) {
    to_extract <- extract %>% dplyr::left_join(gpa$l2_model_setup$fsl, by = c("l1_model", "l2_model"))
  } else if (level == 3L) {
    to_extract <- extract %>% dplyr::left_join(gpa$l3_model_setup$fsl, by = c("l1_model", "l2_model", "l3_model"))
  }

  to_extract <- to_extract %>%
    dplyr::filter(feat_dir_exists == TRUE & feat_complete == TRUE)

  if (nrow(to_extract) == 0L) {
    lg$warn("No valid model outputs to extract for input:")
    lg$warn("%s", capture.output(print(extract)))
    return(NULL)
  }

  extra <- NULL # columns to extract
  if (level == 1L) {
    stat_results <- to_extract %>%
      inner_join(get_l1_cope_df(gpa, extract), by = c("id", "session", "l1_model"))

    for (ww in what) {
      # calculate the expected image location for this contrast and subject based on row values in stat_results data.frame
      stat_results <- stat_results %>%
        mutate("{ww}" := glue_data(., "{feat_dir}/stats/{ww}{l1_cope_number}.nii.gz"))
    }
    extra <- "run_number" # also extract
  } else if (level == 2L) {
    stat_results <- to_extract %>%
      inner_join(get_l2_cope_df(gpa, extract), by = c("id", "session", "l2_model")) %>%
      inner_join(get_l1_cope_df(gpa, extract), by = c("id", "session", "l1_model"))

    for (ww in what) {
      # calculate the expected image location for this contrast and subject based on row values in stat_results data.frame
      stat_results <- stat_results %>%
        mutate("{ww}" := glue_data(., "{feat_dir}/cope{l1_cope_number}.feat/stats/{ww}{l2_cope_number}.nii.gz"))
    }
  } else if (level == 3L) {
    browser()
    stat_results <- to_extract %>%
      inner_join(get_l3_cope_df(gpa, extract), by = c("id", "session", "l3_model")) %>%
      inner_join(get_l2_cope_df(gpa, extract), by = c("id", "session", "l2_model")) %>%
      inner_join(get_l1_cope_df(gpa, extract), by = c("id", "session", "l1_model"))
  }

  stat_results <- stat_results %>%
    dplyr::select(
      id, session, all_of(extra), feat_dir,
      matches("l[1-3]_model"), matches("l[1-3]_cope_number"), matches("l[1-3]_cope_name"), all_of(what)
    )

  mask_img <- oro.nifti::readNIfTI(mask_file, reorient = FALSE)
  if (mask_img@dim_[1L] == 4) {
    lg$warn("Mask is a 4D image. Only the first volume will be used")
    mask_img@.Data <- mask_img[, , , 1]
    mask_img@dim_[1] <- 3 # 3d image
    mask_img@dim_[5] <- 1 # 1 volume
  }

  mask_name <- sub("(\\.nii)*(\\.gz)*$", "", basename(mask_file))

  # find indices of mask that are not close to zero (use 2* to support single floating point precision)
  m_indices <- which(abs(mask_img) > 2 * .Machine$double.eps, arr.ind = TRUE)

  m_vals <- unique(mask_img[m_indices])
  is_int_mask <- all(sapply(m_vals, checkmate::test_integerish))
  m_coordinates <- cbind(m_indices, t(apply(m_indices, 1, function(r) { translateCoordinate(i=r, nim=mask_img, verbose=FALSE) })))
  m_coordinates <- as.data.frame(m_coordinates) %>%
    setNames(c("i", "j", "k", "x", "y", "z")) %>%
    dplyr::mutate(vnum = 1:n(), mask_value = mask_img[m_indices], mask_name = mask_name) %>%
    dplyr::select(vnum, mask_value, everything())
  mask_dim <- dim(mask_img)

  # extract each statistic requested
  stat_results <- stat_results %>%
    tidyr::pivot_longer(cols = all_of(what), names_to = "statistic", values_to = "img") %>%
    dplyr::mutate(img_exists = file.exists(img))

  miss_imgs <- stat_results %>% dplyr::filter(img_exists == FALSE)
  if (nrow(miss_imgs) > 0L) {
    lg$warn("Missing expected images in extract_fsl_betas. These will be omitted in outputs.")
    lg$warn("  Missing image: %s", miss_imgs$img)
  }

  stat_results <- stat_results %>%
    dplyr::filter(img_exists == TRUE)

  # The slow part of beta extraction is the reading of images and calculation of summaries. Use mclapply to help
  # Pre-compute ROI-level center of gravity, rather than computing the means repeatedly for each image
  if (isTRUE(is_int_mask) && isTRUE(aggregate)) {
    agg_coordinates <- m_coordinates %>%
      group_by(mask_value) %>%
      dplyr::summarize(mask_name = mask_name[1], x=mean(x), y=mean(y), z=mean(z), .groups="drop") %>%
      setDT(key = "mask_value")
  } else {
    agg_coordinates <- NULL
  }

  # peek at the dimensions of the first 5 images and compare against the mask
  dim_info <- bind_rows(lapply(stat_results$img[1:min(5, nrow(stat_results))], lookup_dim))
  all_eq <- sapply(dim_info, function(x) all(x == x[1]))
  if (isFALSE(all(all_eq, na.rm = TRUE))) {
    msg <- "Dimensions of first images in extract_glm_betas_in_mask do not match. Cannot continue!"
    lg$error(msg)
    stop(msg)
  }

  img_dim <- dim_info[1, seq_along(mask_dim)] %>% unlist(use.names=FALSE)
  if (!isTRUE(all.equal(img_dim, mask_dim))) {
    msg <- "Cannot proceed with extract_fsl_betas due to a mismatch between stat image dimensions and mask dimensions."
    lg$error(msg)
    lg$error(glue("Mask dimensions: {paste(mask_dim, collapse=', ')}"))
    lg$error(glue("Stat image dimensions: {paste(img_dim, collapse=', ')}"))
    stop(msg)
  }

  if (scheduler == "local") {
    stat_results$img_stats <- mclapply(stat_results$img, function(img) {
      lg$debug("Processing image: %s", img)
      tryCatch(
        get_img_stats(img,
          mask = list(indices = m_indices, coordinates = m_coordinates, agg_coordinates = agg_coordinates, dim = mask_dim),
          is_int_mask = is_int_mask, aggregate = aggregate, aggFUN = aggFUN, remove_zeros = remove_zeros
        ),
        error = function(e) {
          lg$error("Problem extracting statistics from image %s. Error: %s", img, as.character(e))
          return(NULL)
        }
      )
    }, mc.cores = ncores)
  } else {
    # 500 images per job, or divide images into 16 jobs, whichever is smaller.
    # If dividing into 16 jobs yields fewer than 56 images per job, revert to 56 -- otherwise jobs are exceedingly brief.
    chunk_size <- min(500, max(56, nrow(stat_results) / 16))
    cores_per_job <- 8 # use mclapply within job to accelerate
    registerDoFuture() # tell dopar to use future compute mechanism

    # options(future.debug = FALSE, future.progress = FALSE)
    # put cache in scratch directory and reduce polling frequency to every 2 seconds since 0.2 was generating errors
    options(future.wait.interval = 2, future.wait.alpha = 1.05) # , future.progress = TRUE, future.debug=TRUE)

    if (scheduler == "torque") {
      future::plan(
        future.batchtools::batchtools_torque,
        template = system.file("templates/torque-simple.tmpl", package="fmri.pipeline"),
        resources = list(
          pbs_directives=list(
            nodes = glue("1:ppn={cores_per_job}"),
            walltime = hours_to_dhms((30 * chunk_size) / (60*60)), # 30-second request per image (upper bound) -- convert from hours to period
            pmem = "1gb" # 1 GB per core
          ),
          sched_args = gpa$parallel$sched_args, # pass through user scheduler settings for cluster
          compute_environment = gpa$parallel$compute_environment # pass through user compute settings for cluster
        )
      )
    } else if (scheduler == "slurm") {
      future::plan(
        future.batchtools::batchtools_slurm,
        template = "slurm-simple",
        resources = list(
          walltime = 30 * chunk_size, # 30-second request per image (upper bound)
          memory = 1024, # 1 GB per core
          ncpus = cores_per_job, # just needs one CPU within each chunk
          chunks.as.arrayjobs = FALSE
        )
      )
    }

    stat_results$img_stats <- foreach(
      img_chunk = iter(stat_results$img), .combine = c, # since each worker returns a list, concatenate into a bigger list
      .options.future = list(chunk.size = chunk_size), .inorder = TRUE,
      .packages = c("parallel", "data.table", "RNifti"),
      .noexport = "gpa" # avoid large object being sent to workers.
    ) %dopar% {

      # use within-job parallelism
      mclapply(img_chunk, function(img) {
        lg$debug("Processing image: %s", img)
        tryCatch(
          get_img_stats(img,
            mask = list(indices = m_indices, coordinates = m_coordinates, agg_coordinates = agg_coordinates, dim = mask_dim),
            is_int_mask = is_int_mask, aggregate = aggregate, aggFUN = aggFUN, remove_zeros = remove_zeros
          ),
          error = function(e) {
            lg$error("Problem extracting statistics from image %s. Error: %s", img, as.character(e))
            return(NULL)
          }
        )
      }, mc.cores = cores_per_job)
    }
  }

  # unnest statistics for each image
  stat_expand <- stat_results %>%
    unnest(img_stats) %>%
    dplyr::select(-img_exists)

  return(stat_expand)
  #save(stat_expand, file="stat_expand.RData")
}

#' helper function to prompt user to choose combinations of l1, l2, and l3 models for beta extraction.
#' @param gpa a \code{glm_pipeline_arguments} object
#' @param lg a logger object for messages
#' @return a list of data frames (l1, l2, l3) consisting of the chosen combinations of l1, l2, and l3 models
#' @keywords internal
build_beta_extraction <- function(gpa, extract_l1="prompt", extract_l2="prompt", extract_l3="prompt", lg=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_class(lg, "Logger")
  checkmate::assert_string(extract_l1)
  checkmate::assert_string(extract_l2)
  checkmate::assert_string(extract_l3)

  # subfunction for choosing l2 and l1 models for a given
  build_lower_set <- function(gpa, l2_chosen="prompt") {
    lower_extract <- list()

    # N.B. When using expand.grid, stringsAsFactors defaults to TRUE. If you then pass one element of a factor in as a subsetting expression
    # of a named vector, it will evaluate to a number (the as.numeric on the factor), rather than looking up based on name. For example,
    # lower_extract$l2_model[1] will become 1, not its character value in an expression like gpa$l1_cope_names[lower_extract$l2_models[1]].
    # Thus, always use stringsAsFactors=FALSE!

    if (isTRUE(gpa$multi_run)) {
      l2_set <- choose_glm_models(gpa, l2_chosen, level = 2L, lg = lg)
        if (is.null(l2_set)) {
          lg$warn("No models chosen at l2.")
        } else {
          for (l2_model in l2_set) {
            cat(sprintf(
              "\n------\nFor level 2 model %s, choose all level 1 models from which to extract statistics.\n", l2_model
            ))
            l1_set <- choose_glm_models(gpa, "prompt", level = 1L, lg = lg)
            if (is.null(l1_set)) {
              lg$warn("No models chosen at l1 for l2 model %s", l2_model)
            } else {
              lower_extract[[l2_model]] <- expand.grid(l2_model = l2_model, l1_model = l1_set, stringsAsFactors=FALSE)
            }
          }
        }
    } else {
      l1_set <- choose_glm_models(gpa, "prompt", level = 1L, lg = lg)
      if (is.null(l1_set)) {
        lg$warn("No models chosen at l1.")
      } else {
        lower_extract <- list(expand.grid(l1_model = l1_set, stringsAsFactors = FALSE))
      }
    }

    do.call(rbind, lower_extract)
  }


  # placeholder lists for extraction selections at each level
  l3 <- list()
  l2 <- list()
  l1 <- list()

  if (extract_l3 != "none") {
    cat("Specify the combination of models for which you want group-level (level 3) statistics.\n")
    cat("This requires you to specify the combination of level 1, 2, and 3 models that define a given group map.\n")
    l3_set <- choose_glm_models(gpa, extract_l3, level = 3L, lg = lg)

    # work our way down to l2 and l1 for each l3 model
    if (is.null(l3_set)) {
      lg$debug("No level 3 models were chosen in build_beta_extraction. We will not extract group-level coefficients.")
    } else {
      for (l3_model in l3_set) {
        cat(sprintf(
          "\n------\nFor level 3 model %s, choose all %slevel 1 (run) models from which to extract group-level betas.\n",
          l3_model, ifelse(isTRUE(gpa$multi_run), "level 2 (subject) and ", "")
        ))
        lower_choices <- build_lower_set(gpa)
        if (!is.null(lower_choices)) lower_choices$l3_model <- l3_model
        l3[[l3_model]] <- lower_choices
      }

      l3 <- do.call(rbind, l3)
      rownames(l3) <- NULL
    }
  }

  if (isTRUE(gpa$multi_run) && extract_l2 != "none") {
    cat("\n------\nNext, specify the combination of models for which you want subject-level (level 2) statistics.\n")
    cat("This requires you to specify the combination of level 1 and 2 models that define a given subject map.\n")

    l2 <- build_lower_set(gpa, l2_chosen = extract_l2)
    rownames(l2) <- NULL
  }

  if (extract_l1 != "none") {
    if (extract_l1 == "prompt") {
      cat("\n------\nNext, choose models from which you want to extract run-level (level 1) statistics.\n")
    }
    l1 <- expand.grid(l1_model = choose_glm_models(gpa, extract_l1, level = 1L, lg = lg), stringsAsFactors = FALSE)
  }

  return(list(l1=l1, l2=l2, l3=l3))

}