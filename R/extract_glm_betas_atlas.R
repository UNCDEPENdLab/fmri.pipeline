#' Function to extract subject-level coefficients from an FSL group analysis
#' 
#' @param gpa a \code{glm_pipeline_arguments} object for which results of GLM analyses are already available
#' @param mask_files vector of filenames for NIfTI-based atlases. The function will loop over each and extract coefficients from
#'   each unique value in each mask.
#' @param what which statistics to extract from each parcel. Default is cope (aka 'beta'), varcope, and zstat.
#' @param out_dir the directory to which statistics are written as .csv.gz files.
#' @param extract_l1 a character vector of l1 models from which to extract run-level (level 1) statistics. If "none",
#'   l1 statistic extraction will be skipped. If "all", then all l1 models will be extracted.
#' @param extract_l2 a character vector of l2 models from which to extract subject-level (level 2) statistics. If "none",
#'   l2 statistic extraction will be skipped. If "all", then all l2 models will be extracted (you will still be prompted
#'   for which l1 models you want to extract from each l2 model).
#' @param extract_l3 a character vector of l3 models from which to extract group-level (level 3) statistics. If "none",
#'   l3 statistic extraction will be skipped. If "all", then all l3 models will be extracted (you will still be prompted
#'   for which l1 and l2 models you want to extract from each l3 model).#' 
#' @param aggregate whether to take the average (or other central tendency measure) of voxels within a given mask value. This
#'   only pertains to integer-valued masks, not continuous ones.
#' @param aggFUN the function used to aggregate statistics for voxels within a given mask value. Default is mean.
#' 
#' @importFrom checkmate test_integerish assert_class assert_data_frame assert_character assert_logical
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom data.table fwrite rbindlist
#' @importFrom glue glue_data
#' @export
extract_glm_betas_atlas <- function(gpa, mask_files, what=c("cope", "varcope", "zstat"), out_dir=getwd(),
  extract_l1="all", extract_l2="all", extract_l3="all",
  ncores=NULL, aggregate=TRUE, aggFUN=mean) {

  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_data_frame(gpa$subject_data)
  checkmate::assert_character(mask_files)
  checkmate::assert_file_exists(mask_files)
  checkmate::assert_logical(extract_l1)
  checkmate::assert_logical(extract_l2)
  checkmate::assert_logical(extract_l3)

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

  # enforce available outputs across GLM levels
  enforce_glms_complete(gpa, level = 1L, lg)
  enforce_glms_complete(gpa, level = 3L, lg)

  # in multi-run setup, an l2_model_setup must be present
  if (isTRUE(gpa$multi_run)) enforce_glms_complete(gpa, level = 2L, lg)

  # ask user to specify the models from which betas should be extracted
  extract_model_spec <- build_beta_extraction(gpa, extract_l1, extract_l2, extract_l3, lg)

  # double loop over mask_files and extract_spec -- make extract_fsl_betas a single mask image

  l1_list <- list()
  l2_list <- list()
  l3_list <- list()

  for (aa in mask_files) {
    lg$info("Extracting statistics from mask: %s", aa)

    if (isTRUE(extract_l3)) {
      lg$info("Extracting group-level statistics for model combinations:")
      lg$info("%s", capture.output(print(extract_model_spec$l3)))
      l3_list[[aa]] <- extract_fsl_betas(gpa,
        extract = extract_model_spec$l3, level = 3L, what = what,
        aggregate = aggregate, aggFUN = aggFUN, mask_file = aa, lg=lg
      )
    }

    if (isTRUE(extract_l2)) {
      lg$info("Extracting subject-level statistics for model combinations:")
      lg$info("%s", capture.output(print(extract_model_spec$l2)))
      l2_list[[aa]] <- extract_fsl_betas(gpa,
        extract = extract_model_spec$l2, level = 2L, what = what,
        aggregate = aggregate, aggFUN = aggFUN, mask_file = aa, lg=lg
      )
    }

    if (isTRUE(extract_l1)) {
      lg$info("Extracting run-level statistics for model combinations:")
      lg$info("%s", capture.output(print(extract_model_spec$l1)))
      l1_list[[aa]] <- extract_fsl_betas(gpa,
        extract = extract_model_spec$l1, level = 1L, what = what,
        aggregate = aggregate, aggFUN = aggFUN, mask_file = aa, lg=lg
      )
    }
  }

  # sort out files based on glue expressions
  if (isTRUE(extract_l3)) {
    l3_list <- rbindlist(l3_list)
    l3_split <- l3_list %>%
      mutate(out_file = glue_data(., !!l3_expression)) %>%
      split(by = "out_file", keep.by = FALSE)
    for (ff in seq_along(l3_split)) data.table::fwrite(l3_split[[ff]], file = names(l3_split)[ff])
  }

  if (isTRUE(extract_l2)) {
    l2_list <- rbindlist(l2_list)
    l2_split <- l2_list %>%
      mutate(out_file = glue_data(., !!l2_expression)) %>%
      split(by = "out_file", keep.by = FALSE)
    for (ff in seq_along(l2_split)) data.table::fwrite(l2_split[[ff]], file = names(l2_split)[ff])
  }
  
  if (isTRUE(extract_l1)) {
    l1_list <- rbindlist(l1_list)
    l1_split <- l1_list %>%
      mutate(out_file = glue_data(., !!l1_expression)) %>%
      split(by = "out_file", keep.by = FALSE)
    for (ff in seq_along(l1_split)) data.table::fwrite(l1_split[[ff]], file=names(l1_split)[ff])
  }

  return(list(l1 = l1_list, l2 = l2_list, l3 = l3_list))
  
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
#' @param mask_file a nifti image containing integer values for each parcel from which we should extract statistics
#' @keywords internal
#' @importFrom oro.nifti readNIfTI
#' @importFrom checkmate assert_class assert_integerish assert_subset assert_logical assert_file_exists
extract_fsl_betas <- function(gpa, extract=NULL, level=NULL, what = c("cope", "zstat"), aggregate = TRUE, aggFUN = mean, mask_file = NULL, lg=NULL) {
  checkmate::assert_class(gpa, "glm_pipeline_arguments")
  checkmate::assert_integerish(level, lower = 1, upper = 3, len = 1)
  checkmate::assert_subset(what, c("cope", "varcope", "zstat", "tstat"))
  checkmate::assert_logical(aggregate, len = 1)
  checkmate::assert_file_exists(mask_file)
  checkmate::assert_class(lg, "Logger")

  if (is.null(extract)) {
    lg$debug("Nothing to extract in extract_fsl_betas.")
    return(NULL)
  } else {
    checkmate::assert_data_frame(extract)
  }

  # helper subfunction to extract voxel statistics from a 3D image within a mask
  get_img_stats <- function(img_name, mask, aggregate=TRUE, aggFUN=mean) {
    img <- readNIfTI(img_name, reorient = FALSE)
    stopifnot(img@dim_[1] == 3L)
    stopifnot(identical(mask$dim, dim(img))) # enforce identical dimensions for mask and target image
    # note that the direct matrix indexing is about 500x faster than the apply approach. leaving here for record-keeping
    # microbenchmark(
    #   coef_vec = img[m_indices],
    #   coef_mat = apply(m_coordinates[, c("i", "j", "k")], 1, function(r) { img[r["i"], r["j"], r["k"]] })
    # )
    coef_vec <- img[mask$indices]
    coef_df <- mask$coordinates %>% bind_cols(value = coef_vec)
    if (isTRUE(is_int_mask) && isTRUE(aggregate)) {
      coef_df <- coef_df %>%
        group_by(mask_value) %>%
        summarize(mask_name = mask_name[1], x=mean(x), y=mean(y), z=mean(z), value = aggFUN(value)) %>% # get center of gravity for parcels
        ungroup()
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

  if (level == 1L) {
    stat_results <- to_extract %>%
      left_join(get_l1_cope_df(gpa, extract), by = c("id", "session", "l1_model"))

    for (ww in what) {
      # calculate the expected image location for this contrast and subject based on row values in stat_results data.frame
      stat_results <- stat_results %>%
        mutate("{ww}" := glue_data(., "{feat_dir}/stats/{ww}{l1_cope_number}.nii.gz"))
    }
  } else if (level == 2L) {
    stat_results <- to_extract %>%
      left_join(get_l2_cope_df(gpa, extract), by = c("id", "session", "l2_model")) %>%
      left_join(get_l1_cope_df(gpa, extract), by = c("id", "session", "l1_model"))
    
    for (ww in what) {
      # calculate the expected image location for this contrast and subject based on row values in stat_results data.frame
      stat_results <- stat_results %>%
        mutate("{ww}" := glue_data(., "{feat_dir}/cope{l1_cope_number}.feat/stats/{ww}{l2_cope_number}.nii.gz"))
    }
  } else if (level == 3L) {
    browser()
    stat_results <- to_extract %>%
      left_join(get_l3_cope_df(gpa, extract), by = c("id", "session", "l3_model")) %>%
      left_join(get_l2_cope_df(gpa, extract), by = c("id", "session", "l2_model")) %>%
      left_join(get_l1_cope_df(gpa, extract), by = c("id", "session", "l1_model"))
  }

  stat_results <- stat_results %>%
    dplyr::select(id, session, feat_dir, matches("l[1-3]_model"), matches("l[1-3]_cope_number"), matches("l[1-3]_cope_name"), all_of(what))

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
    mutate(vnum = 1:n(), mask_value = mask_img[m_indices], mask_name = mask_name) %>%
    select(vnum, mask_value, everything())

  # extract each statistic requested
  stat_results <- stat_results %>%
    pivot_longer(cols = all_of(what), names_to = "statistic", values_to = "img") %>%
    mutate(img_exists = file.exists(img))

  miss_imgs <- stat_results %>% dplyr::filter(img_exists == FALSE)
  if (nrow(miss_imgs) > 0L) {
    lg$warn("Missing expected images in extract_fsl_betas. These will be omitted in outputs.")
    lg$warn("Image: %s", miss_imgs$img)
  }

  stat_results <- stat_results %>%
      dplyr::filter(img_exists == TRUE)

  stat_results$img_stats <- lapply(stat_results$img, function(img) {
    tryCatch(
      get_img_stats(img,
        mask = list(indices = m_indices, coordinates = m_coordinates, dim = dim(mask_img)),
        aggregate = aggregate, aggFUN = aggFUN
      ),
      error = function(e) {
        lg$error("Problem extracting statistics from image %s. Error: %s", img, as.character(e))
        return(NULL)
      }
    )
  })

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
              lower_extract[[l2_model]] <- expand.grid(l2_model = l2_model, l1_model = l1_set)
            }
          }
        }
    } else {
      l1_set <- choose_glm_models(gpa, "prompt", level = 1L, lg = lg)
      if (is.null(l1_set)) {
        lg$warn("No models chosen at l1.")
      } else {
        lower_extract <- list(expand.grid(l1_model = l1_set))
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
    cat("\n------\nNext, choose models from which you want to extract run-level (level 1) statistics.\n")
    l1 <- expand.grid(l1_model = choose_glm_models(gpa, extract_l1, level = 1L, lg = lg))
  }

  return(list(l1=l1, l2=l2, l3=l3))

}