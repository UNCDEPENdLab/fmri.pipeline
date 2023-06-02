# allow a set of .feat folders as input?

#' Generate a set of diagnostic images and plots for brain masks relative to one another and a template
#' @param input Can be a character vector of mask files, a \code{gpa} object, or a \code{data.frame} containing
#'   a \code{brain_mask} field, or at least a \code{$run_nifti} field. If a \code{gpa} object is provided, we will
#'   use the \code{$run_data} \code{data.frame}.
#' @param reference_mask a string pointing to the reference mask used for run brain mask comparisons. Must be in the same
#'   stereotaxic space and have the same dimensions as the run brain masks.
#' @param reference_t1 a string pointing to the T1w image of the reference brain. Used to aid interpretability of
#'   brain mask visualizations.
#' @param output_directory a string pointing to a directory for mask diagnostics. Defaults to "mask_diagnostics" in 
#'   the current directory. This directory will be created if it does not exist.
#' @param generate_automask a logical indicating whether to use AFNI 3dAutomask to generate brain masks if they are
#'   not provided in the $brain_mask field. If NULL (default), the user will be prompted for whether they wish to
#'   generate missing brain masks. If \code{FALSE}, masks will not be generate, and if \code{TRUE}, they will be
#'   generated.   
#'
#' @importFrom RNifti readNifti niftiHeader
#' @importFrom glue glue
#' @importFrom ggplot2 aes scale_fill_viridis_c
#' @importFrom patchwork plot_layout plot_annotation
#' @importFrom abind abind
#' @import ggbrain
#' @export
generate_mask_diagnostics <- function(input = NULL, reference_mask = NULL, reference_t1 = NULL, 
                                      output_directory = NULL, generate_automask = NULL, generate_run_plots = TRUE,
                                      ncores = 8L) {

  if (inherits(input, "glm_pipeline_arguments")) input <- input$run_data
  
  if (inherits(input, "data.frame")) {
    mask_files <- input$brain_mask
    if (is.null(mask_files) && is.character(input$run_nifti)) {
      # attempt fmriprep-friendly search for masks
      message("Attempting to search for fmriprep-friendly brain_mask NIfTI files based on $run_nifti")
      mask_files <- sub("(.*sub-.*_space.*)-preproc_.*\\.nii\\.gz", "\\1-brain_mask.nii.gz", mask_files, perl=TRUE)
    }
  } else if (checkmate::test_character(input)) {
    mask_files <- input
  }
  
  if (is.null(output_directory)) output_directory <- file.path(getwd(), "mask_diagnostics")
  if (!dir.exists(output_directory)) dir.create(output_directory)
  
  intersect_file <- file.path(output_directory, "intersect_mask.nii.gz")
  union_file <- file.path(output_directory, "union_mask.nii.gz")
  proportion_file <- file.path(output_directory, "proportion_mask.nii.gz")

  if (is.null(reference_mask)) {
    reference_mask <- intersect_file # use intersection mask as reference if none provided
  } else {
    checkmate::assert_string(reference_mask)
    checkmate::assert_file_exists(reference_mask)
  }
  
  # for plotting
  if (is.null(reference_t1)) reference_t1 <- reference_mask
  else checkmate::assert_file_exists(reference_t1)
  
  f_exists <- file.exists(mask_files)
  if (all(!f_exists)) {
    stop("None of the brain masks provided exist")
  } else if (any(!f_exists)) {
    miss_files <- mask_files[!f_exists]
    mask_files <- mask_files[f_exists]
    warning("Could not find mask files for some inputs")
    cat(miss_files, sep="\n")
  }
  
  feat_dirs <- grepl("\\.feat$", mask_files)
  if (any(feat_dirs)) {
    # pull out masks from .feat structure
    mask_files[feat_dirs] <- sapply(mask_files[feat_dirs], function(ff) {
      finfo <- read_feat_dir(ff)$mask_file
    })
  }
  
  # generate 4D mask
  # m4d <- tempfile(pattern="mask4d")
  # mcmd <- glue("fslmerge -t {m4d} {paste(mask_files, collapse=' ')}")
  # run_fsl_command(mcmd)
  # m4d <- RNifti::readNifti(m4d)
  
  # this is much faster than calling fslmerge
  m4d <- do.call(abind, list(along=4, lapply(mask_files, function(x) { readNifti(x) })))
  
  # take means along first 3 dimensions to get proportions
  prop_img <- rowMeans(m4d, dims=3) #this is way faster than apply
  #prop_img <- apply(m4d, c(1,2,3), mean))
  
  # rather than taking min/max in an apply (which is slow), use the sums image
  sum_img <- rowSums(m4d, dims=3)
  intersect_img <- 1L*(sum_img == dim(m4d)[4L])
  union_img <- 1L*(sum_img > 0) # 1L*logical converts to integer
  
  ni3d <- niftiHeader(mask_files[1])
  # ni3d$dim <- c(3, ni3d$dim[2:4], 1, 1, 1, 1) # only necessary if converting from 4d
  
  union_nim <- asNifti(union_img, ni3d)
  writeNifti(union_nim, file = union_file)
  
  intersect_nim <- asNifti(intersect_img, ni3d)
  writeNifti(intersect_nim, file = intersect_file)
  
  prop_nim <- asNifti(prop_img, ni3d)
  writeNifti(prop_nim, file = proportion_file)
  
  # Colors: These derive from ColorBrewer palette "Dark2"
  reference_color <- "#7570b3" # slate blue
  intersect_color <- "#D95F02" # orange
  union_color <-     "#1B9E77" # green
  #mask_color <-      "#E7298A" # magenta
  mask_color <-      "#FFFF99" # yellowish
  
  # plot subject versus reference
  intersect_gg <- ggbrain() +
    images(c(underlay=reference_t1)) +
    images(c(intersect=intersect_file)) +
    images(c(union=union_file)) +
    images(c(reference=reference_mask)) +
    geom_brain("underlay") +
    geom_brain(definition="conj := union = union > 0; intersection = intersect > 0",
               mapping = aes(fill = label), alpha=0.7,
               fill_scale = scale_fill_brewer("Group mask", palette="Dark2")) +
    #geom_outline(definition="ref := reference == 1L", outline="blue", show_legend = TRUE) +
    geom_outline(definition="reference", outline=reference_color, alpha=0.8) +
    #geom_brain(definition="aa := reference == 1L") +
    slices(montage("z", n=5, min = 0.2, max=0.8)) +
    slices(montage("y", n=5, min = 0.2, max=0.8)) +
    slices(montage("x", n=5, min = 0.2, max=0.8)) +
    annotate_coordinates(x="q90") + render() +
    plot_annotation(
      title="Group intersection and union masks",
      subtitle = "Outlines: slate blue = reference",
      theme = theme(plot.title=element_text(color="white"), plot.subtitle = element_text(color="white"))
    ) +
    plot_layout(nrow=3)
    
  ggsave(file.path(output_directory, "group_mask.png"), plot = intersect_gg, width = 12, height=9)
  
  # proportion plot
  proportion_gg <- ggbrain() +
    images(c(underlay=reference_t1)) +
    images(c(intersect=intersect_file)) +
    images(c(union=union_file)) +
    images(c(reference=reference_mask)) +
    images(c(prop = proportion_file)) +
    geom_brain("underlay") +
    geom_brain("Proportion := prop") + #, fill_scale = scale_fill_viridis_c("Proportion", breaks=c(0, 0.5, 1))) +
    geom_outline(definition="reference", outline=reference_color, alpha=0.8) +
    geom_outline(definition="union", outline=union_color, alpha=0.8) +
    #geom_brain(definition="aa := reference == 1L") +
    slices(montage("z", n=5, min = 0.2, max=0.8)) +
    slices(montage("y", n=5, min = 0.2, max=0.8)) +
    slices(montage("x", n=5, min = 0.2, max=0.8)) +
    annotate_coordinates(x="q90") +
    render() +
    plot_annotation(
      title = "Proportion of runs with voxel",
      subtitle = "Outlines: slate blue = reference; green = union",
      #caption = "caption",
      theme = theme(plot.title=element_text(color="white"), plot.subtitle = element_text(color="white"))
    ) +
    plot_layout(nrow=3)

  ggsave(file.path(output_directory, "proportion_mask.png"), plot = proportion_gg, width = 12, height=9)
  
  ## Calculate statistics
  ref_img <- readNifti(reference_mask)
  imstats <- apply(m4d, 4, function(im) {
    c(
      mask_voxels = sum(im),
      reference_overlap = sum(im & ref_img),
      intersect_overlap = sum(im & intersect_img),
      reference_only = sum(ref_img[im < .9]),
      run_not_in_reference = sum(im[ref_img < .9]),
      run_not_in_intersect = sum(im[intersect_img < .9])
    )
  }) %>% 
    t() %>% 
    data.frame() %>%
    dplyr::mutate(mask = mask_files) %>%
    dplyr::select(mask, everything())
    
  # now compare each subject against the group and the reference
  if (isTRUE(generate_run_plots)) {
    run_dir <- file.path(output_directory, "by_run")
    if (!dir.exists(run_dir)) dir.create(run_dir)
    parallel::mclapply(seq_along(mask_files), function(ii) {
      run_gg <- ggbrain() +
        images(c(underlay = reference_t1)) +
        images(c(intersect = intersect_file)) +
        images(c(union = union_file)) +
        images(c(reference = reference_mask)) +
        images(c(run = mask_files[ii])) +
        geom_brain("underlay") +
        # geom_brain(definition="conj := union = union > 0; intersection = intersect > 0",
        #            mapping = aes(fill = label), alpha=0.7,
        #            fill_scale = scale_fill_brewer("Mask", palette="Dark2")) +
        #geom_outline(definition="ref := reference == 1L", outline="blue", show_legend = TRUE) +
        geom_brain(definition="run", fill=mask_color, alpha=0.5) +
        geom_outline(definition="reference", outline=reference_color, alpha=0.8) +
        geom_outline(definition="union", outline=union_color, alpha=0.8) +
        geom_outline(definition="intersect", outline=intersect_color, alpha=0.8) +
        slices(montage("z", n=5, min = 0.2, max=0.8)) +
        slices(montage("y", n=5, min = 0.2, max=0.8)) +
        slices(montage("x", n=5, min = 0.2, max=0.8)) +
        annotate_coordinates(x="q90") +
        render() +
        plot_annotation(
          title = glue::glue("Dataset: {basename(mask_files[ii])}"),
          subtitle = "Outlines: slate blue = reference; green = union; orange = intersection",
          theme = theme(plot.title=element_text(color="white"), plot.subtitle = element_text(color="white"))
        ) +
        plot_layout(nrow=3)
      
      ggsave(file.path(run_dir, paste0(fmri.pipeline:::file_sans_ext(basename(mask_files[ii])), ".png")), plot = run_gg, width = 12, height=9)
    }, mc.cores = ncores)
  }
  
  return(imstats)
  
}

# flexible function when last dim is not the one to average over
# means.along <- function(a, i) {
#   n <- length(dim(a))
#   b <- aperm(a, c(seq_len(n)[-i], i))
#   rowMeans(b, dims = n - 1)
# }
# system.time(prop_img2 <- means.along(m4d, 4))
