ggbrain_images_r6 <- R6::R6Class(
  classname = "ggbrain_images_r6",
  private = list(
    pvt_imgs = list(), # image data
    pvt_dims = c(), # x, y, z extent
    set_images = function(images, divide_overlay_posneg = TRUE) {
      checkmate::assert_character(images)
      checkmate::assert_file_exists(images)
      if (is.null(names(images))) {
        warning(
          "The images vector does not contain any names. ",
          "This may lead to weird behaviors downstream if 'underlay' and 'overlay' are requested."
        )
        names(images) <- paste0("img", seq_along(images))
      } else if (!"underlay" %in% names(images)) {
        warning("'underlay' is not among the images provided. This may lead to weirdness downstream.")
      }

      if (isTRUE(divide_overlay_posneg) && !"overlay" %in% names(images)) {
        stop("divide_overlay_posneg is TRUE, but 'overlay' is not among the images provided")
      }

      private$pvt_imgs <- sapply(images, RNifti::readNifti, simplify = FALSE)

      if (isTRUE(divide_overlay_posneg)) {
        pos <- private$pvt_imgs[["overlay"]]
        neg <- private$pvt_imgs[["overlay"]]
        pos[pos < 0] <- 0
        neg[neg > 0] <- 0
        message("Adding positive overlay values to $images$overlay_pos")
        message("Adding negative overlay values to $images$overlay_neg")
        private$pvt_imgs[["overlay_pos"]] <- pos
        private$pvt_imgs[["overlay_neg"]] <- neg
      }
    }
  ),
  active = list(
    #' @field sided Whether to clusterize 1-sided ('one'), two-sided ('two'), or bi-sided ('bi')
    images = function(val) {
      if (missing(val)) {
        return(private$pvt_imgs)
      } else {
        private$set_images(val)
      }
    }
  ),

  public = list(
    initialize = function(images, divide_overlay_posneg = TRUE) {
      private$set_images(images, divide_overlay_posneg)
    },
    get_all_slices = function(slices) {

    },
    get_slices = function(imgs, slice_numbers, plane, drop=FALSE) {
      if (!checkmate::test_subset(imgs, names(private$pvt_imgs))) {
        stop(glue("The img input to $get_slice() must be one of: {paste(names(private$pvt_imgs), collapse=', ')}"))
      }

      checkmate::assert_integerish(slice_numbers, lower = 1)
      checkmate::assert_subset(plane, c("sagittal", "coronal", "axial"))
      #return named list of slices for the images requested
      sapply(imgs, function(iname) {
        dat <- private$pvt_imgs[[iname]]
        if (plane == "sagittal") {
          slc_mat <- aperm(dat[slice_numbers, , , drop = FALSE], c(1, 2, 3))
        } else if (plane == "coronal") {
          slc_mat <- aperm(dat[, slice_numbers, , drop = FALSE], c(2, 1, 3))
        } else if (plane == "axial") {
          slc_mat <- aperm(dat[, , slice_numbers, drop = FALSE], c(3, 1, 2))
        }

        attr(slc_mat, "slice_numbers") <- slice_numbers
        attr(slc_mat, "plane") <- plane

        if (isTRUE(drop)) slc_mat <- drop(slc_mat)
        return(slc_mat)
      }, simplify = FALSE)

    }
  )

)


xx <- ggbrain_images_r6$new(c(underlay = "template_brain.nii", overlay = "zstat6_ptfce_fwep_0.05_1mm.nii.gz"))


ggbrain_r6 <- R6::R6Class(
  classname = "ggbrain_r6",
  private = list(
    layer_imgs = list(), # keep original data?
    pvt_panels = list(),
    composite_plot = NULL,
    set_panels = function(panels) {
      if (checkmate::assert_class(panels, "gg")) {

      }
      checkmate::assert_list(panels)
      sapply(panels, function(x) { checkmate::assert_class(x, "ggbrain_panel_r6") })
      private$pvt_panels <- panels
    }
  ),
  active = list(
    #' @field sided Whether to clusterize 1-sided ('one'), two-sided ('two'), or bi-sided ('bi')
    panels = function(val) {
      if (missing(val)) {
        return(private$pvt_panels)
      } else {
        private$set_panels(val)
      }
    }
  ),
  public = list(
    initialize = function(panels = NULL) {
      if (is.null(panels)) {
        stop("Cannot create a ggbrain object without panels!")
      } else {
        self$set_panels(panels)
      }
      
    },
    plot = function() {
      plot(private$composite_plot)
    },
    # future idea? -- multiple views based on cached data
    create_view = function(slices) {
      private$views <- c(private$views, "VIEW HERE")
    }
  )
)

# allow for gg + theme() type stuff
`+.ggbrain_r6` <- function(gg, args) {
  gg_new <- gg$clone(deep=TRUE) # need a new object to modify the panels in memory (not by reference)
  gg_new$panels <- lapply(gg_new$get_panels(), function(gg) { gg + args }) # add to each panel
  return(gg_new)
}

# allow for gg + theme() type stuff
`+.ggbrain_panel_r6` <-  function(gg, args) {
  gg_new <- gg$clone(deep=TRUE) # need a new object to modify the panels in memory (not by reference)
  gg_new$ggobj <- gg_new$get() + args # add args to panel
  return(gg_new)
}

ggbrain_panel_r6 <- R6::R6Class(
  classname = "ggbrain_panel_r6",
  private = list(
    pos_limits = c(),
    neg_limits = c(),
    img_df = NULL,
    ggobj = NULL
  ),
  public = list(
    initialize = function(ggobj = NULL, neg_limits = NULL, pos_limits = NULL) {
      checkmate::assert_class(ggobj, "gg")
      private$ggobj <- ggobj
    },
    plot = function(use_global_limits = TRUE) {
      # add enforcement of limits
      plot(private$ggobj)
    },
    get = function() {
      private$ggobj
    },
    set = function(ggobj) {
      checkmate::assert_class(ggobj, "gg")
      private$ggobj <- ggobj
    }
  )
)

plot.ggbrain_r6 <- function(obj) {
  obj$plot()
}




#' Plot fMRI data on an underlay image
#' 
#' @param underlay a 3D nifti image used for the image underlay (default b/w)
#' @param overlay a 4D nifti image used for plotting stats on underlay (color)
#' @param color_col a position in the 4th dim of overlay use to color plots
#' @param alpha_col a position in the 4th dim of overlay use to set alpha transparency of plots
#' @param underlay_colorscale A ggplot scale_fill_* function call used for coloration of underlay
#' @param overlay_colorscale A ggplot scale_fill_* function call used for coloration of overlay
#' @param remove_null_space If TRUE, quantiles are computed on the non-zero slices and plots are trimmed for
#'   empty space.
#' @param pos_thresh The positive threshold to be applied to the \code{overlay} image. Any voxel that exceeds
#'   the 
#' @param zero_underlay Any voxels in the underlay image whose absolute values are less than \code{zero_underlay}
#'   are set to precisely zero. This helps with small color variation in black/empty space.
#' @param trim_underlay Winsorize the extreme values of the underlay based on the low and high quantiles provided
#' @importFrom checkmate assert_numeric
#' @importFrom ggnewscale new_scale_fill
#' @importFrom cowplot plot_grid add_sub
#' @importFrom ggplot2 geom_raster scale_fill_gradient theme
#' @export
ggbrain <- function(underlay=NULL, overlay=NULL, 
                    color_col=NULL, alpha_col=NULL,
                    underlay_colorscale = scale_fill_gradient(low="grey8", high="grey92", na.value = "transparent"),
                    negative_colorscale = scale_fill_distiller(palette="Blues", direction = 1),
                    positive_colorscale = scale_fill_distiller(palette="Reds"),
                    slices = data.frame(
                      coord = c("x = 25%", "x = 50%", "x = 75%",
                                "y = 25%", "y = 50%", "y = 75%",
                                "z = 25%", "z = 50%", "z = 75%")
                    ),
                    remove_null_space=TRUE,
                    pos_thresh = 1,
                    neg_thresh = -1, legend_label = expression(italic("z")),
                    background_color = "gray10", text_color = "white",
                    trim_underlay = c(.01, .99), zero_underlay = 1e-3, symmetric_legend = TRUE,
                    panel_labels = NULL, underlay_contrast = "none", panel_borders = TRUE,
                    theme_custom = NULL, base_size = 14
) {
  
  # slices <- data.frame(
  #   coord = c("x = 55.4", "x = 50%", "y = 10", "k = 20"),
  #   panel_title = c("entropy", "value", "test", "test1"),
  #   xlab = c("xx1", "xx2", "xx3", "xx4"),
  #   ylab = c("gg", "gg", "gg", "gg"),
  #   coord_label = c(TRUE, TRUE, TRUE, FALSE)
  # )
  
  # handle data.frame
  if (is.list(slices)) {
    slices <- bind_rows(slices)
  }
  
  if (is.null(slices$coord_label)) {
    slices$coord_label <- TRUE # default to labeling slices
  }
  
  max_bg <- "gray95" # brightest value on underlay
  
  require(ggnewscale)
  require(cowplot)
  checkmate::assert_file_exists(underlay)
  checkmate::assert_file_exists(overlay)
  checkmate::assert_class(theme_custom, "theme", null.ok = TRUE)
  
  stat_decimals <- 2 # rounding for overlay

  underlay <- RNifti::readNifti(underlay)
  overlay <- RNifti::readNifti(overlay)
  
  # round very small values to zero
  if (!is.null(zero_underlay) && zero_underlay > 0) {
    underlay[underlay > -1*zero_underlay & underlay < zero_underlay] <- 0
  }
  
  # winsorize extreme values in underlay based on quantiles of non-zero voxels
  checkmate::assert_numeric(trim_underlay, lower=0, upper=1, len = 2)
  if (trim_underlay[1] > 0) {
    lthresh <- quantile(underlay[underlay > 0], trim_underlay[1])
    underlay[underlay < lthresh & underlay > 0] <- lthresh
  }
  
  if (trim_underlay[2] < 1) {
    uthresh <- quantile(underlay[underlay > 0], trim_underlay[2])
    underlay[underlay > uthresh] <- uthresh
  }
  
  underlay[abs(underlay) < 1e-8] <- NA
  
  # sigmoid transform
  #underlay2 <- underlay + underlay*underlay_contrast*(1/(1+exp(-underlay)))
  #underlay <- underlay*5*(1/(1+exp(-underlay)))
  #browser()
  
  # Beta of 1 is very weak contrast enhancement, 10 is very steep
  sigmoid <- function(x, beta=1) {
    # make beta (slope) scale invariant by standardizing (note that scale is very slow for array-level normalization)
    m_x <- mean(x, na.rm=T)
    s_x <- sd(x, na.rm=T)
    z_x <- (x - m_x)/s_x
    
    y <- 1/(1+exp(-beta*(z_x)))
    return(y*x)
  }
  
  if (underlay_contrast == "low") {
    beta <- .1
  } else if (underlay_contrast == "medium") {
    beta <- .8
  } else if (underlay_contrast == "high") {
    beta <- 1.6
  }
  
  if (underlay_contrast != "none") {
    underlay <- sigmoid(underlay, beta)  
  }
  
  
  #verify that i,j,k (1,2,3) dimensions of underlay match dimensions of overlay
  stopifnot(identical(dim(underlay)[1:3], dim(overlay)[1:3]))
  
  # need grouping variables to be factors to maintain order in group_by %>% group_split
  coords_df <- lookup_slices(slices$coord, underlay, remove_null_space)
  
  # split into row-wise list for lapply inside slice lookup
  coords <- coords_df %>% group_by(slice_index) %>% group_split()
  
  # helper subfunction to extract specific slices from 3d image
  get_slice <- function(img, slice_number, plane) {
    checkmate::assert_number(slice_number)
    if (plane == "sagittal") {
      slc_mat <- aperm(img[slice_number, , , drop=FALSE], c(1, 2, 3))
    } else if (plane == "coronal") {
      slc_mat <- aperm(img[, slice_number, , drop=FALSE], c(2, 1, 3))
    } else if (plane == "axial") {
      slc_mat <- aperm(img[, , slice_number, drop=FALSE], c(3, 1, 2))
    }
    
    return(slc_mat)
  }
  
  get_all_slices <- function(img) {
    arg_name <- as.character(match.call())[-1]
    slc <- lapply(coords, function(slc) { get_slice(img, slc$slice_number, slc$plane)})
    # add numeric slice number on the front so that the levels of the factor that are created below follow the proper slice order
    names(slc) <- paste0(sprintf("%03d", seq_along(slc)), coords_df$plane, coords_df$slice_number)

    # use the max x and y sizes across slices to create a uniform plot size
    size_x <- max(sapply(slc, function(s) dim(s)[2]))
    size_y <- max(sapply(slc, function(s) dim(s)[3]))
    square_mat <- matrix(NA_real_, nrow=size_x, ncol=size_y)
    square_melt <- reshape2::melt(square_mat, varnames=c("dim1", "dim2"), value.name="dummy")
    
    slc <- dplyr::bind_rows(lapply(seq_along(slc), function(ll) {
      df <- reshape2::melt(slc[[ll]], varnames=c("slice", "dim1", "dim2"))
      dim1_diff <- size_x - max(df$dim1)
      dim2_diff <- size_y - max(df$dim2)

      if (dim1_diff > 1) { df$dim1 <- df$dim1 + round(dim1_diff/2) } # center in x
      if (dim2_diff > 1) { df$dim2 <- df$dim2 + round(dim2_diff/2) } # center in y

      # place image onto shared square matrix
      df <- square_melt %>% left_join(df, by=c("dim1", "dim2"))

      #uniquely identify slices by number and type
      df$slice <- names(slc)[ll]
      return(df)
    }))
    
    slc <- slc %>% mutate(
      image = arg_name,
      slice = factor(slice) # store as factor so that if a slice drops in the filter, we still get an empty list
    )
    return(slc)
  }

  anat_slices <- get_all_slices(underlay) #%>% rename(uval = value) # for fill scale to work properly with overlay, need different variable names
  overlay_slices <- get_all_slices(overlay) #%>% rename(oval = value)

  # separate layers for pos and neg
  pos_plot <- overlay_slices %>% filter(value >= !!pos_thresh) #%>% rename(pos = oval)
  neg_plot <- overlay_slices %>% filter(value <= !!neg_thresh) #%>% rename(neg = oval)

  has_pos <- nrow(pos_plot) > 0L
  has_neg <- nrow(neg_plot) > 0L

  h_stat <- ifelse(has_pos, max(pos_plot$value), 0)
  l_stat <- ifelse(has_neg, min(neg_plot$value), 0)
  if (isTRUE(symmetric_legend)) {
    biggest <- max(abs(c(h_stat, l_stat)))
    h_stat <- biggest
    l_stat <- -1*biggest
  }

  # legend color scale limits
  pos_limits <- round(c(pos_thresh, h_stat), stat_decimals)
  neg_limits <- round(c(l_stat, neg_thresh), stat_decimals)

  # remove empty space across images, if requested

  comb_df <- bind_rows(anat_slices, pos_plot, neg_plot) 
  nz_df <- comb_df %>% filter(abs(value) > 1e-5) # find voxels on any layer that are not zero
  empty_d1 <- setdiff(comb_df$dim1, nz_df$dim1)
  empty_d2 <- setdiff(comb_df$dim2, nz_df$dim2)
  
  rm(comb_df)
  if (length(empty_d1) > 0L && isTRUE(remove_null_space)) {
    anat_slices <- anat_slices %>% filter(!dim1 %in% !!empty_d1)
    pos_plot <- pos_plot %>% filter(!dim1 %in% !!empty_d1)
    neg_plot <- neg_plot %>% filter(!dim1 %in% !!empty_d1)
  }
  
  if (length(empty_d2) > 0L && isTRUE(remove_null_space)) {
    anat_slices <- anat_slices %>% filter(!dim2 %in% !!empty_d2)
    pos_plot <- pos_plot %>% filter(!dim2 %in% !!empty_d2)
    neg_plot <- neg_plot %>% filter(!dim2 %in% !!empty_d2)
  }

  # go to split approach, rather than dealing with facet_wrap issues -- need drop = FALSE to keep lengths the same if some slices blank
  anat_slices <- anat_slices %>% group_by(slice, .drop = FALSE) %>% group_split()
  pos_plot <- pos_plot %>% group_by(slice, .drop = FALSE) %>% group_split()
  neg_plot <- neg_plot %>% group_by(slice, .drop = FALSE) %>% group_split()

  # from here: https://joshuacook.netlify.app/post/integer-values-ggplot-axis/
  integer_breaks <- function(n = 5, ...) {
    fxn <- function(x) {
      breaks <- floor(pretty(x, n, ...))
      names(breaks) <- attr(breaks, "labels")
      breaks
    }
    return(fxn)
  }

  # breaks function for including min + max with labels, and a few unlabeled ticks in between
  range_breaks <- function(n=3, ...) {
    fxn <- function(x) {
      breaks <- round(seq(from = min(x, na.rm = T), to = max(x, na.rm = T), length.out = n + 2), ...)
      # breaks <- signif(c(min(x, na.rm = T), max(x, na.rm = T)), 2)

      # names(breaks) <- attr(breaks, "labels")
      bnames <- as.character(breaks)
      bnames[2:(length(bnames) - 1)] <- "" # don't label interior breaks
      names(breaks) <- bnames
      breaks
      print(breaks)
    }
    return(fxn)
  }
  
  make_slice_plot <- function(i) {
    a_df <- anat_slices[[i]]
    p_df <- pos_plot[[i]]
    n_df <- neg_plot[[i]]

    # combine voxels to plot and look for empty/zero rows/columns to trim out
    # comb_df <- bind_rows(anat_slices[[i]], pos_plot[[i]], neg_plot[[i]]) 
    # nz_df <- comb_df %>% filter(abs(value) > 1e-5)
    # empty_d1 <- setdiff(comb_df$dim1, nz_df$dim1)
    # empty_d2 <- setdiff(comb_df$dim2, nz_df$dim2)
    
    # rename columns to ensure that the multiple fill scales in the plot have unique variable names
    a_df <- a_df %>% dplyr::rename(uval = value)
    p_df <- p_df %>% dplyr::rename(pos = value)
    n_df <- n_df %>% dplyr::rename(neg = value)
    
    # if (length(empty_d1) > 0L && isTRUE(remove_null_space)) {
    #   a_df <- a_df %>% filter(!dim1 %in% !!empty_d1)
    #   p_df <- p_df %>% filter(!dim1 %in% !!empty_d1)
    #   n_df <- n_df %>% filter(!dim1 %in% !!empty_d1)
    # }
    
    # if (length(empty_d2) > 0L && isTRUE(remove_null_space)) {
    #   a_df <- a_df %>% filter(!dim2 %in% !!empty_d2)
    #   p_df <- p_df %>% filter(!dim2 %in% !!empty_d2)
    #   n_df <- n_df %>% filter(!dim2 %in% !!empty_d2)
    # }
    
    g <- ggplot(mapping = aes(x=dim1, y=dim2)) +
      #geom_tile(data = anat_slices, mapping = aes(fill=uval), show.legend = FALSE, color = NA, height = 1.01, width = 1.01) +
      #geom_tile(data = anat_slices, mapping = aes(fill=uval, color=uval), show.legend = FALSE) +
      geom_raster(data = a_df, mapping = aes(fill=uval), show.legend = FALSE, interpolate = FALSE) +
      underlay_colorscale # scales[[1]]
    
    # negative overlay values
    if (has_neg) {
      g <- g +
        new_scale_fill() +
        geom_raster(data = n_df, mapping = aes(fill=neg), interpolate = FALSE) +
        negative_colorscale # scales [[2]]
    
      # examples
      # scale_fill_viridis_c(begin = .48, end=0) +
      # scale_fill_distiller("", palette="Blues", direction = 1,  breaks = integer_breaks(), limits = neg_limits,
      #                      guide = guide_colorbar(order = 2))
    }
    
    # positive overlay values
    if (has_pos) {
      g <- g +
        new_scale_fill() + 
        geom_raster(data = p_df, mapping = aes(fill=pos), interpolate = FALSE) +
        positive_colorscale # scales[[3]]

      #scale_fill_viridis_c(begin = 0.52, end=1) +
      # scale_fill_distiller(legend_label, palette="Reds", breaks = integer_breaks(), limits = pos_limits, 
      #                       guide = guide_colorbar(order = 1))  
    }
    
    # modify scale breaks, limits, and guide order -- need to hack this from the existing scale since adding a new scale
    # overrides all of the information, rather than keeping what we want
    if (has_neg) {
      #g$scales$scales[[2]]$breaks <- integer_breaks()
      #g$scales$scales[[2]]$breaks <- pretty_breaks()
      g$scales$scales[[2]]$breaks <- range_breaks(digits = stat_decimals)
      g$scales$scales[[2]]$limits <- neg_limits
      g$scales$scales[[2]]$name <- ifelse(has_pos, "", legend_label) # only add label to neg scale if there are no positive values
      g$scales$scales[[2]]$guide <- guide_colorbar(order = 2, available_aes = c("fill", "fill_new"), ticks.colour = text_color)
    }

    # position of positive and negative color scales in g$scales$scales list
    ppos <- ifelse(has_pos && has_neg, 3, 2)

    if (has_pos) {
      #g$scales$scales[[ppos]]$breaks <- integer_breaks()
      #g$scales$scales[[ppos]]$breaks <- pretty_breaks() # ggplot default for now
      g$scales$scales[[ppos]]$breaks <- range_breaks(digits = stat_decimals) # ggplot default for now
      g$scales$scales[[ppos]]$limits <- pos_limits
      g$scales$scales[[ppos]]$name <- legend_label # label above positive extent
      g$scales$scales[[ppos]]$guide <- guide_colorbar(order = 1, available_aes = c("fill", "fill_new"), ticks.colour = text_color)
    }

    # annotation
    # g <- g +
    #   annotate(geom="text", x = Inf, y = -Inf, label=coords_df$coord_label[i], color = text_color, hjust = 1.5, vjust = -1.5)

    # theme refinements
    g <- g +
      # geom_text(data = coords_df, mapping = aes(label = coord_label, x = Inf, y = -Inf), 
      #           color="white", vjust = 0.5, hjust = 1, size=1) +
      theme_void(base_size = base_size) + coord_fixed() + #+ facet_wrap(~slice, drop=T) + 
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        plot.background = element_rect(fill=background_color, color=ifelse(isTRUE(panel_borders), text_color, NA)),
        panel.background = element_rect(fill=background_color, color=ifelse(isTRUE(panel_borders), text_color, NA)),
        text = element_text(color = text_color),
        legend.spacing.y = unit(0.1, "lines"),
        legend.position = "right",
        plot.margin = unit(c(0.0, 0.5, 0.0, 0.5), "lines") # space on L and R, but not T and B
      )
    
    # add subcaption with cowplot -- note that this makes the object a gtable, not a ggplot object
    # g <- g %>% add_sub(coords_df$coord_label[i], x = 0.9, hjust = 1, color = text_color) # right justified
    
    # new approach in ggplot 3+: use caption for label since that maintains this as a ggplot object
    if (isTRUE(slices$coord_label[i])) {
      g <- g + theme(plot.caption = element_text(hjust = 0.8, vjust = 8)) + labs(caption = coords_df$coord_label[i])  
    }
    
    # add x axis label if requested
    if (!is.null(slices$xlab[i]))  {
      g <- g + theme(axis.title.x = element_text()) + xlab(slices$xlab[i])
    }
    
    # add y axis label if requested
    if (!is.null(slices$ylab[i]))  {
      g <- g + theme(axis.title.y = element_text()) + ylab(slices$ylab[i])
    }
    
    # add title if requested
    if (!is.null(slices$panel_title[i]))  {
      g <- g + ggtitle(slices$panel_title[i])
    }

    # add custom theme elements to each panel, if requested
    if (!is.null(theme_custom)) {
      g <- g + theme_custom
    }

    # cache legend
    leg <- get_legend(g)

    # remove legend from subplots
    g <- g + theme(legend.position = "none")

    # cache legend for extraction
    attr(g, "legend") <- leg

    # ggp <- ggbrain_panel_r6$new(g)

    return(g)
  }
  
  # geom_raster is grumpy about blank space
  #glist <- suppressWarnings(lapply(seq_along(anat_slices), make_slice_plot))
  glist <- lapply(seq_along(anat_slices), make_slice_plot)
  #leg <- get_legend(glist[[1]] + theme(legend.position="right"))
  
  g_all <- plot_grid(plotlist = glist, labels = panel_labels) #, labels = "AUTO", label_colour = text_color)
  g_all <- plot_grid(g_all, attr(glist[[1]], "legend"), rel_widths=c(1, 0.2)) +
    theme(plot.background = element_rect(fill=background_color, color = NA))
  #+theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
  plot(g_all)
}

#' internal function to lookup which slices to display along each axis based on their quantile,
#'   xyz coordinate, or ijk coordinate
#' @param slices A character vector of coordinates for slices to display
#' @param underlay the RNifti object of the underlay that defines the coordinate system
#' @param remove_null_space If TRUE, plots are trimmed for empty space.
#' @keywords internal
lookup_slices <- function(slices, underlay, remove_null_space = TRUE) {
  checkmate::assert_character(slices)
  
  if (isTRUE(remove_null_space)) {
    nzpos <- which(abs(underlay) > 1e-6, arr.ind=TRUE)
    minx <- min(nzpos[,1])
    maxx <- max(nzpos[,1])
    miny <- min(nzpos[,2])
    maxy <- max(nzpos[,2])
    minz <- min(nzpos[,3])
    maxz <- max(nzpos[,3])
    
    xrange <- min(nzpos[,1]):max(nzpos[,1])
    yrange <- min(nzpos[,2]):max(nzpos[,2])
    zrange <- min(nzpos[,3]):max(nzpos[,3])
    
  } else {
    xrange <- seq_len(dim(underlay)[1])
    yrange <- seq_len(dim(underlay)[2])
    zrange <- seq_len(dim(underlay)[3])
  }
  
  # translate ijk to xyz for each axis
  xcoords <- RNifti::voxelToWorld(cbind(seq_len(dim(underlay)[1]), 1, 1), underlay)[,1]
  ycoords <- RNifti::voxelToWorld(cbind(1, seq_len(dim(underlay)[2]), 1), underlay)[,2]
  zcoords <- RNifti::voxelToWorld(cbind(1, 1, seq_len(dim(underlay)[3])), underlay)[,3]
  
  # helper subfunction to lookup slice number, plane, and label for any ijk, xyz, or % input
  get_slice_num <- function(coord_str) {
    res <- tolower(trimws(strsplit(coord_str, "\\s*=\\s*", perl = TRUE)[[1]]))
    axis <- res[1]
    number <- res[2]
    is_pct <- grepl("[\\d.]+%", number, perl = TRUE)
    if (isTRUE(is_pct)) {
      number <- as.numeric(sub("%", "", number, fixed=TRUE))
      checkmate::assert_number(number, lower = 0, upper = 100)
      number <- number/100 # convert to quantile
    } else {
      number <- as.numeric(number)
    }
      
    # determine plane of slice to display
    plane <- switch(
      axis,
      i = "sagittal",
      j = "coronal",
      k = "axial",
      x = "sagittal",
      y = "coronal",
      z = "axial",
      stop(sprintf("Cannot interpret input: %s", coord_str))
    )
    
    # determine world or voxel coordinate system
    coord <- switch(
      axis,
      i = "voxel",
      j = "voxel",
      k = "voxel",
      x = "world",
      y = "world",
      z = "world",
      stop(sprintf("Cannot interpret input: %s", coord_str))
    )
    
    axis_label <- switch(
      plane,
      sagittal = "x",
      coronal = "y",
      axial = "z"
    )
    
    # validate input and lookup slice
    if (isTRUE(is_pct)) {
      if (plane == "sagittal") {
        rr <- xrange
      } else if (plane == "coronal") {
        rr <- yrange
      } else if (plane == "axial") {
        rr <- zrange
      }
      
      slc_num <- round(quantile(rr, number))
    } else {
      if (coord == "world") { # xyz
        if (plane == "sagittal") {
          coords <- xcoords
        } else if (plane == "coronal") {
          coords <- ycoords
        } else if (plane == "axial") {
          coords <- zcoords
        }
        
        checkmate::assert_number(number, lower=min(coords), upper=max(coords))
        slc_num <- which.min(abs(number - coords))
      } else if (coord == "voxel") { #ijk
        if (plane == "sagittal") {
          coords <- seq_len(dim(underlay)[1])
        } else if (plane == "coronal") {
          coords <- seq_len(dim(underlay)[2])
        } else if (plane == "axial") {
          coords <- seq_len(dim(underlay)[3])
        }
        
        checkmate::assert_integerish(number, lower=min(coords), upper=max(coords), len=1L)
        slc_num <- number
      }
    }
    
    # slc_num is the slice number in the plane of interest
    if (axis_label == "x") {
      slc_coords <- xcoords[slc_num]
    } else if (axis_label == "y") {
      slc_coords <- ycoords[slc_num]
    } else if (axis_label == "z") {
      slc_coords <- zcoords[slc_num]
    }
    
    slc_coords <- round(slc_coords, 1) # for display
    df <- data.frame(plane = plane, slice_number = slc_num, coord_label = paste(axis_label, "=", slc_coords))
    
    # we want to respect the user arrangement order
    # but for extracti
    
    return(df)
  
  }
  
  slice_df <- lapply(slices, get_slice_num) %>% 
    bind_rows() %>%
    distinct() %>% # remove any dupes
    tibble::remove_rownames() %>% # unneeded labels
    mutate(slice_index = 1:n())
  
  return(slice_df)
}




## layer(
#   geom = "raster",
#   data = anat_slices,
#   stat = "identity",
#   position = "identity",
#   
# )
# 
# pal_extreme <- max(abs(overlay_slices$oval))
# pal <- c(viridis_pal(begin=0, end=0.5)(2), "grey90", "grey90", viridis_pal(begin=0.5, end=1)(2))
# val <- scales::rescale(c(-1*pal_extreme, neg_thresh+.01, neg_thresh+.011, pos_thresh-.011, pos_thresh-.01, pal_extreme))
# 
# slice_plots <- ggplot(mapping = aes(x=dim1, y=dim2)) +
#   geom_tile(data = anat_slices, mapping = aes(fill=uval), show.legend = FALSE) +
#   scale_fill_gradient(low="grey10", high="grey90") +
#   #guides(fill="none") +
#   new_scale_fill() + 
#   geom_tile(data = overlay_slices %>% filter(oval >= !!pos_thresh | oval <= !!neg_thresh), mapping = aes(fill=oval)) + 
#   #scale_fill_viridis_c() +
#   #scale_fill_gradientn(colors = c(viridis_pal(begin=0, end=0.5)(2), "grey90", viridis_pal(begin=0.5, end=1)(2)), 
#   #                                 values=scales::rescale(c(-1, -0.2, 0, 0.2, 1))) +
#                        
#   scale_fill_gradientn(colors = pal, values=val) +
#   
#   coord_fixed() + theme_void() + facet_wrap(~slice)



# slice_plots <- ggplot(mapping = aes(x=dim1, y=dim2)) +
#   #geom_tile(data = anat_slices, mapping = aes(fill=uval), show.legend = FALSE, color = NA, height = 1.01, width = 1.01) +
#   #geom_tile(data = anat_slices, mapping = aes(fill=uval, color=uval), show.legend = FALSE) +
#   geom_raster(data = anat_slices, mapping = aes(fill=uval), show.legend = FALSE, interpolate = FALSE) +
#   scale_fill_gradient(low=background_color, high="grey90") +
#   scale_color_gradient(low=background_color, high="grey90") +
#   #guides(fill="none") +
#   
#   # negative overlay values
#   new_scale_fill() + 
#   new_scale_color() + 
#   #geom_tile(data = neg_plot, mapping = aes(fill=neg, color = neg)) + 
#   geom_raster(data = neg_plot, mapping = aes(fill=neg), interpolate = FALSE) +
#   #scale_fill_viridis_c(begin = .48, end=0) +
#   scale_fill_distiller("", palette="Blues", direction = 1, guide = guide_colorbar(order = 2)) + 
#   scale_color_distiller("", palette="Blues", direction = 1, guide = guide_colorbar(order = 2)) + 
#   
#   # positive overlay values
#   new_scale_fill() + 
#   new_scale_color() + 
#   #geom_tile(data = pos_plot, mapping = aes(fill=pos, color = pos)) + 
#   geom_raster(data = pos_plot, mapping = aes(fill=pos), interpolate = FALSE) +
#   #scale_fill_viridis_c(begin = 0.52, end=1) +
#   scale_fill_distiller(legend_label, palette="Reds", guide = guide_colorbar(order = 1)) +
#   scale_color_distiller(legend_label, palette="Reds", guide = guide_colorbar(order = 1)) +
#   geom_text(data = coords_df, mapping = aes(label = coord_label, x = Inf, y = -Inf), 
#             color="white", vjust = 0.5, hjust = 1, size=1) +
#   theme_void() + facet_wrap(~slice, drop=T) + coord_fixed() + 
#   theme(panel.spacing = unit(-1, "lines"), strip.background = element_blank(),
#         strip.text.x = element_blank(), plot.background = element_rect(fill=background_color, color=NA),
#         text = element_text(color = text_color)
#   )
# 
# # handle the need to extend the background color beyond the coord_fixed borders using cowplot
# # https://stackoverflow.com/questions/23349181/how-to-color-entire-background-in-ggplot2-when-using-coord-fixed
# cowplot::ggdraw(slice_plots) + 
#   theme(plot.background = element_rect(fill=background_color, color = NA))
# 
#plot(slice_plots)
