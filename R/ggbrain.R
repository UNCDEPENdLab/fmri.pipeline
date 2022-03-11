#' Plot fMRI data on an underlay image
#' 
#' @param underlay a 3D nifti image used for the image underlay (default b/w)
#' @param overlay a 4D nifti image used for plotting stats on underlay (color)
#' @param color_col a position in the 4th dim of overlay use to color plots
#' @param alpha_col a position in the 4th dim of overlay use to set alpha transparency of plots
#' @param underlay_colorscale A ggplot scale_fill_* function call used for coloration of underlay
#' @param overlay_colorscale A ggplot scale_fill_* function call used for coloration of overlay
#' @param coronal_slices A list of coronal slices to be displayed. Elements can be \code{xyz}, \code{ijk},
#'   and/or \code{quantiles}.
#' @param sagittal_slices A list of sagittal slices to be displayed. Elements can be \code{xyz}, \code{ijk},
#'   and/or \code{quantiles}.
#' @param axial_slices A list of axial slices to be displayed. Elements can be \code{xyz}, \code{ijk},
#'   and/or \code{quantiles}
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
#' @importFrom ggplot geom_raster scale_fill_gradient theme
#' @export
ggbrain <- function(underlay=NULL, overlay=NULL, 
                    color_col=NULL, alpha_col=NULL,
                    underlay_colorscale=scale_fill_gradient(low="grey10", high="grey90"),
                    overlay_colorscale=scale_fill_gradient2(midpoint = 0, low = "blue", mid="grey90", high="red"),
                    coronal_slices=list(quantiles = c(.25, .50, .75), ijk = NULL, xyz = NULL),
                    sagittal_slices=list(quantiles = c(.25, .50, .75), ijk = NULL, xyz = NULL),
                    axial_slices=list(quantiles = c(.25, .50, .75), ijk = NULL, xyz = NULL),
                    remove_null_space=TRUE,
                    pos_thresh = 3,
                    neg_thresh = -3, legend_label = "z",
                    background_color = "gray10", text_color = "white",
                    trim_underlay = c(.01, .99), zero_underlay = 1e-3, symmetric_legend = TRUE,
                    add_panel_labels = FALSE
) {
  
  max_bg <- "gray95" # brightest value on underlay
  
  require(ggnewscale)
  require(cowplot)
  checkmate::assert_file_exists(underlay)
  checkmate::assert_file_exists(overlay)
  
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
  
  #verify that i,j,k (1,2,3) dimensions of underlay match dimensions of overlay
  stopifnot(identical(dim(underlay)[1:3], dim(overlay)[1:3]))
  
  coords <- lookup_slices(sagittal_slices, coronal_slices, axial_slices, underlay, remove_null_space)
  coords_df <- bind_rows(coords) %>% group_by(slice) %>% group_split()
  
  get_slices <- function(img) {
    slc <- list()
    arg_name <- as.character(match.call())[-1]
    if (!is.null(sagittal_slices)) {
      slc[["sagittal"]] <- aperm(img[coords$sagittal$slice_number,,,drop=FALSE], c(1,2,3))
    }
    
    if (!is.null(coronal_slices)) {
      slc[["coronal"]] <- aperm(img[,coords$coronal$slice_number,,drop=FALSE], c(2,1,3))
    }
    
    if (!is.null(axial_slices)) {
      slc[["axial"]] <- aperm(img[,,coords$axial$slice_number,drop=FALSE], c(3,1,2))
    }

    slc <- dplyr::bind_rows(lapply(seq_along(slc), function(ll) {
      df <- reshape2::melt(slc[[ll]], varnames=c("slice", "dim1", "dim2"))
      df$type <- names(slc)[ll]
      
      #uniquely identify slices by number and type
      df$slice <- paste0(substr(names(slc)[ll],1,1), df$slice)
      return(df)
    }))
    
    # slc <- lapply(seq_along(slc), function(ll) {
    #   df <- reshape2::melt(slc[[ll]], varnames=c("slice", "dim1", "dim2"))
    #   df$type <- names(slc)[ll]
    #   
    #   #uniquely identify slices by number and type
    #   df$slice <- paste0(substr(names(slc)[ll],1,1), df$slice)
    #   df$image <- arg_name
    #   return(df)
    # }))
    
    slc <- slc %>% mutate(
      image = arg_name,
      slice = factor(slice) # store as factor so that if a slice drops in the filter, we still get an empty list
    )
    return(slc)
  }
  
  #validate 0-1 bounds on arguments: axial_slices, sagittal_slices, coronal_slices
  anat_slices <- get_slices(underlay) #%>% rename(uval = value) # for fill scale to work properly with overlay, need different variable names
  overlay_slices <- get_slices(overlay) #%>% rename(oval = value)
  
  # new solution: separate layers for pos and neg
  pos_plot <- overlay_slices %>% filter(value >= !!pos_thresh) #%>% rename(pos = oval)
  neg_plot <- overlay_slices %>% filter(value <= !!neg_thresh) #%>% rename(neg = oval)
  
  h_stat <- max(pos_plot$value)
  l_stat <- min(neg_plot$value)
  if (isTRUE(symmetric_legend)) {
    biggest <- max(abs(c(h_stat, l_stat)))
    h_stat <- biggest
    l_stat <- -1*biggest
  }
  
  # legend color scale limits
  pos_breaks <- c(pos_thresh, h_stat)
  neg_breaks <- c(l_stat, neg_thresh)
  
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
  
  make_slice_plot <- function(i) {
    a_df <- anat_slices[[i]]
    p_df <- pos_plot[[i]]
    n_df <- neg_plot[[i]]
    
    # combine voxels to plot and look for empty/zero rows/columns to trim out
    comb_df <- bind_rows(anat_slices[[i]], pos_plot[[i]], neg_plot[[i]]) 
    nz_df <- comb_df %>% filter(abs(value) > 1e-5)
    empty_d1 <- setdiff(comb_df$dim1, nz_df$dim1)
    empty_d2 <- setdiff(comb_df$dim2, nz_df$dim2)
    
    # rename columns to ensure that the multiple fill scales in the plot have unique variable names
    a_df <- a_df %>% dplyr::rename(uval = value)
    p_df <- p_df %>% dplyr::rename(pos = value)
    n_df <- n_df %>% dplyr::rename(neg = value)
    
    if (length(empty_d1) > 0L && isTRUE(remove_null_space)) {
      a_df <- a_df %>% filter(!dim1 %in% !!empty_d1)
      p_df <- p_df %>% filter(!dim1 %in% !!empty_d1)
      n_df <- n_df %>% filter(!dim1 %in% !!empty_d1)
    }
    
    if (length(empty_d2) > 0L && isTRUE(remove_null_space)) {
      a_df <- a_df %>% filter(!dim2 %in% !!empty_d2)
      p_df <- p_df %>% filter(!dim2 %in% !!empty_d2)
      n_df <- n_df %>% filter(!dim2 %in% !!empty_d2)
    }
    
    g <- ggplot(mapping = aes(x=dim1, y=dim2)) +
      #geom_tile(data = anat_slices, mapping = aes(fill=uval), show.legend = FALSE, color = NA, height = 1.01, width = 1.01) +
      #geom_tile(data = anat_slices, mapping = aes(fill=uval, color=uval), show.legend = FALSE) +
      geom_raster(data = a_df, mapping = aes(fill=uval), show.legend = FALSE, interpolate = FALSE) +
      scale_fill_gradient(low=background_color, high=max_bg)
    
    # negative overlay values
    g <- g + 
      new_scale_fill() + 
      #new_scale_color() + 
      #geom_tile(data = neg_plot, mapping = aes(fill=neg, color = neg)) + 
      geom_raster(data = n_df, mapping = aes(fill=neg), interpolate = FALSE) +
      #scale_fill_viridis_c(begin = .48, end=0) +
      scale_fill_distiller("", palette="Blues", direction = 1,  breaks = integer_breaks(), limits = neg_breaks,
                           guide = guide_colorbar(order = 2))
    
    # positive overlay values
    g <- g +
      new_scale_fill() + 
      #new_scale_color() + 
      #geom_tile(data = pos_plot, mapping = aes(fill=pos, color = pos)) + 
      geom_raster(data = p_df, mapping = aes(fill=pos), interpolate = FALSE) +
      #scale_fill_viridis_c(begin = 0.52, end=1) +
      scale_fill_distiller(legend_label, palette="Reds", breaks = integer_breaks(), limits = pos_breaks, 
                           guide = guide_colorbar(order = 1))
    
    # annotation
    # g <- g +
    #   annotate(geom="text", x = Inf, y = -Inf, label=coords_df[[i]]$coord_label, color = text_color, hjust = 1.5, vjust = -1.5)
    
    # theme refinements
    g <- g +
      # geom_text(data = coords_df, mapping = aes(label = coord_label, x = Inf, y = -Inf), 
      #           color="white", vjust = 0.5, hjust = 1, size=1) +
      theme_void() + coord_fixed() + #+ facet_wrap(~slice, drop=T) + 
      theme(
        panel.spacing = unit(-1, "lines"), 
        strip.background = element_blank(),
        strip.text.x = element_blank(), 
        plot.background = element_rect(fill=background_color, color=NA),
        text = element_text(color = text_color),
        legend.spacing.y = unit(0.1, "lines"),
        legend.position = "right",
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")
      )
    
    # cache legend
    leg <- get_legend(g)
    
    # remove legend from subplots
    g <- g + theme(legend.position = "none")
    
    # add subcaption with cowplot -- note that this makes the object a gtable, not a ggplot object
    g <- g %>% add_sub(coords_df[[i]]$coord_label, x = 0.9, hjust = 1, color = text_color) # right justified
    
    # cache legend for extraction
    attr(g, "legend") <- leg
    
    return(g)
  }
  
  # geom_raster is grumpy about blank space
  glist <- suppressWarnings(lapply(seq_along(anat_slices), make_slice_plot))
  #leg <- get_legend(glist[[1]] + theme(legend.position="right"))
  
  g_all <- plot_grid(plotlist = glist) #, labels = "AUTO", label_colour = text_color)
  g_all <- plot_grid(g_all, attr(glist[[1]], "legend"), rel_widths=c(1, 0.1)) + 
    theme(plot.background = element_rect(fill=background_color, color = NA))
  
  plot(g_all)
}

#' internal function to lookup which slices to display along each axis based on their quantile,
#'   xyz coordinate, or ijk coordinate
#' @param sagittal_slices List of sagittal slices to display
#' @param coronal_slices List of coronal slices to display
#' @param axial_slices List of axial slices to display
#' @param underlay the RNifti object of the underlay that defines the coordinate system
#' @keywords internal
lookup_slices <- function(sagittal_slices = NULL, coronal_slices = NULL, axial_slices = NULL, underlay, remove_null_space = TRUE) {
  checkmate::assert_list(sagittal_slices)
  checkmate::assert_list(coronal_slices)
  checkmate::assert_list(axial_slices)
  
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
  
  # look up coordinates and slice numbers for requested slice positions
  coords <- list()
  sag_nums <- cor_nums <- ax_nums <- NULL # slice numbers
  
  if (!is.null(sagittal_slices)) {
    sag_nums <- c()
    if (!is.null(sagittal_slices$quantiles)) {
      sag_nums <- round(quantile(xrange, sagittal_slices$quantiles))
    }
    
    if (!is.null(sagittal_slices$ijk)) {
      checkmate::assert_integerish(sagittal_slices$ijk, lower = 1, upper = dim(underlay)[1])
      sag_nums <- c(sag_nums, sagittal_slices$ijk)
    }
    
    if (!is.null(sagittal_slices$xyz)) {
      checkmate::assert_integerish(sagittal_slices$xyz, lower = min(xcoords), upper = max(xcoords))
      sag_nums <- c(sag_nums, sapply(sagittal_slices$xyz, function(p) which.min(abs(p - xcoords))))
    }
    
    sag_nums <- sort(unique(sag_nums)) # don't display the same slice twice
    #sag_coords <- round(RNifti::voxelToWorld(cbind(sag_nums, 1, 1), underlay)[,1], 1)
    sag_coords <- round(xcoords[sag_nums], 1)
    coords[["sagittal"]] <- data.frame(slice = paste0("s", seq_along(sag_nums)), slice_number = sag_nums, coord_label = paste("x =", sag_coords))
  }
  
  if (!is.null(coronal_slices)) {
    cor_nums <- c()
    if (!is.null(coronal_slices$quantiles)) {
      cor_nums <- unique(round(quantile(yrange, coronal_slices$quantiles)))
    }
    
    if (!is.null(coronal_slices$ijk)) {
      checkmate::assert_integerish(coronal_slices$ijk, lower = 1, upper = dim(underlay)[2])
      cor_nums <- c(cor_nums, coronal_slices$ijk)
    }
    
    if (!is.null(coronal_slices$xyz)) {
      checkmate::assert_integerish(coronal_slices$xyz, lower = min(ycoords), upper = max(ycoords))
      cor_nums <- c(cor_nums, sapply(coronal_slices$xyz, function(p) which.min(abs(p - ycoords))))
    }
    
    cor_nums <- sort(unique(cor_nums))
    cor_coords <- round(ycoords[cor_nums], 1)
    coords[["coronal"]] <- data.frame(slice = paste0("c", seq_along(cor_nums)), slice_number = cor_nums, coord_label = paste("y =", cor_coords))
  }
  
  if (!is.null(axial_slices)) {
    ax_nums <- c()
    if (!is.null(axial_slices$quantiles)) {
      checkmate::assert_numeric(axial_slices$quantiles, lower=0, upper=1)
      ax_nums <- c(ax_nums, round(quantile(zrange, axial_slices$quantiles)))
    }
    
    if (!is.null(axial_slices$ijk)) {
      checkmate::assert_integerish(axial_slices$ijk, lower = 1, upper = dim(underlay)[3])
      ax_nums <- c(ax_nums, axial_slices$ijk)
    }
    
    if (!is.null(axial_slices$xyz)) {
      checkmate::assert_integerish(axial_slices$xyz, lower = min(zcoords), upper = max(zcoords))
      ax_nums <- c(ax_nums, sapply(axial_slices$xyz, function(p) which.min(abs(p - zcoords))))
      
      # skipping this for now because it returns matrix for multiple slices, and vector for single slice (indexing code)
      #ax_nums <- c(ax_nums, RNifti::worldToVoxel(cbind(0, rep(axial_slices$xyz, 2), 0), underlay)[,3])
    }
    
    ax_nums <- sort(unique(ax_nums))
    ax_coords <- round(zcoords[ax_nums], 1)
    coords[["axial"]] <- data.frame(slice = paste0("a", seq_along(ax_nums)), slice_number = ax_nums, coord_label = paste("z =", ax_coords))
  }
  
  return(coords)
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
