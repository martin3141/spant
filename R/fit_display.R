#' Plot the fitting results of an object of class \code{fit_result}.
#' @param x fit_result object.
#' @param dyn the dynamic index to plot.
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param coil the coil element number to plot.
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param data_only display only the processed data (logical).
#' @param label character string to add to the top left of the plot window.
#' @param plot_sigs a character vector of signal names to add to the plot.
#' @param n single index element to plot (overrides other indices when given).
#' @param sub_bl subtract the baseline from the data and fit (logical).
#' @param mar option to adjust the plot margins. See ?par.
#' @param restore_def_par restore default plotting par values after the plot has 
#' been made.
#' @param ylim range of values to display on the y-axis, eg ylim = c(0,10).
#' @param y_scale option to display the y-axis values (logical).
#' @param show_grid plot gridlines behind the data (logical). Defaults to TRUE.
#' @param grid_nx number of cells of the grid in x and y direction. When NULL
#' the grid aligns with the tick marks on the corresponding default axis (i.e.,
#' tickmarks as computed by axTicks). When NA, no grid lines are drawn in the
#' corresponding direction.
#' @param grid_ny as above.
#' @param ... further arguments to plot method.
#' @export
plot.fit_result <- function(x, dyn = 1, x_pos = 1, y_pos = 1, z_pos = 1,
                            coil = 1, xlim = NULL, data_only = FALSE,
                            label = NULL, plot_sigs = NULL, n = NULL,
                            sub_bl = FALSE, mar = NULL, restore_def_par = TRUE, 
                            ylim = NULL, y_scale = FALSE, show_grid = TRUE,
                            grid_nx = NULL, grid_ny = NA, ...) {
  
  .pardefault <- graphics::par(no.readonly = T)
  
  if (is.null(n)) {
    ind <- (x$res_tab$X == x_pos) & (x$res_tab$Y == y_pos) & 
           (x$res_tab$Z == z_pos) & (x$res_tab$Dynamic == dyn) &
           (x$res_tab$Coil == coil) 
    
    n <- which(ind)
  }
  
  opts <- x$opts
  x <- x$fits[[n]]
  
  if (anyNA(x)) { 
    graphics::plot.new()
    return(NULL)
  }
  
  if (is.null(xlim)) {
    if (class(opts) == "list") {
      if ((!is.null(opts$ppm_left)) & (!is.null(opts$ppm_right))) {
        xlim <- c(opts$ppm_left, opts$ppm_right)
      } else {
        xlim <- rev(range(x$PPMScale))
      }
    } else {
      xlim <- rev(range(x$PPMScale))
    }
  }
  
  graphics::par("xaxs" = "i", "yaxs" = "i") # tight axes limits
  
  if (!is.null(mar)) {
    graphics::par(mar = mar)
  } else {
    if (y_scale) {
      graphics::par(mar = c(3.5, 3.5, 1.2, 1.2)) # space around the plot
    } else {
      graphics::par(mar = c(3.5, 1.2, 1.2, 1.2)) # space around the plot
    }
  }
  
  if (y_scale) {
    ylab <- "Intensity (au)"
    yaxt <- "s"
  } else {
    ylab <- ""
    yaxt <- "n"
  }
  
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  
  if (xlim[1] > xlim[2]) {
    ind <- x$PPMScale < xlim[1] & x$PPMScale > xlim[2]
  }
  else {
    ind <- x$PPMScale < xlim[2] & x$PPMScale > xlim[1]
  }
  
  if (data_only) {
    max_dp <- max(x$Data[ind])
    min_dp <- min(x$Data[ind])
    marg = (max_dp - min_dp) * 0.02
    
    graphics::plot(x$PPMScale, x$Data, type = 'l', xlim = xlim, 
         ylim = c(min_dp - marg, max_dp + marg), yaxt = yaxt, ylab = ylab,
         xlab = "Chemical shift (ppm)",
         panel.first = {if (show_grid) graphics::grid(nx = grid_nx,
                                                      ny = grid_ny)}, ...)
    
    if (!is.null(label)) {
      graphics::par(xpd = TRUE)
      graphics::text(xlim[1], max_dp, label, cex = 2.5)
      graphics::par(xpd = FALSE) 
    }
    
  } else {
    if (sub_bl) {
      x$Data <- x$Data - x$Baseline
      x$Baseline <- rep(0,length(x$Baseline))
    }
    
    fit_line <- x$Fit + x$Baseline
    
    if (is.null(ylim)) {
      max_dp <- max(x$Data[ind], fit_line[ind])
      min_dp <- min(x$Data[ind], fit_line[ind], x$Baseline[ind])
    } else {
      max_dp <- ylim[2]
      min_dp <- ylim[1]
    }
    
    res <- x$Data - fit_line
    res_range <- max(res[ind]) - min(res[ind])
    offset <- max_dp - min(res[ind])
    
    if (!is.null(ylim)) {
      x$Data[ind][x$Data[ind] > max_dp] <- max_dp
      fit_line[fit_line > max_dp] <- max_dp
    }
    
    graphics::plot(x$PPMScale, x$Data, type = 'l', xlim = xlim, 
         ylim = c(min_dp,max_dp + res_range), yaxt = yaxt, ylab = ylab,
         xlab = "Chemical shift (ppm)",
         panel.first = {if (show_grid) graphics::grid(nx = grid_nx,
                                                      ny = grid_ny)}, ...)
    
    graphics::lines(x$PPMScale, fit_line, col = 'Red', lw = 2)
    if (!sub_bl) {
      graphics::lines(x$PPMScale, x$Baseline)
    }
    graphics::lines(x$PPMScale, res + offset)
    graphics::abline(h = max_dp)
  }
  
  for (sig in plot_sigs) {
    graphics::lines(x$PPMScale, x[sig][[1]] + x$Baseline, col = "blue")
  }
  
  if (restore_def_par) graphics::par(.pardefault)
}

#' @export
summary.fit_result <- function(object, ...) {
  x <- object
  cat("Summary of data quality\n", sep = "")
  cat("-----------------------\n", sep = "")
  spectra <- nrow(stats::na.omit(x$res_tab))
  cat("Spectra analysed : ", spectra, "\n\n", sep = "")
  
  mean_lw <- sprintf("%.4f", mean(x$res_tab$tNAA_lw, na.rm = TRUE))
  sd_lw <-   sprintf("%.4f", sd(x$res_tab$tNAA_lw, na.rm = TRUE))
  max_lw <-  sprintf("%.4f", max(x$res_tab$tNAA_lw, na.rm = TRUE))
  min_lw <-  sprintf("%.4f", min(x$res_tab$tNAA_lw, na.rm = TRUE))
  
  cat("Mean FWHM : ", mean_lw, " ppm\n", sep = "")
  cat("SD   FWHM : ", sd_lw, " ppm\n", sep = "")
  cat("Max  FWHM : ", max_lw, " ppm\n", sep = "")
  cat("Min  FWHM : ", min_lw, " ppm\n\n", sep = "")
  
  mean_snr <- sprintf("%.0f", mean(x$res_tab$SNR, na.rm = TRUE))
  sd_snr <-   sprintf("%.0f", sd(x$res_tab$SNR, na.rm = TRUE))
  max_snr <-  sprintf("%.0f", max(x$res_tab$SNR, na.rm = TRUE))
  min_snr <-  sprintf("%.0f", min(x$res_tab$SNR, na.rm = TRUE))
  
  cat("Mean SNR  : ", mean_snr, "\n", sep = "")
  cat("SD   SNR  : ", sd_snr, "\n", sep = "")
  cat("Max  SNR  : ", max_snr, "\n", sep = "")
  cat("Min  SNR  : ", min_snr, "\n\n", sep = "")
  
  mean_fqn <- sprintf("%.2f", mean(x$res_tab$FQN, na.rm = TRUE))
  sd_fqn <-   sprintf("%.2f", sd(x$res_tab$FQN, na.rm = TRUE))
  max_fqn <-  sprintf("%.2f", max(x$res_tab$FQN, na.rm = TRUE))
  min_fqn <-  sprintf("%.2f", min(x$res_tab$FQN, na.rm = TRUE))
  
  cat("Mean FQN  : ", mean_fqn, "\n", sep = "")
  cat("SD   FQN  : ", sd_fqn, "\n", sep = "")
  cat("Max  FQN  : ", max_fqn, "\n", sep = "")
  cat("Min  FQN  : ", min_fqn, "\n", sep = "")
}

#' Plot the fitting results of an object of class \code{fit_result} with 
#' individual basis set components shown.
#' @param x fit_result object.
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param y_offset separate basis signals in the y-axis direction by this value.
#' @param dyn the dynamic index to plot.
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param coil the coil element number to plot.
#' @param n single index element to plot (overrides other indices when given).
#' @param sub_bl subtract the baseline from the data and fit (logical).
#' @param labels print signal labels at the right side of the plot.
#' @param label_names provide a character vector of signal names to replace the
#' defaults determined from the basis set.
#' @param sig_col colour of individual signal components.
#' @param restore_def_par restore default plotting par values after the plot has 
#' been made.
#' @param omit_signals a character vector of basis signal names to be removed 
#' from the plot.
#' @param combine_lipmm combine all basis signals with names starting with "Lip"
#' or "MM".
#' @param combine_metab combine all basis signals with names not starting with
#' "Lip" or "MM".
#' @param mar option to adjust the plot margins. See ?par.
#' @param show_grid plot gridlines behind the data (logical). Defaults to TRUE.
#' @param grid_nx number of cells of the grid in x and y direction. When NULL
#' the grid aligns with the tick marks on the corresponding default axis (i.e.,
#' tickmarks as computed by axTicks). When NA, no grid lines are drawn in the
#' corresponding direction.
#' @param grid_ny as above.
#' @param ... further arguments to plot method.
#' @export
stackplot.fit_result <- function(x, xlim = NULL, y_offset = 0, dyn = 1, 
                                 x_pos = 1, y_pos = 1, z_pos = 1, coil = 1,
                                 n = NULL, sub_bl = FALSE, labels = FALSE,
                                 label_names = NULL, sig_col = "black",
                                 restore_def_par = TRUE, omit_signals = NULL,
                                 combine_lipmm = FALSE, combine_metab = FALSE,
                                 mar = NULL, show_grid = TRUE, grid_nx = NULL,
                                 grid_ny = NA,...) {
  
  .pardefault <- graphics::par(no.readonly = T)
  
  if (is.null(n)) {
    ind <- (x$res_tab$X == x_pos) & (x$res_tab$Y == y_pos) & 
           (x$res_tab$Z == z_pos) & (x$res_tab$Dynamic == dyn) &
           (x$res_tab$Coil == coil) 
    
    n <- which(ind)
  }
  
  opts <- x$opts
  x <- x$fits[[n]]
  
  # remove any signals that were requested
  if (!is.null(omit_signals)) x <- x[, !(names(x) %in% omit_signals)]
  
  if (combine_lipmm) {
    # find lip/mm indices
    indices <- c(grep("^Lip",colnames(x)), grep("^MM",colnames(x)))
    new_col <- rowSums(x[indices])
    x <- x[, -indices]
    x$LipMM <- new_col
    cols <- length(colnames(x))
    reorder <- c(1:4, cols, 5:(cols-1))
    x <- x[,reorder]
  }
  
  if (combine_metab) {
    # find lip/mm indices
    matches <- grepl("^Lip", colnames(x)) | grepl("^MM", colnames(x))
    matches <- !matches
    indices <- which(matches)[-c(1:4)]
    new_col <- rowSums(x[indices])
    x <- x[, -indices]
    x$Metabs <- new_col
    cols <- length(colnames(x))
    reorder <- c(1:4, cols, 5:(cols-1))
    x <- x[,reorder]
  }
  
  # label names 
  if (is.null(label_names)) {
    names <- colnames(x)[5:length(colnames(x))]
  } else {
    names <- label_names
  }
  
  if (is.null(xlim)) {
    if ((!is.null(opts$ppm_left)) & (!is.null(opts$ppm_right))) {
      xlim <- c(opts$ppm_left, opts$ppm_right)
    } else {
      xlim <- rev(range(x$PPMScale))
    }
  }
  
  graphics::par("xaxs" = "i", "yaxs" = "i") # tight axes limits
  
  right_mar <- ifelse(labels, 4, 1.2)
  
  graphics::par(mar = c(3.5, 1.2, 1.2, right_mar)) # space around the plot
  
  if (!is.null(mar)) graphics::par(mar = mar)
  
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  
  if (xlim[1] > xlim[2]) {
    ind <- x$PPMScale < xlim[1] & x$PPMScale > xlim[2]
  }
  else {
    ind <- x$PPMScale < xlim[2] & x$PPMScale > xlim[1]
  }
  
  if (sub_bl) {
    x$Data <- x$Data - x$Baseline
    x$Baseline <- rep(0,length(x$Baseline))
  }
  
  fit_line <- x$Fit + x$Baseline
  max_dp <- max(x$Data[ind],fit_line[ind])
  min_dp <- min(x$Data[ind],fit_line[ind],x$Baseline[ind])
  
  res <- x$Data - fit_line
  res_range <- max(res[ind]) - min(res[ind])
  offset <- max_dp - min(res[ind])
  
  basis_yoff <- (max_dp - min_dp) * y_offset / length(names)
  
  min_basis <- Inf
  for (p in 5:ncol(x)) {
    x[,p] <- x[,p] - (p - 4) * basis_yoff + min_dp
    if (min(x[ind, p]) < min_basis) min_basis <- min(x[ind, p])
  }
  
  graphics::plot(x$PPMScale, x$Data, type = 'l', xlim = xlim, 
       ylim = c(min_basis - basis_yoff, max_dp + res_range), yaxt = "n",
       ylab = "", xlab = "Chemical shift (ppm)", panel.first = {if (show_grid)
                                                    graphics::grid(nx = grid_nx,
                                                    ny = grid_ny)}, ...)
  
  graphics::lines(x$PPMScale, fit_line, col = 'Red', lw = 2)
  
  if (!sub_bl) {
    graphics::lines(x$PPMScale, x$Baseline)
  }
  
  graphics::lines(x$PPMScale, res + offset)
  graphics::abline(h = max_dp)
  
  for (p in 5:ncol(x)) {
    if (substr(colnames(x)[p], 1, 3) == "SP_") {
      graphics::lines(x$PPMScale, x[,p], col = "blue")
    } else {
      graphics::lines(x$PPMScale, x[,p], col = sig_col)
    }
  }
  
  if (labels) {
    for (p in 5:ncol(x)) {
      graphics::text(graphics::par("usr")[2], -(p - 4) * basis_yoff + min_dp, 
                     names[p - 4], xpd = TRUE, pos = 4, offset = 0.25)
    }
  }
  
  if (restore_def_par) graphics::par(.pardefault)
}

#' Print a summary of an object of class \code{fit_result}.
#' @param x \code{fit_result} object.
#' @param ... further arguments.
#' @export
print.fit_result <- function(x, ...) {
  cat("Fitting results\n", sep = "")
  cat("-------------------------------\n", sep = "")
  cat("Analysis duration : ", x$proc_time[3],"s\n", sep = "")
  cat("Number of spectra : ", nrow(stats::na.omit(x$res_tab)),"\n", sep = "")
  cat("Basis elements    : ", dim(x$basis$data)[2], "\n\n", sep = "")
  cat("Basis names\n", sep = "")
  cat("-------------------------------\n")
  cat(x$basis$names, sep = ",", fill = 31)
}

#' Print fit coordinates from a single index.
#' @param n fit index.
#' @param fit_res \code{fit_result} object.
#' @export
n2coord <- function(n, fit_res) {
  print(fit_res$res_tab[n, 1:5])
}

#' Write fit results table to a csv file.
#' @param fit_res fit result object.
#' @param fname filename of csv file.
#' @param unscaled output the unscaled result table (default = FALSE).
#' @export
fit_res2csv <- function(fit_res, fname, unscaled = FALSE) {
  if (unscaled) {
    out_tab <- fit_res$res_tab
  } else {
    out_tab <- fit_res$res_tab_unscaled
  }
  
  utils::write.csv(out_tab, fname, quote = FALSE, row.names = FALSE)
}

#' Plot a 2D slice from an MRSI fit result object.
#' @param fit_res \code{fit_result} object.
#' @param map fit result values to display as a colour map. Can be specified as
#' a character string or array of numeric values. Defaults to "tNAA".
#' @param map_denom fit result values to divide the map argument by. Can be
#' specified as a character string (eg "tCr") or array of numeric values.
#' @param slice slice to plot in the z direction.
#' @param zlim range of values to plot.
#' @param interp interpolation factor.
#' @export
plot_slice_fit <- function(fit_res, map, map_denom = NULL, slice = 1,
                           zlim = NULL, interp = 1) {
  
  if (class(map) == "character") map <- get_fit_map(fit_res, map)
  
  if (class(map_denom) == "character") map_denom <- get_fit_map(fit_res,
                                                                map_denom)
  
  if (is.null(map)) map <- get_fit_map(fit_res, "tNAA") 
  
  if (!is.null(map_denom)) map <- map / map_denom
  
  plot_map <- map[1,,, slice, 1, 1]
  plot_map <- pracma::fliplr(plot_map)
  
  col <- viridisLite::viridis(64)
  
  if (interp != 1) {
    plot_map <- mmand::rescale(plot_map, interp, mmand::mnKernel())
  }
  
  if (is.null(zlim)) { 
    fields::image.plot(plot_map, col = col, useRaster = TRUE, 
                       asp = 1, axes = FALSE, legend.shrink = 0.8)
  } else {
    plot_map <- crop_range(plot_map, zlim[1], zlim[2])
    breaks <- seq(from = zlim[1], to = zlim[2], length.out = 65)
    fields::image.plot(plot_map, col = col, useRaster = TRUE, 
                       asp = 1, axes = FALSE, legend.shrink = 0.8,
                       breaks = breaks)
  }
}

#' Get a data array from a fit result.
#' @param fit_res \code{fit_result} object.
#' @param name name of the quantity to plot, eg "tNAA".
#' @export
get_fit_map <- function(fit_res, name) {
  
  # check name is valid
  if (!(name %in% colnames(fit_res$res_tab))) {
    stop("Following column not found in fit result : ", name)
  }
  
  result_map <- fit_res$res_tab[[name]]
  dim(result_map) <- c(1, dim(fit_res$data$data)[2:6])
  result_map
}