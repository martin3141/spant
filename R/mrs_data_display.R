#' Print a summary of mrs_data parameters.
#' @param x mrs_data object.
#' @param full print all parameters (default FALSE).
#' @param ... further arguments.
#' @export
print.mrs_data <- function(x, full = FALSE, ...) {
  
  # if (class(x)[1] == "list") { # should be the same as below
  if (inherits(x, "list", which = TRUE) == 1) {
    class(x) <- "list"
    print(x)
    return(NULL)
  }
  
  cat("MRS Data Parameters\n")
  cat("----------------------------------\n")
  cat(paste(c("Trans. freq (MHz)       : ", round(x$ft * 1e-6, 4), "\n")),
      sep = "")
  cat(paste(c("FID data points         : ", dim(x$data)[7], "\n")), sep = "")
  cat(paste(c("X,Y,Z dimensions        : ", dim(x$data)[2], "x", dim(x$data)[3],
              "x", dim(x$data)[4], "\n")), sep = "")
  cat(paste(c("Dynamics                : ", dim(x$data)[5], "\n")), sep = "")
  cat(paste(c("Coils                   : ", dim(x$data)[6], "\n")) ,sep = "")
  cat(paste(c("Voxel resolution (mm)   : ", round(x$resolution[2], 2),
              "x", round(x$resolution[3], 2),
              "x", round(x$resolution[4], 2), "\n")), sep = "")
  cat(paste(c("Sampling frequency (Hz) : ",
              1 / x$resolution[7], "\n")), sep = "")
  cat(paste(c("Reference freq. (ppm)   : ", x$ref, "\n")), sep = "")
  cat(paste(c("Nucleus                 : ", x$nuc, "\n")), sep = "")
  cat(paste(c("Spectral domain         : ", x$freq_domain[7], "\n")), sep = "")
  
  if (!is.null(x$meta$EchoTime)) {
    cat(paste(c("Echo time (s)           :", x$meta$EchoTime, "\n")),
        sep = " ")
  }
  
  if (!is.null(x$meta$RepetitionTime)) {
    cat(paste(c("Repetition time (s)     :", x$meta$RepetitionTime, "\n")),
        sep = " ")
  }
  
  if (!is.null(x$meta$Manufacturer)) {
    cat(paste(c("Manufacturer            :", x$meta$Manufacturer, "\n")),
        sep = " ")
  }
  
  if (!is.null(x$meta$PulseSequenceType)) {
    cat(paste(c("Pulse sequence type     :", x$meta$PulseSequenceType, "\n")),
        sep = " ")
  }
  
  if (!is.null(x$meta$SequenceName)) {
    cat(paste(c("Sequence name           :", x$meta$SequenceName, "\n")),
        sep = " ")
  }
  
  if (!is.null(x$meta$ProtocolName)) {
    cat(paste(c("Protocol name           :", x$meta$ProtocolName, "\n")),
        sep = " ")
  }
  
  if (full) {
    #cat(paste(c("Row vector              :", x$row_vec, "\n")), sep = " ")
    #cat(paste(c("Column vector           :", x$col_vec, "\n")), sep = " ")
    #cat(paste(c("Slice vector            :", x$sli_vec, "\n")), sep = " ")
    #cat(paste(c("Position vector         :", x$pos_vec, "\n")), sep = " ")
  }
}

#' Plotting method for objects of class mrs_data.
#' @param x object of class mrs_data.
#' @param dyn the dynamic index to plot.
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param coil the coil element number to plot.
#' @param fd display data in the frequency-domain (default), or time-domain 
#' (logical).
#' @param x_units the units to use for the x-axis, can be one of: "ppm", "hz", 
#' "points" or "seconds".
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param y_scale option to display the y-axis values (logical).
#' @param x_ax option to display the x-axis values (logical).
#' @param mode representation of the complex numbers to be plotted, can be one
#' of: "re", "im", "mod" or "arg".
#' @param lwd plot linewidth.
#' @param bty option to draw a box around the plot. See ?par.
#' @param label character string to add to the top left of the plot window.
#' @param restore_def_par restore default plotting par values after the plot has 
#' been made.
#' @param mar option to adjust the plot margins. See ?par.
#' @param xaxis_lab x-axis label.
#' @param yaxis_lab y-axis label.
#' @param xat x-axis tick label values.
#' @param xlabs x-axis tick labels.
#' @param yat y-axis tick label values.
#' @param ylabs y-axis tick labels.
#' @param show_grid plot gridlines behind the data (logical). Defaults to TRUE.
#' @param grid_nx number of cells of the grid in x and y direction. When NULL
#' the grid aligns with the tick marks on the corresponding default axis (i.e.,
#' tickmarks as computed by axTicks). When NA, no grid lines are drawn in the
#' corresponding direction.
#' @param grid_ny as above.
#' @param col set the line colour, eg col = rgb(0.5, 0.5, 0.5).
#' @param alpha set the line transparency, eg alpha = 0.5 is 50% transparency.
#' Overrides any transparency levels set by col.
#' @param bl_lty linetype for the y = 0 baseline trace. A default value NULL
#' results in no baseline being plotted.
#' @param ... other arguments to pass to the plot method.
#' @export
plot.mrs_data <- function(x, dyn = 1, x_pos = 1, y_pos = 1, z_pos = 1, coil = 1,
                          fd = TRUE, x_units = NULL, xlim = NULL, 
                          y_scale = FALSE, x_ax = TRUE, mode = "re",
                          lwd = NULL, bty = NULL, label = "",
                          restore_def_par = TRUE, mar = NULL,
                          xaxis_lab = NULL, yaxis_lab = NULL, xat = NULL,
                          xlabs = TRUE, yat = NULL, ylabs = TRUE,
                          show_grid = TRUE, grid_nx = NULL, grid_ny = NA,
                          col = NULL, alpha = NULL, bl_lty = NULL, ...) {
  
  .pardefault <- graphics::par(no.readonly = T)
 
  # remove data we don't need 
  x <- get_subset(x, x_set = x_pos, y_set = y_pos, z_set = z_pos, dyn_set = dyn,
                  coil_set = coil) 
 
  # has this data element been masked? 
  if (anyNA(x$data)) {
    graphics::par(mar = c(0, 0, 0 ,0))
    graphics::plot.new()
    # plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
    return(NULL)
  }
  
  # default to blue
  if (is.null(col)) col <- grDevices::rgb(0, 0.45, 0.7, 1)
  
  # if alpha is specified, override the default value in col
  if (!is.null(alpha)) col <- add_alpha(col, alpha)
  
  # convert to the correct domain for plotting
  if (fd & !is_fd(x)) {
    x <- td2fd(x)
  } else if (!fd & is_fd(x)) {
    x <- fd2td(x)
  }
  
  if (is.null(lwd)) lwd <- 2.0
  
  if (is.null(bty) && !y_scale) bty <- "n"
  if (is.null(bty) && y_scale) bty <- "l"
  
  if (fd) {
    xlab <- "Chemical shift"  
  } else {
    xlab <- "Time"  
  }
  
  if (is.null(x_units) & fd) {
    x_units = "ppm"
  } else if (is.null(x_units) & !fd) {
    x_units = "seconds"
  }
  
  if ( x_units == "ppm" ) {
    x_scale <- ppm(x)
    xlab <- paste(xlab, "(ppm)")
  } else if (x_units == "hz") {
    x_scale <- hz(x)
    xlab <- paste(xlab, "(Hz)")
  } else if (x_units == "points") {
    x_scale <- pts(x)
    xlab <- paste(xlab, "(Data Points)")
  } else if (x_units == "seconds") {
    x_scale <- seconds(x)
    xlab <- paste(xlab, "(s)")
  } else {
    stop("Invalid x_units option, should be one of : 'ppm', 'hz', 'points' or 'seconds'") 
  }
  
  if (!is.null(xaxis_lab)) xlab <- xaxis_lab
  
  if (is.null(yaxis_lab)) yaxis_lab <- "Intensity (au)"
  
  if (is.null(xlim)) xlim <- c(x_scale[1], x_scale[Npts(x)])
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  #graphics::par("xaxs" = "i", "yaxs"="i") # tight axes limits
  graphics::par("xaxs" = "i") # tight axes limits
  
  plot_data <- x$data[1, 1, 1, 1, 1, 1,]
  
  if (mode == "re") {
    plot_data <- Re(plot_data)
  } else if (mode == "im") {
    plot_data <- Im(plot_data)
  } else if (mode == "mod") {
    plot_data <- Mod(plot_data)
  } else if (mode == "arg") {
    plot_data <- Arg(plot_data)
  } else {
    stop("Invalid mode option, should be one of : 're', 'im', 'mod' or 'arg'") 
  }
  
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  
  if (!is.null(mar)) graphics::par(mar = mar)
  
  if (y_scale) {
    if (is.null(mar)) graphics::par(mar = c(3.5, 3.5, 1, 1))
    graphics::plot(x_scale[subset], plot_data[subset], type = 'l', xlim = xlim, 
                   xlab = xlab, ylab = yaxis_lab, lwd = lwd, bty = bty, 
                   xaxt = "n", yaxt = "n", col = col,
                   panel.first = {if (show_grid) graphics::grid(nx = grid_nx,
                                                            ny = grid_ny)}, ...)
    graphics::axis(2, lwd = 0, lwd.ticks = 1, at = yat, labels = ylabs)
  } else {
    if (is.null(mar)) graphics::par(mar = c(3.5, 1, 1, 1))
    graphics::plot(x_scale[subset], plot_data[subset], type = 'l', xlim = xlim,
         xlab = xlab, yaxt = "n", xaxt = "n", ylab = "", lwd = lwd, bty = bty,
         col = col, panel.first = {if (show_grid) graphics::grid(nx = grid_nx,
                                                      ny = grid_ny)}, ...)
  }
  
  if (x_ax) graphics::axis(1, lwd = 0, lwd.ticks = 1, at = xat, labels = xlabs)
  
  if (bty == "n") graphics::abline(h = graphics::par("usr")[3]) 
  
  if (!is.null(label)) {
    max_dp <- max(plot_data[subset])
    graphics::par(xpd = T)
    graphics::text(xlim[1], max_dp * 1.03, label, cex = 2.5)
    graphics::par(xpd = F) 
  }
  
  # draw baseline(s)
  if (!is.null(bl_lty)) {
      graphics::abline(h = 0, lty = bl_lty, lwd = 0.5)
  }
  
  if (restore_def_par) graphics::par(.pardefault)
}

#' Image plot method for objects of class mrs_data.
#' @param x object of class mrs_data.
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param mode representation of the complex numbers to be plotted, can be one
#' of: "re", "im", "mod" or "arg".
#' @param col Colour map to use, defaults to viridis.
#' @param plot_dim the dimension to display on the y-axis, can be one of: "dyn", 
#' "x", "y", "z", "coil" or NULL. If NULL (the default) all spectra will be
#' collapsed into the dynamic dimension and displayed.
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param dyn the dynamic index to plot.
#' @param coil the coil element number to plot.
#' @param restore_def_par restore default plotting par values after the plot has 
#' been made.
#' @param y_ticks a vector of indices specifying where to place tick marks.
#' @param vline draw a vertical line at the value of vline.
#' @param hline draw a horizontal line at the value of hline.
#' @param ... other arguments to pass to the plot method.
#' @export
image.mrs_data <- function(x, xlim = NULL, mode = "re", col = NULL, 
                           plot_dim = NULL, x_pos = NULL, y_pos = NULL,
                           z_pos = NULL, dyn = 1, coil = 1,
                           restore_def_par = TRUE, y_ticks = NULL, 
                           vline = NULL, hline = NULL, ...) { 
  
  .pardefault <- graphics::par(no.readonly = T)
  
  if (!is_fd(x)) x <- td2fd(x)
  
  x_scale <- ppm(x)
  
  if (is.null(xlim)) xlim <- c(x_scale[1], x_scale[Npts(x)])
  
  if (is.null(y_ticks)) {
    graphics::par(mar = c(3.5, 3.5, 1, 1)) # margins
  } else {
    graphics::par(mar = c(3.5, 3.5, 1, 1.5)) # margins
  }
  
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  if (is.null(plot_dim)) {
    x <- collapse_to_dyns(x)
    plot_dim = "dyn"
  }
  
  data_dim <- dim(x$data)
  
  if (is.null(x_pos)) x_pos <- as.integer(data_dim[2] / 2) + 1
  
  if (is.null(y_pos)) y_pos <- as.integer(data_dim[3] / 2) + 1
  
  if (is.null(z_pos)) z_pos <- as.integer(data_dim[4] / 2) + 1
  
  if (plot_dim == "dyn") {
    plot_data <- t(x$data[1, x_pos, y_pos, y_pos, , coil, subset])
    yN <- data_dim[5]
    y_title = "Dynamic"
  } else if (plot_dim == "x") {
    plot_data <- t(x$data[1, , y_pos, z_pos, dyn, coil, subset])
    yN <- data_dim[2]
    y_title = "x position"
  } else if (plot_dim == "y") {
    plot_data <- t(x$data[1, x_pos, , z_pos, dyn, coil, subset])
    yN <- data_dim[3]
    y_title = "y position"
  } else if (plot_dim == "z") {
    plot_data <- t(x$data[1, x_pos, y_pos, , dyn, coil, subset])
    yN <- data_dim[4]
    y_title = "z position"
  } else if (plot_dim == "coil") {
    plot_data <- t(x$data[1, x_pos, y_pos, z_pos, dyn, , subset])
    yN <- data_dim[6]
    y_title = "Coil"
  } else {
    stop("Unrecognised dim value. Should be one of: dyn, x, y, z, coil")
  } 
  
  if (nrow(plot_data) == 1) {
    warning("image is designed for plotting multiple spectra but only one has been selected.")
    plot_data <- t(plot_data)
  }
  
  if (mode == "re") {
    plot_data <- Re(plot_data)
  } else if (mode == "im") {
    plot_data <- Im(plot_data)
  } else if (mode == "mod") {
    plot_data <- Mod(plot_data)
  } else if (mode == "arg") {
    plot_data <- Arg(plot_data)
  } else {
    stop("Invalid mode option, should be one of : 're', 'im', 'mod' or 'arg'") 
  }
  
  # remove any columns with NAs
  # plot_data <- t(stats::na.omit(t(plot_data)))
  
  # set masked spectra to zero
  plot_data[is.na(plot_data)] <- 0
  
  
  yN <- ncol(plot_data)
  
  col <- viridisLite::viridis(128)
  
  graphics::image(x_scale[subset][length(subset):1], (1:yN),
                  plot_data[length(subset):1,,drop = F], xlim = xlim,
                  xlab = "Chemical shift (ppm)", ylab = y_title, 
                  col = col, ...)
  
  if (!is.null(y_ticks)) {
    graphics::axis(4, at = y_ticks, labels = F, col = NA, col.ticks = "red")
    graphics::axis(2, at = y_ticks, labels = F, col = NA, col.ticks = "red")
  }
  
  if (!is.null(vline)) graphics::abline(v = vline, col = "white")
  
  if (!is.null(hline)) graphics::abline(h = hline, col = "white")
  
  if (restore_def_par) graphics::par(.pardefault)
}

#' Produce a plot with multiple traces.
#' @param x object for plotting.
#' @param ... arguments to be passed to methods.
#' @export
stackplot <- function(x, ...) {
  UseMethod("stackplot", x)
}

#' @export
stackplot.list <- function(x, ...) {
  # make them all td or fd
  combined <- append_scan(x)
  stackplot(combined, plot_dim = "scan", ...)
}

#' Stackplot plotting method for objects of class mrs_data.
#' @param x object of class mrs_data.
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param mode representation of the complex numbers to be plotted, can be one
#' of: "re", "im", "mod" or "arg".
#' @param fd display data in the frequency-domain (default), or time-domain 
#' (logical).
#' @param x_units the units to use for the x-axis, can be one of: "ppm", "hz", 
#' "points" or "seconds".
#' @param col set the colour of the line, eg col = rgb(1, 0, 0, 0.5).
#' @param alpha set the line transparency, eg alpha = 0.5 is 50% transparency.
#' Overrides any transparency levels set by col.
#' @param x_offset separate plots in the x-axis direction by this value. 
#' Default value is 0.
#' @param y_offset separate plots in the y-axis direction by this value.
#' @param plot_dim the dimension to display on the y-axis, can be one of: "dyn", 
#' "x", "y", "z", "coil" or NULL. If NULL (the default) all spectra will be
#' collapsed into the dynamic dimension and displayed.
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param dyn the dynamic index to plot.
#' @param coil the coil element number to plot.
#' @param bty option to draw a box around the plot. See ?par.
#' @param labels add labels to each data item.
#' @param lab_cex label size.
#' @param right_marg change the size of the right plot margin.
#' @param bl_lty linetype for the y = 0 baseline trace. A default value NULL
#' results in no baseline being plotted.
#' @param restore_def_par restore default plotting par values after the plot has 
#' been made.
#' @param show_grid plot gridlines behind the data (logical). Defaults to TRUE.
#' @param grid_nx number of cells of the grid in x and y direction. When NULL
#' the grid aligns with the tick marks on the corresponding default axis (i.e.,
#' tickmarks as computed by axTicks). When NA, no grid lines are drawn in the
#' corresponding direction.
#' @param grid_ny as above.
#' @param lwd plot linewidth.
#' @param ... other arguments to pass to the matplot method.
#' @export
stackplot.mrs_data <- function(x, xlim = NULL, mode = "re", x_units = NULL,
                               fd = TRUE, col = NULL, alpha = NULL, 
                               x_offset = 0, y_offset = 0, plot_dim = NULL,
                               x_pos = NULL, y_pos = NULL, z_pos = NULL,
                               dyn = 1, coil = 1, bty = NULL, labels = NULL,
                               lab_cex = 1, right_marg = NULL, bl_lty = NULL,
                               restore_def_par = TRUE, show_grid = NULL,
                               grid_nx = NULL, grid_ny = NA, lwd = NULL, ...) {
  
  .pardefault <- graphics::par(no.readonly = T)
  
  # convert to the correct domain for plotting
  if (fd & !is_fd(x)) {
    x <- td2fd(x)
  } else if (!fd & is_fd(x)) {
    x <- fd2td(x)
  }
  
  # default to blue
  if (is.null(col)) col <- grDevices::rgb(0, 0.45, 0.7, 1)
  
  # show the grid by default unless x_offset is non-zero
  if (is.null(show_grid)) {
    if (x_offset == 0) {
      show_grid = TRUE 
    } else {
      show_grid = FALSE
    }
  }
  
  # if alpha is specified, override the default value in col
  if (!is.null(alpha)) col <- add_alpha(col, alpha)
  
  if (is.null(bty)) bty <- "n"
  
  if (is.null(lwd)) lwd <- 2.0
  
  if (is.null(right_marg) && is.null(labels)) right_marg = 1
  if (is.null(right_marg) && !is.null(labels)) right_marg = 4
  
  graphics::par("xaxs" = "i") # tight axes limits
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  graphics::par(mar = c(3.5, 1, 1, right_marg)) # margins
  
  if (fd) {
    xlab <- "Chemical shift"  
  } else {
    xlab <- "Time"  
  }
  
  if (is.null(x_units) & fd) {
    x_units = "ppm"
  } else if (is.null(x_units) & !fd) {
    x_units = "seconds"
  }
  
  if ( x_units == "ppm" ) {
    x_scale <- ppm(x)
    xlab <- paste(xlab, "(ppm)")
  } else if (x_units == "hz") {
    x_scale <- hz(x)
    xlab <- paste(xlab, "(Hz)")
  } else if (x_units == "points") {
    x_scale <- pts(x)
    xlab <- paste(xlab, "(Data Points)")
  } else if (x_units == "seconds") {
    x_scale <- seconds(x)
    xlab <- paste(xlab, "(s)")
  } else {
    stop("Invalid x_units option, should be one of : 'ppm', 'hz', 'points' or 'seconds'") 
  }
  
  if (is.null(xlim)) xlim <- c(x_scale[1], x_scale[Npts(x)])
  
  xlim <- sort(xlim)
  
  if (is.null(plot_dim)) {
    x <- collapse_to_dyns(x)
    plot_dim = "dyn"
  }
  
  data_dim <- dim(x$data)
  
  if (is.null(x_pos)) x_pos <- as.integer(data_dim[2] / 2) + 1
  
  if (is.null(y_pos)) y_pos <- as.integer(data_dim[3] / 2) + 1
  
  if (is.null(z_pos)) z_pos <- as.integer(data_dim[4] / 2) + 1
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  if (plot_dim == "dyn") {
    plot_data <- t(x$data[1, x_pos, y_pos, z_pos, , coil, subset])
    yN <- data_dim[5]
    y_title = "Dynamic"
  } else if (plot_dim == "x") {
    plot_data <- t(x$data[1, , y_pos, z_pos, dyn, coil, subset])
    yN <- data_dim[2]
    y_title = "x position"
  } else if (plot_dim == "y") {
    plot_data <- t(x$data[1, x_pos, , z_pos, dyn, coil, subset])
    yN <- data_dim[3]
    y_title = "y position"
  } else if (plot_dim == "z") {
    plot_data <- t(x$data[1, x_pos, y_pos, , dyn, coil, subset])
    yN <- data_dim[4]
    y_title = "z position"
  } else if (plot_dim == "coil") {
    plot_data <- t(x$data[1, x_pos, y_pos, z_pos, dyn, , subset])
    yN <- data_dim[5]
    y_title = "Coil"
  } else if (plot_dim == "scan") {
    plot_data <- t(x$data[, x_pos, y_pos, z_pos, dyn, coil, subset])
    yN <- data_dim[1]
    y_title = "Scan"
  } else {
    stop("Unrecognised dim value. Should be one of: dyn, x, y, z, coil")
  }
  
  if (nrow(plot_data) == 1) {
    warning("stackplot is designed for plotting multiple spectra but only one has been selected.")
    plot_data <- t(plot_data)
  }
  
  if (mode == "re") {
    plot_data <- Re(plot_data)
  } else if (mode == "im") {
    plot_data <- Im(plot_data)
  } else if (mode == "mod") {
    plot_data <- Mod(plot_data)
  } else if (mode == "arg") {
    plot_data <- Arg(plot_data)
  } else {
    stop("Invalid mode option, should be one of : 're', 'im', 'mod' or 'arg'") 
  }
  
  # remove any columns with NAs
  plot_data <- t(stats::na.omit(t(plot_data)))
  
  max_val <- max(abs(plot_data), na.rm = TRUE)
  y_offset_vec <- 0:(ncol(plot_data) - 1) * max_val * y_offset / 100
  y_offset_mat <- matrix(y_offset_vec, nrow = nrow(plot_data),
                         ncol = ncol(plot_data), byrow = TRUE)
  
  plot_data <- plot_data + y_offset_mat
  
  # linetype for spectral data
  lty <- rep(1, ncol(plot_data))
  
  # add baseline traces to the plot?
  #if (!is.null(bl_lty)) {
    # only need one baseline trace if y_offset is zero
  #  if (y_offset == 0) y_offset_mat <- y_offset_mat[, 1, drop = FALSE]
    
  #  plot_data <- cbind(plot_data, y_offset_mat)
  #  lty <- c(lty, rep(bl_lty, ncol(y_offset_mat)))
    # col <- c(col, rep("black", ncol(y_offset_mat)))
  #}
  
  x_scale_mat <- matrix(x_scale[subset], nrow = nrow(plot_data),
                        ncol = ncol(plot_data), byrow = FALSE)
  
  x_offset_mat <- matrix((0:(ncol(plot_data) - 1) * 
                         (xlim[2] - xlim[1]) * -x_offset / 100), 
                         nrow = nrow(plot_data), ncol = ncol(plot_data),
                         byrow = TRUE)
  
  x_scale_mat <- x_scale_mat + x_offset_mat
  
  xlim_labs <- xlim
  xlim <- range(x_scale_mat)
  
  # bug fix for rounding errors
  if (xlim[1] > xlim_labs[1]) xlim[1] <- xlim_labs[1]
  if (xlim[2] < xlim_labs[2]) xlim[2] <- xlim_labs[2]
  
  if ( x_units == "ppm" ) xlim <- rev(xlim)
  
  graphics::matplot(x_scale_mat[length(subset):1,],
                    plot_data[length(subset):1,], type = "l", 
                    lty = lty, col = col, xlab = xlab, ylab = "",
                    yaxt = "n", xaxt = "n", xlim = xlim,
                    bty = bty, lwd = lwd, panel.first = {if (show_grid)
                                                   graphics::grid(nx = grid_nx,
                                                   ny = grid_ny)}, ...)
  
  graphics::axis(1, lwd = 0, lwd.ticks = 1, at = pretty(xlim_labs))
  
  if (bty == "n") graphics::lines(xlim_labs, c(graphics::par("usr")[3],
                                               graphics::par("usr")[3]))
  
  # if (x_offset != 0) {
  #   # graphics::lines(c(0, utils::tail(as.numeric(x_offset_mat),1)),
  #   #       c(graphics::par("usr")[3], graphics::par("usr")[3] +
  #   #           utils::tail(as.numeric(y_offset_mat),1)))
  #   
  #   graphics::lines(c(0, utils::tail(as.numeric(x_offset_mat),1)),
  #             c(0, utils::tail(as.numeric(y_offset_mat),1)), lty = 2)
  # }
  
  # write text labels if provided
  if (!is.null(labels)) {
    
    # allow text outside axes
    graphics::par(xpd = T)
    for (n in 1:length(labels)) {
      graphics::text(xlim[2], y_offset_vec[n], labels[n], pos = 4,
                     cex = lab_cex)
    }
    graphics::par(xpd = F)
  }
  
  # draw baseline(s)
  if (!is.null(bl_lty)) {
    # only need one baseline trace if y_offset is zero
    if (y_offset == 0) {
      graphics::abline(h = 0, lty = bl_lty, lwd = 0.5)
    } else {
      for (offset in y_offset_vec) {
        graphics::abline(h = offset, lty = bl_lty, lwd = 0.5)
      }
    }
  }
  
  # if (show_grid) graphics::grid(nx = grid_nx, ny = grid_ny)
  
  if (restore_def_par) graphics::par(.pardefault)
}

#' Convenience function to plot a baseline estimate with the original data.
#' @param orig_data the original data.
#' @param bc_data the baseline corrected data.
#' @param ... other arguments to pass to the stackplot function.
#' @export
plot_bc <- function(orig_data, bc_data, ...) {
  bl <- orig_data - bc_data
  combined <- append_scan(orig_data, bl)
  stackplot(combined, dim = "scan", ...)
}

#' Plot a slice from a 7 dimensional array.
#' @param data 7d array of values to be plotted.
#' @param zlim smallest and largest values to be plotted.
#' @param mask_map matching map with logical values to indicate if the 
#' corresponding values should be plotted.
#' @param mask_cutoff minimum values to plot (as a percentage of the maximum).
#' @param interp map interpolation factor.
#' @param slice the slice index to plot.
#' @param dyn the dynamic index to plot.
#' @param coil the coil element number to plot.
#' @param ref reference index to plot.
#' @param denom map to use as a denominator.
#' @param horizontal display the colourbar horizontally (logical).
#' @export
plot_slice_map <- function(data, zlim = NULL, mask_map = NULL,
                           mask_cutoff = 20, interp = 1, slice = 1, dyn = 1,
                           coil = 1, ref = 1, denom = NULL,
                           horizontal = FALSE) {
  
  if (inherits(data, "mrs_data")) {
    data <- get_subset(data, coil_set = coil) # speeds things up
    data <- int_spec(data, mode = "mod")
  }
  
  graphics::par(mar = c(0, 0, 0, 2))
  
  data_mask <- is.na(data)
  if (any(data_mask)) {
    data[data_mask] <- 0 # set NAs to zero
  } else {
    data_mask <- NULL
  }
  
  data <- data[ref,,, slice, dyn, coil]
  data <- pracma::fliplr(data) # ?
  data <- mmand::rescale(data, interp, mmand::mnKernel())
  
  if (!is.null(mask_map)) {
    mask_map <- mask_map[ref,,, slice, dyn, coil]
    mask_map <- pracma::fliplr(mask_map) # ?
    mask_map <- mmand::rescale(mask_map, interp, mmand::mnKernel())
    max_mask <- max(mask_map)
    data <- ifelse(mask_map < (max_mask * mask_cutoff / 100), NA, data)
  }
  
  if (!is.null(data_mask)) {
    data_mask <- data_mask[ref,,, slice, dyn, coil]
    data_mask <- pracma::fliplr(data_mask)
    data_mask <- mmand::rescale(data_mask, interp, mmand::boxKernel())
    data[data_mask == 1] <- NA
  }
  
  if (!is.null(denom)) {
    denom <- denom[ref,,, slice, dyn, coil]
    denom <- pracma::fliplr(denom) # ?
    denom <- mmand::rescale(denom, interp,mmand::mnKernel())
    data <- data / denom
  }
  
  asp <- ncol(data) / nrow(data)
  
  if (!is.null(zlim)) {
    data <- crop_range(data, zlim[1], zlim[2])
    breaks <- seq(from = zlim[1], to = zlim[2], length.out = 129)
    fields::image.plot(data, col = viridisLite::viridis(128), useRaster = T,
                       asp = asp, axes = F, breaks = breaks,
                       horizontal = horizontal, legend.shrink = 0.5,
                       legend.mar = 7)
    
    #image(data, col=viridisLite::viridis(128), useRaster = T, asp = 1, axes = F,
    #      breaks = breaks)
  } else {
    fields::image.plot(data, col = viridisLite::viridis(128), useRaster = T,
                       asp = asp, axes = F, horizontal = horizontal,
                       legend.shrink = 0.5, legend.mar = 7)
    
    #image(data, col = viridis::viridisLite(128), useRaster = T, asp = 1, axes = F)
  }
}

#' Arrange spectral plots in a grid.
#' @param x object for plotting.
#' @param ... arguments to be passed to methods.
#' @export
gridplot <- function(x, ...) {
  UseMethod("gridplot", x)
}

#' @export
gridplot.list <- function(x, ...) {
  # make them all td or fd
  combined <- append_scan(x)
  gridplot(combined, ...)
}

#' Arrange spectral plots in a grid.
#' @param x object of class mrs_data.
#' @param rows number of grid rows.
#' @param cols number of grid columns.
#' @param mar option to adjust the plot margins. See ?par.
#' @param oma outer margin area.
#' @param bty option to draw a box around the plot. See ?par.
#' @param restore_def_par restore default plotting par values after the plot has 
#' been made.
#' @param ... other arguments to pass to the plot method.
#' @export
gridplot.mrs_data <- function(x, rows = NA, cols = NA, mar = c(0, 0, 0, 0),
                              oma = c(3.5, 1, 1, 1), bty = "o",
                              restore_def_par = TRUE, ...) {
  
  .pardefault <- graphics::par(no.readonly = T)
  
  mrs_data_dyns <- collapse_to_dyns(x)
  Nspec <- Ndyns(mrs_data_dyns)
  
  # auto choose rows and cols for MRSI
  if (is.na(rows) && (Nx(x) + Ny(x) > 2)) rows <- Ny(x)
  if (is.na(cols) && (Nx(x) + Ny(x) > 2)) cols <- Nx(x)
  
  # set to rows and cols to be squareish if not specified
  if (is.na(rows) & is.na(cols)) {
    rows <- ceiling(Nspec ^ 0.5)
    cols <- ceiling(Nspec / rows)
  } else if (is.na(rows)) {
    rows <- ceiling(Nspec / cols)
  } else if (is.na(cols)) {
    cols <- ceiling(Nspec / rows)
  }
  
  graphics::par(mfrow = c(rows, cols), oma = oma)
  
  if (Ndyns(mrs_data_dyns) > rows * cols) {
    warning("not enough rows and columns to show all spectra")
    mrs_data_dyns <- get_dyns(mrs_data_dyns, 1:(rows * cols))
  }
  
  for (n in 1:Ndyns(mrs_data_dyns)) {
    if (n > (rows * cols - cols)) {
      x_ax = TRUE
    } else {
      x_ax = FALSE
    }
      
    graphics::plot(mrs_data_dyns, restore_def_par = FALSE, dyn = n, mar = mar,
                   bty = bty, x_ax = x_ax, ...)
  }
  
  graphics::mtext(text="Chemical shift (ppm)", side = 1, line = 1.8, outer=TRUE,
                  cex = 0.8)
  
  if (restore_def_par) graphics::par(.pardefault)
}