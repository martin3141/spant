#' Print a summary of mrs_data parameters.
#' @param x mrs_data object.
#' @param ... further arguments.
#' @export
print.mrs_data <- function(x, ...) {
  cat("MRS Data Parameters\n")
  cat("----------------------------------\n")
  cat(paste(c("Trans. freq (MHz)       : ", round(x$ft * 1e-6, 4), "\n")), sep = "")
  cat(paste(c("FID data points         : ", dim(x$data)[7], "\n")), sep = "")
  cat(paste(c("X,Y,Z dimensions        : ", dim(x$data)[2], "x", dim(x$data)[3],
              "x", dim(x$data)[4], "\n")), sep = "")
  cat(paste(c("Dynamics                : ", dim(x$data)[5], "\n")), sep = "")
  cat(paste(c("Coils                   : ", dim(x$data)[6], "\n")) ,sep = "")
  cat(paste(c("Voxel resolution (mm)   : ", x$resolution[2],
              "x", x$resolution[3], "x", x$resolution[4], "\n")), sep = "")
  cat(paste(c("Sampling frequency (Hz) : ",
              1 / x$resolution[7], "\n")), sep = "")
  cat(paste(c("Reference freq. (ppm)   : ", x$ref, "\n")), sep = "")
  cat(paste(c("Spectral domain         : ", x$freq_domain[7], "\n")), sep = "")
  
  # next line only works for 1H, add nucleus option?
  #cat(paste(c("Field strength (Tesla)  : ",
  #             sprintf("%.1f",x$ft/42.58e6),"\n")),sep="")
  #cat(paste(c("Contains referece data  : ",
  #            dim(x$data)[1] == 2, "\n")), sep = "")
}

#' Plotting method for objects of class mrs_data.
#' @param x object of class mrs_data.
#' @param fd display data in the frequency-domain (default), or time-domain 
#' (logical).
#' @param x_units the units to use for the x-axis, can be one of: "ppm", "hz", 
#' "points" or "seconds".
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param y_scale option to display the y-axis values (logical).
#' @param mode representation of the complex numbers to be plotted, can be one
#' of: "real", "imag" or "abs".
#' @param dyn the dynamic index to plot.
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param coil the coil element number to plot.
#' @param lwd plot linewidth.
#' @param bty option to draw a box around the plot. See ?par.
#' @param label character string to add to the top left of the plot window.
#' @param ... other arguments to pass to the plot method.
#' @export
plot.mrs_data <- function(x, fd = TRUE, x_units = NULL, xlim = NULL,
                          y_scale = FALSE, mode = "real", dyn = 1, x_pos = 1,
                          y_pos = 1, z_pos = 1, coil = 1, lwd = NULL, 
                          bty = NULL, label = "", ...) {
  
  # convert to the correct domain for plotting
  if (fd & !is_fd(x)) {
    x <- td2fd(x)
  } else if (!fd & is_fd(x)) {
    x <- fd2td(x)
  }
  
  if (is.null(lwd)) {
    lwd <- 1.2
  }
  
  if (is.null(bty)) {
    bty <- "o"
  }
  
  if (fd) {
    xlab <- "Chemical Shift"  
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
  }
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(x)])
  }
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  #graphics::par("xaxs" = "i", "yaxs"="i") # tight axes limits
  graphics::par("xaxs" = "i") # tight axes limits
  plot_data <- x$data[1, x_pos, y_pos, z_pos, dyn, coil,]
  
  if (mode == "real") {
    plot_data <- Re(plot_data)
  } else if (mode == "imag") {
    plot_data <- Im(plot_data)
  } else if (mode == "abs") {
    plot_data <- Mod(plot_data)
  }
  
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  if (y_scale) {
    graphics::par(mar = c(3.5, 3.5, 1, 1))
    graphics::plot(x_scale[subset], plot_data[subset],type = 'l',xlim = xlim, 
                   xlab = xlab, ylab = "Intensity (au)", lwd = lwd,
                   bty = bty, ...)
  } else {
    graphics::par(mar = c(3.5, 1, 1, 1))
    graphics::plot(x_scale[subset], plot_data[subset], type = 'l', xlim = xlim,
         xlab = xlab, yaxt = "n", ylab = "", lwd = lwd, bty = bty, ...)
  }
  
  if (bty == "n") {
    graphics::abline(a = graphics::par("usr")[3], b = 0, lwd = 1.0) 
  }
  
  if (!is.null(label)) {
    max_dp <- max(plot_data[subset])
    graphics::par(xpd = T)
    graphics::text(xlim[1],max_dp * 1.03, label, cex = 2.5)
    graphics::par(xpd = F) 
  }
}

#' Image plot method for objects of class mrs_data.
#' @param x object of class mrs_data.
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param mode representation of the complex numbers to be plotted, can be one
#' of: "real", "imag" or "abs".
#' @param col Colour map to use, defaults to viridis if the package is 
#' available.
#' @param dim the dimension to display on the y-axis, can be one of: "dyn", "x",
#' "y", "z" or "coil".
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param dyn the dynamic index to plot.
#' @param coil the coil element number to plot.
#' @param ... other arguments to pass to the plot method.
#' @export
image.mrs_data <- function(x, xlim = NULL, mode = "real", col = NULL, 
                           dim = "dyn", x_pos = NULL, y_pos = NULL,
                           z_pos = NULL, dyn = 1, coil = 1, ...) { 
  if (!is_fd(x)) {
    x <- td2fd(x)
  }
  
  x_scale <- ppm(x)
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(x)])
  }
  
  graphics::par(mar = c(3.5, 3.5, 1, 1)) # margins
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  data_dim <- dim(x$data)
  
  if (is.null(x_pos)) {
    x_pos <- as.integer(data_dim[2] / 2) + 1
  }
  
  if (is.null(y_pos)) {
    y_pos <- as.integer(data_dim[3] / 2) + 1
  }
  
  if (is.null(z_pos)) {
    z_pos <- as.integer(data_dim[4] / 2) + 1
  }
  
  if (dim == "dyn") {
    plot_data <- t(x$data[1, x_pos, y_pos, y_pos, , coil, subset])
    yN <- data_dim[5]
    y_title = "Dynamic"
  } else if (dim == "x") {
    plot_data <- t(x$data[1, , y_pos, z_pos, dyn, coil, subset])
    yN <- data_dim[2]
    y_title = "x position"
  } else if (dim == "y") {
    plot_data <- t(x$data[1, x_pos, , z_pos, dyn, coil, subset])
    yN <- data_dim[3]
    y_title = "y position"
  } else if (dim == "z") {
    plot_data <- t(x$data[1, x_pos, y_pos, , dyn, coil, subset])
    yN <- data_dim[4]
    y_title = "z position"
  } else if (dim == "coil") {
    plot_data <- t(x$data[1, x_pos, y_pos, z_pos, dyn, , subset])
    yN <- data_dim[6]
    y_title = "Coil"
  } else {
    stop("Unrecognised dim value. Should be one of: dyn, x, y, z, coil")
  } 
  
  if (mode == "real") {
    plot_data <- Re(plot_data)
  } else if (mode == "imag") {
    plot_data <- Im(plot_data)
  } else if (mode == "abs") {
    plot_data <- Mod(plot_data)
  }
  
  if (is.null(col)) {
    if (is.installed("viridis")) {
      col <- viridis::viridis(64)
    } else if (is.installed("viridisLite")) {
      col <- viridisLite::viridis(64)
    } else {
     col <- grDevices::heat.colors(64)
    }
  }
  
  graphics::image(x_scale[subset][length(subset):1], (1:yN),
                  plot_data[length(subset):1,], xlim = xlim,
                  xlab = "Frequency (ppm)", ylab = y_title, 
                  col = col, ...)
}

#' Produce a plot with multiple traces.
#' @param x object for plotting.
#' @param ... Arguments to be passed to methods.
#' @export
stackplot <- function(x, ...) {
  UseMethod("stackplot", x)
}

#' Stackplot plotting method for objects of class mrs_data.
#' @param x object of class mrs_data.
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param mode representation of the complex numbers to be plotted, can be one
#' of: "real", "imag" or "abs".
#' @param col set the colour of the line, eg col = rgb(1,0,0,0.5).
#' @param x_offset separate plots in the x-axis direction by this value. 
#' Default value is 0.
#' @param y_offset separate plots in the y-axis direction by this value.
#' @param dim the dimension to stack in the y-axis direction, can be one of: 
#' "dyn", "x", "y", "z" or "coil".
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param dyn the dynamic index to plot.
#' @param coil the coil element number to plot.
#' @param ... other arguments to pass to the matplot method.
#' @export
stackplot.mrs_data <- function(x, xlim = NULL, mode = "real", col = NULL, 
                               x_offset = 0, y_offset = 5, dim = "dyn", 
                               x_pos = NULL, y_pos = NULL, z_pos = NULL, 
                               dyn = 1, coil = 1, ...) {
  
  if (!is_fd(x)) {
    x <- td2fd(x)
  }
  
  if (is.null(col)) {
    col <- 1
  }
  
  #par("xaxs" = "i") # tight axes limits
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  graphics::par(mar = c(3.5, 1, 1, 1)) # margins
  
  x_scale <- ppm(x)
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(x)])
  }
  xlim <- sort(xlim)
  
  data_dim <- dim(x$data)
  
  if (is.null(x_pos)) {
    x_pos <- as.integer(data_dim[2] / 2) + 1
  }
  
  if (is.null(y_pos)) {
    y_pos <- as.integer(data_dim[3] / 2) + 1
  }
  
  if (is.null(z_pos)) {
    z_pos <- as.integer(data_dim[4] / 2) + 1
  }
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  if (dim == "dyn") {
    plot_data <- t(x$data[1, x_pos, y_pos, y_pos, , coil, subset])
    yN <- data_dim[5]
    y_title = "Dynamic"
  } else if (dim == "x") {
    plot_data <- t(x$data[1, , y_pos, z_pos, dyn, coil, subset])
    yN <- data_dim[2]
    y_title = "x position"
  } else if (dim == "y") {
    plot_data <- t(x$data[1, x_pos, , z_pos, dyn, coil, subset])
    yN <- data_dim[3]
    y_title = "y position"
  } else if (dim == "z") {
    plot_data <- t(x$data[1, x_pos, y_pos, , dyn, coil, subset])
    yN <- data_dim[4]
    y_title = "z position"
  } else if (dim == "coil") {
    plot_data <- t(x$data[1, x_pos, y_pos, z_pos, dyn, , subset])
    yN <- data_dim[5]
    y_title = "Coil"
  } else {
    stop("Unrecognised dim value. Should be one of: dyn, x, y, z, coil")
  } 
  
  if (mode == "real") {
    plot_data <- Re(plot_data)
  } else if (mode == "imag") {
    plot_data <- Im(plot_data)
  } else if (mode == "abs") {
    plot_data <- Mod(plot_data)
  }
  
  max_val <- max(abs(plot_data))
  y_offset_vec <- 0:(ncol(plot_data) - 1) * max_val * y_offset / 100
  y_offset_mat <- matrix(y_offset_vec, nrow = nrow(plot_data),
                         ncol = ncol(plot_data), byrow = TRUE)
  
  plot_data <- plot_data + y_offset_mat
  
  x_scale_mat <- matrix(x_scale[subset], nrow = nrow(plot_data),
                        ncol = ncol(plot_data), byrow = FALSE)
  
  x_offset_mat <- matrix((0:(ncol(plot_data) - 1) * 
                          (xlim[2] - xlim[1]) * x_offset / 100), 
                         nrow = nrow(plot_data), ncol = ncol(plot_data),
                         byrow = TRUE)
  
  x_scale_mat <- x_scale_mat + x_offset_mat
  
  graphics::matplot(x_scale_mat[length(subset):1,],
                    plot_data[length(subset):1,], type = "l", 
                    lty = 1, col = col, xlab = "Frequency (PPM)", ylab = "",
                    yaxt = "n", xaxt = "n", xlim = rev(range(x_scale_mat)),
                    bty = "n", ...)
  
  graphics::axis(1, pretty(xlim))
  
  #graphics::matplot(x_scale[subset][length(subset):1],
  #                  plot_data[length(subset):1,], type = "l", xlim = xlim,
  #                  lty = 1, col = 1, xlab = "Frequency (PPM)", ylab = "",
  #                  yaxt = "n", ...)
  
  #abline(a = par("usr")[3], b = 0, lw = 2.0) # looks better for bty="n"
  
  #matplot(x_scale[subset][length(subset):1])
          #, (1:dyns(mrs_data)), plot_data[length(subset):1,],
        #xlim=xlim, xlab="Frequency (ppm)", ylab="Dynamic", 
        #col=gray.colors(64), ...)
}


#' @export
plot_slice_map <- function(data, lower = NULL, upper = NULL, mask_map = NULL,
                           mask_cutoff = 20, interp = 16, slice = 1, dyn = 1,
                           coil = 1, ref = 1, denom = NULL,
                           horizontal = FALSE) {
  
  data <- data[ref,,, slice, dyn, coil]
  data <- pracma::fliplr(data)
  data <- mmand::rescale(data, interp,mmand::mnKernel())
  
  if (!is.null(mask_map)) {
    mask_map <- mask_map[ref,,, slice, dyn, coil]
    mask_map <- pracma::fliplr(mask_map)
    mask_map <- mmand::rescale(mask_map, interp, mmand::mnKernel())
    max_mask <- max(mask_map)
    data <- ifelse(mask_map < (max_mask * mask_cutoff / 100), NA, data)
  }
  
  if (!is.null(denom)) {
    denom <- denom[ref,,, slice, dyn, coil]
    denom <- pracma::fliplr(denom)
    denom <- mmand::rescale(denom, interp,mmand::mnKernel())
    data <- data / denom
  }
  
  if (!is.null(lower) & !is.null(upper)) {
    data <- crop_range(data, lower, upper)
    breaks <- seq(from = lower, to = upper, length.out = 129)
    fields::image.plot(data, col = viridis::viridis(128), useRaster = T,
                       asp = 1, axes = F, breaks = breaks,
                       horizontal = horizontal)
    
    #image(data, col=viridis::viridis(128), useRaster = T, asp = 1, axes = F,
    #      breaks = breaks)
  } else {
    fields::image.plot(data, col = viridis::viridis(128), useRaster = T,
                       asp = 1, axes = F, horizontal = horizontal)
    
    #image(data, col = viridis::viridis(128), useRaster = T, asp = 1, axes = F)
  }
}