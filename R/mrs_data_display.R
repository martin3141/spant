#' @export
print.mrs_data <- function(x, ...) {
  cat("MRS data parameters\n")
  cat("-------------------------------\n")
  # next line only works for 1H, add nucleus option?
  #cat(paste(c("Field strength (Tesla)  : ",
  #             sprintf("%.1f",x$ft/42.58e6),"\n")),sep="")
  cat(paste(c("Trans. freq (MHz)       : ", x$ft * 1e-6, "\n")), sep = "")
  cat(paste(c("FID data points         : ", dim(x$data)[7], "\n")), sep = "")
  cat(paste(c("X,Y,Z dimensions        : ", dim(x$data)[2], "x", dim(x$data)[3],
              "x", dim(x$data)[4], "\n")), sep = "")
  
  cat(paste(c("Dynamics                : ", dim(x$data)[5], "\n")), sep = "")
  cat(paste(c("Coils                   : ", dim(x$data)[6], "\n")) ,sep = "")
  cat(paste(c("Voxel resolution (mm)   : ", x$resolution[2],
              "x", x$resolution[3], "x", x$resolution[4], "\n")), sep = "")
  cat(paste(c("Sampling frequency (Hz) : ",
              1 / x$resolution[7], "\n")), sep = "")
  
  cat(paste(c("Contains referece data  : ",
              dim(x$data)[1] == 2, "\n")), sep = "")
  
  cat(paste(c("Spectral domain         : ", x$freq_domain[7], "\n")), sep = "")
  cat(paste(c("Reference freq. (PPM)   : ", x$ref, "\n")), sep = "")
}

# TODO add dim option
#' @export
image.mrs_data <- function(x, mode = "real", xlim = NULL, col = NULL, ...) { 
  if (!is_fd(x)) {
    x <- td2fd(x)
  }
  
  x_scale <- ppm(x)
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(x)])
  }
  
  par(mar = c(3.5, 3.5, 1, 1)) # margins
  par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  
  x_inds <- get_seg_ind(x_scale, xlim[1], xlim[2])
  subset <- x_inds[1]:x_inds[2]
  plot_data <- t(x$data[1, 1, 1, 1, , 1, subset])
  
  if (mode == "real") {
    plot_data <- Re(plot_data)
  } else if (mode == "imag") {
    plot_data <- Im(plot_data)
  } else if (mode == "abs") {
    plot_data <- Mod(plot_data)
  }
  
  if (is.null(col)) {
    is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
    if (is.installed("viridis")) {
      col <- viridis::viridis(64)
    } else if (is.installed("viridisLite")) {
      col <- viridisLite::viridis(64)
    } else {
     col <- grDevices::heat.colors(64)
    }
  }
  
  graphics::image(x_scale[subset][length(subset):1], (1:dyns(x)),
                  plot_data[length(subset):1,], xlim = xlim,
                  xlab = "Frequency (ppm)", ylab = "Dynamic", 
                  col = col, ...)
}

#' Produce a plot with multiple traces.
#' @param x object for plotting.
#' @param ... Arguments to be passed to methods.
#' @export
stackplot <- function(x, ...) {
  UseMethod("stackplot", x)
}

# TODO make consistant with plot
#' @export
stackplot.mrs_data <- function(x, mode = "real", xlim = NULL,
                               x_offset = 5, ...) {
  if (!is_fd(x)) {
    x <- td2fd(x)
  }
  
  par("xaxs" = "i") # tight axes limits
  par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  par(mar = c(3.5, 1, 1, 1)) # margins
  
  x_scale <- ppm(x)
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(x)])
  }
  
  x_inds <- get_seg_ind(x_scale, xlim[1], xlim[2])
  subset <- x_inds[1]:x_inds[2]
  plot_data <- t(x$data[1, 1, 1, 1,, 1, subset])
  
  if (mode == "real") {
    plot_data <- Re(plot_data)
  } else if (mode == "imag") {
    plot_data <- Im(plot_data)
  } else if (mode == "abs") {
    plot_data <- Mod(plot_data)
  }
  
  max_val <- max(abs(plot_data))
  x_offset_vec <- 0:(ncol(plot_data) - 1) * max_val * x_offset / 100
  x_offset_mat <- matrix(x_offset_vec, nrow = nrow(plot_data),
                         ncol = ncol(plot_data), byrow = TRUE)
  
  plot_data <- plot_data + x_offset_mat
  
  graphics::matplot(x_scale[subset][length(subset):1],
                    plot_data[length(subset):1,], type = "l", xlim = xlim,
                    lty = 1, col = 1, xlab = "Frequency (PPM)", ylab = "",
                    yaxt = "n", ...)
  
  abline(a = par("usr")[3], b = 0, lw = 2.0) # looks better for bty="n"
  
  #matplot(x_scale[subset][length(subset):1])
          #, (1:dyns(mrs_data)), plot_data[length(subset):1,],
        #xlim=xlim, xlab="Frequency (ppm)", ylab="Dynamic", 
        #col=gray.colors(64), ...)
}

#' @importFrom graphics par plot abline text                   
#' @export
plot.mrs_data <- function(x, fd = TRUE, scale = NULL, xlim = NULL,
                          yscale = FALSE, mode = "real", x_pos = 1, y_pos = 1,
                          z_pos = 1, dyn = 1, coil = 1, ref = 1, label = "", 
                          ...) {
  
  # convert to the correct domain for plotting
  if (fd & !is_fd(x)) {
    x <- td2fd(x)
  } else if (!fd & is_fd(x)) {
    x <- fd2td(x)
  }
  
  if ( fd ) {
    xlab <- "Chemical Shift"  
  } else {
    xlab <- "Time"  
  }
  
  if (is.null(scale) & fd) {
    scale = "ppm"
  } else if (is.null(scale) & !fd) {
    scale = "seconds"
  }
  
  if ( scale == "ppm" ) {
    x_scale <- ppm(x)
    xlab <- paste(xlab, "(ppm)")
  } else if (scale == "hz") {
    x_scale <- hz(x)
    xlab <- paste(xlab, "(Hz)")
  } else if (scale == "points") {
    x_scale <- pts(x)
    xlab <- paste(xlab, "(Data Points)")
  } else if (scale == "seconds") {
    x_scale <- seconds(x)
    xlab <- paste(xlab, "(s)")
  }
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(x)])
  }
  
  x_inds <- get_seg_ind(x_scale, xlim[1], xlim[2])
  subset <- x_inds[1]:x_inds[2]
  
  #par("xaxs" = "i", "yaxs"="i") # tight axes limits
  par("xaxs" = "i") # tight axes limits
  plot_data <- x$data[ref, x_pos, y_pos, z_pos, dyn, coil,]
  
  if (mode == "real") {
    plot_data <- Re(plot_data)
  } else if (mode == "imag") {
    plot_data <- Im(plot_data)
  } else if (mode == "abs") {
    plot_data <- Mod(plot_data)
  }
  
  par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  if (yscale) {
    par(mar = c(3.5, 3.5, 1, 1))
    plot(x_scale[subset], plot_data[subset],type = 'l',xlim = xlim, xlab = xlab,
         ylab = "Intensity (au)", ...)
  } else {
    par(mar = c(3.5, 1, 1, 1))
    plot(x_scale[subset], plot_data[subset], type = 'l', xlim = xlim,
         xlab = xlab, yaxt = "n", ylab = "", ...)
  }
  abline(a = par("usr")[3], b = 0, lw = 2.0) # looks better for bty="n"
  
  if (!is.null(label)) {
    max_dp <- max(plot_data[subset])
    par(xpd = T)
    text(xlim[1],max_dp * 1.03, label, cex = 2.5)
    par(xpd = F) 
  }
}

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