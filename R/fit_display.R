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
#' @param ... further arguments to plot method.
#' @export
plot.fit_result <- function(x, dyn = 1, x_pos = 1, y_pos = 1, z_pos = 1,
                            coil = 1,xlim = NULL, data_only = FALSE,
                            label = NULL, plot_sigs = NULL, n = NULL,
                            sub_bl = FALSE, mar = NULL, restore_def_par = TRUE, 
                            ylim = NULL, y_scale = FALSE, ...) {
  
  .pardefault <- graphics::par(no.readonly = T)
  
  if (is.null(n)) {
    ind <- (x$res_tab$X == x_pos) & (x$res_tab$Y == y_pos) & 
           (x$res_tab$Z == z_pos) & (x$res_tab$Dynamic == dyn) &
           (x$res_tab$Coil == coil) 
    
    n <- which(ind)
  }
  
  x <- x$fits[[n]]
  
  if (anyNA(x)) { 
    plot.new()
    return(NULL)
  }
  
  if (is.null(xlim)) {
    xlim <- rev(range(x$PPMScale))
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
         xlab = "Chemical shift (ppm)", ...)
    
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
         xlab = "Chemical shift (ppm)", ...)
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
#' @param ... further arguments to plot method.
#' @export
stackplot.fit_result <- function(x, xlim = NULL, y_offset = 1, dyn = 1, 
                                 x_pos = 1, y_pos = 1, z_pos = 1, coil = 1,
                                 n = NULL, sub_bl = FALSE, labels = FALSE,
                                 label_names = NULL, sig_col = "black",
                                 restore_def_par = TRUE, omit_signals = NULL,
                                 combine_lipmm = FALSE, combine_metab = FALSE,
                                 mar = NULL, ...) {
  
  .pardefault <- graphics::par(no.readonly = T)
  
  if (is.null(n)) {
    ind <- (x$res_tab$X == x_pos) & (x$res_tab$Y == y_pos) & 
           (x$res_tab$Z == z_pos) & (x$res_tab$Dynamic == dyn) &
           (x$res_tab$Coil == coil) 
    
    n <- which(ind)
  }
  
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
    xlim <- rev(range(x$PPMScale))
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
       ylim = c(min_basis - basis_yoff, max_dp + res_range), yaxt = "n", ylab = "",
       xlab = "Chemical shift (ppm)", ...)
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
  cat("--------------------------\n", sep = "")
  cat("Analysis duration : ", x$proc_time[3],"s\n", sep = "")
  #cat("Mean residual     : ", mean(x$res_tab$res.deviance),"\n", sep = "")
  #cat("Mean iterations   : ", mean(x$res_tab$res.niter),"\n", sep = "")
  cat("Number of spectra : ", Nspec(x$data),"\n", sep = "")
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
#' @param x fit results table.
#' @param fname filename of csv file.
#' @param pvc output PVC or raw results (logical).
#' @export
fit_tab2csv <- function(x, fname, pvc = FALSE) {
  utils::write.csv(x, fname, quote = FALSE, row.names = FALSE)
}

#' Plot a 2D slice from an MRSI fit result object.
#' @param fit_res \code{fit_result} object.
#' @param name name of the quantity to plot, eg "TNAA".
#' @param slice slice to plot in the z direction.
#' @param zlim range of values to plot.
#' @param interp interpolation factor.
#' @export
plot_slice_fit <- function(fit_res, name, slice = 1, zlim = NULL, interp = 1) {
  result_map <- fit_res$res_tab[[name]]
  dim(result_map) <- dim(fit_res$data$data)[2:6]
  col <- viridisLite::viridis(64)
  plot_map <- result_map[,, slice, 1, 1]
  plot_map <- pracma::fliplr(plot_map)
  plot_map <- mmand::rescale(plot_map, interp, mmand::mnKernel())
  
  if (is.null(zlim)) { 
    fields::image.plot(plot_map, col = col, useRaster = TRUE, 
                       asp = 1, axes = FALSE, legend.shrink = 0.8)
  } else {
    plot_map <- crop_range(plot_map, zlim[1], zlim[2])
    breaks <- seq(from = zlim[1], to = zlim[2], length.out = 65)
    fields::image.plot(plot_map, col = col, useRaster = TRUE, 
                       asp = 1, axes = FALSE, legend.shrink = 0.8, breaks = breaks)
  }
}

#' Get a data array from a fit result.
#' @param fit_res \code{fit_result} object.
#' @param name name of the quantity to plot, eg "TNAA".
#' @export
get_fit_map <- function(fit_res, name) {
  result_map <- fit_res$res_tab[[name]]
  dim(result_map) <- c(1, dim(fit_res$data$data)[2:6])
  result_map
}