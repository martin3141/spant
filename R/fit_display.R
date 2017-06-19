#' Plot the fitting results of an object of class \code{fit_result}.
#' @param x fit_result object.
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param plt_title title to add to the plot.
#' @param data_only display only the processed data (logical).
#' @param label character string to add to the top left of the plot window.
#' @param plot_sigs a character vector of signal names to add to the plot.
#' @param n index to the fit number to display.
#' @param ... further arguments to plot method.
#' @export
plot.fit_result <- function(x, xlim = c(4, 0.5), plt_title = FALSE,
                           data_only = FALSE, label=NULL, 
                           plot_sigs = NULL, n = NA, ...) {
  
  if ( is.na(n) && length(x$fits) > 1 ) {
    warning("Fit number n not specified, plotting the first one.")
    n = 1
  }
  
  if ( is.na(n)) {n = 1} # SVS case
  
  x <- x$fits[[n]]
  
  graphics::par("xaxs" = "i", "yaxs" = "i") # tight axes limits
  
  if ( plt_title == FALSE ) {
    graphics::par(mar = c(3.5, 1.2, 1.2, 1.2)) # space around the plot
  } else {
    graphics::par(mar = c(3.5, 1.2, 2, 1.2)) # space around the plot
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
         ylim = c(min_dp - marg, max_dp + marg), yaxt = "n", ylab = "",
         xlab = "Chemical Shift (ppm)", ...)
    
    if (!is.null(label)) {
      graphics::par(xpd = T)
      graphics::text(xlim[1], max_dp, label, cex = 2.5)
      graphics::par(xpd = F) 
    }
    
  } else {
    fit_line <- x$Fit + x$Baseline
    max_dp <- max(x$Data[ind],fit_line[ind])
    min_dp <- min(x$Data[ind],fit_line[ind],x$Baseline[ind])
    
    res <- x$Data - fit_line
    res_range <- max(res[ind]) - min(res[ind])
    offset <- max_dp - min(res[ind])
    
    graphics::plot(x$PPMScale, x$Data, type = 'l', xlim = xlim, 
         ylim = c(min_dp,max_dp + res_range), yaxt = "n", ylab = "",
         xlab = "Chemical Shift (ppm)", ...)
    graphics::lines(x$PPMScale, fit_line, col = 'Red', lw = 2)
    graphics::lines(x$PPMScale, x$Baseline)
    graphics::lines(x$PPMScale, res + offset)
    graphics::abline(h = max_dp)
  }
  
  for (sig in plot_sigs) {
    graphics::lines(x$PPMScale, x[sig][[1]] + x$Baseline, col = "blue")
  }
}

#' Plot the fitting results of an object of class \code{fit_result} with 
#' individual basis set components shown.
#' @param x fit_result object.
#' @param xlim the range of values to display on the x-axis, eg xlim = c(4,1).
#' @param y_offset separate basis signals in the y-axis direction by this value.
#' @param plt_title title to add to the plot.
#' @param n index to the fit number to display.
#' @param ... further arguments to plot method.
#' @export
stackplot.fit_result <- function(x, xlim = c(4, 0.5), y_offset = 0.04, 
                                 plt_title = FALSE, n = NA, ...) {
  
  if ( is.na(n) && length(x$fits) > 1 ) {
    warning("Fit number n not specified, plotting the first one.")
    n = 1
  }
  
  if ( is.na(n)) {n = 1} # SVS case
  
  x <- x$fits[[n]]
  
  graphics::par("xaxs" = "i", "yaxs" = "i") # tight axes limits
  
  if ( plt_title == FALSE ) {
    graphics::par(mar = c(3.5, 1.2, 1.2, 1.2)) # space around the plot
  } else {
    graphics::par(mar = c(3.5, 1.2, 2, 1.2)) # space around the plot
  }
  
  graphics::par(mgp = c(1.8, 0.5, 0)) # distance between axes and labels
  
  if (xlim[1] > xlim[2]) {
    ind <- x$PPMScale < xlim[1] & x$PPMScale > xlim[2]
  }
  else {
    ind <- x$PPMScale < xlim[2] & x$PPMScale > xlim[1]
  }
  
  fit_line <- x$Fit + x$Baseline
  max_dp <- max(x$Data[ind],fit_line[ind])
  min_dp <- min(x$Data[ind],fit_line[ind],x$Baseline[ind])
  
  res <- x$Data - fit_line
  res_range <- max(res[ind]) - min(res[ind])
  offset <- max_dp - min(res[ind])
  
  basis_yoff <- (max_dp - min_dp) * y_offset
  
  min_basis <- Inf
  for (p in 5:ncol(x)) {
    x[,p] <- x[,p] - (p - 4) * basis_yoff + min_dp
    if (min(x[ind, p]) < min_basis) min_basis <- min(x[ind, p])
  }
  
  graphics::plot(x$PPMScale, x$Data, type = 'l', xlim = xlim, 
       ylim = c(min_basis - basis_yoff, max_dp + res_range), yaxt = "n", ylab = "",
       xlab = "Chemical Shift (ppm)", ...)
  graphics::lines(x$PPMScale, fit_line, col = 'Red', lw = 2)
  graphics::lines(x$PPMScale, x$Baseline)
  graphics::lines(x$PPMScale, res + offset)
  graphics::abline(h = max_dp)
  
  for (p in 5:ncol(x)) {
    graphics::lines(x$PPMScale, x[,p])
  }
}

#' Print a summary of an object of class \code{fit_result}.
#' @param x \code{fit_result} object.
#' @param ... further arguments.
#' @export
print.fit_result <- function(x, ...) {
  print(x$res_tab)
}

output_csv <- function(analysis, fname, pvc=FALSE) {
  if (pvc == TRUE) {
    utils::write.csv(analysis$results_pvc, fname, quote = FALSE, 
                     row.names = FALSE)
  } else {
    utils::write.csv(analysis$results, fname, quote = FALSE, row.names = FALSE)
  }
}

plot_slice <- function(analysis, name) {
  result_map <- analysis$results[[name]]
  dim(result_map) <- dim(analysis$data$data)[2:6]
  graphics::image(result_map[,, 1, 1, 1])
}