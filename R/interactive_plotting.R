#' Plot a 2D slice from an MRSI fit result object
#' @param fit_res \code{fit_result} object
#' @param name name of the quantity to plot, eg "TNAA"
#' @param slice slice to plot in the z direction
#' @param zlim range of values to plot
#' @param interp interpolation factor
#' @export
plot_fit_slice_inter <- function(fit_res, name, slice = 1, zlim = NULL, 
                                 interp = 1) {
  
  map <- get_fit_map(fit_res, name) 
  plot_slice_map_inter(map = map, mrs_data = fit_res, slice = slice, 
                       interp = interp, zlim = zlim)
}

#' Plot an interactive slice map from a data array where voxels can be selected
#' to display a corresponding spectrum
#' @param map array of values to be plotted
#' @param mrs_data spectral data
#' @param xlim spectral region to plot
#' @param slice the slice index to plot
#' @param zlim smallest and largest values to be plotted
#' @param mask_map matching map with logical values to indicate if the
#' corresponding values should be plotted
#' @param denom map to use as a denominator
#' @param mask_cutoff minimum values to plot (as a percentage of the maximum)
#' @param interp map interpolation factor
#' @export
#' @importFrom tkrplot tkrplot
plot_slice_map_inter <- function(map, mrs_data, xlim = NULL, slice = 1,
                                 zlim = NULL, mask_map = NULL, denom = NULL, 
                                 mask_cutoff = 20, interp = 1) {
  
  # kill any existing plots
  if (exists("plot_env")) {
    if (exists("win1",plot_env)) {
      tcltk::tkdestroy(plot_env$win1)
    }
  }

  assign("plot_env", new.env(hash = TRUE), envir = baseenv())
  #assign("plot_env", new.env(hash = TRUE), envir = globalenv())

  if (class(mrs_data) == "mrs_data") {
    x_scale <- ppm(mrs_data)
  } else {
    x_scale <- mrs_data$fits[[1]]$PPMScale
  }

  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[length(x_scale)])
  }

  plot_env$xlim <- xlim
  plot_env$slice <- slice
  plot_env$mask_map <- mask_map
  plot_env$zlim <- zlim
  plot_env$denom <- denom
  plot_env$mask_cutoff <- mask_cutoff
  plot_env$interp <- interp

  #map_data <- map[1,,,1,1,1]
  #map_data <- pracma::fliplr(map_data)
  #map_data <- pracma::fliplr(map_data)
  mrs_data <- mrs_data

  plot_env$map_data <- map
  plot_env$mrs_data <- mrs_data

  plot_env$x <- 1
  plot_env$y <- 1


  plot_env$win1 <- tcltk::tktoplevel()

  #plot_env$win1$env$plot <- tkrplot::tkrplot(plot_env$win1, fun = plotTk,
  #                                           hscale = 3.0, vscale = 1.5)


  plot_env$win1$env$plot <- tkrplot::tkrplot(plot_env$win1, fun = plotTk,
                                             hscale = 2.5, vscale = 1.25)
  tcltk::tkgrid(plot_env$win1$env$plot)

  tcltk::tkbind(plot_env$win1$env$plot, "<Button-1>", onLeftClick)
  tcltk::tkconfigure(plot_env$win1$env$plot, cursor = "hand2")
}

onLeftClick <- function(x, y) {
  xClick <- x
  yClick <- y
  width  <- as.numeric(tcltk::tclvalue(tcltk::tkwinfo("reqwidth",
                                                      plot_env$win1$env$plot)))

  height <- as.numeric(tcltk::tclvalue(tcltk::tkwinfo("reqheight",
                                                      plot_env$win1$env$plot)))

  xMin <- plot_env$parPlotSize[1] * width
  xMax <- plot_env$parPlotSize[2] * width
  yMin <- plot_env$parPlotSize[3] * height
  yMax <- plot_env$parPlotSize[4] * height

  rangeX <- plot_env$usrCoords[2] - plot_env$usrCoords[1]
  rangeY <- plot_env$usrCoords[4] - plot_env$usrCoords[3]

  #imgXcoords <- (xCoords - usrCoords[1]) * (xMax - xMin) / rangeX + xMin
  #imgYcoords <- (yCoords - usrCoords[3]) * (yMax - yMin) / rangeY + yMin

  xClick <- as.numeric(xClick) + 0.5
  yClick <- as.numeric(yClick) + 0.5
  yClick <- height - yClick

  xPlotCoord <- plot_env$usrCoords[1] + (xClick - xMin) * rangeX / (xMax - xMin)
  yPlotCoord <- plot_env$usrCoords[3] + (yClick - yMin) * rangeY / (yMax - yMin)

  #plot_env$xPlotCoord <- xPlotCoord
  #plot_env$yPlotCoord <- yPlotCoord

  #print(xPlotCoord)
  #print(yPlotCoord)
  x_len <- dim(plot_env$mrs_data)[2]
  y_len <- dim(plot_env$mrs_data)[3]
  plot_env$x <- round(xPlotCoord * (x_len - 1)) + 1
  plot_env$y <- y_len - round(yPlotCoord*(y_len - 1))

  plot_env$xPlotCoord <- (plot_env$x - 1)  / (x_len - 1)
  plot_env$yPlotCoord <- -(plot_env$y - y_len) / (y_len - 1)

  tkrplot::tkrreplot(plot_env$win1$env$plot)

  #msg <- paste0("Label the point closest to these ",
  #              "approximate plot coordinates: \n",
  #              "x = ", format(xPlotCoord, digits = 2),
  #              ", y = ", format(yPlotCoord, digits = 2), "?")
  #mbval <- tkmessageBox(title =
  #                 "Label Point Closest to These Approximate Plot Coordinates",
  #                      message = msg, type = "yesno", icon = "question")

  #if (tclvalue(mbval)== "yes")
  #  labelClosestPoint(xClick, yClick, imgXcoords, imgYcoords)
}

plotTk <- function() {
  graphics::par(mfrow = c(1,2))

  #image(plot_env$map_data, col = viridis::viridis(128), useRaster = T,
  #      asp = 1, axes = FALSE)

  graphics::par(mar = c(1,1,1,3))
  plot_slice_map(plot_env$map_data, slice = plot_env$slice,
                 mask_map = plot_env$mask_map, zlim = plot_env$zlim,
                 denom = plot_env$denom, mask_cutoff = plot_env$mask_cutoff, 
                 interp = plot_env$interp,
                 horizontal = FALSE)

  graphics::points((plot_env$xPlotCoord), (plot_env$yPlotCoord), col = "white",
                   cex = 4, lw = 3)

  plot_env$parPlotSize    <- graphics::par("plt")
  plot_env$parPlotSize[1] <- plot_env$parPlotSize[1] / 2 # correction for subplt
  plot_env$parPlotSize[2] <- plot_env$parPlotSize[2] / 2 # correction for subplt
  plot_env$usrCoords      <- graphics::par("usr")
  text = paste("X=", plot_env$x, ", Y=", plot_env$y, sep = "")
  cat(text, "\n")
  graphics::plot(plot_env$mrs_data, x_pos = plot_env$x, y_pos = plot_env$y,
                 z_pos = plot_env$slice, xlim = plot_env$xlim)
                 #z_pos = plot_env$slice, yscale = TRUE, xlim = plot_env$xlim)
}