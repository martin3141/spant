plot_slice_map_inter <- function(map, mrs_data, xlim = NULL, slice = 1, 
                                 mask_map = NULL, upper = NULL, lower = NULL,
                                 denom = NULL, mask_cutoff = 20, interp = 16) {
  
  assign("plot_env", new.env(hash = T), envir = baseenv())
  
  x_scale <- ppm(mrs_data)
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(mrs_data)])
  }
  
  
  plot_env$xlim <- xlim
  plot_env$slice <- slice
  plot_env$mask_map <- mask_map
  plot_env$upper <- upper
  plot_env$lower <- lower
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
  
  plot_env$xPlotCoord <- xPlotCoord  
  plot_env$yPlotCoord <- yPlotCoord  
  
  #print(xPlotCoord)
  #print(yPlotCoord)
  x_len <- dim(plot_env$mrs_data)[2]
  y_len <- dim(plot_env$mrs_data)[3]
  plot_env$x <- round(xPlotCoord * (x_len - 1)) + 1
  plot_env$y <- y_len - round(yPlotCoord*(y_len - 1))
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
  #      asp = 1, axes = F)
  
  plot_slice_map(plot_env$map_data, slice = plot_env$slice, 
                 mask_map = plot_env$mask_map, lower = plot_env$lower,
                 upper = plot_env$upper, denom = plot_env$denom,
                 mask_cutoff = plot_env$mask_cutoff, interp = plot_env$interp,
                 horizontal = T)
  
  graphics::points(plot_env$xPlotCoord, plot_env$yPlotCoord, col = "white",
                   cex = 3, lw = 3)
  
  plot_env$parPlotSize    <- graphics::par("plt")
  plot_env$parPlotSize[1] <- plot_env$parPlotSize[1] / 2 # correction for subplt
  plot_env$parPlotSize[2] <- plot_env$parPlotSize[2] / 2 # correction for subplt
  plot_env$usrCoords      <- graphics::par("usr")
  text = paste("X=", plot_env$x, ", Y=", plot_env$y, sep = "")
  cat(text, "\n")
  graphics::plot(plot_env$mrs_data, x_pos = plot_env$x, y_pos = plot_env$y, 
                 z_pos = plot_env$slice, yscale = TRUE, xlim = plot_env$xlim)
}