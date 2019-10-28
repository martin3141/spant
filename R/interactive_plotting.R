#' Plot a 2D slice from an MRSI fit result object.
#' @param fit_res \code{fit_result} object.
#' @param map array of values to be plotted, defaults to a "TNAA" map.
#' @param slice slice to plot in the z direction.
#' @param zlim range of values to plot.
#' @param interp interpolation factor.
#' @export
plot_slice_fit_inter <- function(fit_res, map = NULL, slice = 1, zlim = NULL, 
                                 interp = 1) {
  
  if (is.null(map)) map <- get_fit_map(fit_res, "TNAA") 
  
  plot_slice_map_inter(mrs_data = fit_res, map = map, slice = slice, 
                       interp = interp, zlim = zlim)
}

#' Plot an interactive slice map from a data array where voxels can be selected
#' to display a corresponding spectrum.
#' @param mrs_data spectral data.
#' @param map array of values to be plotted, defaults to the integration of the
#' modulus of the full spectral width.
#' @param xlim spectral region to plot.
#' @param slice the slice index to plot.
#' @param zlim smallest and largest values to be plotted.
#' @param mask_map matching map with logical values to indicate if the
#' corresponding values should be plotted.
#' @param denom map to use as a denominator.
#' @param mask_cutoff minimum values to plot (as a percentage of the maximum).
#' @param interp map interpolation factor.
#' @param mode representation of the complex spectrum to be plotted, can be one
#' of: "re", "im", "mod" or "arg".
#' @param y_scale option to display the y-axis values (logical).
#' @param ylim intensity range to plot.
#' @param coil coil element to plot.
#' @export
#' @importFrom tkrplot tkrplot
plot_slice_map_inter <- function(mrs_data, map = NULL, xlim = NULL, slice = 1,
                                 zlim = NULL, mask_map = NULL, denom = NULL, 
                                 mask_cutoff = 20, interp = 1, mode = "re",
                                 y_scale = FALSE, ylim = NULL, coil = 1) {
  
  # kill any existing plots
  if (exists("plot_env")) {
    if (exists("win1",plot_env)) {
      tcltk::tkdestroy(plot_env$win1)
    }
  }

  assign("plot_env", new.env(hash = TRUE), envir = baseenv())
  #assign("plot_env", new.env(hash = TRUE), envir = globalenv())
  
  if (is.null(map)) map <- int_spec(mrs_data, mode = "mod")

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
  plot_env$coil <- coil
  plot_env$mask_map <- mask_map
  plot_env$zlim <- zlim
  plot_env$denom <- denom
  plot_env$mask_cutoff <- mask_cutoff
  plot_env$interp <- interp
  plot_env$mode <- mode
  plot_env$y_scale <- y_scale
  plot_env$ylim <- ylim

  #map_data <- map[1,,,1,1,1]
  #map_data <- pracma::fliplr(map_data)
  #map_data <- pracma::fliplr(map_data)
  mrs_data <- mrs_data

  plot_env$map_data <- map
  plot_env$mrs_data <- mrs_data

  plot_env$x <- 1
  plot_env$y <- 1

  plot_env$win1 <- tcltk::tktoplevel(class = "spant_plot")

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

  graphics::par(mar = c(1,1,1,6))
  plot_slice_map(plot_env$map_data, slice = plot_env$slice,
                 mask_map = plot_env$mask_map, zlim = plot_env$zlim,
                 denom = plot_env$denom, mask_cutoff = plot_env$mask_cutoff, 
                 interp = plot_env$interp,
                 horizontal = FALSE, coil = plot_env$coil)

  graphics::points((plot_env$xPlotCoord), (plot_env$yPlotCoord), col = "white",
                    cex = 4, lw = 3)

  plot_env$parPlotSize    <- graphics::par("plt")
  plot_env$parPlotSize[1] <- plot_env$parPlotSize[1] / 2 # correction for subplt
  plot_env$parPlotSize[2] <- plot_env$parPlotSize[2] / 2 # correction for subplt
  plot_env$usrCoords      <- graphics::par("usr")
  val <- plot_env$map_data[1, plot_env$x, plot_env$y, plot_env$slice, 1, 1]
  text = paste("X=", plot_env$x, "\t Y=", plot_env$y, "\t val=", val, sep = "")
  cat(text, "\n")
  
  if (class(plot_env$mrs_data) == "mrs_data") {
    graphics::plot(plot_env$mrs_data, x_pos = plot_env$x, y_pos = plot_env$y,
                   z_pos = plot_env$slice, xlim = plot_env$xlim,
                   mode = plot_env$mode, y_scale = plot_env$y_scale,
                   ylim = plot_env$ylim, coil = plot_env$coil)
  } else {
    graphics::plot(plot_env$mrs_data, x_pos = plot_env$x, y_pos = plot_env$y,
                   z_pos = plot_env$slice, xlim = plot_env$xlim,
                   coil = plot_env$coil)
  }
}

#' Diaply an orthographic projection plot of a nifti object.
#' @param underlay underlay image to be shown in grayscale.
#' @param overlay optional overlay image.
#' @param xyz x, y, z slice coordinates to display.
#' @param zlim underlay intensity limits.
#' @param zlim_ol overlay intensity limits.
#' @param alpha transparency of overlay.
#' @param col_ol color palette of overlay.
#' @param orient_lab display orientation labels (default TRUE).
#' @param rescale rescale factor for the underlay and overlay images.
#' @param crosshairs display the crosshairs (default TRUE).
#' @export
ortho3 <- function(underlay, overlay = NULL, xyz = NULL, zlim = NULL,
                   zlim_ol = NULL, alpha = 1, col_ol = viridisLite::viridis(64),
                   orient_lab = TRUE, rescale = 1, crosshairs = TRUE) {
  
  if ((RNifti::orientation(underlay) != "RAS") && (orient_lab)) {
    warning("Underlay image is not in RAS format, orientation labels may be incorrect.")
  }
  
  graphics::par(bg = "black", mar = c(0,0,0,0))

  img_dim <- dim(underlay)[1:3]

  if (is.null(xyz)) xyz <- ceiling(img_dim / 2)

  cor <- underlay[img_dim[1]:1, xyz[2],]
  sag <- underlay[xyz[1], img_dim[2]:1,]
  ax <-  underlay[img_dim[1]:1,, xyz[3]]

  xcutoff <- img_dim[1] / (img_dim[1] + img_dim[2])
  ycutoff <- img_dim[2] / (img_dim[2] + img_dim[3])

  dummy <- matrix(0, img_dim[2], img_dim[2])

  full <- cbind(rbind(ax, dummy), rbind(cor, sag))
  
  if (rescale != 1) {
    full <- mmand::rescale(full, rescale, mmand::triangleKernel())
  }

  if (is.null(zlim)) {
    zlim <- range(underlay)
  } else {
    full[full < zlim[1]] <- zlim[1]
    full[full > zlim[2]] <- zlim[2]
  }

  asp <- dim(full)[2] / dim(full)[1]
  graphics::image(full, useRaster = TRUE, col = grDevices::gray(0:64 / 64),
                  axes = FALSE, asp = asp, zlim = zlim)
  
  if (!is.null(overlay)) {
    if ((RNifti::orientation(overlay) != "RAS") && (orient_lab)) {
      warning("Overlay image is not in RAS format, orientation labels may be incorrect.")
    }
    cor_y <- overlay[img_dim[1]:1, xyz[2],]
    sag_y <- overlay[xyz[1], img_dim[2]:1,]
    ax_y <-  overlay[img_dim[1]:1,, xyz[3]]
    dummy_y <- matrix(NA, img_dim[2], img_dim[2])
    full_y <- cbind(rbind(ax_y, dummy_y), rbind(cor_y, sag_y))
    
    if (rescale != 1) {
      full_y <- mmand::rescale(full_y, rescale, mmand::triangleKernel())
    }
    
    if (is.null(zlim_ol)) zlim_ol <- range(full_y, na.rm = TRUE)
      
    full_y[full_y <= zlim_ol[1]] <- NA
    full_y[full_y > zlim_ol[2]] <- zlim_ol[2]
    
    col_ol <- add_alpha(col_ol, alpha)
    
    graphics::image(full_y, useRaster = TRUE, col = col_ol, axes = FALSE, asp = asp,
          add = TRUE, zlim = zlim_ol)
  }

  if (crosshairs) {
    # top right verical
    graphics::lines(rep(xcutoff + (img_dim[2] - xyz[2]) / img_dim[2] * (1 - xcutoff), 2),
          c(ycutoff, 1), col = "red")
    
    # upper left horizonal
    graphics::lines(c(0, 1), rep(ycutoff + xyz[3] / img_dim[3] * (1 - ycutoff), 2),
          col = "red")
    
    # lower left horizontal
    graphics::lines(c(0, xcutoff), rep(xyz[2] / img_dim[2] * ycutoff, 2), col = "red")
    
    # lower left vertical
    graphics::lines(rep((img_dim[1] - xyz[1]) / img_dim[1] * xcutoff, 2), c(0, 1),
          col = "red")
  }


  if (orient_lab) {
    cex_lab <- 0.8
    lab_marg <- 0.5
    lab_font <- 2
    lm_pos <- 1 + lab_marg
    lm_neg <- -lab_marg
    
    graphics::text(0.0, ycutoff / 2, "R", col = "white", cex = cex_lab, adj = lm_neg,
         font = lab_font)
    #text(xcutoff, ycutoff / 2, "L", col = "white", cex = cex_lab, adj = lm_pos,
    #     font = lab_font)
    
    graphics::text(0.0, ycutoff + (1 - ycutoff) / 2, "R", col = "white", cex = cex_lab,
         adj = lm_neg, font = lab_font)
    #text(xcutoff, ycutoff + (1 - ycutoff) / 2, "L", col = "white", cex = cex_lab,
    #     adj = lm_pos, font = lab_font)
    
    graphics::text(xcutoff / 2, 0, "P", col = "white", cex = cex_lab, adj = c(1, lm_neg),
         font = lab_font)
    #text(xcutoff / 2, ycutoff, "A", col = "white", cex = cex_lab,
    #     adj = c(1, lm_pos), font = lab_font)
    
    graphics::text(xcutoff / 2, 1, "S", col = "white", cex = cex_lab, adj = c(1, lm_pos),
         font = lab_font)
    #text(xcutoff / 2, ycutoff, "I", col = "white", cex = cex_lab,
    #     adj = c(1.7, lm_neg), font = lab_font)
    
    graphics::text(1.0, ycutoff + (1 - ycutoff) / 2, "P", col = "white", cex = cex_lab,
         adj = lm_pos, font = lab_font)
    #text(xcutoff, ycutoff + (1 - ycutoff) / 2, "A", col = "white", cex = cex_lab,
    #     adj = lm_neg, font = lab_font)
    
    graphics::text(xcutoff + (1 - xcutoff) / 2, 1, "S", col = "white", cex = cex_lab,
                   adj = c(1, lm_pos), font = lab_font)
    # text(xcutoff + (1 - xcutoff) / 2, ycutoff, "I", col = "white", cex = cex_lab,
    #     adj = c(1.9, lm_neg), font = lab_font)
  }
}

#' Diaply an interactive orthographic projection plot of a nifti object.
#' @param underlay underlay image to be shown in grayscale.
#' @param overlay optional overlay image.
#' @param xyz x, y, z slice coordinates to display.
#' @param zlim underlay intensity limits.
#' @param zlim_ol overlay intensity limits.
#' @param alpha transparency of overlay.
#' @param ... other options to be passed to the ortho3 function.
#' @export
ortho3_int <- function(underlay, overlay = NULL, xyz = NULL, zlim = NULL,
                       zlim_ol = NULL, alpha = 1, ...) {
  
  img_dim <- dim(underlay)[1:3]
  if (is.null(xyz)) xyz <- ceiling(img_dim / 2)
  
  mri_range <- signif(range(underlay), 3)
  if (is.null(zlim)) zlim <- mri_range
  
  if (is.null(overlay)) {
    mri_range_y <- c(0, 1)
    zlim_ol <- c(0, 1)
  } else {
    mri_range_y <- signif(range(overlay), 3)
    if (is.null(zlim_ol)) zlim_ol <- mri_range_y
  }
  
  xstep <- (mri_range[2] - mri_range[1]) / 100
  ystep <- (mri_range_y[2] - mri_range_y[1]) / 100
  
  ui <- miniUI::miniPage(
    
    miniUI::miniContentPanel(
      shiny::fillRow(flex = c(1, NA),
              shiny::plotOutput("plot", height = "100%", click = "plot_click"),
              shiny::fillCol(width = "90%", flex = rep(NA, 7),
                    shiny::sliderInput("zlim_sli", "Range", mri_range[1], mri_range[2], zlim, width = 200, step = xstep),
                    shiny::numericInput("x_in", "x", xyz[1], 1, img_dim[1], width = 80), 
                    shiny::numericInput("y_in", "y", xyz[2], 1, img_dim[2], width = 80),
                    shiny::numericInput("z_in", "z", xyz[3], 1, img_dim[3], width = 80),
                    shiny::sliderInput("zlim_y_sli", "Overlay range", mri_range_y[1], mri_range_y[2], zlim_ol, width = 200, step = ystep),
                    shiny::sliderInput("alpha_sli", "Overlay alpha", 0, 1, alpha, width = 200),
                    shiny::actionButton("done", "Done", icon = shiny::icon("check"), width = 150))
              ), padding = 0))
  
  server <- function(input, output, session) {
    # Render the plot
    output$plot <- shiny::renderPlot({
      ortho3(underlay, overlay, xyz, zlim, zlim_ol, alpha, ...)
    })
    
    shiny::observeEvent(input$zlim_sli, {zlim <<- input$zlim_sli
                 output$plot <- shiny::renderPlot(ortho3(underlay, overlay, xyz, zlim, zlim_ol, alpha, ...))}
                 )
    
    shiny::observeEvent(input$zlim_y_sli, {zlim_ol <<- input$zlim_y_sli
                 output$plot <- shiny::renderPlot(ortho3(underlay, overlay, xyz, zlim, zlim_ol, alpha, ...))}
                 )
    
    shiny::observeEvent(input$alpha_sli, {alpha <<- input$alpha_sli
                 output$plot <- shiny::renderPlot(ortho3(underlay, overlay, xyz, zlim, zlim_ol, alpha, ...))}
                 )
      
    shiny::observeEvent(input$x_in, {xyz[1] <<- input$x_in
                 output$plot <- shiny::renderPlot(ortho3(underlay, overlay, xyz, zlim, zlim_ol, alpha, ...))}
                 )
    
    shiny::observeEvent(input$y_in, {xyz[2] <<- input$y_in
                 output$plot <- shiny::renderPlot(ortho3(underlay, overlay, xyz, zlim, zlim_ol, alpha, ...))}
                 )
    
    shiny::observeEvent(input$z_in, {xyz[3] <<- input$z_in
                 output$plot <- shiny::renderPlot(ortho3(underlay, overlay, xyz, zlim, zlim_ol, alpha, ...))}
                 )
      
    shiny::observeEvent(input$plot_click, {
              
                 xpos <- input$plot_click$x
                 ypos <- input$plot_click$y
                 img_dim <- dim(underlay)
    
                 if (xpos > 1) xpos <- 1
                 if (xpos < 0) xpos <- 0
                 if (ypos > 1) ypos <- 1
                 if (ypos < 0) ypos <- 0
    
                 xcutoff <- img_dim[1] / (img_dim[1] + img_dim[2])
                 ycutoff <- img_dim[2] / (img_dim[2] + img_dim[3])
    
                 if (xpos < xcutoff) {
                    left <- TRUE
                    xpos <- xpos / xcutoff
                 } else {
                    left <- FALSE
                    xpos <- (xpos - xcutoff) / (1 - xcutoff)
                 }
    
                 if (ypos > ycutoff) {
                    top <- TRUE
                    ypos <- (ypos - ycutoff) / (1 - ycutoff)
                 } else {
                    top <- FALSE
                    ypos <- ypos / ycutoff
                 }
    
                 if (left && top) {
                    xyz[1] <<- round((1-xpos) * img_dim[1])
                    xyz[3] <<- round(ypos * img_dim[3])
                 } else if (left && !top) {
                    xyz[1] <<- round((1-xpos) * img_dim[1])
                    xyz[2] <<- round(ypos * img_dim[2])
                 } else if (!left && top) {
                    xyz[2] <<- round((1-xpos) * img_dim[2])
                    xyz[3] <<- round(ypos * img_dim[3])
                 }
    
                 if (xyz[1] == 0) xyz[1] <<- 1
                 if (xyz[2] == 0) xyz[2] <<- 1
                 if (xyz[3] == 0) xyz[3] <<- 1
                     
                 shiny::updateNumericInput(session, "x_in", value = xyz[1])
                 shiny::updateNumericInput(session, "y_in", value = xyz[2])
                 shiny::updateNumericInput(session, "z_in", value = xyz[3])
                 })
    
    # Handle the Done button being pressed.
    shiny::observeEvent(input$done, shiny::stopApp())
  }
  
  shiny::runGadget(ui, server, viewer = shiny::paneViewer())
}
