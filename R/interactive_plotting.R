#' Plot a 2D slice from an MRSI fit result object.
#' @param fit_res \code{fit_result} object.
#' @param map array of values to be plotted, defaults to a "TNAA" map.
#' @param slice slice to plot in the z direction.
#' @param zlim range of values to plot.
#' @param interp interpolation factor.
#' @param xlim spectral plot limits for the x axis.
#' @export
plot_slice_fit_inter <- function(fit_res, map = NULL, slice = 1, zlim = NULL, 
                                 interp = 1, xlim = NULL) {
  
  if (is.null(map)) map <- get_fit_map(fit_res, "TNAA") 
  
  plot_slice_map_inter(mrs_data = fit_res, map = map, slice = slice, 
                       interp = interp, zlim = zlim, xlim = xlim)
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
plot_slice_map_inter <- function(mrs_data, map = NULL, xlim = NULL, slice = 1,
                                  zlim = NULL, mask_map = NULL, denom = NULL,
                                  mask_cutoff = 20, interp = 1, mode = "re",
                                  y_scale = FALSE, ylim = NULL, coil = 1) {
  
  if (is.null(map)) map <- int_spec(mrs_data, mode = "mod")
  
  if (class(mrs_data) == "mrs_data") {
    x_scale <- ppm(mrs_data)
  } else {
    non_na_res <- which(!is.na(mrs_data$fits))[[1]]
    x_scale <- mrs_data$fits[[non_na_res]]$PPMScale
  }
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[length(x_scale)])
  }
  
  ui <- miniUI::miniPage(
    miniUI::gadgetTitleBar("Select point on the map to show spectrum."),
    miniUI::miniContentPanel(
      shiny::fillRow(
        shiny::plotOutput("map", height = "100%", click = "plot_click"),
        shiny::plotOutput("spec", height = "100%")
      )
    )
  )
  
  server <- function(input, output, session) {
    x_min <- - 1 / (Nx(mrs_data) * interp - 1) / 2
    x_max <- 1 + 1 / (Nx(mrs_data) * interp - 1) / 2
    y_min <- - 1 / (Ny(mrs_data) * interp - 1) / 2
    y_max <- 1 + 1 / (Ny(mrs_data) * interp - 1) / 2
    
    x <- round(Nx(mrs_data) / 2)
    y <- round(Ny(mrs_data) / 2)
    xpos_round <- (x - 0.5) * (x_max - x_min) / Nx(mrs_data) + x_min
    ypos_round <- (y - 0.5) * (y_max - y_min) / Ny(mrs_data) + y_min
    
    output$map <- shiny::renderPlot({
      plot_slice_map(map, slice = slice, mask_map = mask_map, zlim = zlim,
                     denom = denom, mask_cutoff = mask_cutoff, interp = interp,
                     horizontal = FALSE, coil = coil)
      
      graphics::points(xpos_round, ypos_round, pch = 1, col = "red", cex = 4,
                       lw = 3)
    })
    output$spec <- shiny::renderPlot({
      if (class(mrs_data) == "mrs_data") {
        graphics::plot(mrs_data, x_pos = x, y_pos = Ny(mrs_data) + 1 - y,
                       z_pos = slice, xlim = xlim, mode = mode,
                       y_scale = y_scale, ylim = ylim, coil = coil)
      } else {
        graphics::plot(mrs_data, x_pos = x, y_pos = Ny(mrs_data) + 1 - y,
                       z_pos = slice, xlim = xlim, coil = coil)
      }
      
    })
    
    shiny::observeEvent(input$plot_click, {
      xpos <- input$plot_click$x
      ypos <- input$plot_click$y
      
      # scale coordinates to be between 0.5 and Nx + 0.5, Ny + 0.5
      x_rescale <- (xpos - x_min) / (x_max - x_min) * Nx(mrs_data) + 0.5
      y_rescale <- (ypos - y_min) / (y_max - y_min) * Ny(mrs_data) + 0.5
      
      # round position to the nearest voxel coord
      x <- round(x_rescale)
      y <- round(y_rescale)
      
      cat("x = ", x, ", y = ", Ny(mrs_data) + 1 - y, "\n", sep = "")
      
      if (x > Nx(mrs_data)) x <- Nx(mrs_data)
      if (y > Ny(mrs_data)) y <- Ny(mrs_data)
      if (x < 1) x <- 1
      if (y < 1) y <- 1
      
      output$spec <- shiny::renderPlot({
        if (class(mrs_data) == "mrs_data") {
          graphics::plot(mrs_data, x_pos = x, y_pos = Ny(mrs_data) + 1 - y,
                         z_pos = slice, xlim = xlim, mode = mode,
                         y_scale = y_scale, ylim = ylim, coil = coil)
        } else {
          graphics::plot(mrs_data, x_pos = x, y_pos = Ny(mrs_data) + 1 - y,
                         z_pos = slice, xlim = xlim, coil = coil)
        }
      })
      
      output$map <- shiny::renderPlot({
        plot_slice_map(map, slice = slice, mask_map = mask_map, zlim = zlim,
                       denom = denom, mask_cutoff = mask_cutoff, 
                       interp = interp, horizontal = FALSE, coil = coil)
        xpos_round <- (x - 0.5) * (x_max - x_min) / Nx(mrs_data) + x_min
        ypos_round <- (y - 0.5) * (y_max - y_min) / Ny(mrs_data) + y_min
        graphics::points(xpos_round, ypos_round, pch = 1, col = "red", 
                         cex = 4, lw = 3)})
    })
    
    # When the Done button is clicked, return a value
    shiny::observeEvent(input$done, {
      shiny::stopApp()
    })
  }
  
  shiny::runGadget(ui, server, viewer = shiny::dialogViewer("Plot slice map",
                                                            width = 1200,
                                                            height = 600))
}

#' Display an orthographic projection plot of a nifti object.
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
#' @param colourbar display a colourbar for the overlay (default TRUE).
#' @export
ortho3 <- function(underlay, overlay = NULL, xyz = NULL, zlim = NULL,
                   zlim_ol = NULL, alpha = 1, col_ol = viridisLite::viridis(64),
                   orient_lab = TRUE, rescale = 1, crosshairs = TRUE,
                   colourbar = TRUE) {
  
  if ((RNifti::orientation(underlay) != "RAS") && (orient_lab)) {
    warning("Underlay image is not in RAS format, orientation labels may be incorrect.")
  }
  
  graphics::par(bg = "black", fg = "white", col.axis = "white",
                mar = c(0,0,0,0))

  img_dim <- dim(underlay)[1:3]

  if (is.null(xyz)) xyz <- ceiling(img_dim / 2)

  cor <- underlay[img_dim[1]:1, xyz[2],]
  sag <- underlay[xyz[1], img_dim[2]:1,]
  ax  <- underlay[img_dim[1]:1,, xyz[3]]

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
    
    if (colourbar) {
      fields::image.plot(full_y, useRaster = TRUE, col = col_ol, axes = FALSE,
                         asp = asp, add = TRUE, zlim = zlim_ol,
                         smallplot = c(0.50, 0.51, 0.1, 0.4))
    } else {
      graphics::image(full_y, useRaster = TRUE, col = col_ol, axes = FALSE,
                      asp = asp, add = TRUE, zlim = zlim_ol)
    }
  }

  if (crosshairs) {
    # top right verical
    graphics::lines(rep(xcutoff + (img_dim[2] - xyz[2]) / img_dim[2] *
                       (1 - xcutoff), 2), c(ycutoff, 1), col = "red")
    
    # upper left horizonal
    graphics::lines(c(0, 1), rep(ycutoff + xyz[3] / img_dim[3] *
                                (1 - ycutoff), 2), col = "red")
    
    # lower left horizontal
    graphics::lines(c(0, xcutoff), rep(xyz[2] / img_dim[2] * ycutoff, 2),
                    col = "red")
    
    # lower left vertical
    graphics::lines(rep((img_dim[1] - xyz[1]) / img_dim[1] * xcutoff, 2),
                    c(0, 1), col = "red")
  }

  if (orient_lab) {
    cex_lab <- 0.8
    lab_marg <- 0.5
    lab_font <- 2
    lm_pos <- 1 + lab_marg
    lm_neg <- -lab_marg
    
    graphics::text(0.0, ycutoff / 2, "R", col = "white", cex = cex_lab,
                   adj = lm_neg, font = lab_font)
    
    #text(xcutoff, ycutoff / 2, "L", col = "white", cex = cex_lab, adj = lm_pos,
    #     font = lab_font)
    
    graphics::text(0.0, ycutoff + (1 - ycutoff) / 2, "R", col = "white",
                   cex = cex_lab, adj = lm_neg, font = lab_font)
    
    #text(xcutoff, ycutoff + (1 - ycutoff) / 2, "L", col = "white", cex = cex_lab,
    #     adj = lm_pos, font = lab_font)
    
    graphics::text(xcutoff / 2, 0, "P", col = "white", cex = cex_lab,
                   adj = c(1, lm_neg), font = lab_font)
    
    #text(xcutoff / 2, ycutoff, "A", col = "white", cex = cex_lab,
    #     adj = c(1, lm_pos), font = lab_font)
    
    graphics::text(xcutoff / 2, 1, "S", col = "white", cex = cex_lab,
                   adj = c(1, lm_pos), font = lab_font)
    
    #text(xcutoff / 2, ycutoff, "I", col = "white", cex = cex_lab,
    #     adj = c(1.7, lm_neg), font = lab_font)
    
    graphics::text(1.0, ycutoff + (1 - ycutoff) / 2, "P", col = "white",
                   cex = cex_lab, adj = lm_pos, font = lab_font)
    
    #text(xcutoff, ycutoff + (1 - ycutoff) / 2, "A", col = "white", cex = cex_lab,
    #     adj = lm_neg, font = lab_font)
    
    graphics::text(xcutoff + (1 - xcutoff) / 2, 1, "S", col = "white",
                   cex = cex_lab, adj = c(1, lm_pos), font = lab_font)
    
    # text(xcutoff + (1 - xcutoff) / 2, ycutoff, "I", col = "white", cex = cex_lab,
    #     adj = c(1.9, lm_neg), font = lab_font)
  }
}

#' Display an interactive orthographic projection plot of a nifti object.
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
  
  mri_range <- signif(range(underlay, na.rm = TRUE), 3)
  if (is.null(zlim)) zlim <- mri_range
  
  if (is.null(overlay)) {
    mri_range_y <- c(0, 1)
    zlim_ol <- c(0, 1)
  } else {
    mri_range_y <- signif(range(overlay, na.rm = TRUE), 3)
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
