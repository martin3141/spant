#' Normalise mrs_data to a spectral region.
#' @param mrs_data MRS data.
#' @param xlim spectral range to be integrated (defaults to full range).
#' @param scale units of xlim, can be : "ppm", "Hz" or "points".
#' @param mode spectral mode, can be : "re", "im", "mod" or "cplx".
#' @param summation can be "sum", "mean" or "l2" (default).
#' @return normalised data.
#' @export
norm_mrs <- function(mrs_data, xlim = NULL, scale = "ppm", mode = "re",
                     summation = "l2") {
  
  warning("norm_mrs is deprecated, use scale_mrs instead")
  
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  amps <- int_spec(mrs_data, xlim, scale, mode, summation)
  amps_full <- array(rep(amps, Npts(mrs_data)), dim = dim(mrs_data$data))
  mrs_data$data <- mrs_data$data / amps_full
  return(mrs_data)
}

#' Integrate a spectral region.
#' @param mrs_data MRS data.
#' @param xlim spectral range to be integrated (defaults to full range).
#' @param scale units of xlim, can be : "ppm", "hz" or "points".
#' @param mode spectral mode, can be : "re", "im", "mod" or "cplx".
#' @param summation can be "sum" (default), "mean" or "l2".
#' @return an array of integral values.
#' @export
int_spec <- function(mrs_data, xlim = NULL, scale = "ppm", mode = "re",
                     summation = "sum") {
  
  warning("int_spec is deprecated, use spec_op instead")
  
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
    
  if ( scale == "ppm" ) {
    x_scale <- ppm(mrs_data)
  } else if (scale == "hz") {
    x_scale <- hz(mrs_data)
  } else if (scale == "points") {
    x_scale <- pts(mrs_data)
  }
  
  if (is.null(xlim)) xlim <- c(x_scale[1], x_scale[Npts(mrs_data)])
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  data_arr <- mrs_data$data[,,,,,, subset, drop = F]
  
  if (mode == "re") {
    data_arr <- Re(data_arr)
  } else if (mode == "im") {
    data_arr <- Im(data_arr)
  } else if (mode == "mod") {
    data_arr <- Mod(data_arr)
  }
 
  if (summation == "l2") {
    data_arr <- data_arr * data_arr
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), sum)
    res <- res ^ 0.5
  } else if (summation == "mean") {
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), mean)
  } else {
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), sum)
  }
  
  return(res) 
}