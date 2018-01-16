#' Apply a function over specified array axes
#' @param x an array
#' @param axes a vector of axes to apply fun over
#' @param fun function to be applied
#' @param ... optional argments to fun
#' @return array
#' @examples
#' z <- array(1:1000, dim = c(10, 10, 10))
#' a <- apply_axes(z, 3, fft)
#' a[1,1,] == fft(z[1,1,])
#' a <- apply_axes(z, 3, sum)
#' a[1,1,] == sum(z[1,1,])
#' @export
apply_axes <- function(x, axes, fun, ...) {
  ndim <- length(dim(x))
  # apply apply function to each axis
  for (axis in axes) {
    margin <- c(1:ndim)[-axis]
    z <- apply(x, margin, fun, ...)
    if (length(dim(z)) != ndim) {
      # insert the dropped dimension
      dim(z) <- append(dim(z), 1, axis - 1)
      x <- z
    } else {
      # permute back to the original order
      perm_vec <- 2:ndim
      perm_vec <- append(perm_vec, 1, axis - 1)
      print(perm_vec)
      x <- aperm(z, perm_vec)
    }
  }
  x
}

# TODO
# see rep dyn for example
#rep_array_dim <- function(x, dim, n) {
#  orig_dim <- dim(x)
#  new_dim <- orig_dim
#  new_dim[dim] <- new_dim[dim] * n
#  # make the dynamic dimension (5th) the last
#  z <- aperm(x, c(1,2,3,4,6,7,5))
#  # duplicate the data
#  z <- rep(z, n)
#  # set the new dimesnions
#  dim(z) <- new_dim[c(1,2,3,4,6,7,5)]
#  # reorder
#  aperm(z, c(1,2,3,4,7,5,6))
#}


#' Covert a beta value in the time-domain to an equivalent linewidth in Hz:
#' x * exp(-i * t * t * beta)
#' @param beta beta damping value 
#' @return linewidth value in Hz
#' @export
beta2lw <- function(beta) {2 * (-beta * log(0.5)) ^ 0.5 / pi}

#' Covert a linewidth in Hz to an equivalent beta value in the time-domain ie:
#' x * exp(-i * t * t * beta)
#' @param lw linewidth in Hz
#' @return beta damping value
#' @export
lw2beta <- function(lw) {(lw * pi / 2) ^ 2 / (-log(0.5))}

alpha2lw <- function(alpha) {alpha / pi}

lw2alpha <- function(lw) {lw * pi}

#' Perform a fft and ffshift on a vector.
#' @param vec_in vector input.
#' @return output vector.
#' @export
ft_shift <- function(vec_in) {pracma::fftshift(stats::fft(vec_in))}

#' Perform an iffshift and ifft on a vector.
#' @param vec_in vector input.
#' @return output vector.
#' @export
ift_shift <- function(vec_in) {pracma::ifft(pracma::ifftshift(vec_in))}

#' Perform a fft and fftshift on a matrix with each column replaced by its 
#' shifted fft.
#' @param mat_in matrix input.
#' @return output matrix.
#' @export
ft_shift_mat <- function(mat_in) {
  mat_in_ft = stats::mvfft(mat_in) # 33s
  #p <- fftw::planFFT(NROW(mat_in),effort=1)
  #mat_in_ft = apply(mat_in,2,fftw::FFT,p) # 34s
  #mat_in_ft = fftwtools::mvfftw(mat_in) # 35s
  mvfftshift(mat_in_ft)
}

#' Perform a fftshift on a matrix, with each column replaced by its shifted 
#' result.
#' @param x matrix input.
#' @return output matrix.
#' @export
mvfftshift <- function(x) {
  m <- NROW(x)
  p <- ceiling(m/2)
  idx <- c((p + 1):m, 1:p)
  x[idx,,drop=FALSE]
}

hilbert <- function(x) {
  x <- Re(x)
  x <- stats::fft(x)
  n <- length(x)
  h <- vector()
  h[1] <- 1
  h[floor(n / 2) + 1] <- 1
  h[2:round(n / 2 + 0.5)] <- 2
  h[floor(n / 2 + 2):n] <- 0
  y <- h * x
  stats::fft(y, inverse = TRUE) / n
}

re_max <- function(x) {
  max(Re(x))
}

re_mean <- function(x) {
  mean(Re(x))
}

re_sd <- function(x) {
  stats::sd(Re(x))
}

is.installed <- function(mypkg) {
  is.element(mypkg, utils::installed.packages()[,1]) 
}

rotate_vec <- function(vec_in, ax, theta) {
  ct <- cos(theta)
  st <- sin(theta)
  rotate_mat <- matrix(nrow = 3, ncol = 3)
  rotate_mat[1, 1] = ct + ax[1] * ax[1] * (1 - ct)
  rotate_mat[1, 2] = ax[1] * ax[2] * (1 - ct) - ax[3] * st
  rotate_mat[1, 3] = ax[1] * ax[3] * (1 - ct) + ax[2] * st
  rotate_mat[2, 1] = ax[2] * ax[1] * (1 - ct) + ax[3] * st
  rotate_mat[2, 2] = ct + ax[1] * ax[1] * (1 - ct)
  rotate_mat[2, 3] = ax[2] * ax[3] * (1 - ct) - ax[1] * st
  rotate_mat[3, 1] = ax[3] * ax[1] * (1 - ct) - ax[2] * st
  rotate_mat[3, 2] = ax[3] * ax[2] * (1 - ct) + ax[1] * st
  rotate_mat[3, 3] = ct + ax[3] * ax[3] * (1 - ct)
  vec_out <- rotate_mat %*% vec_in
  as.numeric(vec_out / sum((vec_out) ^ 2))
}

# From : http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
measure_lorentz_amp <- function(y, t, start_pnt = 10, end_pnt = 50) {
  # crop to time region and take Mod
  y <- Mod(y[start_pnt:end_pnt])
  t <- t[start_pnt:end_pnt]
  
  a <- sum((t ^ 2) * y) * sum(y * log(y)) -  sum(t * y) * sum(t * y * log(y))
  a <- a / (sum(y) * sum(t ^ 2 * y) - sum(t * y) ^ 2)
  A <- exp(a)
  A
}

#' Perform a polynomial fit, subtract and return the standard deviation of the
#' residuals.
#' @param y an array.
#' @param degree polynomial degree.
#' @return the standard deviation of the fit residuals.
#' @export
calc_sd_poly <- function(y, degree = 1) {
  x <- 1:length(y)
  model <- stats::lm(y ~ poly(x, degree))
  stats::sd(model$residual)
}

# Compute the vector cross product between x and y, and return the components
# indexed by i. Stolen from: 
# http://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function
crossprod_3d <- function(x, y, i = 1:3) {
  # Project inputs into 3D, since the cross product only makes sense in 3D.
  To3D <- function(x) utils::head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)
  
  # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
  # so on).  Index3D() lets us do that using R's convention of 1-based (rather
  # than 0-based) arrays.
  Index3D <- function(i) (i - 1) %% 3 + 1
  
  # The i'th component of the cross product is:
  # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat the indices cyclically.
  return(x[Index3D(i + 1)] * y[Index3D(i + 2)] -
            x[Index3D(i + 2)] * y[Index3D(i + 1)])
}

# stolen from the interweb
add_alpha <- function(col, alpha = 1) {
  if (missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, grDevices::col2rgb) / 255, 2, 
        function(x) grDevices::rgb(x[1], x[2], x[3], alpha = alpha))  
}

crop_range <- function(map, lower, upper) {
  # TODO check upper > lower and both are >0<1
  map_range <- range(map, na.rm = TRUE)
  #upper_lim <- map_range[1]+upper/100*(map_range[2] - map_range[1])
  #lower_lim <- map_range[1]+lower/100*(map_range[2] - map_range[1])
  upper_lim <- upper
  lower_lim <- lower
  map <- ifelse(map > upper_lim, upper_lim,map)  
  ifelse(map < lower_lim, lower_lim, map)  
}

# https://stackoverflow.com/questions/13432863/determine-level-of-nesting-in-r
depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)