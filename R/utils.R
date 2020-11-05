#' Apply a function over specified array axes.
#' @param x an array.
#' @param axes a vector of axes to apply fun over.
#' @param fun function to be applied.
#' @param ... optional arguments to fun.
#' @return array.
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

#' Repeat an array over a given dimension.
#' @param x array.
#' @param rep_dim dimension to extend.
#' @param n number of times to repeat.
#' @return extended array.
#' @export
rep_array_dim <- function(x, rep_dim, n) {
  # add a dimension if needed
  if ((rep_dim == 7) && length(dim(x)) == 6) dim(x) <- c(dim(x), 1)
  
  dims <- length(dim(x))
  orig_dim <- dim(x)
  new_dim <- orig_dim
  new_dim[rep_dim] <- new_dim[rep_dim] * n
  # make rep_dim the last dimension
  perm_vec <- (1:dims)[-rep_dim]
  perm_vec <- c(perm_vec, rep_dim)
  z <- aperm(x, perm_vec)
  # duplicate the data
  z <- rep(z, n)
  # set the new dimensions
  dim(z) <- new_dim[perm_vec]
  # reorder
  perm_vec <- append((1:(dims - 1)), dims, rep_dim - 1)
  aperm(z, perm_vec)
}

#rep_array_dim <- function(x, rep_dim, n) {
#  if (dim(x)[rep_dim] != 1) stop("Starting dimension extent does not equal one.")
#  z <- replicate(n, x)
#  dims <- length(dim(z))
#  perm_vec <- 1:dims
#  perm_vec[rep_dim] <- dims
#  perm_vec[dims] <- rep_dim
#  z <- aperm(z, perm_vec)
#  # drop the last dimension
#  abind::adrop(z, dims)
#}

#' Covert a beta value in the time-domain to an equivalent linewidth in Hz:
#' x * exp(-i * t * t * beta).
#' @param beta beta damping value.
#' @return linewidth value in Hz.
#' @export
beta2lw <- function(beta) {2 * (-beta * log(0.5)) ^ 0.5 / pi}

#' Covert a linewidth in Hz to an equivalent beta value in the time-domain ie:
#' x * exp(-t * t * beta).
#' @param lw linewidth in Hz.
#' @return beta damping value.
#' @export
lw2beta <- function(lw) {(lw * pi / 2) ^ 2 / (-log(0.5))}

alpha2lw <- function(alpha) {alpha / pi}

#' Covert a linewidth in Hz to an equivalent alpha value in the time-domain ie:
#' x * exp(-t * alpha).
#' @param lw linewidth in Hz.
#' @return beta damping value.
#' @export
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
  mat_in_ft <- stats::mvfft(mat_in)
  mvfftshift(mat_in_ft)
}

#' Perform an ifft and ifftshift on a matrix with each column replaced by its 
#' shifted ifft.
#' @param mat_in matrix input.
#' @return output matrix.
#' @export
ift_shift_mat <- function(mat_in) {
  mat_in_shift <- mvifftshift(mat_in)
  stats::mvfft(mat_in_shift, inverse = TRUE) / NROW(mat_in)
}

#' Perform a fftshift on a matrix, with each column replaced by its shifted 
#' result.
#' @param x matrix input.
#' @return output matrix.
#' @export
mvfftshift <- function(x) {
  m <- NROW(x)
  p <- ceiling(m / 2)
  idx <- c((p + 1):m, 1:p)
  x[idx,,drop = FALSE]
}

#' Perform an ifftshift on a matrix, with each column replaced by its shifted 
#' result.
#' @param x matrix input.
#' @return output matrix.
#' @export
mvifftshift <- function(x) {
  m <- NROW(x)
  p <- floor(m / 2)
  idx <- c((p + 1):m, 1:p)
  x[idx,,drop = FALSE]
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

re_max_interp <- function(data_pts, interp_f) {
  data_pts <- stats::spline(Re(data_pts), n = interp_f * length(data_pts))
  max(data_pts$y)
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

measure_td_amp <- function(y, start_pnt = 10, end_pnt = 50) {
  
  # check for masked data
  if (is.na(y[1])) return(NA)
  
  # crop to time region and take Mod
  y <- Mod(y[start_pnt:end_pnt])
  stats::spline(start_pnt:end_pnt, y, xout=1, method = "natural")$y
}

# return the sum of squares of a complex vector
cplx_sum_sq <- function(x) sum(Re(x) ^ 2) + sum(Im(x) ^ 2)

#' Perform a polynomial fit, subtract and return the standard deviation of the
#' residuals.
#' @param y array.
#' @param degree polynomial degree.
#' @return standard deviation of the fit residuals.
#' @export
calc_sd_poly <- function(y, degree = 1) {
  x <- 1:length(y)
  model <- stats::lm(y ~ poly(x, degree))
  stats::sd(model$residual)
}

#' Compute the vector cross product between vectors x and y. Adapted from 
#' http://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function
#' @param x vector of length 3.
#' @param y vector of length 3.
#' @return vector cross product of x and y.
#' @export
crossprod_3d <- function(x, y) {
  i <- 1:3
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

#' Simulate MRS data with a similar appearance to normal brain (by default).
#' @param acq_paras list of acquisition parameters or an mrs_data object. See
#' \code{\link{def_acq_paras}}.
#' @param type type of spectrum, only "normal" is implemented currently.
#' @param pul_seq pulse sequence function to use.
#' @param xlim range of frequencies to simulate in ppm.
#' @param full_output when FALSE (default) only output the simulated MRS data.
#' When TRUE output a list containing the MRS data, basis set object and 
#' corresponding amplitudes.
#' @param amps a vector of basis amplitudes may be specified to modify the
#' output spectrum.
#' @param ... extra parameters to pass to the pulse sequence function.
#' @return see full_output option.
#' @export
sim_brain_1h <- function(acq_paras = def_acq_paras(), type = "normal_v1",
                         pul_seq = seq_press_ideal, xlim = c(0.5, 4.2), 
                         full_output = FALSE, amps = NULL,  ...) {
  
  
  
  if (is.null(amps)) {
    if (type == "normal_v1") {
      brain_basis_paras <- get_1h_brain_basis_paras_v1(ft = acq_paras$ft)
      amps <- c(0.000000000, 0.009799548, 0.072152490, 0.077845526, 0.045575002, 
                0.005450371, 0.000000000, 0.028636132, 0.076469056, 0.028382618,
                0.069602483, 0.001763720, 0.042031981, 0.013474549, 0.000000000, 
                0.000000000, 0.015705756, 0.001903173, 0.063409950, 0.051807236,
                0.222599530, 0.110431759, 0.023451957, 0.000000000, 0.032474090,
                0.007033074, 0.000000000)
    } else {
      stop("invalid spectrum type")
    }
  }
  
  basis <- sim_basis(brain_basis_paras, pul_seq, acq_paras, xlim = xlim, ...)

  mrs_data <- basis2mrs_data(basis, sum_elements = TRUE, amp = amps)
  
  if (!full_output) {
    return(mrs_data)
  } else {
    return(list(mrs_data = mrs_data, basis = basis, amps = amps))
  }
}

#' Combine the results from multiple csv format files into a table.
#' @param pattern glob string to match csv files.
#' @param supp_mess suppress messages from the read_csv function.
#' @param ... extra parameters to pass to read_csv.
#' @return results table.
#' @export
comb_csv_results <- function(pattern, supp_mess = TRUE, ...) {
  files <- Sys.glob(pattern)
  n <- length(files)
  
  if (n == 0) stop("No matching files found.")
  
  cat(paste(n, "file(s) found:\n"))
  
  # list the matches
  for (file in files) cat(file, "\n")
  
  # read the first file
  if (supp_mess) {
    res <- suppressMessages(readr::read_csv(files[1], ...))
  } else {
    res <- readr::read_csv(files[1], ...)
  }
  
  if (length(files) > 1) {
    for (n in 2:length(files)) {
      if (supp_mess) {
        res_temp <- suppressMessages(readr::read_csv(files[n], ...))
      } else {
        res_temp <- readr::read_csv(files[n])
      }
      res <- rbind(res, res_temp)
    }
  }
  
  # add the filename
  res <- tibble::add_column(res, files, .before = 1)
  res
}

#' Get the point spread function (PSF) for a 2D phase encoded MRSI scan.
#' @param FOV field of view in mm.
#' @param mat_size acquisition matrix size (not interpolated).
#' @param sampling can be either "circ" for circular or "rect" for rectangular.
#' @param hamming should Hamming k-space weighting be applied (default FALSE).
#' @return A matrix of the PSF with 1mm resolution.
#' @export
get_2d_psf <- function(FOV = 160, mat_size = 16, sampling = "circ",
                       hamming = FALSE) {
  
  if (sampling == "circ") {
    g <- expand.grid(1:FOV, 1:FOV)
    dist <- sqrt((g$Var1 - FOV / 2 - 0.5) ^ 2 + (g$Var2 - FOV / 2 - 0.5) ^ 2)
    dist <- matrix(dist, FOV, FOV)
    kspace <- 1 * (dist <= mat_size / 2)
  } else if (sampling == "rect") {
    kspace <- matrix(0, nrow = FOV, ncol = FOV)
    start_ind <- FOV / 2 - mat_size / 2 + 1
    inds <- seq(from = start_ind, by = 1, length.out = mat_size)
    kspace[inds, inds] <- 1
  } else {
    stop("Incorrect sampling option.")
  }
  
  if (hamming) {
    ham_mat <- pracma::repmat(signal::hamming(mat_size), mat_size, 1)
    ham_mat <- ham_mat * t(ham_mat)
    inds <- seq(from = start_ind, by = 1, length.out = mat_size)
    kspace[inds, inds] <- kspace[inds, inds] * ham_mat
  }
  
  kspace_shift <- apply(kspace, 1, pracma::ifftshift)
  kspace_shift <- apply(kspace_shift, 1, pracma::ifftshift)
  imspace <- apply(kspace_shift, 1, stats::fft, FALSE)
  imspace <- apply(imspace, 1, stats::fft, FALSE)
  imspace <- apply(imspace, 1, pracma::fftshift)
  imspace <- apply(imspace, 1, pracma::fftshift)
  imspace <- imspace / max(Re(imspace))
  
  return(imspace)
}

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom magrittr %$%
#' @export
magrittr::`%$%`

#' @importFrom RNifti readNifti
#' @export
RNifti::readNifti

# prob not as fast as functions in matrixStats but better than
# apply
colSdColMeans <- function(x, na.rm) {
  if (na.rm) {
    n <- colSums(!is.na(x)) # thanks @flodel
  } else {
    n <- nrow(x)
  }
  colVar <- colMeans(x * x, na.rm = na.rm) - (colMeans(x, na.rm = na.rm)) ^ 2
  return(sqrt(colVar * n / (n - 1)))
}

# Rotation matrix from axis and angle
# https://en.wikipedia.org/wiki/Rotation_matrix#cite_note-5
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
  as.numeric(vec_out / sum((vec_out) ^ 2) ^ 0.5)
}

cross <- function(a, b) {
  vec_out <- c(NA, NA, NA)
  vec_out[1] <- a[2] * b[3] - a[3] * b[2]
  vec_out[2] <- a[3] * b[1] - a[1] * b[3]
  vec_out[3] <- a[1] * b[2] - a[2] * b[1]
  vec_out <- vec_out / (sum(vec_out ^ 2) ^ 0.5)
  vec_out
}

#' Create a logical circular mask spanning the full extent of an n x n matrix.
#' @param d diameter of the mask.
#' @param n number of matrix rows and columns.
#' @param offset offset the mask center in matrix dimension units.
#' @return logical n x n mask matrix.
#' @export
circ_mask <- function(d, n, offset = 1) {
  g <- expand.grid(1:n, 1:n)
  dist <- sqrt((g$Var1 - n / 2 - offset) ^ 2 + (g$Var2 - n / 2 - offset) ^ 2)
  dist <- matrix(dist, n, n)
  return(dist <= d / 2)
}