#' Find valid MRS data files recursively from a directory path.
#' @param dir a directory path.
#' @return a vector of valid MRS data files.
#' @export
find_mrs_files <- function(dir) {
  
  # recursively find all files
  all_files <- list.files(dir, recursive = TRUE, full.names = TRUE)
  
  # empty vector of found files
  mrs_files <- c()
  
  # check each file, catching any errors
  for (file in all_files) {
    tryCatch(
      {
        dummy <- read_mrs(file)
        mrs_files <- c(mrs_files, file)
      },
      error = function(cond) {
      }
    )
  }
  return(mrs_files)
}

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

#' Simulate an ideal pulse excitation profile by smoothing a top-hat function 
#' with a Gaussian.
#' @param bw top-hat bandwidth (Hz).
#' @param sigma Gaussian width smoothing parameter (Hz).
#' @param fa intended flip angle of the pulse.
#' @return data frame containing the frequency scale, excitation profile and
#' corresponding flip-angles.
#' @export
sim_th_excit_profile <- function(bw = 1500, sigma = 50, fa = 180) {
  y <- c(rep(0, bw / 2), rep(1, bw), rep(0, bw / 2))
  if (sigma != 0) y <- mmand::gaussianSmooth(y, sigma)
  width <- length(y)
  x <- seq(from = -width / 2, to = width / 2, length.out = width)
  y[y >  1] <-  1 # acos doesn't like values gt 1
  y[y < -1] <- -1 # acos doesn't like values lt -1
  angle <- acos(1 - 2 * y) * fa / pi
  return(data.frame(freq = x, Mxy = y, fa = angle))
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
  if (is.na(data_pts[1])) return(NA)
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
  
  # check for masked data
  if (is.na(y[1])) return(NA)
  
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
  stats::spline(start_pnt:end_pnt, y, xout = 1, method = "natural")$y
}

measure_td_amp_poly <- function(y, start_pnt = 10, end_pnt = 50) {
  
  # check for masked data
  if (is.na(y[1])) return(NA)
  
  x <- start_pnt:end_pnt
  y <- Mod(y[x])
  lm_res <- stats::lm(y ~ poly(x, 3))
  new <- data.frame(x = 1)
  return(stats::predict(lm_res, new))
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

l2_norm_vec <- function(x) {
  return(x / norm(x, "2"))
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
#' @param basis_lb apply additional Gaussian line-broadening to the basis (Hz).
#' @param zero_lip_mm zero the amplitudes of any lipid or macromolecular
#' components based on their name starting with "MM" or "Lip".
#' @param remove_lip_mm remove any lipid or macromolecular basis components
#' based on their name starting with "MM" or "Lip".
#' @param ... extra parameters to pass to the pulse sequence function.
#' @return see full_output option.
#' @export
sim_brain_1h <- function(acq_paras = def_acq_paras(), type = "normal_v2",
                         pul_seq = seq_slaser_ideal, xlim = c(0.5, 4.2), 
                         full_output = FALSE, amps = NULL,
                         basis_lb = NULL, zero_lip_mm = FALSE, 
                         remove_lip_mm = FALSE, ...) {
  
  if (type == "normal_v1") {
    brain_basis_paras <- get_1h_brain_basis_paras_v1(ft = acq_paras$ft)
    if (is.null(amps)) {
      amps <- c(0.000000000, 0.009799548, 0.072152490, 0.077845526, 0.045575002, 
                0.005450371, 0.000000000, 0.028636132, 0.076469056, 0.028382618,
                0.069602483, 0.001763720, 0.042031981, 0.013474549, 0.000000000, 
                0.000000000, 0.015705756, 0.001903173, 0.063409950, 0.051807236,
                0.222599530, 0.110431759, 0.023451957, 0.000000000, 0.032474090,
                0.007033074, 0.000000000)
    }
  } else if (type == "normal_v2") {
    brain_basis_paras <- get_1h_brain_basis_paras_v3(ft = acq_paras$ft)
    if (is.null(amps)) {
      # levels are taken from within the range specified by "In Vivo NMR
      # Spectroscopy 3rd ed" in normal grey matter 
      amps <- c("-CrCH2" = 0,
                "Ala"    = 0.2, 
                "Asp"    = 1.5,
                "Cr"     = 5,
                "GABA"   = 1,
                "Glc"    = 1,
                "Gln"    = 3,
                "Gly"    = 1,
                "GSH"    = 2,
                "Glu"    = 10,
                "GPC"    = 1.2,
                "Ins"    = 5,
                "Lac"    = 0.5,
                "Lip09"  = 0,
                "Lip13a" = 0,
                "Lip13b" = 0,
                "Lip20"  = 0,
                "MM09"   = 4,
                "MM12"   = 4,
                "MM14"   = 4,
                "MM17"   = 4,
                "MM20"   = 4,
                "NAA"    = 10,
                "NAAG"   = 2,
                "PCh"    = 0.5,
                "PCr"    = 4.5,
                "PEth"   = 0.8,
                "sIns"   = 0.4,
                "Tau"    = 1.5)
    }
  } else {
    stop("invalid spectrum type")
  }
  
  basis <- sim_basis(brain_basis_paras, pul_seq, acq_paras, xlim = xlim, ...)
  
  if (zero_lip_mm) {
    if (length(grep("^MM", basis$names)) > 0) {
      zero_idx <- grep("^MM", basis$names)
      amps[zero_idx] <- 0
    }
    
    if (length(grep("^Lip", basis$names)) > 0) {
      zero_idx <- grep("^Lip", basis$names)
      amps[zero_idx] <- 0
    }
  }
  
  if (remove_lip_mm) {
    if (length(grep("^MM", basis$names)) > 0) {
      rem_idx <- grep("^MM", basis$names)
      amps <- amps[-rem_idx]
      basis$names <- basis$names[-rem_idx]
      basis$data  <- basis$data[, -rem_idx]
    }
    
    if (length(grep("^Lip", basis$names)) > 0) {
      rem_idx <- grep("^Lip", basis$names)
      amps <- amps[-rem_idx]
      basis$names <- basis$names[-rem_idx]
      basis$data  <- basis$data[, -rem_idx]
    }
  }
  
  names(amps) <- basis$names
  
  if (!is.null(basis_lb)) basis <- lb(basis, basis_lb)

  mrs_data <- basis2mrs_data(basis, sum_elements = TRUE, amps = amps)
  
  if (!full_output) {
    return(mrs_data)
  } else {
    return(list(mrs_data = mrs_data, basis = basis, amps = amps))
  }
}

#' Get the point spread function (PSF) for a 2D phase encoded MRSI scan.
#' @param FOV field of view in mm.
#' @param mat_size acquisition matrix size (not interpolated).
#' @param sampling can be either "circ" for circular or "rect" for rectangular.
#' @param hamming should Hamming k-space weighting be applied (default FALSE).
#' @param ensure_odd add 1mm to the FOV when required to ensure the output pdf
#' has odd dimensions. Required when using get_mrsi2d_seg.
#' @return A matrix of the PSF with 1mm resolution.
#' @export
get_2d_psf <- function(FOV = 160, mat_size = 16, sampling = "circ",
                       hamming = FALSE, ensure_odd = TRUE) {
  
  # round FOV to ensure int value
  FOV <- round(FOV)
  
  if (ensure_odd) {
    if ((FOV %% 2) == 0) FOV <- FOV + 1
  }
  
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

#' Mask fit result spectra depending on a vector of bool values.
#' @param fit_result fit result object to be masked.
#' @param mask_vec a Boolean vector with the same number of rows as there are
#' rows in the results table.
#' @param amps_only only mask the amplitude and associated error estimate
#' columns.
#' @return a masked fit result object.
#' @export
mask_fit_res <- function(fit_result, mask_vec, amps_only = FALSE) {
  if (length(mask_vec) != nrow(fit_result$res_tab)) {
    stop("mask_vec and fit_res dimensions do not agree")
  }

  # mask any NA values too
  mask_vec[is.na(mask_vec)] <- TRUE
  
  if (amps_only) {
    fit_result$res_tab[mask_vec, 6:(fit_result$amp_cols * 2 + 5)] <- NA
  } else {
    fit_result$res_tab[mask_vec, 6:ncol(fit_result$res_tab)] <- NA
  }
  
  return(fit_result)
}

#' @importFrom RNifti readNifti
#' @export
RNifti::readNifti

#' @importFrom RNifti writeNifti
#' @export
RNifti::writeNifti

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
#' @param offset offset the mask centre in matrix dimension units.
#' @return logical n x n mask matrix.
#' @export
circ_mask <- function(d, n, offset = 1) {
  g <- expand.grid(1:n, 1:n)
  dist <- sqrt((g$Var1 - n / 2 - offset) ^ 2 + (g$Var2 - n / 2 - offset) ^ 2)
  dist <- matrix(dist, n, n)
  return(dist <= d / 2)
}

round_dp <- function(x, dp) format(round(x, dp), nsmall = dp)

# equation for following function taken from:
# https://math.stackexchange.com/questions/2645689/what-is-the-parametric-equation-of-a-rotated-ellipse-given-the-angle-of-rotatio

#' Create an elliptical mask stored as a matrix of logical values.
#' @param xN number of pixels in the x dimension.
#' @param yN number of pixels in the y dimension.
#' @param x0 centre of ellipse in the x direction in units of pixels.
#' @param y0 centre of ellipse in the y direction in units of pixels.
#' @param xr radius in the x direction in units of pixels.
#' @param yr radius in the y direction in units of pixels.
#' @param angle angle of rotation in degrees.
#' @return logical mask matrix with dimensions fov_yN x fov_xN.
#' @export
elliptical_mask <- function(xN, yN, x0, y0, xr, yr, angle) {
  g <- pracma::meshgrid(0:(xN - 1) - xN / 2, 0:(yN - 1) - yN / 2)
  phi <- angle / 180 * pi
  mask <- (((g$X - x0) * cos(phi) + (g$Y - y0) * sin(phi))  ^ 2 / xr ^ 2 + 
          (((g$X - x0) * sin(phi) - (g$Y - y0) * cos(phi))) ^ 2 / yr ^ 2) <= 1
  return(mask)
}

#' Create a rectangular mask stored as a matrix of logical values.
#' @param xN number of pixels in the x dimension.
#' @param yN number of pixels in the y dimension.
#' @param x0 centre of rectangle in the x direction in units of pixels.
#' @param y0 centre of rectangle in the y direction in units of pixels.
#' @param xw width in the x direction in units of pixels.
#' @param yw width in the y direction in units of pixels.
#' @param angle angle of rotation in degrees.
#' @return logical mask matrix with dimensions fov_yN x fov_xN.
#' @export
rectangular_mask <- function(xN, yN, x0, y0, xw, yw, angle) {
  g <- pracma::meshgrid(0:(xN - 1) - xN / 2, 0:(yN - 1) - yN / 2)
  phi <- angle / 180 * pi
  
  x_mask <- (((g$X - x0) * cos(phi) + (g$Y - y0) * sin(phi)) ^ 2 / 
               (xw / 2) ^ 2 <= 1) 
  
  y_mask <- (((g$X - x0) * sin(phi) - (g$Y - y0) * cos(phi)) ^ 2 /
               (yw / 2) ^ 2 <= 1)
  
  mask <- x_mask & y_mask
  return(mask)
}

#' Create a two dimensional Gaussian window function stored as a matrix.
#' @param xN number of pixels in the x dimension.
#' @param yN number of pixels in the y dimension.
#' @param x0 centre of window function in the x direction in units of pixels.
#' Note, only integer values are applied.
#' @param y0 centre of window function in the y direction in units of pixels.
#' Note, only integer values are applied.
#' @param xw the reciprocal of the standard deviation of the Gaussian window in
#' x direction.
#' @param yw the reciprocal of the standard deviation of the Gaussian window in
#' y direction.
#' @return matrix with dimensions fov_yN x fov_xN.
#' @export
gausswin_2d <- function(xN, yN, x0, y0, xw, yw) {
  x_gaus_vec <- signal::gausswin(xN, xw)
  x_gaus_vec <- pracma::circshift(x_gaus_vec, x0)
  x_gaus_mat <- pracma::repmat(x_gaus_vec, yN, 1)
  
  y_gaus_vec <- signal::gausswin(yN, yw)
  y_gaus_vec <- pracma::circshift(y_gaus_vec, y0)
  y_gaus_mat <- pracma::repmat(y_gaus_vec, xN, 1)
  y_gaus_mat <- t(y_gaus_mat)

  return(x_gaus_mat * y_gaus_mat)
}

#' Check if an object is defined, which is the same as being not NULL.
#' @param x object to test for being NULL.
#' @return logical value.
is.def <- function(x) !is.null(x)

convolve_td <- function(x, y) {
  x <- stats::fft(x * Conj(y) * length(x))
  return(x)
}

# Directly from the matrixcalc package
# I don't like copying code like this, but CRAN don't like packages
# that have lots of dependencies so...
hankel.matrix <- function(n, x) {
  if (!is.numeric(n)) 
    stop("argument n is not numeric")
  if (n != trunc(n)) 
    stop("argument n ix not an integer")
  if (n < 2) 
    stop("argument n is less than 2")
  if (!is.vector(x)) 
    stop("argument x is not a vector")
  m <- length(x)
  if (m < n) 
    stop("length of argument x is less than n")
  y <- x
  H <- matrix(0, nrow = n, ncol = n)
  H[1, ] <- y[1:n]
  for (i in 2:n) {
    y <- c(y[2:m], y[1])
    H[i, ] <- y[1:n]
  }
  return(H)
}

# Directly from the MASS package
ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                           t(Xsvd$u[, Positive, drop = FALSE]))
}

#' Matrix exponential function taken from complexplus package to reduce the
#' number of spant dependencies.
#' @param x a square complex matrix.
#' @return the matrix exponential of x.
#' @export
matexp <- function(x) {
  if (utils::compareVersion(as.character(getRversion()), "4.4.0") >= 0) {
    d  <- dim(x)
    n  <- d[1]
    eA <- expm::expm(x) # complex matrices supported since expm v1.0
  } else {
    # older R versions don't support complex norm so have to use this slower 
    # method
    d  <- dim(x)
    n  <- d[1]
    p  <- d[2]
    Ar <- Re(x)
    Ai <- Im(x)
    E  <- rbind(cbind(Ar, -Ai), cbind(Ai, Ar))
    eE <- expm::expm(E)
    eA <- eE[1:n, 1:n] + (0 + 1i) * eE[(1:n) + n, 1:n]
  }
  eA <- matrix(Imzap(eA), ncol = n)
  return(eA)
}

#' Complex rounding function taken from complexplus package to reduce the number
#' of spant dependencies.
#' @param x a scalar or vector, real or complex.
#' @param tol a tolerance, 10^-6 by default. Prevents possible numerical
#' problems. Can be set to 0 if desired.
#' @export
Imzap <- function(x, tol = 1e-06) {
  if (all(abs(Im(z <- zapsmall(x))) <= tol)) 
    suppressWarnings(as.double(x))
  else x
}

# a trick to allow nice interpolation to work with data containing NAs
interpolate_nas <- function(map, factor, sigma) {
  map_box_interp <- mmand::rescale(map, factor, mmand::boxKernel())
  map_gaus       <- mmand::gaussianSmooth(map, sigma)
  map_gaus[is.na(map_gaus)] <- 0 # set any remaing NA values to 0
  map_replace    <- map
  map_replace[is.na(map_replace)] <- map_gaus[is.na(map_replace)]
  map_interp     <- mmand::rescale(map_replace, factor, mmand::mnKernel())
  map_interp[is.na(map_box_interp)] <- NA
  return(map_interp)
}

# is this dataset SVS?
is_svs <- function(mrs_data) identical(dim(mrs_data$data)[2:4], c(1L, 1L, 1L))

#' Calculate the mean of adjacent blocks in a vector.
#' @param x input vector.
#' @param block_size number of adjacent elements to average over.
#' @return vector data averaged in blocks.
#' @export 
mean_vec_blocks <- function(x, block_size) {
  
  if ((length(x) %% block_size) != 0) {
    warning("Block size does not fit into vector length without truncation.")
  }
  
  new_dyns <- floor(length(x) / block_size)
  x_out    <- x[seq(1, new_dyns * block_size, block_size)]
  
  for (n in 2:block_size) {
    x_out <- x_out + x[seq(n, new_dyns * block_size, block_size)]
  }
  
  return(x_out / block_size)
}

#' Output a table of fit amplitudes and error estimates for a single-voxel fit.
#' @param fit_res input vector.
#' @param format_out reduce the accuracy of values to aid table formatting.
#' @return data.frame of values.
#' @export
sv_res_table <- function(fit_res, format_out = FALSE) {
  
  # pretend it's a full fit result if only a table is passed in
  if (is.data.frame(fit_res)) fit_res <- list(res_tab = fit_res)
  
  first_sig        <- names(fit_res$res_tab)[6]
  first_sig_sd     <- paste0(first_sig, ".sd")
  first_sig_sd_idx <- which(names(fit_res$res_tab) == first_sig_sd)
  amp_indx <- 6:(first_sig_sd_idx - 1)
  amps     <- as.numeric(fit_res$res_tab[1, amp_indx])
  names    <- colnames(fit_res$res_tab[1, amp_indx])
  sds      <- as.numeric(fit_res$res_tab[1, amp_indx + 
                                            amp_indx[length(amp_indx)] - 5])
  sds_perc <- sds / amps * 100
  CI95_UB <- amps + sds * 1.96
  CI95_LB <- amps - sds * 1.96
  CI95_LB[CI95_LB < 0] <- 0
  
  if (format_out) {
    df_out <- data.frame(amps, sds, sds_perc, CI95_LB, CI95_UB,
                         row.names = names)
    df_out <- format(df_out, digits = 3)
    df_out$sds_perc <- round(sds_perc)
    df_out$sds_perc[sds_perc > 999] <- "999"
    df_out$sds_perc <- paste0(df_out$sds_perc, "%")
    df_out$CI95 <- paste0("[", gsub(" ", "", df_out$CI95_LB), ", ",
                          gsub(" ", "", df_out$CI95_UB), "]")
  } else {
    df_out <- data.frame(amps, sds, sds_perc, CI95_LB, CI95_UB,
                         row.names = names)
  }
  
  return(df_out)
}

add_fit_res_tab_amp_sd <- function(res_tab, name, amp_col, sd_col) {
  first_sig_idx    <- 6
  first_sig        <- names(res_tab[first_sig_idx])
  first_sig_sd     <- paste0(first_sig, ".sd")
  first_sig_sd_idx <- which(names(res_tab) == first_sig_sd)
  last_sig_idx     <- first_sig_sd_idx - 1
  sig_idxs         <- first_sig_idx:last_sig_idx
  last_sig_sd_idx  <- first_sig_sd_idx + length(sig_idxs) - 1
  sig_sd_idxs      <- first_sig_sd_idx:last_sig_sd_idx
  
  new_amp_col <- data.frame(amp_col)
  colnames(new_amp_col) <- name
  sd_name <- paste0(name, ".sd")
  new_sd_col <- data.frame(sd_col)
  colnames(new_sd_col) <- sd_name
  
  res_tab_out <- cbind(res_tab[1:last_sig_idx], new_amp_col,
                       res_tab[first_sig_sd_idx:last_sig_sd_idx], new_sd_col,
                       res_tab[(last_sig_sd_idx + 1):ncol(res_tab)])
  
  return(res_tab_out)
}

#' Segment T1 weighted MRI data using FSL FAST and write to file. Runs deface
#' and bet as preprocessing steps by default.
#' 
#' This function requires a working installation of FSL and uses the fslr 
#' package. You may need to specify the fsl install directory, eg: 
#' 'options(fsl.path = "/path/to/fsl")'
#' 
#' @param mri_path path to the volumetric T1 data.
#' @param out_dir optional output directory. Defaults to the same directory
#' as mri_path if not specified.
#' @param deface deface the input T1 data before analysis. Defaults to FALSE.
#' @param bet_fit fractional intensity threshold for bet brain extraction.
#' Values should be between 0 and 1. Defaults to 0.5 with smaller values giving
#' larger brain estimates.
#' @export
segment_t1_fsl <- function(mri_path, out_dir = NULL, deface = FALSE,
                           bet_fit = 0.5) {
  
  if (is.null(out_dir)) {
    dir_path <- dirname(mri_path)
  } else {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    dir_path <- out_dir
  }
  
  if (deface) {
    deface_path <- file.path(dir_path, "t1_deface")
    fslr::fsl_deface(mri_path, outfile = deface_path, retimg = FALSE,
                     verbose = FALSE)
    mri_path <- deface_path
  }
  brain_path <- file.path(dir_path, "t1_brain")
  bet_opts <- paste("-B -f", bet_fit)
  fslr::fslbet(mri_path, outfile = brain_path, retimg = FALSE, opts = bet_opts,
               betcmd = "bet", verbose = FALSE)
  seg_prefix <- file.path(dir_path, "t1")
  fslr::fslfast(brain_path, outfile = seg_prefix, retimg = FALSE,
                verbose = FALSE)
}

#' Read the log from Dinesh's MoCo sLASER sequence.
#' @param path path to the log file.
#' @param date scan date, eg "01-01-2025".
#' @param end_time scan time, eg "14:00".
#' @return data.frame of moco parameters.
#' @export
read_dkd_moco_log <- function(path, date, end_time) {
  
  lines <- readLines(path)
  
  for (n in 1:(length(lines) - 1)) {
    date_found <- (lines[n] == paste0(" # Scan date: ", date))
    time_found <- (lines[n + 1] == paste0(" # Scan end time: ", end_time))
    if (date_found & time_found) break
  }
  
  if (n == (length(lines) - 1)) stop("Date and time not found.")
  
  data_start <- n + 5
  n <- 0
  repeat {
    line <- lines[data_start + n]
    if (is.na(line)) break
    if (line == "") break
    n <- n + 1
  }
  
  data_end <- data_start + n - 1
  
  cnames <- trimws(gsub("[]\\[#;]", "", lines[data_start-2]))
  cnames <- gsub("Shim_DACunit:", "", cnames)
  cnames <- strsplit(cnames, " ")[[1]]
  
  data <- utils::read.table(textConnection(lines[data_start:data_end]))
  
  # remove weird characters
  rem_cols <- c(which(data[1,] == ";"), which(data[1,] == "["),
                which(data[1,] == "]"))
  
  data <- data[, -rem_cols]
  colnames(data) <- cnames
  return(data)
}

#' Split paths into parts based on backslash, forwardslash and dot characters
#' and return a data frame of these parts.
#' @param paths a vector of paths.
#' @param col_num when set, only the specified column number of the data frame
#' will be returned.
#' @param extra_regex extra regular expression for splitting the paths. Example,
#' "_|-" could be used to additionally split on underscores and hyphen
#' characters.
#' @return data.frame of separated path elements.
#' @export 
paths2df <- function(paths, col_num = NULL, extra_regex = NULL) {
  
  regex_str <- "/|\\\\|\\."
  
  if (!is.null(extra_regex)) regex_str <- paste0(regex_str, "|", extra_regex)
    
  str_list <- strsplit(paths, regex_str)
  str_list <- lapply(str_list, \(x) x[x!=""])
  out_df <- as.data.frame(do.call("rbind", str_list))
  # out_df <- cbind(paths = paths, out_df)
  if (is.null(col_num)) {
    return(out_df)
  } else {
    return(out_df[,col_num])
  }
}

#' Trim a vector of filesystem paths.
#' @param paths vectors of filesystem paths.
#' @param dir number of times to apply the base dirname function.
#' @param char number of characters to trim from the end of each path. Note this
#' is performed after the dir based trimming.
#' @return a vector of trimmed paths.
#' @export 
trim_paths <- function(paths, dir = 0, char = 0) {
  
  for (n in 1:length(paths)) {
    if (dir > 0) for (m in 1:dir) paths[n] <- dirname(paths[n])
    if (char > 0) paths[n] <- substr(paths[n], 1, nchar(paths[n]) - char)
  }
    
  return(paths)
}

#' Match files based on a vector of input paths and a glob pattern. The glob
#' pattern is appended to each path and should match one file only.
#' @param paths vectors of filesystem paths.
#' @param glob pattern to append to each path before passing to Sys.glob.
#' @return matched files. NA values correspond to either no match or multiple
#' matches.
#' @export 
match_files <- function(paths, glob) {
  paths_out <- rep(NA, length(paths))
  for (n in 1:length(paths)) {
    match <- Sys.glob(paste0(paths[n], glob))
    if (length(match) == 0) {
      warning(paste0("no match for : ", paste0(paths[n], glob)))
    } else if (length(match) > 1) {
      warning(paste0("multiple matches for : ", paste0(paths[n], glob)))
    } else if (length(match) == 1) {
      paths_out[n] <- match
    } else {
      stop("I don't belong here")
    }
  }
  return(paths_out)
}