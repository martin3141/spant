
check_mrs_data <- function(mrs_data) {
  in_class <- class(mrs_data)[1]
  if (in_class == "mrs_data") {
    return()
  } else if (in_class == "list") {
    if (class(mrs_data[[1]]) == "mrs_data") {
      stop("Error, input is a list of mrs_data objects. Please only pass a
           single mrs_data object to this function.")
    } else {
      stop("Error, input is not mrs_data class.")
    }
  } else {
    stop("Error, input is not mrs_data class.")
  }
}

#' Simulate a MRS data object containing a set of simulated resonances.
#' @param freq resonance frequency.
#' @param amp resonance amplitude.
#' @param lw line width in Hz.
#' @param lg Lorentz-Gauss lineshape parameter (between 0 and 1).
#' @param phase phase in degrees.
#' @param freq_ppm frequencies are given in ppm units if set to TRUE, otherwise
#' Hz are assumed.
#' @param acq_paras list of acquisition parameters. See
#' \code{\link{def_acq_paras}}
#' @param fp_scale multiply the first data point by 0.5.
#' @param back_extrap_pts number of data points to back extrapolate.
#' @param sum_resonances sum all resonances (default is TRUE), otherwise return
#' a dynamic mrs_data object. 
#' @return MRS data object.
#' @examples
#' sim_data <- sim_resonances(freq = 2, lw = 5)
#' @export
sim_resonances <- function(freq = 0, amp = 1, lw = 0, lg = 0, phase = 0, 
                           freq_ppm = TRUE, acq_paras = def_acq_paras(),
                           fp_scale = TRUE, back_extrap_pts = 0,
                           sum_resonances = TRUE) {
  
  if (class(acq_paras) == "mrs_data") acq_paras <- get_acq_paras(acq_paras)
  
  sig_n <- length(freq)
  if (sig_n != length(amp)) {
    if (length(amp) != 1) warning("amp length is not the same as freq length")
    amp <- rep_len(amp, sig_n)
  }
  
  if (sig_n != length(lw)) {
    if (length(lw) != 1) warning("lw length is not the same as freq length")
    lw <- rep_len(lw, sig_n)
  }
  
  if (sig_n != length(phase)) {
    if (length(phase) != 1) {
      warning("phase length is not the same as freq length")
    }
    phase <- rep_len(phase, sig_n)
  }
  
  if (sig_n != length(lg)) {
    if (length(lg) != 1) warning("lg length is not the same as freq length")
    lg <- rep_len(lg, sig_n)
  }
  
  # generate data in TD
  t <- seq(from = -back_extrap_pts / acq_paras$fs,
           to = (acq_paras$N - 1) / acq_paras$fs, by = 1 / acq_paras$fs)
  
  # covert freqs to Hz
  if (freq_ppm) {
    f_hz <- (acq_paras$ref - freq) * acq_paras$ft / 1e6
  } else {
    f_hz <- freq
  }
  
  data <- matrix(NA, sig_n, acq_paras$N + back_extrap_pts)
  
  for (n in 1:sig_n) {
    data[n,] <- amp[n] * exp(1i * pi * phase[n] / 180 + 2i * pi * f_hz[n] * t)
    
    # LG peak model
    data[n,] <- data[n,] * ((1 - lg[n]) * exp(-lw[n] * t * pi) + 
                             lg[n] * exp(-lw2beta(lw[n]) * t * t))
  }
  
  # first point correction
  if (fp_scale) data[, 1] <- data[, 1] * 0.5
  
  mrs_data <- mat2mrs_data(data, fs = acq_paras$fs, ft = acq_paras$ft,
                           ref = acq_paras$ref, nuc = acq_paras$nuc,
                           fd = FALSE)
  
  if (sum_resonances) mrs_data <- sum_dyns(mrs_data)
  
  return(mrs_data)
}

sim_resonances_fast <- function(freq = 0, amp = 1, freq_ppm = TRUE,
                                N = def_N(), fs = def_fs(), ft = def_ft(),
                                ref = def_ref(), nuc = def_nuc()) {
  
  sig_n <- length(freq)
  if (sig_n != length(amp)) {
    amp <- rep_len(amp, sig_n)
  }
  
  # generate data in TD
  t <- seq(from = 0, to = (N - 1) / fs, by = 1 / fs)
  
  # covert freqs to Hz
  if (freq_ppm) {
    f_hz <- (ref - freq) * ft / 1e6
  } else {
    f_hz <- freq
  }
  
  data <- rep(0, N)
  for (n in 1:sig_n) {
    temp_data <- amp[n] * exp(2i * pi * f_hz[n] * t)
    data <- data + temp_data
  }
  
  # first point correction
  data[1] <- data[1] * 0.5
  
  data <- array(data,dim = c(1, 1, 1, 1, 1, 1, N))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
    
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = rep(FALSE, 7), affine = NULL,
                       meta = NULL, extra = NULL)
  
  return(mrs_data)
}

sim_resonances_fast2 <- function(freq = 0, amp = 1, freq_ppm = TRUE,
                                 N = def_N(), fs = def_fs(), ft = def_ft(), 
                                 ref = def_ref(), nuc = def_nuc()) {
  
  sig_n <- length(freq)
  if (sig_n != length(amp)) {
    amp <- rep_len(amp, sig_n)
  }
   # covert freqs to Hz
  if (freq_ppm) {
    f_hz <- (ref - freq) * ft / 1e6
  } else {
    f_hz <- freq
  }
  
  # generate data in TD
  t <- seq(from = 0, to = (N - 1) / fs, by = 1 / fs)
  t_i_omega <- t * 2i * pi
  
  # TODO should be able to loose a tranpose here
  t_i_omega_mat <- matrix(t_i_omega, nrow = N, ncol = sig_n)
  f_hz_mat <- matrix(f_hz, nrow = N, ncol = sig_n, byrow = TRUE)
  temp <- t_i_omega_mat * f_hz_mat
  #temp <- t(t(t_i_omega_mat) * f_hz)
  td_sig <- exp(temp) 
  #td_sig <- matrix(exp(as.vector(temp)), nrow = N, ncol = sig_n, byrow = FALSE)
  #expp <- exp(1)
  #e <- matrix(expp, nrow = N, ncol = sig_n)
  #td_sig <- e ^ temp
  data <- td_sig %*% amp 
  
  # first point correction
  data[1] <- data[1] * 0.5
  
  data <- array(data, dim = c(1, 1, 1, 1, 1, 1, N))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = rep(FALSE, 7), affine = NULL,
                       meta = NULL, extra = NULL)
  
  return(mrs_data)
}

#' Convert a vector into a mrs_data object.
#' @param vec the data vector.
#' @param fs sampling frequency in Hz.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param nuc resonant nucleus.
#' @param dyns replicate the data across the dynamic dimension.
#' @param fd flag to indicate if the matrix is in the frequency domain (logical).
#' @return mrs_data object.
#' @export
vec2mrs_data <- function(vec, fs = def_fs(), ft = def_ft(), ref = def_ref(),
                         nuc = def_nuc(), dyns = 1, fd = FALSE) {
  
  data <- array(vec, dim = c(length(vec), dyns))
  data <- aperm(data,c(2, 1))
  dim(data) <- c(1, 1, 1, 1, dyns, 1, length(vec))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = c(rep(FALSE, 6), fd),
                       affine = NULL, meta = NULL, extra = NULL)
  
  return(mrs_data)
}

#' Convert a 7 dimensional array in into a mrs_data object. The array dimensions
#' should be ordered as : dummy, X, Y, Z, dynamic, coil, FID.
#' @param data_array 7d data array.
#' @param fs sampling frequency in Hz.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param nuc nucleus that is resonant at the transmitter frequency.
#' @param fd flag to indicate if the matrix is in the frequency domain (logical).
#' @return mrs_data object.
#' @export
array2mrs_data <- function(data_array, fs = def_fs(), ft = def_ft(),
                           ref = def_ref(), nuc = def_nuc(), fd = FALSE) {
  
  if (length(dim(data_array)) != 7) stop("Incorrect number of dimensions.")
  
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  
  mrs_data <- mrs_data(data = data_array, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = c(rep(FALSE, 6), fd),
                       affine = NULL, meta = NULL, extra = NULL)
  
  return(mrs_data)
}

#' Convert mrs_data object to a matrix, with spectral points in the column
#' dimension and dynamics in the row dimension.
#' @param mrs_data MRS data object or list of MRS data objects.
#' @param collapse collapse all other dimensions along the dynamic dimension, eg
#' a 16x16 MRSI grid would be first collapsed across 256 dynamic scans.
#' @return MRS data matrix.
#' @export
mrs_data2mat <- function(mrs_data, collapse = TRUE) {
  if (class(mrs_data) == "list") mrs_data <- append_dyns(mrs_data)
  
  if (collapse) mrs_data <- collapse_to_dyns(mrs_data)
  
  as.matrix(mrs_data$data[1,1,1,1,,1,])
}

#' Convert mrs_data object to a vector.
#' @param mrs_data MRS data object.
#' @param dyn dynamic index.
#' @param x_pos x index.
#' @param y_pos y index.
#' @param z_pos z index.
#' @param coil coil element index.
#' @return MRS data vector.
#' @export
mrs_data2vec <- function(mrs_data, dyn = 1, x_pos = 1,
                          y_pos = 1, z_pos = 1, coil = 1) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  as.vector(mrs_data$data[1, x_pos, y_pos, z_pos, dyn, coil,])
}

#' Convert a matrix (with spectral points in the column dimension and dynamics
#' in the row dimensions) into a mrs_data object.
#' @param mat data matrix.
#' @param fs sampling frequency in Hz.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param nuc resonant nucleus.
#' @param fd flag to indicate if the matrix is in the frequency domain (logical).
#' @return mrs_data object.
#' @export
mat2mrs_data <- function(mat, fs = def_fs(), ft = def_ft(), ref = def_ref(),
                         nuc = def_nuc(), fd = FALSE) {
  
  data <- array(mat, dim = c(1, 1, 1, 1, nrow(mat), 1, ncol(mat)))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = c(rep(FALSE, 6), fd),
                       affine = NULL, meta = NULL, extra = NULL)
  
  return(mrs_data)
}

#' Simulate an mrs_data object containing simulated Gaussian noise.
#' @param sd standard deviation of the noise.
#' @param fs sampling frequency in Hz.
#' @param ft transmitter frequency in Hz.
#' @param N number of data points in the spectral dimension.
#' @param ref reference value for ppm scale.
#' @param dyns number of dynamic scans to generate.
#' @param fd return data in the frequency-domain (TRUE) or time-domain (FALSE)
#' @return mrs_data object.
#' @export
sim_noise <- function(sd = 0.1, fs = def_fs(), ft = def_ft(), N = def_N(),
                      ref = def_ref(), dyns = 1, fd = TRUE) {
 
  data_pts <- dyns * N 
  vec <- stats::rnorm(data_pts, 0, sd) + 1i*stats::rnorm(data_pts, 0, sd)
  data_array <- array(vec, dim = c(1, 1, 1, 1, dyns, 1, N))
  array2mrs_data(data_array, fs = fs, ft = ft, ref = ref, fd = fd)
}

#' Add noise to an mrs_data object.
#' @param mrs_data data to add noise to.
#' @param sd standard deviation of the noise.
#' @param fd generate the noise samples in the frequency-domain (TRUE) or
#' time-domain (FALSE). This is required since the absolute value of the 
#' standard deviation of noise samples changes when data is Fourier transformed.
#' @return mrs_data object with additive normally distributed noise.
#' @export
add_noise <- function(mrs_data, sd = 0.1, fd = TRUE) {
  
  # covert data to the frequency domain if needed
  if (fd & !is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  # covert data to the time domain if needed
  if (!fd & is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
 
  data_pts <- length(mrs_data$data)
  noise <- stats::rnorm(data_pts, 0, sd) + 1i*stats::rnorm(data_pts, 0, sd)
  mrs_data$data <- mrs_data$data + noise
  mrs_data
}

#' Simulate an mrs_data object containing complex zero valued samples.
#' @param fs sampling frequency in Hz.
#' @param ft transmitter frequency in Hz.
#' @param N number of data points in the spectral dimension.
#' @param ref reference value for ppm scale.
#' @param dyns number of dynamic scans to generate.
#' @return mrs_data object.
#' @export
sim_zero <- function(fs = def_fs(), ft = def_ft(), N = def_N(),
                      ref = def_ref(), dyns = 1) {
  
  data_pts <- dyns * N 
  vec <- rep(0, data_pts) * 1i
  data_array <- array(vec, dim = c(1, 1, 1, 1, dyns, 1, N))
  array2mrs_data(data_array, fs = fs, ft = ft, ref = ref)
}

#' Apply a function across given dimensions of a MRS data object.
#' @param mrs_data MRS data.
#' @param dims dimensions to apply the function.
#' @param fun name of the function.
#' @param ... arguments to the function.
#' @param data_only return an array rather than an MRS data object.
#' @export
apply_mrs <- function(mrs_data, dims, fun, ..., data_only = FALSE) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  dims <- sort(dims)
  margins <- c(1:7)[-dims]
  mrs_data$data <- plyr::aaply(mrs_data$data, margins, fun, ..., .drop = FALSE)
  perm_vec <- 1:(7 - length(dims))
  for (n in (1:length(dims))) {
    perm_vec <- append(perm_vec, (n + 7 - length(dims)), after = dims[n] - 1)
  }
  
  #print(perm_vec)
  mrs_data$data <- aperm(mrs_data$data, perm_vec)
  if (data_only == FALSE) {
    return(mrs_data)
  } else {
    return(mrs_data$data)
  }
}

#' Apply a frequency shift to MRS data.
#' @param mrs_data MRS data.
#' @param shift frequency shift (in ppm by default).
#' @param units of the shift ("ppm" or "hz").
#' @return frequency shifted MRS data.
#' @export
shift <- function(mrs_data, shift, units = "ppm") {
  
  # covert to time-domain
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  if (units == "hz") {
    shift_hz <- shift
  } else if (units == "ppm") {
    shift_hz <- ppm2hz(shift, mrs_data$ft, 0)
  } else {
    stop("Error, did not recognise the units.") 
  }
  
  t_orig <- rep(seconds(mrs_data), each = Nspec(mrs_data))
  t_array <- array(t_orig, dim = dim(mrs_data$data))
  
  if (length(shift_hz) == 1) {
    shift_array <- array(shift_hz, dim = dim(mrs_data$data))
  } else if (length(shift_hz) == Ndyns(mrs_data)) {
    # assume array should be applied in the dynamic dimension
    shift_array <- array(shift_hz, dim = c(1, 1, 1, 1, Ndyns(mrs_data), 1,
                                           Npts(mrs_data)))
  } else if (dim(shift_hz)[1:6] && dim(mrs_data$data)[1:6]) {
    shift_array <- array(shift_hz, dim = dim(mrs_data$data))
  } else {
    stop("Shift vector has an incorrect dimensions.")
  }
  
  shift_array <- exp(2i * pi * t_array * shift_array)
  mrs_data$data <- mrs_data$data * shift_array
  mrs_data
}

#' Apply phasing parameters to MRS data.
#' @param mrs_data MRS data.
#' @param zero_order zero'th order phase term in degrees.
#' @param first_order first order (frequency dependent) phase term in ms.
#' @return MRS data with applied phase parameters.
#' @export
phase <- function(mrs_data, zero_order, first_order = 0) {
  
  if (class(mrs_data) == "list") {
    res <- lapply(mrs_data, phase, zero_order = zero_order,
                  first_order = first_order)
    return(res)
  }
  
  # check the input
  check_mrs_data(mrs_data)
  
  if ((first_order == 0) && (length(zero_order) == 1)) { 
    # single zero order phase term given
    mrs_data$data <- mrs_data$data * exp(1i * zero_order * pi / 180)
  } else if ((first_order == 0) && (is.null(dim(zero_order)))) { 
    # array of zero order phase terms given
    if (length(zero_order) != Ndyns(mrs_data)) {
      stop("Shift vector has an incorrect length.")
    }
    # assume array should be applied in the dynamic dimension
    phase_array <- array(zero_order, dim = c(1, 1, 1, 1, Ndyns(mrs_data), 1,
                                             Npts(mrs_data)))
    mrs_data$data = mrs_data$data * exp(1i * phase_array * pi / 180)
  } else if ((first_order == 0) && (dim(zero_order)[1:6] == dim(mrs_data$data)[1:6])) {
    phase_array <- array(rep(zero_order, Npts(mrs_data)), 
                         dim = dim(mrs_data$data))
    mrs_data$data <- mrs_data$data * exp(1i * phase_array * pi / 180)
  } else if ((length(zero_order) == 1) && (first_order != 0)) {
    freq <- rep(hz(mrs_data), each = Nspec(mrs_data))
    freq_mat <- array(freq, dim = dim(mrs_data$data))
    # needs to be a freq-domain operation
    if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
    }
    mrs_data$data = mrs_data$data * exp(2i * pi * (zero_order / 360 - freq_mat
                                                   * first_order / 1000))
  } else {
    stop("Unsupported input options.")
  }
  return(mrs_data)
}

#' Perform a zeroth order phase correction based on the phase of the first data 
#' point in the time-domain.
#' @param mrs_data MRS data to be corrected.
#' @param ret_phase return phase values (logical).
#' @return corrected data or a list with corrected data and optional phase 
#' values.
#' @export
fp_phase_correct <- function(mrs_data, ret_phase = FALSE) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, fp_phase_correct, ret_phase = ret_phase))  
  }
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  phases <- Arg(mrs_data$data[,,,,,, 1, drop = F])
  mrs_data$data <- mrs_data$data * array(exp(-1i * phases), dim = dim(mrs_data))
  
  if (ret_phase) {
    return(list(mrs_data, 180 / pi * abind::adrop(phases, 7)))
  } else {
    return(mrs_data)
  }
}

#' Return the first time-domain data point.
#' @param mrs_data MRS data.
#' @return first time-domain data point.
#' @export
get_fp <- function(mrs_data) {
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  # drop the chem shift dimension
  mrs_data$data[,,,,,, 1, drop = F]
}

#' Scale the first time-domain data point in an mrs_data object.
#' @param mrs_data MRS data.
#' @param scale scaling value, defaults to 0.5.
#' @return scaled mrs_data object.
#' @export
fp_scale <- function(mrs_data, scale = 0.5) {
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  # drop the chem shift dimension
  mrs_data$data[,,,,,, 1] <- mrs_data$data[,,,,,, 1] * scale
  
  return(mrs_data)
}

fp_mag <- function(mrs_data) {
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  # drop the chem shift dimension
  abind::adrop(Mod(mrs_data$data[,,,,,, 1, drop = F]), 7)
}

#' Convolve two MRS data objects.
#' @param mrs_data MRS data to be convolved.
#' @param conv convolution data stored as an mrs_data object.
#' @return convolved data.
#' @export
conv_mrs <- function(mrs_data, conv) {
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  if (is_fd(conv)) conv <- fd2td(conv)
  
  if (Ndyns(mrs_data) > 1) {
    warning("Repeating convolution data to match mrs_data dynamics.")
    conv <- rep_dyn(conv, Ndyns(mrs_data))
  }
  
  mrs_data * conv 
}

#' Return the phase of the first data point in the time-domain.
#' @param mrs_data MRS data.
#' @return phase values in degrees.
#' @export
fp_phase <- function(mrs_data) {
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  # drop the chem shift dimension
  abind::adrop(Arg(mrs_data$data[,,,,,, 1, drop = F]), 7) * 180 / pi
}

#' Apply line-broadening (apodisation) to MRS data or basis object.
#' @param x input mrs_data or basis_set object.
#' @param lb amount of line-broadening in Hz.
#' @param lg Lorentz-Gauss lineshape parameter (between 0 and 1).
#' @return line-broadened data.
#' @rdname lb
#' @export
lb <- function(x, lb, lg = 1) UseMethod("lb")

#' @rdname lb
#' @export
lb.list <- function(x, lb, lg = 1) {
  res <- lapply(x, lb, lb = lb, lg = lg)
  return(res)
}

#' @rdname lb
#' @export
lb.mrs_data <- function(x, lb, lg = 1) {
  if (lg > 1 | lg < 0) {
    cat("Error, lg values not between 0 and 1.")  
    stop()
  }
  
  # collapse to simplify 
  lb <- as.vector(drop(lb))
  lg <- as.vector(drop(lg))
  orig_dim <- dim(x$data)
  orig_dimnames <- dimnames(x$data)
  x <- collapse_to_dyns(x)
  
  # needs to be a time-domain operation
  if (is_fd(x)) {
    x <- fd2td(x)
  }
  t <- rep(seconds(x), each = Nspec(x))
  
  if (lg < 1)  x$data = x$data * exp(-(1 - lg) * lb * t * pi)
  
  if (lg > 0) {
    sign <- ifelse(lb > 0, 1, -1)
    x$data = x$data * exp((sign * lg * lb ^ 2 * pi ^ 2 / 4 / log(0.5)) * 
                          (t ^ 2))
  }
  
  # revert back to original dims and unname
  dim(x$data) <- orig_dim
  
  return(x)
}

#' @rdname lb
#' @export
lb.basis_set <- function(x, lb, lg = 1) {
  mrs_data2basis(lb(basis2mrs_data(x), lb, lg), x$names)
}

#' Apply a weighting to the FID to enhance spectral resolution.
#' @param mrs_data data to be enhanced.
#' @param re resolution enhancement factor (rising exponential factor).
#' @param alpha alpha factor (Gaussian decay)
#' @return resolution enhanced mrs_data.
#' @export
re_weighting <- function(mrs_data, re, alpha) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, re_weighting, re = re, alpha = alpha))  
  }
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data) 
  
  t <- rep(seconds(mrs_data), each = Nspec(mrs_data))
  
  mrs_data$data = mrs_data$data * exp(re * t) * exp(-alpha * t ^ 2)
  
  return(mrs_data)
}

#' Zero-fill MRS data in the time domain.
#' @param x input mrs_data or basis_set object.
#' @param factor zero-filling factor, factor of 2 returns a dataset with
#' twice the original data points.
#' @return zero-filled data.
#' @rdname zf
#' @export
zf <- function(x, factor = 2) UseMethod("zf")

#' @rdname zf
#' @export
zf.list <- function(x, factor = 2) {
  lapply(x, zf, factor = factor)
}

#' @rdname zf
#' @export
zf.mrs_data <- function(x, factor = 2) {
  set_td_pts(x, factor * Npts(x))
}

#' @rdname zf
#' @export
zf.basis_set <- function(x, factor = 2) {
  x_mrs_data <- basis2mrs_data(x)
  mrs_data2basis(set_td_pts(x_mrs_data, factor * Npts(x_mrs_data)), x$names)
}

#' Set the number of time-domain data points, truncating or zero-filling as
#' appropriate.
#' @param mrs_data MRS data.
#' @param pts number of data points.
#' @return MRS data with pts data points.
#' @export
set_td_pts <- function(mrs_data, pts) {
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  data_dim <- dim(mrs_data$data)
  if (data_dim[7] > pts) {
    data_dim_trunc <- data_dim
    data_dim_trunc[7] <- pts
    mrs_data$data = array(mrs_data$data, data_dim_trunc)
  } else if (data_dim[7] < pts) {
    zero_dim <- data_dim
    zero_dim[7] <- pts - data_dim[7]
    zero_array <- array(0, dim = zero_dim)
    mrs_data$data = abind::abind(mrs_data$data, zero_array, along = 7)
  }
  dimnames(mrs_data$data) <- NULL
  return(mrs_data)
}

#' Set the ppm reference value (eg ppm value at 0Hz).
#' @param mrs_data MRS data.
#' @param ref reference value for ppm scale.
#' @export
set_ref <- function(mrs_data, ref) {
  
  # check the input
  check_mrs_data(mrs_data)
  
  mrs_data$ref = ref
  return(mrs_data)
}

#' Check if the chemical shift dimension of an MRS data object is in the
#' frequency domain.
#' @param mrs_data MRS data.
#' @return logical value.
#' @export
is_fd <- function(mrs_data) {
  
  # check the input is an mrs_data object
  check_mrs_data(mrs_data)
  
  mrs_data$freq_domain[7]
}

#' Apply the Fourier transform over the dynamic dimension.
#' @param mrs_data MRS data where the dynamic dimension is in the time-domain.
#' @param ft_shift apply FT shift to the output, default is FALSE.
#' @param ret_mod return the modulus out the transform, default is FALSE.
#' @param fd transform the chemical shift axis to the frequency domain first,
#' default is TRUE.
#' @return transformed MRS data.
#' @export
ft_dyns <- function(mrs_data, ft_shift = FALSE, ret_mod = FALSE, fd = TRUE) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, ft_dyns, ft_shift = ft_shift, ret_mod = ret_mod,
                  fd = fd))
  }
  
  # convert to the correct domain
  if (fd & !is_fd(mrs_data)) {
    mrs_data <- td2fd(mrs_data)
  } else if (!fd & is_fd(mrs_data)) {
    mrs_data <- fd2td(mrs_data)
  }
  
  data <- mrs_data$data
  
  # permute the dynamic dimension to be the last (7th)
  data <- aperm(data, c(1, 2, 3, 4, 6, 7, 5)) 
  data_dim <- dim(data)
  
  # convert to a matrix and transform
  data <- matrix(data, ncol = Ndyns(mrs_data))
  
  # fft and shift if requested
  if (ft_shift) {
    data <- t(ft_shift_mat(t(data)))
  } else {
    data <- t(stats::mvfft(t(data)))
  }
  
  # revert back to an array
  dim(data) <- data_dim
  
  # permute the dynamic dimension to be 5th
  data <- aperm(data, c(1, 2, 3, 4, 7, 5, 6)) 
  
  if (ret_mod) data <- Mod(data)
  
  mrs_data$data <- data
  
  return(mrs_data)
}

#' Transform time-domain data to the frequency-domain.
#' @param mrs_data MRS data in time-domain representation.
#' @return MRS data in frequency-domain representation.
#' @export
td2fd <- function(mrs_data) {
  
  if (class(mrs_data) == "list") return(lapply(mrs_data, td2fd))
  
  if (mrs_data$freq_domain[7] == TRUE) {
    warning("Data is alread in the frequency-domain.")
  }
  
  freq <- matrix(mrs_data$data, ncol = Npts(mrs_data))
  freq <- t(ft_shift_mat(t(freq)))
  dim(freq) <- dim(mrs_data$data)
  mrs_data$data <- freq
  
  mrs_data$freq_domain[7] <- TRUE
  unname(mrs_data$data)
  return(mrs_data)
}

#' Transform frequency-domain data to the time-domain.
#' @param mrs_data MRS data in frequency-domain representation.
#' @return MRS data in time-domain representation.
#' @export
fd2td <- function(mrs_data) {
  
  if (class(mrs_data) == "list") return(lapply(mrs_data, fd2td))
  
  if (mrs_data$freq_domain[7] == FALSE) {
    warning("Data is alread in the time-domain.")
  }
  
  time <- matrix(mrs_data$data, ncol = Npts(mrs_data))
  time <- t(ift_shift_mat(t(time)))
  dim(time) <- dim(mrs_data$data)
  mrs_data$data <- time
  
  mrs_data$freq_domain[7] <- FALSE
  dimnames(mrs_data$data) <- NULL
  return(mrs_data)
}

# recon complex td data from real part of fd data
recon_imag <- function(mrs_data) {
  # data needs to be in the FD
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  mrs_data <- apply_mrs(mrs_data, 7, recon_imag_vec)
  mrs_data$freq_domain[7] <- FALSE
  return(mrs_data)
}

#' Return acquisition parameters from a MRS data object.
#' @param mrs_data MRS data.
#' @return list of acquisition parameters.
#' @export
get_acq_paras <- function(mrs_data) {
  
  # check the input is an mrs_data object
  check_mrs_data(mrs_data)
  
  list(ft = mrs_data$ft, fs = fs(mrs_data), N = Npts(mrs_data),
       ref = mrs_data$ref, nuc = mrs_data$nuc)
}

ft <- function(mrs_data, dims) {
  apply_mrs(mrs_data, dims, ft_shift)
}

#' Apply the diff operator to an MRS dataset in the FID/spectral dimension.
#' @param mrs_data MRS data.
#' @param ... additional arguments to the diff function.
#' @return MRS data following diff operator.
#' @export
diff_mrs <- function(mrs_data, ...) {
  apply_mrs(mrs_data, 7, fun = diff, ...)
}

#' Apply the max operator to an MRS dataset.
#' @param mrs_data MRS data.
#' @return MRS data following max operator.
#' @export
max_mrs <- function(mrs_data) {
  apply_mrs(mrs_data, 7, max, data_only = TRUE)
}

#' Apply the max operator to an interpolated MRS dataset.
#' @param mrs_data MRS data.
#' @param interp_f interpolation factor.
#' @return Array of maximum values (real only).
#' @export
max_mrs_interp <- function(mrs_data, interp_f = 4) {
  apply_mrs(mrs_data, 7, re_max_interp, interp_f, data_only = TRUE)
}

#' Apply Re operator to an MRS dataset.
#' @param z MRS data.
#' @return MRS data following Re operator.
#' @export
Re.mrs_data <- function(z) {
  z$data <- Re(z$data)
  z
}

#' Apply Im operator to an MRS dataset.
#' @param z MRS data.
#' @return MRS data following Im operator.
#' @export
Im.mrs_data <- function(z) {
  z$data <- Im(z$data)
  z
}

#' Apply Mod operator to an MRS dataset.
#' @param z MRS data.
#' @return MRS data following Mod operator.
#' @export
Mod.mrs_data <- function(z) {
  z$data <- Mod(z$data)
  z
}

#' Apply Arg operator to an MRS dataset.
#' @param z MRS data.
#' @return MRS data following Arg operator.
#' @export
Arg.mrs_data <- function(z) {
  z$data <- Arg(z$data)
  z
}

#' Apply Conj operator to an MRS dataset.
#' @param z MRS data.
#' @return MRS data following Conj operator.
#' @export
Conj.mrs_data <- function(z) {
  z$data <- Conj(z$data)
  z
}

#' Downsample an MRS signal by a factor of 2 using an FFT "brick-wall" filter.
#' @param mrs_data MRS data object.
#' @return downsampled data.
#' @export
downsample_mrs_fd <- function(mrs_data) {
  # needs to be a FD operation
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  N <- Npts(mrs_data)
  mrs_data <- crop_spec(mrs_data, xlim = c(N / 2 - N / 4 + 1, N / 2 + N / 4),
                        scale = "points")
  mrs_data
}

#' Downsample an MRS signal by a factor of 2 by removing every other data point
#' in the time-domain. Note, signals outside the new sampling frequency will be
#' aliased.
#' @param mrs_data MRS data object.
#' @return downsampled data.
#' @export
downsample_mrs_td <- function(mrs_data) {
  # needs to be a TD operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  mrs_data <- get_subset(mrs_data, td_set = seq(1, Npts(mrs_data), 2))
  mrs_data$resolution[7] <- mrs_data$resolution[7] * 2
  return(mrs_data)
}

#' Decimate an MRS signal by filtering in the time domain before downsampling.
#' @param mrs_data MRS data object.
#' @param q integer factor to downsample by (default = 2).
#' @param n filter order used in the downsampling.
#' @param ftype	filter type, "iir" or "fir".
#' @return decimated data.
#' @export
decimate_mrs_td <- function(mrs_data, q = 2, n = 4, ftype = "iir") {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, decimate_mrs_td, q = q, n = n, ftype = ftype))
  }
  
  # needs to be a TD operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  mrs_data_re <- apply_mrs(Re(mrs_data), 7, fun = signal::decimate, q, n, ftype)
  mrs_data_im <- apply_mrs(Im(mrs_data), 7, fun = signal::decimate, q, n, ftype)
  mrs_data$data <- (mrs_data_re$data + 1i * mrs_data_im$data) * q
  mrs_data$resolution[7] <- mrs_data$resolution[7] * q
  dimnames(mrs_data$data) <- NULL
  mrs_data
}

#' Decimate an MRS signal to half the original sampling frequency by filtering
#' in the frequency domain before down sampling.
#' @param mrs_data MRS data object.
#' @return decimated data at half the original sampling frequency.
#' @export
decimate_mrs_fd <- function(mrs_data) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, decimate_mrs_fd))
  }
  
  # needs to be a FD operation initially
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  mrs_data <- apply_mrs(mrs_data, 7, smooth_high_freq_vec)
  dimnames(mrs_data$data) <- NULL
  mrs_data <- fd2td(mrs_data)
  mrs_data <- downsample_mrs_td(mrs_data)
  return(mrs_data)
}

smooth_high_freq_vec <- function(vec) {
  vec <- pracma::ifftshift(vec) 
  # smooth the central portion
  N <- length(vec)
  inds <- (N / 4 + 1):(N - N / 4)
  vec_smooth_re <- stats::predict(stats::smooth.spline(Re(vec[inds])),
                                  1:length(inds))$y
  
  vec_smooth_im <- stats::predict(stats::smooth.spline(Im(vec[inds])),
                                  1:length(inds))$y
  
  vec[inds] <- vec_smooth_re + 1i * vec_smooth_im
  vec <- pracma::fftshift(vec) 
  return(vec)
}

# alternate TD method that might cause a bit of phase distortion
# downsample_mrs <- function(mrs_data) {
#  # needs to be a TD operation
#  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
#  mrs_data$data <- mrs_data$data[,,,,,,c(TRUE, FALSE), drop = FALSE] +
#                   mrs_data$data[,,,,,,c(FALSE, TRUE), drop = FALSE]
#  
#  # apply half a data points worth of frequency dep phase correction
#  # TODO
#  
#   mrs_data$resolution <- mrs_data$resolution * 2
#   mrs_data
# }

ift <- function(mrs_data, dims) {
  apply_mrs(mrs_data, dims, ift_shift)
}

dim.mrs_data <- function(x) {
  dim(x$data)
}

#' Return the total number of spectra in an MRS dataset.
#' @param mrs_data MRS data.
#' @export
Nspec <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_dims <- dim(mrs_data$data)
  (mrs_dims[1] * mrs_dims[2] * mrs_dims[3] * mrs_dims[4] * mrs_dims[5] *
   mrs_dims[6])
}

#' Return the total number of x locations in an MRS dataset.
#' @param mrs_data MRS data.
#' @export
Nx <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  dim(mrs_data$data)[2]
}

#' Return the total number of y locations in an MRS dataset.
#' @param mrs_data MRS data.
#' @export
Ny <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  dim(mrs_data$data)[3]
}

#' Return the total number of z locations in an MRS dataset.
#' @param mrs_data MRS data.
#' @export
Nz <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data)
  
  dim(mrs_data$data)[4]
}

#' Return the total number of dynamic scans in an MRS dataset.
#' @param mrs_data MRS data.
#' @export
Ndyns <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data)
  
  dim(mrs_data$data)[5]
}

#' Return the total number of coil elements in an MRS dataset.
#' @param mrs_data MRS data.
#' @export
Ncoils <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data)
  
  dim(mrs_data$data)[6]
}

#' Return the number of data points in an MRS dataset.
#' @param mrs_data MRS data.
#' @return number of data points.
#' @export
Npts <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data)
  
  dim(mrs_data$data)[7]
}

#' Return the sampling frequency in Hz of an MRS dataset.
#' @param mrs_data MRS data.
#' @return sampling frequency in Hz.
#' @export
fs <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data)
  
  1 / mrs_data$resolution[7]
}

#' Return the frequency scale of an MRS dataset in Hz.
#' @param mrs_data MRS data.
#' @param fs sampling frequency in Hz.
#' @param N number of data points in the spectral dimension.
#' @return frequency scale.
#' @export
hz <- function(mrs_data, fs = NULL, N = NULL) {
  
  # check the input
  if (!missing(mrs_data)) check_mrs_data(mrs_data)
  
  if (is.null(fs)) fs <- fs(mrs_data)
  
  if (is.null(N)) N <- Npts(mrs_data)
     
  seq(from = -fs / 2, to = fs / 2 - fs / N, length.out = N)
}

#' Return the ppm scale of an MRS dataset or fit result.
#' @param x MRS dataset of fit result.
#' @param ft transmitter frequency in Hz, does not apply when the object is a
#' fit result.
#' @param ref reference value for ppm scale, does not apply when the object is a
#' fit result.
#' @param fs sampling frequency in Hz, does not apply when the object is a
#' fit result.
#' @param N number of data points in the spectral dimension, does not apply when the object is a
#' fit result.
#' @return ppm scale.
#' @export
ppm <- function(x, ft = NULL, ref = NULL, fs= NULL, N = NULL) UseMethod("ppm")

#' @rdname ppm
#' @export 
ppm.mrs_data <- function(x, ft = NULL, ref = NULL, fs= NULL, N = NULL) {
  
  # check the input
  check_mrs_data(x)
  
  if (is.null(ft)) ft <- x$ft
  
  if (is.null(ref)) ref <- x$ref
  
  if (is.null(fs)) fs <- fs(x)
  
  if (is.null(N)) N <- Npts(x)
  
  -hz(fs = fs, N = N) / ft * 1e6 + ref
}

#' @rdname ppm
#' @export 
ppm.fit_result <- function(x, ft = NULL, ref = NULL, fs= NULL, N = NULL) {
  fit_is_na <- is.na(x$fits)
  if (sum(fit_is_na) > 0) {
    first_non_na <- which(!fit_is_na)[[1]]
    return(x$fits[[first_non_na]]$PPMScale)
  } else {
    return(x$fits[[1]]$PPMScale)
  }
}

n2hz <- function(n, N, fs) {
  -fs / 2 + (fs / N) * (n - 1)
}
  
hz2ppm <- function(hz_in, ft, ref) {
  ref - hz_in / ft * 1e6
}

ppm2hz <- function(ppm_in, ft, ref) {
  (ref - ppm_in) * ft / 1e6
}

pts <- function(mrs_data) {
  seq(from = 1, to = Npts(mrs_data))  
}

#' Return a time scale vector to match the FID of an MRS data object.
#' @param mrs_data MRS data.
#' @return time scale vector in units of seconds.
#' @export
seconds <- function(mrs_data) {
  fs <- fs(mrs_data)
  seq(from = 0, to = (Npts(mrs_data) - 1) / fs, by = 1 / fs)
}

#' Get the indices of data points lying between two values (end > x > start).
#' @param scale full list of values.
#' @param start smallest value in the subset.
#' @param end largest value in the subset.
#' @return set of indices.
#' @export
get_seg_ind <- function(scale, start, end) {
  if (start > end) {
    tmp <- end
    end <- start
    start <- tmp
  }
  which(scale >= start & scale <= end)
}

#' Crop \code{mrs_data} object data points in the time-domain.
#' @param mrs_data MRS data.
#' @param start starting data point (defaults to 1).
#' @param end ending data point (defaults to the last saved point).
#' @return cropped \code{mrs_data} object.
#' @export
crop_td_pts <- function(mrs_data, start = NULL, end = NULL) {
  
  # needs to be a TD operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  if (is.null(start)) start <- 1
  
  if (is.null(end)) end <- Npts(mrs_data)
  
  mrs_data$data <- mrs_data$data[,,,,,,start:end, drop = F]
  
  mrs_data
}

#' Crop \code{mrs_data} object based on a frequency range.
#' @param mrs_data MRS data.
#' @param xlim range of values to crop in the spectral dimension eg 
#' xlim = c(4, 0.2).
#' @param scale the units to use for the frequency scale, can be one of: "ppm", 
#' "hz" or "points".
#' @return cropped \code{mrs_data} object.
#' @export
crop_spec <- function(mrs_data, xlim = c(4, 0.2), scale = "ppm") {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, crop_spec, xlim = xlim, scale = scale))
  }
  
  # needs to be a FD operation
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  if (scale == "ppm") {
    x_scale <- ppm(mrs_data)
  } else if (scale == "hz") {
    x_scale <- hz(mrs_data)
  } else if (scale == "points") {
    x_scale <- pts(mrs_data)
  } else {
    stop("Error, scale not recognised.")
  }
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[Npts(mrs_data)])
  }
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  old_ppm <- ppm(mrs_data)
  #old_ref <- mrs_data$ref
  
  # update fs
  mrs_data$resolution[7] <- (mrs_data$resolution[7] / length(subset) *
                             Npts(mrs_data))
  
  mrs_data$data <- mrs_data$data[,,,,,, subset, drop = F]
  #print(length(subset))
  
  # not sure why subset[2] works better than subset[1]
  new_ppm = (old_ppm[subset[length(subset)]] + old_ppm[subset[2]])/2
  mrs_data$ref <- new_ppm
  
  dimnames(mrs_data$data) <- NULL
    
  mrs_data
}

#' Align spectra to a reference frequency using a convolution based method.
#' @param mrs_data data to be aligned.
#' @param ref_freq reference frequency in ppm units. More than one frequency
#' may be specified.
#' @param zf_factor zero filling factor to increase alignment resolution.
#' @param lb line broadening to apply to the reference signal.
#' @param max_shift maximum allowable shift in Hz.
#' @param ret_df return frequency shifts in addition to aligned data (logical).
#' @param mean_dyns align the mean spectrum and apply the same shift to each
#' dynamic.
#' @return aligned data object.
#' @export
align <- function(mrs_data, ref_freq = 4.65, zf_factor = 2, lb = 2,
                  max_shift = 20, ret_df = FALSE, mean_dyns = FALSE) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, align, ref_freq = ref_freq, zf_factor = zf_factor,
                  lb = lb, max_shift = max_shift, ret_df = ret_df,
                  mean_dyns = mean_dyns))
  }
  
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  if (mean_dyns) {
    mrs_data_orig <- mrs_data
    mrs_data <- mean_dyns(mrs_data)
  }
  
  mrs_data_zf <- zf(mrs_data, zf_factor)
  
  freq <- ppm2hz(ref_freq, mrs_data$ft, mrs_data$ref)
  t_zf <- seconds(mrs_data_zf)
  
  freq_mat <- matrix(freq, length(freq), length(t_zf), byrow = FALSE)
  t_zf_mat <- matrix(t_zf, length(freq), length(t_zf), byrow = TRUE)
  
  ref_data <- colSums(exp(2i * t_zf_mat * pi * freq_mat - lb * t_zf_mat * pi))
  
  window <- floor(max_shift * Npts(mrs_data_zf) * mrs_data$resolution[7])
  
  shifts <- apply_mrs(mrs_data_zf, 7, conv_align, ref_data, window,
                      1 / mrs_data$resolution[7], FALSE, data_only = TRUE)
  
  if (mean_dyns) mrs_data <- mrs_data_orig
  
  t_orig <- rep(seconds(mrs_data), each = Nspec(mrs_data))
  t_array <- array(t_orig, dim = dim(mrs_data$data))
  shift_array <- array(shifts, dim = dim(mrs_data$data))
  shift_array <- exp(2i * pi * t_array * shift_array)
  mrs_data$data <- mrs_data$data * shift_array
  
  if (ret_df) {
    return(list(data = mrs_data, shifts = abind::adrop(shifts, 7)))
  } else {
    return(mrs_data)
  }
}

#' Return an array of amplitudes derived from fitting the initial points in the
#' time domain and extrapolating back to t=0.
#' @param mrs_data MRS data.
#' @param nstart first data point to fit.
#' @param nend last data point to fit.
#' @param method method for measuring the amplitude, one of "spline" or "exp".
#' @return array of amplitudes.
#' @export
get_td_amp <- function(mrs_data, nstart = 10, nend = 50, method = "spline") {
  
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  if (method == "spline") {
    amps <- apply_mrs(mrs_data, 7, measure_td_amp, nstart, nend)$data
  } else if (method == "exp") {
    t <- seconds(mrs_data)
    amps <- apply_mrs(mrs_data, 7, measure_lorentz_amp, t, nstart, nend)$data
  } else {
    stop("incorrect method for get_td_amp")
  }
 
  amps <- abind::adrop(amps, 7)
  amps
}

conv_align <- function(acq, ref, window, fs, fd) {
  if (fd) {
    conv <- pracma::fftshift(Mod(stats::convolve(acq, ref)))
  } else {
    conv <- pracma::fftshift(Mod(convolve_td(acq, ref)))
  }
  conv_crop <- array(conv[(length(acq) / 2 - window + 1):(length(acq) / 2 +
                     window + 1)])
  #plot(conv_crop)
  #print(-which.max(conv_crop))
  shift_pts <- (-which.max(conv_crop) + window + 1)
  shift_hz <- shift_pts * fs / length(acq)
  return(shift_hz)
}

shift_hz <- function(fid_in, shifts, t) {
  return(fid_in * exp(2i * pi * t * shift_hz))
}

#' Extract a subset of dynamic scans.
#' @param mrs_data dynamic MRS data.
#' @param subset vector containing indices to the dynamic scans to be 
#' returned.
#' @return MRS data containing the subset of requested dynamics.
#' @export
get_dyns <- function(mrs_data, subset) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- mrs_data$data[,,,, subset,,, drop = FALSE]
  return(mrs_data)
}

#' Interleave the first and second half of a dynamic series.
#' @param mrs_data dynamic MRS data.
#' @return interleaved data.
#' @export
interleave_dyns <- function(mrs_data) {
  total <- Ndyns(mrs_data)
  fh <- 1:(total / 2)
  sh <- (total / 2 + 1):total
  new_idx <- c(rbind(fh, sh))
  get_dyns(mrs_data, new_idx)
}

set_dyns <- function(mrs_data, subset, mrs_data_in) {
  mrs_data$data[,,,, subset,,] = mrs_data_in$data[,,,, 1,,]
  return(mrs_data)
}

#' Remove a subset of dynamic scans.
#' @param mrs_data dynamic MRS data.
#' @param subset vector containing indices to the dynamic scans to be 
#' removed.
#' @return MRS data without the specified dynamic scans.
#' @export
rm_dyns <- function(mrs_data, subset) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- mrs_data$data[,,,, -subset,,, drop = F]
  mrs_data
}

#' Return a single voxel from a larger mrs dataset.
#' @param mrs_data MRS data.
#' @param x_pos the x index to plot.
#' @param y_pos the y index to plot.
#' @param z_pos the z index to plot.
#' @param dyn the dynamic index to plot.
#' @param coil the coil element number to plot.
#' @return MRS data.
#' @export
get_voxel <- function(mrs_data, x_pos = 1, y_pos = 1, z_pos = 1, dyn = 1,
                      coil = 1) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- mrs_data$data[1, x_pos, y_pos, z_pos, dyn, coil, , 
                                 drop = FALSE]
  
  return(mrs_data)
}

#' Return a single slice from a larger MRSI dataset.
#' @param mrs_data MRSI data.
#' @param z_pos the z index to extract.
#' @return MRS data.
#' @export
get_slice <- function(mrs_data, z_pos) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- mrs_data$data[,,,z_pos,,,, drop = FALSE]
  return(mrs_data)
}

#' Extract a subset of MRS data.
#' @param mrs_data MRS data object.
#' @param x_set x indices to include in the output (default all).
#' @param y_set y indices to include in the output (default all).
#' @param z_set z indices to include in the output (default all). 
#' @param dyn_set dynamic indices to include in the output (default all). 
#' @param coil_set coil indices to include in the output (default all). 
#' @param fd_set frequency domain data indices to include in the output (default
#' all). 
#' @param td_set time-domain indices to include in the output (default all). 
#' @return selected subset of MRS data.
#' @export
get_subset <- function(mrs_data, x_set = NULL, y_set = NULL, z_set = NULL,
                       dyn_set = NULL, coil_set = NULL, fd_set = NULL,
                       td_set = NULL) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  if (!is.null(td_set) && !is.null(fd_set)) {
    stop("td_set OR fd_set should be specified - not both.")
  }
  
  orig_dims <- dim(mrs_data$data)
  
  if (is.null(x_set)) x_set       <- 1:orig_dims[2]
  if (is.null(y_set)) y_set       <- 1:orig_dims[3]
  if (is.null(z_set)) z_set       <- 1:orig_dims[4]
  if (is.null(dyn_set)) dyn_set   <- 1:orig_dims[5]
  if (is.null(coil_set)) coil_set <- 1:orig_dims[6]
  n_set <- 1:orig_dims[7]
  
  if (!is.null(fd_set)) {
    # needs to be a FD operation
    if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
    n_set <- fd_set
  }
  
  if (!is.null(td_set)) {
    # needs to be a TD operation
    if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
    n_set <- td_set
  }
  
  mrs_data$data <- mrs_data$data[, x_set, y_set, z_set, dyn_set, coil_set,
                                 n_set, drop = FALSE]
  
  dimnames(mrs_data$data) <- NULL
  
  return(mrs_data)
}

#' Crop an MRSI dataset in the x-y direction
#' @param mrs_data MRS data object.
#' @param x_dim x dimension output length.
#' @param y_dim y dimension output length.
#' @return selected subset of MRS data.
#' @export
crop_xy <- function(mrs_data, x_dim, y_dim) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mid_pt_x <- Nx(mrs_data) / 2
  mid_pt_y <- Ny(mrs_data) / 2
  x_set <- seq(from = mid_pt_x - x_dim / 2 + 1, by = 1, length.out = x_dim)
  y_set <- seq(from = mid_pt_y - y_dim / 2 + 1, by = 1, length.out = y_dim)
  x_set <- floor(x_set) # could be floor or ceil, need to test
  y_set <- floor(y_set) # could be floor or ceil, need to test
  x_offset <- x_set[1] - 1
  y_offset <- y_set[1] - 1
  
  affine <- mrs_data$affine
  
  affine[1:3, 4] <- affine[1:3, 4] + mrs_data$affine[1:3, 1] * x_offset +
                                     mrs_data$affine[1:3, 2] * y_offset
  
  mrs_data$affine <- affine
  
  return(get_subset(mrs_data, x_set = x_set, y_set = y_set))
}

#' Mask an MRSI dataset in the x-y direction
#' @param mrs_data MRS data object.
#' @param x_dim x dimension output length.
#' @param y_dim y dimension output length.
#' @return masked MRS data.
#' @export
mask_xy <- function(mrs_data, x_dim, y_dim) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mid_pt_x <- Nx(mrs_data) / 2
  mid_pt_y <- Ny(mrs_data) / 2
  x_set <- seq(from = mid_pt_x - x_dim / 2 + 1, by = 1, length.out = x_dim)
  y_set <- seq(from = mid_pt_y - y_dim / 2 + 1, by = 1, length.out = y_dim)
  x_set <- floor(x_set) # could be floor or ceil, need to test
  y_set <- floor(y_set) # could be floor or ceil, need to test
  mask_mat <- matrix(TRUE, Nx(mrs_data), Ny(mrs_data))
  mask_mat[x_set, y_set] <- FALSE
  mrs_data <- mask_xy_mat(mrs_data, mask_mat)
  return(mrs_data)
}

#' Mask a 2D MRSI dataset in the x-y dimension.
#' @param mrs_data MRS data object.
#' @param mask matrix of boolean values specifying the voxels to mask, where a
#' value of TRUE indicates the voxel should be removed.
#' @param value the value to set masked data to (usually NA or 0).
#' @return masked dataset.
#' @export
mask_xy_mat <- function(mrs_data, mask, value = NA) {
  
  # check the input
  check_mrs_data(mrs_data)
  
  dim(mask) <- c(1, nrow(mask), ncol(mask), 1, 1, 1, 1)
  mask <- rep_array_dim(mask, 7, Npts(mrs_data))
  mrs_data$data[mask] <- value
  return(mrs_data)
}

#' Set the masked voxels in a 2D MRSI dataset to given spectrum.
#' @param mrs_data MRSI data object.
#' @param mask matrix of boolean values specifying the voxels to set, where a
#' value of TRUE indicates the voxel should be set to mask_mrs_data.
#' @param mask_mrs_data the spectral data to be assigned to the masked voxels.
#' @return updated dataset.
#' @export
set_mask_xy_mat <- function(mrs_data, mask, mask_mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data)
  check_mrs_data(mask_mrs_data)
  
  if (Nspec(mask_mrs_data) > 1) {
    stop("mask_mrs_data should only contain one spectrum")
  }
  
  if (Npts(mrs_data) != Npts(mask_mrs_data)) {
    stop("mrs_data and mask_mrs_data don't have the same number of data points")
  }
  
  # make sure mask_mrs_data is in the same domain as mrs_data 
  if (is_fd(mrs_data) & !is_fd(mask_mrs_data)) {
    mask_mrs_data <- td2fd(mask_mrs_data)
  }
  
  if (!is_fd(mrs_data) & is_fd(mask_mrs_data)) {
    mask_mrs_data <- fd2td(mask_mrs_data)
  }
  
  voxels <- sum(mask)
  
  dim(mask) <- c(1, nrow(mask), ncol(mask), 1, 1, 1, 1)
  mask <- rep_array_dim(mask, 7, Npts(mrs_data))
  mrs_data$data[mask] <- rep(drop(mask_mrs_data$data), each = voxels)
  return(mrs_data)
}

#' Mask an MRS dataset in the dynamic dimension.
#' @param mrs_data MRS data object.
#' @param mask vector of boolean values specifying the dynamics to mask, where a
#' value of TRUE indicates the spectrum should be removed.
#' @return masked dataset.
#' @export
mask_dyns <- function(mrs_data, mask) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  dim(mask) <- c(1, 1, 1, 1, length(mask), 1, 1)
  mask <- rep_array_dim(mask, 7, Npts(mrs_data))
  mrs_data$data[mask] <- NA
  return(mrs_data)
}

#' Return the first scans of a dynamic series.
#' @param mrs_data dynamic MRS data.
#' @param n the number of dynamic scans to return.
#' @return first scans of a dynamic series.
#' @export
get_head_dyns <- function(mrs_data, n = 1) {
  index <- 1:n
  get_dyns(mrs_data, index)
}

#' Return the last scans of a dynamic series.
#' @param mrs_data dynamic MRS data.
#' @param n the number of dynamic scans to return.
#' @return last scans of a dynamic series.
#' @export
get_tail_dyns <- function(mrs_data, n = 1) {
  index <- (Ndyns(mrs_data) - n + 1):Ndyns(mrs_data)
  get_dyns(mrs_data, index)
}

#' Return the first half of a dynamic series.
#' @param mrs_data dynamic MRS data.
#' @return first half of the dynamic series.
#' @export
get_fh_dyns <- function(mrs_data) {
  fh <- 1:(Ndyns(mrs_data) / 2)
  get_dyns(mrs_data, fh)
}

#' Return the second half of a dynamic series.
#' @param mrs_data dynamic MRS data.
#' @return second half of the dynamic series.
#' @export
get_sh_dyns <- function(mrs_data) {
  sh <- (Ndyns(mrs_data) / 2 + 1):Ndyns(mrs_data)
  get_dyns(mrs_data, sh)
}
  
#' Return odd numbered dynamic scans starting from 1 (1,3,5...).
#' @param mrs_data dynamic MRS data.
#' @return dynamic MRS data containing odd numbered scans.
#' @export
get_odd_dyns <- function(mrs_data) {
  subset <- seq(1, Ndyns(mrs_data), 2)
  get_dyns(mrs_data, subset)
}

#' Return even numbered dynamic scans starting from 1 (2,4,6...).
#' @param mrs_data dynamic MRS data.
#' @return dynamic MRS data containing even numbered scans.
#' @export
get_even_dyns <- function(mrs_data) {
  subset <- seq(2, Ndyns(mrs_data), 2)
  get_dyns(mrs_data, subset)
}

#' Invert odd numbered dynamic scans starting from 1 (1,3,5...).
#' @param mrs_data dynamic MRS data.
#' @return dynamic MRS data with inverted odd numbered scans.
#' @export
inv_odd_dyns <- function(mrs_data) {
  subset <- seq(1, Ndyns(mrs_data), 2)
  mrs_data$data[,,,, subset,,] <- -1 * mrs_data$data[,,,, subset,,]
  return(mrs_data)
}

#' Invert even numbered dynamic scans starting from 1 (2,4,6...).
#' @param mrs_data dynamic MRS data.
#' @return dynamic MRS data with inverted even numbered scans.
#' @export
inv_even_dyns <- function(mrs_data) {
  subset <- seq(2, Ndyns(mrs_data), 2)
  mrs_data$data[,,,, subset,,] <- -1 * mrs_data$data[,,,, subset,,]
  return(mrs_data)
}

#' Combine a reference and metabolite mrs_data object.
#' @param metab metabolite mrs_data object.
#' @param ref reference mrs_data object.
#' @return combined metabolite and reference mrs_data object.
#' @export
comb_metab_ref <- function(metab, ref) {
  
  # check the input
  check_mrs_data(metab) 
  check_mrs_data(ref) 
  
  metab$data <- abind::abind(metab$data, ref$data, along = 1)
  metab
}

#' Extract the reference component from an mrs_data object.
#' @param mrs_data MRS data.
#' @return reference component.
#' @export
get_ref <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- mrs_data$data[2,,,,,,,drop = FALSE]
  mrs_data
}

#' Extract the metabolite component from an mrs_data object.
#' @param mrs_data MRS data.
#' @return metabolite component.
#' @export
get_metab <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- mrs_data$data[1,,,,,,,drop = FALSE]
  mrs_data
}

#' Append MRS data across the coil dimension, assumes they matched across the
#' other dimensions.
#' @param ... MRS data objects as arguments, or a list of MRS data objects.
#' @return a single MRS data object with the input objects concatenated together.
#' @export
append_coils <- function(...) {
  x <- list(...)
  
  # were the arguments a list already? 
  if (depth(x) == 3) x <- x[[1]]
  
  first_dataset <- x[[1]]
  
  # data needs to be in the same domain
  if (is_fd(first_dataset)) {
    for (n in 1:length(x)) {
      if (!is_fd(x[[n]])) {
        x[[n]] <- td2fd(x[[n]])
      }
      x[[n]] <- x[[n]]$data
    }
  } else {
    for (n in 1:length(x)) {
      if (is_fd(x[[n]])) {
        x[[n]] <- fd2td(x[[n]])
      }
      x[[n]] <- x[[n]]$data
    }
  }
  
  new_data <- abind::abind(x, along = 6)
  first_dataset$data <- unname(new_data)
  first_dataset
}

#' Append MRS data across the dynamic dimension, assumes they matched across the
#' other dimensions.
#' @param ... MRS data objects as arguments, or a list of MRS data objects.
#' @return a single MRS data object with the input objects concatenated together.
#' @export
append_dyns <- function(...) {
  x <- list(...)
  
  # were the arguments a list already? 
  if (class(x[[1]])[1] == "list") x <- x[[1]]
  
  # were the arguments a list already? 
  # this doesn't work now we have the meta field :(
  # if (depth(x) == 3) x <- x[[1]]
  
  first_dataset <- x[[1]]
  
  # data needs to be in the same domain
  if (is_fd(first_dataset)) {
    for (n in 1:length(x)) {
      if (!is_fd(x[[n]])) {
        x[[n]] <- td2fd(x[[n]])
      }
      x[[n]] <- x[[n]]$data
    }
  } else {
    for (n in 1:length(x)) {
      if (is_fd(x[[n]])) {
        x[[n]] <- fd2td(x[[n]])
      }
      x[[n]] <- x[[n]]$data
    }
  }
  
  new_data <- abind::abind(x, along = 5)
  first_dataset$data <- unname(new_data)
  first_dataset
}

append_scan <- function(...) {
  x <- list(...)
  
  # were the arguments a list already? 
  if (class(x[[1]])[1] == "list") x <- x[[1]]
  
  # this doesn't work now we have the meta field :(
  # if (depth(x) == 3) x <- x[[1]]
  
  first_dataset <- x[[1]]
  
  if (is_fd(first_dataset)) {
    first_dataset <- fd2td(first_dataset)
  }
  
  dim_first_dataset <- dim(first_dataset$data)
  
  # data needs to be in the same domain
  for (n in 1:length(x)) {
    if (!identical(dim(x[[n]]$data), dim_first_dataset)) {
      stop("Data dimensions passed to stackplot do not match.")
    }
    if (is_fd(x[[n]])) {
        x[[n]] <- fd2td(x[[n]])
    }
    x[[n]] <- x[[n]]$data
  }
  
  new_data <- abind::abind(x, along = 1)
  first_dataset$data <- unname(new_data)
  first_dataset
}

split_metab_ref <- function(mrs_data) {
  metab <- mrs_data
  ref <- mrs_data
  metab$data = metab$data[1,,,,,,, drop = FALSE]
  ref$data = ref$data[2,,,,,,, drop = FALSE]
  return(list(metab, ref))
}

bc <- function(mrs_data, lambda = 1e3, p = 0.1) {
  if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
  }
  # extract real part
  mrs_data$data <- Re(mrs_data$data)
  apply_mrs(mrs_data, 7, ptw::baseline.corr, lambda, p)
}

#' Sum two mrs_data objects.
#' @param a first mrs_data object to be summed.
#' @param b second mrs_data object to be summed.
#' @param force set to TRUE to force mrs_data objects to be summed, even if they
#' are in different time/frequency domains.
#' @return a + b 
#' @export
sum_mrs <- function(a, b, force = FALSE) {
 
  if (class(a) != "mrs_data") stop("a argument is not an mrs_data object")
  
  if (class(b) != "mrs_data") stop("b argument is not an mrs_data object")
  
  # if they are not in the same domain and we are not forcing, then swap b to 
  # match a
  if ((is_fd(a) != is_fd(b)) & !force) {
    if (is_fd(b)) {
      b <- fd2td(b)
    } else {
      b <- td2fd(b)
    }
  }
  
  a$data <- a$data + b$data
  
  return(a)
}

#' Perform a mathematical operation on a spectral region.
#' @param mrs_data MRS data.
#' @param xlim spectral range to be integrated (defaults to full range).
#' @param operator can be "sum" (default), "mean", "l2", "max", "min" or
#' "max-min".
#' @param freq_scale units of xlim, can be : "ppm", "hz" or "points".
#' @param mode spectral mode, can be : "re", "im", "mod" or "cplx".
#' @return an array of integral values.
#' @export
spec_op <- function(mrs_data, xlim = NULL, operator = "sum", freq_scale = "ppm",
                    mode = "re") {
  
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
    
  if (freq_scale == "ppm") {
    x_scale <- ppm(mrs_data)
  } else if (freq_scale == "hz") {
    x_scale <- hz(mrs_data)
  } else if (freq_scale == "points") {
    x_scale <- pts(mrs_data)
  } else {
    stop("unknown freq_scale for spec_op function")
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
  } else if (mode == "cplx") {
    data_arr <- data_arr
  } else {
    stop("unknown mode for spec_op function")
  }
 
  if (operator == "l2") {
    data_arr <- data_arr * data_arr
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), sum)
    res <- res ^ 0.5
  } else if (operator == "mean") {
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), mean)
  } else if (operator == "sum") {
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), sum)
  } else if (operator == "max") {
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), max)
  } else if (operator == "min") {
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), min)
  } else if (operator == "max-min") {
    res <- apply(data_arr, c(1, 2, 3, 4, 5, 6), max) - 
           apply(data_arr, c(1, 2, 3, 4, 5, 6), min)
  } else {
    stop("unknown operator argument for spec_op function") 
  }
  
  return(res) 
}

#' Integrate a spectral region.
#' 
#' See spec_op function for a more complete set of spectral operations.
#' 
#' @param mrs_data MRS data.
#' @param xlim spectral range to be integrated (defaults to full range).
#' @param freq_scale units of xlim, can be : "ppm", "hz" or "points".
#' @param mode spectral mode, can be : "re", "im", "mod" or "cplx".
#' @return an array of integral values.
#' @export
int_spec <- function(mrs_data, xlim = NULL, freq_scale = "ppm", mode = "re") {
  
  res <- spec_op(mrs_data = mrs_data, xlim = xlim, operator = "sum",
                 freq_scale = freq_scale, mode = mode)
  
  return(res) 
}

#' Scale an mrs_data object by a scalar or vector or amplitudes.
#' @param mrs_data data to be scaled.
#' @param amp multiplicative factor, must have length equal to 1 or
#' Nspec(mrs_data).
#' @return mrs_data object multiplied by the amplitude scale factor.
#' @export
scale_mrs_amp <- function(mrs_data, amp) {
  
  if (class(mrs_data) == "list") {
    res <- lapply(mrs_data, scale_mrs_amp, amp = amp)
    return(res)
  }
 
  if (class(mrs_data) != "mrs_data") {
    stop("first argument is not an mrs_data object")
  }
  
  if (length(amp) == 1) {
    mrs_data$data <- mrs_data$data * amp
  } else if (length(amp) == Nspec(mrs_data)) {
    amps_full <- array(rep(amp, Npts(mrs_data)), dim = dim(mrs_data$data))
    mrs_data$data <- mrs_data$data * amps_full
  } else {
    stop("scale should either have length equal to 1 or Nspec(mrs_data)")
  }
  
  return(mrs_data)
}

#' Scale mrs_data to a spectral region.
#' @param mrs_data MRS data.
#' @param xlim spectral range to be integrated (defaults to full range).
#' @param operator can be "sum" (default), "mean", "l2", "max", "min" or
#' "max-min".
#' @param freq_scale units of xlim, can be : "ppm", "Hz" or "points".
#' @param mode spectral mode, can be : "re", "im", "mod" or "cplx".
#' @param mean_dyns mean the dynamic scans before applying the operator. The
#' same scaling value will be applied to each individual dynamic.
#' @return normalised data.
#' @export
scale_spec <- function(mrs_data, xlim = NULL, operator = "sum",
                       freq_scale = "ppm", mode = "re", mean_dyns = TRUE) {
  
  if (class(mrs_data) == "list") {
    res <- lapply(mrs_data, scale_spec, xlim = xlim, operator = operator,
                  freq_scale = freq_scale, mode = mode, mean_dyns = mean_dyns)
    return(res)
  }
  
  if (mean_dyns) {
    amp <- spec_op(mean_dyns(mrs_data), xlim, operator, freq_scale, mode)
    amp <- as.numeric(amp)
  } else {
    amp <- spec_op(mrs_data, xlim, operator, freq_scale, mode)
  }
  
  mrs_data <- scale_mrs_amp(mrs_data, 1 / amp)
  
  return(mrs_data)
}

#' @export
`+.mrs_data` <- function(a, b) {
  if (class(b) == "mrs_data" ) {
    if (is_fd(a) != is_fd(b)) {
      warning("sum warning, spectral domains do not match")
    }
    a$data <- a$data + b$data
  } else if (class(b) %in% c("numeric", "integer", "complex")) {
    a$data <- a$data + b
  } else {
    stop("unrecognised type")
  }
  return(a)
}

#' @export
`-.mrs_data` <- function(a, b = NULL) {
  if (class(b) == "mrs_data" ) {
    if (is_fd(a) != is_fd(b)) {
      warning("subtract warning, time/frequency domains do not match")
    }
    a$data = a$data - b$data
  } else if (is.null(b)) {
    a$data = -a$data
  } else if (class(b) %in% c("numeric", "integer", "complex")) {
    a$data = a$data - b
  } else {
    stop("unrecognised type")
  }
  return(a)
}

#' @export
`*.mrs_data` <- function(a, b) {
  if (class(b) == "mrs_data" ) {
    if (is_fd(a) != is_fd(b)) {
      warning("multiply warning, time/frequency domains do not match")
    }
    a$data <- a$data * b$data
  } else if (class(b) %in% c("numeric", "integer", "complex")) {
    a$data <- a$data * b
  } else {
    stop("unrecognised type")
  }
  return(a)
}

#' @export
`/.mrs_data` <- function(a, b) {
  if (class(b) == "mrs_data" ) {
    if (is_fd(a) != is_fd(b)) {
      warning("divide warning, time/frequency domains do not match")
    }
    a$data <- a$data / b$data
  } else if (class(b) %in% c("numeric", "integer", "complex")) {
    a$data <- a$data / b
  } else {
    stop("unrecognised type")
  }
  return(a)
}

#' Calculate the mean spectrum from an mrs_data object.
#' @param x object of class mrs_data.
#' @param ... other arguments to pass to the colMeans function.
#' @return mean mrs_data object.
#' @export
mean.mrs_data <- function(x, ...) {
  data_pts <- x$data
  data_N <- Npts(x)
  dim(data_pts) <- c(length(data_pts) / data_N, data_N)
  x$data <- colMeans(data_pts, ...)
  dim(x$data) <- c(1, 1, 1, 1, 1, 1, data_N)
  x
}

#' Calculate the standard deviation spectrum from an mrs_data object.
#' @param x object of class mrs_data.
#' @param na.rm remove NA values.
#' @return sd mrs_data object.
#' @export
sd.mrs_data <- function(x, na.rm = FALSE) {
  data_pts <- x$data
  data_N <- Npts(x)
  dim(data_pts) <- c(length(data_pts) / data_N, data_N)
  x$data <- colSdColMeans(data_pts, na.rm)
  dim(x$data) <- c(1, 1, 1, 1, 1, 1, data_N)
  x
}

## make an S3 generic for sd (cos R Core don't do this for some reason!)
## see https://cran.r-project.org/doc/manuals/R-exts.html#Adding-new-generics

#' Calculate the standard deviation spectrum from an mrs_data object.
#' @param x object of class mrs_data.
#' @param na.rm remove NA values.
#' @return sd mrs_data object.
#' @export
sd <- function(x, na.rm) UseMethod("sd")

## take the usual definition of sd,
## and set it to be the default method
#' @export
sd.default <- function(x, na.rm = FALSE) stats::sd(x, na.rm)

#' Collapse MRS data by concatenating spectra along the dynamic dimension.
#' @param x data object to be collapsed (mrs_data or fit_result object).
#' @param rm_masked remove masked dynamics from the output.
#' @return collapsed data with spectra or fits concatenated along the dynamic
#' dimension.
#' @rdname collapse_to_dyns
#' @export
collapse_to_dyns <- function(x, rm_masked = FALSE) UseMethod("collapse_to_dyns")

#' @rdname collapse_to_dyns
#' @export
collapse_to_dyns.mrs_data <- function(x, rm_masked = FALSE) {
  data_pts <- x$data
  data_N <- Npts(x)
  dim(data_pts) <- c(1, 1, 1, 1, length(data_pts) / data_N, 1, data_N)
  x$data <- data_pts
  
  if (rm_masked) {
    keepers <- !is.na(x$data[,,,,,,1])
    x$data <- x$data[,,,,keepers,,,drop = FALSE]
  }
  
  return(x)
}

#' @rdname collapse_to_dyns
#' @export
collapse_to_dyns.fit_result <- function(x, rm_masked = FALSE) {
  x$res_tab[c(1, 2, 3, 5)] <- 1
  dyns <- nrow(x$res_tab)
  x$res_tab[4] <- 1:dyns
  x
}

#' Calculate the mean dynamic data.
#' @param mrs_data dynamic MRS data.
#' @return mean dynamic data.
#' @export
mean_dyns <- function(mrs_data) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, mean_dyns))
  }
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- aperm(mrs_data$data, c(5, 1, 2, 3, 4, 6, 7))
  mrs_data$data <- colMeans(mrs_data$data, na.rm = TRUE)
  new_dim <- dim(mrs_data$data)
  dim(mrs_data$data) <- c(new_dim[1:4], 1, new_dim[5:6])
  mrs_data
}

#' Subtract the mean dynamic spectrum from a dynamic series.
#' @param mrs_data dynamic MRS data.
#' @return subtracted data.
#' @export
sub_mean_dyns <- function(mrs_data) {
  mean_mrs_data <- mean_dyns(mrs_data)
  mrs_data_mean_sub <- mrs_data - rep_dyn(mean_mrs_data, Ndyns(mrs_data))
  return(mrs_data_mean_sub)
}

#' Calculate the mean of adjacent dynamic scans.
#' @param mrs_data dynamic MRS data.
#' @param block_size number of adjacent dynamics scans to average over.
#' @return dynamic data averaged in blocks.
#' @export
mean_dyn_blocks <- function(mrs_data, block_size) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, mean_dyn_blocks, block_size = block_size))
  }
  
  if ((Ndyns(mrs_data) %% block_size) != 0) {
    warning("Block size does not fit into the number of dynamics without truncation.")
  }
  
  new_dyns <-  floor(Ndyns(mrs_data) / block_size)
  mrs_out <- get_dyns(mrs_data, seq(1, new_dyns * block_size, block_size))
  for (n in 2:block_size) {
    mrs_out <- mrs_out + get_dyns(mrs_data, seq(n, new_dyns * block_size, block_size))
  }
  
  mrs_out / block_size
}

#' Calculate the pairwise means across a dynamic data set.
#' @param mrs_data dynamic MRS data.
#' @return mean dynamic data of adjacent dynamic pairs.
#' @export
mean_dyn_pairs <- function(mrs_data) {
  pairs <- get_odd_dyns(mrs_data) + get_even_dyns(mrs_data)
  pairs / 2
}

#' Calculate the sum of data dynamics.
#' @param mrs_data dynamic MRS data.
#' @return sum of data dynamics.
#' @export
sum_dyns <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- aperm(mrs_data$data, c(5,1,2,3,4,6,7))
  mrs_data$data <- colSums(mrs_data$data, na.rm = TRUE)
  new_dim <- dim(mrs_data$data)
  dim(mrs_data$data) <- c(new_dim[1:4],1,new_dim[5:6])
  mrs_data
}

#' Calculate the sum across receiver coil elements.
#' @param mrs_data MRS data split across receiver coil elements.
#' @return sum across coil elements.
#' @export
sum_coils <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- aperm(mrs_data$data, c(6,1,2,3,4,5,7))
  mrs_data$data <- colSums(mrs_data$data)
  new_dim <- dim(mrs_data$data)
  dim(mrs_data$data) <- c(new_dim[1:5],1,new_dim[6])
  mrs_data
}

cplx_median <- function(input) {
  stats::median(Re(input), na.rm = TRUE) +
    stats::median(Im(input), na.rm = TRUE) * 1i
}

#' Calculate the median dynamic data.
#' @param mrs_data dynamic MRS data.
#' @return median dynamic data.
#' @export
median_dyns <- function(mrs_data) {
  return(apply_mrs(mrs_data, 5, cplx_median))
}

# TODO correct first imaginary data point?
recon_imag_vec <- function(data) {
  data <- Conj(hilbert(Re(data)))
  data <- ift_shift(data)
  fh <- data[1:(length(data) / 2)]
  sh <- c(data[(length(data) / 2 + 2):length(data)], 0)
  data <- fh + Conj(rev(sh))
}

conv_filt_vec <- function(fid, K = 25, ext = 1) {
  k = -K:K
  filt_fun = exp(-4 * k ^ 2 / K ^ 2)
  filt_fun = filt_fun / sum(filt_fun)
  N = 2 * K + 1
  # filter this signal
  filt_sig <- (stats::filter(Re(fid), filt_fun) +
              1i * stats::filter(Im(fid), filt_fun))
  
  # extrapolate at the left edge
  # real
  x1 <- (N - 1) / 2 + 1
  x2 <- x1 + ext
  y1 <- Re(filt_sig)[x1]
  y2 <- Re(filt_sig)[x2]
  m <- (y2 - y1) / (x2 - x1)
  c = y1 - m * x1
  st <- 1
  end <- (N - 1) / 2
  filt_sig[st:end] <- m * (st:end) + c
  # imag
  y1 <- Im(filt_sig)[x1]
  y2 <- Im(filt_sig)[x2]
  m <- (y2 - y1) / (x2 - x1)
  c = y1 - m * x1
  st <- 1
  end <- (N - 1) / 2
  filt_sig[st:end] <- filt_sig[st:end] + 1i * (m * (st:end) + c)
  
  # extrapolate at the right edge
  # real
  x1 <- length(fid) - ((N - 1) / 2 + 1)
  x2 <- x1 - ext
  y1 <- Re(filt_sig)[x1]
  y2 <- Re(filt_sig)[x2]
  m <- (y2 - y1) / (x2 - x1)
  c = y1 - m * x1
  st <- length(fid) - (N - 1) / 2
  end <- length(fid)
  filt_sig[st:end] <- m*(st:end) + c
  # imag
  y1 <- Im(filt_sig)[x1]
  y2 <- Im(filt_sig)[x2]
  m <- (y2 - y1) / (x2 - x1)
  c = y1 - m * x1
  st <- length(fid) - (N - 1) / 2
  end <- length(fid)
  filt_sig[st:end] <- filt_sig[st:end] + 1i * (m * (st:end) + c)
  
  # subtract from the data and ret
  array(fid - filt_sig)
}

#' Time-domain convolution based filter.
#' 
#' Time-domain convolution based filter described by:
#' Marion D, Ikura M, Bax A. Improved solvent suppression in one-dimensional and
#' twodimensional NMR spectra by convolution of time-domain data. J Magn Reson 
#' 1989;84:425-430.
#' 
#' @param mrs_data MRS data to be filtered.
#' @param K window width in data points.
#' @param ext point separation for linear extrapolation.
#' @export
td_conv_filt <- function(mrs_data, K = 25, ext = 1) {
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  apply_mrs(mrs_data, 7, conv_filt_vec, K, ext)
}

#' Frequency-domain convolution based filter.
#' @param mrs_data MRS data to be filtered.
#' @param K window width in data points.
#' @param ext point separation for linear extrapolation.
#' @export
fd_conv_filt <- function(mrs_data, K = 25, ext = 1) {
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  apply_mrs(mrs_data, 7, conv_filt_vec, K, ext)
}

#' HSVD based signal filter.
#' 
#' HSVD based signal filter described in:
#' Barkhuijsen H, de Beer R, van Ormondt D. Improved algorithm for noniterative 
#' and timedomain model fitting to exponentially damped magnetic resonance
#' signals. J Magn Reson 1987;73:553-557.
#' 
#' @param mrs_data MRS data to be filtered.
#' @param xlim frequency range to filter, default units are Hz which can be
#' changed to ppm using the "scale" argument.
#' @param comps number of Lorentzian components to use for modelling.
#' @param irlba option to use irlba SVD (logical).
#' @param max_damp maximum allowable damping factor.
#' @param scale either "hz" or "ppm" to set the frequency units of xlim.
#' @param return_model by default the filtered spectrum is returned. Set
#' return_model to TRUE to return the HSVD model of the data.
#' @return filtered data or model depending on the return_model argument.
#' @export
hsvd_filt <- function(mrs_data, xlim = c(-30, 30), comps = 40, irlba = TRUE,
                      max_damp = 10, scale = "hz", return_model = FALSE) {
  
  if (class(mrs_data) == "list") {
    res <- lapply(mrs_data, hsvd_filt, xlim = xlim, comps = comps,
                  irlba = irlba, max_damp = max_damp, scale = scale,
                  return_model = return_model)
    return(res)
  }
  
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  if ( scale == "ppm" ) {
    xlim <- ppm2hz(xlim, mrs_data$ft, mrs_data$ref) 
  } else if (scale != "hz") {
    stop("Invalid scale option, should be 'ppm' or 'hz'.")
  }
  
  xlim <- sort(xlim)
  
  filt <- apply_mrs(mrs_data, 7, hsvd_filt_vec, fs = fs(mrs_data),
                    region = xlim, comps = comps, irlba, max_damp = max_damp)
  
  if (return_model) {
    model <- sum_mrs(mrs_data, -filt)
    return(model)
  } else {
    return(filt)
  }
}

hsvd_filt_vec <- function(fid, fs, region = c(-30, 30), comps = 40, 
                          irlba = TRUE, max_damp = 10) {
  
  hsvd_res <- hsvd_vec(fid, fs, comps = comps, irlba, max_damp)
  if (is.null(region)) {
    model <- rowSums(hsvd_res$basis)
  } else {
    idx <- (hsvd_res$reson_table$frequency_hz < region[2]) &
           (hsvd_res$reson_table$frequency_hz > region[1] )
    model <- rowSums(hsvd_res$basis[,idx])
  }
  fid - model
}

#' HSVD of an mrs_data object.
#' 
#' HSVD method as described in:
#' Barkhuijsen H, de Beer R, van Ormondt D. Improved algorithm for noniterative 
#' and timedomain model fitting to exponentially damped magnetic resonance
#' signals. J Magn Reson 1987;73:553-557.
#' 
#' @param mrs_data mrs_data object to be decomposed.
#' @param comps number of Lorentzian components to use for modelling.
#' @param irlba option to use irlba SVD (logical).
#' @param max_damp maximum allowable damping factor.
#' @return basis matrix and signal table.
#' @export
hsvd <- function(mrs_data, comps = 40, irlba = TRUE, max_damp = 10) {
  
  if (Nspec(mrs_data) > 1) {
    warning("Data contains multiple spectra, but only one has been processed.")
  }
  
  res <- hsvd_vec(mrs_data2vec(mrs_data), fs = fs(mrs_data), comps = comps,
                  irlba = irlba, max_damp = max_damp)
  
  res$reson_table$frequency_ppm <- hz2ppm(res$reson_table$frequency,
                                          mrs_data$ft, mrs_data$ref)
  
  res$reson_table$lw_hz <- -res$reson_table$damping / pi
  
  res$basis <- mat2mrs_data(t(res$basis), fs = fs(mrs_data),
                            ref = mrs_data$ref, nuc = mrs_data$nuc,
                            ft = mrs_data$ft, fd = FALSE)
  
  res$model <- sum_dyns(res$basis)
  
  return(res)
}

#' HSVD of a complex vector.
#' 
#' HSVD method as described in:
#' Barkhuijsen H, de Beer R, van Ormondt D. Improved algorithm for noniterative 
#' and timedomain model fitting to exponentially damped magnetic resonance
#' signals. J Magn Reson 1987;73:553-557.
#' 
#' @param y time domain signal to be filtered as a vector.
#' @param fs sampling frequency of y.
#' @param comps number of Lorentzian components to use for modelling.
#' @param irlba option to use irlba SVD (logical).
#' @param max_damp maximum allowable damping factor. Default value of 0 ensures
#' resultant model is damped.
#' @return basis matrix and signal table.
#' @export
hsvd_vec <- function(y, fs, comps = 40, irlba = TRUE, max_damp = 0) {
  N <- length(y)
  L <- floor(0.5 * N)
  
  # scale the input vector to keep things stable
  sc_factor <- max(Mod(y))
  y <- y / sc_factor

  # H is the LxM Hankel LP matrix
  H <- hankel.matrix(L + 1, y)
  H <- H[1:L,]
  
  if (irlba)  {
    svd_res <- irlba::irlba(H, comps)
  } else {
    svd_res <- svd(H)
  }
  
  # construct H of rank K
  Uk <- svd_res$u[,1:comps]
  rows <- nrow(Uk)
  Ukt <- Uk[2:rows,]
  Ukb <- Uk[1:(rows - 1),]
  Zp  <- ginv(Ukb) %*% Ukt

  # find the poles
  q <- pracma::eig(Zp)
  q <- log(q)
  dt <- 1 / fs
  dampings <- Re(q) / dt
  frequencies <- Im(q)/(2 * pi) / dt
  
  # large +ve dampings can cause stability issues where the basis signals
  # can have values 1e88 at the end of the FID causing ginv to fail
  # cap these positive dampings to max_damp
  if (!is.null(max_damp)) dampings[dampings > max_damp] <- max_damp
  
  t <- seq(from = 0, to = (N - 1) / fs, by = 1 / fs)
  t_mat <- matrix(t, ncol = comps, nrow = N)
  
  # TODO not sure if the next line is faster
  #basis <- t(exp(t(t_mat) * (dampings + 2i * pi * frequencies)))
  
  freq_damp <- matrix(dampings + 2i * pi * frequencies, ncol = comps, nrow = N,
                      byrow = TRUE)
  
  basis <- exp(t_mat * freq_damp)
  
  ahat <- ginv(basis) %*% y
  
  # Undo scaling
  ahat <- ahat * sc_factor
  # yhat <- basis %*% ahat 
  
  # scale basis by ahat
  ahat_mat <- matrix(ahat, ncol = comps, nrow = N, byrow = TRUE)
  basis <- basis * ahat_mat
  
  # generate a table of resonances
  reson_table <- data.frame(amplitude = Mod(ahat), phase = Arg(ahat) * 180 / pi,
                            frequency_hz = frequencies, damping = dampings)
  
  list(basis = basis, reson_table = reson_table)
}

#' Perform zeroth-order phase correction based on the minimisation of the
#' squared difference between the real and magnitude components of the
#' spectrum.
#' @param mrs_data an object of class \code{mrs_data}.
#' @param xlim frequency range (default units of PPM) to including in the phase.
#' @param ret_phase return phase values (logical).
#' @return MRS data object and phase values (optional).
#' @export
auto_phase <- function(mrs_data, xlim = NULL, ret_phase = FALSE) {
  
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  mrs_data_proc <- mrs_data
  
  if (!is.null(xlim)) mrs_data_proc <- crop_spec(mrs_data_proc, xlim)
    
  phases <- apply_mrs(mrs_data_proc, 7, auto_phase_vec, data_only = TRUE)
  
  if (length(phases) == 1) phases <- as.numeric(phases)
  
  # TODO update phase function and remove drop
  mrs_data <- phase(mrs_data, phases)
  
  if (ret_phase) {
    return(list(mrs_data = mrs_data, phase = abind::adrop(phases, 7)))
  } else {
    return(mrs_data)
  }
}

auto_phase_vec <- function(vec) {
  res <- stats::optim(0, phase_obj_fn, gr = NULL, vec, method = "Brent",
                      lower = -180, upper = 180)
  #vec * exp(1i * res$par / 180 * pi)
  res$par
}

phase_obj_fn <- function(phi, vec) {
  vec_adj <- vec * exp(1i * phi / 180 * pi)
  sum((Mod(vec_adj) - Re(vec_adj)) ^ 2)
}

ecc_2d_array <- function(array) {
  array <- drop(array)
  metab <- array[1,]
  ref <- array[2,]
  ref_phase <- Arg(ref)
  metab_ecc <- metab * exp(-1i * ref_phase)
  #ref_ecc <- ref * exp(-1i*ref_phase)
  aperm(abind::abind(metab_ecc, ref, along = 2), c(2,1))
}

ecc_1d_array <- function(array) {
  ref <- drop(array)
  ref_phase <- Arg(ref)
  ref * exp(-1i * ref_phase)
}

ecc_ref <- function(mrs_data) {
  if (is_fd(mrs_data)) {
      mrs_data <- fd2td(mrs_data)
  }
  apply_mrs(mrs_data, 7, ecc_1d_array)
}

#' Eddy current correction.
#' 
#' Apply eddy current correction using the Klose method.
#' 
#' In vivo proton spectroscopy in presence of eddy currents.
#' Klose U.
#' Magn Reson Med. 1990 Apr;14(1):26-30.
#' 
#' @param metab MRS data to be corrected.
#' @param ref reference dataset.
#' @param rev reverse the correction.
#' @return corrected data in the time domain.
#' @export
ecc <- function(metab, ref, rev = FALSE) {
  
  if (class(metab) == "list") {
    if (class(ref) != "list") stop("metab is a list but ref is not")
    if (length(metab) != length(ref)) {
      stop("metab and ref must have the same length")
    }
    return(mapply(ecc, metab = metab, ref = ref, MoreArgs = list(rev = rev),
                  SIMPLIFY = FALSE))
  } 
  
  if (is_fd(metab)) metab <- fd2td(metab)
  if (is_fd(ref))   ref <- fd2td(ref)
  
  if (rev) ref <- Conj(ref)
  
  if (Ndyns(ref) > 1) {
    ref <- mean_dyns(ref)
    warning("Using the mean reference signal for ECC.")
  }
  
  # repeat the refernce signal to match the number of dynamics
  if (Ndyns(metab) > 1) {
    ref <- rep_dyn(ref, Ndyns(metab))
  }
  
  mrs_data <- comb_metab_ref(metab, ref)
  ecc_data <- apply_mrs(mrs_data, c(1,7), ecc_2d_array)
  get_metab(ecc_data)
}

#' Apodise MRSI data in the x-y direction with a k-space filter.
#' @param mrs_data MRSI data.
#' @param func must be "hamming" or "gaussian".
#' @param w the reciprocal of the standard deviation for the Gaussian function.
#' @return apodised data.
#' @export
apodise_xy <- function(mrs_data, func = "hamming", w = 2.5) {
  
  # check the input
  check_mrs_data(mrs_data)
  
  x_dim <- Nx(mrs_data)
  y_dim <- Ny(mrs_data)
  z_dim <- Nz(mrs_data)
  dyns  <- Ndyns(mrs_data)
  coils <- Ncoils(mrs_data)
  N     <- Npts(mrs_data)
  
  # transform to k-space
  mrs_data <- img2kspace_xy(mrs_data)
  
  # create the k-space weighting matrix
  weight_mat <- matrix(1, nrow = x_dim, ncol = y_dim)
  
  if (func == "hamming") {
    x_fun <- signal::hamming(x_dim)
    y_fun <- signal::hamming(y_dim)
  } else if (func == "gaussian") {
    x_fun <- signal::gausswin(x_dim, w = w)
    y_fun <- signal::gausswin(y_dim, w = w)
  } else{
    stop("error func not recognised")
  }
    
  weight_mat   <- t(t(weight_mat * x_fun) * y_fun)
  weight_array <- array(weight_mat, c(1, x_dim, y_dim, 1, 1, 1, 1))
  
  if (z_dim != 1) weight_array <- rep_array_dim(weight_array, 4, z_dim)
  if (dyns != 1)  weight_array <- rep_array_dim(weight_array, 5, dyns)
  if (coils != 1) weight_array <- rep_array_dim(weight_array, 6, coils)
  if (N != 1)     weight_array <- rep_array_dim(weight_array, 7, N)
  
  # apply the weighting to the dataset
  mrs_data$data <- mrs_data$data * weight_array
  
  # transform back to image space
  mrs_data <- kspace2img_xy(mrs_data)
  return(mrs_data)
}

#' Grid shift MRSI data in the x/y dimension.
#' @param mrs_data MRSI data in the spatial domain.
#' @param x_shift shift to apply in the x-direction in units of voxels.
#' @param y_shift shift to apply in the y-direction in units of voxels.
#' @return shifted data.
#' @export
grid_shift_xy <- function(mrs_data, x_shift, y_shift) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  # TODO adjust pos vec to match
  mrsi_dims <- dim(mrs_data$data)
  x_dim <- mrsi_dims[2]
  y_dim <- mrsi_dims[3]
  N <- mrsi_dims[7]
  
  mrs_data <- img2kspace_xy(mrs_data)
  
  mat <- mrs_data$data
  mat <- drop(mat)
  dim(mat) <- c(x_dim, y_dim * N)
  
  mat <- mat * exp(1i * (seq(from = 0, to = (x_dim - 1) / x_dim,
                            length.out = x_dim) - 0.5) * x_shift * 2 * pi)
  
  dim(mat) <- c(x_dim, y_dim, N)
  mat <- aperm(mat, c(2, 1, 3))
  dim(mat) <- c(y_dim, x_dim * N)
  
  mat <- mat * exp(1i * (seq(from = 0, to = (y_dim - 1) / y_dim,
                            length.out = y_dim) - 0.5) * y_shift * 2 * pi)
  
  dim(mat) <- c(y_dim, x_dim, N)
  mat <- aperm(mat, c(2, 1, 3))
  dim(mat) <- mrsi_dims
  mrs_data$data <- mat
  
  # put xy dims back to space
  mrs_data <- kspace2img_xy(mrs_data)
  return(mrs_data)
}

#' Zero-fill MRSI data in the k-space x-y direction.
#' @param mrs_data MRSI data.
#' @param factor zero-filling factor, a factor of 2 returns a dataset with
#' twice the original points in the x-y directions. Factors smaller than one
#' are permitted, such that a factor of 0.5 returns half the k-space points in
#' the x-y directions.
#' @return zero-filled data.
#' @export
zf_xy <- function(mrs_data, factor = 2) {
  # TODO check this works for odd numbers of rows and cols...
  
  # check the input
  check_mrs_data(mrs_data)
  
  if (factor == 1) return(mrs_data)
  
  # put xy dims into k-space
  mrs_data <- img2kspace_xy(mrs_data)
  
  orig_dims <- dim(mrs_data$data)
  new_dims  <- orig_dims
  new_dims[2] <- new_dims[2] * factor
  new_dims[3] <- new_dims[3] * factor
  
  if (factor > 1) {
    new_data <- array(0, dim = new_dims)
    x_inds <- 1:orig_dims[2] + new_dims[2] / 2 - orig_dims[2] / 2
    y_inds <- 1:orig_dims[3] + new_dims[3] / 2 - orig_dims[3] / 2
    new_data[,x_inds, y_inds,,,,] <- mrs_data$data
    mrs_data$data <- new_data 
  } else {
    x_inds <- 1:new_dims[2] + orig_dims[2] / 2 - new_dims[2] / 2
    y_inds <- 1:new_dims[3] + orig_dims[3] / 2 - new_dims[3] / 2
    mrs_data$data <- mrs_data$data[, x_inds, y_inds,,,,,drop = FALSE]
  }
  
  # put xy dims back to spatial domain
  mrs_data <- kspace2img_xy(mrs_data)
}

hamming_vec <- function(vector) {
  vector * signal::hamming(length(vector))
}

# zero pad vector
zp_vec <- function(vector, n) {
  zp_vec <- rep(0, n)
  start_pt <- pracma::ceil((n - length(vector)) / 2) + 1
  zp_vec[start_pt:(start_pt + length(vector) - 1)] <- vector
  zp_vec
}

#' Combine coil data based on the first data point of a reference signal.
#' 
#' By default, elements are phased and scaled prior to summation. Where a 
#' reference signal is not given, the mean dynamic signal will be used
#' instead.
#' @param metab MRS data containing metabolite data.
#' @param ref MRS data containing reference data (optional).
#' @param noise MRS data from a noise scan (optional).
#' @param scale option to rescale coil elements based on the first data point
#' (logical).
#' @param scale_method one of "sig_noise_sq", "sig_noise" or "sig".
#' @param sum_coils sum the coil elements as a final step (logical).
#' @param noise_region the spectral region (in ppm) to estimate the noise.
#' @param average_ref_dyns take the mean of the reference scans in the dynamic
#' dimension before use.
#' @param ref_pt_index time-domain point to use for estimating phase and scaling 
#' values.
#' @return MRS data.
#' @export
comb_coils <- function(metab, ref = NULL, noise = NULL, scale = TRUE,
                       scale_method = "sig_noise_sq", sum_coils = TRUE,
                       noise_region = c(-0.5, -2.5), average_ref_dyns = TRUE,
                       ref_pt_index = 1) {
  
  if (class(metab) == "list") {
    # if (class(ref) != "list") stop("metab is a list but ref is not")
    # if (length(metab) != length(ref)) {
    #   stop("metab and ref must have the same length")
    # }
    
    if (!is.null(ref)) {
      if (class(ref) != "list") stop("ref is not a list but metab is")
      
      if (length(metab) != length(ref)) {
        stop("metab and w_ref must have the same length")
      }
    } else {
      ref <- vector("list", length(metab))
    }
    
    if (!is.null(noise)) {
      if (class(noise) != "list") stop("noise is not a list but metab is")
      
      if (length(metab) != length(noise)) {
        stop("metab and noise must have the same length")
      }
    } else {
      noise <- vector("list", length(metab))
    }
    
    more_args <- list(scale = scale, scale_method = scale_method,
                      sum_coils = sum_coils, noise_region = noise_region,
                      average_ref_dyns = average_ref_dyns,
                      ref_pt_index = ref_pt_index)
    
    res <- mapply(comb_coils, metab = metab, ref = ref, noise = noise,
                  MoreArgs = more_args, SIMPLIFY = FALSE)
    
    if (is.null(ref[[1]])) {
      return(res)
    } else {
      metab_list <- lapply(res, '[[', 1)
      ref_list   <- lapply(res, '[[', 2)
      out        <- list(metab = metab_list, ref = ref_list)
      class(out) <- c("list", "mrs_data")
      return(out)
    }
  } 
  
  metab_only <- FALSE
  if (is.null(ref)) {
    ref <- metab
    metab_only <- TRUE
  }
  
  if (!(scale_method %in% c("sig_noise_sq", "sig_noise", "sig"))) {
    stop("incorrect scale_method requested")
  }
  
  if (is_fd(metab)) metab <- fd2td(metab)
  
  if (is_fd(ref)) ref <- fd2td(ref)
  
  # get the dynamic mean of the ref data (better to take the first dynamic in
  # some cases?)
  if (average_ref_dyns) ref <- mean_dyns(ref)
  
  # ref_pt <- get_fp(ref)
  ref_pt <- get_subset(ref, td_set = ref_pt_index)$data
  phi <- Arg(ref_pt)
  amp <- Mod(ref_pt)
  
  # maintain original spatial scaling
  mean_amps <- apply(amp, c(1,2,3,4,5), mean)
  dim(mean_amps) <- c(dim(mean_amps), 1, 1)
  mean_amps <- rep_array_dim(mean_amps, 6, Ncoils(ref))
  amp <- amp / mean_amps
  
  if (scale & (scale_method != "sig")) {
    if (!is.null(noise)) {
      # estimate noise from noise data
      if (scale_method == "sig_noise_sq") {
        amp <- amp / (calc_coil_noise_sd(noise) ^ 2)
      } else {
        amp <- amp / calc_coil_noise_sd(noise)
      }
      
    } else {
      # estimate noise from first FID of the metab data
      metab_first <- get_dyns(metab, 1)
      noise_data <- crop_spec(metab_first, noise_region)
      noise_sd <- est_noise_sd(noise_data, offset = 0, n = Npts(noise_data),
                               p_order = 2)
      
      if (scale_method == "sig_noise_sq") {
        amp <- amp / (noise_sd ^ 2)
      } else {
        amp <- amp / noise_sd
      }
    }
  }
  
  # phase and scale ref data
  
  # repeat across the FID and dynamic dimensions
  ang <- rep_array_dim(phi, 7, Npts(ref))
  if (Ndyns(ref) > 1) ang <- rep_array_dim(ang, 5, Ndyns(ref))
  
  if (scale) {
    # repeat across the FID and dynamic dimensions
    scale_f <- rep_array_dim(amp, 7, Npts(ref))
    if (Ndyns(ref) > 1) scale_f <- rep_array_dim(scale_f, 5, Ndyns(ref))
    
    ref_ps <- ref
    ref_ps$data <- ref$data * exp(-1i * ang) * scale_f
  } else {
    ref_ps <- ref
    ref_ps$data <- ref$data * exp(-1i * ang)
  }
  
  if (sum_coils) ref_ps <- sum_coils(ref_ps)
  
  # phase and scale metab data
  
  # repeat across the FID and dynamic dimensions
  ang <- rep_array_dim(phi, 7, Npts(metab))
  if (Ndyns(metab) > 1) ang <- rep_array_dim(ang, 5, Ndyns(metab))
  
  if (scale) {
    # repeat across the FID and dynamic dimensions
    scale_f <- rep_array_dim(amp, 7, Npts(metab))
    if (Ndyns(metab) > 1) scale_f <- rep_array_dim(scale_f, 5, Ndyns(metab))
    
    metab_ps <- metab
    metab_ps$data <- metab$data * exp(-1i * ang) * scale_f
  } else {
    metab_ps <- metab
    metab_ps$data <- metab$data * exp(-1i * ang)
  }
  
  if (sum_coils) metab_ps <- sum_coils(metab_ps)
  
  if (metab_only) {
    return(metab_ps)
  } else {
    out <- list(metab = metab_ps, ref = ref_ps)
    return(out)
  }
}

#' Replicate a scan in the dynamic dimension.
#' @param mrs_data MRS data to be replicated.
#' @param times number of times to replicate.
#' @return replicated data object.
#' @export
rep_dyn <- function(mrs_data, times) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrs_data$data <- rep_array_dim(mrs_data$data, 5, times)
  mrs_data
}

#' Replicate a scan over a given dimension.
#' @param mrs_data MRS data to be replicated.
#' @param x_rep number of x replications.
#' @param y_rep number of y replications.
#' @param z_rep number of z replications.
#' @param dyn_rep number of dynamic replications.
#' @param coil_rep number of coil replications.
#' @param warn print a warning when the data dimensions do not change.
#' @return replicated data object.
#' @export
rep_mrs <- function(mrs_data, x_rep = 1, y_rep = 1, z_rep = 1, dyn_rep = 1,
                    coil_rep = 1, warn = TRUE) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  old_dims <- dim(mrs_data$data) 
  
  if (x_rep != 1) mrs_data$data <- rep_array_dim(mrs_data$data, 2, x_rep)
  if (y_rep != 1) mrs_data$data <- rep_array_dim(mrs_data$data, 3, y_rep)
  if (z_rep != 1) mrs_data$data <- rep_array_dim(mrs_data$data, 4, z_rep)
  if (dyn_rep != 1) mrs_data$data <- rep_array_dim(mrs_data$data, 5, dyn_rep)
  if (coil_rep != 1) mrs_data$data <- rep_array_dim(mrs_data$data, 6, coil_rep)
  
  if (warn & (identical(old_dims, dim(mrs_data$data)))) {
    warning("Data dimensions not changed.")
  }
  
  mrs_data
}

#' Estimate the standard deviation of the noise from a segment of an mrs_data
#' object.
#' @param mrs_data MRS data object.
#' @param n number of data points (taken from the end of array) to use in the
#' estimation.
#' @param offset number of final points to exclude from the calculation.
#' @param p_order polynomial order to fit to the data before estimating the
#' standard deviation.
#' @return standard deviation array.
#' @export
est_noise_sd <- function(mrs_data, n = 100, offset = 100, p_order = 2) {
  apply_mrs(mrs_data, 7, est_noise_sd_vec, n, offset, p_order, data_only = TRUE)
}

est_noise_sd_vec <- function(x, n = 100, offset = 100, p_order = 2) {
  N <- length(x)
  seg <- Re(x[(N - offset - n + 1):(N - offset)])
  lm_res <- stats::lm(seg ~ stats::poly(1:n, p_order))
  stats::sd(lm_res$residual)
}

#' Calculate the noise correlation between coil elements.
#' @param noise_data \code{mrs_data} object with one FID for each coil element.
#' @return correlation matrix.
#' @export
calc_coil_noise_cor <- function(noise_data) {
  
  # check the input
  check_mrs_data(noise_data) 
  
  cplx_data <- drop(noise_data$data)
  # concat real and imag parts
  real_data <- cbind(Re(cplx_data), Im(cplx_data))
  stats::cor(t(real_data))
}

#' Calculate the noise standard deviation for each coil element.
#' @param noise_data \code{mrs_data} object with one FID for each coil element.
#' @return array of standard deviations.
#' @export
calc_coil_noise_sd <- function(noise_data) {
  
  # check the input
  check_mrs_data(noise_data) 
  
  cplx_data <- drop(noise_data$data)
  # concat real and imag parts
  real_data <- cbind(Re(cplx_data), Im(cplx_data))
  apply(real_data, 1, stats::sd)
}

#' Calculate the spectral SNR.
#' 
#' SNR is defined as the maximum signal value divided by the standard deviation 
#' of the noise.
#' 
#' The mean noise value is subtracted from the maximum signal value to reduce DC
#' offset bias. A polynomial detrending fit (second order by default) is applied 
#' to the noise region before the noise standard deviation is estimated.
#' 
#' @param mrs_data an object of class \code{mrs_data}.
#' @param sig_region a ppm region to define where the maximum signal value
#' should be estimated.
#' @param noise_region a ppm region to defined where the noise level should be 
#' estimated.
#' @param p_order polynomial order to fit to the noise region before estimating 
#' the standard deviation.
#' @param interp_f interpolation factor to improve detection of the highest
#' signal value.
#' @param full_output output signal, noise and SNR values separately.
#' @return an array of SNR values.
#' @export
calc_spec_snr <- function(mrs_data, sig_region = c(4,0.5), 
                          noise_region = c(-0.5,-2.5), p_order = 2,
                          interp_f = 4, full_output = FALSE) {
  
  sig_data <- crop_spec(mrs_data, sig_region)
  noise_data <- crop_spec(mrs_data, noise_region)
  
  #max_sig <- apply_mrs(sig_data, 7, re_max, data_only = TRUE)
  max_sig <- apply_mrs(sig_data, 7, re_max_interp, interp_f, data_only = TRUE)
  noise_mean <- apply_mrs(noise_data, 7, re_mean, data_only = TRUE)
  max_sig <- max_sig - noise_mean
  
  #noise_sd <- apply_mrs(noise_data, 7, re_sd, data_only = TRUE)
  
  noise_sd <- est_noise_sd(noise_data, offset = 0, n = Npts(noise_data), 
                           p_order = p_order)
  
  snr <- max_sig / noise_sd
  
  # drop the last dimension for plotting functions
  snr <- abind::adrop(snr, 7)
  
  if (full_output) {
    max_sig  <- abind::adrop(max_sig, 7)
    noise_sd <- abind::adrop(noise_sd, 7)
    return(list(snr = snr, max_sig = max_sig, noise_sd = noise_sd))
  } else {
    return(snr)
  }
}

#' Search for the highest peak in a spectral region and return the frequency,
#' height and FWHM.
#' @param mrs_data an object of class \code{mrs_data}.
#' @param xlim frequency range (default units of PPM) to search for the highest 
#' peak.
#' @param interp_f interpolation factor, defaults to 4x.
#' @param scale the units to use for the frequency scale, can be one of: "ppm", 
#' "hz" or "points".
#' @param mode spectral mode, can be : "real", "imag" or "mod".
#' @return list of arrays containing the highest peak frequency, height and FWHM
#' in units of PPM and Hz.
#' @export
peak_info <- function(mrs_data, xlim = c(4,0.5), interp_f = 4, 
                           scale = "ppm", mode = "real") {
  
  mrs_data_crop <- crop_spec(mrs_data, xlim, scale)
  
  if (mode == "real") {
    mrs_data_crop$data <- Re(mrs_data_crop$data)
  } else if (mode == "imag") {
    mrs_data_crop$data <- Im(mrs_data_crop$data)
  } else if (mode == "mod") {
    mrs_data_crop$data <- Mod(mrs_data_crop$data)
  }
  
  res <- apply_mrs(mrs_data_crop, 7, calc_peak_info_vec, interp_f,
                   data_only = TRUE)
  
  pos_n <- res[,,,,,,1, drop = FALSE]
  pos_hz <- n2hz(pos_n, Npts(mrs_data_crop), fs(mrs_data_crop))
  pos_ppm <- hz2ppm(pos_hz, mrs_data_crop$ft, mrs_data_crop$ref)
  height <- res[,,,,,,2, drop = FALSE]
  fwhm_n <- res[,,,,,,3, drop = FALSE]
  fwhm_hz <- fwhm_n * fs(mrs_data_crop) / Npts(mrs_data_crop)
  fwhm_ppm <- fwhm_hz / mrs_data_crop$ft * 1e6 
  pos_ppm <- abind::adrop(pos_ppm, 7)
  pos_hz <- abind::adrop(pos_hz, 7)
  height <- abind::adrop(height, 7)
  fwhm_ppm <- abind::adrop(fwhm_ppm, 7)
  fwhm_hz <- abind::adrop(fwhm_hz, 7)
  list(freq_ppm = pos_ppm, freq_hz = pos_hz, height = height,
       fwhm_ppm = fwhm_ppm, fwhm_hz = fwhm_hz)
}

#' Calculate the FWHM of a peak from a vector of intensity values.
#' @param data_pts input vector.
#' @param interp_f interpolation factor to improve the FWHM estimate.
#' @return a vector of: x position of the highest data point, maximum peak
#' value in the y axis, FWHM in the units of data points.
#' @export
calc_peak_info_vec <- function(data_pts, interp_f) {
  data_pts <- stats::spline(data_pts, n = interp_f * length(data_pts))
  data_pts_x <- data_pts$x
  data_pts <- data_pts$y
  peak_pos_n <- which.max(data_pts)
  peak_height <- data_pts[peak_pos_n]
  hh <- peak_height / 2
  
  # right side of peak
  rs <- peak_pos_n + 
        min(which((data_pts < hh)[peak_pos_n:length(data_pts)])) - 1
  
  rs_slope <- (data_pts[rs] - data_pts[rs - 1])
  rs_intercept <- data_pts[rs] - rs_slope * rs
  rs_x_hh <- (hh - rs_intercept) / rs_slope
  
  # left side of peak
  ls <- peak_pos_n - min(which((data_pts < hh)[peak_pos_n:1])) + 1
  ls_slope <- (data_pts[ls + 1] - data_pts[ls])
  ls_intercept <- data_pts[ls] - ls_slope * ls
  ls_x_hh <- (hh - ls_intercept) / ls_slope
  
  fwhm <- (rs_x_hh - ls_x_hh) / interp_f
  
  array(c(data_pts_x[peak_pos_n], peak_height, fwhm))
}

#' Remove a constant baseline offset based on a reference spectral region.
#' @param mrs_data MRS data.
#' @param xlim spectral range containing a flat baseline region to measure the 
#' offset.
#' @return baseline corrected data.
#' @export
bc_constant <- function(mrs_data, xlim = NULL) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, bc_constant, xlim = xlim))
  }
  
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  offsets <- spec_op(mrs_data, xlim = xlim, mode = "cplx", operator = "mean")
  offsets_rep <- array(rep(offsets, Npts(mrs_data)), dim = dim(mrs_data$data))
  mrs_data$data <- mrs_data$data - offsets_rep
  return(mrs_data)
}



#' Baseline correction using the ALS method.
#' @param mrs_data mrs_data object.
#' @param lambda lambda parameter.
#' @param p p parameter.
#' @return baseline corrected data.
#' @export
bc_als <- function(mrs_data, lambda = 1e4, p = 0.001) {
  
  if (class(mrs_data) == "list") {
    return(lapply(mrs_data, bc_als, lambda = lambda, p = p))
  }
  
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  apply_mrs(mrs_data, 7, bc_als_vec, lambda, p)
}

bc_als_vec <- function(vec, lambda, p) {
  if (is.na(vec[1]))
    return(vec) 
  else {
    return(ptw::baseline.corr(Re(vec), lambda = lambda, p = p))
  }
}

#' Back extrapolate time-domain data points using an autoregressive model.
#' @param mrs_data mrs_data object.
#' @param extrap_pts number of points to extrapolate.
#' @param pred_pts number of points to base the extrapolation on.
#' @param method character string specifying the method to fit the model. Must
#' be one of the strings in the default argument (the first few characters are
#' sufficient). Defaults to "burg".
#' @param rem_add remove additional points from the end of the FID to maintain
#' the original length of the dataset. Default to TRUE.
#' @param ... additional arguments to specific methods, see ?ar.
#' @return back extrapolated data.
#' @export
back_extrap_ar <- function(mrs_data, extrap_pts, pred_pts = NULL,
                           method = "burg", rem_add = TRUE, ...) {
  
  # a time-domain operation
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  Np <- Npts(mrs_data) 
  
  mrs_data$data <- mrs_data$data[,,,,,,Np:1, drop = FALSE]
  mrs_data <- apply_mrs(mrs_data, 7, back_extrap_vec, extrap_pts, pred_pts,
                        method, ...)
  mrs_data$data <- mrs_data$data[,,,,,,(Np + extrap_pts):1, drop = FALSE]
  
  if (rem_add) mrs_data <- get_subset(mrs_data, td_set = 1:Np)
  
  dimnames(mrs_data$data) <- NULL
  
  return(mrs_data)
}

back_extrap_vec <- function(vec, extrap_pts, pred_pts, method, ...) {
  
  # use the full vector unless instructed otherwise
  if (is.null(pred_pts)) {
    pred_vec <- vec 
  } else {
    pred_vec <- utils::tail(vec, pred_pts)
  }
  
  new_pts_re <- as.numeric(stats::predict(stats::ar(Re(pred_vec),
                                                    method = method, ...),
                                          n.ahead = extrap_pts,
                                          se.fit = FALSE))
  
  new_pts_im <- as.numeric(stats::predict(stats::ar(Im(pred_vec),
                                                    method = method, ...),
                                          n.ahead = extrap_pts,
                                          se.fit = FALSE))
  
  out <- c(vec, new_pts_re + new_pts_im * 1i)
  return(out)
}

#' Calculate the sum of squares differences between two mrs_data objects.
#' @param mrs_data mrs_data object.
#' @param ref reference mrs_data object to calculate differences.
#' @param xlim spectral limits to perform calculation.
#' @return an array of the sum of squared difference values.
#' @export
calc_spec_diff <- function(mrs_data, ref = NULL, xlim = c(4, 0.5)) {
  
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  # diff from mean dynamic if ref not given
  if (is.null(ref)) ref <- mean_dyns(mrs_data)
  
  mrs_data_crop <- crop_spec(mrs_data, xlim)
  ref_crop <- crop_spec(ref, xlim)
  ref_crop <- rep_dyn(ref_crop, Ndyns(mrs_data))
  res <- mrs_data_crop - ref_crop
  apply_mrs(res, 7, cplx_sum_sq, data_only = TRUE)
}

#' Transform 2D MRSI data to k-space in the x-y direction.
#' @param mrs_data 2D MRSI data.
#' @return k-space data.
#' @export
img2kspace_xy <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrsi_dims <- dim(mrs_data$data) 
  x_dim <- mrsi_dims[2]
  y_dim <- mrsi_dims[3]
  coils <- mrsi_dims[6]
  N <- mrsi_dims[7]
  mat <- mrs_data$data
  mat <- drop(mat)
  dim(mat) <- c(x_dim, y_dim * N * coils)
  mat <- ft_shift_mat(mat)
  dim(mat) <- c(x_dim, y_dim, N * coils)
  mat <- aperm(mat, c(2, 1, 3))
  dim(mat) <- c(y_dim, x_dim * N * coils)
  mat <- ft_shift_mat(mat)
  dim(mat) <- c(y_dim, x_dim, N * coils)
  mat <- aperm(mat, c(2, 1, 3))
  dim(mat) <- mrsi_dims
  mrs_data$data <- mat
  return(mrs_data)
}

#' Transform 2D MRSI data from k-space to image space in the x-y direction.
#' @param mrs_data 2D MRSI data.
#' @return MRSI data in image space.
#' @export
kspace2img_xy <- function(mrs_data) {
  
  # check the input
  check_mrs_data(mrs_data) 
  
  mrsi_dims <- dim(mrs_data$data) 
  x_dim <- mrsi_dims[2]
  y_dim <- mrsi_dims[3]
  coils <- mrsi_dims[6]
  N <- mrsi_dims[7]
  mat <- mrs_data$data
  mat <- drop(mat)
  dim(mat) <- c(x_dim, y_dim * N * coils)
  mat <- ift_shift_mat(mat)
  dim(mat) <- c(x_dim, y_dim, N * coils)
  mat <- aperm(mat, c(2, 1, 3))
  dim(mat) <- c(y_dim, x_dim * N * coils)
  mat <- ift_shift_mat(mat)
  dim(mat) <- c(y_dim, x_dim, N * coils)
  mat <- aperm(mat, c(2, 1, 3))
  dim(mat) <- mrsi_dims
  mrs_data$data <- mat
  return(mrs_data)
}

#' Apply line-broadening to an mrs_data object to achieve a specified linewidth.
#' @param mrs_data data in.
#' @param lw target linewidth in units of ppm.
#' @param xlim region to search for peaks to obtain a linewidth estimate.
#' @return line-broadened data.
#' @export
set_lw <- function(mrs_data, lw, xlim = c(4, 0.5)) {
  
  if (class(mrs_data) == "list") {
    res <- lapply(mrs_data, set_lw, lw = lw, xlim = xlim)
    return(res)
  }
  
  # check the input
  check_mrs_data(mrs_data) 
  
  # start in the frequency-domain
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  # get an example spectrum for the data parameters 
  single_mrs <- get_subset(mrs_data, x_set = 1, y_set = 1, z_set = 1,
                           dyn_set = 1, coil_set = 1)
  
  lb_res <- apply_mrs(mrs_data, 7, optim_set_lw, lw, xlim, single_mrs,
                   data_only = TRUE)
  
  # apply the lb parameter to the full dataset
  res <- lb(mrs_data, lb_res)
  
  return(res)
}

optim_set_lw <- function(x, lw, xlim, single_mrs) {
  single_mrs$data[1,1,1,1,1,1,] <- x
  
  # return NA if the spectrum has been masked
  if (is.na(x[1])) return(NA)
  
  # measure current lw and check it is narrower than requested
  init_lw <- peak_info(single_mrs, xlim)$fwhm_ppm[1]
  if (init_lw > lw) {
    warning("Target linewidth is too narrow, masking spectrum.")
    return(NA)
  }
  
  # convert lw to Hz to get the upper value for 1D search
  upper_lw <- lw * single_mrs$ft / 1e6
  
  res <- stats::optim(0, lw_obj_fn, NULL, single_mrs, lw, xlim, lower = 0,
                      upper = upper_lw, method = "Brent")
  
  return(res$par[1])
}

lw_obj_fn <- function(lb_val, mrs_data, lw, xlim) {
  mrs_data <- lb(mrs_data, lb_val)
  new_lw   <- peak_info(mrs_data, xlim)$fwhm_ppm[1]
  Mod(new_lw - lw)
}

#' Perform l2 regularisation artefact suppression.
#' 
#' Perform l2 regularisation artefact suppression using the method proposed by
#' Bilgic et al. JMRI 40(1):181-91 2014.
#' @param mrs_data input data for artefact suppression.
#' @param thresh threshold parameter to extract lipid signals from mrs_data
#' based on the spectral integration of the thresh_xlim region in magnitude
#' mode.
#' @param b regularisation parameter.
#' @param A set of spectra containing the artefact basis signals. The thresh
#' parameter is ignored when A is specified.
#' @param xlim spectral limits in ppm to restrict the reconstruction range.
#' Defaults to the full spectral width.
#' @param thresh_xlim spectral limits in ppm to integrate for the threshold map.
#' @param A_append additional spectra to append to the A basis.
#' @param ret_norms return the residual norm and solution norms.
#' @return l2 reconstructed mrs_data object.
#' @export
l2_reg <- function(mrs_data, thresh = 0.05, b = 1e-11, A = NA, xlim = NA,
                   thresh_xlim = NULL, A_append = NULL, ret_norms = FALSE) {
  
  # generally done as a FD operation
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  # crop mrs_data to xlim if specified
  if (!anyNA(xlim)) mrs_data <- crop_spec(mrs_data, xlim)
  
  # get the data dimensions per coil
  res_dim <- dim(mrs_data$data)
  res_dim[6] <- 1
  
  if (!anyNA(A)) {
    if (!is_fd(A)) A <- td2fd(A)
    # crop A to xlim if specified
    if (!anyNA(xlim)) A <- crop_spec(A, xlim)
    A_coil <- t(stats::na.omit(mrs_data2mat(A)))
  }
  
  # prepare any extra signals and convert to matrix
  if (!is.null(A_append)) {
    if (!anyNA(xlim)) A_append <- crop_spec(A_append, xlim)
    A_append_mat <- t(stats::na.omit(mrs_data2mat(A_append)))
  }
  
  if (ret_norms) {
    resid_norm <- array(dim = Ncoils(mrs_data))
    soln_norm  <- array(dim = Ncoils(mrs_data))
  }
  
  for (coil in 1:Ncoils(mrs_data)) {
    
    mrs_data_coil <- get_subset(mrs_data, coil_set = coil)
    
    if (anyNA(A)) {
      map <- drop(spec_op(mrs_data_coil, mode = "mod", xlim = thresh_xlim))
      map_bool <- map > (max(map) * thresh)
      mrsi_mask <- mask_xy_mat(mrs_data_coil, mask = !map_bool)
      A_coil <- t(stats::na.omit(mrs_data2mat(mrsi_mask)))
    }
    
    # Append any extra signals
    if (!is.null(A_append)) A_coil <- cbind(A_coil, A_append_mat)
    
    # original data as a matrix
    x0 <- t(mrs_data2mat(mrs_data_coil))
    
    # recon. matrix 
    recon_mat <- solve(diag(nrow(A_coil)) + b * A_coil %*% Conj(t(A_coil)))
    
    # recon data
    x <- recon_mat %*% x0
    
    if (ret_norms) {
      # mrsi_mask <- mask_xy_mat(mrs_data_coil, mask = map_bool)
      # x0_mask <- t(stats::na.omit(mrs_data2mat(mrsi_mask)))
      # x_mask <-  recon_mat %*% x0_mask
      
      # residual norm
      resid_norm[coil] <- norm(x - x0, "2")
      # resid_norm[coil] <- sum(Mod(x - x0) ^ 2) ^ 0.5
      # resid_norm[coil] <- sum(rowSums(Mod(x - x0) ^ 2) ^ 0.5)
      
      # solution norm
      soln_norm[coil] <- norm(Conj(t(A_coil)) %*% x, "2")
      # soln_norm[coil] <- sum(Mod(Conj(t(A_coil)) %*% x) ^ 2) ^ 0.5
      # soln_norm[coil] <- sum(rowSums(Mod(Conj(t(A_coil)) %*% x) ^ 2) ^ 0.5)
    }
    
    x <- t(x)
    dim(x) <- res_dim
    mrs_data$data[,,,,,coil,] <- x
  }
    
  if (ret_norms) {
    return(list(mrs_data = mrs_data, resid_norm = resid_norm,
                soln_norm = soln_norm))
  } else {
    return(mrs_data)
  }
}

#' Signal space projection method for lipid suppression.
#' 
#' Signal space projection method as described in:
#' Tsai SY, Lin YR, Lin HY, Lin FH. Reduction of lipid contamination in MR
#' spectroscopy imaging using signal space projection. Magn Reson Med 2019
#' Mar;81(3):1486-1498.
#' 
#' @param mrs_data MRS data object.
#' @param comps the number of spatial components to use.
#' @param xlim spectral range (in ppm) covering the lipid signals.
#' @return lipid suppressed \code{mrs_data} object.
#' @export
ssp <- function(mrs_data, comps = 5, xlim = c(1.5, 0.8)) {
  
  # normally a FD operation
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  # chop out the lipid region
  D_ol_mrs <- crop_spec(mrs_data, xlim)
  
  # get the data dimensions per coil
  res_dim <- dim(mrs_data$data)
  res_dim[6] <- 1
  
  for (coil in 1:Ncoils(mrs_data)) {
    # extract the lipid region
    D_ol <- mrs_data2mat(get_subset(D_ol_mrs, coil_set = coil))
    
    # SVD
    svd_res <- svd(D_ol)
    
    # get the top spatial components
    U_m <- svd_res$u[,1:comps]
    
    # remove the lipids from the input
    D_ori <- mrs_data2mat(get_subset(mrs_data, coil_set = coil))
    
    if (is.matrix(U_m)) {
      P <- diag(nrow(U_m)) - U_m %*% Conj(t(U_m))
    } else{
      P <- diag(U_m) - U_m %*% Conj(t(U_m))
    }
    
    D_supp <- P %*% D_ori
    
    # restructure back into an mrs_data object
    dim(D_supp) <- res_dim
    
    # replace the input data points
    mrs_data$data[,,,,,coil,] <- D_supp
  }
  dimnames(mrs_data$data) <- NULL
  
  return(mrs_data)
}

#' Generate mrs_data from a table of single Lorentzian resonances.
#' @param reson_table as produced by the hsvd function.
#' @param acq_paras list of acquisition parameters. See
#' @param back_extrap_pts number of data points to back extrapolate
#' \code{\link{def_acq_paras}}
#' @return mrs_data object.
#' @export
reson_table2mrs_data <- function(reson_table, acq_paras = def_acq_paras(),
                                 back_extrap_pts = 0) {
  
  if (class(acq_paras) == "mrs_data") acq_paras <- get_acq_paras(acq_paras)
  
  sim_resonances(freq = reson_table$frequency_ppm, amp = reson_table$amplitude,
                 phase = reson_table$phase, lw = reson_table$lw_hz,
                 acq_paras = acq_paras, fp_scale = FALSE,
                 back_extrap_pts = back_extrap_pts)
}

#' Papoulis-Gerchberg (PG) algorithm method for k-space extrapolation.
#' 
#' PG method as described in: Haupt CI, Schuff N, Weiner MW, Maudsley AA.
#' Removal of lipid artifacts in 1H spectroscopic imaging by data extrapolation.
#' Magn Reson Med. 1996 May;35(5):678-87. Extrapolation is performed to expand
#' k-space coverage by a factor of 2, with the aim to reduce Gibbs ringing.
#' 
#' @param mrs_data MRS data object.
#' @param img_mask a boolean matrix of voxels with strong signals to be
#' extrapolated. Must be twice the dimensions of the input data.
#' @param kspace_mask a boolean matrix of kspace points that have been sampled.
#' Typically a circle for MRSI, but defaults to the full rectangular area of
#' k-space covered by the input data. Must match the x-y dimensions of the input
#' data.
#' @param intensity_thresh used to define img_mask based on the strength of the
#' signal in each voxel. Defaults to intensities greater than 15% of the
#' maximum. Ignored if img_mask is specified as argument.
#' @param iters number of iterations to perform.
#' @return extrapolated \code{mrs_data} object.
#' @export
pg_extrap_xy <- function(mrs_data, img_mask = NULL, kspace_mask = NULL,
                         intensity_thresh = 0.15, iters = 50) {
  
  # zero fill kspace by a factor of two
  mrs_data_zf <- zf_xy(mrs_data)
  
  # central mask of original k-space values
  # would typically be circular for MRSI but we assume rectangular unless told
  # otherwise
  kspace_mask_full <- matrix(FALSE, Nx(mrs_data_zf), Ny(mrs_data_zf)) 
  x_inds <- 1:Nx(mrs_data) + Nx(mrs_data_zf) / 2 - Nx(mrs_data) / 2
  y_inds <- 1:Ny(mrs_data) + Ny(mrs_data_zf) / 2 - Nx(mrs_data) / 2
  if (is.null(kspace_mask)) {
    # assume rectangular
    kspace_mask_full[x_inds, y_inds] <- TRUE
  } else {
    kspace_mask_full[x_inds, y_inds] <- kspace_mask
  }
  
  # make sure points outside the kspace mask are zero at the start
  mrs_data_ksp_zf <- img2kspace_xy(mrs_data_zf)
  mrs_data_ksp_zf_zerod <- mask_xy_mat(mrs_data_ksp_zf,
                                       mask = !kspace_mask_full, value = 0)
  mrs_data_zf <- kspace2img_xy(mrs_data_ksp_zf_zerod)
 
  # save a copy 
  mrs_data_zf_orig <- mrs_data_zf
  
  # get an image mask from the interpolated data, generally a scalp lipid mask
  if (is.null(img_mask)) {
    # img_map  <- drop(int_spec(mrs_data_zf, mode = "mod"))
    img_map  <- drop(spec_op(mrs_data_zf, mode = "mod"))
    img_mask <- img_map > (max(img_map) * intensity_thresh)
  }
  
  for (n in 1:iters) {
    # set all voxels not included in the image mask to zero
    mask_only <- mask_xy_mat(mrs_data_zf, mask = !img_mask, value = 0)
    
    # transform to kspace
    mask_only_ksp <- img2kspace_xy(mask_only)
    
    # set inner kspace values (in the kspace_mask_full) to zero to maintain the
    # original values in the following step
    mask_only_ksp <- mask_xy_mat(mask_only_ksp, kspace_mask_full, value = 0)
    
    # add new peripheral k-space values to the original zero-filled data
    mrs_data_zf <- kspace2img_xy(img2kspace_xy(mrs_data_zf_orig) +
                                     mask_only_ksp)
  }
  
  return(mrs_data_zf)
}

#' Return the mean of a list of mrs_data objects.
#' @param mrs_list list of mrs_data objects.
#' @return mean \code{mrs_data} object.
#' @export
mean_mrs_list <- function(mrs_list) {
  mean_mrs  <- sum_mrs_list(mrs_list) / length(mrs_list)
  return(mean_mrs)
}

#' Return the sum of a list of mrs_data objects.
#' @param mrs_list list of mrs_data objects.
#' @return sum \code{mrs_data} object.
#' @export
sum_mrs_list <- function(mrs_list) {
  sum_mrs <- mrs_list[[1]]
  for (n in (2:length(mrs_list))) sum_mrs <- sum_mrs(sum_mrs, mrs_list[[n]])
  return(sum_mrs)
}

#' Reconstruct 2D MRSI data from a twix file loaded with read_mrs.
#' @param twix_mrs raw dynamic data.
#' @return reconstructed data.
#' @export
recon_twix_2d_mrsi <- function(twix_mrs) {
  
  # figure out the required output array size
  max_lin <- max(twix_mrs$twix_inds$Lin) + 1
  max_seg <- max(twix_mrs$twix_inds$Seg) + 1
  
  twix_recon      <- twix_mrs 
  twix_recon$data <- array(0, c(1, max_lin, max_seg, 1, 1, Ncoils(twix_mrs),
                                Npts(twix_mrs)))
  
  # map the twix indices to a full data array
  for (n in 1:nrow(twix_mrs$twix_inds)) {
    x_ind <- twix_mrs$twix_inds$Seg[n] + 1
    y_ind <- twix_mrs$twix_inds$Lin[n] + 1
    twix_recon$data[1, x_ind, y_ind, 1, 1, , ] <- twix_mrs$data[1, 1, 1, 1, n, , ]
  }
  
  # invert every other k-space line (for some reason)
  k_sp_corr_x <- array(1, dim(twix_recon$data))
  k_sp_corr_y <- k_sp_corr_x
  k_sp_corr_x[,c(T, F),,,,,] <- -1
  k_sp_corr_y[,,c(F, T),,,,] <- -1
  twix_recon$data <- twix_recon$data * k_sp_corr_x * k_sp_corr_y
  
  # kspace to image space
  twix_recon <- kspace2img_xy(twix_recon)
  
  # mirror in the x-direction  
  twix_recon <- get_subset(twix_recon, x_set = Nx(twix_recon):1)
  
  return(twix_recon)
}

#' Fade a spectrum to zero by frequency domain multiplication with a tanh
#' function. Note this operation distorts data points at the end of the FID.
#' @param mrs_data data to be faded.
#' @param start_ppm start point of the fade in ppm units.
#' @param end_ppm end point of the fade in ppm units.
#' @return modified mrs_data object.
#' @export
zero_fade_spec <- function(mrs_data, start_ppm, end_ppm) {
  
  # fd operation
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  fade_spec <- get_voxel(mrs_data)
  fade_spec$data[,,,,,,] <- 0
  
  ppm_scale <- ppm(mrs_data)
  inds      <- get_seg_ind(ppm_scale, start_ppm, end_ppm)
  
  if (start_ppm < end_ppm) {
    fade <- (tanh(seq(from = -3, to = 3, length.out = length(inds))) + 1) / 2
    fade_spec$data[,,,,,,(inds[length(inds)] + 1):Npts(mrs_data)] <- 1
  } else {
    fade <- (tanh(seq(from = 3, to = -3, length.out = length(inds))) + 1) / 2
    fade_spec$data[,,,,,,0:(inds[1] - 1)] <- 1
  }

  fade_spec$data[,,,,,,inds] <- fade
  
  fade_spec <- rep_mrs(fade_spec, x_rep = Nx(mrs_data), y_rep = Ny(mrs_data),
                       z_rep = Nz(mrs_data), dyn_rep = Ndyns(mrs_data),
                       coil_rep = Ncoils(mrs_data), warn = FALSE)
  
  return(mrs_data * fade_spec)
}