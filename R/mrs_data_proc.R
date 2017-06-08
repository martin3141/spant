#' Simulate a MRS data object containing a set of simulated resonances.
#' @param freq Resonance frequency.
#' @param amp Resonance amplitude.
#' @param lw Line width in Hz.
#' @param lg Lorentz-Gauss lineshape parameter (between 0 and 1).
#' @param phase Phase in degrees.
#' @param freq_ppm Frequencies are given in ppm units if set to TRUE, otherwise
#' Hz are assumed.
#' @param acq_paras List of acquisition parameters. See
#' \code{\link{def_acq_paras}}
#' @return MRS data object.
#' @examples
#' sim_data <- sim_resonances(freq = 2, lw = 5)
#' @export
sim_resonances <- function(freq = 0, amp = 1, lw = 0, lg = 0, phase = 0, 
                           freq_ppm = TRUE, acq_paras = def_acq_paras()) {
  
  # TODO check this works for vectors
  #if ((sum(lg > 1) + sum(lg < 0)) > 0) {
  #  cat("Error, lg values not between 0 and 1.")  
  #  stop()
  #}
  
  sig_n <- length(freq)
  if (sig_n != length(amp)) {
    amp <- rep_len(amp, sig_n)
  }
  
  if (sig_n != length(lw)) {
    lw <- rep_len(lw, sig_n)
  }
  
  if (sig_n != length(phase)) {
    phase <- rep_len(phase, sig_n)
  }
  
  if (sig_n != length(lg)) {
    lg <- rep_len(lg, sig_n)
  }
  
  # generate data in TD
  t <- seq(from = 0, to = (acq_paras$N - 1) / acq_paras$fs,
           by = 1 / acq_paras$fs)
  
  # covert freqs to Hz
  if (freq_ppm) {
    f_hz <- (acq_paras$ref - freq) * acq_paras$ft / 1e6
  } else {
    f_hz <- freq
  }
  
  data <- rep(0, acq_paras$N)
  for (n in 1:sig_n) {
    temp_data <- amp[n] * exp(1i * pi * phase[n] / 180 + 2i * pi * f_hz[n] *
                          t - (1 - lg[n]) * lw[n] * t * pi)
    
    if (lg[n] > 0) {
      temp_data <- temp_data * exp(((lg[n] * lw[n]) ^ 2 * pi ^ 2 / 4 /
                               log(0.5)) * (t ^ 2))
    }
    
    data <- data + temp_data
  }
  
  #if (lg < 1) {
  #  mrs_data$data = mrs_data$data*exp(-(1-lg)*lb*t*pi)
  #}
  
  #if (lg > 0) {
  #  mrs_data$data = mrs_data$data*exp(((lg*lb)^2*pi^2/4/log(0.5))*(t^2))
  #}
  
  data <- array(data,dim = c(1, 1, 1, 1, 1, 1, acq_paras$N))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / acq_paras$fs)
  mrs_data <- list(ft = acq_paras$ft, data = data, resolution = res, te = 0, 
                   ref = acq_paras$ref, row_vec = c(1,0,0), col_vec = c(0,1,0),
                   pos_vec = c(0,0,0), freq_domain = rep(FALSE, 7))
  
  class(mrs_data) <- "mrs_data"
  return(mrs_data)
}

sim_resonances_fast <- function(freq = 0, amp = 1, freq_ppm = TRUE,
                                N = def_N(), fs = def_fs(), ft = def_ft(),
                                ref = def_ref()) {
  
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
  
  data <- array(data,dim = c(1, 1, 1, 1, 1, 1, N))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  mrs_data <- list(ft = ft, data = data, resolution = res, te = 0, ref = ref, 
                  row_vec = c(1,0,0), col_vec = c(0,1,0), pos_vec = c(0,0,0), 
                  freq_domain = rep(FALSE, 7))
  
  class(mrs_data) <- "mrs_data"
  
  return(mrs_data)
}

sim_resonances_fast2 <- function(freq = 0, amp = 1, freq_ppm = TRUE,
                                 N = def_N(), fs = def_fs(), ft = def_ft(), 
                                 ref = def_ref()) {
  
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
  
  data <- array(data, dim = c(1, 1, 1, 1, 1, 1, N))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  mrs_data <- list(ft = ft, data = data, resolution = res, te = 0, ref = ref, 
                   row_vec = c(1,0,0), col_vec = c(0,1,0), pos_vec = c(0,0,0), 
                   freq_domain = rep(FALSE, 7))
  
  class(mrs_data) <- "mrs_data"
  return(mrs_data)
}

vec2mrs_data <- function(vec, fs = def_fs(), ft = def_ft(), ref = def_ref(),
                         dyns = 1, fd = FALSE) {
  
  data <- array(vec, dim = c(length(vec), dyns))
  data <- aperm(data,c(2, 1))
  dim(data) <- c(1, 1, 1, 1, dyns, 1, length(vec))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  mrs_data <- list(ft = ft, data = data, resolution = res, te = 0, ref = ref, 
                  row_vec = c(1,0,0), col_vec = c(0,1,0), pos_vec = c(0,0,0), 
                  freq_domain = c(rep(FALSE, 6), fd))
  
  class(mrs_data) <- "mrs_data"
  return(mrs_data)
}

array2mrs_data <- function(data_array, fs = def_fs(), ft = def_ft(),
                           ref = def_ref(), fd = FALSE) {
  
  if (length(dim(data_array)) != 7) stop("Incorrect number of dimensions.")
  
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  mrs_data <- list(ft = ft, data = data_array, resolution = res, te = 0,
                   ref = ref, row_vec = c(1,0,0), col_vec = c(0,1,0),
                   pos_vec = c(0,0,0), freq_domain = c(rep(FALSE, 6), fd))
  
  class(mrs_data) <- "mrs_data"
  return(mrs_data)
}

#' Convert mrs_data object to a matrix, with spectral points in the row dimension
#' and dynamics in the column dimension.
#' @param mrs_data MRS data object.
#' @return A matrix.
#' @export
mrs_data2mat <- function(mrs_data) {
  t(as.matrix(mrs_data$data[1,1,1,1,,1,]))
}

#' Convert a matrix (with spectral points in the row dimension and dynamics in
#' the column dimensions) into a mrs_data object.
#' @param mat The data matrix.
#' @param fs Sampling frequency in Hz.
#' @param ft Transmitter frequency in Hz.
#' @param ref Reference value for ppm scale.
#' @param fd Flag to indicate if the matrix is in the frequency domain (logical).
#' @export
mat2mrs_data <- function(mat, fs = def_fs(), ft = def_ft(), ref = def_ref(),
                         fd = FALSE) {
  
  data <- array(t(mat), dim = c(1, 1, 1, 1, ncol(mat), 1, nrow(mat)))
  res <- c(NA, 1, 1, 1, 1, NA, 1 / fs)
  mrs_data <- list(ft = ft, data = data, resolution = res, te = 0, ref = ref, 
                   row_vec = c(1,0,0), col_vec = c(0,1,0), pos_vec = c(0,0,0), 
                   freq_domain = c(rep(FALSE, 6), fd))
  
  class(mrs_data) <- "mrs_data"
  return(mrs_data)
}

#' Simulate a time-domain mrs_data object containing simulated Gaussian noise.
#' @param sd standard deviation of the noise.
#' @param fs Sampling frequency in Hz.
#' @param ft Transmitter frequency in Hz.
#' @param N Number of data points in the spectral dimension.
#' @param ref Reference value for ppm scale.
#' @param dyns Number of dynamic scans to generate.
#' @return mrs_data object.
#' @export
sim_noise <- function(sd = 0.1, fs = def_fs(), ft = def_ft(), N = def_N(),
                      ref = def_ref(), dyns = 1) {
 
  data_pts <- dyns * N 
  # generate data in TD
  vec <- stats::rnorm(data_pts, 0, sd) + 1i*stats::rnorm(data_pts, 0, sd)
  data_array <- array(vec, dim = c(1, 1, 1, 1, dyns, 1, N))
  array2mrs_data(data_array, fs = fs, ft = ft, ref = ref)
}

sim_zeros <- function(fs = def_fs(), ft = def_ft(), N = def_N(),
                      ref = def_ref(), dyns = 1) {
  
  data_pts <- dyns * N 
  vec <- rep(0, data_pts) * 1i
  data_array <- array(vec, dim = c(1, 1, 1, 1, dyns, 1, N))
  array2mrs_data(data_array, fs = fs, ft = ft, ref = ref)
}

apply_mrs <- function(mrs_data, dims, fun, ..., data_only = FALSE) {
  dims <- sort(dims)
  margins <- c(1:7)[-dims]
  mrs_data$data <- plyr::aaply(mrs_data$data, margins, fun, ..., .drop = FALSE)
  perm_vec <- 1:(7 - length(dims))
  for (n in (1:length(dims))) {
    perm_vec <- append(perm_vec, (n + 7 - length(dims)), after = dims[n] - 1)
  }
  
  #print(perm_vec)
  mrs_data$data <- aperm(mrs_data$data, perm_vec)
  if ( data_only == FALSE ) {
    return(mrs_data)
  } else {
    return(mrs_data$data)
  }
}

#' Apply phasing parameters to MRS data.
#' @param mrs_data MRS data.
#' @param zero_order Zero'th order phase term in degrees.
#' @param first_order First order (frequency dependent) phase term in ms.
#' @return MRS data with applied phase parameters.
#' @export
phase <- function(mrs_data, zero_order, first_order = 0) {
  if (first_order == 0) {
    mrs_data$data = mrs_data$data * exp(1i * zero_order * pi / 180)
  } else {
    freq <- rep(hz(mrs_data), each = Nspec(mrs_data))
    freq_mat <- array(freq, dim = dim(mrs_data$data))
    # needs to be a freq-domain operation
    if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
    }
    mrs_data$data = mrs_data$data * exp(2i * pi * (zero_order / 360 - freq_mat
                                                   * first_order / 1000))
  }
  return(mrs_data)
}

#' Perform a zeroth order phase correction based on the phase of the first data 
#' point in the time-domain.
#' @param mrs_data MRS data to be corrected.
#' @param ret_phase Return phase values (logical).
#' @return Corrected data or a list with corrected data and optional phase 
#' values.
#' @export
fp_phase_correct <- function(mrs_data, ret_phase = FALSE) {
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) {
    mrs_data <- fd2td(mrs_data)
  }
  
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
#' @return First time-domain data point.
#' @export
get_fp <- function(mrs_data) {
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) {
    mrs_data <- fd2td(mrs_data)
  }
  
  # drop the chem shift dimension
  mrs_data$data[,,,,,, 1, drop = F]
}

fp_mag <- function(mrs_data) {
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) {
    mrs_data <- fd2td(mrs_data)
  }
  
  # drop the chem shift dimension
  abind::adrop(Mod(mrs_data$data[,,,,,, 1, drop = F]), 7)
}

#' Return the phase of the first data point in the time-domain.
#' @param mrs_data MRS data.
#' @return Phase values in degrees.
#' @export
fp_phase <- function(mrs_data) {
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) {
    mrs_data <- fd2td(mrs_data)
  }
  
  # drop the chem shift dimension
  abind::adrop(Arg(mrs_data$data[,,,,,, 1, drop = F]), 7) * 180 / pi
}

#' Conjugate MRS data.
#' @param mrs_data input data
#' @return Conjugated data
#' @export
conj <- function(mrs_data) {
    mrs_data$data = Re(mrs_data$data) - Im(mrs_data$data) * 1i
    mrs_data
}

#' Apply line-broadening (apodisation) to MRS data.
#' @param mrs_data input data
#' @param lb amount on line-broadening in Hz
#' @param lg Lorentz-Gauss lineshape parameter (between 0 and 1)
#' @return line-broadened data
#' @export
lb <- function(mrs_data, lb, lg = 1) {
  if ((sum(lg > 1) + sum(lg < 0)) > 0 ) {
    cat("Error, lg values not between 0 and 1.")  
    stop()
  }
  
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) {
    mrs_data <- fd2td(mrs_data)
  }
  t <- rep(seconds(mrs_data), each = Nspec(mrs_data))
  
  if (lg < 1) {
    mrs_data$data = mrs_data$data * exp(-(1 - lg) * lb * t * pi)
  }
  
  if (lg > 0) {
    mrs_data$data = mrs_data$data * exp(((lg * lb) ^ 2 * pi ^ 2 / 4 /
                                          log(0.5)) * (t ^ 2))
  }
  return(mrs_data)
}

#' Zero-fill MRS data in the time domain.
#' @param mrs_data MRS data.
#' @param factor Zero-filling factor, factor of 2 returns a dataset with
#' twice the original data points.
#' @return Zero-filled data.
#' @export
zf <- function(mrs_data, factor = 2) {
  set_td_pts(mrs_data, factor * N(mrs_data))
}

#' Set the number of time-domain data points, truncating or zero-filling as
#' appropriate.
#' @param mrs_data MRS data.
#' @param pts Number of data points.
#' @return MRS data with pts data points.
#' @export
set_td_pts <- function(mrs_data, pts) {
  # needs to be a time-domain operation
  if (is_fd(mrs_data)) {
    mrs_data <- fd2td(mrs_data)
  }
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

# TODO should be set_ref really
ref <- function(mrs_data, ref) {
  mrs_data$ref = ref
  return(mrs_data)
}

is_fd <- function(mrs_data) {
  mrs_data$freq_domain[7]
}

#' Transform time-domain data to the frequency-domain.
#' @param mrs_data MRS data in time-domain representation.
#' @return MRS data in frequency-domain representation.
#' @export
td2fd <- function(mrs_data) {
  mrs_data <- ft(mrs_data, 7)
  mrs_data$freq_domain[7] = TRUE
  return(mrs_data)
}

#' Transform frequency-domain data to the time-domain.
#' @param mrs_data MRS data in frequency-domain representation.
#' @return MRS data in time-domain representation.
#' @export
fd2td <- function(mrs_data) {
  mrs_data <- ift(mrs_data, 7)
  mrs_data$freq_domain[7] = FALSE
  return(mrs_data)
}

# recon complex td data from real part of fd data
recon_imag <- function(mrs_data) {
  # data needs to be in the FD
  if (!is_fd(mrs_data)) {
    mrs_data <- td2fd(mrs_data)
  }
  mrs_data <- apply_mrs(mrs_data, 7, recon_imag_vec)
  mrs_data$freq_domain[7] = FALSE
  return(mrs_data)
}

#' Return acquisition parameters from a MRS data object.
#' @param mrs_data MRS data.
#' @return A list of acquisition parameters.
#' @export
get_acq_paras <- function(mrs_data) {
  list(ft = mrs_data$ft, fs = fs(mrs_data), N = N(mrs_data), ref = mrs_data$ref)
}

ft <- function(mrs_data, dims) {
  apply_mrs(mrs_data, dims, ft_shift)
}

ift <- function(mrs_data, dims) {
  apply_mrs(mrs_data, dims, ift_shift)
}

dim.mrs_data <- function(x) {
  dim(x$data)
}

#' Return the number of data points in an MRS dataset.
#' @param mrs_data MRS data.
#' @return Number of data points.
#' @export
N <- function(mrs_data) {
  dim(mrs_data$data)[7]
}

#' Return the number of dynamic scans in an MRS dataset.
#' @param mrs_data MRS data.
#' @return Number of dynamic scans
#' @export
dyns <- function(mrs_data) {
  dim(mrs_data$data)[5]
}

Nspec <- function(mrs_data) {
  mrs_dims <- dim(mrs_data$data)
  (mrs_dims[1] * mrs_dims[2] * mrs_dims[3] * mrs_dims[4] * mrs_dims[5] *
   mrs_dims[6])
}

#' Return the sampling frequency in Hz of an MRS dataset.
#' @param mrs_data MRS data.
#' @return Sampling frequency in Hz.
#' @export
fs <- function(mrs_data) {
  1 / mrs_data$resolution[7]
}

hz <- function(mrs_data, fs = NULL, N = NULL) {
  if (is.null(fs)) {
    fs <- fs(mrs_data)
  }
  
  if (is.null(N)) {
    N <- N(mrs_data)
  }
     
  seq(from = -fs / 2, to = fs / 2 - fs / N, length.out = N)
}

ppm <- function(mrs_data, ft = NULL, ref = NULL, fs= NULL, N = NULL) {
  if (is.null(ft)) {
    ft <- mrs_data$ft
  }
  
  if (is.null(ref)) {
    ref <- mrs_data$ref
  }
  
  if (is.null(fs)) {
    fs <- fs(mrs_data)
  }
  
  if (is.null(N)) {
    N <- N(mrs_data)
  }
  
  -hz(fs = fs, N = N) / mrs_data$ft * 1e6 + mrs_data$ref
}

hz2ppm <- function(hz_in, ft, ref) {
  ref - hz_in / ft * 1e6
}

ppm2hz <- function(ppm_in, ft, ref) {
  (ref - ppm_in) * ft / 1e6
}

pts <- function(mrs_data) {
  seq(from = 1, to = N(mrs_data))  
}

#' Return a time scale vector to match the FID of an MRS data object.
#' @param mrs_data MRS data.
#' @return A time scale vector in units of seconds.
seconds <- function(mrs_data) {
  fs <- fs(mrs_data)
  seq(from = 0, to = (N(mrs_data) - 1)/fs, by = 1 / fs)
}

get_seg_ind <- function(scale, start, end) {
  st_ind  <- sum(scale <= start)  
  end_ind <- sum(scale <= end)  
  if ( scale[1] > scale[2] ) {
    st_ind  <- length(scale) - st_ind + 1
    end_ind <- length(scale) - end_ind + 1
  }
  c(st_ind, end_ind)
}

crop_spec <- function(mrs_data, xlim = c(4,0.5), scale = "ppm") {
  # needs to be a fd operation
  if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
  }
  
  if (scale == "ppm") {
    x_scale <- ppm(mrs_data)
  } else if (scale == "hz") {
    x_scale <- hz(mrs_data)
  } else if (scale == "points") {
    x_scale <- pts(mrs_data)
  } else if (scale == "seconds") {
    x_scale <- seconds(mrs_data)
  }
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(mrs_data)])
  }
  
  x_inds <- get_seg_ind(x_scale, xlim[1], xlim[2])
  subset <- x_inds[1]:x_inds[2]
  
  old_ppm <- ppm(mrs_data)
  
  # update fs
  mrs_data$resolution[7] <- (mrs_data$resolution[7] / length(subset) *
                             N(mrs_data))
  
  mrs_data$data <- mrs_data$data[,,,,,, subset, drop = F]
  
  # update ref TODO +1 may only be needed in some cases?
  new_ppm = old_ppm[which.min(abs(hz(mrs_data))) + x_inds[1] + 1]
  mrs_data$ref <- new_ppm
    
  mrs_data
}

#' Align spectra to a reference frequency using a convolution based method.
#' @param mrs_data Data to be aligned.
#' @param ref_peak A reference frequency in ppm units.
#' @param zf_factor Zero filling factor to increase alignment resolution.
#' @param lb Line broadening to apply to the reference signal.
#' @param max_shift Maximum allowable shift in Hz.
#' @param ret_df Return frequency shifts in addition to aligned data (logical).
#' @return Aligned data object.
#' @export
align <- function(mrs_data, ref_peak = 4.65, zf_factor = 2, lb = 2,
                  max_shift = 20, ret_df = FALSE) {
  
  if (is_fd(mrs_data)) {
    mrs_data <- fd2td(mrs_data)
  }
  
  # zero fill
  #mrs_data_zf <- set_td_pts(mrs_data,zf*N(mrs_data))
  mrs_data_zf <- zf(mrs_data, zf_factor)
  mrs_data_zf <- td2fd(mrs_data_zf)
  freq <- -(ref_peak - mrs_data$ref) * mrs_data$ft * 1.0e-6 # TODO use ppm2hz
  #print(freq)
  t_zf <- seconds(mrs_data_zf)
  ref_data <- ft_shift(exp(2i * t_zf * pi * freq - lb * t_zf * pi))
  #plot(Re(ref_data_td),type="l")
  window <- floor(max_shift * N(mrs_data_zf) * mrs_data$resolution[7])
  
  #plot(Re(ref_data[1000:1500])/1000,type="l")
  #lines(Re(mrs_data_zf$data[1,1,1,1,1,1,1000:1500]),type="l")
  #shift_hz <- conv_align(mrs_data_zf$data[1,1,1,1,1,1,], ref_data, window,
  # 1/mrs_data$resolution[7])
  #print(shift_hz)
  #t_orig <- seconds(mrs_data)
  #mrs_data$data <- mrs_data$data * exp(2i*pi*t_orig*shift_hz)
  
  shifts <- apply_mrs(mrs_data_zf, 7, conv_align, ref_data, window,
                      1/mrs_data$resolution[7], data_only = TRUE)
  
  #return(shifts)
  t_orig <- rep(seconds(mrs_data), each = Nspec(mrs_data))
  t_array <- array(t_orig, dim = dim(mrs_data$data))
  shift_array <- array(shifts, dim = dim(mrs_data$data))
  shift_array <- exp(2i * pi * t_array * shift_array)
  mrs_data$data <- mrs_data$data * shift_array
  if (ret_df) {
    return(list(mrs_data, abind::adrop(shifts, 7)))
  } else {
    return(mrs_data)
  }
}

get_fwhm <- function(mrs_data) {
  
  if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
  }
  
  fwhm <- apply_mrs(mrs_data, 7, calc_fwhm)
  
  # convert points to ppm
  fwhm_ppm <- fwhm$data * fs(mrs_data) / N(mrs_data) / mrs_data$ft * 1e6
  
  abind::adrop(fwhm_ppm, 7)
}

#' Return an array of amplitudes derived from fitting the initial points in the
#' time domain and extrapolating back to t=0.
#' @param mrs_data MRS data.
#' @param nstart The first data point to fit.
#' @param nend The last data point to fit.
#' @return An array of amplitudes.
#' @export
get_td_amp <- function(mrs_data, nstart = 10, nend = 50) {
  
  if (is_fd(mrs_data)) {
      mrs_data <- fd2td(mrs_data)
  }
  t <- seconds(mrs_data)
  amps <- apply_mrs(mrs_data, 7, measure_lorentz_amp, t, nstart, nend)$data
  
  abind::adrop(amps, 7)
  amps
}

calc_fwhm <- function(mrs_data) {
  mrs_data <- Re(mrs_data)
  max_pt <- which.max(mrs_data)
  max_val <- mrs_data[max_pt]
  gthm <- mrs_data > (max_val / 2) # greater than half max.
  sum(gthm)
}

conv_align <- function(acq, ref, window, fs) {
  conv <- pracma::fftshift(Mod(stats::convolve(acq, ref)))
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
#' @param mrs_data Dynamic MRS data.
#' @param subset A vector containing indices to the dynamic scans to be 
#' returned.
#' @return MRS data containing the subset of requested dynamics.
#' @export
get_dyns <- function(mrs_data, subset) {
  mrs_data$data <- mrs_data$data[,,,, subset,,, drop = FALSE]
  return(mrs_data)
}

set_dyns <- function(mrs_data, subset, mrs_data_in) {
  mrs_data$data[,,,, subset,,] = mrs_data_in$data[,,,, 1,,]
  return(mrs_data)
}

#' Remove a subset of dynamic scans.
#' @param mrs_data Dynamic MRS data.
#' @param subset A vector containing indices to the dynamic scans to be 
#' removed.
#' @return MRS data without the specified dynamic scans.
#' @export
rm_dyns <- function(mrs_data, subset) {
  mrs_data$data <- mrs_data$data[,,,, -subset,,, drop = F]
  mrs_data
}

get_voxel <- function(mrs_data, x = 1, y = 1, z = 1, ref = 1, coil = 1) {
  mrs_data$data <- mrs_data$data[ref, x, y, z, coil,,, drop = FALSE]
  return(mrs_data)
}

get_slice <- function(mrs_data, slice = 1, ref = 1, coil = 1) {
  mrs_data$data <- mrs_data$data[ref,,, slice, coil,,, drop = FALSE]
  return(mrs_data)
}
  
#' Return odd numbered dynamic scans starting from 1 (1,3,5...).
#' @param mrs_data Dynamic MRS data.
#' @return Dynamic MRS data containing odd numbered scans.
#' @export
get_odd_dyns <- function(mrs_data) {
  subset <- seq(1, dyns(mrs_data), 2)
  get_dyns(mrs_data, subset)
}

#' Return even numbered dynamic scans starting from 1 (2,4,6...).
#' @param mrs_data Dynamic MRS data.
#' @return Dynamic MRS data containing even numbered scans.
#' @export
get_even_dyns <- function(mrs_data) {
  subset <- seq(2, dyns(mrs_data), 2)
  get_dyns(mrs_data, subset)
}

#' Invert odd numbered dynamic scans starting from 1 (1,3,5...).
#' @param mrs_data Dynamic MRS data.
#' @return Dynamic MRS data with inverted odd numbered scans.
#' @export
inv_odd_dyns <- function(mrs_data) {
  subset <- seq(1, dyns(mrs_data), 2)
  mrs_data$data[,,,, subset,,] <- -1 * mrs_data$data[,,,, subset,,]
  return(mrs_data)
}

#' Invert even numbered dynamic scans starting from 1 (2,4,6...).
#' @param mrs_data Dynamic MRS data.
#' @return Dynamic MRS data with inverted even numbered scans.
#' @export
inv_even_dyns <- function(mrs_data) {
  subset <- seq(2, dyns(mrs_data), 2)
  mrs_data$data[,,,, subset,,] <- -1 * mrs_data$data[,,,, subset,,]
  return(mrs_data)
}


#' Combine a reference and metabolite mrs_data object.
#' @param metab metabolite mrs_data object.
#' @param ref reference mrs_data object.
#' @return Combined metabolite and reference mrs_data object.
#' @export
combine_metab_ref <- function(metab, ref) {
  metab$data <- abind::abind(metab$data, ref$data, along = 1)
  metab
}

#' Extract the reference component from an mrs_data object.
#' @param mrs_data MRS data.
#' @return Reference component.
#' @export
get_ref <- function(mrs_data) {
  mrs_data$data <- mrs_data$data[2,,,,,,,drop = FALSE]
  mrs_data
}

#' Extract the metabolite component from an mrs_data object.
#' @param mrs_data MRS data.
#' @return Metabolite component.
#' @export
get_metab <- function(mrs_data) {
  mrs_data$data <- mrs_data$data[1,,,,,,,drop = FALSE]
  mrs_data
}


#' Append MRS data across the dynamic dimension.
#' @param ... MRS data objects.
#' @return A single MRS data object with the input objects concatenated together.
#' @export
append_dyns <- function(...) {
  x <- list(...)
  
  first_dataset <- x[[1]]
  
  if (is_fd(first_dataset)) {
        first_dataset <- fd2td(first_dataset)
  }
  
  # data needs to be in the same domain
  for (n in 1:length(x)) {
    if (is_fd(x[[n]])) {
        x[[n]] <- fd2td(x[[n]])
    }
    x[[n]] <- x[[n]]$data
  }
  
  new_data <- abind::abind(x, along = 5)
  first_dataset$data <- new_data
  first_dataset
}

split_metab_ref <- function(mrs_data) {
  metab <- mrs_data
  ref <- mrs_data
  metab$data = metab$data[1,,,,,,, drop = FALSE]
  ref$data = ref$data[2,,,,,,, drop = FALSE]
  return(list(metab, ref))
}

bc <- function(mrs_data, lambda, p) {
  if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
  }
  # extract real part
  mrs_data$data <- Re(mrs_data$data)
  apply_mrs(mrs_data, 7, ptw::baseline.corr, lambda, p)
}

#' @export
`+.mrs_data` <- function(a, b) {
  if ( class(b) == "mrs_data" ) {
    a$data <- a$data + b$data
  } else if (class(b) == "numeric") {
    a$data <- a$data + b
  }
  return(a)
}

#' @export
`-.mrs_data` <- function(a, b = NULL) {
  if ( class(b) == "mrs_data" ) {
    a$data = a$data - b$data
  } else if (is.null(b)) {
    a$data = -a$data
  } else if ( class(b) == "numeric") {
    a$data = a$data - b
  }
  return(a)
}

#' Calculate the mean dynamic data.
#' @param mrs_data Dynamic MRS data.
#' @return Mean dynamic data.
#' @export
mean_dyns <- function(mrs_data) {
  mrs_data$data <- aperm(mrs_data$data, c(5,1,2,3,4,6,7))
  mrs_data$data <- colMeans(mrs_data$data)
  new_dim <- dim(mrs_data$data)
  dim(mrs_data$data) <- c(new_dim[1:4],1,new_dim[5:6])
  mrs_data
}

#' Calculate the mean of adjacent dynamic scans.
#' @param mrs_data Dynamic MRS data.
#' @param block_size Number of adjacent dynamics scans to average over.
#' @return Dynamic data averaged in blocks.
#' @export
mean_dyn_blocks <- function(mrs_data, block_size) {
  
  if ((dyns(mrs_data) %% block_size) != 0) {
    warning("Block size does not fit into the number of dynamics without truncation.")
  }
  
  new_dyns <-  floor(dyns(mrs_data) / block_size)
  mrs_out <- get_dyns(mrs_data, seq(1,new_dyns*block_size, block_size))
  for (n in 2:block_size) {
    mrs_out <- mrs_out + get_dyns(mrs_data, seq(n,new_dyns*block_size, block_size))
  }
  
  mrs_out / block_size
}

#' Calculate the pairwise means across a dynamic data set.
#' @param mrs_data Dynamic MRS data.
#' @return Mean dynamic data of adjacent dynamic pairs.
#' @export
mean_dyn_pairs <- function(mrs_data) {
  pairs <- get_odd_dyns(mrs_data) + get_even_dyns(mrs_data)
  pairs / 2
}

#' Calculate the sum of data dynamics.
#' @param mrs_data Dynamic MRS data.
#' @return Sum of data dynamics.
#' @export
sum_dyns <- function(mrs_data) {
  mrs_data$data <- aperm(mrs_data$data, c(5,1,2,3,4,6,7))
  mrs_data$data <- colSums(mrs_data$data)
  new_dim <- dim(mrs_data$data)
  dim(mrs_data$data) <- c(new_dim[1:4],1,new_dim[5:6])
  mrs_data
}

#' Calculate the sum across receiver coil elements.
#' @param mrs_data MRS data split across receiver coil elements.
#' @return Sum across coil elements.
#' @export
sum_coils <- function(mrs_data) {
  mrs_data$data <- aperm(mrs_data$data, c(6,1,2,3,4,5,7))
  mrs_data$data <- colSums(mrs_data$data)
  new_dim <- dim(mrs_data$data)
  dim(mrs_data$data) <- c(new_dim[1:5],1,new_dim[6])
  mrs_data
}

cplx_median <- function(input) {
  stats::median(Re(input)) + stats::median(Im(input)) * 1i
}

median_dyns <- function(mrs_data) {
  return(apply_mrs(mrs_data, 5, cplx_median))
}

#' @export
`*.mrs_data` <- function(a, b) {
  if ( class(b) == "mrs_data" ) {
    a$data <- a$data * b$data
  } else if ( class(b) == "numeric") {
    a$data <- a$data * b
  }
  return(a)
}

#' @export
`/.mrs_data` <- function(a,b) {
  if ( class(b) == "mrs_data" ) {
    a$data <- a$data / b$data
  } else if (class(b) == "numeric") {
    a$data <- a$data / b
  }
  return(a)
}

# TODO correct first imaginary data point?
recon_imag_vec <- function(data) {
  data <- Conj(hilbert(Re(data)))
  data <- ift_shift(data)
  fh <- data[1:(length(data) / 2)]
  sh <- c(data[(length(data) / 2 + 2):length(data)], 0)
  data <- fh + Conj(rev(sh))
}

td_conv_filt_vec <- function(fid, K = 25, ext = 1)
{
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
  if (is_fd(mrs_data)) {
      mrs_data <- fd2td(mrs_data)
  }
  apply_mrs(mrs_data, 7, td_conv_filt_vec, K, ext)
}

#' HSVD based signal filter.
#' 
#' HSVD based signal filter described in:
#' Barkhuijsen H, de Beer R, van Ormondt D. Improved algorithm for noniterative 
#' and timedomain model fitting to exponentially damped magnetic resonance
#' signals. J Magn Reson 1987;73:553-557.
#' 
#' @param mrs_data MRS data to be filtered.
#' @param xlim Frequency range in Hz to filter.
#' @param comps Number of Lorentzian components to use for modelling.
#' @param propack Option to use PROPACK SVD (logical).
#' @export
hsvd_filt <- function(mrs_data, xlim = c(-30, 30), comps = 50, propack = FALSE) {
# TODO, region defn in ppm and allow multiple region selection
  if (is_fd(mrs_data)) {
      mrs_data <- fd2td(mrs_data)
  }
  apply_mrs(mrs_data, 7, hsvd_filt_vec, fs = fs(mrs_data), region = xlim,
            comps = comps, propack)
}

hsvd_filt_vec <- function(fid, fs, region = c(-30,30), comps = 50, 
                          propack = TRUE) {
  
  hsvd_res <- hsvd(fid, fs, K = comps, propack)  
  idx <- (hsvd_res$reson_table$frequency < region[2]) &
         (hsvd_res$reson_table$frequency > region[1] )
  model <- rowSums(hsvd_res$basis[,idx])
  fid - model
}

# TODO add option to remove -ve dampings
hsvd <- function(y, fs, K = 50, propack = TRUE, fast_hank = TRUE) {
  N <- length(y)
  L <- floor(0.5 * N)
  # M <- N + 1 - L
  
  # scale the input vector to keep things stable
  sc_factor <- max(Mod(y))
  y <- y / sc_factor

  # H is the LxM Hankel LP matrix
  if (fast_hank) {
    H <- matrixcalc::hankel.matrix(L + 1, y)
    H <- H[1:L,]
  } else {
    H <- pracma::hankel(y[1:L], y[L:N])
  }
  
  if (propack)  {
    # K1 formulation
    H_real <- rbind(cbind(Re(H), -Im(H)), cbind(Im(H), Re(H)))
    #H_real <- matrix(, 2*nrow(H), 2*ncol(H))
    #H_real[seq(1, 2*nrow(H), 2), seq(1, 2*ncol(H), 2)] <- Re(H)
    #H_real[seq(2, 2*nrow(H), 2), seq(1, 2*ncol(H), 2)] <- -Im(H)
    #H_real[seq(2, 2*nrow(H), 2), seq(2, 2*ncol(H), 2)] <- Re(H)
    #H_real[seq(1, 2*nrow(H), 2), seq(2, 2*ncol(H), 2)] <- Im(H)
    svd_res <- svd::propack.svd(H_real, neig = (K))
    #svd_res <- svd::trlan.svd(H_real, neig = K)
    #svd_res$u <- svd_res$u[1:(nrow(H)/2),1:(ncol(H)/2)] - 1i*svd_res$u[(nrow(H)/2):nrow(H),1:(ncol(H)/2)]
    #svd_res <- svd(H_real)
    #svd_res$u <- svd_res$u
  } else {
    svd_res <- svd(H)
  }
  
  # construct H of rank K
  Uk <- svd_res$u[,1:K]
  rows <- nrow(Uk)
  Ukt <- Uk[2:rows,]
  Ukb <- Uk[1:(rows - 1),]
  Zp = MASS::ginv(Ukb) %*% Ukt

  # find the poles
  q <- pracma::eig(Zp)
  q <- log(q)
  dt <- 1 / fs
  dampings <- Re(q) / dt;
  frequencies <- Im(q)/(2 * pi) / dt;
  
  t <- seq(from = 0, to = (N - 1) / fs, by = 1 / fs)
  t_mat <- matrix(t, ncol = K, nrow = N)
  
  # TODO not sure if the next line is faster
  #basis <- t(exp(t(t_mat) * (dampings + 2i * pi * frequencies)))
  
  freq_damp <- matrix(dampings + 2i * pi * frequencies, ncol = K, nrow = N,
                      byrow = TRUE)
  
  basis <- exp(t_mat * freq_damp)
  
  ahat <- MASS::ginv(basis) %*% y
  
  # Undo scaling
  ahat <- ahat * sc_factor
  #yhat <- basis%*%ahat 
  
  # scale basis by ahat
  ahat_mat <- matrix(ahat, ncol = K, nrow = N, byrow = TRUE)
  basis <- basis * ahat_mat
  
  # generate a table of resonances
  reson_table <- data.frame(amplitude = Mod(ahat), phase = Arg(ahat) * 180 / pi,
                            frequency = frequencies, damping = dampings)
  
  list(basis = basis, reson_table = reson_table)
}

auto_phase <- function(mrs_data) {
  if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
  }
  apply_mrs(mrs_data, 7, auto_phase_vec)
}

auto_phase_vec <- function(vec) {
  res <- stats::optim(0, phase_obj_fn, gr = NULL, vec, method = "Brent",
                      lower = -180, upper = 180)
  vec * exp(1i * (res$par + 180) / 180 * pi)
}

phase_obj_fn <- function(phi, vec) {
  sum(Re(vec * exp(1i * phi / 180 * pi)))
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

#' Apply eddy current correction using the Klose method:
#' 
#' In vivo proton spectroscopy in presence of eddy currents.
#' Klose U.
#' Magn Reson Med. 1990 Apr;14(1):26-30.
#' 
#' @param metab MRS data to be corrected.
#' @param ref Reference dataset.
#' @return Corrected data in the time domain.
#' @export
ecc <- function(metab, ref) {
  if (is_fd(metab)) {
      metab <- fd2td(metab)
  }
  
  if (is_fd(ref)) {
      ref <- fd2td(ref)
  }
  
  if (dyns(ref) > 1) {
    ref <- mean_dyns(ref)
    warning("Using the mean reference signal for ECC.")
  }
  
  # repeat the refernce signal to match the number of dynamics
  if (dyns(metab) > 1) {
    ref <- rep_dyn(ref, dyns(metab))
  }
  
  mrs_data <- combine_metab_ref(metab, ref)
  ecc_data <- apply_mrs(mrs_data, c(1,7), ecc_2d_array)
  get_metab(ecc_data)
}

#' Apodise MRSI data in the x-y direction with a k-space hamming filter.
#' @param mrs_data MRSI data.
#' @return Apodised data.
#' @export
apodise_xy <- function(mrs_data) {
  # TODO check data is 2D in xy dirn and make faster...
  
  # put xy dims into k-space
  mrs_data <- apply_mrs(mrs_data, 2, ft_shift)
  mrs_data <- apply_mrs(mrs_data, 3, ft_shift)
  # apply filter
  mrs_data <- apply_mrs(mrs_data, 2, hamming_vec)
  mrs_data <- apply_mrs(mrs_data, 3, hamming_vec)
  # put xy dims back to space
  mrs_data <- apply_mrs(mrs_data, 2, ift_shift)
  apply_mrs(mrs_data, 3, ift_shift)
}

#' Zero-fill MRSI data in the k-space x-y direction.
#' @param mrs_data MRSI data.
#' @param factor Zero-filling factor, factor of 2 returns a dataset with
#' twice the original points in the x-y directions.
#' @return Zero-filled data.
#' @export
zf_xy <- function(mrs_data, factor = 2) {
  # TODO check data is 2D in xy dirn and make (much) faster by using afill
  # TODO check this works for even numbers of rows and cols...
  
  # put xy dims into k-space
  mrs_data <- apply_mrs(mrs_data, 2, ft_shift)
  mrs_data <- apply_mrs(mrs_data, 3, ft_shift)
  # apply filter
  mrs_data <- apply_mrs(mrs_data, 2, zp_vec, dim(mrs_data)[2] * factor)
  mrs_data <- apply_mrs(mrs_data, 3, zp_vec, dim(mrs_data)[3] * factor)
  # put xy dims back to space
  mrs_data <- apply_mrs(mrs_data, 2, ift_shift)
  apply_mrs(mrs_data, 3, ift_shift)
}

hamming_vec <- function(vector) {
  vector*signal::hamming(length(vector))
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
#' @param metab_mrs MRS data containing metabolite data.
#' @param ref_mrs MRS data containing reference data (optional).
#' @param noise_mrs MRS data from a noise scan (optional).
#' @param scale option to rescale coil elements based on the first data point
#' (logical).
#' @param sum_coils sum the coil elements as a final step (logical).
#' @return MRS data.
#' @export
comb_coils <- function(metab_mrs, ref_mrs = NULL, noise_mrs = NULL, 
                       scale = TRUE, sum_coils = TRUE) {
  
  metab_only <- FALSE
  if (is.null(ref_mrs)) {
    ref_mrs <- mean_dyns(metab_mrs)
    metab_only <- TRUE
  }
  
  if (is_fd(metab_mrs)) {
      metab <- fd2td(metab)
  }
  
  if (is_fd(ref_mrs)) {
      ref <- fd2td(ref)
  }
  
  # get the first dynamic of the ref data
  # first_ref <- get_dyns(ref_mrs, 1)
  # fp <- get_fp(first_ref)
  
  # get the dynamic mean of the ref data
  mean_ref <- mean_dyns(ref_mrs)
  fp <- get_fp(mean_ref)
  
  phi <- Arg(fp)
  amp <- Mod(fp)
  
  if (!is.null(noise_mrs)) {
    amp <- amp / (calc_coil_noise_sd(noise_mrs) ^ 2)
  }
  
  # phase and scale ref data
  ref_dims <- dim(ref_mrs$data)
  
  ang <- rep(phi, prod(ref_dims[-6]))
  dim(ang) <- c(ref_dims[c(6, 2, 3, 4, 5, 1, 7)])
  ang <- aperm(ang, c(6, 2, 3, 4, 5, 1, 7))
  
  if (scale) {
    scale_f <- rep(amp, prod(ref_dims[-6]))
    dim(scale_f) <- c(ref_dims[c(6, 2, 3, 4, 5, 1, 7)])
    scale_f <- aperm(scale_f, c(6, 2, 3, 4, 5, 1, 7))
    
    ref_mrs_ps <- ref_mrs
    ref_mrs_ps$data <- ref_mrs$data * exp(-1i * ang) * scale_f
  } else {
    ref_mrs_ps <- ref_mrs
    ref_mrs_ps$data <- ref_mrs$data * exp(-1i * ang)
  }
  
  ref_mrs_ps <- sum_coils(ref_mrs_ps)
  
  # phase and scale metab data
  metab_dims <- dim(metab_mrs$data)
  
  ang <- rep(phi, prod(metab_dims[-6]))
  dim(ang) <- c(metab_dims[c(6, 2, 3, 4, 5, 1, 7)])
  ang <- aperm(ang, c(6, 2, 3, 4, 5, 1, 7))
  
  if (scale) {
    scale_f <- rep(amp, prod(metab_dims[-6]))
    dim(scale_f) <- c(metab_dims[c(6, 2, 3, 4, 5, 1, 7)])
    scale_f <- aperm(scale_f, c(6, 2, 3, 4, 5, 1, 7))
    
    metab_mrs_ps <- metab_mrs
    metab_mrs_ps$data <- metab_mrs$data * exp(-1i * ang) * scale_f
  } else {
    metab_mrs_ps <- metab_mrs
    metab_mrs_ps$data <- metab_mrs$data * exp(-1i * ang)
  }
  
  if (sum_coils == TRUE) {
    metab_mrs_ps <- sum_coils(metab_mrs_ps)
  }
  
  if (metab_only) {
    return(metab_mrs_ps)
  } else {
    return(list(metab = metab_mrs_ps, ref = ref_mrs_ps))
  }
}

#' Replicate a scan in the dynamic dimension.
#' @param mrs_data MRS data to be replicated.
#' @param times Number of times to replicate.
#' @return Replicated data object.
#' @export
rep_dyn <- function(mrs_data, times) {
  orig_dim <- dim(mrs_data$data)
  new_dim <- orig_dim
  new_dim[5] <- new_dim[5] * times
  # make the dynamic dimension (5th) the last
  rep_data <- aperm(mrs_data$data, c(1,2,3,4,6,7,5))
  # duplicate the data
  rep_data <- rep(rep_data, times)
  # set the new dimesnions
  dim(rep_data) <- new_dim[c(1,2,3,4,6,7,5)]
  # reorder
  rep_data <- aperm(rep_data, c(1,2,3,4,7,5,6))
  mrs_data$data <- rep_data
  mrs_data
}

#' Estimate the standard deviation of the noise from a segment of an mrs_data object.
#' @param mrs_data MRS data object.
#' @param n The number of data points (taken from the end of array) to use in the estimation.
#' @param offset The number of final points to exclude from the calculation.
#' @param p_order Polynomial order to fit to the data before estimating the standard deviation.
#' @return Standard deviation array.
#' @export
est_noise_sd <- function(mrs_data, n = 100, offset = 100, p_order = 2) {
  apply_mrs(mrs_data, 7, est_noise_sd_vec, n, offset, p_order, data_only = TRUE)
}

est_noise_sd_vec <- function(x, n = 100, offset = 100, p_order = 2) {
  N <- length(x)
  seg <- Re(x[(N - offset - n + 1):(N - offset)])
  lm_res <- stats::lm(seg ~ poly(1:n, p_order))
  stats::sd(lm_res$residual)
}

#' Calculate the noise correlation between coil elements.
#' @param noise_data \code{mrs_data} object with one FID for each coil element.
#' @return Correlation matrix.
#' @export
calc_coil_noise_cor <- function(noise_data) {
  cplx_data <- drop(noise_data$data)
  # concat real and imag parts
  real_data <- cbind(Re(cplx_data), Im(cplx_data))
  stats::cor(t(real_data))
}

#' Calculate the noise stanard deviation for each coil element.
#' @param noise_data \code{mrs_data} object with one FID for each coil element.
#' @return array of standard deviations.
#' @export
calc_coil_noise_sd <- function(noise_data) {
  cplx_data <- drop(noise_data$data)
  # concat real and imag parts
  real_data <- cbind(Re(cplx_data), Im(cplx_data))
  apply(real_data, 1, stats::sd)
}