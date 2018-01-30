#' Frequency and phase corretion of MRS data in the time-domain.
#' 
#' An implementation of the method published by Near et al MRM 73:44-50 (2015).
#' 
#' @param mrs_data MRS data to be corrected
#' @param ref optional MRS data to use as a reference, the first dynamic of 
#' mrs_data is used if this argument is not supplied.
#' @param xlim optional frequency range to perform optimsiation.
#' @param max_t truncate the FID when longer than max_t to reduce time taken
#' @return a list containing the corrected data; phase and shift values in units
#' of degrees and Hz respectivly.
#' @export
fp_corr_td <- function(mrs_data, ref = NULL, xlim = NULL, max_t = 0.2) {
  # align to first dynamic if ref is not given
  if (is.null(ref)) ref <- get_dyns(mrs_data, 1)
  
  if (is.null(xlim)) {
    if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
    if (is_fd(ref)) ref <- fd2td(ref)
    ref_mod <- ref
    mrs_data_mod <- mrs_data
  } else {
    if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
    if (!is_fd(ref)) ref <- td2fd(ref)
    mrs_data <- fd2td(mrs_data)
    mrs_data_mod <- fd2td(crop_spec(mrs_data, xlim))
    ref_mod <- fd2td(crop_spec(ref, xlim))
  }
  
  # truncate the FID to improve speed 
  t <- seconds(ref_mod)
  pts <- sum(t < max_t)
  t <- t[1:pts]
  ref_data <- as.complex(ref_mod$data)
  ref_data <- ref_data[1:pts]
  mrs_data_mod$data <- mrs_data_mod$data[,,,,,,1:pts, drop = FALSE]
  
  # optimisation step
  res <- apply_mrs(mrs_data_mod, 7, optim_fp, ref_data, t, data_only = TRUE)
  phases <- res[,,,,,,1, drop = FALSE]
  shifts <- res[,,,,,,2, drop = FALSE]
  
  # apply to original data
  t_orig <- rep(seconds(mrs_data), each = Nspec(mrs_data))
  t_array <- array(t_orig, dim = dim(mrs_data$data))
  shift_array <- array(shifts, dim = dim(mrs_data$data))
  phase_array <- array(phases, dim = dim(mrs_data$data))
  mod_array <- exp(2i * pi * t_array * shift_array + 1i * phase_array * pi / 180)
  mrs_data$data <- mrs_data$data * mod_array
 
  # results 
  list(corrected = mrs_data, phases = phases, shifts = shifts)
}

phase_drift_obj_fn_td <- function(par, x, ref, t) {
  phase <- par[1]
  shift <- par[2]
  x <- x * exp(1i * (phase / 180 * pi + 2 * pi * shift * t))
  c(Re(x), Im(x)) - c(Re(ref), Im(ref))
}

optim_fp <- function(x, ref, t) {
  ctrl <- minpack.lm::nls.lm.control()
  ctrl$maxiter <- 100
  res <- minpack.lm::nls.lm(c(0, 0), NULL, NULL, phase_drift_obj_fn_td, 
                            jac = NULL, control = ctrl,
                            x, ref, t)
  
  res$par
}
