#' Robust Alignment to a Target Spectrum (RATS).
#' 
#' @param mrs_data MRS data to be corrected.
#' @param ref optional MRS data to use as a reference, the mean of all dynamics
#' is used if this argument is not supplied.
#' @param xlim optional frequency range to perform optimisation, set to NULL
#' to use the full range.
#' @param max_shift maximum allowable frequency shift in Hz.
#' @param p_deg polynomial degree used for baseline modelling. Negative values
#' disable baseline modelling.
#' @param sp_N number of spline functions, note the true number will be sp_N +
#' sp_deg.
#' @param sp_deg degree of spline functions.
#' @param max_t truncate the FID when longer than max_t to reduce time taken,
#' set to NULL to use the entire FID.
#' @param basis_type may be one of "poly" or "spline".
#' @param rescale_output rescale the bl_matched_spec and bl output to improve
#' consistency between dynamic scans.
#' @param phase_corr apply phase correction (in addition to frequency). TRUE by
#' default.
#' @param ret_corr_only return the corrected mrs_data object only.
#' @return a list containing the corrected data; phase and shift values in units
#' of degrees and Hz respectively.
#' @export
rats <- function(mrs_data, ref = NULL, xlim = c(4, 0.5), max_shift = 20,
                 p_deg = 2, sp_N = 2, sp_deg = 3, max_t = 0.2,
                 basis_type = "poly", rescale_output = TRUE,
                 phase_corr = TRUE, ret_corr_only = TRUE) {
  
  if (inherits(mrs_data, "list")) {
    res <- lapply(mrs_data, rats, ref = ref, xlim = xlim, max_shift = max_shift,
                  p_deg = p_deg, sp_N = sp_N, sp_deg = sp_deg, max_t = max_t,
                  basis_type = basis_type, rescale_output = rescale_output,
                  phase_corr = phase_corr, ret_corr_only = ret_corr_only) 
    return(res)
  }
  
  # move mrs_data to the time-domain
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
   
  t <- seconds(mrs_data)
  
  # truncate the FID to improve speed 
  mrs_data_mod <- mrs_data
  if (!is.null(max_t)) {
    pts <- sum(t < max_t)
    t <- t[1:pts]
    mrs_data_mod$data <- mrs_data_mod$data[,,,,,,1:pts, drop = FALSE]
  }
  
  # align to mean if ref is not given
  if (is.null(ref)) ref <- mean(mrs_data, na.rm = TRUE)
  
  if (!is.null(max_t)) {
    # move ref to the time-domain
    if (is_fd(ref)) ref <- fd2td(ref)
    ref$data <- ref$data[,,,,,,1:pts, drop = FALSE]
  }
  
  # move ref back to the freq-domain
  if (!is_fd(ref)) ref <- td2fd(ref)
  
  # ref_mod is in the fd 
  ref_mod <- crop_spec(ref, xlim)
  ref_data <- as.complex(ref_mod$data)
  inds <- get_seg_ind(ppm(mrs_data_mod), xlim[1], xlim[2]) 
  
  if (basis_type == "poly") {
    if (p_deg == 0) {
      basis <- rep(1, length(inds))
    } else if (p_deg > 0) {
      basis <- cbind(rep(1, length(inds)), stats::polym(1:length(inds),
                                                        degree = p_deg))
    } else {
      basis <- NULL
    }
  } else if (basis_type == "spline") {
    basis <- bbase(length(inds), sp_N, sp_deg)
  } else{
    stop("I don't belong here.")
  }
  
  # optimisation step
  res <- apply_mrs(mrs_data_mod, 7, optim_rats, ref_data, t, inds, basis, 
                   max_shift, data_only = TRUE)
  
  phases <- Re(res[,,,,,,1, drop = FALSE])
  shifts <- Re(res[,,,,,,2, drop = FALSE])
  amps   <- Re(res[,,,,,,3, drop = FALSE])
  
  corr_spec <- ref_mod
  corr_spec$data <- res[,,,,,,4:(length(inds) + 3), drop = FALSE]
  bl_spec <- ref_mod
  bl_spec$data <- res[,,,,,,(length(inds) + 4):(2 * length(inds) + 3),
                      drop = FALSE]
  
  # apply to original data
  t_orig <- rep(seconds(mrs_data), each = Nspec(mrs_data))
  t_array <- array(t_orig, dim = dim(mrs_data$data))
  shift_array <- array(shifts, dim = dim(mrs_data$data))
  
  if (!phase_corr) phases <- 0
  
  phase_array <- array(phases, dim = dim(mrs_data$data))
  mod_array <- exp(2i * pi * t_array * shift_array + 1i * phase_array * pi / 180)
  mrs_data$data <- mrs_data$data * mod_array
 
  # maintain original intensities for bl_matched_spec and bl output
  if (!rescale_output) {
    corr_spec <- scale_mrs_amp(corr_spec, 1 / amps)
    bl_spec   <- scale_mrs_amp(bl_spec,   1 / amps)
  }
  
  # results
  res <- list(corrected = mrs_data, phases = -phases, shifts = -shifts,
              amps = amps, bl_matched_spec = corr_spec, bl = -bl_spec)
  
  if (ret_corr_only) {
    return(res$corrected)
  } else {
    return(res)
  }
}

optim_rats <- function(x, ref, t, inds, basis, max_shift) {
  
  # masked spectra are special case
  if (is.na(x[1])) {
    res <- c(NA, NA, NA, rep(NA, 2 * length(inds)))
    return(res)
  }
  
  # optim step
  res <- stats::optim(c(0), rats_obj_fn, gr = NULL, x, ref, t, inds, basis, 
               method = "Brent", lower = -max_shift, upper = max_shift)
  
  # find the phase
  shift <- res$par[1]
  x <- x * exp(2i * pi * shift * t)
  x <- ft_shift(x)
  x <- x[inds]
  
  if (is.null(basis)) {
    basis_mod <- x
    ahat <- unname(qr.solve(basis_mod, ref))
    yhat <- basis_mod * ahat
  } else {
    basis_mod <- cbind(x, basis)
    ahat <- unname(qr.solve(basis_mod, ref))
    yhat <- basis_mod %*% ahat
    bl   <- basis_mod %*% c(0, ahat[2:length(ahat)])
  }
  
  res <- c(Arg(ahat[1]) * 180 / pi, res$par, Mod(ahat[1]), yhat, bl)
  return(res)
}

rats_obj_fn <- function(par, x, ref, t, inds, basis) {
  shift <- par[1]
  x <- x * exp(2i * pi * shift * t)
  x <- ft_shift(x)
  x <- x[inds]
  
  if (is.null(basis)) {
    basis_mod <- x
  } else {
    basis_mod <- cbind(x, basis)
  }
  
  # use ginv
  #inv_basis <- ginv(basis)
  #ahat <- inv_basis %*% ref
  
  # use qr
  ahat <- qr.solve(basis_mod, ref)
  
  if (is.null(basis)) {
    yhat <- basis_mod * ahat
  } else {
    yhat <- basis_mod %*% ahat
  }
  
  res <- c(Re(yhat), Im(yhat)) - c(Re(ref), Im(ref))
  
  sum(res ^ 2)
}

quick_phase_ref <- function(mrs_data) {
  ref <- sim_resonances(freq = c(2.01, 3.03, 3.22),
                        acq_paras = get_acq_paras(mrs_data))
  
  # correct the first dynamic
  res <- rats(get_dyns(mrs_data,1), ref)
  
  # apply to full dataset
  mrs_data <- phase(mrs_data, res$phases[1])
  mrs_data <- shift(mrs_data, res$shifts[1], units = "hz")
  mrs_data
}
