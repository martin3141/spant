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
#' @param max_t truncate the FID when longer than max_t to reduce time taken
#' @return a list containing the corrected data; phase and shift values in units
#' of degrees and Hz respectively.
#' @export
rats <- function(mrs_data, ref = NULL, xlim = c(4, 0.5), max_shift = 20,
                 p_deg = 2, max_t = 0.2) {
  
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
  if (is.null(ref)) ref <- mean(mrs_data)
  
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
  
  if (p_deg == 0) {
    poly_basis <- rep(1, length(inds))
  } else if (p_deg > 0) {
    poly_basis <- cbind(rep(1, length(inds)), stats::polym(1:length(inds),
                                                           degree = p_deg))
  } else {
    poly_basis <- NULL
  }
  
  # optimisation step
  res <- apply_mrs(mrs_data_mod, 7, optim_rats, ref_data, t, inds, poly_basis, 
                   max_shift, data_only = TRUE)
  
  phases <- Re(res[,,,,,,1, drop = FALSE])
  shifts <- Re(res[,,,,,,2, drop = FALSE])
  amps   <- Re(res[,,,,,,3, drop = FALSE])
  
  corr_spec <- ref_mod
  corr_spec$data <- res[,,,,,,4:(length(inds) + 2), drop = FALSE]
  
  # apply to original data
  t_orig <- rep(seconds(mrs_data), each = Nspec(mrs_data))
  t_array <- array(t_orig, dim = dim(mrs_data$data))
  shift_array <- array(shifts, dim = dim(mrs_data$data))
  phase_array <- array(phases, dim = dim(mrs_data$data))
  mod_array <- exp(2i * pi * t_array * shift_array + 1i * phase_array * pi / 180)
  mrs_data$data <- mrs_data$data * mod_array
 
  # results 
  list(corrected = mrs_data, phases = -phases, shifts = -shifts, amps = amps,
       bl_matched_spec = corr_spec)
}

optim_rats <- function(x, ref, t, inds, poly_basis, max_shift) {
  # optim step
  res <- stats::optim(c(0), rats_obj_fn, gr = NULL, x, ref, t, inds, poly_basis, 
               method = "Brent", lower = -max_shift, upper = max_shift)
  
  # find the phase
  shift <- res$par[1]
  x <- x * exp(2i * pi * shift * t)
  x <- ft_shift(x)
  x <- x[inds]
  
  if (is.null(poly_basis)) {
    basis <- x
    ahat <- unname(qr.solve(basis, ref))
    yhat <- basis * ahat
  } else {
    basis <- cbind(x, poly_basis)
    ahat <- unname(qr.solve(basis, ref))
    yhat <- basis %*% ahat
  }
  
  c(Arg(ahat[1]) * 180 / pi, res$par, Mod(ahat[1]), yhat)
}

rats_obj_fn <- function(par, x, ref, t, inds, poly_basis) {
  shift <- par[1]
  x <- x * exp(2i * pi * shift * t)
  x <- ft_shift(x)
  x <- x[inds]
  
  if (is.null(poly_basis)) {
    basis <- x
  } else {
    basis <- cbind(x, poly_basis)
  }
  
  # use ginv
  #inv_basis <- MASS::ginv(basis)
  #ahat <- inv_basis %*% ref
  
  # use qr
  ahat <- qr.solve(basis, ref)
  
  if (is.null(poly_basis)) {
    yhat <- basis * ahat
  } else {
    yhat <- basis %*% ahat
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
