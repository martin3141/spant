# ABFit (Adaptive Baseline Fitting)
abfit <- function(y, acq_paras, basis, opts = NULL) {
  
  # has this data been masked?
  if (is.na(y[1])) return(list(amps = NA, crlbs = NA, diags = NA, fit = NA))
  
  #### init ####
  # generate spline basis and prepare objects for the fitting routine
  
  # use default fitting opts if not specified 
  if (is.null(opts)) opts <- abfit_opts()
  
  # check the noise_region is within the spectral range
  ppm_sc <- -hz(fs = acq_paras$fs, N = acq_paras$N) / acq_paras$ft * 1e6 +
             acq_paras$ref
  ppm_lim <- sort(range(ppm_sc))
  noise_reg_lim <- sort(opts$noise_region)
  if ((noise_reg_lim[1] < ppm_lim[1]) | (noise_reg_lim[2] > ppm_lim[2])) {
    stop(paste0("The spectral range for the noise region estimate is outside",
                " the acquired spectral width. Change the noise_region",
                " parameter in the fitting options."))
  }
  
  # check the basis has the correct number of data points
  if (length(y) != nrow(basis$data)) {
    stop("Data points in the basis and acquired data do not match.")
  }
 
  # zero pad input to twice length
  if (opts$zp) y <- c(y, rep(0, length(y))) 
  
  # convert y vec to mrs_data object to use convenience functions
  mrs_data <- vec2mrs_data(y, fs = acq_paras$fs, ft = acq_paras$ft, 
                           ref = acq_paras$ref)
  
  #### 1 coarse freq align ####
  if (opts$pre_align) {
    max_init_shift_hz <- acq_paras$ft * 1e-6 * opts$max_pre_align_shift
    init_shift <- as.numeric(align(mrs_data, opts$pre_align_ref_freqs,
                                   max_shift = max_init_shift_hz, ret_df = TRUE,
                                   lb = 0)$shifts)
  } else {
    init_shift <- 0
  }
  
  # time scale in seconds
  t <- seconds(mrs_data)
  
  # freq scale in Hz
  f <- seq(from = -acq_paras$fs / 2,
           to = acq_paras$fs / 2 - acq_paras$fs / length(y),
           length.out = length(y))
  
  Nbasis <- ncol(basis$data)
  
  # tranform metab basis to time-domain
  raw_metab_basis <- apply(basis$data, 2, ift_shift)
  
  # zero pad metab basis matrix to twice length
  if (opts$zp) {
    zero_mat <- matrix(0, nrow = nrow(raw_metab_basis),
                       ncol = ncol(raw_metab_basis))
    
    raw_metab_basis <- rbind(raw_metab_basis, zero_mat)
  }
  
  # set phase, damping and shift limits for the pre-fit
  par   <- c(0, opts$init_damping, init_shift)
  upper <- c(opts$max_phase * pi / 180, opts$max_damping,
             init_shift + opts$max_shift * acq_paras$ft * 1e-6)
  
  lower <- c(-opts$max_phase * pi / 180, 0,
             init_shift - opts$max_shift * acq_paras$ft * 1e-6)
  
  #### 2 approx iter fit ####
  # 3 para pre-fit with flexable bl and no broad signals in basis
  # to get good starting values
  if (opts$maxiters_pre > 0) {
    
    # generate spline basis
    sp_bas_pf <- generate_sp_basis(mrs_data, opts$pre_fit_ppm_right,
                                   opts$pre_fit_ppm_left,
                                   opts$bl_comps_pppm)
    
    # 10 ED per ppm is a good start for 3T 1H brain MRS
    pre_fit_ed_pppm <- opts$pre_fit_bl_ed_pppm
    
    # calculate the right lambda
    lambda <- calc_lambda_from_ed(sp_bas_pf$bl_bas, sp_bas_pf$deriv_mat,
                                  pre_fit_ed_pppm * sp_bas_pf$ppm_range)
    
    # augment spline basis with penalty matrix
    bl_basis_pre_fit <- rbind(sp_bas_pf$bl_bas, sp_bas_pf$deriv_mat *
                             (lambda ^ 0.5))
    
    # remove broad signal components for the prefit
    if (opts$remove_lip_mm_prefit) {
      broad_indices <- c(grep("^Lip", basis$names), grep("^MM", basis$names))
      metab_basis_pre <- raw_metab_basis[, -broad_indices]
    } else {
      metab_basis_pre <- raw_metab_basis
    }
    
    # Do a 1D search to improve the starting value of the phase estimate.
    # Note, this is not part of the published method, but was added in Jan 2021
    # to avoid issues with the simplex "prefit" getting stuck in a local minima
    # in rare cases.
    if (opts$prefit_phase_search) {
      # 5 degree increments between -175 to 180
      phi_zero  <- seq(-pi + 2 * pi / 72, pi, length.out = 72)
      par_frame <- data.frame(phi_zero, opts$init_damping, init_shift)
      
      phase_res <- apply(par_frame, 1, abfit_3p_obj, y = y,
                         raw_metab_basis = metab_basis_pre,
                         bl_basis = bl_basis_pre_fit, t = t, f = f,
                         inds = sp_bas_pf$inds, bl_comps = sp_bas_pf$bl_comps,
                         sum_sq = TRUE, ret_full = FALSE,
                         ahat_calc_method = opts$ahat_calc_method)
      
      # find the optimum value and update the starting value for the next step
      par[1] <- phi_zero[which.min(phase_res)]
    }
    
    if (opts$algo_pre == "levmar") {
      # get the default nnls control paras
      ctrl <- minpack.lm::nls.lm.control()
      ctrl$maxiter <- opts$maxiters_pre
      
      prelim_res <- minpack.lm::nls.lm(par, lower, upper, abfit_3p_obj, NULL,
                                       ctrl, y, metab_basis_pre,
                                       bl_basis_pre_fit, t, f, sp_bas_pf$inds,
                                       sp_bas_pf$bl_comps, FALSE, FALSE,
                                       opts$ahat_calc_method)
      
      res        <- prelim_res
    } else {
      nloptr_opts <- list("algorithm" = opts$algo_pre,
                          "maxeval" = opts$maxiters_pre)
      
      prelim_res <- nloptr::nloptr(x0 = par, eval_f = abfit_3p_obj,
                                   eval_grad_f = NULL, lb = lower, ub = upper,
                                   opts = nloptr_opts, y = y,
                                   raw_metab_basis = metab_basis_pre,
                                   bl_basis = bl_basis_pre_fit, t = t, f = f,
                                   inds = sp_bas_pf$inds,
                                   bl_comps = sp_bas_pf$bl_comps, sum_sq = TRUE,
                                   ret_full = FALSE,
                                   ahat_calc_method = opts$ahat_calc_method)
      
      res          <- prelim_res
      res$par      <- prelim_res$solution
      res$deviance <- prelim_res$objective
      res$niter    <- prelim_res$iterations
      res$info     <- prelim_res$status
    }
  } else {
    # starting values for the full fit if prefit hasn't been done
    res <- NULL
    res$par <- c(par[1], par[2], par[3])
  }
  
  #### 3 est bl smoothness ####
  if (opts$auto_bl_flex) {
    
    # generate spline basis
    sp_bas_ab <- generate_sp_basis(mrs_data, opts$ppm_right, opts$ppm_left,
                                   opts$bl_comps_pppm)
    
    # array of bl flex values to try
    bl_n <- opts$auto_bl_flex_n
    
    # calculate the right lambda
    lambda_start <- calc_lambda_from_ed(sp_bas_ab$bl_bas, sp_bas_ab$deriv_mat,
                                        opts$max_bl_ed_pppm *
                                        sp_bas_ab$ppm_range)
    
    if (is.null(opts$min_bl_ed_pppm)) {
      lambda_end   <- calc_lambda_from_ed(sp_bas_ab$bl_bas, sp_bas_ab$deriv_mat,
                                          2.01)
    } else {
      lambda_end   <- calc_lambda_from_ed(sp_bas_ab$bl_bas, sp_bas_ab$deriv_mat,
                                          opts$min_bl_ed_pppm *
                                          sp_bas_ab$ppm_range)
    }
    
    lambda_vec   <- 10 ^ (seq(log10(lambda_start), log10(lambda_end),
                              length.out = bl_n))
    
    optim_model_vec <- rep(NULL, bl_n)
    ed_vec          <- rep(NULL, bl_n)
    
    # loop over candidate bl fwhm's
    for (n in 1:bl_n) {
      ed_vec[n] <- calc_ed_from_lambda(sp_bas_ab$bl_bas, sp_bas_ab$deriv_mat,
                                       lambda_vec[n])
      
      # augment spline basis with penalty matrix
      bl_basis_ab <- rbind(sp_bas_ab$bl_bas, sp_bas_ab$deriv_mat *
                     (lambda_vec[n] ^ 0.5))
      
      # apply par to metabolite basis and find the fit residual
      ab_res <- abfit_3p_obj(res$par, y = y, raw_metab_basis = raw_metab_basis,
                             bl_basis = bl_basis_ab, t = t, f = f,
                             inds = sp_bas_ab$inds,
                             bl_comps = sp_bas_ab$bl_comps, sum_sq = FALSE,
                             ret_full = TRUE, opts$ahat_calc_method)
      
      resid_spec <- ab_res$resid
      
      if (opts$optimal_smooth_criterion == "maic") {
        optim_model_vec[n] <- log(sum(resid_spec ^ 2)) + 
                              opts$aic_smoothing_factor * 2 * ed_vec[n] / 
                              length(sp_bas_ab$inds)
      } else if (opts$optimal_smooth_criterion == "aic") {
        optim_model_vec[n] <- log(sum(resid_spec ^ 2)) + 2 * ed_vec[n] / 
                              length(sp_bas_ab$inds)
      } else if (opts$optimal_smooth_criterion == "bic") {
        optim_model_vec[n] <- length(sp_bas_ab$inds) * 
                              log(sum(resid_spec ^ 2)) + ed_vec[n] *
                              log(length(sp_bas_ab$inds))
      } else if (opts$optimal_smooth_criterion == "cv") {
        full_hat_mat <- sp_bas_ab$bl_bas %*% pracma::inv(t(sp_bas_ab$bl_bas) %*%
                        sp_bas_ab$bl_bas + lambda_vec[n] *
                        t(sp_bas_ab$deriv_mat) %*% sp_bas_ab$deriv_mat) %*%
                        t(sp_bas_ab$bl_bas)
        optim_model_vec[n] <- sum((resid_spec / (1 - diag(full_hat_mat))) ^ 2)
      } else {
        stop("unknown optimal smoothing criterion method")
      }
    }
    
    # find the optimal flexability 
    optim_idx <- which.min(optim_model_vec)
    
    if (optim_idx == 1) {
      # the most flexible baseline available was used - could indicate a problem
      max_bl_flex_used <- "TRUE"
    } else {
      max_bl_flex_used <- "FALSE"
    }
    
    opts$bl_ed_pppm <- ed_vec[optim_idx] / sp_bas_ab$ppm_range
    lambda          <- lambda_vec[optim_idx]
    
    resid_vec_frame <- data.frame(t(optim_model_vec))
    colnames(resid_vec_frame) <- paste("auto_bl_crit_", 
                                       as.character(round(ed_vec / 
                                             sp_bas_ab$ppm_range, 3)), sep = "")
    
    resid_vec_frame <- cbind(resid_vec_frame)
  }
  
  #### 4 detailed fit ####
  if (opts$maxiters > 0) {
    
    # estimate the required spine functions based on the auto bl flex
    if (is.null(opts$bl_comps_pppm)) {
      opts$bl_comps_pppm <- opts$bl_ed_pppm * 1.25
    }
    
    # generate the spline basis
    sp_bas_full <- generate_sp_basis(mrs_data, opts$ppm_right, opts$ppm_left,
                                     opts$bl_comps_pppm)
    
    if (is.null(opts$lambda)) {
      lambda <- calc_lambda_from_ed(sp_bas_full$bl_bas, sp_bas_full$deriv_mat,
                                    opts$bl_ed_pppm * sp_bas_full$ppm_range)
    } else {
      lambda <- opts$lambda
    }
    
    # augment spline basis with penalty matrix
    bl_basis_full <- rbind(sp_bas_full$bl_bas, sp_bas_full$deriv_mat *
                     (lambda ^ 0.5))
  
    # set phase, damping, shift, asym, basis_shift and basis_damping limits
    
    lb_init   <- 0.001 # make small but not zero to avoid being equal to the
                       # lower limit
    
    asym_init <- 0
    
    par   <- c(res$par[1], res$par[2], res$par[3], asym_init, rep(0, Nbasis),
               rep(lb_init, Nbasis))
   
    # find any signals with names starting with Lip or MM as they may have
    # different parameter limits
    broad_indices <- c(grep("^Lip", basis$names), grep("^MM", basis$names))
    
    max_basis_shifts <- rep(opts$max_basis_shift, Nbasis) * acq_paras$ft * 1e-6
    
    max_basis_shifts[broad_indices] <- opts$max_basis_shift_broad * 
                                       acq_paras$ft * 1e-6
    
    max_basis_dampings <- rep(opts$max_basis_damping, Nbasis)
    max_basis_dampings[broad_indices] <- opts$max_basis_damping_broad
    
    upper <- c(opts$max_phase * pi / 180, opts$max_damping,
               res$par[3] + opts$max_shift, opts$max_asym,
               max_basis_shifts, max_basis_dampings)
    
    lower <- c(-opts$max_phase * pi / 180, 0, res$par[3] -opts$max_shift,
               -opts$max_asym, -max_basis_shifts, rep(0, Nbasis))
    
    if (opts$phi1_optim) {
      par   <- append(par,    opts$phi1_init, 4)
      upper <- append(upper,  opts$phi1_init + opts$max_dphi1, 4)
      lower <- append(lower,  opts$phi1_init - opts$max_dphi1, 4)
    }
    
    # get the default nnls control paras
    ctrl <- minpack.lm::nls.lm.control()
    ctrl$maxiter <- opts$maxiters
    
    if (opts$anal_jac) {
      jac_fn <- abfit_full_anal_jac
    } else {
      jac_fn <- abfit_full_num_jac
    }
   
    if (is.def(opts$freq_reg)) { 
      # apply best guess for phase parameter to data
      y_mod_noise_est <- y * exp(1i * res$par[1])
      
      # apply best guess shift parameter to data
      y_mod_noise_est <- y_mod_noise_est * exp(2i * pi * t * res$par[3])
      mrs_data_corr_noise_est <- vec2mrs_data(y_mod_noise_est,
                                              fs = acq_paras$fs,
                                              ft = acq_paras$ft, 
                                              ref = acq_paras$ref)
      
      # estimate the noise sd 
      noise_sd_est <- as.numeric(calc_spec_snr(mrs_data_corr_noise_est,
                                               noise_region = opts$noise_region,
                                               full_output = TRUE)$noise_sd)
      # convert ppm sd to Hz
      freq_reg_scaled <- noise_sd_est / (opts$freq_reg * acq_paras$ft * 1e-6)
    } else {
      freq_reg_scaled <- NULL
    }
    
    res <- minpack.lm::nls.lm(par, lower, upper, abfit_full_obj, jac_fn,
                              ctrl, y, raw_metab_basis, bl_basis_full, t,
                              f, sp_bas_full$inds, sp_bas_full$bl_comps, FALSE,
                              NULL, opts$phi1_optim, opts$ahat_calc_method,
                              freq_reg_scaled)
  } 
  
  if (opts$maxiters == 0) {
    res$par <- c(res$par[1], res$par[2], res$par[3], 0, rep(0, Nbasis),
                 rep(0, Nbasis))
    
    if (opts$maxiters_pre == 0) { # don't overwrite if a prefit was done
      res$deviance <- NA
      res$niter <- NA
      res$info <- NA
      res$message <- NA
    }
  }
  
  if (is.null(opts$lambda)) {
    sp_bas_temp <- generate_sp_basis(mrs_data, opts$ppm_right, opts$ppm_left,
                                     opts$bl_comps_pppm)
    
    lambda <- calc_lambda_from_ed(sp_bas_temp$bl_bas, sp_bas_temp$deriv_mat,
                                  opts$bl_ed_pppm * sp_bas_temp$ppm_range)
  } else {
    lambda <- opts$lambda
  }
  
  #### fit resut ####
  # apply fit parameters to data and basis and construct fit result object
  
  final_par <- res$par
  
  # diagnostic info frame 
  diags <- data.frame(phase = res$par[1] * 180 / pi, lw = res$par[2],
                      shift = -res$par[3] / acq_paras$ft * 1e6,
                      asym = res$par[4], res$deviance, res$niter, res$info,
                      res$message, bl_ed_pppm = opts$bl_ed_pppm,
                      stringsAsFactors = TRUE)
  
  if (opts$auto_bl_flex) diags$max_bl_flex_used <- max_bl_flex_used
  
  if (opts$phi1_optim) diags$phi1 <- res$par[5]
  
  # apply phase parameter to data
  y_mod <- y * exp(1i * res$par[1])
  
  # apply shift parameter to data
  y_mod <- y_mod * exp(2i * pi * t * res$par[3])
  Y_mod <- ft_shift(y_mod)
  
  # apply phi1 phase parameter if adjustment has been requested
  if (opts$phi1_optim) Y_mod <- Y_mod * exp(-2i * pi * f * res$par[5] / 1000)
  
  # apply lineshape parameter to metab basis
  fs <- 1 / t[2]
  N  <- length(t)
  damp_vec <- asy_pvoigt_ls(fs, N, res$par[2], 1, res$par[4], TRUE)
  
  damp_mat <- matrix(damp_vec, nrow = nrow(raw_metab_basis), 
                     ncol = ncol(raw_metab_basis), byrow = F)
  metab_basis_damped <- raw_metab_basis * damp_mat
  
  # apply shift and lb terms to individual basis signals
  Nbasis <- dim(raw_metab_basis)[2]
  t_mat <- matrix(t, nrow = N, ncol = Nbasis)
  
  if (opts$phi1_optim) {
    freq_vec_hz <- res$par[6:(5 + Nbasis)]
    lb_vec_hz   <- res$par[(6 + Nbasis):(5 + 2 * Nbasis)]
  } else {
    freq_vec_hz <- res$par[5:(4 + Nbasis)]
    lb_vec_hz   <- res$par[(5 + Nbasis):(4 + 2 * Nbasis)]
  }
  
  freq_vec      <- 2i * pi * freq_vec_hz
  lb_vec        <- lw2alpha(lb_vec_hz)
  
  freq_lb_mat <- matrix(freq_vec - lb_vec, nrow = N, ncol = Nbasis,
                        byrow = TRUE) 
  
  metab_basis_damped <- metab_basis_damped * exp(t_mat * freq_lb_mat)
  
  # back to fd
  raw_metab_basis_fd <- ft_shift_mat(metab_basis_damped)
  
  # generate spline basis
  sp_bas_final <- generate_sp_basis(mrs_data, opts$ppm_right, opts$ppm_left,
                                    opts$bl_comps_pppm)
  
  # cut out fit region 
  metab_basis_fd_cut <- raw_metab_basis_fd[sp_bas_final$inds,,drop = FALSE]
  metab_comps <- dim(raw_metab_basis_fd)[2]
  
  # augment metabolite basis with zeros to match spline basis dimensions 
  metab_basis_fd_aug <- rbind(metab_basis_fd_cut,
                              matrix(0, sp_bas_final$bl_comps - 2, metab_comps))
  
  bl_basis_final <- rbind(sp_bas_final$bl_bas, sp_bas_final$deriv_mat *
                    (lambda ^ 0.5))
  
  # combine metabolite and spline basis 
  full_bas <- cbind(bl_basis_final, metab_basis_fd_aug)
  
  # augment spectrum with zeros to match full basis
  fit_seg <- c(Y_mod[sp_bas_final$inds], rep(0, sp_bas_final$bl_comps - 2))
  
  # estimate amplitudes
  ahat <- calc_ahat(Re(full_bas), Re(fit_seg), k = sp_bas_final$bl_comps,
                    opts$ahat_calc_method)
  ahat[is.na(ahat)] <- 0
  
  # signal estimate (fit)
  Y_hat <- Re(full_bas) %*% ahat
  
  # residual
  res <- Re(fit_seg) - Y_hat
  full_res <- sum(res ^ 2) # full residual inc. penalty
  diags$full_res <- full_res
  
  # if maxiters = 0 and maxiters_pre = 0 then we need to calculate the spectral
  # residual in the final phase
  if (is.na(diags$res.deviance)) {
    diags$res.deviance = sum(res[1:length(sp_bas_final$inds)] ^ 2)
  }
  
  # number of data points used in the fit
  diags$fit_pts <- length(sp_bas_final$inds)
  
  # ppm range in the fit
  diags$ppm_range <- sp_bas_final$ppm_range
  
  # baseline
  bl <- Re(sp_bas_final$bl_bas) %*% ahat[1:sp_bas_final$bl_comps]
  
  # turn amplitudes into a matrix to scale invividual signals
  metab_amp_mat <- matrix(ahat[(sp_bas_final$bl_comps + 1):length(ahat)], 
                          nrow = nrow(metab_basis_fd_cut),
                          ncol = ncol(metab_basis_fd_cut), byrow = TRUE)
  
  # scaled metabolite basis signals
  basis_frame <- as.data.frame(Re(metab_basis_fd_cut) * metab_amp_mat, 
                               row.names = NA)
  
  colnames(basis_frame) <- basis$names
  
  # turn amplitudes into a matrix to scale invividual splines
  spline_amp_mat <- matrix(ahat[(1:sp_bas_final$bl_comps)], 
                          nrow = nrow(sp_bas_final$bl_bas),
                          ncol = ncol(sp_bas_final$bl_bas), byrow = TRUE)
  
  # remove augmented data points for plotting
  bl_plot <- as.vector(bl[1:length(sp_bas_final$inds)])
  Y_plot <- as.vector(Re(Y_mod[sp_bas_final$inds]))
  Y_hat_plot <- as.vector(Y_hat[1:length(sp_bas_final$inds)])
  
  # fit frame for plotting
  fit_frame <- data.frame(PPMScale = sp_bas_final$x_scale, Data = Y_plot,
                          Fit = Y_hat_plot - bl_plot, Baseline = bl_plot)
  
  # SNR calc
  mrs_data_corr <- vec2mrs_data(y_mod, fs = acq_paras$fs, ft = acq_paras$ft, 
                                ref = acq_paras$ref)
  
  noise_sd  <- as.numeric(calc_spec_snr(mrs_data_corr,
                                        noise_region = opts$noise_region,
                                        full_output = TRUE)$noise_sd)
  
  res_sd <- sd(res[1:length(sp_bas_final$inds)])
  diags$SNR <- max(fit_frame$Fit) / noise_sd
  diags$SRR <- max(fit_frame$Fit) / res_sd
  diags$FQN <- stats::var(res[1:length(sp_bas_final$inds)]) / noise_sd ^ 2 
  
  # combine fit and basis signals 
  fit_frame <- cbind(fit_frame, basis_frame)
  class(fit_frame) <- c("fit_table", "data.frame")
  
  if (opts$export_sp_fit) {
    # scaled spline basis signals
    spline_frame <- as.data.frame(Re(sp_bas_final$bl_bas) * spline_amp_mat, 
                                  row.names = NA)
    
    colnames(spline_frame) <- paste("SP_",
                                    as.character(1:ncol(sp_bas_final$bl_bas)),
                                    sep = "")
    
    fit_frame <- cbind(fit_frame, spline_frame)
  }
  
  # metabolite amplitudes
  amps <- data.frame(t(ahat[(sp_bas_final$bl_comps + 1):length(ahat)]))
  colnames(amps) <- basis$names
  
  # tNAA_lw calc 
  if (("NAA" %in% colnames(amps)) & ("NAAG" %in% colnames(amps))) {
    tnaa_sig_pts <- basis_frame$NAA + basis_frame$NAAG
    diags$tNAA_lw <- calc_peak_info_vec(tnaa_sig_pts, 2)[3] * 
                     (sp_bas_final$x_scale[1] - sp_bas_final$x_scale[2])
  } else if ("NAA" %in% colnames(amps)) {
    diags$NAA_lw <- calc_peak_info_vec(basis_frame$NAA, 2)[3] * 
                    (sp_bas_final$x_scale[1] - sp_bas_final$x_scale[2])
  }
  
  #### crlb calc ####
  
  # calculate the analytical jacobian for the non-linear parameters
  para_crlb     <- abfit_full_anal_jac(final_par, y, raw_metab_basis,
                                       bl_basis_final, t, f, sp_bas_final$inds,
                                       sp_bas_final$bl_comps, FALSE, NULL,
                                       opts$phi1_optim, opts$ahat_calc_method,
                                       NULL) # nb reg_freq not included in the
                                             # crlb calc
   
  bl_comps_crlb <- sp_bas_final$bl_comps
  para_crlb     <- rbind(Re(para_crlb), matrix(0, nrow = bl_comps_crlb - 2,
                                               ncol = ncol(para_crlb)))
  amp_crlb      <- Re(full_bas)
  
  ## alternate method by generating spline basis with ED components in place of
  ## the penatly matrix. Gives similar numbers to above method, so kept here as
  ## a reference example for validation.
  ##
  ## ppm_range   <- opts$ppm_left - opts$ppm_right
  ## crlb_comps  <- round(opts$bl_ed_pppm * ppm_range) - 2
  ## sp_bas_crlb <- generate_sp_basis(mrs_data, opts$ppm_right, opts$ppm_left,
  ##                                  crlb_comps / ppm_range)
  ## 
  ## bl_comps_crlb <- sp_bas_crlb$bl_comps
  ## amp_crlb      <- cbind(sp_bas_crlb$bl_bas, Re(metab_basis_fd_cut))
  
  # remove any zero columns from para_crlb
  para_crlb <- para_crlb[,!(colSums(abs(para_crlb)) == 0)] 
  
  # non-zero paras in jacobian
  nz_paras <- dim(para_crlb)[2]
  
  D <- cbind(para_crlb, amp_crlb)
  F <- t(D) %*% D / (res_sd ^ 2)
  #F_inv <- pracma::inv(F) # sometimes F becomes singular and inv fails
  F_inv <- pracma::pinv(F)
  crlbs <- diag(F_inv) ^ 0.5
  crlbs_out <- crlbs[(bl_comps_crlb + nz_paras + 1):length(crlbs)]
 
  comb_sigs <- 0
  
  D_cut <- D
  colnames(D_cut) <- c(rep("para", nz_paras), rep("sp", bl_comps_crlb),
                       colnames(amps))
  
  D_cut <- as.data.frame(D_cut)
  
  # create some common metabolite combinations
  if (("NAA" %in% colnames(amps)) & ("NAAG" %in% colnames(amps))) {
    amps['tNAA'] <- amps['NAA'] + amps['NAAG']
    comb_sigs <- comb_sigs + 1
    D_cut$tNAA <- D_cut$NAA + D_cut$NAAG
    D_cut$NAA <- NULL
    D_cut$NAAG <- NULL
  }
  
  if (("Cr" %in% colnames(amps)) & ("PCr" %in% colnames(amps))) {
    amps['tCr'] <- amps['Cr'] + amps['PCr']
    comb_sigs <- comb_sigs + 1
    D_cut$tCr <- D_cut$Cr + D_cut$PCr
    D_cut$Cr <- NULL
    D_cut$PCr <- NULL
  }
  
  if (("PCh" %in% colnames(amps)) & ("GPC" %in% colnames(amps))) {
    amps['tCho'] <- amps['PCh'] + amps['GPC']
    comb_sigs <- comb_sigs + 1
    D_cut$tCho <- D_cut$PCh + D_cut$GPC
    D_cut$PCh <- NULL
    D_cut$GPC <- NULL
  }
  
  if (("Glu" %in% colnames(amps)) & ("Gln" %in% colnames(amps))) {
    amps['Glx'] <- amps['Glu'] + amps['Gln']
    comb_sigs <- comb_sigs + 1
    D_cut$Glx <- D_cut$Glu + D_cut$Gln
    D_cut$Glu <- NULL
    D_cut$Gln <- NULL
  }
  
  if (("Lip09" %in% colnames(amps)) & ("MM09" %in% colnames(amps))) {
    amps['tLM09'] <- amps['Lip09'] + amps['MM09']
    comb_sigs <- comb_sigs + 1
    D_cut$tLM09 <- D_cut$Lip09 + D_cut$MM09
    D_cut$Lip09 <- NULL
    D_cut$MM09 <- NULL
  }
  
  if (("Lip13a" %in% colnames(amps)) & ("Lip13b" %in% colnames(amps)) & 
        ("MM12" %in% colnames(amps)) & ("MM14" %in% colnames(amps))) {
    amps["tLM13"] <- amps["Lip13a"] + amps["Lip13b"] + amps["MM12"] + 
                     amps["MM14"]
    comb_sigs <- comb_sigs + 1
    D_cut$tLM13 <- D_cut$Lip13a + D_cut$Lip13b + D_cut$MM12 + D_cut$MM14
    D_cut$Lip13a <- NULL
    D_cut$Lip13b <- NULL
    D_cut$MM12 <- NULL
    D_cut$MM14 <- NULL
  }
  
  if (("Lip20" %in% colnames(amps)) & ("MM20" %in% colnames(amps))) {
    amps['tLM20'] <- amps['Lip20'] + amps['MM20']
    comb_sigs <- comb_sigs + 1
    D_cut$tLM20 <- D_cut$Lip20 + D_cut$MM20
    D_cut$Lip20 <- NULL
    D_cut$MM20  <- NULL
  }
  
  if (comb_sigs != 0) {
    D_cut <- as.matrix(D_cut)
    F_cut <- t(D_cut) %*% D_cut / (res_sd ^ 2)
    # F_cut_inv <- inv(F_cut) # sometimes F becomes singular and inv fails
    F_cut_inv <- pracma::pinv(F_cut)
    crlbs_cut <- diag(F_cut_inv) ^ 0.5
    crlb_n <- length(crlbs_cut)
    # append combined CRLB estimates
    crlbs_out <- as.numeric(c(crlbs_out,
                              crlbs_cut[(crlb_n - comb_sigs + 1):crlb_n]))
  }
  
  # vector of coarse fit residuals at differing levels of baseline flexibility
  if (opts$auto_bl_flex) diags <- cbind(diags, resid_vec_frame)
  
  if (opts$output_all_paras) {
    # vector of shifts
    freq_vec_ppm <- -freq_vec_hz / acq_paras$ft * 1e6
    names(freq_vec_ppm) <- paste(basis$names, "shift.ppm", sep = ".")
    names(lb_vec_hz)    <- paste(basis$names, "lb.hz", sep = ".")
    diags <- cbind(diags, t(freq_vec_ppm))
    diags <- cbind(diags, t(lb_vec_hz))
  }
  
  # construct output
  list(amps = amps, crlbs = t(crlbs_out), diags = diags, fit = fit_frame)
}

#' Return a list of options for an ABfit analysis.
#' 
#' @param init_damping initial value of the Gaussian global damping parameter
#' (Hz). Very poorly shimmed or high field data may benefit from a larger value.
#' @param maxiters The maximum number of iterations to run for the detailed fit.
#' @param max_shift The maximum allowable shift to be applied in the
#' optimisation phase of fitting (ppm).
#' @param max_damping maximum permitted value of the global damping parameter
#' (Hz).
#' @param max_phase the maximum absolute permitted value of the global
#' zero-order phase term (degrees). Note, the prefit_phase_search option is not
#' constrained by this term.
#' @param lambda manually set the the baseline smoothness parameter.
#' @param ppm_left downfield frequency limit for the fitting range (ppm).
#' @param ppm_right upfield frequency limit for the fitting range (ppm).
#' @param zp zero pad the data to twice the original length before fitting.
#' @param bl_ed_pppm manually set the the baseline smoothness parameter (ED per
#' ppm).
#' @param auto_bl_flex automatically determine the level of baseline smoothness.
#' @param bl_comps_pppm spline basis density (signals per ppm).
#' @param export_sp_fit add the fitted spline functions to the fit result.
#' @param max_asym maximum allowable value of the asymmetry parameter.
#' @param max_basis_shift maximum allowable frequency shift for individual basis
#' signals (ppm).
#' @param max_basis_damping maximum allowable Lorentzian damping factor for
#' individual basis signals (Hz).
#' @param maxiters_pre maximum iterations for the coarse (pre-)fit.
#' @param algo_pre optimisation method for the coarse (pre-)fit.
#' @param min_bl_ed_pppm minimum value for the candidate baseline flexibility
#' analyses (ED per ppm).
#' @param max_bl_ed_pppm minimum value for the candidate baseline flexibility
#' analyses (ED per ppm).
#' @param auto_bl_flex_n number of candidate baseline analyses to perform.
#' @param pre_fit_bl_ed_pppm level of baseline flexibility to use in the coarse
#' fitting stage of the algorithm (ED per ppm).
#' @param remove_lip_mm_prefit remove broad signals in the coarse fitting stage
#' of the algorithm.
#' @param pre_align perform a pre-alignment step before coarse fitting.
#' @param max_pre_align_shift maximum allowable shift in the pre-alignment step
#' (ppm).
#' @param pre_align_ref_freqs a vector of prominent spectral frequencies used in
#' the pre-alignment step (ppm).
#' @param noise_region spectral region to estimate the noise level (ppm).
#' @param optimal_smooth_criterion method to determine the optimal smoothness.
#' @param aic_smoothing_factor modification factor for the AIC calculation.
#' @param anal_jac use a analytical approximation to the jacobian in the 
#' detailed fitting stage.
#' @param pre_fit_ppm_left downfield frequency limit for the fitting range in
#' the coarse fitting stage of the algorithm (ppm).
#' @param pre_fit_ppm_right upfield frequency limit for the fitting range in the
#' coarse fitting stage of the algorithm (ppm).
#' @param phi1_optim apply and optimise a frequency dependant phase term.
#' @param phi1_init initial value for the frequency dependant phase term (ms).
#' @param max_dphi1 maximum allowable change from the initial frequency
#' dependant phase term (ms).
#' @param max_basis_shift_broad maximum allowable shift for broad signals in the
#' basis (ppm). Determined based on their name beginning with Lip or MM.
#' @param max_basis_damping_broad maximum allowable Lorentzian damping for broad
#' signals in the basis (Hz). Determined based on their name beginning with Lip
#' or MM.
#' @param ahat_calc_method method to calculate the metabolite amplitudes. May be
#' one of: "lh_pnnls" or "ls".
#' @param prefit_phase_search perform a 1D search for the optimal phase in the
#' prefit stage of the algorithm.
#' @param freq_reg frequency shift parameter.
#' @param output_all_paras include more fitting parameters in the fit table,
#' e.g. individual shift and damping factors for each basis set element.
#' @return full list of options.
#' @examples
#' opts <- abfit_opts(ppm_left = 4.2, noise_region = c(-1, -3))
#' @export
abfit_opts <- function(init_damping = 5, maxiters = 1024, max_shift = 0.078, 
                       max_damping = 15, max_phase = 360, lambda = NULL, 
                       ppm_left = 4, ppm_right = 0.2, zp = TRUE,
                       bl_ed_pppm = 2.0, auto_bl_flex = TRUE,
                       bl_comps_pppm = 15, export_sp_fit = FALSE,
                       max_asym = 0.25, max_basis_shift = 0.0078,
                       max_basis_damping = 2, maxiters_pre = 1000,
                       algo_pre = "NLOPT_LN_NELDERMEAD", min_bl_ed_pppm = NULL,
                       max_bl_ed_pppm = 7, auto_bl_flex_n = 20, 
                       pre_fit_bl_ed_pppm = 1, remove_lip_mm_prefit = FALSE,
                       pre_align = TRUE, max_pre_align_shift = 0.1,
                       pre_align_ref_freqs = c(2.01, 3.03, 3.22),
                       noise_region = c(-0.5, -2.5),
                       optimal_smooth_criterion = "maic",
                       aic_smoothing_factor = 5, anal_jac = TRUE,
                       pre_fit_ppm_left = 4, pre_fit_ppm_right = 1.8,
                       phi1_optim = FALSE, phi1_init = 0, max_dphi1 = 0.2,
                       max_basis_shift_broad = 0.0078,
                       max_basis_damping_broad = 2,
                       ahat_calc_method = "lh_pnnls",
                       prefit_phase_search = TRUE, freq_reg = NULL,
                       output_all_paras = FALSE) {
                         
  list(init_damping = init_damping, maxiters = maxiters,
       max_shift = max_shift, max_damping = max_damping, max_phase = max_phase,
       lambda = lambda, ppm_left = ppm_left, ppm_right = ppm_right, zp = zp,
       bl_ed_pppm = bl_ed_pppm, auto_bl_flex = auto_bl_flex,
       bl_comps_pppm = bl_comps_pppm, export_sp_fit = export_sp_fit,
       max_asym = max_asym, max_basis_shift = max_basis_shift, 
       max_basis_damping = max_basis_damping, maxiters_pre = maxiters_pre,
       algo_pre = algo_pre, min_bl_ed_pppm = min_bl_ed_pppm,
       max_bl_ed_pppm = max_bl_ed_pppm, auto_bl_flex_n = auto_bl_flex_n,
       pre_fit_bl_ed_pppm = pre_fit_bl_ed_pppm,
       remove_lip_mm_prefit = remove_lip_mm_prefit, pre_align = pre_align,
       max_pre_align_shift = max_pre_align_shift,
       pre_align_ref_freqs = pre_align_ref_freqs, noise_region = noise_region,
       optimal_smooth_criterion = optimal_smooth_criterion,
       aic_smoothing_factor = aic_smoothing_factor, anal_jac = anal_jac,
       pre_fit_ppm_left = pre_fit_ppm_left,
       pre_fit_ppm_right = pre_fit_ppm_right, phi1_optim = phi1_optim,
       phi1_init = phi1_init, max_dphi1 = max_dphi1,
       max_basis_shift_broad = max_basis_shift_broad,
       max_basis_damping_broad = max_basis_damping_broad,
       ahat_calc_method = ahat_calc_method,
       prefit_phase_search = prefit_phase_search, freq_reg = freq_reg,
       output_all_paras = output_all_paras)
}

#' Return a list of options for an ABfit analysis to maintain comparability with
#' analyses performed with version 1.9.0 (and earlier) of spant.
#' @param ... arguments passed to [spant::abfit_opts].
#' @return full list of options.
#' @export
abfit_opts_v1_9_0 <- function(...) {
  opts <- abfit_opts(prefit_phase_search = FALSE, ...)
  return(opts)
}

# objective function for 4 parameter full spine fitting method
abfit_full_obj <- function(par, y, raw_metab_basis, bl_basis, t, f, inds,
                           bl_comps, sum_sq, basis_paras, phi1_optim,
                           ahat_calc_method, freq_reg) {
  
  if (!is.null(basis_paras)) par <- c(par, basis_paras)
  
  # apply phase parameter to data
  y <- y * exp(1i * par[1])
  
  # apply shift parameter to data
  y <- y * exp(2i * pi * t * par[3])
  Y <- ft_shift(y)
  
  # apply phi1 phase parameter if adjustment has been requested
  if (phi1_optim) Y <- Y * exp(-2i * pi * f * par[5] / 1000)
  
  # apply lineshape parameters to basis
  fs <- 1 / t[2]
  N <- length(t)
  damp_vec <- asy_pvoigt_ls(fs, N, par[2], 1, par[4], TRUE)
  
  damp_mat <- matrix(damp_vec, nrow = nrow(raw_metab_basis), 
                     ncol = ncol(raw_metab_basis), byrow = F)
  raw_metab_basis <- raw_metab_basis * damp_mat
  
  # apply shift and lb terms to individual basis signals
  Nbasis <- dim(raw_metab_basis)[2]
  t_mat <- matrix(t, nrow = N, ncol = Nbasis)
  
  if (phi1_optim) {
    freq_shifts <- par[6:(5 + Nbasis)]
    freq_vec <- 2i * pi * freq_shifts
    lb_vec <- lw2alpha(par[(6 + Nbasis):(5 + 2 * Nbasis)])
  } else {
    freq_shifts <- par[5:(4 + Nbasis)]
    freq_vec <- 2i * pi * freq_shifts
    lb_vec <- lw2alpha(par[(5 + Nbasis):(4 + 2 * Nbasis)])
  }
  
  freq_lb_mat <- matrix(freq_vec - lb_vec, nrow = N, ncol = Nbasis,
                        byrow = TRUE) 
  
  raw_metab_basis <- raw_metab_basis * exp(t_mat * freq_lb_mat)
  
  # back to fd
  raw_metab_basis <- ft_shift_mat(raw_metab_basis)
  
  # cut out fitting region
  raw_metab_basis <- raw_metab_basis[inds,,drop = FALSE]
  
  # augment metabolite basis with zeros to match spline basis 
  metab_comps <- dim(raw_metab_basis)[2]
  raw_metab_basis <- rbind(raw_metab_basis, 
                           matrix(0, bl_comps - 2, metab_comps))
  
  # combine metabolite and spline basis
  full_bas <- cbind(bl_basis, raw_metab_basis)
  
  # augment signal with zeros to match basis dimensions
  fit_seg <- c(Y[inds], rep(0, bl_comps - 2))
  
  # estimtate amplitudes
  ahat <- calc_ahat(Re(full_bas), Re(fit_seg), k = bl_comps,
                    ahat_calc_method)
  
  # signal estimate 
  Y_hat <- Re(full_bas) %*% ahat
  
  if (sum_sq) {
    return(sum((Re(Y[inds]) - Y_hat[1:length(inds)]) ^ 2))
  } else {
    if (is.null(freq_reg)) {
      return(Re(Y[inds]) - Y_hat[1:length(inds)])
    } else {
      return(c(Re(Y[inds]) - Y_hat[1:length(inds)], freq_reg * freq_shifts))
    }
  }
}

abfit_full_num_jac <- function(par, y, raw_metab_basis, bl_basis, t, f,
                               inds, bl_comps, sum_sq, basis_paras,
                               phi1_optim, ahat_calc_method, freq_reg) {
  
  numDeriv::jacobian(func = abfit_full_obj, x = par, method = "simple",
                     y = y, raw_metab_basis = raw_metab_basis,
                     bl_basis = bl_basis, t = t, f = f, inds = inds,
                     bl_comps = bl_comps, sum_sq = sum_sq, basis_paras = NULL,
                     phi1_optim = phi1_optim,
                     ahat_calc_method = ahat_calc_method, freq_reg = freq_reg)
}

abfit_partial_num_jac <- function(par, y, raw_metab_basis, bl_basis, t, f,
                                  inds, bl_comps, sum_sq, basis_paras,
                                  phi1_optim, ahat_calc_method, freq_reg) {
  
  if (phi1_optim) {
    global_paras <- par[1:5]
    basis_paras_fixed <- par[6:length(par)]
  } else {
    global_paras <- par[1:4]
    basis_paras_fixed <- par[5:length(par)]
  }
  
  numDeriv::jacobian(func = abfit_full_obj, x = global_paras, method = "simple",
                     y = y, raw_metab_basis = raw_metab_basis,
                     bl_basis = bl_basis, t = t, f = f, inds = inds,
                     bl_comps = bl_comps, sum_sq = sum_sq,
                     basis_paras = basis_paras_fixed, phi1_optim = phi1_optim,
                     ahat_calc_method = ahat_calc_method, freq_reg = freq_reg)
}

# attempt to calc approx Jacobian for some of the global parameters
# seems to be less accuarate
abfit_full_anal_jac_test <- function(par, y, raw_metab_basis, bl_basis, t,
                                     f, inds, bl_comps, sum_sq,
                                     basis_paras, phi1_optim,
                                     ahat_calc_method, freq_reg) {
  
  # calculate the first 4 paras numerically
  global_paras_jac <- abfit_partial_num_jac(par, y, raw_metab_basis, bl_basis,
                                            t, f, inds, bl_comps, sum_sq,
                                            basis_paras, phi1_optim,
                                            ahat_calc_method, freq_reg)
  
  # apply phase parameter to data
  y <- y * exp(1i * par[1])
  
  # apply shift parameter to data
  y <- y * exp(2i * pi * t * par[3])
  Y <- ft_shift(y)
  
  # apply lineshape parameters to basis
  fs <- 1 / t[2]
  N <- length(t)
  damp_vec <- asy_pvoigt_ls(fs, N, par[2], 1, par[4], TRUE)
  
  damp_mat <- matrix(damp_vec, nrow = nrow(raw_metab_basis), 
                     ncol = ncol(raw_metab_basis), byrow = F)
  raw_metab_basis <- raw_metab_basis * damp_mat
  raw_metab_basis_orig <- raw_metab_basis
  
  # apply shift and lb terms to individual basis signals
  Nbasis <- dim(raw_metab_basis)[2]
  t_mat <- matrix(t, nrow = N, ncol = Nbasis)
  freq_vec <- 2i * pi * par[5:(4 + Nbasis)]
  lb_vec <- lw2alpha(par[(5 + Nbasis):(4 + 2 * Nbasis)])
  freq_lb_mat <- matrix(freq_vec - lb_vec, nrow = N, ncol = Nbasis,
                        byrow = TRUE) 
  
  raw_metab_basis <- raw_metab_basis * exp(t_mat * freq_lb_mat)
  
  # back to fd
  raw_metab_basis <- ft_shift_mat(raw_metab_basis)
  
  # cut out fitting region
  raw_metab_basis <- raw_metab_basis[inds,,drop = FALSE]
  
  # augment metabolite basis with zeros to match spline basis 
  metab_comps <- dim(raw_metab_basis)[2]
  raw_metab_basis <- rbind(raw_metab_basis, 
                           matrix(0, bl_comps - 2, metab_comps))
  
  # combine metabolite and spline basis
  full_bas <- cbind(bl_basis, raw_metab_basis)
  
  # augment signal with zeros to match basis dimensions
  fit_seg <- c(Y[inds], rep(0, bl_comps - 2))
  
  # estimtate amplitudes
  ahat <- calc_ahat(Re(full_bas), Re(fit_seg), k = bl_comps, ahat_calc_method)
  
  # multiply by a * 2i * pi * f * t for basis freq shifts
  freq_vec <- ahat[(bl_comps + 1):length(ahat)] * 2i * pi
  freq_mat <- matrix(freq_vec, nrow = N, ncol = Nbasis, byrow = TRUE)
  freq_jac <- raw_metab_basis_orig * t_mat * freq_mat
  
  # back to fd
  freq_jac <- ft_shift_mat(freq_jac)
  
  # cut out fitting region
  freq_jac <- -Re(freq_jac[inds,,drop = FALSE])
  
  # multiply by -a * t for basis lw adjustment
  lb_vec <- -ahat[(bl_comps + 1):length(ahat)]
  lb_mat <- matrix(lb_vec, nrow = N, ncol = Nbasis, byrow = TRUE)
  lb_jac <- raw_metab_basis_orig * t_mat * lb_mat * pi
  
  # back to fd
  lb_jac <- ft_shift_mat(lb_jac)
  
  # cut out fitting region
  lb_jac <- -Re(lb_jac[inds,,drop = FALSE])
  
  # global phase term
  global_paras_jac[, 1] <- Re(1i * Y[inds])
  
  # global freq term
  global_paras_jac[, 3] <- Re(ft_shift(2i * pi * t * y)[inds])
  
  return(cbind(global_paras_jac, freq_jac, lb_jac))
}

# jacobian function for the full spine fitting method
abfit_full_anal_jac <- function(par, y, raw_metab_basis, bl_basis, t, f, inds,
                                bl_comps, sum_sq, basis_paras, phi1_optim,
                                ahat_calc_method, freq_reg) {
  
  # calculate the first 4 paras numerically
  global_paras_jac <- abfit_partial_num_jac(par, y, raw_metab_basis, bl_basis,
                                            t, f, inds, bl_comps, sum_sq,
                                            basis_paras, phi1_optim,
                                            ahat_calc_method, freq_reg)
  
  # apply phase parameter to data
  y <- y * exp(1i * par[1])
  
  # apply shift parameter to data
  y <- y * exp(2i * pi * t * par[3])
  Y <- ft_shift(y)
  
  # apply phi1 phase parameter if adjustment has been requested
  if (phi1_optim) Y <- Y * exp(-2i * pi * f * par[5] / 1000)
  
  # apply lineshape parameters to basis
  fs <- 1 / t[2]
  N <- length(t)
  damp_vec <- asy_pvoigt_ls(fs, N, par[2], 1, par[4], TRUE)
  
  damp_mat <- matrix(damp_vec, nrow = nrow(raw_metab_basis), 
                     ncol = ncol(raw_metab_basis), byrow = F)
  raw_metab_basis <- raw_metab_basis * damp_mat
  raw_metab_basis_orig <- raw_metab_basis
  
  # apply shift and lb terms to individual basis signals
  Nbasis <- dim(raw_metab_basis)[2]
  t_mat <- matrix(t, nrow = N, ncol = Nbasis)
  
  if (phi1_optim) {
    freq_vec <- 2i * pi * par[6:(5 + Nbasis)]
    lb_vec <- lw2alpha(par[(6 + Nbasis):(5 + 2 * Nbasis)])
  } else {
    freq_vec <- 2i * pi * par[5:(4 + Nbasis)]
    lb_vec <- lw2alpha(par[(5 + Nbasis):(4 + 2 * Nbasis)])
  }
  
  freq_lb_mat <- matrix(freq_vec - lb_vec, nrow = N, ncol = Nbasis,
                        byrow = TRUE) 
  
  raw_metab_basis <- raw_metab_basis * exp(t_mat * freq_lb_mat)
  raw_metab_basis_mod <- raw_metab_basis
  
  # back to fd
  raw_metab_basis <- ft_shift_mat(raw_metab_basis)
  
  # cut out fitting region
  raw_metab_basis <- raw_metab_basis[inds,,drop = FALSE]
  
  # augment metabolite basis with zeros to match spline basis 
  metab_comps <- dim(raw_metab_basis)[2]
  raw_metab_basis <- rbind(raw_metab_basis, 
                           matrix(0, bl_comps - 2, metab_comps))
  
  # combine metabolite and spline basis
  full_bas <- cbind(bl_basis, raw_metab_basis)
  
  # augment signal with zeros to match basis dimensions
  fit_seg <- c(Y[inds], rep(0, bl_comps - 2))
  
  # estimtate amplitudes
  ahat <- calc_ahat(Re(full_bas), Re(fit_seg), k = bl_comps, ahat_calc_method)
  
  # multiply by a * 2i * pi * t for basis freq shifts
  freq_vec <- ahat[(bl_comps + 1):length(ahat)] * 2i * pi
  freq_mat <- matrix(freq_vec, nrow = N, ncol = Nbasis, byrow = TRUE)
  freq_jac <- raw_metab_basis_mod * t_mat * freq_mat
  
  # back to fd
  freq_jac <- ft_shift_mat(freq_jac)
  
  # cut out fitting region
  freq_jac <- -Re(freq_jac[inds,,drop = FALSE])
  
  # multiply by -a * t * pi for basis lw adjustment
  lb_vec <- ahat[(bl_comps + 1):length(ahat)]
  lb_mat <- matrix(lb_vec, nrow = N, ncol = Nbasis, byrow = TRUE)
  lb_jac <- -raw_metab_basis_mod * t_mat * lb_mat * pi
  
  # back to fd
  lb_jac <- ft_shift_mat(lb_jac)
  
  # cut out fitting region
  lb_jac <- -Re(lb_jac[inds,,drop = FALSE])
  
  if (is.null(freq_reg)) {
    ret_mat <- cbind(global_paras_jac, freq_jac, lb_jac)
  } else {
    ret_mat <- cbind(global_paras_jac,
                     rbind(freq_jac, matrix(0, ncol = Nbasis, nrow = Nbasis)),
                     rbind(lb_jac,   matrix(0, ncol = Nbasis, nrow = Nbasis)))
  }
  
  return(ret_mat)
}

# objective function for 3 parameter spine fitting method
abfit_3p_obj <- function(par, y, raw_metab_basis, bl_basis, t, f, inds,
                           bl_comps, sum_sq, ret_full, ahat_calc_method) {
  
  # apply phase parameter to data
  y <- y * exp(1i * par[1])
  
  # apply shift parameter to data
  y <- y * exp(2i * pi * t * par[3])
  Y <- ft_shift(y)
  
  # apply damping parameter to basis
  damp_vec <- exp(-t * t * lw2beta(par[2]))
  damp_mat <- matrix(damp_vec, nrow = nrow(raw_metab_basis), 
                     ncol = ncol(raw_metab_basis), byrow = F)
  raw_metab_basis <- raw_metab_basis * damp_mat
  
  # back to fd
  raw_metab_basis <- ft_shift_mat(raw_metab_basis)
  
  # cut out fitting region
  raw_metab_basis <- raw_metab_basis[inds,,drop = FALSE]
  
  # augment metabolite basis with zeros to match spline basis 
  metab_comps <- dim(raw_metab_basis)[2]
  raw_metab_basis <- rbind(raw_metab_basis, 
                           matrix(0, bl_comps - 2, metab_comps))
  
  # combine metabolite and spline basis
  full_bas <- cbind(bl_basis, raw_metab_basis)
  
  # augment signal with zeros to match basis dimensions
  fit_seg <- c(Y[inds], rep(0, bl_comps - 2))
  
  # estimtate amplitudes
  ahat <- calc_ahat(Re(full_bas), Re(fit_seg), k = bl_comps, ahat_calc_method)
  
  # signal estimate 
  Y_hat <- Re(full_bas) %*% ahat
  
  if (ret_full) {
    res <- list(resid = (Re(Y[inds]) - Y_hat[1:length(inds)]), amp = ahat)
    return(res)
  }
  
  if (sum_sq) {
    return(sum((Re(Y[inds]) - Y_hat[1:length(inds)]) ^ 2))
  } else {
    return(Re(Y[inds]) - Y_hat[1:length(inds)])
  }
}

# Directly from "Splines, knots, and penalties", Eilers 2010
tpower <- function(x, t, p) {
  # Truncated p-th power function
  (x - t) ^ p * (x > t)
}

#' Generate a spline basis, slightly adapted from : "Splines, knots, and
#' penalties", Eilers 2010.
#'
#' @param N number of data points.
#' @param number number of spline functions.
#' @param deg spline degree : deg = 1 linear, deg = 2 quadratic, deg = 3 cubic.
#' @return spline basis as a matrix.
#' @export
bbase <- function(N, number, deg = 3) {
  # Construct a B-spline basis of degree
  x <- 1:N
  xl <- min(x)
  xr <- max(x)
  dx <- (xr - xl) / number
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
  B <- (-1) ^ (deg + 1) * P %*% t(D)
  B
}

# Gaussian lineshape
G <- function(freq, fwhm) {
  sqrt(4 * log(2) / pi / fwhm ^ 2) * exp(-4 * log(2) * (freq / fwhm) ^ 2)
}

# Lorentzian lineshape
L <- function(freq, fwhm) {
  fwhm / (2 * pi) / ((fwhm / 2) ^ 2 + freq ^ 2)
}

# Asymmetric lineshape function
asy_fwhm_fn <- function(freq, fwhm, asy) {
  asy_fwhm <- 2 * fwhm / (1 + exp(-asy * freq))
  asy_fwhm[asy_fwhm < .Machine$double.eps] <- .Machine$double.eps
  asy_fwhm
}

# Asymmetric pseudo-Voigt lineshape function
asy_pvoigt_ls <- function(fs, N, fwhm, lg = 0, asy = 0, td_only = FALSE) {
  freq_oversamp <- seq(from = -fs / 2, to = fs / 2 - fs / (N * 2),
                       length.out = N * 2)
  
  asy_fwhm <- asy_fwhm_fn(freq_oversamp, fwhm, asy)
  fd_oversamp <- lg * G(freq_oversamp, asy_fwhm) +
                 (1 - lg) * L(freq_oversamp, asy_fwhm)
  
  td <- pracma::ifft(pracma::ifftshift(fd_oversamp))[1:N]
  td <- td / Re(td[1])
  
  # for small fwhm values the td curve can remain constant when using the ifft
  # method above - this is a fix
  if (identical(unique(Re(td)), 1) && (fwhm != 0)) {
    t <- seq(from = 0, to = (N - 1) / fs, by = 1 / fs)
    beta <- lw2beta(fwhm)
    td <- exp(-(t ^ 2) * beta) + 0i
  }
    
  if (td_only) {
    return(td)
  } else {
    t <- seq(from = 0, to = (N - 1) / fs, by = 1 / fs)
    freq <- seq(from = -fs / 2, to = fs / 2 - fs / N, length.out = N)
    asy_fwhm <- asy_fwhm_fn(freq, fwhm, asy)
    fd <- lg * G(freq, asy_fwhm) + (1 - lg) * L(freq, asy_fwhm)
    return(list(fd = fd, td = td, freq = freq, t = t))
  }
}

auto_pspline_smoother <- function(y, spline_basis, deriv_mat, maxiters = 25,
                                  verbose = FALSE) {
  
  # set so that the initial lambda = 1
  sig_sq <- 1
  tau_sq <- 1
  
  # ED of 1 is not possible - so at least one iteration is required
  ED <- 1
  
  for (iter in 1:maxiters) {
    lambda <- sig_sq / tau_sq
    
    inv_mat <- MASS::ginv(rbind(spline_basis, lambda ^ 0.5 * deriv_mat))
    alpha <-   inv_mat %*% c(y, rep(0, dim(deriv_mat)[1]))
    
    yhat <- spline_basis %*% alpha
    
    zero_mat <- matrix(0, dim(deriv_mat)[1], dim(deriv_mat)[2])
    G        <- inv_mat %*% rbind(spline_basis, zero_mat)
    ED_new   <- sum(diag(G))
    
    sig_sq <- sum((y - yhat) ^ 2) / (length(y) - ED_new)
    
    # - 2 below comes from the d = difference pen matrix
    tau_sq <- sum(((deriv_mat) %*% alpha) ^ 2) / (ED_new - 2) 
    lambda_new <- sig_sq / (tau_sq + sig_sq * 1e-8)
   
    if (verbose) cat("iter :", iter, "lambda =", lambda_new,
                     "ED =", ED_new, "\n")
    
    if (ED_new < 2) {
      warning("instability in auto pspline smoother")
      ED_new <- 2
    }
    
    # check for convergence
    if (Mod(ED - ED_new) < 0.01 | ED_new < 2.01) {
      return(list(yhat = yhat, ED = ED_new, lambda = lambda, iters = iter,
                  converged = TRUE))
    }
     
    ED <- ED_new
  }
  return(list(yhat = yhat, ED = ED, lambda = lambda, iters = iter,
              converged = FALSE))
}

calc_ed_from_lambda_stable <- function(spline_basis, deriv_mat, lambda) {
  inv_mat  <- MASS::ginv(rbind(spline_basis, (lambda ^ 0.5) * deriv_mat))
  zero_mat <- matrix(0, dim(deriv_mat)[1], dim(deriv_mat)[2])
  G        <- inv_mat %*% rbind(spline_basis, zero_mat)
  return(sum(diag(G)))
}

#' Calculate the effective dimensions of a spline smoother from lambda.
#' @param spline_basis spline basis.
#' @param deriv_mat derivative matrix.
#' @param lambda smoothing parameter.
#' @return the effective dimension value.
#' @export
calc_ed_from_lambda <- function(spline_basis, deriv_mat, lambda) {
  inv_mat <- solve(t(spline_basis) %*% spline_basis +
                     lambda * (t(deriv_mat) %*% deriv_mat))
  H       <- inv_mat %*% (t(spline_basis) %*% spline_basis)
  return(sum(diag(H)))
}

ed_obj_fn <- function(par, spline_basis, deriv_mat, target_ed) {
  ed <- calc_ed_from_lambda(spline_basis, deriv_mat, par[1]) 
  return((ed - target_ed) ^ 2)
}

calc_lambda_from_ed <- function(spline_basis, deriv_mat, target_ed,
                                upper_lim = 1e10, lower_lim = 1e-6,
                                start_val = 1.0) {
  
  res <- stats::optim(start_val, ed_obj_fn, method = "Brent", lower = lower_lim,
               upper = upper_lim, spline_basis = spline_basis,
               deriv_mat = deriv_mat, target_ed = target_ed)
  
  if (res$value > 1e-6) warning("correct lambda not found")
  
  return(res$par)
}

generate_sp_basis <- function(mrs_data, ppm_right, ppm_left, bl_comps_pppm) {
  inds      <- get_seg_ind(ppm(mrs_data), ppm_right, ppm_left)
  x_scale   <- ppm(mrs_data)[inds]
  ppm_range <- ppm_left - ppm_right
  comps     <- round(bl_comps_pppm * ppm_range)
  bl_bas    <- bbase(length(inds), comps)
  bl_comps  <- comps + 3
  deriv_mat <- diff(diag(bl_comps), lag = 1, differences = 2)
  list(inds = inds, x_scale = x_scale, bl_bas = bl_bas, bl_comps = bl_comps,
       deriv_mat = deriv_mat, ppm_range = ppm_range)
}

# a - basis matrix (eg simulated metabolite signals)
# b - response vector (eg acquired data)
# k - the first k coefficients are not NN-restricted (eg number of spline fn's)
# method - one of: "lh_pnnls", "glmnet_pnnls", "ls"
calc_ahat <- function(a, b, k, ahat_calc_method) {
  if (ahat_calc_method == "lh_pnnls") {
    ahat <- pnnls(a, b, k = k)$x
  } else if (ahat_calc_method == "ls") {
    ahat <- stats::.lm.fit(a, b)$coefficients
  # } else if (ahat_calc_method == "bvls") {
  #   lower  <- rep(0, ncol(a))
  #   if (k > 0) lower[1:k] <- -Inf
  #   ahat <- bvls::bvls(a, b, bl = lower, bu = rep(Inf, ncol(a)))$x
  # } else if (ahat_calc_method == "glmnet_pnnls") {
  #   lower  <- rep(0, ncol(a))
  #   if (k > 0) lower[1:k] <- -Inf
  #   fit <- glmnet::glmnet(a, b, lambda = 0, lower.limits = lower,
  #                         intercept = FALSE)
  #   ahat <- glmnet::coef.glmnet(fit)[-1] # -1 to remove the intercept (always 0)
  } else {
    stop("Invalid method for calc_ahat.")
  }
  return(ahat)
}