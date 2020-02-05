varpro_3_para <- function(y, acq_paras, basis, opts = NULL) {
  mrs_data <- vec2mrs_data(y, fs = acq_paras$fs, ft = acq_paras$ft, 
                           ref = acq_paras$ref)
  
  # has this data been masked?
  if (is.na(y[1])) return(list(amps = NA, crlbs = NA, diags = NA, fit = NA))
  
  # use default fitting opts if not specified 
  if (is.null(opts)) {
      opts <- varpro_3_para_opts()
  }
  
  # convert basis from FD to TD
  basis_td <- apply(basis$data, 2, ift_shift)
  
  y <- drop(mrs_data$data)
  Npts <- length(y)
  Nbasis <- dim(basis$data)[2]
  
  # phase, global damping, global shift
  par <- c(0, opts$init_damping, 0)
           
  t <- seconds(mrs_data)  
  # lm control options
  ctrl <- minpack.lm::nls.lm.control()
  ctrl$maxiter = opts$maxiters
  # do the fit
  lower <- c(-pi, 0, -opts$max_shift)
  upper <- c(pi, opts$max_damping, opts$max_shift)
  
  if (opts$anal_jac) {
    res <- minpack.lm::nls.lm(par, lower, upper, varpro_3_para_obj,
                              varpro_3_para_anal_jac, ctrl, y, basis_td, t,
                              opts$nstart)
  } else {
    res <- minpack.lm::nls.lm(par, lower, upper, varpro_3_para_obj, NULL, ctrl,
                              y, basis_td, t, opts$nstart)
  }
  
  # apply phase to y
  y <- y * exp(1i * (res$par[1]))
  
  # apply global broadening term to basis
  basis_mod <- basis_td * matrix(exp(-t * t * lw2beta(res$par[2])),
                                 ncol = ncol(basis_td), nrow = nrow(basis_td),
                                 byrow = FALSE)
  
  # apply shift terms to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * rep(res$par[3], Nbasis)
  freq_mat <- matrix(freq_vec, nrow = Npts, ncol = Nbasis, byrow = TRUE)
 
  basis_mod <- basis_mod * exp(t_mat * freq_mat)
  
  # get ahat
  y_real <- c(Re(y[opts$nstart:Npts]), Im(y[opts$nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[opts$nstart:Npts,]),
                      Im(basis_mod[opts$nstart:Npts,]))
  
  ahat <- nnls::nnls(basis_real, y_real)$x
  
  yhat <- basis_mod %*% ahat
  amat <- matrix(ahat, nrow = Npts, ncol = Nbasis, byrow = TRUE)
  basis_sc <- basis_mod * amat
  zero_mat <- matrix(0, nrow = Npts, ncol = Nbasis)
  basis_sc <- rbind(basis_sc, zero_mat)
  BASIS_SC <- apply(basis_sc, 2, ft_shift)
  basis_frame <- as.data.frame(Re(BASIS_SC), row.names = NA)
  colnames(basis_frame) <- basis$names
  
  # zero pad
  yhat <- c(yhat, rep(0, Npts))
  YHAT <- ft_shift(as.vector(yhat))
  
  # zero pad
  y <- c(y, rep(0, Npts))
  Y <- ft_shift(y)
  resid <- Y - YHAT
  
  suppressPackageStartupMessages(
  BL <- smoother::smth.gaussian(Re(resid), opts$bl_smth_pts, tails = TRUE) + 
        1i * smoother::smth.gaussian(Im(resid), opts$bl_smth_pts, tails = TRUE)
  )
  
  RESID <- Y - YHAT
  
  offset <- max(Re(Y)) - min(Re(RESID))
  resid <- vec2mrs_data(RESID + offset, fd = TRUE)
  
  amps <- data.frame(t(ahat))
  colnames(amps) <- basis$names
  
  # create some common metabolite combinations
  if (("NAA" %in% colnames(amps)) & ("NAAG" %in% colnames(amps))) {
    amps['TNAA'] <- amps['NAA'] + amps['NAAG']
  }
  
  if (("PCh" %in% colnames(amps)) & ("GPC" %in% colnames(amps))) {
    amps['TCho'] <- amps['PCh'] + amps['GPC']
  }
  
  if (("Cr" %in% colnames(amps)) & ("PCr" %in% colnames(amps))) {
    amps['TCr'] <- amps['Cr'] + amps['PCr']
  }
  
  if (("Glu" %in% colnames(amps)) & ("Gln" %in% colnames(amps))) {
    amps['Glx'] <- amps['Glu'] + amps['Gln']
  }
  
  if (("Lip09" %in% colnames(amps)) & ("MM09" %in% colnames(amps))) {
    amps['TLM09'] <- amps['Lip09'] + amps['MM09']
  }
  
  if (("Lip13a" %in% colnames(amps)) & ("Lip13b" %in% colnames(amps)) & 
        ("MM12" %in% colnames(amps)) & ("MM14" %in% colnames(amps))) {
    amps["TLM13"] <- amps["Lip13a"] + amps["Lip13b"] + amps["MM12"] + 
                     amps["MM14"]
  }
  
  if (("Lip20" %in% colnames(amps)) & ("MM20" %in% colnames(amps))) {
    amps['TLM20'] <- amps['Lip20'] + amps['MM20']
  }
  
  fit <- data.frame(PPMScale = ppm(mrs_data, N = Npts * 2), Data = Re(Y),
                    Fit = Re(YHAT), Baseline = Re(BL))
  
  fit <- cbind(fit, basis_frame)
  
  class(fit) <- c("fit_table", "data.frame")
  
  diags <- data.frame(phase = res$par[1] * 180 / pi, damping = res$par[2],
                      shift = res$par[3], res$deviance, res$niter, res$info, 
                      res$message)
  
  list(amps = amps, crlbs = t(rep(NA, length(amps))), diags = diags, fit = fit)
}

varpro_3_para_obj <- function(par, y, basis, t, nstart, sc_res = FALSE) {
  Npts <- length(y)
  Nbasis <- dim(basis)[2]
  
  # apply phase to y
  y <- y * exp(1i * (par[1]))
  
  # apply global broadening term to basis
  basis_mod <- basis * matrix(exp(-t * t * lw2beta(par[2])), ncol = ncol(basis),
                              nrow = nrow(basis), byrow = F)
  
  # apply global shift to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * rep(par[3], Nbasis)
  freq_mat <- matrix(freq_vec, nrow = Npts, ncol = Nbasis, byrow = TRUE)
 
  basis_mod <- basis_mod * exp(t_mat * freq_mat)
  
  y_real <- c(Re(y[nstart:Npts]), Im(y[nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[nstart:Npts,]), Im(basis_mod[nstart:Npts,]))
  
  ahat <- nnls::nnls(basis_real, y_real)$x
  res <- y_real - basis_real %*% ahat
  
  if ( sc_res ) {
    res <- sum(res ^ 2)
  }
  res
}

# following doc example gives an error:
# @examples
# varpro_3_para_opts(nstart = 10)

#' Return a list of options for VARPRO based fitting with 3 free parameters:
#' * zero'th order phase correction
#' * global damping
#' * global frequency shift.
#' @param nstart position in the time-domain to start fitting, units of data
#' points.
#' @param init_damping starting value for the global Gaussian line-broadening
#' term - measured in Hz.
#' @param maxiters maximum number of levmar iterations to perform.
#' @param max_shift maximum global shift allowed, measured in Hz.
#' @param max_damping maximum damping allowed, FWHM measured in Hz.
#' @param anal_jac option to use the analytic or numerical Jacobian (logical).
#' @param bl_smth_pts number of data points to use in the baseline smoothing
#' calculation.
#' @return list of options.
#' @export
varpro_3_para_opts <- function(nstart = 20, init_damping = 2, maxiters = 200,
                        max_shift = 5, max_damping = 5, anal_jac = FALSE,
                        bl_smth_pts = 80) {
  
  list(nstart = nstart, init_damping = init_damping, maxiters = maxiters,
       max_shift = max_shift, max_damping = max_damping, anal_jac = anal_jac,
       bl_smth_pts = bl_smth_pts)
}

varpro_3_para_anal_jac <- function(par, y, basis, t, nstart) {
  Npts <- length(y)
  Nbasis <- dim(basis)[2]

  # apply phase to y
  y <- y*exp(1i * (par[1]))

  # apply global broadening term to basis
  basis_mod <- basis * matrix(exp(-t * t * lw2beta(par[2])), ncol = ncol(basis),
                              nrow = nrow(basis), byrow = F)

  # apply global shift to basis
  t_mat <- matrix(t, nrow = Npts, ncol = Nbasis)
  freq_vec <- 2i * pi * rep(par[3], Nbasis)
  freq_mat <- matrix(freq_vec, nrow = Npts, ncol = Nbasis, byrow = TRUE)
  basis_mod <- basis_mod * exp(t_mat * freq_mat)

  t_cut <- c(t[nstart:Npts], t[nstart:Npts])
  y_real <- c(Re(y[nstart:Npts]), Im(y[nstart:Npts]))
  basis_real <- rbind(Re(basis_mod[nstart:Npts,]), Im(basis_mod[nstart:Npts,]))
  ahat <- nnls::nnls(basis_real, y_real)$x

  unmod_basis_real <- basis_real %*% ahat
  unmod_basis_cplx <- basis_mod %*% ahat

  phase_jac <- unmod_basis_cplx * 1i
  phase_jac_real <- c(Re(phase_jac[nstart:Npts]), Im(phase_jac[nstart:Npts]))
  g_lw_jac_real <- (unmod_basis_real * t_cut * t_cut * 2 * par[2] * (pi ^ 2) / 4
                    / (-log(0.5)))
  
  shift_jac <- -unmod_basis_cplx * 2i * pi * t
  shift_jac_real <- c(Re(shift_jac[nstart:Npts]), Im(shift_jac[nstart:Npts]))

  c(phase_jac_real, g_lw_jac_real, shift_jac_real)
}