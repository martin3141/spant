varpro_basic <- function(y, acq_paras, basis, opts = NULL) {
  mrs_data <- vec2mrs_data(y, fs = acq_paras$fs, ft = acq_paras$ft, 
                           ref = acq_paras$ref)
  
  # has this data been masked?
  if (is.na(y[1])) return(list(amps = NA, crlbs = NA, diags = NA, fit = NA))
  
  # use default fitting opts if not specified 
  if (is.null(opts)) opts <- varpro_basic_opts()
  
  y <- drop(mrs_data$data)
  
  
  if (opts$method == "td") {
    proj <- td_projection(y, basis, opts$nnls)
  } else if (opts$method == "fd") {
    proj <- fd_projection(y, basis, opts$nnls)
  } else if (opts$method == "fd_re") {
    proj <- fd_re_projection(y, basis, opts$nnls)
  } else {
    stop("unrecognised varpro method")
  }

  # create some common metabolite combinations
  amps <- append_metab_combs(proj$amps)
  
  fit <- data.frame(PPMScale = ppm(mrs_data, N = 2 * length(y)),
                    Data = Re(proj$Y), Fit = Re(proj$YHAT),
                    Baseline = Re(proj$BL))
  
  fit <- cbind(fit, proj$basis_frame)
  
  class(fit) <- c("fit_table", "data.frame")
  
  diags <- data.frame(dummy = 1, stringsAsFactors = TRUE)
  
  list(amps = amps, crlbs = t(rep(NA, length(amps))), diags = diags, fit = fit)
}

fd_projection <- function(y, basis, nnls) {
  
  Npts   <- length(y)
  Nbasis <- dim(basis$data)[2]
  
  Y <- ft_shift(y)
  Y_real <- c(Re(Y), Im(Y))
  
  basis_fd   <- basis$data
  basis_real <- rbind(Re(basis_fd), Im(basis_fd))
  
  if (nnls) {
    ahat <- nnls(basis_real, Y_real)$x
  } else {
    ahat <- stats::.lm.fit(basis_real, Y_real)$coefficients
  }
  
  YHAT <- basis_fd %*% ahat
  yhat <- ift_shift(as.vector(YHAT))
  amat <- matrix(ahat, nrow = Npts, ncol = Nbasis, byrow = TRUE)
  basis_sc <- basis_fd * amat
  basis_sc <- apply(basis_sc, 2, ift_shift)
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
  BL <- rep(0, length(Y))
  
  amps <- data.frame(t(ahat))
  colnames(amps) <- basis$names
  
  return(list(amps = amps, Y = Y, YHAT = YHAT, BL = BL,
              basis_frame = basis_frame))
}

fd_re_projection <- function(y, basis, nnls) {
  
  Npts   <- length(y)
  Nbasis <- dim(basis$data)[2]
  
  # zf y 
  y <- c(y, rep(0, Npts))
  
  Y <- ft_shift(y)
  Y_real <- Re(Y)
  
  basis_td <- apply(basis$data, 2, ift_shift)
  zero_mat <- matrix(0, nrow = Npts, ncol = Nbasis)
  basis_td <- rbind(basis_td, zero_mat)
  basis_fd <- apply(basis_td, 2, ft_shift)
  
  basis_real <- Re(basis_fd)
  
  if (nnls) {
    ahat <- nnls(basis_real, Y_real)$x
  } else {
    ahat <- stats::.lm.fit(basis_real, Y_real)$coefficients
  }
  
  YHAT <- basis_real %*% ahat
  amat <- matrix(ahat, nrow = Npts * 2, ncol = Nbasis, byrow = TRUE)
  BASIS_SC <- basis_real * amat
  
  basis_frame <- as.data.frame(Re(BASIS_SC), row.names = NA)
  colnames(basis_frame) <- basis$names
  
  BL <- rep(0, length(Y))
  
  amps <- data.frame(t(ahat))
  colnames(amps) <- basis$names
  
  return(list(amps = amps, Y = Y, YHAT = YHAT, BL = BL,
              basis_frame = basis_frame))
}

td_projection <- function(y, basis, nnls) {
  
  Npts   <- length(y)
  Nbasis <- dim(basis$data)[2]
  
  y_real <- c(Re(y), Im(y))
  
  basis_td   <- apply(basis$data, 2, ift_shift)
  basis_real <- rbind(Re(basis_td), Im(basis_td))
  
  if (nnls) {
    ahat <- nnls(basis_real, y_real)$x
  } else {
    ahat <- stats::.lm.fit(basis_real, y_real)$coefficients
  }
  
  yhat <- basis_td %*% ahat
  amat <- matrix(ahat, nrow = Npts, ncol = Nbasis, byrow = TRUE)
  basis_sc <- basis_td * amat
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
  BL <- rep(0, length(Y))
  
  amps <- data.frame(t(ahat))
  colnames(amps) <- basis$names
  
  return(list(amps = amps, Y = Y, YHAT = YHAT, BL = BL,
              basis_frame = basis_frame))
}

append_metab_combs <- function(amps) {
    
  if (("NAA" %in% colnames(amps)) & ("NAAG" %in% colnames(amps))) {
    amps['tNAA'] <- amps['NAA'] + amps['NAAG']
  }
  
  if (("PCh" %in% colnames(amps)) & ("GPC" %in% colnames(amps))) {
    amps['tCho'] <- amps['PCh'] + amps['GPC']
  }
  
  if (("Cr" %in% colnames(amps)) & ("PCr" %in% colnames(amps))) {
    amps['tCr'] <- amps['Cr'] + amps['PCr']
  }
  
  if (("Glu" %in% colnames(amps)) & ("Gln" %in% colnames(amps))) {
    amps['Glx'] <- amps['Glu'] + amps['Gln']
  }
  
  if (("Lip09" %in% colnames(amps)) & ("MM09" %in% colnames(amps))) {
    amps['tLM09'] <- amps['Lip09'] + amps['MM09']
  }
  
  if (("Lip13a" %in% colnames(amps)) & ("Lip13b" %in% colnames(amps)) & 
        ("MM12" %in% colnames(amps)) & ("MM14" %in% colnames(amps))) {
    amps["tLM13"] <- amps["Lip13a"] + amps["Lip13b"] + amps["MM12"] + 
                     amps["MM14"]
  }
  
  if (("Lip20" %in% colnames(amps)) & ("MM20" %in% colnames(amps))) {
    amps['tLM20'] <- amps['Lip20'] + amps['MM20']
  }
  
  return(amps)
}

#' Return a list of options for a basic VARPRO analysis.
#' 
#' @param method one of "td", "fd", "fd_re".
#' @param nnls restrict basis amplitudes to non-negative values.
#' @return full list of options.
#' @export
varpro_basic_opts <- function(method = "fd_re", nnls = TRUE) {
  list(method = method, nnls = nnls)
}