#' Correct linear frequency drift.
#' @param mrs_data MRS data to be corrected.
#' @param max_hz_s the maximum drift rate to search over.
#' @param tr mrs_data repetition time.
#' @param ret_corr_only return the corrected mrs_data object only.
#' @param outlier_thresh threshold to remove outliers.
#' @param xlim spectral width (in ppm) to evaluate outliers.
#' @param order correction order.
#' @return drift corrected mrs_data object.
#' @export
lofdc <- function(mrs_data, max_hz_s = 0.1, tr = NULL, ret_corr_only = TRUE, 
                  outlier_thresh = 3, xlim = c(4, 0.5), order = 1) {
  
  if (inherits(mrs_data, "list")) {
    
    res <- lapply(mrs_data, lofdc, max_hz_s = max_hz_s, tr = tr,
                  ret_corr_only = ret_corr_only,
                  outlier_thresh = outlier_thresh, xlim = xlim, order = order)
                  
    return(res)
  }
  
  # covert to time-domain
  if (is_fd(mrs_data)) mrs_data <- fd2td(mrs_data)
  
  if (is.null(tr)) {
    if (is.null(mrs_data$meta$RepetitionTime)) {
      stop("drift_corr function requires the TR.")
    } else {
      tr <- mrs_data$meta$RepetitionTime
    }
  }
  
  max_hz <- max_hz_s * tr * Ndyns(mrs_data)
  
  # generate the shift matrix
  t_pts <- Ndyns(mrs_data)
  drift_mat <- stats::poly(0:(t_pts - 1), degree = order, simple = TRUE)
  
  # make sure t = 0 has a shift of 0 Hz
  drift_mat <- drift_mat - matrix(drift_mat[1,], nrow = nrow(drift_mat),
                                  ncol = ncol(drift_mat), byrow = TRUE)
  
  # scale matrix columns to have max val of 1
  for (n in 1:ncol(drift_mat)) {
    drift_mat[, n] <- drift_mat[, n] / max(Mod(drift_mat[, n]))
  }
  
  if (order == 1) {
    optim_res <- stats::optim(0, lofdc_obj_fn, method = "Brent",
                              lower = -max_hz, upper = max_hz,
                              mrs_data = mrs_data, drift_mat = drift_mat,
                              control = list(fnscale = -1))
  } else {
    start_vals <- rep(0, order)
    optim_res  <- stats::optim(start_vals, lofdc_obj_fn, method = "Nelder-Mead",
                               mrs_data = mrs_data, drift_mat = drift_mat,
                               control = list(fnscale = -1))
  }
  
  hz_s       <- optim_res$par[1] / (tr * Ndyns(mrs_data))
  
  drift_vec     <- drift_mat %*% optim_res$par
  mrs_data_corr <- shift(mrs_data, drift_vec, units = "hz")
  
  drift_range_hz <- max(drift_vec) - min(drift_vec)
  
  if (is.null(outlier_thresh)) {
    mrs_data_corr_masked <- mrs_data_corr
  } else {
    # detect and mask outliers
    mean_spec <- mean_dyns(mrs_data_corr)
    
    resids <- mrs_data_corr - rep_mrs(mean_spec, dyn_rep = Ndyns(mrs_data_corr))
    
    int <- int_spec(resids, xlim = xlim, mode = "mod")
    
    resids_med <- stats::median(int)
    resids_mad <- stats::mad(int)
    
    upper_lim <- resids_med + resids_mad * outlier_thresh
    lower_lim <- resids_med - resids_mad * outlier_thresh
    
    good <- (int < upper_lim) & (int > lower_lim)
    
    mrs_data_corr_masked <- mask_dyns(mrs_data_corr, drop(!good))
  }
  
  if (ret_corr_only) {
    return(mrs_data_corr_masked)
  } else {
    perc_outlier <- sum(!good) / Ndyns(mrs_data_corr) * 100
    return(list(mrs_data_corr_masked = mrs_data_corr_masked,
                mrs_data_corr = mrs_data_corr, hz_s = hz_s,
                perc_outlier = perc_outlier, drift_mat = drift_mat,
                drift_vec = drift_vec, drift_range_hz = drift_range_hz))
  }
}

lofdc_obj_fn <- function(par, mrs_data, drift_mat) {
  drift_vec <- drift_mat %*% par
  mrs_data  <- shift(mrs_data, drift_vec, units = "hz")
  val       <- sum(Mod(mean_dyns(mrs_data)$data))
  return(val)
}