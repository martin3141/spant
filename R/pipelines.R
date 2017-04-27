#' @export
pipe_mpress_gaba <- function(mrs_data, TE1 = 35, TE2 = 35) {
  dynamics <- dyns(mrs_data)
  
  # guess the phasing of the acquisition scheme
  phases <- as.numeric(fp_phase(mrs_data))
  odd_phase <- mean(phases[seq(1, to = dynamics, 2)])
  even_phase <- mean(phases[seq(2, to = dynamics, 2)])
  phase_diff <- abs(odd_phase - even_phase)
  
  # stop if phases don't look right 
  phase_tol <- 10
  if (phase_diff < phase_tol) {
    inv_dyns <- FALSE
  } else if (abs(phase_diff - 180) < phase_tol) {
    inv_dyns <- TRUE
  } else {
    stop("Unexpected data phases, check for problems.")
  }
  
  # lets have the same phase for each dynamic
  if (inv_dyns)  {
    mrs_data <- inv_odd_dyns(mrs_data)
  }
  
  # freq_align dynamics using the residual water peak
  align_res <- align(mrs_data, ret_df = TRUE)
  mrs_data_aligned <- align_res[[1]]
  shifts <- as.numeric(align_res[[2]])
  
  # take an average
  mean_aligned <- mean_dyns(mrs_data_aligned)
  
  # remove residual water
  mean_aligned <- hsvd_filt(mean_aligned)
  
  # simulate a brain basis set
  acq_paras <- get_acq_paras(mrs_data)
  basis <- sim_basis_1h_brain_press(acq_paras = acq_paras, TE1 = TE1, TE2 = TE2)
  fit <- fit_mrs(mean_aligned, basis, method = "VARPRO")
  
  list(shifts = shifts, fit = fit)
}