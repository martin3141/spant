#' @export
gen_mpress_gaba_report <- function(...) {
  anal_results <- pipe_mpress_gaba(...)
  rep_file <- system.file("reports", "mpress_gaba.Rmd", package = "spant")
  rmarkdown::render(rep_file, output_dir = "~")
}

#' @export
pipe_mpress_gaba <- function(mrs_data, TE1 = 30, TE2 = 40) {
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
  fit <- fit_mrs(mean_aligned, basis, method = "VARPRO_3P")
  
  mean_phase <- fit$results$phase
  
  mean_aligned <- phase(mean_aligned, mean_phase)
  
  mrs_data_aligned <- phase(mrs_data_aligned, mean_phase)
  
  mean_odd_dyns <- hsvd_filt(mean_dyns(get_odd_dyns(mrs_data_aligned)))
  mean_even_dyns <- hsvd_filt(mean_dyns(get_even_dyns(mrs_data_aligned)))
  
  odd_int <- as.numeric(int_spec(mean_odd_dyns, xlim = c(1.9, 2.1)))
  even_int <- as.numeric(int_spec(mean_even_dyns, xlim = c(1.9, 2.1)))
  
  if (odd_int > even_int) {
    ed_on  <- mean_even_dyns
    ed_off <- mean_odd_dyns 
  } else {
    ed_on  <- mean_odd_dyns
    ed_off <- mean_even_dyns
  }
  
  ed_gaba <- ed_on - ed_off
  
  ed_off_fit <- fit_mrs(ed_off, basis, method = "VARPRO",
                        opts = varpro_opts(anal_jac = TRUE))
  
  list(shifts = shifts, ed_off_fit = ed_off_fit, ed_on = ed_on, ed_off = ed_off,
       ed_gaba = ed_gaba)
}