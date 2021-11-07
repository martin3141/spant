#' Simulate and fit some spectra with ABfit for benchmarking purposes. Basic
#' timing and performance metrics will be printed.
#' @param noise_reps number of spectra to fit with differing noise samples.
#' @param return_res return a list of fit_result objects.
#' @param opts ABfit options structure.
#' @export
spant_abfit_benchmark <- function(noise_reps = 10, return_res = FALSE,
                                  opts = abfit_opts()) {
  
  # simulate data and basis
  sim_res <- sim_brain_1h(full_output = TRUE)
  
  # apply 5Hz line broadening
  metab <- lb(sim_res$mrs_data, 5)
  
  # add a residual water signal
  metab <- metab + sim_resonances(4.65, 3, 5, 90)
  
  # apply 120 deg. phase error
  metab <- phase(metab, 120)
 
  # shift by 10 Hz
  metab <- shift(metab, 10, units = "hz")
  
  # make multiple versions with different noise samples
  metab_rep  <- rep_mrs(metab, dyn_rep = noise_reps)
  set.seed(101)
  metab_rep  <- metab_rep + sim_noise(fd = FALSE, dyns = noise_reps, sd = 0.01)
  
  # time the fitting duration
  start_time <- Sys.time()
  fit_res    <- fit_mrs(metab_rep, sim_res$basis, opts = opts)
  end_time   <- Sys.time()
  diff_time  <- as.double(end_time - start_time, units = "secs")
  
  # output results
  cat("Total fit time    :", diff_time, "seconds.\n")
  cat("Time per fit      :", diff_time / noise_reps, "seconds.\n")
  cat("Mean of residuals : ", mean(fit_res$res_tab$res.deviance), ".\n",
      sep = "")
  
  if (return_res) return(fit_res)
}

#' Simulate a typical metabolite basis set for benchmarking. Timing metrics will
#' be printed on completion.
#' @param sim_reps number of times to simulate the basis set.
#' @param N number of FID data points to simulate.
#' @export
spant_simulation_benchmark <- function(sim_reps = 10, N = 1024) {
  
  start_time <- Sys.time()
  
  for (n in 1:sim_reps) {
    sim_basis_1h_brain_press(acq_paras = def_acq_paras(N = N)) 
  }
  
  end_time <- Sys.time()
  
  diff_time  <- as.double(end_time - start_time, units = "secs")
  cat("Total simulation time :", diff_time, "seconds.\n")
}