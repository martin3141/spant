#' @export
spant_benchmark <- function(noise_reps = 10, return_res = FALSE) {
  
  # simulate data and basis
  sim_res <- sim_brain_1h(full_output = TRUE)
  
  # apply 5Hz line broadening
  metab <- lb(sim_res$mrs_data, 5)
  
  # make multiple versions with different noise samples
  metab_rep  <- rep_mrs(metab, dyn_rep = noise_reps)
  set.seed(101)
  metab_rep  <- metab_rep + sim_noise(fd = FALSE, dyns = noise_reps, sd = 0.01)
  
  # time the fitting duration
  start_time <- Sys.time()
  fit_res    <- fit_mrs(metab_rep, sim_res$basis)
  end_time   <- Sys.time()
  diff_time  <- as.double(end_time - start_time)
  
  # output results
  cat("Total fit time    :", diff_time, "seconds.\n")
  cat("Time per fit      :", diff_time / noise_reps, "seconds.\n")
  cat("Mean of residuals : ", mean(fit_res$res_tab$res.deviance), ".\n",
      sep = "")
  
  if (return_res) return(fit_res)
}