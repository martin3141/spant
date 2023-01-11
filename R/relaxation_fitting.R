#' Fit a T2 relaxation curve, from multiple TEs, to a set of amplitudes.
#' @param te_vec vector of TE values in seconds.
#' @param amp_vec vector of amplitudes.
#' @param lower minimum allowable T2 value.
#' @param upper maximum allowable T2 value.
#' @param output_fit_res temporal resolution (seconds) of the ideal output
#' relaxation curve.
#' @param ret_full return full fitting information including ideal relaxation
#' curve.
#' @return a list containing relaxation parameters and an ideal curve for fit
#' evaluation.
#' @export
fit_t2_te_array <- function(te_vec, amp_vec, lower = 0, upper = 10,
                            output_fit_res = 0.01, ret_full = TRUE) {
  
  fit_res <- stats::optimize(t2_te_array_obj, te_vec = te_vec,
                             amp_vec = amp_vec, lower = lower, upper = upper)
  
  t2 <- fit_res$minimum
  
  model <- exp(-te_vec / t2)
  
  c <- as.numeric(stats::lm(amps ~ model + 0)$coefficients)
 
  if (ret_full) {
    t <- seq(from = min(te_vec), to = max(te_vec), by = output_fit_res)
    fit_amps <- c * exp(-t / t2)
    fit_table <- data.frame(t = t, amps = fit_amps)
    return(list(t2 = t2, c = c, fit_table = fit_table))
  } else {
    return(list(t2 = t2, c = c))
  }
}

t2_te_array_obj <- function(par, te_vec, amp_vec) {
  model <- exp(-te_vec / par)
  c <- as.numeric(stats::lm(amp_vec ~ model + 0)$coefficients)
  return(sum((amp_vec - c * model) ^ 2))
}

#' Fit a T1 recovery curve, from multiple TRs, to a set of amplitudes.
#' @param tr_vec vector of TR values in seconds.
#' @param amp_vec vector of amplitudes.
#' @param lower minimum allowable T1 value.
#' @param upper maximum allowable T1 value.
#' @param output_fit_res temporal resolution (seconds) of the ideal output
#' relaxation curve.
#' @param ret_full return full fitting information including ideal relaxation
#' curve.
#' @return a list containing relaxation parameters and an ideal curve for fit
#' evaluation.
#' @export
fit_t1_tr_array <- function(tr_vec, amp_vec, lower = 0, upper = 10,
                            output_fit_res = 0.01, ret_full = TRUE) {
  
  fit_res <- stats::optimize(t1_tr_array_obj, tr_vec = tr_vec,
                             amp_vec = amp_vec, lower = lower, upper = upper)
  
  t1 <- fit_res$minimum
  
  model <- (1 - exp(-tr_vec / t1))
  
  c <- as.numeric(stats::lm(amps ~ model + 0)$coefficients)
 
  if (ret_full) {
    t <- seq(from = min(tr_vec), to = max(tr_vec), by = output_fit_res)
    fit_amps <- c * (1 - exp(-t / t1))
    fit_table <- data.frame(t = t, amps = fit_amps)
    return(list(t1 = t1, c = c, fit_table = fit_table))
  } else {
    return(list(t1 = t1, c = c))
  }
}

t1_tr_array_obj <- function(par, tr_vec, amp_vec) {
  model <- 1 - exp(-tr_vec / par)
  c <- as.numeric(stats::lm(amp_vec ~ model + 0)$coefficients)
  return(sum((amp_vec - c * model) ^ 2))
}

#' Fit a T1 recovery curve, from multiple TIs, to a set of amplitudes.
#' @param ti_vec vector of TI values in seconds.
#' @param amp_vec vector of amplitudes.
#' @param lower minimum allowable T1 value.
#' @param upper maximum allowable T1 value.
#' @param output_fit_res temporal resolution (seconds) of the ideal output
#' relaxation curve.
#' @param ret_full return full fitting information including ideal relaxation
#' curve.
#' @return a list containing relaxation parameters and an ideal curve for fit
#' evaluation.
#' @export
fit_t1_ti_array <- function(ti_vec, amp_vec, lower = 0, upper = 10,
                            output_fit_res = 0.01, ret_full = TRUE) {
  
  fit_res <- stats::optimize(t1_ti_array_obj, ti_vec = ti_vec,
                             amp_vec = amp_vec, lower = lower, upper = upper)
  
  t1 <- fit_res$minimum
  
  model <- (1 - 2 * exp(-ti_vec / t1))
  
  c <- as.numeric(stats::lm(amps ~ model + 0)$coefficients)
 
  if (ret_full) {
    t <- seq(from = min(ti_vec), to = max(ti_vec), by = output_fit_res)
    fit_amps <- c * (1 - 2 * exp(-t / t1))
    fit_table <- data.frame(t = t, amps = fit_amps)
    return(list(t1 = t1, c = c, fit_table = fit_table))
  } else {
    return(list(t1 = t1, c = c))
  }
}

t1_ti_array_obj <- function(par, ti_vec, amp_vec) {
  model <- 1 - 2 * exp(-ti_vec / par)
  c <- as.numeric(stats::lm(amp_vec ~ model + 0)$coefficients)
  return(sum((amp_vec - c * model) ^ 2))
}