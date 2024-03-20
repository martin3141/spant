#' Generate a trapezoidal response function.
#' @param onset stimulus onset in seconds.
#' @param duration stimulus duration in seconds.
#' @param trial_type string label for the stimulus.
#' @param mrs_data mrs_data object for timing information.
#' @param rise_t rise time of the trapezoid in seconds.
#' @param fall_t fall time of the trapezoid in seconds
#' @param smo_sigma standard deviation of Gaussian smoothing kernel in seconds.
#' Set to NULL to disable (default behavior).
#' @param match_tr match the output to the input mrs_data.
#' @param dt timing resolution for internal calculations.
#' @param normalise normalise the response function to have a maximum value of 
#' one.
#' @return trapezoidal response function.
#' @export
gen_trap_rf <- function(onset, duration, trial_type, mrs_data, rise_t = 0,
                        fall_t = 0, smo_sigma = NULL, match_tr = TRUE,
                        dt = 0.01, normalise = TRUE) {
                         
  if (is.na(tr(mrs_data)) | is.null(tr(mrs_data))) {
    stop("TR not set, use set_tr function to set the repetition time.")
  }
  
  if (is.na(Ntrans(mrs_data)) | is.null(Ntrans(mrs_data))) {
    stop("Number of transients not set, use set_Ntrans function to set the 
         number of transients.")
  }
  
  # check everything is the right length 
  input_lengths <- c(length(onset), length(duration), length(trial_type))
  if (length(unique(input_lengths)) != 1) {
    stop("Stim length input error.")
  }
  
  # make a time scale with dt seconds resolution for the duration of the scan
  # time
  n_trans   <- mrs_data$meta$NumberOfTransients
  TR        <- tr(mrs_data)
  n_dyns    <- Ndyns(mrs_data)
  t_fine    <- seq(from = 0, to = n_trans * TR - TR, by = dt)
  end       <- onset + duration
  
  stim_frame <- data.frame(onset, end, trial_type)
  
  trial_types  <- unique(trial_type)
  trial_type_n <- length(trial_types)
  
  if (match_tr) {
    empty_mat <- matrix(NA, nrow = n_dyns, ncol = trial_type_n + 1)
  } else {
    empty_mat <- matrix(NA, nrow = length(t_fine), ncol = trial_type_n + 1)
  }
  
  output_frame <- data.frame(empty_mat)
  colnames(output_frame) <- c("time", trial_types)
    
  # loop over trial types
  for (m in 1:trial_type_n) {
    stim_fine <- rep(0, length(t_fine))
    
    # filter out the stim of interest
    stim_frame_trial <- stim_frame[(stim_frame$trial_type == trial_types[m]),]
    
    # loop over stims of the same trial type
    for (n in 1:length(stim_frame_trial$onset)) {
      
      if (rise_t > (stim_frame$end[n] - stim_frame$onset[n])) {
        stop("Rise time cannot be shorter than duration.")
      }
      
      rise_inds <- (t_fine > stim_frame$onset[n]) &
                   (t_fine <= (stim_frame$onset[n] + rise_t))
      
      rise_pts <- seq(from = 0, to = 1, length.out = 1 + sum(rise_inds))
      
      if (sum(stim_fine[rise_inds]) > 0) {
        stop("Overlapping response functions (rise)")
      }
      
      stim_fine[rise_inds] <- rise_pts[2:length(rise_pts)]
      
      plateau_inds <- (t_fine > (stim_frame$onset[n] + rise_t)) &
                      (t_fine <= stim_frame$end[n])
      
      if (sum(stim_fine[plateau_inds]) > 0) {
        stop("Overlapping response functions (plateau)")
      }
      
      stim_fine[plateau_inds] <- rep(1, sum(plateau_inds))
      
      fall_inds <- (t_fine > (stim_frame$end[n])) &
                   (t_fine <= (stim_frame$end[n] + fall_t))
      
      if (sum(stim_fine[fall_inds]) > 0) {
        stop("Overlapping response functions (fall)")
      }
      
      fall_pts <- seq(from = 1, to = 0, length.out = 1 + sum(fall_inds))
      stim_fine[fall_inds] <- fall_pts[2:length(fall_pts)]
    }
    
    if (!is.null(smo_sigma)) {
      # generate a 1D Gaussian kernel 
      gaus_ker  <- mmand::gaussianKernel(smo_sigma / dt)
      stim_fine <- mmand::morph(stim_fine, gaus_ker, operator = "*",
                                merge = "sum")
    }
    
    if (normalise) stim_fine <- stim_fine / max(stim_fine)
    
    t_acq    <- seq(from = 0, by = TR, length.out = n_trans)
    stim_acq <- stats::approx(t_fine, stim_fine, t_acq, method='linear')$y
    
    if (n_trans != n_dyns) {
      if (n_trans%%n_dyns != 0) stop("Dynamics and transients do not match")
      
      block_size <- n_trans / n_dyns
      
      t_acq    <- colMeans(matrix(t_acq, nrow = block_size))
      stim_acq <- colMeans(matrix(stim_acq, nrow = block_size))
    }
    
    if (match_tr) {
      if (m == 1) output_frame[, 1] <- t_acq
      output_frame[, (1 + m)] <- stim_acq
    } else {
      if (m == 1) output_frame[, 1] <- t_fine
      output_frame[, (1 + m)] <- stim_fine
    }
  }
  
  return(output_frame)
}

#' Generate a BOLD response function.
#' @param onset stimulus onset in seconds.
#' @param duration stimulus duration in seconds.
#' @param trial_type string label for the stimulus.
#' @param mrs_data mrs_data object for timing information.
#' @param match_tr match the output to the input mrs_data.
#' @param dt timing resolution for internal calculations.
#' @return BOLD response function.
#' @export
gen_bold_rf <- function(onset, duration, trial_type, mrs_data, match_tr = TRUE,
                        dt = 0.1) {
  
  if (is.na(tr(mrs_data)) | is.null(tr(mrs_data))) {
    stop("TR not set, use set_tr function to set the repetition time.")
  }
  
  if (is.na(Ntrans(mrs_data)) | is.null(Ntrans(mrs_data))) {
    stop("Number of transients not set, use set_Ntrans function to set the 
         number of transients.")
  }
  
  # check everything is the right length 
  input_lengths <- c(length(onset), length(duration), length(trial_type))
  if (length(unique(input_lengths)) != 1) {
    stop("Stim length input error.")
  }
  
  # make a time scale with dt seconds resolution for the duration of the scan
  # time
  n_trans   <- mrs_data$meta$NumberOfTransients
  TR        <- tr(mrs_data)
  n_dyns    <- Ndyns(mrs_data)
  t_fine    <- seq(from = 0, to = n_trans * TR - TR, by = dt)
  end       <- onset + duration
  
  stim_frame <- data.frame(onset, end, trial_type)
  
  trial_types  <- unique(trial_type)
  trial_type_n <- length(trial_types)
  
  if (match_tr) {
    empty_mat <- matrix(NA, nrow = n_dyns, ncol = trial_type_n + 1)
  } else {
    empty_mat <- matrix(NA, nrow = length(t_fine), ncol = trial_type_n + 1)
  }
  
  output_frame <- data.frame(empty_mat)
  colnames(output_frame) <- c("time", trial_types)
  
  for (m in 1:trial_type_n) {
    stim_fine <- rep(0, length(t_fine))
    
    # filter out the stim of interest
    stim_frame_trial <- stim_frame[(stim_frame$trial_type == trial_types[m]),]
    
    for (n in 1:length(stim_frame_trial$onset)) {
      stim_fine[t_fine >= stim_frame_trial$onset[n] & 
                t_fine < stim_frame_trial$end[n]] <- 1
    }
    
    resp_fn   <- gen_hrf(res_t = dt)$hrf
    stim_fine <- stats::convolve(stim_fine, rev(resp_fn), type = 'open')
    stim_fine <- stim_fine[1:length(t_fine)]
    
    t_acq    <- seq(from = 0, by = TR, length.out = n_trans)
    stim_acq <- stats::approx(t_fine, stim_fine, t_acq, method='linear')$y
    
    if (n_trans != n_dyns) {
      if (n_trans%%n_dyns != 0) stop("Dynamics and transients do not match")
      
      block_size <- n_trans / n_dyns
      
      t_acq    <- colMeans(matrix(t_acq, nrow = block_size))
      stim_acq <- colMeans(matrix(stim_acq, nrow = block_size))
    }
    
    if (match_tr) {
      if (m == 1) output_frame[, 1] <- t_acq
      output_frame[, (1 + m)] <- stim_acq
    } else {
      if (m == 1) output_frame[, 1] <- t_fine
      output_frame[, (1 + m)] <- stim_fine
    }
    
  }
  
  return(output_frame)
}

# gen double gamma model of hrf (as used in SPM) with 10ms resolution
# https://github.com/spm/spm12/blob/main/spm_hrf.m
gen_hrf <- function(end_t = 30, res_t = 0.01) {
  t_hrf <- seq(from = 0, to = end_t, by = res_t)
  a1 <- 6; a2 <- 16; b1 <- 1; b2 <- 1; c <- 1/6
  hrf <-     t_hrf ^ (a1 - 1) * b1 ^ a1 * exp(-b1 * t_hrf) / gamma(a1) -
         c * t_hrf ^ (a2 - 1) * b2 ^ a2 * exp(-b2 * t_hrf) / gamma(a2)
  hrf <- hrf / sum(hrf)
  return(list(hrf = hrf, t = t_hrf))
}