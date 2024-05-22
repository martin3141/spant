#' Generate trapezoidal regressors.
#' @param onset stimulus onset in seconds.
#' @param duration stimulus duration in seconds.
#' @param trial_type string label for the stimulus.
#' @param mrs_data mrs_data object for timing information.
#' @param rise_t time to reach a plateau from baseline in seconds.
#' @param fall_t time to fall from plateau level back to baseline in seconds.
#' @param exp_fall model an exponential fall instead of linear.
#' @param exp_fall_power exponential fall power.
#' @param smo_sigma standard deviation of Gaussian smoothing kernel in seconds.
#' Set to NULL to disable (default behavior).
#' @param match_tr match the output to the input mrs_data.
#' @param dt timing resolution for internal calculations.
#' @param normalise normalise the response function to have a maximum value of 
#' one.
#' @return trapezoidal regressor data frame.
#' @export
gen_trap_reg <- function(onset, duration, trial_type = NULL, mrs_data,
                         rise_t = 0, fall_t = 0, exp_fall = FALSE,
                         exp_fall_power = 1, smo_sigma = NULL, match_tr = TRUE,
                         dt = 0.01, normalise = FALSE) {
                         
  if (is.na(tr(mrs_data)) | is.null(tr(mrs_data))) {
    stop("TR not set, use set_tr function to set the repetition time.")
  }
  
  if (is.na(Ntrans(mrs_data)) | is.null(Ntrans(mrs_data))) {
    stop("Number of transients not set, use set_Ntrans function to set the 
         number of transients.")
  }
  
  if (is.null(trial_type)) trial_type <- rep("stim", length(onset))
  
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
  
  # time for a 95% reduction
  if (exp_fall) lambda <- -fall_t ^ exp_fall_power / log(0.05)
    
  # loop over trial types
  for (m in 1:trial_type_n) {
    stim_fine <- rep(0, length(t_fine))
    stim_bool <- rep(FALSE, length(t_fine))
    
    # filter out the stim of interest
    stim_frame_trial <- stim_frame[(stim_frame$trial_type == trial_types[m]),]
    
    # loop over stims of the same trial type
    for (n in 1:length(stim_frame_trial$onset)) {
      stim_seg <- t_fine > stim_frame$onset[n] & t_fine <= stim_frame$end[n]
      stim_bool[stim_seg] <- TRUE
    }
    
    last_val <- 0
    for (n in 1:length(stim_fine)) {
      
      if (stim_bool[n]) {
        new_val <- last_val + dt / rise_t
      } else {
        if (exp_fall) {
          t <- (-lambda * log(last_val)) ^ (1 / exp_fall_power) + dt
          new_val <- exp(-(t ^ exp_fall_power / lambda))
          # new_val <-  last_val - dt / lambda * last_val ^ 2
          # new_val <-  last_val - dt / lambda * last_val
      } else {
          new_val <- last_val - dt / fall_t
        }
      }
      
      if (new_val > 1) new_val <- 1
      if (new_val < 0) new_val <- 0
      
      stim_fine[n] <- new_val
      
      last_val <- new_val
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
   
    # correct for missmatch between n_trans and n_dyns due to temporal averaging 
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

#' Generate BOLD regressors.
#' @param onset stimulus onset in seconds.
#' @param duration stimulus duration in seconds.
#' @param trial_type string label for the stimulus.
#' @param mrs_data mrs_data object for timing information.
#' @param match_tr match the output to the input mrs_data.
#' @param dt timing resolution for internal calculations.
#' @return BOLD regressor data frame.
#' @export
gen_bold_reg <- function(onset, duration, trial_type = NULL, mrs_data,
                         match_tr = TRUE, dt = 0.1) {
  
  if (is.na(tr(mrs_data)) | is.null(tr(mrs_data))) {
    stop("TR not set, use set_tr function to set the repetition time.")
  }
  
  if (is.na(Ntrans(mrs_data)) | is.null(Ntrans(mrs_data))) {
    stop("Number of transients not set, use set_Ntrans function to set the 
         number of transients.")
  }
  
  if (is.null(trial_type)) trial_type <- rep("stim_bold", length(onset))
  
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

#' Perform a GLM analysis of dynamic MRS data in the spectral domain.
#' @param mrs_data single-voxel dynamics MRS data.
#' @param regressor_df a data frame containing temporal regressors to be applied
#' to each spectral datapoint.
#' @return list of statistical results.
#' @export
glm_spec <- function(mrs_data, regressor_df) {
  
  # warning, any column named time in regressor_df will be removed
  
  # needs to be a FD operation
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  mrs_mat <- Re(mrs_data2mat(mrs_data))
  
  # drop the time column if present
  regressor_df<- regressor_df[, !names(regressor_df) %in% c("time"),
                              drop = FALSE]
  
  lm_res_list <- vector("list", ncol(mrs_mat))
  for (n in 1:ncol(mrs_mat)) {
    lm_res_list[[n]] <- summary(stats::lm(mrs_mat[, n] ~ ., regressor_df))
  }
  
  # extract stats
  get_glm_stat <- function(x, name) x$coefficients[-1, name, drop = FALSE]
  
  beta_weight <- as.data.frame(t(as.data.frame(sapply(lm_res_list, get_glm_stat,
                                               "Estimate", simplify = FALSE))))
  p_value <- as.data.frame(t(as.data.frame(sapply(lm_res_list, get_glm_stat,
                                           "Pr(>|t|)", simplify = FALSE))))
  
  row.names(beta_weight) <- NULL
  row.names(p_value)     <- NULL
  
  ppm_sc      <-  ppm(mrs_data)
  beta_weight <-  cbind(ppm = ppm_sc, beta_weight)
  p_value_log <- -log10(p_value)
  p_value     <-  cbind(ppm = ppm_sc, p_value)
  p_value_log <-  cbind(ppm = ppm_sc, p_value_log)
  
  p_value_log_mrs <- mat2mrs_data(t(p_value_log[, -1]), fs = fs(mrs_data),
                                  ft = mrs_data$ft, ref = mrs_data$ref,
                                  nuc = mrs_data$nuc, fd = TRUE)
  
  p_value_mrs     <- mat2mrs_data(t(p_value[, -1]), fs = fs(mrs_data),
                                  ft = mrs_data$ft, ref = mrs_data$ref,
                                  nuc = mrs_data$nuc, fd = TRUE)
  
  beta_weight_mrs <- mat2mrs_data(t(beta_weight[, -1]), fs = fs(mrs_data),
                                  ft = mrs_data$ft, ref = mrs_data$ref,
                                  nuc = mrs_data$nuc, fd = TRUE)
  
  return(list(beta_weight = beta_weight, p_value = p_value,
              p_value_log = p_value_log, p_value_log_mrs = p_value_log_mrs,
              p_value_mrs = p_value_mrs, beta_weight_mrs = beta_weight_mrs,
              lm_res_list = lm_res_list))
}

#' Calculate the efficiency of a regressor data frame.
#' @param regressor_df input regressor data frame.
#' @param contrasts a vector of contrast values.
#' @export
calc_design_efficiency <- function(regressor_df, contrasts) {
  X   <- as.matrix(regressor_df[, -1])
  eff <- 1 / sum(diag(t(contrasts) %*% ginv(t(X) %*% X) %*% contrasts))
  return(eff)
}

#' Plot regressors as an image.
#' @param regressor_df input regressor data frame.
#' @export
plot_reg <- function(regressor_df) {
  time  <- regressor_df$time
  names <- colnames(regressor_df)
  X     <- t(regressor_df[, -1])
  graphics::image(y = time, z = X, col = viridisLite::viridis(128),
                  ylab = "Time (s)")
}

#' Append multiple regressor data frames into a single data frame.
#' @param ... input regressor data frames.
#' @return output regressor data frame.
#' @export
append_regs <- function(...) {
  df_list    <- list(...)
  time       <- df_list[[1]]$time
  df_list_nt <- lapply(df_list, subset, select = -time)
  output     <- do.call("cbind", df_list_nt)
  return(cbind(time, output)) 
}

#' Generate baseline regressor.
#' @param mrs_data mrs_data object for timing information.
#' @return a single baseline regressor with value of 1.
#' @export
gen_baseline_reg <- function(mrs_data) {
  time   <- dyn_acq_times(mrs_data)
  reg_df <- data.frame(time = time, baseline = rep(1, length(t)))
  return(reg_df)
}

#' Generate polynomial regressors.
#' @param mrs_data mrs_data object for timing information.
#' @param degree the degree of the polynomial.
#' @return polynomial regressors.
#' @export
gen_poly_reg <- function(mrs_data, degree) {
  time       <- dyn_acq_times(mrs_data)
  poly_mat   <- stats::poly(time, degree)
  scale_vals <- apply(Mod(poly_mat), 2, max)
  poly_mat   <- scale(poly_mat, center = FALSE, scale = scale_vals)
  reg_df     <- data.frame(time = time, poly = poly_mat)
  return(reg_df)
}

#' Generate impulse regressors.
#' @param onset stimulus onset in seconds.
#' @param trial_type string label for the stimulus.
#' @param mrs_data mrs_data object for timing information.
#' @return impulse regressors data frame.
#' @export
gen_impulse_reg <- function(onset, trial_type = NULL, mrs_data) {
  
  if (is.na(tr(mrs_data)) | is.null(tr(mrs_data))) {
    stop("TR not set, use set_tr function to set the repetition time.")
  }
  
  if (is.na(Ntrans(mrs_data)) | is.null(Ntrans(mrs_data))) {
    stop("Number of transients not set, use set_Ntrans function to set the 
         number of transients.")
  }
  
  if (is.null(trial_type)) trial_type <- rep("stim_imp", length(onset))
  
  trial_types  <- unique(trial_type)
  trial_type_n <- length(trial_types)
  
  stim_frame <- data.frame(onset, trial_type)
 
  n_dyns    <- Ndyns(mrs_data)
  empty_mat <- matrix(NA, nrow = n_dyns, ncol = trial_type_n)
  
  output_frame <- data.frame(empty_mat)
  colnames(output_frame) <- c(trial_types)
    
  time <- dyn_acq_times(mrs_data)
  
  for (m in 1:trial_type_n) {
    stim <- rep(0, length(time))
    stim_frame_trial <- stim_frame[(stim_frame$trial_type == trial_types[m]),]
    for (n in 1:length(stim_frame_trial$onset)) {
      ind <- which.min(Mod(time - stim_frame_trial$onset[n]))
      stim[ind] <- 1
      if (Mod(stim_frame_trial$onset[n] - time[ind]) > 0.01) {
        warning("onset and output impulse differ by more than 10 ms")
      }
    }
    output_frame[, m] <- stim
  }
  output_frame <- cbind(time, output_frame)
  
  return(output_frame)
}