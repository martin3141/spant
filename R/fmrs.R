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
gen_trap_reg <- function(onset, duration, trial_type = NULL, mrs_data = NULL,
                         rise_t = 0, fall_t = 0, exp_fall = FALSE,
                         exp_fall_power = 1, smo_sigma = NULL, match_tr = TRUE,
                         dt = 0.01, normalise = FALSE) {
  
  if (is.null(mrs_data)) {
    seq_tr   <- 2
    N_scans  <- 800
    mrs_data <- sim_resonances()
    mrs_data <- set_tr(mrs_data, seq_tr)
    mrs_data <- set_Ntrans(mrs_data, N_scans)
    mrs_data <- rep_dyn(mrs_data, N_scans)
  }
                         
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
    print(input_lengths)
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
    
    if (normalise) stim_acq <- stim_acq / max(stim_acq)
   
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
#' @param normalise normalise the response function to have a maximum value of 
#' one.
#' @return BOLD regressor data frame.
#' @export
gen_bold_reg <- function(onset, duration = NULL, trial_type = NULL,
                         mrs_data = NULL, match_tr = TRUE, dt = 0.1,
                         normalise = FALSE) {
  
  # create a dummy dataset if not specified
  if (is.null(mrs_data)) {
    seq_tr   <- 2
    N_scans  <- 800
    mrs_data <- sim_resonances()
    mrs_data <- set_tr(mrs_data, seq_tr)
    mrs_data <- set_Ntrans(mrs_data, N_scans)
    mrs_data <- rep_dyn(mrs_data, N_scans)
  }
  
  if (is.null(duration)) duration <- rep(dt, length(onset))
  
  # set the minimum duration to dt * 1.1
  min_dur <- dt * 1.1
  duration[duration < min_dur] <- min_dur
  
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
  if (length(unique(input_lengths)) != 1) stop("Stim length input error.")
  
  # make a time scale with dt seconds resolution for the duration of the scan
  # time
  n_trans   <- mrs_data$meta$NumberOfTransients
  TR        <- tr(mrs_data)
  n_dyns    <- Ndyns(mrs_data)
  t_fine    <- seq(from = 0, to = n_trans * TR, by = dt)
  end       <- onset + duration
  
  stim_frame   <- data.frame(onset, end, trial_type, duration)
  
  trial_types  <- unique(trial_type)
  trial_type_n <- length(trial_types)
  
  if (match_tr) {
    empty_mat <- matrix(NA, nrow = n_dyns, ncol = trial_type_n + 1)
  } else {
    empty_mat <- matrix(NA, nrow = length(t_fine), ncol = trial_type_n + 1)
  }
  
  output_frame <- data.frame(empty_mat)
  colnames(output_frame) <- c("time", trial_types)
  
  resp_fn   <- gen_hrf(res_t = dt)$hrf
    
  for (m in 1:trial_type_n) {
    stim_fine <- rep(0, length(t_fine))
    
    # filter out the stim of interest
    stim_frame_trial <- stim_frame[(stim_frame$trial_type == trial_types[m]),]
    
    for (n in 1:length(stim_frame_trial$onset)) {
      index_bool <- t_fine >= stim_frame_trial$onset[n] & 
                    t_fine < stim_frame_trial$end[n]
      index <- which(index_bool)
      
      # only use one point if an impulse
      if (stim_frame_trial$duration[n] == min_dur) index <- index[1]
      
      stim_fine[index] <- 1
    }
    
    stim_fine <- stats::convolve(stim_fine, rev(resp_fn), type = 'open')
    stim_fine <- stim_fine[1:length(t_fine)]
    
    if (normalise) stim_fine <- stim_fine / max(stim_fine)
    
    t_acq    <- seq(from = 0, by = TR, length.out = n_trans)
    stim_acq <- stats::approx(t_fine, stim_fine, t_acq, method='linear')$y
    
    # if (normalise) stim_acq <- stim_acq / max(stim_acq)
    
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

#' Generate regressors by convolving a specified response function with a
#' stimulus.
#' @param onset stimulus onset in seconds.
#' @param duration stimulus duration in seconds.
#' @param trial_type string label for the stimulus.
#' @param mrs_data mrs_data object for timing information.
#' @param resp_fn a data frame specifying the response function to be convolved.
#' @param match_tr match the output to the input mrs_data.
#' @param normalise normalise the response function to have a maximum value of 
#' one.
#' @return BOLD regressor data frame.
#' @export
gen_conv_reg <- function(onset, duration = NULL, trial_type = NULL,
                         mrs_data = NULL, resp_fn, match_tr = TRUE,
                         normalise = FALSE) {
  
  # create a dummy dataset if not specified
  if (is.null(mrs_data)) {
    seq_tr   <- 2
    N_scans  <- 800
    mrs_data <- sim_resonances()
    mrs_data <- set_tr(mrs_data, seq_tr)
    mrs_data <- set_Ntrans(mrs_data, N_scans)
    mrs_data <- rep_dyn(mrs_data, N_scans)
  }
  
  dt <- resp_fn[2, 1] - resp_fn[1, 1]
  
  if (is.null(duration)) duration <- rep(0, length(onset))
  
  # set the minimum duration to dt * 1.1
  min_dur <- dt * 1.1
  duration[duration < min_dur] <- min_dur
  
  if (is.na(tr(mrs_data)) | is.null(tr(mrs_data))) {
    stop("TR not set, use set_tr function to set the repetition time.")
  }
  
  if (is.na(Ntrans(mrs_data)) | is.null(Ntrans(mrs_data))) {
    stop("Number of transients not set, use set_Ntrans function to set the 
         number of transients.")
  }
  
  if (is.null(trial_type)) trial_type <- rep("stim_conv", length(onset))
  
  # check everything is the right length 
  input_lengths <- c(length(onset), length(duration), length(trial_type))
  if (length(unique(input_lengths)) != 1) {
    print(input_lengths)
    stop("Stim length input error.")
  }
  
  # make a time scale with dt seconds resolution for the duration of the scan
  # time
  n_trans   <- mrs_data$meta$NumberOfTransients
  TR        <- tr(mrs_data)
  n_dyns    <- Ndyns(mrs_data)
  t_fine    <- seq(from = 0, to = n_trans * TR, by = dt)
  end       <- onset + duration
  
  stim_frame   <- data.frame(onset, end, trial_type, duration)
  
  trial_types  <- unique(trial_type)
  trial_type_n <- length(trial_types)
  
  if (match_tr) {
    empty_mat <- matrix(NA, nrow = n_dyns, ncol = trial_type_n + 1)
  } else {
    empty_mat <- matrix(NA, nrow = length(t_fine), ncol = trial_type_n + 1)
  }
  
  output_frame <- data.frame(empty_mat)
  colnames(output_frame) <- c("time", trial_types)
  
  resp_fn <- resp_fn[,2] 
    
  for (m in 1:trial_type_n) {
    stim_fine <- rep(0, length(t_fine))
    
    # filter out the stim of interest
    stim_frame_trial <- stim_frame[(stim_frame$trial_type == trial_types[m]),]
    
    for (n in 1:length(stim_frame_trial$onset)) {
      index_bool <- t_fine >= stim_frame_trial$onset[n] & 
                    t_fine < stim_frame_trial$end[n]
      index <- which(index_bool)
      
      # only use one point if an impulse
      if (stim_frame_trial$duration[n] == min_dur) index <- index[1]
      stim_fine[index] <- 1
    }
    
    stim_fine <- stats::convolve(stim_fine, rev(resp_fn), type = 'open')
    stim_fine <- stim_fine[1:length(t_fine)]
    
    if (normalise) stim_fine <- stim_fine / max(stim_fine)
    
    t_acq    <- seq(from = 0, by = TR, length.out = n_trans)
    stim_acq <- stats::approx(t_fine, stim_fine, t_acq, method='linear')$y
    
    # if (normalise) stim_acq <- stim_acq / max(stim_acq)
    
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

#' Generate impulse regressors.
#' @param onset stimulus onset in seconds.
#' @param trial_type string label for the stimulus.
#' @param mrs_data mrs_data object for timing information.
#' @return impulse regressors data frame.
#' @export
gen_impulse_reg <- function(onset, trial_type = NULL, mrs_data = NULL) {
  
  if (is.null(mrs_data)) {
    seq_tr   <- 2
    N_scans  <- 800
    mrs_data <- sim_resonances()
    mrs_data <- set_tr(mrs_data, seq_tr)
    mrs_data <- set_Ntrans(mrs_data, N_scans)
    mrs_data <- rep_dyn(mrs_data, N_scans)
  }
  
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

# gen double gamma model of hrf (as used in SPM) with 10ms resolution
# https://github.com/spm/spm12/blob/main/spm_hrf.m
gen_hrf <- function(end_t = 30, res_t = 0.01) {
  t_hrf <- seq(from = 0, to = end_t, by = res_t)
  a1 <- 6; a2 <- 16; b1 <- 1; b2 <- 1; c <- 1 / 6
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
  names <- colnames(regressor_df)[-1]
  X     <- t(regressor_df[, -1])
  graphics::image(y = time, z = X, col = viridisLite::viridis(128),
                  ylab = "Time (s)", axes = FALSE)
  graphics::axis(1, at=seq(0, 1, length = length(names)), labels = names)
  graphics::axis(2)
  graphics::box()
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
  colnames(reg_df) <- gsub("\\.", "_", colnames(reg_df)) # dots in names = bad
  return(reg_df)
}

#' Create a BIDS directory and file structure from a list of mrs_data objects.
#' @param mrs_data_list list of mrs_data objects.
#' @param output_dir the base directory to create the BIDS structure.
#' @param runs number of runs per subject and session.
#' @param sessions number of sessions.
#' @param sub_labels optional labels for subject level identification.
#' @export
mrs_data_list2bids <- function(mrs_data_list, output_dir, runs = 1,
                               sessions = 1, sub_labels = NULL) {
  
  Nmrs <- length(mrs_data_list)
  
  # check the number of datasets can be cleanly divided by the number of runs
  # and sessions
  if ((Nmrs %% (runs * sessions)) != 0) stop("inconsistent number of datasets")
  
  Nsubs <- Nmrs / runs / sessions
  
  if (is.null(sub_labels)) sub_labels <- auto_pad_seq(1:Nsubs)
  
  if (Nsubs != length(sub_labels)) stop("inconsistent subject labels")
  
  sub_labels <- paste0("sub-", sub_labels)
  sub_labels <- rep(sub_labels, each = runs * sessions)
  
  if (sessions != 1) {
    ses_labels <- auto_pad_seq(1:sessions)
    ses_labels <- paste0("ses-", ses_labels)
    ses_labels <- rep(ses_labels, each = runs)
    ses_labels <- rep(ses_labels, Nsubs * sessions)
  }
  
  if (runs != 1) {
    run_labels <- auto_pad_seq(1:runs)
    run_labels <- paste0("run-", run_labels)
    run_labels <- rep(run_labels, runs * sessions * Nsubs)
  }
  
  # construct the directories
  if (sessions == 1) {
    dirs <- file.path(output_dir, sub_labels, "mrs")
  } else {
    dirs <- file.path(output_dir, sub_labels, ses_labels, "mrs")
  }
  
  # create the directories
  for (dir in dirs) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  # construct the filenames
  fnames <- paste0(sub_labels)
  
  if (sessions != 1) fnames <- paste0(fnames, "_", ses_labels)
  
  if (runs != 1) fnames <- paste0(fnames, "_", run_labels)
  
  fnames <- paste0(fnames, "_svs.nii.gz")
  
  # construct the full path
  full_path <- file.path(dirs, fnames)
  
  # write the data
  for (n in 1:Nmrs) {
    write_mrs(mrs_data_list[[n]], fname = full_path[n], format = "nifti",
              force = TRUE) 
  }
}

auto_pad_seq <- function(x, min_pad = 2) {
  x_int <- sprintf("%d", x)
  pad_n <- max(nchar(x_int))
  if (pad_n < min_pad) pad_n <- min_pad
  fmt_string <- paste0("%0", pad_n, "d")
  return(sprintf(fmt_string, x))
}

#' Search for MRS data files in a BIDS filesystem structure.
#' @param path path to the directory containing the BIDS structure.
#' @return data frame containing full paths and information on each MRS file.
#' @export
find_bids_mrs <- function(path) {
  
  # find the "mrs" directories
  mrs_dirs <- dir(path, recursive = TRUE, include.dirs = TRUE, pattern = "mrs",
                  full.names = TRUE)
  
  # list all files in "mrs" directories
  mrs_paths <- list.files(mrs_dirs, full.names = TRUE)
  
  # remove any .json files
  mrs_paths <- grep(".json$", mrs_paths, invert = TRUE, value = TRUE)
  
  mrs_names <- basename(mrs_paths)
  mrs_names <- tools::file_path_sans_ext(tools::file_path_sans_ext(mrs_names))
  
  tags    <- strsplit(mrs_names, "_")
  tags_ul <- unlist(tags)
  
  sub <- grep("sub-", tags_ul, value = TRUE)
  sub <- substring(sub, 5)
  
  mrs_info <- data.frame(path = mrs_paths, sub = as.factor(sub))
  
  ses <- grep("ses-", tags_ul, value = TRUE)
  if (length(ses) != 0) {
    ses <- substring(ses, 5)
    mrs_info <- cbind(mrs_info, ses = as.factor(ses))
  }
  
  run <- grep("run-", tags_ul, value = TRUE)
  if (length(run) != 0) {
    run <- substring(run, 5)
    mrs_info <- cbind(mrs_info, run = as.factor(run))
  }
  
  suffix <- grep("-", tags_ul, value = TRUE, invert = TRUE)
  
  mrs_info <- cbind(mrs_info, suffix = as.factor(suffix))
 
  return(mrs_info) 
}

#' Preprocess and perform quality assessment of a single SVS data set.
#' @param path path to the fMRS data file or IMA directory.
#' @param label a label to describe the data set.
#' @param output_dir output directory.
#' @export
preproc_svs <- function(path, label = NULL, output_dir = NULL) {
  
  # TODO combine coils if needed, make the noise region a parameter
  # TODO deal with GE style data with wref included in the same file
  
  if (dir.exists(path)) {
    mrs_data <- read_ima_dyn_dir(path)
  } else {
    mrs_data <- read_mrs(path)
  }
  
  if (is.null(output_dir)) {
    output_dir <- getwd()
  } else {
    if (!dir.exists(output_dir)) dir.create(output_dir)
  }
  
  if (is.null(label)) {
    label <- basename(path)
    label <- tools::file_path_sans_ext(label)
    label <- tools::file_path_sans_ext(label)
  }
  
  mrs_rats  <- rats(mrs_data, xlim = c(4, 1.9), zero_freq_shift_t0 = TRUE,
                    ret_corr_only = FALSE)
  
  # perform simple baseline offset corrected based on a noisy spectral region
  mrs_rats$corrected <- bc_constant(mrs_rats$corrected, xlim = c(-0.5, -2.5))
  
  mean_mrs  <- mean_dyns(mrs_rats$corrected)
    
  # frequency and phase correct the mean spectrum
  ref <- sim_resonances(acq_paras = mrs_data, freq = c(2.01, 3.03, 3.22),
                        amp = 1, lw = 4, lg = 0)
  
  res <- rats(mean_mrs, ref, xlim = c(4, 1.9), p_deg = 4, ret_corr_only = FALSE)
  
  # apply mean spectrum phase and shift to the single shots
  mrs_proc <- phase(mrs_rats$corrected, -as.numeric(res$phases))
  mrs_proc <- shift(mrs_proc, -as.numeric(res$shifts), units = "hz")
  
  mrs_uncorr <-  phase(mrs_data,   -as.numeric(res$phases))
  mrs_uncorr <-  shift(mrs_uncorr,
                       -as.numeric(res$shifts) - mean(mrs_rats$shifts),
                       units = "hz")
  
  mean_uncorr <- mean_dyns(mrs_uncorr)
  
  snr         <- as.numeric(calc_spec_snr(mrs_proc))
  
  # single shot SNR
  ss_median_snr <- stats::median(snr)
  
  dyn_peak_info <- peak_info(mrs_proc, xlim = c(1.8, 2.2))
  lw_ppm        <- as.numeric(dyn_peak_info$fwhm_ppm)
  tnaa_height   <- as.numeric(dyn_peak_info$height)
  
  if (anyNA(lw_ppm)) {
    lw_ppm_smo <- rep(NA, length(lw_ppm))
  } else {
    lw_ppm_smo <- stats::smooth.spline(lw_ppm, spar = 0.8)$y
  }
  
  diag_table <- data.frame(dynamics = 1:length(snr),
                           time_sec = dyn_acq_times(mrs_data),
                           shifts_hz = as.numeric(mrs_rats$shifts),
                           phases = as.numeric(mrs_rats$phases),
                           snr = snr, lw_ppm = lw_ppm,
                           lw_ppm_smo = lw_ppm_smo, tnaa_height = tnaa_height)
  
  # scale data to the tCr peak
  amp <- spec_op(zf(res$corrected), xlim = c(2.9, 3.1), operator = "max-min")
  amp <- as.numeric(amp)
  mrs_proc      <- scale_mrs_amp(mrs_proc, 1 / amp)
  mrs_uncorr    <- scale_mrs_amp(mrs_uncorr, 1 / amp)
  res$corrected <- scale_mrs_amp(res$corrected, 1 / amp)
  mean_uncorr   <- scale_mrs_amp(mean_uncorr, 1 / amp)
  
  # mean spec SNR and LW
  mean_corr_spec_snr   <- calc_spec_snr(res$corrected)
  corr_peak_info       <- peak_info(res$corrected, xlim = c(1.8, 2.2))
  mean_corr_spec_lw    <- corr_peak_info$fwhm_ppm 
  mean_uncorr_spec_snr <- calc_spec_snr(mean_uncorr)
  uncorr_peak_info     <- peak_info(bc_constant(mean_uncorr,
                                                xlim = c(-0.5, -2.5)),
                                    xlim = c(1.8, 2.2))
  mean_uncorr_spec_lw  <- uncorr_peak_info$fwhm_ppm
  
  # measure the dynamic fluctuation range
  mrs_proc_smoothed <- smooth_dyns(crop_spec(lb(mrs_proc, 5)), 10)
  mrs_mean_sub      <- sub_mean_dyns(mrs_proc_smoothed)
  mrs_mean_sub_bc   <- bc_poly(mrs_mean_sub, 2)
  dfr               <- diff(range(Re(mrs_data2mat(mrs_mean_sub_bc))))
  
  # measure the dynamic lipid fluctuation range
  mrs_proc_smoothed_lip <- crop_spec(mrs_proc_smoothed, xlim = c(1.8, 0.4))
  mrs_mean_sub          <- sub_mean_dyns(mrs_proc_smoothed)
  mrs_mean_sub_bc       <- bc_poly(mrs_mean_sub, 2)
  dlfr                  <- diff(range(Re(mrs_data2mat(mrs_mean_sub_bc))))
  
  summary_diags <- c(mean_corr_spec_snr = mean_corr_spec_snr,
                     mean_corr_spec_lw = mean_corr_spec_lw,
                     mean_uncorr_spec_snr = mean_uncorr_spec_snr,
                     mean_uncorr_spec_lw = mean_uncorr_spec_lw,
                     ss_median_spec_snr = ss_median_snr,
                     lw_ppm_smo_range = diff(range(lw_ppm_smo)),
                     shift_hz_range = diff(range(diag_table$shifts_hz)),
                     dfr = dfr, dlfr = dlfr)
  
  res <- list(corrected = mrs_proc, uncorrected = mrs_uncorr,
              mean_corr = res$corrected, mean_uncorr = mean_uncorr,
              diag_table = diag_table, summary_diags = summary_diags,
              mrs_mean_sub = mrs_mean_sub, mrs_mean_sub_bc = mrs_mean_sub_bc)
  
  rmd_file <- system.file("rmd", "single_scan_svs_qa.Rmd",
                          package = "spant")
  
  rmd_out_f <- file.path(tools::file_path_as_absolute(output_dir), label)
  
  rmarkdown::render(rmd_file, params = list(data = res, label = label),
                    output_file = rmd_out_f)
  
  return(res)
}

#' Preprocess and perform quality assessment of one or more SVS data sets.
#' @param paths paths to the fMRS data file or IMA directory.
#' @param labels labels to describe each data set.
#' @param output_dir output directory.
#' @param exclude_labels vector of labels of scans to exclude, eg poor quality
#' data.
#' @export
preproc_svs_dataset <- function(paths, labels = NULL,
                                output_dir = "spant_analysis",
                                exclude_labels = NULL) {
  
  # TODO print warning if there are more preprocessed results than input
  # paths
  
  if (is.null(labels)) {
    labels <- basename(paths)
    labels <- tools::file_path_sans_ext(labels)
    labels <- tools::file_path_sans_ext(labels)
  }
  
  # lock in factor level order for ggplot output
  labels <- factor(labels, levels = labels)
  
  # detect non unique labels and quit if found
  if (any(table(labels) > 1)) stop("Labels are non-unique.")
  
  # find the indices of any excluded scans
  if (!is.null(exclude_labels)) {
    if (any(table(exclude_labels) > 1)) stop("Exclude labels are non-unique.")
    N_excl <- length(exclude_labels)
    exclude_inds <- rep(NA, N_excl)
    for (n in 1:N_excl) {
      found_bool <- (labels == exclude_labels[n])
      if (sum(found_bool) == 1) exclude_inds[n] <- which(found_bool)
    }
    if (anyNA(exclude_inds)) {
      print(exclude_labels[is.na(exclude_inds)])
      stop("Above exclude labels not found.")
    }
  }
  
  # root directory for all analysis results
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # directory for qa html reports 
  qa_dir <- file.path(output_dir, "qa")
  if (!dir.exists(qa_dir)) dir.create(qa_dir)
  
  # directory for pre-processing results
  preproc_dir <- file.path(output_dir, "preproc")
  if (!dir.exists(preproc_dir)) dir.create(preproc_dir)

  tot_num <- length(paths)
  preproc_res_list <- vector(mode = "list", length = tot_num)
  for (n in 1:tot_num) {
    cat(c("Processing ", n, " of ", tot_num, " : ",
          as.character(labels[n]), "\n"), sep = "")
    
    preproc_rds <- file.path(output_dir, "preproc", paste0(labels[n], ".rds"))
    
    if (file.exists(preproc_rds)) {
      cat("Preprocessing already performed, using saved results.\n")
      preproc_res_list[[n]] <- readRDS(preproc_rds)
      next
    }
    
    preproc_res_list[[n]] <- preproc_svs(paths[n], labels[n],
                                         file.path(output_dir, "qa"))
    
    saveRDS(preproc_res_list[[n]], preproc_rds)
  }
  
  preproc_summary <- data.frame(t(sapply(preproc_res_list,
                                         (\(x) x$summary_diags))))
  
  preproc_summary <- cbind(labels, preproc_summary)
  
  corrected_list <- lapply(preproc_res_list, \(x) x$corrected)
  
  if (tot_num > 1) {
    mean_dataset <- mean_mrs_list(corrected_list)
  } else {
    mean_dataset <- corrected_list[[1]]
  }
  
  res <- list(res_list = preproc_res_list, summary = preproc_summary,
              mean_dataset = mean_dataset, exclude_labels = NULL)
  
  rmd_file <- system.file("rmd", "dataset_summary_svs_qa.Rmd",
                          package = "spant")
  
  rmd_out_f <- file.path(tools::file_path_as_absolute(output_dir), "qa",
                         "dataset_summary_full")
  
  rmarkdown::render(rmd_file, params = list(data = res),
                    output_file = rmd_out_f)
  
  if (!is.null(exclude_labels)) {
    # exclude unwanted scans
    preproc_res_list <- preproc_res_list[-exclude_inds]
    
    preproc_summary <- data.frame(t(sapply(preproc_res_list,
                                           (\(x) x$summary_diags))))
    
    preproc_summary <- cbind(labels = labels[-exclude_inds], preproc_summary)
    
    corrected_list <- lapply(preproc_res_list, \(x) x$corrected)
    
    if (tot_num > 1) {
      mean_dataset <- mean_mrs_list(corrected_list)
    } else {
      mean_dataset <- corrected_list
    }
    
    res <- list(res_list = preproc_res_list, summary = preproc_summary,
                mean_dataset = mean_dataset, exclude_labels = exclude_labels)
    
    rmd_file <- system.file("rmd", "dataset_summary_svs_qa.Rmd",
                            package = "spant")
    
    rmd_out_f <- file.path(tools::file_path_as_absolute(output_dir), "qa",
                           "dataset_summary_subset")
    
    rmarkdown::render(rmd_file, params = list(data = res),
                      output_file = rmd_out_f)
  }
  
  return(mean_dataset)
  
  # saveRDS(res, file.path(output_dir, "preproc_full.rds"))
  
  #cut_res <- list(corrected = corrected_list, mean_dataset = mean_dataset,
  #                labels = labels)
  
  #saveRDS(cut_res, file.path(output_dir, "preproc_corrected.rds"))
}

#' Perform spectral GLM analysis of an fMRS dataset.
#' @param regressor_df a data frame containing temporal regressors to be applied
#' to each spectral datapoint.
#' @param analysis_dir directory containing preprocessed data generated by
#' the preproc_svs_dataset function.
#' @param exclude_labels vector of labels of scans to exclude, eg poor quality
#' @param labels labels to describe each data set.
#' @export
fmrs_dataset_spec_glm <- function(regressor_df,
                                  analysis_dir = "spant_analysis",
                                  exclude_labels = NULL, labels = NULL) {
  
  # TODO add optional arguments for datasets to preserve original
  # ordering if needed
  
  # check preproc_dir exists
  if (!dir.exists(analysis_dir)) stop("analysis_dir not found.")
  
  if (is.null(labels)) {
    # discover data labels from file names
    labels <- tools::file_path_sans_ext(basename(dir(file.path(analysis_dir,
                                                                   "preproc"))))
  }
  
  # detect non unique labels and quit if found
  if (any(table(labels) > 1)) stop("Labels are non-unique.")
  
  # remove any excluded labels
  if (!is.null(exclude_labels)) {
    if (any(table(exclude_labels) > 1)) stop("Exclude labels are non-unique.")
    N_excl <- length(exclude_labels)
    exclude_inds <- rep(NA, N_excl)
    for (n in 1:N_excl) {
      found_bool <- (labels == exclude_labels[n])
      if (sum(found_bool) == 1) exclude_inds[n] <- which(found_bool)
    }
    if (anyNA(exclude_inds)) {
      print(exclude_labels[is.na(exclude_inds)])
      stop("Above exclude labels not found.")
    }
    labels <- labels[-exclude_inds]
  }
  
  # directory for spec glm html reports 
  spec_glm_dir <- file.path(analysis_dir, "spec_glm")
  if (!dir.exists(spec_glm_dir)) dir.create(spec_glm_dir)
  
  # read preprocessed results
  Nscans <- length(labels)
  preproc_res_list <- vector(mode = "list", length = Nscans)
  for (n in 1:Nscans) {
    cat(c("Reading ", n, " of ", Nscans, " : ",
          as.character(labels[n]), "\n"), sep = "")
    preproc_path <- file.path(analysis_dir, "preproc",
                              paste0(labels[n], ".rds"))
    
    if (!file.exists(preproc_path)) stop("Preprocessing result not found.")
    
    preproc_res_list[[n]] <- readRDS(preproc_path) 
  }
  
  # extract just the corrected data
  corrected_list <- lapply(preproc_res_list, \(x) x$corrected)
  
  # calculate the mean spectrum
  if (Nscans > 1) {
    mean_dataset <- mean_mrs_list(corrected_list)
  } else {
    mean_dataset <- corrected_list[[1]]
  }
  
  # process the data
  mean_spec_glm <- lb(mean_dataset, 3)
  mean_spec_glm <- zf(mean_spec_glm)
  mean_spec_glm <- crop_spec(mean_spec_glm)
  mean_spec_glm <- bc_poly(mean_spec_glm, 1)
  
  glm_spec_res <- glm_spec(mean_spec_glm, regressor_df)
  
  # write output
  rmd_file <- system.file("rmd", "spec_glm_results.Rmd", package = "spant")
  
  rmd_out_f <- file.path(tools::file_path_as_absolute(analysis_dir), "spec_glm",
                         "dataset_mean")
  
  rmarkdown::render(rmd_file, params = list(data = glm_spec_res, 
                                            regressor_df = regressor_df,
                                            label = "dataset mean"),
                    output_file = rmd_out_f)
  
  return(glm_spec_res)
}
