#' Generate trapezoidal regressors.
#' @param onset stimulus onset in seconds.
#' @param duration stimulus duration in seconds.
#' @param trial_type string label for the stimulus.
#' @param mrs_data mrs_data object for timing information.
#' @param tr repetition time.
#' @param Ndyns number of dynamic scans stored, potentially less than Ntrans
#' if block averaging has been performed.
#' @param Ntrans number of dynamic scans acquired.
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
                         tr = NULL, Ndyns = NULL, Ntrans = NULL, rise_t = 0,
                         fall_t = 0, exp_fall = FALSE, exp_fall_power = 1,
                         smo_sigma = NULL, match_tr = TRUE, dt = 0.01,
                         normalise = FALSE) {
  
  res <- check_dyn_input(mrs_data, tr, Ndyns, Ntrans)  
  
  if (is.null(trial_type)) trial_type <- rep("stim", length(onset))
  
  # check everything is the right length 
  input_lengths <- c(length(onset), length(duration), length(trial_type))
  if (length(unique(input_lengths)) != 1) {
    print(input_lengths)
    stop("Stim length input error.")
  }
  
  # make a time scale with dt seconds resolution for the duration of the scan
  # time
  n_trans <- res$Ntrans
  TR      <- res$tr
  n_dyns  <- res$Ndyns
  t_fine  <- seq(from = 0, to = n_trans * TR - TR, by = dt)
  end     <- onset + duration
  
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
      stim_seg <- t_fine > stim_frame_trial$onset[n] & t_fine <= stim_frame_trial$end[n]
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
#' @param tr repetition time.
#' @param Ndyns number of dynamic scans stored, potentially less than Ntrans
#' if block averaging has been performed.
#' @param Ntrans number of dynamic scans acquired.
#' @param match_tr match the output to the input mrs_data.
#' @param dt timing resolution for internal calculations.
#' @param normalise normalise the response function to have a maximum value of 
#' one.
#' @return BOLD regressor data frame.
#' @export
gen_bold_reg <- function(onset, duration = NULL, trial_type = NULL,
                         mrs_data = NULL, tr = NULL, Ndyns = NULL,
                         Ntrans = NULL, match_tr = TRUE, dt = 0.1,
                         normalise = FALSE) {
  
  res <- check_dyn_input(mrs_data, tr, Ndyns, Ntrans)  
  
  if (is.null(duration)) duration <- rep(dt, length(onset))
  
  # set the minimum duration to dt * 1.1
  min_dur <- dt * 1.1
  duration[duration < min_dur] <- min_dur
  
  if (is.null(trial_type)) trial_type <- rep("stim_bold", length(onset))
  
  # check everything is the right length 
  input_lengths <- c(length(onset), length(duration), length(trial_type))
  if (length(unique(input_lengths)) != 1) stop("Stim length input error.")
  
  # make a time scale with dt seconds resolution for the duration of the scan
  # time
  n_trans <- res$Ntrans
  TR      <- res$tr
  n_dyns  <- res$Ndyns
  t_fine  <- seq(from = 0, to = n_trans * TR, by = dt)
  end     <- onset + duration
  
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
  
  resp_fn   <- get_hrf(res_t = dt)$resp_fn
    
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
#' @param tr repetition time.
#' @param Ndyns number of dynamic scans stored, potentially less than Ntrans
#' if block averaging has been performed.
#' @param Ntrans number of dynamic scans acquired.
#' @param resp_fn a data frame specifying the response function to be convolved.
#' @param match_tr match the output to the input mrs_data.
#' @param normalise normalise the response function to have a maximum value of 
#' one.
#' @return BOLD regressor data frame.
#' @export
gen_conv_reg <- function(onset, duration = NULL, trial_type = NULL,
                         mrs_data = NULL, tr = NULL, Ndyns = NULL,
                         Ntrans = NULL, resp_fn, match_tr = TRUE,
                         normalise = FALSE) {
  
  res <- check_dyn_input(mrs_data, tr, Ndyns, Ntrans)  
  
  dt <- resp_fn[2, 1] - resp_fn[1, 1]
  
  if (is.null(duration)) duration <- rep(0, length(onset))
  
  # set the minimum duration to dt * 1.1
  min_dur <- dt * 1.1
  duration[duration < min_dur] <- min_dur
  
  if (is.null(trial_type)) trial_type <- rep("stim_conv", length(onset))
  
  # check everything is the right length 
  input_lengths <- c(length(onset), length(duration), length(trial_type))
  if (length(unique(input_lengths)) != 1) {
    print(input_lengths)
    stop("Stim length input error.")
  }
  
  # make a time scale with dt seconds resolution for the duration of the scan
  # time
  n_trans <- res$Ntrans
  TR      <- res$tr
  n_dyns  <- res$Ndyns
  t_fine  <- seq(from = 0, to = n_trans * TR, by = dt)
  end     <- onset + duration
  
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
#' @param tr repetition time.
#' @param Ndyns number of dynamic scans stored, potentially less than Ntrans
#' if block averaging has been performed.
#' @param Ntrans number of dynamic scans acquired.
#' @return impulse regressors data frame.
#' @export
gen_impulse_reg <- function(onset, trial_type = NULL, mrs_data = NULL,
                            tr = NULL, Ndyns = NULL, Ntrans = NULL) {
  
  time <- dyn_acq_times(mrs_data, tr, Ndyns, Ntrans)
  
  if (is.null(trial_type)) trial_type <- rep("stim_imp", length(onset))
  
  trial_types  <- unique(trial_type)
  trial_type_n <- length(trial_types)
  
  stim_frame <- data.frame(onset, trial_type)
 
  n_dyns    <- length(time)
  empty_mat <- matrix(NA, nrow = n_dyns, ncol = trial_type_n)
  
  output_frame <- data.frame(empty_mat)
  colnames(output_frame) <- c(trial_types)
  
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

#' Generate baseline regressor.
#' @param mrs_data mrs_data object for timing information.
#' @param tr repetition time.
#' @param Ndyns number of dynamic scans stored, potentially less than Ntrans
#' if block averaging has been performed.
#' @param Ntrans number of dynamic scans acquired.
#' @return a single baseline regressor with value of 1.
#' @export
gen_baseline_reg <- function(mrs_data = NULL, tr = NULL, Ndyns = NULL,
                               Ntrans = NULL) {
    
  time   <- dyn_acq_times(mrs_data, tr, Ndyns, Ntrans)
  reg_df <- data.frame(time = time, baseline = rep(1, length(t)))
  return(reg_df)
}

#' Generate polynomial regressors.
#' @param degree the degree of the polynomial.
#' @param mrs_data mrs_data object for timing information.
#' @param tr repetition time.
#' @param Ndyns number of dynamic scans stored, potentially less than Ntrans
#' if block averaging has been performed.
#' @param Ntrans number of dynamic scans acquired.
#' @return polynomial regressors.
#' @export
gen_poly_reg <- function(degree, mrs_data = NULL, tr = NULL, Ndyns = NULL,
                         Ntrans = NULL) {
  
  time       <- dyn_acq_times(mrs_data, tr, Ndyns, Ntrans)
  poly_mat   <- stats::poly(time, degree)
  scale_vals <- apply(Mod(poly_mat), 2, max)
  poly_mat   <- scale(poly_mat, center = FALSE, scale = scale_vals)
  reg_df     <- data.frame(time = time, poly = poly_mat)
  colnames(reg_df) <- gsub("\\.", "_", colnames(reg_df)) # dots in names = bad
 
  if (degree == 1) colnames(reg_df) <- c("time", "poly_1")
  
  return(reg_df)
}

#' Generate a double gamma model of the HRF as used in SPM.
#' @param end_t last time point to generate in seconds.
#' @param res_t temporal resolution in seconds, defaults to 10ms.
#' @return a data.frame of time and HRF vectors.
#' @export
# https://github.com/spm/spm12/blob/main/spm_hrf.m
get_hrf <- function(end_t = 30, res_t = 0.01) {
  t_hrf <- seq(from = 0, to = end_t, by = res_t)
  a1 <- 6; a2 <- 16; b1 <- 1; b2 <- 1; c <- 1 / 6
  hrf <-     t_hrf ^ (a1 - 1) * b1 ^ a1 * exp(-b1 * t_hrf) / gamma(a1) -
         c * t_hrf ^ (a2 - 1) * b2 ^ a2 * exp(-b2 * t_hrf) / gamma(a2)
  hrf <- hrf / sum(hrf)
  return(data.frame(t = t_hrf, resp_fn = hrf))
}

#' Perform a GLM analysis of dynamic MRS data in the spectral domain.
#' @param mrs_data single-voxel dynamics MRS data.
#' @param regressor_df a data frame containing temporal regressors to be applied
#' to each spectral datapoint.
#' @param full_output append mrs_data and regressor_df to the output list.
#' @return list of statistical results.
#' @export
glm_spec <- function(mrs_data, regressor_df, full_output = FALSE) {
  
  # warning, any column named time in regressor_df will be removed
  
  # needs to be a FD operation
  if (!is_fd(mrs_data)) mrs_data <- td2fd(mrs_data)
  
  mrs_mat <- Re(mrs_data2mat(mrs_data))
  
  regressor_df_in <- regressor_df
  
  # drop the time column if present
  regressor_df<- regressor_df[, !names(regressor_df) %in% c("time"),
                              drop = FALSE]
  
  lm_res_list <- vector("list", ncol(mrs_mat))
  for (n in 1:ncol(mrs_mat)) {
    # no intercept
    lm_res_list[[n]] <- summary(stats::lm(mrs_mat[, n] ~ . + 0, regressor_df))
    
    # with intercept
    # lm_res_list[[n]] <- summary(stats::lm(mrs_mat[, n] ~ ., regressor_df))
  }
  
  # extract stats with intercept
  # get_glm_stat <- function(x, name) x$coefficients[-1, name, drop = FALSE]
  
  # extract stats no intercept
  get_glm_stat <- function(x, name) x$coefficients[, name, drop = FALSE]
  
  beta_weight <- as.data.frame(t(as.data.frame(sapply(lm_res_list, get_glm_stat,
                                               "Estimate", simplify = FALSE))))
  #beta_weight <- as.data.frame(t(as.data.frame(sapply(lm_res_list, get_glm_stat,
  #                                             "t value", simplify = FALSE))))
  p_value <- as.data.frame(t(as.data.frame(sapply(lm_res_list, get_glm_stat,
                                           "Pr(>|t|)", simplify = FALSE))))
  
  row.names(beta_weight) <- NULL
  row.names(p_value)     <- NULL
  
  ppm_sc      <-  ppm(mrs_data)
  beta_weight <-  cbind(ppm = ppm_sc, beta_weight)
  p_value_log <- -log10(p_value)
  p_value_log[p_value_log > 300] <- 300
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
  
  output_list <- list(beta_weight = beta_weight, p_value = p_value,
              p_value_log = p_value_log, p_value_log_mrs = p_value_log_mrs,
              p_value_mrs = p_value_mrs, beta_weight_mrs = beta_weight_mrs,
              lm_res_list = lm_res_list)
  
  if (full_output) output_list <- c(output_list, list(mrs_data = mrs_data,
                                               regressor_df = regressor_df_in))
  
  return(output_list)
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
  
  graphics::par(mar = c(2.8, 4.3, 1, 1))
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

#' Expand a regressor matrix for a group analysis.
#' @param regressor_df input regressor data frame.
#' @param n number of datasets n the group.
#' @export
gen_group_reg <- function(regressor_df, n) {
  # remove the time column and convert to a matrix
  tr <- regressor_df$time[2] - regressor_df$time[1]
  Nreg <- ncol(regressor_df) - 1
  offset <- utils::tail(regressor_df$time, 1) + tr
  reg_mat <- as.matrix(regressor_df[-1])
  reg_mat_group <- kronecker(diag(n), reg_mat)
  group_regressor_df <- as.data.frame(reg_mat_group)
  time <- rep(regressor_df$time, n) + 
          rep(offset * (0:(n - 1)), each = nrow(regressor_df))
  col_inds <- seq(from = 1, by = Nreg, length.out = n) + 
              rep(0:(Nreg - 1), each = n)
  
  group_regressor_df <- group_regressor_df[col_inds]
  
  colnames(group_regressor_df) <- paste0(rep(colnames(regressor_df)[-1], 
                                         each = n), "_scan_", 1:n)
  
  group_regressor_df <- cbind(time = time, group_regressor_df)
  return(group_regressor_df)
}

#' Create a BIDS file structure from a vector of MRS data paths or list of
#' mrs_data objects.
#' @param mrs_data vector of MRS data paths or list of mrs_data objects.
#' @param output_dir the base directory to create the BIDS structure.
#' @param suffix optional vector of file suffixes. Default behaviour is to
#' automatically determine these from the input data, however it is recommended
#' that they are specified to allow more efficient skipping of existing data.
#' @param sub optional vector of subject labels. If not specified, these will be
#' automatically generated as a series of increasing zero-padded integer values
#' corresponding to the mrs_data input indices.
#' @param ses optional vector of session labels.
#' @param task optional vector of task labels.
#' @param acq optional vector of acquisition labels.
#' @param nuc optional vector of nucleus labels.
#' @param voi optional vector of volume of interest labels.
#' @param rec optional vector of reconstruction labels.
#' @param run optional vector of run indices.
#' @param echo optional vector of echo time indices.
#' @param inv optional vector of inversion indices.
#' @param skip_existing skip any data files that have already been converted.
#' Defaults to TRUE, set to FALSE to force an overwrite of any existing data
#' files.
#' @export
mrs_data2bids <- function(mrs_data, output_dir, suffix = NULL, sub = NULL,
                          ses = NULL, task = NULL, acq = NULL, nuc = NULL,
                          voi = NULL, rec = NULL, run = NULL, echo = NULL,
                          inv = NULL, skip_existing = TRUE) {
  
  warning("mrs_data2bids is deprecated, use mr_data2bids instead")
  
  Nscans <- length(mrs_data)
  
  if (!is.null(suffix)) {
    
    if (length(suffix) != Nscans) {
      if (length(suffix) == 1) {
        suffix <- rep(suffix, Nscans)
      } else {
        stop("suffix length does not match.")
      }
    } 
    
    allowed <- c("svs", "mrsi", "unloc", "mrsref")
    wrong <- !(suffix %in% allowed)
    if (any(wrong)) stop(paste0("one or more suffix labels are invalid,",
                          " must be one of 'svs', 'mrsi', 'unloc', 'mrsref'"))
  }
  
  if (is.null(sub)) sub <- auto_pad_seq(1:Nscans)
  
  if (length(sub) != Nscans) stop("sub length does not match.")
  
  if (!is.null(ses)) {
    if (length(ses) == 1) ses <- rep(ses, Nscans)
    if (length(ses) != Nscans) stop("ses length does not match.")
  }
  
  if (!is.null(acq)) {
    if (length(acq) == 1) acq <- rep(acq, Nscans)
    if (length(acq) != Nscans) stop("acq length does not match.")
  }
  
  if (!is.null(nuc)) {
    if (length(nuc) == 1) nuc <- rep(nuc, Nscans)
    if (length(nuc) != Nscans) stop("nuc length does not match.")
  }
  
  if (!is.null(voi)) {
    if (length(voi) == 1) voi <- rep(voi, Nscans)
    if (length(voi) != Nscans) stop("voi length does not match.")
  }
  
  if (!is.null(rec)) {
    if (length(rec) == 1) rec <- rep(rec, Nscans)
    if (length(rec) != Nscans) stop("rec length does not match.")
  }
  
  if (!is.null(run)) {
    if (length(run) == 1) run <- rep(run, Nscans)
    if (length(run) != Nscans) stop("run length does not match.")
    run <- as.integer(run)
    if (any(is.na(run))) stop("non integer run")
  }
  
  if (!is.null(echo)) {
    if (length(echo) == 1) echo <- rep(echo, Nscans)
    if (length(echo) != Nscans) stop("echo length does not match.")
    echo <- as.integer(echo)
    if (any(is.na(echo))) stop("non integer echo")
  }
  
  if (!is.null(inv)) {
    if (length(inv) == 1) inv <- rep(inv, Nscans)
    if (length(inv) != Nscans) stop("inv length does not match.")
    inv <- as.integer(inv)
    if (any(is.na(inv))) stop("non integer inv")
  }
  
  for (n in 1:Nscans) {
    
    sub_lab <- paste0("sub-", sub[n])
    
    if (!is.null(ses)) ses_lab <- paste0("ses-", ses[n])
    
    # generate the directory structure
    if (is.null(ses)) {
      dir <- file.path(output_dir, sub_lab, "mrs")
    } else {
      dir <- file.path(output_dir, sub_lab, ses_lab, "mrs")
    }
    
    if (is.character(mrs_data[n])) {
      
      # skip if possible and we already know the suffix
      if (skip_existing & !is.null(suffix)) {
          
        # construct the filename
        fname <- paste0(sub_lab)
        if (!is.null(ses))  fname <- paste0(fname, "_", ses_lab)
        if (!is.null(task)) fname <- paste0(fname, "_task-", task[n])
        if (!is.null(acq))  fname <- paste0(fname, "_acq-", acq[n])
        if (!is.null(nuc))  fname <- paste0(fname, "_nuc-", nuc[n])
        if (!is.null(voi))  fname <- paste0(fname, "_voi-", voi[n])
        if (!is.null(rec))  fname <- paste0(fname, "_rec-", rec[n])
        if (!is.null(run))  fname <- paste0(fname, "_run-", run[n])
        if (!is.null(echo)) fname <- paste0(fname, "_echo-", echo[n])
        if (!is.null(inv))  fname <- paste0(fname, "_inv-", inv[n])
  
        # suffix 
        fname_main <- paste0(fname, "_", suffix[n], ".nii.gz")
  
        # construct the full path
        full_path_main <- file.path(dir, fname_main)
        
        if (file.exists(full_path_main)) {
          cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_main,
              "\n", sep = "")
          next
        }
      }
      
      mrs_n <- read_mrs(mrs_data[n])
    } else {
      mrs_n <- mrs_data[[n]]
    }
    
    if (identical(class(mrs_n), c("list", "mrs_data"))) {
      main    <- mrs_n$metab
      ref     <- mrs_n$ref
      ref_ecc <- mrs_n$ref_ecc
    } else {
      main    <- mrs_n
      ref     <- NULL
      ref_ecc <- NULL
    }
    
    if (is.null(suffix)) {
      voxels <- Nx(main) * Ny(main) * Nz(main)
      if (voxels == 1) {
        suffix_main <- "svs" 
      } else {
        suffix_main <- "mrsi" 
      }
    } else {
      suffix_main <- suffix[n] 
    }
  
    # create the directory
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    
    # construct the filename
    fname <- paste0(sub_lab)
    if (!is.null(ses))  fname <- paste0(fname, "_", ses_lab)
    if (!is.null(task)) fname <- paste0(fname, "_task-", task[n])
    if (!is.null(acq))  fname <- paste0(fname, "_acq-", acq[n])
    if (!is.null(nuc))  fname <- paste0(fname, "_nuc-", nuc[n])
    if (!is.null(voi))  fname <- paste0(fname, "_voi-", voi[n])
    if (!is.null(rec))  fname <- paste0(fname, "_rec-", rec[n])
    if (!is.null(run))  fname <- paste0(fname, "_run-", run[n])
    if (!is.null(echo)) fname <- paste0(fname, "_echo-", echo[n])
    if (!is.null(inv))  fname <- paste0(fname, "_inv-", inv[n])
  
    # suffix 
    fname_main <- paste0(fname, "_", suffix_main, ".nii.gz")
  
    # construct the full path
    full_path_main <- file.path(dir, fname_main)
    
    if (skip_existing & file.exists(full_path_main)) {
      cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_main, "\n",
          sep = "")
    } else {
      cat("Writing dataset ", n," of ", Nscans, " : ", full_path_main, "\n",
          sep = "")
      # write the data
      write_mrs(main, fname = full_path_main, format = "nifti", force = TRUE)
    }
    
    # just one ref file
    if (!is.null(ref) & is.null(ref_ecc)) {
      fname_ref <- paste0(fname, "_", "mrsref", ".nii.gz")
  
      # construct the full path
      full_path_ref <- file.path(dir, fname_ref)
      
      if (skip_existing & file.exists(full_path_ref)) {
        cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_ref, "\n",
            sep = "")
      } else {
        cat("Writing dataset ", n," of ", Nscans, " : ", full_path_ref, "\n",
            sep = "")
        write_mrs(ref, fname = full_path_ref, format = "nifti", force = TRUE)
      }
    }
    
    # both conc and ecc ref files
    if (!is.null(ref) & !is.null(ref_ecc)) {
      
      # conc filename
      if (is.null(acq)) {
        acq_lab <- paste0("_acq-conc")
      } else {
        acq_lab <- paste0("_acq-", acq[n], "conc")
      }
      fname <- paste0(sub_lab)
      if (!is.null(ses))  fname <- paste0(fname, "_", ses_lab)
      if (!is.null(task)) fname <- paste0(fname, "_task-", task[n])
      fname <- paste0(fname, acq_lab)
      if (!is.null(nuc))  fname <- paste0(fname, "_nuc-", nuc[n])
      if (!is.null(voi))  fname <- paste0(fname, "_voi-", voi[n])
      if (!is.null(rec))  fname <- paste0(fname, "_rec-", rec[n])
      if (!is.null(run))  fname <- paste0(fname, "_run-", run[n])
      if (!is.null(echo)) fname <- paste0(fname, "_echo-", echo[n])
      if (!is.null(inv))  fname <- paste0(fname, "_inv-", inv[n])
      fname_ref_conc <- paste0(fname, "_", "mrsref", ".nii.gz")
      full_path_conc <- file.path(dir, fname_ref_conc)
      
      if (skip_existing & file.exists(full_path_conc)) {
        cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_conc, "\n",
            sep = "")
      } else {
        cat("Writing dataset ", n," of ", Nscans, " : ", full_path_conc, "\n",
            sep = "")
        write_mrs(ref, fname = full_path_conc, format = "nifti", force = TRUE)
      }
      
      # ecc filename
      if (is.null(acq)) {
        acq_lab <- paste0("_acq-ecc")
      } else {
        acq_lab <- paste0("_acq-", acq[n], "ecc")
      }
      fname <- paste0(sub_lab)
      if (!is.null(ses))  fname <- paste0(fname, "_", ses_lab)
      if (!is.null(task)) fname <- paste0(fname, "_task-", task[n])
      fname <- paste0(fname, acq_lab)
      if (!is.null(nuc))  fname <- paste0(fname, "_nuc-", nuc[n])
      if (!is.null(voi))  fname <- paste0(fname, "_voi-", voi[n])
      if (!is.null(rec))  fname <- paste0(fname, "_rec-", rec[n])
      if (!is.null(run))  fname <- paste0(fname, "_run-", run[n])
      if (!is.null(echo)) fname <- paste0(fname, "_echo-", echo[n])
      if (!is.null(inv))  fname <- paste0(fname, "_inv-", inv[n])
      fname_ref_ecc <- paste0(fname, "_", "mrsref", ".nii.gz")
      full_path_ecc <- file.path(dir, fname_ref_ecc)
      
      if (skip_existing & file.exists(full_path_ecc)) {
        cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_ecc, "\n",
            sep = "")
      } else {
        cat("Writing dataset ", n," of ", Nscans, " : ", full_path_ecc, "\n",
            sep = "")
        write_mrs(ref_ecc, fname = full_path_ecc, format = "nifti",
                  force = TRUE)
      }
    }
  }
}

#' Create a BIDS file structure from a vector of data paths or list of
#' mri/mrs data objects.
#' @param mr_data vector of data paths or list of mri/mrs objects.
#' @param suffix vector of file suffixes, eg : c("svs", "mrsi", "T1w).
#' @param output_dir the base directory to create the BIDS structure.
#' @param sub optional vector of subject labels. If not specified, these will be
#' automatically generated as a series of increasing zero-padded integer values
#' corresponding to the mrs_data input indices.
#' @param ses optional vector of session labels.
#' @param task optional vector of task labels.
#' @param acq optional vector of acquisition labels.
#' @param nuc optional vector of nucleus labels.
#' @param voi optional vector of volume of interest labels.
#' @param rec optional vector of reconstruction labels.
#' @param run optional vector of run indices.
#' @param echo optional vector of echo time indices.
#' @param inv optional vector of inversion indices.
#' @param skip_existing skip any data files that have already been converted.
#' Defaults to TRUE, set to FALSE to force an overwrite of any existing data
#' files.
#' @param mri_format defaults to "nifti", can also be "dicom" provided the 
#' divest packages is installed.
#' @param deface_mri option to apply fsl_deface to the mri as a preprocessing
#' step. Defaults to FALSE, requires the fslr package to be installed when TRUE.
#' @export
mr_data2bids <- function(mr_data, suffix, output_dir, sub = NULL,
                         ses = NULL, task = NULL, acq = NULL, nuc = NULL,
                         voi = NULL, rec = NULL, run = NULL, echo = NULL,
                         inv = NULL, skip_existing = TRUE, 
                         mri_format = "nifti", deface_mri = FALSE) {
  
  if (!identical(class(mr_data), "list")) mr_data <- list(mr_data)
  
  Nscans <- length(mr_data)
  
  if (length(suffix) != Nscans) {
    if (length(suffix) == 1) {
      suffix <- rep(suffix, Nscans)
    } else {
      stop("suffix length does not match.")
    }
  } 
  
  mrs_suffix  <- c("svs", "mrsi", "unloc", "mrsref")
  func_suffix <- c("bold", "cbv")
  anat_suffix <- c("FLAIR", "PDT2", "PDw", "T1w", "T2starw", "T2w", "UNIT1",
                   "angio", "inplaneT1", "inplaneT2")
  
  allowed <- c(mrs_suffix, func_suffix, anat_suffix)
  
  wrong <- !(suffix %in% allowed)
  if (any(wrong)) stop(paste0("one or more suffix labels are invalid,",
                              " must be one of : ", 
                              paste(allowed, collapse = ", ")))
  
  if (is.null(sub)) sub <- auto_pad_seq(1:Nscans)
  
  if (length(sub) != Nscans) stop("sub length does not match.")
  
  if (!is.null(ses)) {
    if (length(ses) == 1) ses <- rep(ses, Nscans)
    if (length(ses) != Nscans) stop("ses length does not match.")
  }
  
  if (!is.null(acq)) {
    if (length(acq) == 1) acq <- rep(acq, Nscans)
    if (length(acq) != Nscans) stop("acq length does not match.")
  }
  
  if (!is.null(nuc)) {
    if (length(nuc) == 1) nuc <- rep(nuc, Nscans)
    if (length(nuc) != Nscans) stop("nuc length does not match.")
  }
  
  if (!is.null(voi)) {
    if (length(voi) == 1) voi <- rep(voi, Nscans)
    if (length(voi) != Nscans) stop("voi length does not match.")
  }
  
  if (!is.null(rec)) {
    if (length(rec) == 1) rec <- rep(rec, Nscans)
    if (length(rec) != Nscans) stop("rec length does not match.")
  }
  
  if (!is.null(run)) {
    if (length(run) == 1) run <- rep(run, Nscans)
    if (length(run) != Nscans) stop("run length does not match.")
    run <- as.integer(run)
    if (any(is.na(run))) stop("non integer run")
  }
  
  if (!is.null(echo)) {
    if (length(echo) == 1) echo <- rep(echo, Nscans)
    if (length(echo) != Nscans) stop("echo length does not match.")
    echo <- as.integer(echo)
    if (any(is.na(echo))) stop("non integer echo")
  }
  
  if (!is.null(inv)) {
    if (length(inv) == 1) inv <- rep(inv, Nscans)
    if (length(inv) != Nscans) stop("inv length does not match.")
    inv <- as.integer(inv)
    if (any(is.na(inv))) stop("non integer inv")
  }
  
  for (n in 1:Nscans) {
    
    sub_lab <- paste0("sub-", sub[n])
    
    if (!is.null(ses)) {
      ses_lab <- paste0("ses-", ses[n])
    } else {
      ses_lab <- NULL
    }
    
    # determine the image type
    
    if (suffix[n] %in% mrs_suffix) {
      image_type <- "mrs"
    } else if (suffix[n] %in% anat_suffix) {
      image_type <- "anat"
    } else if (suffix[n] %in% func_suffix) {
      image_type <- "func"
    } else {
      stop("image type not found for given suffix")
    }
    
    # generate the directory structure
    if (is.null(ses)) {
      dir <- file.path(output_dir, sub_lab, image_type)
    } else {
      dir <- file.path(output_dir, sub_lab, ses_lab, image_type)
    }
    
    if (identical(class(mr_data[[n]]), "character")) {
      
      # skip if possible and we already know the suffix
      if (skip_existing & !is.null(suffix)) {
        
        # construct the filename
        fname <- paste0(sub_lab)
        fname <- paste0(fname, build_bids_fname(ses_lab, task, acq, nuc, voi,
                                                rec, run, echo, inv, n))
        
        # suffix 
        fname_main <- paste0(fname, "_", suffix[n], ".nii.gz")
        
        # construct the full path
        full_path_main <- file.path(dir, fname_main)
        
        if (file.exists(full_path_main)) {
          cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_main,
              "\n", sep = "")
          next
        }
      }
      
      if (image_type == "mrs") {
        mr_n <- read_mrs(mr_data[[n]])
      } else {
        if (mri_format == "nifti") {
          mr_n <- readNifti(mr_data[[n]], json = "read")  
        } else if (mri_format == "dicom") {
          mr_n <- divest::readDicom(mr_data[[n]], interactive = FALSE,
                                    verbosity = -1)[[1]]
        } else {
          stop("Incorrect mri_format.")
        }
        
        if (deface_mri) {
          deface_fsl <- fslr::fsl_deface(mr_n, verbose = FALSE)
          deface     <- RNifti::asNifti(deface_fsl$outfile)
          # copy the image meta data from the original mri
          RNifti::imageAttributes(deface) <- RNifti::imageAttributes(mr_n)
          mr_n <- deface
        }
      }
       
    } else {
      mr_n <- mr_data[[n]]
    }
    
    if (identical(class(mr_n), c("list", "mrs_data"))) {
      main    <- mr_n$metab
      ref     <- mr_n$ref
      ref_ecc <- mr_n$ref_ecc
    } else if (identical(class(mr_n), c("mrs_data"))) {
      main    <- mr_n
      ref     <- NULL
      ref_ecc <- NULL
    } else {
      main    <- mr_n
      ref     <- NULL
      ref_ecc <- NULL
    }
    
    # create the directory
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    
    # construct the filename
    fname <- paste0(sub_lab)
    fname <- paste0(fname, build_bids_fname(ses_lab, task, acq, nuc, voi, rec, 
                                            run, echo, inv, n))
    
    # suffix 
    fname_main <- paste0(fname, "_", suffix[n], ".nii.gz")
    
    # construct the full path
    full_path_main <- file.path(dir, fname_main)
    
    if (skip_existing & file.exists(full_path_main)) {
      cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_main, "\n",
          sep = "")
    } else {
      cat("Writing dataset ", n," of ", Nscans, " : ", full_path_main, "\n",
          sep = "")
      # write the data
      if (image_type == "mrs") {
        write_mrs(main, fname = full_path_main, format = "nifti", force = TRUE)
      } else {
        writeNifti(main, file = full_path_main, json = TRUE)
      }
    }
    
    # the following only applies to MRS data
    
    # just one ref file
    if (!is.null(ref) & is.null(ref_ecc)) {
      fname_ref <- paste0(fname, "_", "mrsref", ".nii.gz")
      
      # construct the full path
      full_path_ref <- file.path(dir, fname_ref)
      
      if (skip_existing & file.exists(full_path_ref)) {
        cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_ref, "\n",
            sep = "")
      } else {
        cat("Writing dataset ", n," of ", Nscans, " : ", full_path_ref, "\n",
            sep = "")
        write_mrs(ref, fname = full_path_ref, format = "nifti", force = TRUE)  
      }
    }
        
    # both conc and ecc ref files
    if (!is.null(ref) & !is.null(ref_ecc)) {
      
      # conc filename
      if (is.null(acq)) {
        acq_lab <- paste0("conc")
      } else {
        acq_lab <- paste0(acq[n], "conc")
      }
      
      fname <- paste0(sub_lab)
      fname <- paste0(fname, build_bids_fname(ses_lab, task, acq_lab, nuc, voi,
                                              rec, run, echo, inv, n))                    
                          
      fname_ref_conc <- paste0(fname, "_", "mrsref", ".nii.gz")
      full_path_conc <- file.path(dir, fname_ref_conc)
      
      if (skip_existing & file.exists(full_path_conc)) {
        cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_conc, "\n",
            sep = "")
      } else {
        cat("Writing dataset ", n," of ", Nscans, " : ", full_path_conc, "\n",
            sep = "")
        write_mrs(ref, fname = full_path_conc, format = "nifti", force = TRUE)
      }
      
      # ecc filename
      if (is.null(acq)) {
        acq_lab <- paste0("ecc")
      } else {
        acq_lab <- paste0(acq[n], "ecc")
      }
      
      fname <- paste0(sub_lab)
      fname <- paste0(fname, build_bids_fname(ses_lab, task, acq_lab, nuc, voi,
                                              rec, run, echo, inv, n))               
                          
      fname_ref_ecc <- paste0(fname, "_", "mrsref", ".nii.gz")
      full_path_ecc <- file.path(dir, fname_ref_ecc)
      
      if (skip_existing & file.exists(full_path_ecc)) {
        cat("Skipping dataset ", n," of ", Nscans, " : ", full_path_ecc, "\n",
            sep = "")
      } else {
        cat("Writing dataset ", n," of ", Nscans, " : ", full_path_ecc, "\n",
            sep = "")
        write_mrs(ref_ecc, fname = full_path_ecc, format = "nifti",
                  force = TRUE)
      }
    }
  }
}

build_bids_fname <- function(ses, task, acq, nuc, voi, rec, run, echo, inv, n) {
  
  if (!is.null(ses))  ses  <- ifelse(length(ses)  == n, ses[n],  ses)
  if (!is.null(task)) task <- ifelse(length(task) == n, task[n], task)
  if (!is.null(acq))  acq  <- ifelse(length(acq)  == n, acq[n],  acq)
  if (!is.null(nuc))  nuc  <- ifelse(length(nuc)  == n, nuc[n],  nuc)
  if (!is.null(voi))  voi  <- ifelse(length(voi)  == n, voi[n],  voi)
  if (!is.null(rec))  rec  <- ifelse(length(rec)  == n, rec[n],  rec)
  if (!is.null(run))  run  <- ifelse(length(run)  == n, run[n],  run)
  if (!is.null(echo)) echo <- ifelse(length(echo) == n, echo[n], echo)
  if (!is.null(inv))  inv  <- ifelse(length(inv)  == n, inv[n],  inv)
  
  fname <- ""
  if (!is.null(ses))  fname <- paste0(fname, "_",      ses)
  if (!is.null(task)) fname <- paste0(fname, "_task-", task)
  if (!is.null(acq))  fname <- paste0(fname, "_acq-",  acq)
  if (!is.null(nuc))  fname <- paste0(fname, "_nuc-",  nuc)
  if (!is.null(voi))  fname <- paste0(fname, "_voi-",  voi)
  if (!is.null(rec))  fname <- paste0(fname, "_rec-",  rec)
  if (!is.null(run))  fname <- paste0(fname, "_run-",  run)
  if (!is.null(echo)) fname <- paste0(fname, "_echo-", echo)
  if (!is.null(inv))  fname <- paste0(fname, "_inv-",  inv)
  
  return(fname)
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
#' @param output_full_path output the full normalised data paths.
#' @return data frame containing full paths and information on each MRS file.
#' @export
find_bids_mrs <- function(path, output_full_path = FALSE) {
  
  # find the "mrs" directories
  mrs_dirs <- dir(path, recursive = TRUE, include.dirs = TRUE, pattern = "mrs$",
                  full.names = TRUE)
  
  # list all files in "mrs" directories
  mrs_paths <- list.files(mrs_dirs, full.names = TRUE)
  
  # remove any .json files
  mrs_paths <- grep(".json$", mrs_paths, invert = TRUE, value = TRUE)
  
  if (output_full_path) mrs_paths <- normalizePath(mrs_paths)
  
  mrs_names <- basename(mrs_paths)
  mrs_names <- tools::file_path_sans_ext(tools::file_path_sans_ext(mrs_names))
  
  Nscans <- length(mrs_paths)
  
  if (Nscans == 0) stop("No data files were found.")
  
  # create an empty data frame
  na_vec  <- rep(NA, Nscans) 
  bids_df <- data.frame(path = mrs_paths, sub = na_vec, ses = na_vec,
                        task = na_vec, acq = na_vec, nuc = na_vec,
                        voi = na_vec, rec = na_vec, run = na_vec,
                        echo = na_vec, inv = na_vec, suffix = na_vec)
  
  # for each scan
  for (n in 1:Nscans) {
    
    tags    <- strsplit(mrs_names[n], "_")
    tags_ul <- unlist(tags)
    
    sub <- grep("sub-", tags_ul, value = TRUE)
    bids_df$sub[n] <- substring(sub, 5)
    
    ses <- grep("ses-", tags_ul, value = TRUE)
    if (length(ses) != 0) bids_df$ses[n] <- substring(ses, 5)
    
    task <- grep("task-", tags_ul, value = TRUE)
    if (length(task) != 0) bids_df$task[n] <- substring(task, 6)
    
    acq <- grep("acq-", tags_ul, value = TRUE)
    if (length(acq) != 0) bids_df$acq[n] <- substring(acq, 5)
    
    nuc <- grep("nuc-", tags_ul, value = TRUE)
    if (length(nuc) != 0) bids_df$nuc[n] <- substring(nuc, 5)
    
    voi <- grep("voi-", tags_ul, value = TRUE)
    if (length(voi) != 0) bids_df$voi[n] <- substring(voi, 5)
    
    rec <- grep("rec-", tags_ul, value = TRUE)
    if (length(rec) != 0) bid_df$rec[n] <- substring(rec, 5)
    
    run <- grep("run-", tags_ul, value = TRUE)
    if (length(run) != 0) bids_df$run[n] <- substring(run, 5)
    
    echo <- grep("echo-", tags_ul, value = TRUE)
    if (length(echo) != 0) bids_df$echo[n] <- substring(echo, 6)
    
    inv <- grep("inv-", tags_ul, value = TRUE)
    if (length(inv) != 0) bids_df$inv[n] <- substring(inv, 5)
    
    suffix <- grep("-", tags_ul, value = TRUE, invert = TRUE)
    bids_df$suffix[n] <- suffix
  }
 
  return(bids_df) 
}

#' Preprocess and perform quality assessment of a single SVS data set.
#' @param path path to the fMRS data file or IMA directory.
#' @param label a label to describe the data set.
#' @param output_dir output directory.
#' @param ref_inds a vector of 1-based indices for any water reference dynamic
#' scans.
#' @export
preproc_svs <- function(path, label = NULL, output_dir = NULL,
                        ref_inds = NULL) {
  
  # TODO deal with GE style data with wref included in the same file
  
  if (dir.exists(path)) {
    mrs_data <- read_ima_dyn_dir(path)
  } else {
    mrs_data <- read_mrs(path)
  }
  
  if (!is.null(ref_inds)) {
    mrs_data <- get_dyns(mrs_data, -ref_inds)
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
  
  # combine coils if needed
  if (Ncoils(mrs_data) > 1) mrs_data <- comb_coils_svs_gls(mrs_data)
  
  # if the spectral width exceeds 20 PPM, then it should be ok
  # to decimate the signal to improve analysis speed
  decimate <- NULL
  if (is.null(decimate) & ((fs(mrs_data) / mrs_data$ft * 1e6) > 20)) {
    decimate <- TRUE
  } else {
    decimate <- FALSE
  }
  
  if (decimate) {
    metab <- decimate_mrs_fd(mrs_data)
    # if (!is.null(w_ref)) w_ref <- decimate_mrs_fd(w_ref)
  }
  
  mrs_rats  <- rats(mrs_data, xlim = c(4, 1.9), zero_freq_shift_t0 = TRUE,
                    ret_corr_only = FALSE)
  
  # perform simple baseline offset corrected based on a noisy spectral region
  mrs_rats$corrected <- bc_constant(mrs_rats$corrected, xlim = c(-0.5, -2.5))
  
  mean_mrs  <- mean_dyns(mrs_rats$corrected)
    
  # frequency and phase correct the mean spectrum
  res <- phase_ref_1h_brain(mean_mrs, ret_corr_only = FALSE)
  
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
#' @param overwrite overwrite saved results, defaults to FALSE.
#' @param ref_inds a vector of 1-based indices for any water reference dynamic
#' scans.
#' @param return_results function will return key outputs, defaults to FALSE.
#' @export
preproc_svs_dataset <- function(paths, labels = NULL,
                                output_dir = "spant_analysis",
                                exclude_labels = NULL, overwrite = FALSE,
                                ref_inds = NULL, return_results = FALSE) {
  
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
  rds_dir <- file.path(output_dir, "preproc", "rds")
  if (!dir.exists(rds_dir)) dir.create(rds_dir)
  metab_dir <- file.path(output_dir, "preproc", "metab_dyn")
  if (!dir.exists(metab_dir)) dir.create(metab_dir)
  metab_dir <- file.path(output_dir, "preproc", "metab_av_dyn")
  if (!dir.exists(metab_dir)) dir.create(metab_dir)

  tot_num <- length(paths)
  preproc_res_list <- vector(mode = "list", length = tot_num)
  for (n in 1:tot_num) {
    cat(c("Processing ", n, " of ", tot_num, " : ",
          as.character(labels[n]), "\n"), sep = "")
    
    preproc_rds <- file.path(output_dir, "preproc", "rds",
                             paste0(labels[n], ".rds"))
    
    if (file.exists(preproc_rds)) {
      if (overwrite) {
        cat("Overwriting saved results.\n")
      }  else {
        cat("Preprocessing already performed, using saved results.\n")
        preproc_res_list[[n]] <- readRDS(preproc_rds)
        next
      }
    }
    
    preproc_res_list[[n]] <- preproc_svs(paths[n], labels[n],
                                         file.path(output_dir, "qa"), ref_inds)
    
    saveRDS(preproc_res_list[[n]], preproc_rds)
    
    # write dynamic MRS data if available
    if (Ndyns(preproc_res_list[[n]]$corrected) > 1) {
      preproc_metab <- file.path(output_dir, "preproc", "metab_dyn",
                                 paste0(labels[n], ".nii.gz"))
    
      write_mrs(preproc_res_list[[n]]$corrected, preproc_metab, force = TRUE)
      
      # write temporally averaged data 
      preproc_metab <- file.path(output_dir, "preproc", "metab_av_dyn",
                                 paste0(labels[n], ".nii.gz"))
    
      write_mrs(mean_dyns(preproc_res_list[[n]]$corrected), preproc_metab,
                force = TRUE)
      
  } else {
      # write temporally averaged data 
      preproc_metab <- file.path(output_dir, "preproc", "metab_av_dyn",
                                 paste0(labels[n], ".nii.gz"))
    
      write_mrs(preproc_res_list[[n]]$corrected, preproc_metab, force = TRUE)
    }
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
  
  csv_out_f <- file.path(tools::file_path_as_absolute(output_dir), "qa",
                         "qa_summary_full.csv")
  
  utils::write.csv(preproc_summary, csv_out_f)
  
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
    
    csv_out_f <- file.path(tools::file_path_as_absolute(output_dir), "qa",
                           "qa_summary_subset.csv")
    
    utils::write.csv(preproc_summary, csv_out_f)
  }
  
  if (return_results) return(list(preproc_metab_list = corrected_list,
                             preproc_metab_mean = mean_dataset))
}

#' Perform first-level spectral GLM analysis of an fMRS dataset.
#' @param regressor_df a data frame containing temporal regressors to be applied
#' to each spectral datapoint.
#' @param analysis_dir directory containing preprocessed data generated by
#' the preproc_svs_dataset function.
#' @param exclude_labels vector of labels of scans to exclude, eg poor quality
#' data.
#' @param labels labels to describe each data set.
#' @param xlim spectral range to include in the analysis.
#' @param vline vertical lines to add to the plot.
#' @param lb linebroading to add in Hz before GLM analysis.
#' @param return_results function will return key outputs, defaults to FALSE.
#' @export
glm_spec_fmrs_fl <- function(regressor_df, analysis_dir = "spant_analysis",
                             exclude_labels = NULL, labels = NULL,
                             xlim = c(4, 0.2),
                             vline = c(1.35, 1.28, 2.35, 2.29),
                             lb = 4,
                             return_results = FALSE) {
  
  # TODO add optional arguments for datasets to preserve original
  # ordering if needed
  
  # check preproc_dir exists
  if (!dir.exists(analysis_dir)) stop("analysis_dir not found.")
  
  if (is.null(labels)) {
    # discover data labels from file names
    labels <- tools::file_path_sans_ext(basename(dir(file.path(analysis_dir,
                                                               "preproc",
                                                               "rds"))))
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
  spec_glm_dir <- file.path(analysis_dir, "spec_glm", "first_level")
  if (!dir.exists(spec_glm_dir)) dir.create(spec_glm_dir, recursive = TRUE)
  
  # directory for processed spectra
  spec_glm_mrs_dir <- file.path(analysis_dir, "spec_glm", "proc_fmrs")
  if (!dir.exists(spec_glm_mrs_dir)) dir.create(spec_glm_mrs_dir)
  
  # read preprocessed results
  Nscans <- length(labels)
  preproc_res_list <- vector(mode = "list", length = Nscans)
  for (n in 1:Nscans) {
    cat(c("Reading ", n, " of ", Nscans, " : ",
          as.character(labels[n]), "\n"), sep = "")
    preproc_path <- file.path(analysis_dir, "preproc", "rds",
                              paste0(labels[n], ".rds"))
    
    if (!file.exists(preproc_path)) stop("Preprocessing result not found.")
    
    preproc_res_list[[n]] <- readRDS(preproc_path) 
  }
  
  # run glm spec on each result separately
  glm_spec_res_list <- vector(mode = "list", length = Nscans)
  for (n in 1:Nscans) {
    mrs_data_glm <- preproc_res_list[[n]]$corrected
    glm_spec_res_list[[n]] <- gen_glm_spec_report(mrs_data_glm, regressor_df,
                                                  labels[n], analysis_dir, xlim,
                                                  vline, lb = lb)
    
    # write processed data 
    preproc_metab <- file.path(spec_glm_mrs_dir, paste0(labels[n], ".nii.gz"))
    
    write_mrs(glm_spec_res_list[[n]]$mrs_data, preproc_metab,
              force = TRUE)
    
  }
  
  # extract just the corrected data
  corrected_list <- lapply(preproc_res_list, \(x) x$corrected)
  
  # calculate the mean spectrum
  if (Nscans > 1) {
    mean_dataset <- mean_mrs_list(corrected_list)
  } else {
    mean_dataset <- corrected_list[[1]]
  }
  
  # run glm spec on the mean dataset
  glm_spec_mean_res <- gen_glm_spec_report(mean_dataset, regressor_df,
                          "dataset_mean", analysis_dir,
                           xlim, vline, exclude_labels, lb = lb)
  
  if (return_results) {
    return(list(indiv_res = glm_spec_res_list, mean_res = glm_spec_mean_res))
  }
}

gen_glm_spec_report <- function(mrs_data, regressor_df, label, analysis_dir,
                                xlim, vline, exclude_labels = NULL, lb) {
  
  # process the data
  mrs_data_glm <- lb(mrs_data, lb)
  mrs_data_glm <- zf(mrs_data_glm)
  mrs_data_glm <- crop_spec(mrs_data_glm, xlim = xlim)
  
  # mrs_data_glm <- bc_poly(mrs_data_glm, 0)
  # mrs_data_glm <- bc_als(mrs_data_glm, lambda = 10)
  
  mrs_data_plot <- zf(mean_dyns(mrs_data))
  
  # write the processed spectra
  
  
  # regress
  glm_spec_res <- glm_spec(mrs_data_glm, regressor_df, full_output = TRUE)
  
  # write output
  rmd_file <- system.file("rmd", "spec_glm_results.Rmd", package = "spant")
  
  rmd_out_f <- file.path(tools::file_path_as_absolute(analysis_dir), "spec_glm",
                         "first_level", label)
  
  rmarkdown::render(rmd_file, params = list(data = glm_spec_res, 
                    regressor_df = regressor_df, label = label,
                    mrs_data_plot = mrs_data_plot, xlim = xlim, vline = vline,
                    exclude_labels = exclude_labels),
                    output_file = rmd_out_f)
  
  return(glm_spec_res)
}

# Doesn't work...
# glm_spec_group_analysis <- function(glm_spec_dataset) {
#   Nreg <- ncol(glm_spec_dataset[[1]]$beta_weight) - 1
#   Npts <- nrow(glm_spec_dataset[[1]]$beta_weight)
#   out_colnames <- colnames(glm_spec_dataset[[1]]$beta_weight)
#   ppm_sc <- glm_spec_dataset[[1]]$beta_weight$ppm
#   mrs_data <- glm_spec_dataset[[1]]$beta_weight_mrs
#   
#   get_betas <- function(x, n) x$beta_weight[, n + 1]
#   get_p_val <- function(x) x$p.value
#   get_beta_mean <- function(x) x$estimate
#   
#   p_value   <- matrix(nrow = Npts, ncol = Nreg)
#   beta_mean <- matrix(nrow = Npts, ncol = Nreg)
#   for (n in 1:Nreg) {
#     betas <- lapply(glm_spec_dataset, get_betas, n = n)
#     betas <- do.call(rbind, betas)
#     stat_res <- apply(betas, 2, stats::t.test)
#     p_value[, n] <- sapply(stat_res, get_p_val)
#     beta_mean[, n] <- sapply(stat_res, get_beta_mean)
#   }
#   
#   beta_mean             <-  as.data.frame(beta_mean)
#   p_value               <-  as.data.frame(p_value)
#   p_value_log           <- -log10(p_value)
#   beta_mean             <-  cbind(ppm = ppm_sc, beta_mean)
#   p_value               <-  cbind(ppm = ppm_sc, p_value)
#   p_value_log           <-  cbind(ppm = ppm_sc, p_value_log)
#   colnames(beta_mean)   <-  out_colnames
#   colnames(p_value)     <-  out_colnames
#   colnames(p_value_log) <-  out_colnames
#   
#   beta_mean_mrs <- mat2mrs_data(t(beta_mean[, -1]), fs = fs(mrs_data),
#                                   ft = mrs_data$ft, ref = mrs_data$ref,
#                                   nuc = mrs_data$nuc, fd = TRUE)
#   
#   p_value_log_mrs <- mat2mrs_data(t(p_value_log[, -1]), fs = fs(mrs_data),
#                                   ft = mrs_data$ft, ref = mrs_data$ref,
#                                   nuc = mrs_data$nuc, fd = TRUE)
#   
#   p_value_mrs <- mat2mrs_data(t(p_value[, -1]), fs = fs(mrs_data),
#                               ft = mrs_data$ft, ref = mrs_data$ref,
#                               nuc = mrs_data$nuc, fd = TRUE)
#   
#   return(list(p_value_log_mrs = p_value_log_mrs,
#               p_value_mrs = p_value_mrs,
#               beta_weight_mrs = beta_mean_mrs,
#               p_value_log = p_value_log,
#               p_value = p_value,
#               beta_weight = beta_mean))
# }

#' Perform a t-test on spectral data points.
#' @param mrs_data an mrs_data object with spectra in the dynamic dimension.
#' @param group vector describing the group membership of each dynamic spectrum.
#' @return a list of statistical results.
#' @export
t_test_spec <- function(mrs_data, group) {
  
  if (missing(group)) stop("group argument is missing.")
  
  mrs_data_mat <- mrs_data2spec_mat(mrs_data)
  
  # catches errors when passing identical values which can crop up when
  # normalising to a maximum spectral data point
  t_test_calc <- function(x, group, y) {
    tryCatch(stats::t.test(x ~ group), error = function(e) list(p.value = 1,
                                                                statistic = 0))
  }

  t_test_res_list <- apply(mrs_data_mat, 2, t_test_calc, group)
  
  pvals     <- sapply(t_test_res_list, \(x) x$p.value)
  pvals_mrs <- vec2mrs_data(pvals, mrs_data = mrs_data, fd = TRUE)
  pvals_log <- -log10(pvals)
  pvals_log_mrs <- vec2mrs_data(pvals_log, mrs_data = mrs_data, fd = TRUE)
  
  tvals     <- sapply(t_test_res_list, \(x) x$statistic)
  tvals_mrs <- vec2mrs_data(tvals, mrs_data = mrs_data, fd = TRUE)
  
  return(list(p_value = pvals, p_value_mrs = pvals_mrs,
              p_value_log = pvals_log, p_value_log_mrs = pvals_log_mrs,
              t_stat = tvals, t_stat_mrs = tvals_mrs,
              t_test_res_list = t_test_res_list))
}

#' Simulate an example fMRS dataset for a block design fMRS experiment and
#' export a BIDS structure.
#' @param output_dir output directory for the BIDS data. Defaults to : 
#' "HOME/sim_fmrs_dataset/data".
#' @export
spant_sim_fmrs_dataset <- function(output_dir = NULL) {
  
  if (is.null(output_dir)) {
    output_dir <- file.path("~", "sim_fmrs_dataset", "data")
  }
  
  seq_tr      <- 2    # set the sequence TR
  N_scans     <- 448  # just under 15 mins with TR = 2s
  bz_inhom_lb <- 4    # Gaussian line-broadening to simulate B0 inhomogeneity
  bold_lb_hz  <- 0.0  # linewidth differences in Hz from BOLD T2* effect
  ss_spec_snr <- 50   # single shot spectral SNR
  subjects    <- 5    # number of subject scans to generate in BIDS format
  set.seed(1)         # random number generator seed used for noise samples
  
  # Make a data frame containing a single row of basis signal amplitudes.
  # Metabolite values are for visual cortex from Bednarik et al 2015 Table 1.
  # Note Alanine and Glycine are not listed in the table and therefore set to 0.
  basis_amps <- data.frame("ala"    = 0.00, "asc"    = 0.96, "asp"   = 3.58,
                           "cr"     = 4.22, "gaba"   = 1.03, "glc"   = 0.62,
                           "gln"    = 2.79, "gly"    = 0.00, "gsh"   = 1.09,
                           "glu"    = 8.59, "gpc"    = 0.54, "ins"   = 6.08,
                           "lac"    = 1.01, "naa"    = 11.9, "naag"  = 1.32,
                           "pch"    = 0.40, "pcr"    = 3.34, "peth"  = 0.93,
                           "sins"   = 0.27, "tau"    = 1.27, "lip09" = 0.00,
                           "lip13a" = 0.00, "lip13b" = 0.00, "lip20" = 0.00,
                           "mm09"   = 4.00, "mm12"   = 4.00, "mm14"  = 4.00,
                           "mm17"   = 4.00, "mm20"   = 4.00)
  
  # Duplicate the row N_scans times to make a table of values
  basis_amps <- basis_amps[rep(1, N_scans),]
  
  # simulate two 120 second blocks of stimulation starting at 100 and 500 seconds
  onsets    <- c(100, 500)
  durations <- rep(120, 2)
  
  # generate a dummy mrs_data object for generating the regressors
  mrs_data_dummy <- set_Ntrans(set_tr(sim_zero(dyns = N_scans), seq_tr),
                               N_scans)
  
  # generate metabolite response functions assuming simple trapezoidal shapes
  glu_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                         rise_t = 120, fall_t = 150)
  
  lac_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                         rise_t = 120, fall_t = 150)
  
  asp_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                         rise_t = 120, fall_t = 150)
  
  glc_rf <- gen_trap_reg(onsets, durations, mrs_data = mrs_data_dummy,
                         rise_t = 120, fall_t = 150)
  
  bold_rf <- gen_bold_reg(onsets, durations, mrs_data = mrs_data_dummy)
  
  # update the amplitudes according to predicted changes from Bednarik et at 2015
  # Table 1
  lac_perc_change <-  29.6
  glu_perc_change <-  3.3
  glc_perc_change <- -16.0
  asp_perc_change <- -5.4
  
  # update metabolite data frame to have dynamic metabolite values
  basis_amps$glu <- basis_amps$glu * (glu_rf$stim * glu_perc_change / 100 + 1)
  basis_amps$lac <- basis_amps$lac * (lac_rf$stim * lac_perc_change / 100 + 1)
  basis_amps$asp <- basis_amps$asp * (asp_rf$stim * asp_perc_change / 100 + 1)
  basis_amps$glc <- basis_amps$glc * (glc_rf$stim * glc_perc_change / 100 + 1)
  
  # simulate a typical basis for TE=28ms semi-LASER acquisition at 3T
  acq_paras <- def_acq_paras(ft = 127.8e6)
  basis     <- sim_basis(names(basis_amps), pul_seq = seq_slaser_ideal,
                         TE1 = 0.008, TE2 = 0.011, TE3 = 0.009)
  
  # apply basis amplitudes to the basis set to generate a simulated fMRS dataset
  mrs_dyn_orig <- basis2dyn_mrs_data(basis, basis_amps, seq_tr)
  
  # broaden basis to simulate B0 inhomogeneity, apply any addition BOLD related 
  # narrowing
  bold_lb_dyn <- (1 - bold_rf$stim_bold) * bold_lb_hz
  mrs_dyn     <- lb(lb(mrs_dyn_orig, bz_inhom_lb), bold_lb_dyn, 0)
  
  # duplicate the data to generate multiple subjects with different noise samples
  mrs_dyn_list <- rep(list(mrs_dyn), subjects)
  
  # apply additional 10 Hz line-broadening to sub-02
  mrs_dyn_list[[2]] <- lb(mrs_dyn_list[[2]], 10)
  
  # add noise
  mrs_dyn_list <- add_noise_spec_snr(mrs_dyn_list, ss_spec_snr,
                                     ref_data = mrs_dyn_list[[1]])
  
  # export to BIDS structure
  mrs_data2bids(mrs_dyn_list, output_dir, skip_existing = FALSE)
}

#' Perform group-level spectral GLM analysis of an fMRS dataset.
#' @param regressor_df a data frame containing temporal regressors to be applied
#' to each spectral datapoint.
#' @param analysis_dir directory containing preprocessed data generated by
#' the preproc_svs_dataset function.
#' @param exclude_labels vector of labels of scans to exclude, eg poor quality
#' data.
#' @param labels labels to describe each data set.
#' @export
glm_spec_fmrs_group <- function(regressor_df, analysis_dir = "spant_analysis",
                                exclude_labels = NULL, labels = NULL) {
  
  # check preproc_dir exists
  if (!dir.exists(analysis_dir)) stop("analysis_dir not found.")
  
  # directory for output
  spec_glm_group_dir <- file.path(analysis_dir, "spec_glm", "group_level")
  if (!dir.exists(spec_glm_group_dir)) dir.create(spec_glm_group_dir)
  
  if (is.null(labels)) {
    # discover data labels from file names
    fnames <- Sys.glob(file.path(analysis_dir, "spec_glm", "proc_fmrs",
                                               "*.nii.gz"))
    labels <- basename(fnames)
    labels <- tools::file_path_sans_ext(tools::file_path_sans_ext(labels))
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
    # remove any unfound labels
    exclude_inds <- exclude_inds[!is.na(exclude_inds)]
    if (length(exclude_inds) > 0) labels <- labels[-exclude_inds]
  }
  
  # read the processed data
  Nscans      <- length(labels)
  preproc_mrs <- vector(mode = "list", length = Nscans)
  for (n in 1:Nscans) {
    fname <- paste0(labels[n], ".nii.gz")
    fpath <- file.path(analysis_dir, "spec_glm", "proc_fmrs", fname)
    preproc_mrs[[n]] <- read_mrs(fpath)
  }
  
  # modelling
  glm_spec_group_res <- glm_spec_group_level(preproc_mrs, regressor_df)
  
  # write the results
  output <- list(glm_spec_group_res = glm_spec_group_res, labels = labels,
                 exclude_labels = exclude_labels, regressor_df = regressor_df,
                 acq_paras = get_acq_paras(preproc_mrs[[1]]))
  
  saveRDS(output, file.path(spec_glm_group_dir, "group_fit.rds"))
}

#' Test a group-level spectral GLM linear hypothesis.
#' @param hmat linear hypothesis matrix.
#' @param analysis_dir directory containing preprocessed data generated by
#' the preproc_svs_dataset function.
#' @export
glm_spec_group_linhyp <- function(hmat, analysis_dir = "spant_analysis") {
  
  group_dir   <- file.path(analysis_dir, "spec_glm", "group_level")
  group_fit_f <- file.path(group_dir, "group_fit.rds")
  group_fit   <- readRDS(group_fit_f)
  
  lin_hyp_res <- lapply(group_fit$glm_spec_group_res, car::linearHypothesis,
                        hmat)
  
  F_vals      <- sapply(lin_hyp_res, \(x) x$`F`[2])
  p_vals      <- sapply(lin_hyp_res, \(x) x$`Pr(>F)`[2])
  p_vals_log  <- -log10(p_vals)
  
  acq_paras <- group_fit$acq_paras
  
  p_vals_mrs <- vec2mrs_data(p_vals, fs = acq_paras$fs,
                             ft = acq_paras$ft, ref = acq_paras$ref,
                             nuc = acq_paras$nuc, fd = TRUE)
  
  p_vals_log_mrs <- vec2mrs_data(p_vals_log, fs = acq_paras$fs,
                                 ft = acq_paras$ft, ref = acq_paras$ref,
                                 nuc = acq_paras$nuc, fd = TRUE)
  
  F_vals_mrs <- vec2mrs_data(F_vals, fs = acq_paras$fs, ft = acq_paras$ft,
                             ref = acq_paras$ref, nuc = acq_paras$nuc,
                             fd = TRUE)
  
  return(list(p_vals_log_mrs = p_vals_log_mrs, p_vals_mrs = p_vals_mrs, 
              F_vals_mrs = F_vals_mrs))
}

glm_spec_group_level <- function(mrs_data_list, regressor_df) {
  Ngroup <- length(mrs_data_list)
  group_regressor_df <- gen_group_reg(regressor_df, Ngroup)
  mrs_data_dyn <- append_dyns(mrs_data_list) 
  mrs_mat <- mrs_data2spec_mat(mrs_data_dyn)
  
  # drop the time column if present
  group_regressor_df <- group_regressor_df[, !names(group_regressor_df) %in% 
                                             c("time"), drop = FALSE]
  lm_res_list <- vector("list", ncol(mrs_mat))
  for (n in 1:ncol(mrs_mat)) {
    lm_res_list[[n]] <- stats::lm(mrs_mat[, n] ~ . + 0, group_regressor_df)
  }
  
  return(lm_res_list)
}