#' Simple pulse and acquire sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param nuc acquisition nucleus.
#' @param acq_delay delay between excitation and acquisition.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_pulse_acquire <- function(spin_params, ft, ref, nuc = "1H",
                              acq_delay = 0) {
  
  sys <- spin_sys(spin_params, ft, ref)
  
  sys$rho <- gen_F(sys, "z", nuc)
   
  angle <- 90
  Fx <- gen_F(sys, "x",  nuc)
  lhs_pulse <- matexp(-Fx * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fx * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  if (acq_delay > 0) {
    # apply delay
    t <- acq_delay
    # find the inverse of the eigenvector matrix
    eig_vec_inv <- solve(sys$H_eig_vecs)
    lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
      eig_vec_inv
    rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
      eig_vec_inv
    sys$rho <- lhs %*% sys$rho %*% rhs
  }
  
  # acquire
  acquire(sys, detect = nuc)
}

#' PRESS sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE1 TE1 sequence parameter in seconds (TE=TE1+TE2).
#' @param TE2 TE2 sequence parameter in seconds.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_press_ideal <- function(spin_params, ft, ref, TE1 = 0.01, TE2 = 0.02) {
  
  sys <- spin_sys(spin_params, ft, ref)
  sys$rho <- -gen_F(sys, "y", "1H")
  
  # apply delay
  t <- TE1 / 2
  
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  angle <- 180
  Fy <- gen_F(sys, "y", "1H")
  lhs_pulse <- matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t <- (TE1 + TE2) / 2
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
         eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
         eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t <- TE2 / 2
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # acquire
  acquire(sys, detect = "1H")
}

#' Spin echo sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param nuc acquisition nucleus.
#' @param TE echo time in seconds.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_spin_echo_ideal <- function(spin_params, ft, ref, nuc = "1H", TE = 0.03) {
  sys <- spin_sys(spin_params, ft, ref)
  
  sys$rho <- gen_F(sys, "z")
   
  angle <- 90
  Fx <- gen_F(sys, "x", nuc)
  lhs_pulse <- matexp(-Fx * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fx * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t = TE / 2
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
         eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  angle <- 180
  Fy <- gen_F(sys, "y", nuc)
  lhs_pulse <- matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # acquire
  acquire(sys, detect = nuc)
}

#' CPMG style sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE echo time in seconds.
#' @param echoes number of echoes.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_cpmg_ideal <- function(spin_params, ft, ref, TE = 0.03, echoes = 4) {
  sys <- spin_sys(spin_params, ft, ref)
  sys$rho <- -gen_F(sys, "y", "1H")
  
  # calc delay operator
  t <- TE / (echoes * 2)
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  # calc pulse operator
  angle <- 180
  Fy <- gen_F(sys, "y", "1H")
  lhs_pulse <- matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fy * 1i * angle * pi / 180)
    
  sys$rho <- lhs %*% sys$rho %*% rhs
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  if (echoes > 1) {
    for (n in 1:(echoes - 1)) {
      sys$rho <- lhs %*% sys$rho %*% rhs
      sys$rho <- lhs %*% sys$rho %*% rhs
      sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
    }
  }
    
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # acquire
  acq_res <- acquire(sys, detect = "1H")
}

#' MEGA-PRESS sequence with ideal localisation pulses and Gaussian shaped
#' editing pulse.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param ed_freq editing pulse frequency in ppm.
#' @param TE1 TE1 sequence parameter in seconds (TE=TE1+TE2).
#' @param TE2 TE2 sequence parameter in seconds.
#' @param BW editing pulse bandwidth in Hz.
#' @param steps number of hard pulses used to approximate the editing pulse.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_mega_press_ideal <- function(spin_params, ft, ref, ed_freq = 1.89,
                                 TE1 = 0.015, TE2 = 0.053, BW = 110,
                                 steps = 50) {
  
  # ed pulse duration 
  duration <- 1.53 / BW # 180 deg 1% Guass (p245 de Graff)
  
  sys <- spin_sys(spin_params, ft, ed_freq)
  sys$rho <- -gen_F(sys, "y", "1H")
  
  # filter -1
  sys$rho <- coherence_filter(sys, sys$rho, -1)
  
  # apply delay
  t <- TE1 / 2
  
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  angle <- 180
  Fy <- gen_F(sys, "y", "1H")
  lhs_pulse <- matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # filter 1
  sys$rho <- coherence_filter(sys, sys$rho, 1)
  
  # apply delay
  t <- (TE1 + TE2 / 2) / 2 - duration / 2
  if (t < 0) stop("Error, negative delay duration required.")
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply selective 180
  pulse <- get_guassian_pulse(180, steps)
  dt <- duration / steps
  
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs_dt <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * dt)) %*%
    eig_vec_inv
  
  rhs_dt <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * dt)) %*%
         eig_vec_inv
  
  for (x in pulse) {
    angle <- x
    Fx <- gen_F(sys, "y", "1H")
    lhs_pulse_gaus <- matexp(-Fx * 1i * angle)
    rhs_pulse_gaus <- matexp( Fx * 1i * angle)
    sys$rho <- lhs_pulse_gaus %*% sys$rho %*% rhs_pulse_gaus
    sys$rho <- lhs_dt %*% sys$rho %*% rhs_dt
  }
  
  # filter 1
  sys$rho <- coherence_filter(sys, sys$rho, 1)
  
  # apply TE2/4 delay
  t <- TE2 / 4 - duration / 2
  if (t < 0) stop("Error, negative delay duration required.")
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # filter -1
  sys$rho <- coherence_filter(sys, sys$rho, -1)
  
  # apply TE2/4 delay
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply selective 180 
  for (x in pulse) {
    angle <- x
    Fx <- gen_F(sys, "x", "1H")
    lhs_pulse_gaus <- matexp(-Fx * 1i * angle)
    rhs_pulse_gaus <- matexp( Fx * 1i * angle)
    sys$rho <- lhs_pulse_gaus %*% sys$rho %*% rhs_pulse_gaus
    sys$rho <- lhs_dt %*% sys$rho %*% rhs_dt
  }
  
  # filter -1
  sys$rho <- coherence_filter(sys, sys$rho, -1)
  
  # apply TE2/4 delay
  sys$rho <- lhs %*% sys$rho %*% rhs

  # acquire
  acq_res <- acquire(sys, detect = "1H")
  
  # shift back to requested ref
  acq_res$freqs <- acq_res$freqs + (ref - ed_freq) * ft / 1e6
  acq_res
}

#' STEAM sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE sequence parameter in seconds.
#' @param TM sequence parameter in seconds.
#' @param amp_scale amplitude scaling factor. Set to 2 (default) to ensure
#' correct scaling for water reference scaling. Set to 1 to maintain the
#' inherent loss of signal associated with STEAM.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_steam_ideal <- function(spin_params, ft, ref, TE = 0.03, TM = 0.02,
                            amp_scale = 2) {
  
  sys <- spin_sys(spin_params, ft, ref)
  sys$rho <- gen_F(sys, "z")
  
  # TE/2 delay operator 
  eig_vec_inv <- solve(sys$H_eig_vecs)
  t <- (TE / 2)
  lhs_half_te <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  rhs_half_te <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  # TM delay operator 
  t <- TM
  lhs_tm <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  rhs_tm <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  # 90 degree pulse operator
  angle <- 90
  Fx <- gen_F(sys, "x", "1H")
  lhs_x_pulse <- matexp(-Fx * 1i * angle * pi / 180)
  rhs_x_pulse <- matexp( Fx * 1i * angle * pi / 180)
  
  basis_size   <- prod(sys$spin_n * 2 + 1)
  rho_combined <- matrix(0, basis_size, basis_size)
  
  # phase cycling loop
  for (n in 0:3) {
    phase <- n * 90
    
    # first and third 90 pulse operator
    angle <- 90
    Fxy <- gen_F_xy(sys, phase, "1H")
    
    # XY rotation operator
    lhs_xy_pulse <- matexp(-Fxy * 1i * angle * pi / 180)
    rhs_xy_pulse <- matexp( Fxy * 1i * angle * pi / 180)
    
    # first 90 plus rotation
    rho <- lhs_xy_pulse %*% sys$rho %*% rhs_xy_pulse
    
    # evolve TE/2
    rho <- lhs_half_te %*% rho %*% rhs_half_te
    
    # second 90
    rho <- lhs_x_pulse %*% rho %*% rhs_x_pulse
    
    # zero non-zero-order coherences
    rho <- coherence_filter(sys, rho)
    
    # evolve TM
    rho <- lhs_tm %*% rho %*% rhs_tm
    
    # third 90 plus rotation
    rho <- lhs_xy_pulse %*% rho %*% rhs_xy_pulse
    
    # evolve TE/2
    rho <- lhs_half_te %*% rho %*% rhs_half_te
    
    rho_combined <- rho_combined + rho / 4
  }
  
  sys$rho <- rho_combined
  
  # acquire, and double the output intensity to ensure consistent concentration
  # scaling
  acquire(sys, rec_phase = 180, detect = "1H", amp_scale = amp_scale)
}

#' STEAM sequence with ideal pulses using the z-rotation gradient simulation
#' method described by Young et al JMR 140, 146-152 (1999).
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE sequence parameter in seconds.
#' @param TM sequence parameter in seconds.
#' @param amp_scale amplitude scaling factor. Set to 2 (default) to ensure
#' correct scaling for water reference scaling. Set to 1 to maintain the
#' inherent loss of signal associated with STEAM.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_steam_ideal_young <- function(spin_params, ft, ref, TE = 0.03, TM = 0.02,
                                  amp_scale = 2) {
  
  sys <- spin_sys(spin_params, ft, ref)
  Fz  <- gen_F(sys, "z")
  sys$rho <- Fz
  
  # TE/2 delay operator 
  eig_vec_inv <- solve(sys$H_eig_vecs)
  t <- (TE / 2)
  lhs_half_te <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  rhs_half_te <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  # TM delay operator 
  t <- TM
  lhs_tm <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  rhs_tm <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  # 90 degree x pulse operator
  angle <- 90
  Fx <- gen_F(sys, "x", "1H")
  lhs_x_pulse <- matexp(-Fx * 1i * angle * pi / 180)
  rhs_x_pulse <- matexp( Fx * 1i * angle * pi / 180)
  
  basis_size   <- prod(sys$spin_n * 2 + 1)
  rho_combined <- matrix(0, basis_size, basis_size)
  
  # apply first 90
  sys$rho <- lhs_x_pulse %*% sys$rho %*% rhs_x_pulse
  
  # evolve TE/2
  sys$rho <- lhs_half_te %*% sys$rho %*% rhs_half_te
  
  # phase cycling loop
  for (n in 0:3) {
    phase_ang <- n * 90
    
    # rotate about z
    lhs_z_rot <- matexp(-Fz * 1i * phase_ang * pi / 180)
    rhs_z_rot <- matexp( Fz * 1i * phase_ang * pi / 180)
    rho <- lhs_z_rot %*% sys$rho %*% rhs_z_rot
    
    # second 90
    rho <- lhs_x_pulse %*% rho %*% rhs_x_pulse
    
    # zero non-zero-order coherences
    rho <- coherence_filter(sys, rho)
    
    # evolve TM
    rho <- lhs_tm %*% rho %*% rhs_tm
    
    # third 90
    rho <- lhs_x_pulse %*% rho %*% rhs_x_pulse
    
    # rotate about z
    rho <- lhs_z_rot %*% rho %*% rhs_z_rot

    rho_combined <- rho_combined + rho / 4
  }
  
  # evolve TE/2
  sys$rho <- lhs_half_te %*% rho_combined %*% rhs_half_te
  
  # acquire, and double the output intensity to ensure consistent concentration
  # scaling
  acquire(sys, detect = "1H", rec_phase = 180, amp_scale = amp_scale)
}

#' STEAM sequence with ideal pulses and coherence order filtering to simulate
#' gradient crushers.
#' 
#' See Landheer et al NMR Biomed 2021 34(5):e4129 and Landheer et al MRM 2019
#' Apr;81(4):2209-2222 for more details on the coherence order filtering method.
#' 
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE sequence parameter in seconds.
#' @param TM sequence parameter in seconds.
#' @param amp_scale amplitude scaling factor. Set to 2 (default) to ensure
#' correct scaling for water reference scaling. Set to 1 to maintain the
#' inherent loss of signal associated with STEAM.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_steam_ideal_cof <- function(spin_params, ft, ref, TE = 0.03, TM = 0.02,
                                amp_scale = 2) {
  
  sys <- spin_sys(spin_params, ft, ref)
  Fz  <- gen_F(sys, "z")
  sys$rho <- Fz
  
  # TE/2 delay operator 
  eig_vec_inv <- solve(sys$H_eig_vecs)
  t <- (TE / 2)
  lhs_half_te <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  rhs_half_te <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  # TM delay operator 
  t <- TM
  lhs_tm <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  rhs_tm <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  # 90 degree x pulse operator
  angle <- 90
  Fx <- gen_F(sys, "x", "1H")
  lhs_x_pulse <- matexp(-Fx * 1i * angle * pi / 180)
  rhs_x_pulse <- matexp( Fx * 1i * angle * pi / 180)
  
  basis_size   <- prod(sys$spin_n * 2 + 1)
  rho_combined <- matrix(0, basis_size, basis_size)
  
  # apply first 90
  sys$rho <- lhs_x_pulse %*% sys$rho %*% rhs_x_pulse
  
  # filter +1
  sys$rho <- coherence_filter(sys, sys$rho, 1)
  
  # evolve TE/2
  sys$rho <- lhs_half_te %*% sys$rho %*% rhs_half_te
  
  # second 90
  sys$rho <- lhs_x_pulse %*% sys$rho %*% rhs_x_pulse
  
  # filter 0
  sys$rho <- coherence_filter(sys, sys$rho, 0)
  
  # evolve TM
  sys$rho <- lhs_tm %*% sys$rho %*% rhs_tm
  
  # third 90
  sys$rho <- lhs_x_pulse %*% sys$rho %*% rhs_x_pulse
  
  # filter -1
  sys$rho <- coherence_filter(sys, sys$rho, -1)
  
  # evolve TE/2
  sys$rho <- lhs_half_te %*% sys$rho %*% rhs_half_te
  
  # acquire, and double the output intensity to ensure consistent concentration
  # scaling
  acquire(sys, rec_phase = 180, detect = "1H", amp_scale = amp_scale)
}

#' sLASER sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE1 first echo time (between exc. and 1st echo) in seconds.
#' @param TE2 second echo time (between 2nd echo and 4th echo) in seconds.
#' @param TE3 third echo time (between 4th echo and 5th echo) in seconds.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_slaser_ideal <- function(spin_params, ft, ref, TE1 = 0.008, TE2 = 0.011,
                             TE3 = 0.009) {
  
  sys <- spin_sys(spin_params, ft, ref)
  
  sys$rho <- gen_F(sys, "z")
  
  angle <- 90
  Fx <- gen_F(sys, "x", "1H")
  lhs_pulse <- matexp(-Fx * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fx * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t <- TE1 / 2
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  angle <- 180
  Fy <- gen_F(sys, "y", "1H")
  lhs_pulse <- matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t = TE1 / 2 + TE2 / 4
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  angle <- 180
  Fy <- gen_F(sys, "y", "1H")
  lhs_pulse <- matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t <- 2 * TE2 / 4
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  angle <- 180
  Fy <- gen_F(sys, "y", "1H")
  lhs_pulse <- matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t <- TE2 / 4 + TE3 /2
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  angle <- 180
  Fy <- gen_F(sys, "y", "1H")
  lhs_pulse <- matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- matexp( Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t <- TE3 /2
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # acquire
  acquire(sys, detect = "1H")
}

#' PRESS sequence with shaped refocusing pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE1 TE1 sequence parameter in seconds (TE=TE1+TE2).
#' @param TE2 TE2 sequence parameter in seconds.
#' @param pulse_file path to refocusing pulse file.
#' @param pulse_dur refocusing pulse duration.
#' @param pulse_file_format file format for the refocusing pulse.
#' @param refoc_flip_angle refocusing pulse flip angle in degrees (defaults to
#' 180).
#' @param xy_pulse_ppm a vector of ppm values for the offset of each
#' sub-simulation.
#' @param resamp option to resample the pulse shape.
#' @param fs_resamp sampling frequency (Hz) to resample.
#' @return list of resonance amplitudes and frequencies.
#' @export
seq_press_2d_shaped <- function(spin_params, ft, ref, TE1 = 0.01, TE2 = 0.02,
                                pulse_file, pulse_dur, pulse_file_format,
                                refoc_flip_angle = 180, xy_pulse_ppm = NULL,
                                resamp = TRUE, fs_resamp = 1e-4) {
  
  if (is.null(xy_pulse_ppm)) xy_pulse_ppm <- seq(from = -4, to = 10,
                                                 length.out = 10)
  
  # read in a 180 degree refocusing pulse
  if (pulse_file_format == "pta") {
    pulse <- read_pulse_pta(pulse_file)
  } else if (pulse_file_format == "bruker") {
    pulse <- read_pulse_bruker(pulse_file)
  } else if (pulse_file_format == "ascii") {
    pulse <- read_pulse_ascii(pulse_file)
  } else {
    stop("Pulse file format not recognised.")
  }
  
  if (resamp) {
    n_resamp   <- round(pulse_dur / fs_resamp)
    mag_interp <- stats::approx(pulse$data$mag, n = n_resamp)$y
    pha_interp <- stats::approx(pulse$data$pha, n = n_resamp)$y
    pulse$data <- data.frame(mag = mag_interp, pha = pha_interp)
  }
  
  # initialize spin system
  sys      <- spin_sys(spin_params, ft, ref)
  sys$rho  <- -gen_F(sys, "y", "1H")
  
  # precomputed for speed
  H_mat_jc <- sys$H_mat_jc
 
  # precomputed for speed
  spin_num <- get_spin_num(spin_params$nucleus)
  precomp_Iz <- vector(mode = "list", length = length(spin_num))
  for (n in 1:length(spin_num)) precomp_Iz[[n]] <- gen_I(n, spin_num, "z")
  
  # filter -1
  sys$rho <- coherence_filter(sys, sys$rho, -1)
  
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  
  # prepare first delay operators
  t <- TE1 / 2 - pulse_dur / 2
  lhs_d1 <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  rhs_d1 <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  # prepare second delay operators
  t <- (TE1 + TE2) / 2 - pulse_dur
  lhs_d2 <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  rhs_d2 <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*% 
    eig_vec_inv
  
  # prepare third delay operators
  t <- TE2 / 2 - pulse_dur / 2
  lhs_d3 <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  rhs_d3 <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
    eig_vec_inv
  
  # apply first delay
  sys$rho <- lhs_d1 %*% sys$rho %*% rhs_d1
    
  basis_size <- prod(sys$spin_n * 2 + 1)
  
  dt <- pulse_dur / nrow(pulse$data) 

  pulse_cplx <- pulse$data$mag * exp(1i * pulse$data$pha)
  pulse$data$mag <- pulse$data$mag / sum(Re(pulse_cplx)) * pi * 
                    refoc_flip_angle / 180
  
  rho_spatial <- matrix(0, basis_size, basis_size)
  start_rho <- sys$rho
  
  # precalculate pulse operators 
  lhs_pulse <- vector(mode = "list", length = nrow(pulse$data))
  rhs_pulse <- vector(mode = "list", length = nrow(pulse$data))
  for (n in 1:nrow(pulse$data)) {
    Fxy <- gen_F_xy(sys, pulse$data$pha[n] * 180 / pi + 90, "1H")
    lhs_pulse[[n]] <- matexp(-Fxy * 1i * pulse$data$mag[n])
    rhs_pulse[[n]] <- matexp( Fxy * 1i * pulse$data$mag[n])
  }
  
  lhs_dt <- vector(mode = "list", length = length(xy_pulse_ppm))
  rhs_dt <- vector(mode = "list", length = length(xy_pulse_ppm))
 
  for (m in 1:length(xy_pulse_ppm)) {
    sys <- spin_sys(spin_params, ft, xy_pulse_ppm[m], precomp_jc_H = H_mat_jc,
                    precomp_Iz = precomp_Iz)
    eig_vec_inv <- solve(sys$H_eig_vecs)
    
    lhs_dt[[m]] <- sys$H_eig_vecs %*%
                   diag(exp(sys$H_eig_vals * 2i * pi * dt)) %*% eig_vec_inv
    rhs_dt[[m]] <- sys$H_eig_vecs %*%
                   diag(exp(-sys$H_eig_vals * 2i * pi * dt)) %*% eig_vec_inv
  }
  
  for (m in 1:length(xy_pulse_ppm)) {
      
      rho <- start_rho
      # first 180
      for (n in 1:nrow(pulse$data)) {
        rho <- lhs_pulse[[n]] %*% rho %*% rhs_pulse[[n]]
        rho <- lhs_dt[[m]]    %*% rho %*% rhs_dt[[m]]
      }
      rho_spatial <- rho_spatial + rho / length(xy_pulse_ppm)
    }
      
    # filter +1
    rho <- coherence_filter(sys, rho_spatial, 1)
    
    # apply second delay
    rho <- lhs_d2 %*% rho %*% rhs_d2
    
    rho_spatial <- matrix(0, basis_size, basis_size)
    start_rho   <- rho
    for (m in 1:length(xy_pulse_ppm)) {
      
      rho <- start_rho
      
      # second 180
      for (n in 1:nrow(pulse$data)) {
        rho <- lhs_pulse[[n]] %*% rho %*% rhs_pulse[[n]]
        rho <- lhs_dt[[m]]    %*% rho %*% rhs_dt[[m]]
      }
      rho_spatial <- rho_spatial + rho / length(xy_pulse_ppm)
    }
    
    # filter -1
    rho <- coherence_filter(sys, rho_spatial, -1)
    
    # reset ref 
    sys <- spin_sys(spin_params, ft, ref)
    
    # apply third delay
    sys$rho <- lhs_d3 %*% rho %*% rhs_d3 
    
    # acquire
    acquire(sys, detect = "1H")
}