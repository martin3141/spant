#' Simple pulse and acquire sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @return a list of resonance amplitudes and frequencies.
#' @export
seq_pulse_acquire <- function(spin_params, ft, ref) {
  sys <- spin_sys(spin_params, ft, ref)
  
  sys$rho <- gen_F(sys, "z")
   
  angle <- 90
  Fx <- gen_F(sys, "x", "1H")
  lhs_pulse <- complexplus::matexp(-Fx * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fx * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # acquire
  acquire(sys, detect = "1H")
}

#' Simple pulse and acquire sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @return a list of resonance amplitudes and frequencies.
#' @export
seq_pulse_acquire_31p <- function(spin_params, ft, ref) {
  sys <- spin_sys(spin_params, ft, ref)
  
  sys$rho <- gen_F(sys, "z", "31P")
   
  angle <- 90
  Fx <- gen_F(sys, "x", "31P")
  lhs_pulse <- complexplus::matexp(-Fx * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fx * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # acquire
  acquire(sys, detect = "31P")
}

#' PRESS sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE1 TE1 sequence parameter in seconds (TE=TE1+TE2).
#' @param TE2 TE2 sequence parameter in seconds.
#' @return a list of resonance amplitudes and frequencies.
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
  lhs_pulse <- complexplus::matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fy * 1i * angle * pi / 180)
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
  t = TE2/2
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # acquire
  acquire(sys,detect = "1H")
}

#' Spin echo sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE echo time in seconds.
#' @return a list of resonance amplitudes and frequencies.
#' @export
seq_spin_echo_ideal <- function(spin_params, ft, ref, TE = 0.03) {
  sys <- spin_sys(spin_params, ft, ref)
  
  sys$rho <- gen_F(sys, "z")
   
  angle <- 90
  Fx <- gen_F(sys, "x", "1H")
  lhs_pulse <- complexplus::matexp(-Fx * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fx * 1i * angle * pi / 180)
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
  Fy <- gen_F(sys, "y", "1H")
  lhs_pulse <- complexplus::matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # acquire
  acquire(sys, detect = "1H")
}

#' Spin echo sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE echo time in seconds.
#' @return a list of resonance amplitudes and frequencies.
#' @export
seq_spin_echo_ideal_31p <- function(spin_params, ft, ref, TE = 0.03) {
  sys <- spin_sys(spin_params, ft, ref)
  
  sys$rho <- gen_F(sys, "z", "31P")
   
  angle <- 90
  Fx <- gen_F(sys, "x", "31P")
  lhs_pulse <- complexplus::matexp(-Fx * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fx * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  t <- TE / 2
  # find the inverse of the eigenvector matrix
  eig_vec_inv <- solve(sys$H_eig_vecs)
  lhs <- sys$H_eig_vecs %*% diag(exp(sys$H_eig_vals * 2i * pi * t)) %*% 
         eig_vec_inv
  
  rhs <- sys$H_eig_vecs %*% diag(exp(-sys$H_eig_vals * 2i * pi * t)) %*%
         eig_vec_inv
  
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply pulse
  angle <- 180
  Fy <- gen_F(sys,"y","31P")
  lhs_pulse <- complexplus::matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # acquire
  acquire(sys, detect = "31P")
}

#' CPMG style sequence with ideal pulses.
#' @param spin_params spin system definition.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param TE echo time in seconds.
#' @param echoes number of echoes.
#' @return a list of resonance amplitudes and frequencies.
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
  lhs_pulse <- complexplus::matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fy * 1i * angle * pi / 180)
    
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
#' @return a list of resonance amplitudes and frequencies.
#' @export
seq_mega_press_ideal <- function(spin_params, ft, ref, ed_freq = 1.89, TE1 = 0.015, 
                             TE2 = 0.053, BW = 110, steps = 50) {
 
  # ed pulse duration 
  duration <- 1.53 / BW # 180 deg 1% Guass (p245 de Graff)
  #print(duration)
  
  sys <- spin_sys(spin_params, ft, ed_freq)
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
  lhs_pulse <- complexplus::matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
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
    lhs_pulse_gaus <- complexplus::matexp(-Fx * 1i * angle)
    rhs_pulse_gaus <- complexplus::matexp(Fx * 1i * angle)
    sys$rho <- lhs_pulse_gaus %*% sys$rho %*% rhs_pulse_gaus
    sys$rho <- lhs_dt %*% sys$rho %*% rhs_dt
  }
  
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
  
  # apply TE2/4 delay
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # apply selective 180 
  for (x in pulse) {
    angle <- x
    Fx <- gen_F(sys, "x", "1H")
    lhs_pulse_gaus <- complexplus::matexp(-Fx * 1i * angle)
    rhs_pulse_gaus <- complexplus::matexp(Fx * 1i * angle)
    sys$rho <- lhs_pulse_gaus %*% sys$rho %*% rhs_pulse_gaus
    sys$rho <- lhs_dt %*% sys$rho %*% rhs_dt
  }
  
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
#' @return a list of resonance amplitudes and frequencies.
#' @export
seq_steam_ideal <- function(spin_params, ft, ref, TE = 0.03, TM = 0.02) {
  
  # TODO, compare results with the method described by Young et al JMR 140,
  # 146-152 (1999) where a rotation about the z-axis is used instead of 
  # the Fxy operator
  
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
  lhs_x_pulse <- complexplus::matexp(-Fx * 1i * angle * pi / 180)
  rhs_x_pulse <- complexplus::matexp(Fx * 1i * angle * pi / 180)
  
  basis_size <- prod(sys$spin_n * 2 + 1)
  rho_combined <- matrix(0, basis_size, basis_size)
  
  # phase cycling loop 
  for (n in 0:3) {
    phase <- n * 360 / 4
    #phase <- 270
    
    # first and third 90 pulse operator
    angle <- 90
    Fxy <- gen_F_xy(sys, phase, "1H")
    # these lines give warnings because cos(90*pi/180) isn't exactly zero
    # but that's ok - so suppressWarnings
    lhs_xy_pulse <- suppressWarnings(complexplus::matexp(-Fxy * 1i * angle * pi / 180))
    rhs_xy_pulse <- suppressWarnings(complexplus::matexp(Fxy * 1i * angle * pi / 180))
    
    # first 90
    rho <- lhs_xy_pulse %*% sys$rho %*% rhs_xy_pulse
    # evole TE/2
    rho <- lhs_half_te %*% rho %*% rhs_half_te
    # second 90
    rho <- lhs_x_pulse %*% rho %*% rhs_x_pulse
    # zero non-zero-order coherences
    rho <- zero_nzoc(sys, rho)
    # evole TM
    rho <- lhs_tm %*% rho %*% rhs_tm
    # third 90
    rho <- lhs_xy_pulse %*% rho %*% rhs_xy_pulse
    # evole TE/2
    rho <- lhs_half_te %*% rho %*% rhs_half_te
    
    rho_combined <- rho_combined + rho / 4
  }
  
  sys$rho <- rho_combined
  
  # acquire
  acquire(sys, detect = "1H", rec_phase = 0)
}