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
#' @param TE1 TE1 sequence parameter in seconds.
#' @param TE2 TE2 sequence parameter in seconds.
#' @return a list of resonance amplitudes and frequencies.
#' @export
seq_press_ideal <- function(spin_params, ft, ref, TE1 = 0.01, TE2 = 0.02) {
  
  sys <- spin_sys(spin_params, ft, ref)
  sys$rho <- -gen_F(sys, "y", "1H")
  
  # apply delay
  t = TE1 / 2
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
  t = (TE1 + TE2) / 2
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
  Fx <- gen_F(sys,"x","1H")
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
  Fy <- gen_F(sys,"y","1H")
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
  
  sys$rho <- gen_F(sys, "z","31P")
   
  angle <- 90
  Fx <- gen_F(sys,"x","31P")
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
  Fy <- gen_F(sys,"y","31P")
  lhs_pulse <- complexplus::matexp(-Fy * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fy * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # apply delay
  sys$rho <- lhs %*% sys$rho %*% rhs
  
  # acquire
  acquire(sys, detect = "31P")
}