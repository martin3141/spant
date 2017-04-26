pulse_acquire <- function(spin_params, B0, offset) {
  sys <- spin_sys(spin_params, B0, offset)
  
  sys$rho <- gen_F(sys, "z")
   
  angle <- 90
  Fx <- gen_F(sys, "x", "1H")
  lhs_pulse <- complexplus::matexp(-Fx * 1i * angle * pi / 180)
  rhs_pulse <- complexplus::matexp(Fx * 1i * angle * pi / 180)
  sys$rho <- lhs_pulse %*% sys$rho %*% rhs_pulse
  
  # acquire
  acquire(sys, detect = "1H")
}
  
press_ideal <- function(spin_params, B0, offset, TE1, TE2) {
  
  sys <- spin_sys(spin_params, B0, offset)
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

spin_echo_ideal <- function(spin_params, B0, offset, TE) {
  sys <- spin_sys(spin_params, B0, offset)
  
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