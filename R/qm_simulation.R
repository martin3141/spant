#' Create a spin system object for pulse sequence simulation.
#' @param spin_params an object describing the spin system properties.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param precomp_jc_H use a precomputed J-coupling H matrix to save time.
#' @param precomp_Iz use precomputed Iz matrices to save time.
#' @return spin system object.
#' @export
spin_sys <- function(spin_params, ft, ref, precomp_jc_H = NULL,
                     precomp_Iz = NULL) {
  
  # force uppercase
  spin_params$nucleus <- toupper(spin_params$nucleus)
  
  spin_num <- get_spin_num(spin_params$nucleus)
  
  # calculate the Hamiltonian
  if (!is.null(precomp_jc_H)) {
    omit_jc <- TRUE
  } else {
    omit_jc <- FALSE
  }
      
  H_mat_list <- H(spin_num, spin_params$nucleus, spin_params$chem_shift, 
                  spin_params$j_coupling_mat, ft, ref, omit_jc, precomp_Iz)
  
  if (!is.null(precomp_jc_H)) {
    H_mat_list$H_mat_jc <- precomp_jc_H
    H_mat_list$H_mat    <- H_mat_list$H_mat + H_mat_list$H_mat_jc
  }
  
  # perform a symmetric eigenvalue decomposition
  H_mat_res <- eigen(H_mat_list$H_mat, symmetric = TRUE)
  
  res <- list(spin_num = spin_num, nucleus = spin_params$nucleus,
              H_mat = H_mat_list$H_mat, H_eig_vals = H_mat_res$values,
              H_eig_vecs = H_mat_res$vectors, H_mat_jc = H_mat_list$H_mat_jc)
  
  return(res)
}

H <- function(spin_n, nucleus, chem_shift, j_coupling_mat, ft, ref,
              omit_jc = FALSE, precomp_Iz = NULL) {
  
  basis_size <- prod(spin_n * 2 + 1)
  
  # chemical shift part
  H_mat_cs <- matrix(0, basis_size, basis_size)
  
  # J-coupling part
  H_mat_jc <- matrix(0, basis_size, basis_size)
  
  # Combined chemical shift and J-coupling
  H_mat <- matrix(0, basis_size, basis_size)
  
  # chemical shift part
  if (is.null(precomp_Iz)) {
    for (n in (1:length(spin_n))) {
      # Convert chem shift to angular freq and apply to Iz
      H_mat_cs <- H_mat_cs + gen_I(n, spin_n, "z") * 
                  ((chem_shift[n] - ref) * ft * 1e-6)
    }
  } else {
    for (n in (1:length(spin_n))) {
      # Convert chem shift to angular freq and apply to Iz
      H_mat_cs <- H_mat_cs + precomp_Iz[[n]] * 
                  ((chem_shift[n] - ref) * ft * 1e-6)
    }
  }
  
  if (omit_jc) {
    H_mat <- H_mat_cs
    return(list(H_mat = H_mat, H_mat_cs = H_mat_cs, H_mat_jc = H_mat_jc))
  }
  
  # Find non-zero elements of j_coupling_mat
  inds <- which((j_coupling_mat != 0), arr.ind = TRUE)
  
  if ((dim(inds)[1]) > 0)  {
    # j-coupling part
    for (n in 1:dim(inds)[1]) {
      j <- j_coupling_mat[inds[n, 1],inds[n, 2]]
      H_mat_jc <- H_mat_jc + j * gen_I(inds[n, 1], spin_n, "z") %*%
                  gen_I(inds[n, 2], spin_n, "z")
     
      # strong coupling for homonuclear spins
      if (nucleus[inds[n, 1]] == nucleus[inds[n, 2]]) {
        H_mat_jc <- H_mat_jc + j * gen_I(inds[n, 1], spin_n, "x") %*%
                    gen_I(inds[n, 2], spin_n, "x")
        
        H_mat_jc <- H_mat_jc + j * gen_I(inds[n, 1], spin_n, "y") %*%
                    gen_I(inds[n, 2], spin_n, "y")
      }
    }
  }
  
  H_mat <- H_mat_cs + H_mat_jc
  
  return(list(H_mat = H_mat, H_mat_cs = H_mat_cs, H_mat_jc = H_mat_jc))
}

#' Simulate pulse sequence acquisition.
#' @param sys spin system object.
#' @param rec_phase receiver phase in degrees.
#' @param tol ignore resonance amplitudes below this threshold.
#' @param detect detection nuclei.
#' @param amp_scale scaling factor for the output amplitudes.
#' @return a list of resonance amplitudes and frequencies.
#' @export
acquire <- function(sys, rec_phase = 0, tol = 1e-4, detect = NULL,
                    amp_scale = 1) {
  
  if (is.null(detect)) {
    Fm <- gen_F(sys, "m")
  } else {
    Fm <- gen_F(sys, "m", detect)
  }
    
  coherence <- Conj(t(sys$H_eig_vecs)) %*% sys$rho %*% sys$H_eig_vecs
  coupled_coherence <- Conj(t(sys$H_eig_vecs)) %*% Fm %*% sys$H_eig_vecs
  amp_mat <- coherence * coupled_coherence
  amp_scaling_factor <- 2i / nrow(coherence)
  
  # find resonances
  sig_amps <- (Mod(amp_mat) > tol)
  indx <- which(sig_amps, arr.ind = TRUE)
  amps <- amp_mat[indx] * (exp(1i * rec_phase * pi / 180) * amp_scaling_factor)
  freqs <- sys$H_eig_vals[indx[,1]] - sys$H_eig_vals[indx[,2]]
  list(amps = amps * amp_scale, freqs = freqs)
}

#' Generate the F product operator.
#' @param sys spin system object.
#' @param op operator, one of "x", "y", "z", "p", "m".
#' @param detect detection nuclei.
#' @return F product operator matrix.
#' @export
gen_F <- function(sys, op, detect = NULL) {
  basis_size <- prod(sys$spin_num * 2 + 1)
  F_mat <- matrix(0, basis_size, basis_size)
  if (is.null(detect)) {
    spin_indices <- 1:length(sys$spin_num)
  } else {
    spin_indices <- which(toupper(sys$nucleus) == toupper(detect))
  }
    
  for (n in spin_indices) {
    F_mat <- F_mat + gen_I(n, sys$spin_num, op)
  }
  F_mat
}

#' Generate the Fxy product operator with a specified phase.
#' @param sys spin system object.
#' @param phase phase angle in degrees.
#' @param detect detection nuclei.
#' @return product operator matrix.
#' @export
gen_F_xy <- function(sys, phase, detect = NULL) {
  F_mat <- cos(phase * pi / 180) * gen_F(sys, "x", detect) +
           sin(phase * pi / 180) * gen_F(sys, "y", detect)
  F_mat
}

#' Get the quantum coherence matrix for a spin system.
#' @param sys spin system object.
#' @return quantum coherence number matrix.
#' @export
qn_states <- function(sys) {
  Fz <- gen_F(sys, "z")
  states_vec <- diag(Fz)
  states_mat <- outer(states_vec, states_vec, '-')
  states_mat <- Re(states_mat)
  mode(states_mat) <- "integer"
  return(states_mat)
}

#' Zero all coherence orders other than the one supplied as an argument.
#' @param sys spin system object.
#' @param rho density matrix.
#' @param order coherence order to keep (default is 0).
#' @return density matrix.
#' @export
coherence_filter <- function(sys, rho, order = 0) {
  qn_states <- qn_states(sys)
  rho[qn_states != order] <- 0
  return(rho)
}

#' Zero all coherences including and above a given order.
#' @param sys spin system object.
#' @param rho density matrix.
#' @param order states higher than or equal to this argument will be set to
#' zero.
#' @return density matrix.
#' @export
zero_higher_orders <- function(sys, rho, order) {
  order <- as.integer(order)
  qn_states <- qn_states(sys)
  rho[Mod(qn_states) >= order] <- 0
  return(rho)
}

#' Generate the I product operator for a single spin.
#' @param n spin index number for the required operator.
#' @param spin_num vector of spin numbers in the system.
#' @param op operator, one of "x", "y", "z", "p", "m".
#' @return I product operator matrix.
#' @export
gen_I <- function(n, spin_num, op) {
  N <- length(spin_num)
  I <- spin_num[n]
  basis_size <- prod(spin_num * 2 + 1)
  
  if (op == "z") {
    op_mat <- Iz_pauli(I) 
  } else if (op == "x") {
    op_mat <- Ix_pauli(I) 
  } else if (op == "y") {
    op_mat <- Iy_pauli(I) 
  } else if (op == "p") {
    op_mat <- Ip_pauli(I) 
  } else if (op == "m") {
    op_mat <- Im_pauli(I) 
  } else {
    stop("Unrecognised operator requested.")
  }
  
  if (n == 1) {
    I_mat_size <- basis_size / (2 * I + 1)
    out <- op_mat %x% diag(I_mat_size)
  } else if (n == N) {
    I_mat_size <- basis_size / (2 * I + 1)
    out <- diag(I_mat_size) %x% op_mat
  } else {
    lsize <- prod(spin_num[1:(n - 1)] * 2 + 1)
    rsize <- prod(spin_num[(1 + n):N] * 2 + 1)
    out <- diag(lsize) %x% op_mat %x% diag(rsize)
  }
  out
}

Iz_pauli <- function(I) {
  size <- 2 * I + 1
  diag(I - seq(0, size - 1))
}

Ip_pauli <- function(I) {
  size <- 2 * I + 1
  mat <- matrix(0, size, size)
  for (n in seq(0, size - 2)) {
    mat[n + 1, n + 2] <- (I * (I + 1) - (I - n - 1) * (I - n)) ^ 0.5
  }
  mat
}

Im_pauli <- function(I) {
  size <- 2 * I + 1
  mat <- matrix(0, size, size)
  for (n in seq(0, size - 2)) {
    mat[n + 2, n + 1] <- (I * (I + 1) - (I - n - 1) * (I - n)) ^ 0.5
  }
  mat
}

Ix_pauli <- function(I) {
  0.5 * (Ip_pauli(I) + Im_pauli(I))
}

Iy_pauli <- function(I) {
  -0.5i * (Ip_pauli(I) - Im_pauli(I))
}

#' Simulate an RF pulse on a single spin.
#' @param sys spin system object.
#' @param rho density matrix.
#' @param spin_n spin index.
#' @param angle RF flip angle in degrees.
#' @param nuc nucleus influenced by the pulse.
#' @param xy x or y pulse.
#' @return density matrix.
#' @export
apply_pulse <- function(sys, rho, spin_n, angle, nuc, xy) {
  
  if (nuc != sys$nucleus[spin_n]) return(rho)
  
  if (is.na(angle)) return(rho)
  
  if (angle < 1) return(rho)
  
  F <- gen_I(spin_n, sys$spin_num, xy)
  
  lhs_pulse <- matexp(-F * 1i * angle * pi / 180)
  rhs_pulse <- matexp( F * 1i * angle * pi / 180)
  
  rho_out <- lhs_pulse %*% rho %*% rhs_pulse
  
  return(rho_out)
}

#' Return the spin number for a given nucleus.
#' @param nucleus nucleus name, eg "1H".
#' @return spin number.
#' @export
get_spin_num <- function(nucleus) {
  spin_lookup <- data.frame(nucleus = c("1H", "31P", "14N", "13C"),
                            spin = c(0.5, 0.5, 1.0, 0.5))
  matches <- match(toupper(nucleus), spin_lookup$nucleus)
  spin_lookup$spin[matches]
}

#' Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
#' analyses.
#' @param ft transmitter frequency in Hz.
#' @param metab_lw linewidth of metabolite signals (Hz).
#' @param lcm_compat when TRUE, lipid, MM and -CrCH molecules will be excluded
#' from the output.
#' @return list of \code{mol_parameter} objects.
#' @export
get_1h_brain_basis_paras_v1 <- function(ft, metab_lw = NULL,
                                        lcm_compat = FALSE) {
  if (!lcm_compat) {
    m_cr_ch2 <- get_m_cr_ch2_paras(metab_lw)
  }
  ala <- get_ala_paras(metab_lw)
  asp <- get_asp_paras(metab_lw)
  cr <- get_cr_paras(metab_lw)
  gaba <- get_gaba_paras(metab_lw)
  glc <- get_a_glc_paras(metab_lw)
  gln <- get_gln_paras(metab_lw)
  gsh <- get_gsh_paras(metab_lw)
  glu <- get_glu_paras(metab_lw)
  gpc <- get_gpc_paras(metab_lw)
  ins <- get_ins_paras(metab_lw)
  lac <- get_lac_paras(metab_lw)
  if (!lcm_compat) {
    lip09 <- get_lip09_paras(ft)
    lip13a <- get_lip13a_paras(ft)
    lip13b <- get_lip13b_paras(ft)
    lip20 <- get_lip20_paras(ft)
    mm09 <- get_mm09_paras(ft)
    mm12 <- get_mm12_paras(ft)
    mm14 <- get_mm14_paras(ft)
    mm17 <- get_mm17_paras(ft)
    mm20 <- get_mm20_paras(ft)
  }
  naa <- get_naa_paras(metab_lw)
  naag <- get_naag_ch3_paras(metab_lw)
  pch <- get_pch_paras(metab_lw)
  pcr <- get_pcr_paras(metab_lw)
  sins <- get_sins_paras(metab_lw)
  tau <- get_tau_paras(metab_lw)
  
  if (!lcm_compat) {
    basis_list <- list(m_cr_ch2, ala, asp, cr, gaba, glc, gln, gsh, glu, gpc,
                       ins, lac, lip09, lip13a, lip13b, lip20, mm09, mm12, mm14,
                       mm17, mm20, naa, naag, pch, pcr, sins, tau)
  } else {
    basis_list <- list(ala, asp, cr, gaba, glc, gln, gsh, glu, gpc, ins, lac,
                       naa, naag, pch, pcr, sins, tau)
  }
  
  basis_list
}

#' Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
#' analyses.
#' @param ft transmitter frequency in Hz.
#' @param metab_lw linewidth of metabolite signals (Hz).
#' @param lcm_compat when TRUE, lipid, MM and -CrCH molecules will be excluded
#' from the output.
#' @return list of \code{mol_parameter} objects.
#' @export
get_1h_brain_basis_paras_v2 <- function(ft, metab_lw = NULL, lcm_compat = FALSE) {
  if (!lcm_compat) {
    m_cr_ch2 <- get_m_cr_ch2_paras(metab_lw)
  }
  ala <- get_ala_paras(metab_lw)
  asp <- get_asp_paras(metab_lw)
  cr <- get_cr_paras(metab_lw)
  gaba <- get_gaba_paras(metab_lw)
  glc <- get_glc_paras(metab_lw)
  gln <- get_gln_paras(metab_lw)
  gsh <- get_gsh_paras(metab_lw)
  glu <- get_glu_paras(metab_lw)
  gpc <- get_gpc_paras(metab_lw)
  ins <- get_ins_paras(metab_lw)
  lac <- get_lac_paras(metab_lw)
  if (!lcm_compat) {
    lip09 <- get_lip09_paras(ft)
    lip13a <- get_lip13a_paras(ft)
    lip13b <- get_lip13b_paras(ft)
    lip20 <- get_lip20_paras(ft)
    mm09 <- get_mm09_paras(ft)
    mm12 <- get_mm12_paras(ft)
    mm14 <- get_mm14_paras(ft)
    mm17 <- get_mm17_paras(ft)
    mm20 <- get_mm20_paras(ft)
  }
  naa <- get_naa_paras(metab_lw)
  naag <- get_naag_paras(metab_lw)
  pch <- get_pch_paras(metab_lw)
  pcr <- get_pcr_paras(metab_lw)
  sins <- get_sins_paras(metab_lw)
  tau <- get_tau_paras(metab_lw)
  
  if (!lcm_compat) {
    basis_list <- list(m_cr_ch2, ala, asp, cr, gaba, glc, gln, gsh, glu, gpc,
                       ins, lac, lip09, lip13a, lip13b, lip20, mm09, mm12, mm14,
                       mm17, mm20, naa, naag, pch, pcr, sins, tau)
  } else {
    basis_list <- list(ala, asp, cr, gaba, glc, gln, gsh, glu, gpc, ins, lac,
                       naa, naag, pch, pcr, sins, tau)
  }
  basis_list
}

#' Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
#' analyses.
#' @param ft transmitter frequency in Hz.
#' @param metab_lw linewidth of metabolite signals (Hz).
#' @param lcm_compat when TRUE, lipid, MM and -CrCH molecules will be excluded
#' from the output.
#' @return list of \code{mol_parameter} objects.
#' @export
get_1h_brain_basis_paras_v3 <- function(ft, metab_lw = NULL, lcm_compat = FALSE) {
  if (!lcm_compat) {
    m_cr_ch2 <- get_m_cr_ch2_paras(metab_lw)
  }
  ala <- get_ala_paras(metab_lw)
  asp <- get_asp_paras(metab_lw)
  cr <- get_cr_paras(metab_lw)
  gaba <- get_gaba_paras(metab_lw)
  glc <- get_glc_paras(metab_lw)
  gln <- get_gln_paras(metab_lw)
  gly <- get_gly_paras(metab_lw)
  gsh <- get_gsh_paras(metab_lw)
  glu <- get_glu_paras(metab_lw)
  gpc <- get_gpc_paras(metab_lw)
  ins <- get_ins_paras(metab_lw)
  lac <- get_lac_paras(metab_lw)
  if (!lcm_compat) {
    lip09 <- get_lip09_paras(ft)
    lip13a <- get_lip13a_paras(ft)
    lip13b <- get_lip13b_paras(ft)
    lip20 <- get_lip20_paras(ft)
    mm09 <- get_mm09_paras(ft)
    mm12 <- get_mm12_paras(ft)
    mm14 <- get_mm14_paras(ft)
    mm17 <- get_mm17_paras(ft)
    mm20 <- get_mm20_paras(ft)
  }
  naa <- get_naa_paras(metab_lw)
  naag <- get_naag_paras(metab_lw)
  pch <- get_pch_paras(metab_lw)
  pcr <- get_pcr_paras(metab_lw)
  peth <- get_peth_paras(metab_lw)
  sins <- get_sins_paras(metab_lw)
  tau <- get_tau_paras(metab_lw)
  
  if (!lcm_compat) {
    basis_list <- list(m_cr_ch2, ala, asp, cr, gaba, glc, gln, gly, gsh, glu,
                       gpc, ins, lac, lip09, lip13a, lip13b, lip20, mm09, mm12,
                       mm14, mm17, mm20, naa, naag, pch, pcr, peth, sins, tau)
  } else {
    basis_list <- list(ala, asp, cr, gaba, glc, gln, gly, gsh, glu, gpc, ins,
                       lac, naa, naag, pch, pcr, peth, sins, tau)
  }
  basis_list
}

#' Return a list of \code{mol_parameter} objects suitable for 1H brain MRS
#' analyses.
#' @param ft transmitter frequency in Hz.
#' @param metab_lw linewidth of metabolite signals (Hz).
#' @param lcm_compat when TRUE, lipid, MM and -CrCH molecules will be excluded
#' from the output.
#' @return list of \code{mol_parameter} objects.
#' @export
get_1h_brain_basis_paras <- function(ft, metab_lw = NULL, lcm_compat = FALSE) {
  get_1h_brain_basis_paras_v2(ft, metab_lw, lcm_compat)
}

#' Return a character vector of molecules included in the Gold Star Phantoms 
#' SPECTRE phantom.
#' @return a character vector of molecule names.
#' @export
get_1h_spectre_basis_names <- function() {
  names <- c("cho_rt", "cr_ch2_rt", "cr_ch3_rt", "gaba", "glu_rt", "ins_rt", 
             "lac_rt", "naa_rt")
  return(names)
}

#' Return a character vector of molecules included in the GE BRAINO phantom.
#' @return a character vector of molecule names.
#' @export
get_1h_braino_basis_names <- function() {
  names <- c("cho_rt", "cr_ch2_rt", "cr_ch3_rt", "glu_rt", "ins_rt", "lac_rt",
             "naa_rt")
  return(names)
}

#' Return a character vector of common 1H molecules found in healthy human
#' brain.
#' 
#' Note, this is a basic set and it may be appropriate to also include Asc, 
#' Gly and PEth for high quality MRS data.
#' 
#' @param add optional character vector of additional molecular names. Eg
#' c("asc", "gly", "peth").
#' @param remove optional character vector of molecular names to remove from the
#' set. Eg c("m_cr_ch2").
#' @param inc_lip_mm include Lipid and MM basis signals.
#' @return a character vector of molecule names.
#' @export
get_1h_brain_basis_names <- function(add = NULL, remove = NULL,
                                     inc_lip_mm = TRUE) {
  
  names <- c("m_cr_ch2", "ala", "asp", "cr", "gaba", "glc", "gln", "gsh", "glu",
             "gpc", "ins", "lac", "naa", "naag", "pch", "pcr", "sins", "tau")
  
  if (!is.null(add)) names <- sort(c(names, add))
  
  if (inc_lip_mm) {
    names <- c(names, "lip09", "lip13a", "lip13b", "lip20", "mm09", "mm12",
               "mm14", "mm17", "mm20")
  }
  
  if (!is.null(remove)) names <- names[!(names %in% remove)]
  
  return(names)
}

#' Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS 
#' sequence. Note, ideal pulses are assumed.
#' @param pul_seq pulse sequence function to use.
#' @param acq_paras list of acquisition parameters or an mrs_data object. See
#' \code{\link{def_acq_paras}}.
#' @param xlim range of frequencies to simulate in ppm.
#' @param lcm_compat exclude lipid and MM signals for use with default LCModel
#' options.
#' @param ... extra parameters to pass to the pulse sequence function.
#' @return basis object.
#' @export
sim_basis_1h_brain <- function(pul_seq = seq_press_ideal, 
                               acq_paras = def_acq_paras(), xlim = c(0.5, 4.2), 
                               lcm_compat = FALSE, ...) {
  
  sim_basis(get_1h_brain_basis_paras(ft = acq_paras$ft, lcm_compat = lcm_compat), 
                                     pul_seq = pul_seq, acq_paras = acq_paras,
                                     xlim = xlim, ...)
}

#' Simulate a macromolecular and lipid basis-set suitable for 1H brain MRS
#' analysis.
#' 
#' @param acq_paras list of acquisition parameters or an mrs_data object. See
#' \code{\link{def_acq_paras}}
#' @return basis object.
#' @export
sim_basis_mm_lip_lcm <- function(acq_paras = def_acq_paras()){
  
  ft <- acq_paras$ft 
  
  mol_list <- list(
    get_lip09_paras(ft),
    get_lip13a_paras(ft),
    get_lip13b_paras(ft),
    get_lip20_paras(ft),
    get_mm09_paras(ft),
    get_mm12_paras(ft),
    get_mm14_paras(ft),
    get_mm17_paras(ft),
    get_mm20_paras(ft)
  )
  
  return(sim_basis(mol_list, acq_paras = acq_paras))
}

#' Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS 
#' sequence. Note, ideal pulses are assumed.
#' @param acq_paras list of acquisition parameters or an mrs_data object. See
#' \code{\link{def_acq_paras}}
#' @param xlim range of frequencies to simulate in ppm.
#' @param lcm_compat exclude lipid and MM signals for use with default LCModel
#' options.
#' @param TE1 TE1 of PRESS sequence (TE = TE1 + TE2).
#' @param TE2 TE2 of PRESS sequence.
#' @return basis object.
#' @export
sim_basis_1h_brain_press <- function(acq_paras = def_acq_paras(),
                                     xlim = c(0.5, 4.2), lcm_compat = FALSE, 
                                     TE1 = 0.01, TE2 = 0.02) {
  
  sim_basis(get_1h_brain_basis_paras(ft = acq_paras$ft, lcm_compat = lcm_compat), 
                                     seq_press_ideal, acq_paras = acq_paras,
                                     xlim = xlim, TE1 = TE1, TE2 = TE2)
}

get_mol_para_list_names <- function(mol_para_list) {
  names <- rep(NA, length(mol_para_list))
  for (n in 1:length(mol_para_list)) {
    names[n] <- mol_para_list[[n]]$name
  }
  names 
}

#' Simulate a basis set object.
#' @param mol_list list of \code{mol_parameter} objects. Alternatively, a 
#' character vector matching molecules may also be provided. Use the 
#' get_mol_names function for a full list of molecules.
#' @param pul_seq pulse sequence function to use.
#' @param acq_paras list of acquisition parameters or an mrs_data object. See
#' \code{\link{def_acq_paras}}
#' @param xlim ppm range limiting signals to be simulated.
#' @param auto_scale scale the basis based on the intensity of a singlet
#' resonance. Needed for sequences with spatial simulation.
#' @param use_basis_cache create and use a cache of simulated basis sets stored
#' in the "spant_basis_cache" folder in the users home directory. Defaults to
#' FALSE.
#' @param verbose output simulation progress and timings.
#' @param ... extra parameters to pass to the pulse sequence function.
#' @return basis object.
#' @export
sim_basis <- function(mol_list, pul_seq = seq_pulse_acquire,
                      acq_paras = def_acq_paras(), xlim = NULL,
                      auto_scale = FALSE, use_basis_cache = FALSE,
                      verbose = FALSE, ...) {
  
  if (inherits(acq_paras, "mrs_data")) acq_paras <- get_acq_paras(acq_paras)
  
  ft  <- acq_paras$ft
  ref <- acq_paras$ref
  fs  <- acq_paras$fs
  N   <- acq_paras$N
  
  if (use_basis_cache) {
    # make sure the cache directory exists
    basis_cache_path <- file.path(path.expand('~'), "spant_basis_cache")
    dir.create(basis_cache_path, showWarnings = FALSE)
    hash_obj <- list(pul_seq = as.character(substitute(pul_seq)), 
                     ft = round(ft, 1), ref = ref, fs = round(fs, 1), N = N,
                     xlim = xlim, auto_scale = auto_scale, ...,
                     mol_list = mol_list)
    
    hash <- digest::digest(hash_obj)
    
    if (verbose) cat("Basis hash :", hash, "\n")
    
    hash_base_path <- file.path(basis_cache_path, hash)
    basis_path     <- file.path(hash_base_path, "basis.rds")
    
    if (file.exists(basis_path)) {
      cat("Using precomputated basis :", basis_path, "\n")
      basis <- readRDS(basis_path) 
      return(basis)
    } else {
      dir.create(hash_base_path)
      if (verbose) cat("Precomputated basis not found, calculating...\n")
      cat("Cached basis not found, simulation may take a few minutes.\n")
      cat("Subseqent requests with matched parameters will use precomputed\n")
      cat("results to save computation time.\n")
    }
  }
  
  if (inherits(mol_list[[1]], "character")) {
    if (length(mol_list) == 1) {
      mol_list <- list(get_mol_paras(mol_list, ft = ft))
    } else {
      mol_list <- get_mol_paras(mol_list, ft = ft)
    }
  }
   
  basis_mrs_data <- sim_zero(ft = ft, ref = ref, fs = fs, N = N,
                             dyns = length(mol_list))
  
  if (verbose) {
    cat("Simulation started.\n")
    start_time_full <- Sys.time()
  }
  
  for (n in 1:length(mol_list)) {
    if (verbose) {
      cat(paste0("Simuating ", n, " of ", length(mol_list), " : ",
                 mol_list[[n]]$full_name, "\n"))
      start_time <- Sys.time()
    }
    mrs_data <- sim_mol(mol_list[[n]], pul_seq, ft, ref, fs, N, xlim, ...) 
    basis_mrs_data <- set_dyns(basis_mrs_data, n, mrs_data)
    if (verbose) {
      end_time <- Sys.time()
      print(round(end_time - start_time, 2))
    }
  }
  
  if (verbose) {
    cat("Simulation finished.\n")
    end_time_full <- Sys.time()
    print(round(end_time_full - start_time_full, 2))
  }
  
  names <- get_mol_para_list_names(mol_list)
  basis <- mrs_data2basis(basis_mrs_data, names = names)
 
  # auto scale the basis
  if (auto_scale) {
    if (verbose) cat("Auto scaling the basis set.\n")
    # simulate a singlet for scaling basis sets with CSD
    protons     <- 1
    singlet     <- get_uncoupled_mol("scaling_singlet", 2, "1H", 1, 2, 0)
    mrs_singlet <- sim_mol(singlet, pul_seq, ft, ref, fs, N, xlim, ...) 
    sc_factor   <- get_td_amp(mrs_singlet, nstart = 2) / protons * 2
    if (verbose) cat("Rescale value :", sc_factor, "\n")
    basis$data  <- basis$data / as.numeric(sc_factor)
  }
  
  if (use_basis_cache) {
    # if we got here then we need to store the basis for future use
    if (verbose) cat("Saving precomputated basis :", basis_path, "\n")
    saveRDS(basis, basis_path)
    jsonlite::write_json(hash_obj, file.path(hash_base_path,
                                             "basis_paras.json"), pretty = TRUE,
                         auto_unbox = TRUE, force = TRUE)
  }
    
  return(basis)
}

#' Simulate a \code{mol_parameter} object.
#' @param mol \code{mol_parameter} object.
#' @param pul_seq pulse sequence function to use.
#' @param ft transmitter frequency in Hz.
#' @param ref reference value for ppm scale.
#' @param fs sampling frequency in Hz.
#' @param N number of data points in the spectral dimension.
#' @param xlim ppm range limiting signals to be simulated.
#' @param ... extra parameters to pass to the pulse sequence function.
#' @return \code{mrs_data} object.
#' @export
sim_mol <- function(mol, pul_seq = seq_pulse_acquire, ft = def_ft(), 
                    ref = def_ref(), fs = def_fs(), N = def_N(),
                    xlim = NULL, ...) {
  # create empty fid
  mrs_data <- sim_zero(fs = fs, N = N, ft = ft, ref = ref)
  for (group in mol$spin_groups) {
    res <- pul_seq(group, ft, ref, ...)
    
    if (!is.null(xlim)) {
      # filter frequences based on xlim
      xlim_hz <- sort(ppm2hz(xlim, ft = ft, ref = ref))
      inds <- c(which(res$freqs < xlim_hz[1]), which(res$freqs > xlim_hz[2]))
      if ( length(inds) > 0 ) {
        res$amps <- res$amps[-inds]
        res$freqs <- res$freqs[-inds]
      }
    }
    
    if (length(res$freqs > 0)) {
      group_data <- sim_resonances_fast(res$freqs, res$amps, freq_ppm = FALSE, 
                                        ft = ft, fs = fs, ref = ref, N = N)
      group_data <- lb(group_data,group$lw, group$lg)
      mrs_data <- mrs_data + (group_data * group$scale_factor)
    }
  }
  
  return(mrs_data)
}