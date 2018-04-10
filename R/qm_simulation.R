#' Create a spin system object for pulse sequence simulation
#' @param spin_params an object describing the spin system properties
#' @param ft transmitter frequency in Hz
#' @param ref reference value for ppm scale
#' @return spin system object
#' @export
spin_sys <- function(spin_params, ft, ref) {
  # TODO checks on input
  
  # force uppercase
  spin_params$nucleus <- toupper(spin_params$nucleus)
  
  spin_num <- get_spin_num(spin_params$nucleus)
  
  # calculate the Hamiltonian
  H_mat <- H(spin_num, spin_params$nucleus, spin_params$chem_shift, 
             spin_params$j_coupling_mat, ft, ref)
  
  # perform a symmetric eigenvalue decomposition
  res <- eigen(H_mat, symmetric = TRUE) 
  
  list(spin_num = spin_num, nucleus = spin_params$nucleus,
       H_mat = H_mat, H_eig_vals = res$values, H_eig_vecs = res$vectors)
}

H <- function(spin_n, nucleus, chem_shift, j_coupling_mat, ft, ref) {
  basis_size <- prod(spin_n * 2 + 1)
  H_mat <- matrix(0, basis_size, basis_size)
  
  # chemical shift part
  for (n in (1:length(spin_n))) {
    # Convert chem shift to angular freq and apply to Iz
    H_mat <- H_mat + gen_I(n, spin_n, "z") * 
             ((-chem_shift[n] + ref) * ft * 1e-6)
  }
  
  # Find non-zero elements of j_coupling_mat
  inds <- which((j_coupling_mat != 0), arr.ind = TRUE)
  if ((dim(inds)[1]) > 0)  {
    # j-coupling part
    for (n in 1:dim(inds)[1]) {
      j <- j_coupling_mat[inds[n,1],inds[n,2]]
      H_mat <- H_mat + j * gen_I(inds[n,1], spin_n, "z") %*%
               gen_I(inds[n,2], spin_n, "z")
     
      # strong coupling for homonuclear spins
      if ( nucleus[inds[n,1]] == nucleus[inds[n,2]] ) {
        H_mat <- H_mat + j * gen_I(inds[n,1], spin_n, "x") %*%
                 gen_I(inds[n,2], spin_n, "x")
        
        H_mat <- H_mat + j * gen_I(inds[n,1], spin_n, "y") %*%
                 gen_I(inds[n,2], spin_n, "y")
      }
    }
  }
  H_mat
}

#' Simulate pulse sequence acquisition.
#' @param sys spin system object
#' @param rec_phase reciever phase in degrees
#' @param tol ignore resonance amplitudes below this threshold
#' @param detect detection nuclie
#' @return a list of resonance amplitudes and frequencies
#' @export
acquire <- function(sys, rec_phase = 180, tol = 1e-4, detect = NULL) {
  if (is.null(detect)) {
    Fp <- gen_F(sys, "p")
  } else {
    Fp <- gen_F(sys, "p", detect)
  }
    
  coherence <- Conj(t(sys$H_eig_vecs)) %*% sys$rho %*% sys$H_eig_vecs
  coupled_coherence <- Conj(t(sys$H_eig_vecs)) %*% Fp %*% sys$H_eig_vecs
  amp_mat <- coherence * coupled_coherence
  amp_scaling_factor <- 2i / nrow(coherence)
  # find resonances
  sig_amps <- (Mod(amp_mat) > tol)
  indx <- which(sig_amps, arr.ind = TRUE)
  amps <- amp_mat[indx] * (exp(1i * rec_phase * pi / 180) * amp_scaling_factor)
  freqs <- sys$H_eig_vals[indx[,1]] - sys$H_eig_vals[indx[,2]]
  list(amps = amps, freqs = freqs)
}

#' Generate the F product operator
#' @param sys spin system object
#' @param op operator, one of "x", "y", "z", "p", "m"
#' @param detect detection nuclei
#' @return F product operator matrix
#' @export
gen_F <- function(sys, op, detect=NULL) {
  basis_size <- prod(sys$spin_num * 2 + 1)
  F_mat <- matrix(0, basis_size, basis_size)
  if (is.null(detect)) {
    spin_indices <- 1:length(sys$spin_num)
  } else {
    spin_indices <- which(toupper(sys$nucleus) == toupper(detect))
  }
    
  for (n in spin_indices ) {
    F_mat = F_mat + gen_I(n, sys$spin_num, op)
  }
  F_mat
}

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
get_1h_brain_basis_paras <- function(ft, metab_lw = NULL, lcm_compat = FALSE) {
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

#' Simulate a basis-set suitable for 1H brain MRS analysis acquired with a PRESS 
#' sequence. Note, ideal pulses are assumed.
#' @param acq_paras List of acquisition parameters or an mrs_data object. See
#' \code{\link{def_acq_paras}}
#' @param xlim Range of frequencies to simulate in ppm.
#' @param lcm_compat Exclude lipid and MM signals for use with default LCModel
#' options.
#' @param TE1 TE1 of PRESS sequence (TE = TE1 + TE2).
#' @param TE2 TE2 of PRESS sequence.
#' @return Basis object.
#' @export
sim_basis_1h_brain_press <- function(acq_paras = def_acq_paras(),
                                     xlim = c(0.5, 4.2), lcm_compat = FALSE, 
                                     TE1 = 0.01, TE2 = 0.02) {
  
  if (class(acq_paras) == "mrs_data") {
    acq_paras <- get_acq_paras(acq_paras)
  }
  
  sim_basis(get_1h_brain_basis_paras(ft = acq_paras$ft, lcm_compat = lcm_compat), 
                                     seq_press_ideal, fs = acq_paras$fs, 
                                     N = acq_paras$N, ref = acq_paras$ref,
                                     ft = acq_paras$ft, xlim = xlim, TE1 = TE1,
                                     TE2 = TE2)
}

get_mol_para_list_names <- function(mol_para_list) {
  names <- rep(NA, length(mol_para_list))
  for (n in 1:length(mol_para_list)) {
    names[n] <- mol_para_list[[n]]$name
  }
  names 
}

#' Simulate a basis set object.
#' @param mol_list a list of \code{mol_parameter} objects.
#' @param pul_seq A pulse sequence function to use.
#' @param ft Transmitter frequency in Hz.
#' @param ref Reference value for ppm scale.
#' @param fs Sampling frequency in Hz.
#' @param N Number of data points in the spectral dimension.
#' @param xlim A ppm range limiting signals to be simulated.
#' @param ... Extra parameters to pass to the pulse sequence function.
#' @return A basis object.
#' @export
sim_basis <- function(mol_list, pul_seq = seq_pulse_acquire, ft = def_ft(),
                      ref = def_ref(), fs = def_fs(), N = def_N(),
                      xlim = NULL, ...) {
  
  basis_mrs_data <- sim_zeros(ft = ft, ref = ref, fs = fs, N = N,
                              dyns = length(mol_list))
  
  for (n in 1:length(mol_list)) {
    mrs_data <- sim_mol(mol_list[[n]], pul_seq, ft, ref, fs, N, xlim, ...) 
    basis_mrs_data <- set_dyns(basis_mrs_data,n,mrs_data)
  }
  names <- get_mol_para_list_names(mol_list)
  mrs_data2basis(basis_mrs_data, names = names)
}

#' Simulate a \code{mol_parameter} object.
#' @param mol a \code{mol_parameter} object.
#' @param pul_seq A pulse sequence function to use.
#' @param ft Transmitter frequency in Hz.
#' @param ref Reference value for ppm scale.
#' @param fs Sampling frequency in Hz.
#' @param N Number of data points in the spectral dimension.
#' @param xlim A ppm range limiting signals to be simulated.
#' @param ... Extra parameters to pass to the pulse sequence function.
#' @return An \code{mrs_data} object.
#' @export
sim_mol <- function(mol, pul_seq = seq_pulse_acquire, ft = def_ft(), 
                    ref = def_ref(), fs = def_fs(), N = def_N(),
                    xlim = NULL, ...) {
  # create empty fid
  mrs_data <- sim_zeros(fs = fs, N = N, ft = ft, ref = ref)
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
  
  # first pt correction - shouldn't need this because already done in 
  # sim_resonances_fast function
  # mrs_data$data[,,,,,,1] <- 0.5 * mrs_data$data[,,,,,,1]
  
  mrs_data
}