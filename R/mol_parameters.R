#' Get a \code{mol_parameters} object for a named molecule.
#' @param name the name of the molecule.
#' @param ... arguments to pass to molecule definition function.
#' @export
get_mol_paras <- function(name, ...) {
  get(paste("get_", tolower(name), "_paras", sep = ""))(...)
}

#' Return a character array of names that may be used with the 
#' \code{get_mol_paras} function.
#' @return a character array of names.
#' @export
get_mol_names <- function() {
  funs <- ls(getNamespace("spant"), all.names = TRUE)
  funs <- funs[!funs %in% c("get_mol_paras", "get_acq_paras",
                            "get_1h_brain_basis_paras")]
  
  sub("_paras", "", sub("get_", "", funs[grep("get_.*_paras", funs)]))
}

#' Generate a \code{mol_parameters} object for a simple spin system with one resonance.
#' @param name name of the molecule.
#' @param chem_shift chemical shift of the resonance (PPM).
#' @param nucleus nucleus (1H, 31P...).
#' @param scale_factor multiplicative scaling factor.
#' @param lw linewidth in Hz.
#' @param lg Lorentz-Gauss lineshape parameter (between 0 and 1).
#' @return mol_parameters object.
#' @export
get_uncoupled_mol <- function(name, chem_shift, nucleus, scale_factor, lw, lg) {
  j_coupling_mat <- matrix(0, 1, 1)
  spin_groups <- vector("list", length(chem_shift))
  for (n in 1:length(chem_shift)) {
    spin_group <- list(nucleus = nucleus[n], chem_shift = chem_shift[n], 
                       j_coupling_mat = j_coupling_mat,
                       scale_factor = scale_factor[n], lw = lw[n], lg = lg[n])
    spin_groups[[n]] <- spin_group
  }
  paras <- list(spin_groups = spin_groups, name = name)
  class(paras) <- "mol_parameters"
  paras
}

#' @export
print.mol_parameters <- function(x, ...) {
  cat(paste(c("Name        : ", x$name, "\n")), sep = "")
  cat(paste(c("Source      : ", x$source, "\n")), sep = "")
  cat(paste(c("Spin groups : ", length(x$spin_groups), "\n")), sep = "")
  for (n in 1:length(x$spin_groups)) {
    cat("\n")
    cat(paste(c("Spin group ", n, "\n")), sep = "")
    cat("------------\n")
    print(data.frame(nucleus = x$spin_groups[[n]]$nucleus,
                     chem_shift = x$spin_groups[[n]]$chem_shift))
    
    if (nrow(x$spin_groups[[n]]$j_coupling_mat) > 1) {
      cat("\n")
      cat("j-coupling matrix\n")
      rownames(x$spin_groups[[n]]$j_coupling_mat) <- x$spin_groups[[n]]$chem_shift
      colnames(x$spin_groups[[n]]$j_coupling_mat) <- x$spin_groups[[n]]$chem_shift
      j_mat <- x$spin_groups[[n]]$j_coupling_ma
      j_mat[upper.tri(j_mat, diag = TRUE)] <- NA
      j_mat[j_mat == 0] <- NA
      #j_mat <- t(j_mat)t
      print(j_mat, na.print = "-")
    }
  }
}

get_m_cr_ch2_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("-CrCH2", 3.913, "1H", -2, lw, lg)
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  paras$source = source
  paras
}

get_ala_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 4)
  chem_shift <- c(3.7746, 1.4667, 1.4667, 1.4667)
  j_coupling_mat <- matrix(0, 4, 4)
  j_coupling_mat[2,1] <- 7.234
  j_coupling_mat[3,1] <- 7.234
  j_coupling_mat[4,1] <- 7.234
  j_coupling_mat[3,2] <- -14.366
  j_coupling_mat[4,2] <- -14.366
  j_coupling_mat[4,3] <- -14.366
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13: 129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Ala", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_asp_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 3)
  chem_shift <- c(3.8914, 2.8011, 2.6533)
  j_coupling_mat <- matrix(0, 3, 3)
  j_coupling_mat[2,1] <- 3.647
  j_coupling_mat[3,1] <- 9.107
  j_coupling_mat[3,2] <- -17.426
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Asp", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_cr_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("Cr", 3.027, "1H", 3, lw, lg)
  paras_b <- get_uncoupled_mol("Cr", 3.913, "1H", 2, lw, lg)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras$source = source
  paras
}

get_gaba_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 6)
  chem_shift <- c(3.0128, 3.0128, 1.889, 1.889, 2.284, 2.284)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- -12.021
  j_coupling_mat[3,1] <- 5.372
  j_coupling_mat[4,1] <- 7.127
  j_coupling_mat[3,2] <- 10.578
  j_coupling_mat[4,2] <- 6.982
  j_coupling_mat[4,3] <- -13.121
  j_coupling_mat[5,3] <- 7.755
  j_coupling_mat[6,3] <- 7.432
  j_coupling_mat[5,4] <- 6.173
  j_coupling_mat[6,4] <- 7.933
  j_coupling_mat[6,5] <- -10.744
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13: 129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "GABA", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_gaba_jn_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 6)
  chem_shift <- c(2.2840, 2.2840, 1.8880, 1.8880, 3.0130, 3.0130)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- -15.938
  j_coupling_mat[3,1] <- 7.678
  j_coupling_mat[4,1] <- 6.980
  j_coupling_mat[3,2] <- 6.980
  j_coupling_mat[4,2] <- 7.678
  j_coupling_mat[4,3] <- -15.000
  j_coupling_mat[5,3] <- 8.510
  j_coupling_mat[6,3] <- 6.503
  j_coupling_mat[5,4] <- 6.503
  j_coupling_mat[6,4] <- 8.510
  j_coupling_mat[6,5] <- -14.062
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "J-Difference Editing of Gamma-Aminobutyric Acid (GABA): 
              Simulated and Experimental Multiplet Patterns. MRM 2013; 
              70:1183-1191."
  
  paras <- list(spin_groups = list(spin_group_a), name = "GABA", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_gln_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 5)
  chem_shift <- c(3.753, 2.129, 2.109, 2.432, 2.454)
  j_coupling_mat <- matrix(0, 5, 5)
  j_coupling_mat[2,1] <- 5.847
  j_coupling_mat[3,1] <- 6.5
  j_coupling_mat[3,2] <- -14.504
  j_coupling_mat[4,2] <- 9.165
  j_coupling_mat[5,2] <- 6.347
  j_coupling_mat[4,3] <- 6.324
  j_coupling_mat[5,3] <- 9.209
  j_coupling_mat[5,4] <- -15.371
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Gln", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_gsh_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- c("1H")
  chem_shift <- c(3.769)
  j_coupling_mat <- matrix(0, 1, 1)
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 2,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H")
  chem_shift <- c(4.5608, 2.9264, 2.9747)
  j_coupling_mat <- matrix(0, 3, 3)
  j_coupling_mat[2,1] <- 7.09
  j_coupling_mat[3,1] <- 4.71
  j_coupling_mat[3,2] <- -14.06
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H", "1H", "1H")
  chem_shift <- c(3.769, 2.159, 2.146, 2.510, 2.560)
  j_coupling_mat <- matrix(0, 5, 5)
  j_coupling_mat[2,1] <- 6.34
  j_coupling_mat[3,1] <- 6.36
  j_coupling_mat[3,2] <- -15.48
  j_coupling_mat[4,2] <- 6.7
  j_coupling_mat[5,2] <- 7.6
  j_coupling_mat[4,3] <- 7.6
  j_coupling_mat[5,3] <- 6.7
  j_coupling_mat[5,4] <- -15.92
  spin_group_c <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b, spin_group_c),
                name = "GSH", source = source)
  
  class(paras) <- "mol_parameters"
  paras
}

get_gly_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("Gly", 3.548, "1H", 2, lw, lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras$source <- source
  paras
}

get_ins_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 0.5
  nucleus <- rep("1H", 6)
  chem_shift <- c(3.5217, 4.0538, 3.5217, 3.6144, 3.269, 3.6144)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- 2.889
  j_coupling_mat[6,1] <- 9.998
  j_coupling_mat[3,2] <- 3.006
  j_coupling_mat[4,3] <- 9.997
  j_coupling_mat[5,4] <- 9.485
  j_coupling_mat[6,5] <- 9.482
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Ins", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_lip09_paras <- function(ft) {
  paras <- get_uncoupled_mol("Lip09", 0.89, "1H", 3, 0.14 * ft / 1e6, 1)
  paras$source <- "LCModel manual."
  paras
}

get_lip13a_paras <- function(ft) {
  paras <- get_uncoupled_mol("Lip13a", 1.28, "1H", 2, 0.15 * ft / 1e6, 1)
  paras$source <- "LCModel manual."
  paras
}

get_lip13b_paras <- function(ft) {
  paras <- get_uncoupled_mol("Lip13b", 1.28, "1H", 2, 0.089 * ft / 1e6, 1)
  paras$source <- "LCModel manual."
  paras
}

get_lip20_paras <- function(ft) {
  paras <- get_uncoupled_mol("Lip20", 2.04, "1H", 1.33, 0.15 * ft / 1e6, 1)
  paras_b <- get_uncoupled_mol("Lip20", 2.25, "1H", 0.67, 0.15 * ft / 1e6, 1)
  paras_c <- get_uncoupled_mol("Lip20", 2.8, "1H", 0.87, 0.2 * ft / 1e6, 1)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  paras$spin_groups[[3]] <- paras_c$spin_groups[[1]]
  paras$source <- "LCModel manual."
  paras
}

get_mm09_paras <- function(ft) {
  paras <- get_uncoupled_mol("MM09", 0.91, "1H", 3, 0.14 * ft / 1e6, 1)
  paras$source <- "LCModel manual."
  paras
}

get_mm12_paras <- function(ft) {
  paras <- get_uncoupled_mol("MM12", 1.21, "1H", 2, 0.15 * ft / 1e6, 1)
  paras$source <- "LCModel manual."
  paras
}

get_mm14_paras <- function(ft) {
  paras <- get_uncoupled_mol("MM14", 1.43, "1H", 2, 0.17 * ft / 1e6, 1)
  paras$source <- "LCModel manual."
  paras
}

get_mm17_paras <- function(ft) {
  paras <- get_uncoupled_mol("MM17", 1.67, "1H", 2, 0.15 * ft / 1e6, 1)
  paras$source <- "LCModel manual."
  paras
}

get_mm20_paras <- function(ft) {
  paras <- get_uncoupled_mol("MM20", 2.08, "1H", 1.33, 0.15 * ft / 1e6, 1)
  paras_b <- get_uncoupled_mol("MM20", 2.25, "1H", 0.33, 0.2 * ft / 1e6, 1)
  paras_c <- get_uncoupled_mol("MM20", 1.95, "1H", 0.33, 0.15 * ft / 1e6, 1)
  paras_d <- get_uncoupled_mol("MM20", 3.0, "1H", 0.4, 0.2 * ft / 1e6, 1)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  paras$spin_groups[[3]] <- paras_c$spin_groups[[1]]
  paras$spin_groups[[4]] <- paras_d$spin_groups[[1]]
  paras$source <- "LCModel manual."
  paras
}

get_naa_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 0.5
  nucleus <- c("1H")
  chem_shift <- c(2.008)
  j_coupling_mat <- matrix(0, 1, 1)
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 3,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H")
  chem_shift <- c(4.3817, 2.6727, 2.4863)
  j_coupling_mat <- matrix(0, 3, 3)
  j_coupling_mat[2,1] <- 3.861
  j_coupling_mat[3,1] <- 9.821
  j_coupling_mat[3,2] <- -15.592
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a,spin_group_b), name = "NAA",
                source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_naag_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 0.5
  paras <- get_uncoupled_mol("NAAG", 2.042, "1H", 3, lw, lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153. Note, only the acetyl moiety is simualted."
  
  paras$source <- source
  paras
}

get_pch_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- c("1H")
  chem_shift <- c(3.208)
  j_coupling_mat <- matrix(0, 1, 1)
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 9,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H", "1H", "14N", "31P")
  chem_shift <- c(4.2805, 4.2805, 3.641, 3.641, 0, 0)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- -14.89
  j_coupling_mat[3,1] <- 2.284
  j_coupling_mat[4,1] <- 7.231
  j_coupling_mat[5,1] <- 2.68
  j_coupling_mat[6,1] <- 6.298
  j_coupling_mat[3,2] <- 7.326
  j_coupling_mat[4,2] <- 2.235
  j_coupling_mat[5,2] <- 2.772
  j_coupling_mat[6,2] <- 6.249
  j_coupling_mat[4,3] <- -14.19
  
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b), name = "PCh",
                source = source)
  
  class(paras) <- "mol_parameters"
  paras
}

get_pcr_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("PCr", 3.029, "1H", 3, lw, lg)
  paras_b <- get_uncoupled_mol("PCr", 3.930, "1H", 2, lw, lg)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras$source <- source
  paras
}

get_sins_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("sIns", 3.34, "1H", 6, lw, lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras$source <- source
  paras
}

get_suc_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("Suc", 2.3920, "1H", 4, lw, lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras$source <- source
  paras
}

get_tau_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 4)
  chem_shift <- c(3.4206, 3.4206, 3.2459, 3.2459)
  j_coupling_mat <- matrix(0, 4, 4)
  j_coupling_mat[2,1] <- -12.438
  j_coupling_mat[3,1] <- 6.742
  j_coupling_mat[4,1] <- 6.464
  j_coupling_mat[3,2] <- 6.403
  j_coupling_mat[4,2] <- 6.792
  j_coupling_mat[4,3] <- -12.93
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Tau", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_lac_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 4)
  chem_shift <- c(4.0974, 1.3142, 1.3142, 1.3142)
  j_coupling_mat <- matrix(0, 4, 4)
  j_coupling_mat[2,1] <- 6.933
  j_coupling_mat[3,1] <- 6.933
  j_coupling_mat[4,1] <- 6.933
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Lac", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_glu_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 5)
  chem_shift <- c(3.7433, 2.0375, 2.1200, 2.3378, 2.352)
  j_coupling_mat <- matrix(0,5,5)
  j_coupling_mat[2,1] <- 7.331 
  j_coupling_mat[3,1] <- 4.651
  j_coupling_mat[3,2] <- -14.849
  j_coupling_mat[4,2] <- 6.413
  j_coupling_mat[5,2] <- 8.406
  j_coupling_mat[4,3] <- 8.478
  j_coupling_mat[5,3] <- 6.875
  j_coupling_mat[5,4] <- -15.915
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Glu", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_glc_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 7)
  chem_shift <- c(5.216, 3.519, 3.698, 3.395, 3.822, 3.826, 3.749)
  j_coupling_mat <- matrix(0, 7, 7)
  j_coupling_mat[2,1] <- 3.8
  j_coupling_mat[3,2] <- 9.6
  j_coupling_mat[4,3] <- 9.4
  j_coupling_mat[5,4] <- 9.9
  j_coupling_mat[6,5] <- 1.5
  j_coupling_mat[7,5] <- 6
  j_coupling_mat[7,6] <- -12.1
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153. Only the alpha-anomer is simulated as the 
              beta-anomer is eliminated by water supression according to the above reference."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Glc", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_gpc_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- c("1H")
  chem_shift <- c(3.212)
  j_coupling_mat <- matrix(0, 1, 1)
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 9,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H", "1H", "1H", "31P")
  chem_shift <- c(3.605, 3.672, 3.903, 3.871, 3.946, 0)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- -14.78
  j_coupling_mat[3,1] <- 5.77
  j_coupling_mat[3,2] <- 4.53
  j_coupling_mat[4,2] <- 5.77
  j_coupling_mat[5,2] <- 4.53
  j_coupling_mat[5,4] <- -14.78
  j_coupling_mat[6,4] <- 6.03
  j_coupling_mat[6,5] <- 6.03
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H", "1H", "14N", "31P")
  chem_shift <- c(4.312, 4.312, 3.659, 3.659, 0, 0)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- -9.32
  j_coupling_mat[3,1] <- 3.1
  j_coupling_mat[4,1] <- 5.9
  j_coupling_mat[5,1] <- 2.67
  j_coupling_mat[6,1] <- 6.03
  j_coupling_mat[3,2] <- 5.9
  j_coupling_mat[4,2] <- 3.1
  j_coupling_mat[5,2] <- 2.67
  j_coupling_mat[6,2] <- 6.03
  j_coupling_mat[4,3] <- -9.32
  spin_group_c <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for brain metabolites.
              NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b, spin_group_c),
                name = "GPC", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_10spin_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 10)
  chem_shift <- c(5.216, 3.519, 3.698, 3.395, 3.822, 3.826, 3.749, 1, 2, 3)
  j_coupling_mat <- matrix(0, 10, 10)
  j_coupling_mat[2,1] <- 3.8
  j_coupling_mat[3,2] <- 9.6
  j_coupling_mat[4,3] <- 9.4
  j_coupling_mat[5,4] <- 9.9
  j_coupling_mat[6,5] <- 1.5
  j_coupling_mat[7,5] <- 6
  j_coupling_mat[7,6] <- -12.1
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "made up molecule"
  
  paras <- list(spin_groups = list(spin_group_a), name = "Glc", source = source)
  class(paras) <- "mol_parameters"
  paras
}

get_9spin_paras <- function(lw = NULL, lg = 0) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 9)
  chem_shift <- c(5.216, 3.519, 3.698, 3.395, 3.822, 3.826, 3.749, 1, 2)
  j_coupling_mat <- matrix(0, 9, 9)
  j_coupling_mat[2,1] <- 3.8
  j_coupling_mat[3,2] <- 9.6
  j_coupling_mat[4,3] <- 9.4
  j_coupling_mat[5,4] <- 9.9
  j_coupling_mat[6,5] <- 1.5
  j_coupling_mat[7,5] <- 6
  j_coupling_mat[7,6] <- -12.1
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "made up molecule"
  
  paras <- list(spin_groups = list(spin_group_a), name = "Glc", source = source)
  class(paras) <- "mol_parameters"
  paras
}