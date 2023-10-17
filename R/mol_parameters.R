#' Get a \code{mol_parameters} object for a named molecule.
#' @param name the name of the molecule.
#' @param ... arguments to pass to molecule definition function.
#' @export
get_mol_paras <- function(name, ...) {
  if (length(name) == 1) {
    return(get(paste("get_", tolower(name), "_paras", sep = ""))(...))
  } else {
    out <- vector(mode = "list", length = length(name))
    for (n in 1:length(name)) {
      out[[n]] <- get(paste("get_", tolower(name[n]), "_paras", sep = ""))(...)
    }
    return(out)
  }
}

#' Return a character array of names that may be used with the 
#' \code{get_mol_paras} function.
#' @return a character array of names.
#' @export
get_mol_names <- function() {
  funs <- ls(getNamespace("spant"), all.names = TRUE)
  funs <- funs[!funs %in% c("get_mol_paras", "get_acq_paras",
                            "get_1h_brain_basis_paras",
                            "get_1h_brain_basis_paras_v1",
                            "get_1h_brain_basis_paras_v2",
                            "get_1h_brain_basis_paras_v3")]
  
  sub("_paras", "", sub("get_", "", funs[grep("get_.*_paras", funs)]))
}

#' Generate a \code{mol_parameters} object for a simple spin system with one
#' resonance.
#' @param name abbreviated name of the molecule.
#' @param chem_shift chemical shift of the resonance (PPM).
#' @param nucleus nucleus (1H, 31P...).
#' @param scale_factor multiplicative scaling factor. Note, this value can be
#' made complex to adjust the phase of the resonance.
#' @param lw linewidth in Hz.
#' @param lg Lorentz-Gauss lineshape parameter (between 0 and 1).
#' @param full_name long name of the molecule (optional).
#' @return mol_parameters object.
#' @export
get_uncoupled_mol <- function(name, chem_shift, nucleus, scale_factor, lw, lg,
                              full_name = NULL) {
  
  if (is.null(full_name)) full_name <- name
  
  j_coupling_mat <- matrix(0, 1, 1)
  spin_groups <- vector("list", length(chem_shift))
  for (n in 1:length(chem_shift)) {
    spin_group <- list(nucleus = nucleus[n], chem_shift = chem_shift[n], 
                       j_coupling_mat = j_coupling_mat,
                       scale_factor = scale_factor[n], lw = lw[n], lg = lg[n])
    spin_groups[[n]] <- spin_group
  }
  paras <- list(spin_groups = spin_groups, name = name, full_name = full_name)
  class(paras) <- "mol_parameters"
  paras
}

#' @export
print.mol_parameters <- function(x, ...) {
  cat(paste(c("Name        : ", x$name, "\n")), sep = "")
  cat(paste(c("Full name   : ", x$full_name, "\n")), sep = "")
  cat(paste(c("Spin groups : ", length(x$spin_groups), "\n")), sep = "")
  cat(paste(c("Source      : ", x$source, "\n")), sep = "")
  for (n in seq_len(length(x$spin_groups))) {
    cat("\n")
    cat(paste(c("Spin group ", n, "\n")), sep = "")
    cat("------------\n")
    cat(paste(c("Scaling factor : ", x$spin_groups[[n]]$scale_factor, "\n")), sep = "")
    cat(paste(c("Linewidth (Hz) : ", x$spin_groups[[n]]$lw, "\n")), sep = "")
    cat(paste(c("L/G lineshape  : ", x$spin_groups[[n]]$lg, "\n\n")), sep = "")
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
      print(j_mat, na.print = "-")
    }
  }
}

get_m_cr_ch2_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("-CrCH2", 3.913, "1H", -2, lw, lg,
                             "Inverted creatine CH2 group")
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  paras$source = source
  paras
}

get_ala_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13: 129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Ala", source = source,
                full_name = "Alanine")
  class(paras) <- "mol_parameters"
  paras
}

get_asc_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 4)
  chem_shift <- c(4.4965, 4.0072, 3.7469, 3.7194)
  j_coupling_mat <- matrix(0, 4, 4)
  j_coupling_mat[2,1] <- 2.055
  j_coupling_mat[3,2] <- 5.78
  j_coupling_mat[4,2] <- 7.373
  j_coupling_mat[4,3] <- -11.585
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Detection of an antioxidant profile in the humain brain in vivo via
              double editing with MEGA-PRESS. MRM. 2006; 56(6):1192-1199."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Asc", source = source,
                full_name = "Ascorbate")
  class(paras) <- "mol_parameters"
  paras
}

get_asp_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Asp", source = source,
                full_name = "Aspartate")
  class(paras) <- "mol_parameters"
  paras
}

get_cho_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- c("1H", "14N")
  chem_shift <- c(3.185, 0)
  j_coupling_mat <- matrix(0, 2, 2)
  j_coupling_mat[2,1] <- 0.57
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 9,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H", "1H", "14N")
  chem_shift <- c(4.054, 4.054, 3.501, 3.501, 0)
  j_coupling_mat <- matrix(0, 5, 5)
  j_coupling_mat[2,1] <- -14.1
  j_coupling_mat[3,1] <- 3.14
  j_coupling_mat[4,1] <- 6.979
  j_coupling_mat[5,1] <- 2.572
  j_coupling_mat[3,2] <- 7.011
  j_coupling_mat[4,2] <- 3.168
  j_coupling_mat[5,2] <- 2.681
  j_coupling_mat[4,3] <- -14.07
  
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for
              brain metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b), name = "Cho",
                source = source, full_name = "Choline")
  
  class(paras) <- "mol_parameters"
  paras
}

get_cr_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("Cr", 3.027, "1H", 3, lw, lg, "Creatine")
  paras_b <- get_uncoupled_mol("Cr", 3.913, "1H", 2, lw, lg)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras$source = source
  paras
}

get_cr_ch2_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("CrCH2", 3.913, "1H", 2, lw, lg,
                             "Creatine CH2 group")
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras$source = source
  paras
}

get_cr_ch3_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("CrCH3", 3.027, "1H", 3, lw, lg,
                             "Creatine CH2 group")
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras$source = source
  paras
}

get_gaba_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for
              brain metabolites. NMR Biomed. 2000; 13: 129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "GABA",
                source = source, full_name = "gamma-Aminobutyric acid")
  class(paras) <- "mol_parameters"
  paras
}

get_gaba_rt_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5 
  nucleus <- rep("1H", 6)
  chem_shift <- c(3.005, 3.005, 1.889, 1.889, 2.284, 2.284)
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153. Modified by MW for
              room temperature phantom scans."
  
  paras <- list(spin_groups = list(spin_group_a), name = "GABA",
                source = source, full_name = "gamma-Aminobutyric acid")
  class(paras) <- "mol_parameters"
  paras
}

get_gaba_jn_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "J-Difference Editing of Gamma-Aminobutyric Acid (GABA): Simulated
              and Experimental Multiplet Patterns. MRM 2013; 70:1183-1191."
  
  paras <- list(spin_groups = list(spin_group_a), name = "GABA",
                source = source, full_name = "gamma-Aminobutyric acid")
  class(paras) <- "mol_parameters"
  paras
}

get_gln_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Gln",
                source = source, full_name = "Glutamine")
  class(paras) <- "mol_parameters"
  paras
}

get_gsh_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b, spin_group_c),
                name = "GSH", source = source, full_name = "Glutathione")
  
  class(paras) <- "mol_parameters"
  paras
}

get_gly_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("Gly", 3.548, "1H", 2, lw, lg, "Glycine")
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras$source <- source
  paras
}

get_bhb_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5
  nucleus <- rep("1H", 6)
  chem_shift <- c(2.388, 2.294, 4.133, 1.186, 1.186, 1.186)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- -14.5
  j_coupling_mat[3,1] <- 7.3
  j_coupling_mat[3,2] <- 6.3
  j_coupling_mat[4,3] <- 6.3
  j_coupling_mat[5,3] <- 6.3
  j_coupling_mat[6,3] <- 6.3
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "In Vivo NMR Spectroscopy: Principles and Techniques,
              Robin A. de Graaf"
  
  paras <- list(spin_groups = list(spin_group_a), name = "BHB",
                source = source, full_name = "beta-Hydroxybutyrate")
  class(paras) <- "mol_parameters"
  paras
}

get_glyc_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5
  nucleus <- rep("1H", 5)
  chem_shift <- c(3.552, 3.640, 3.770, 3.640, 3.552)
  j_coupling_mat <- matrix(0, 5, 5)
  
  # these are possible also:
  #j_coupling_mat[2,1] <- -11.72
  #j_coupling_mat[3,1] <- 4.43
  #j_coupling_mat[3,2] <- 6.49
  #j_coupling_mat[4,3] <- 4.43
  #j_coupling_mat[5,3] <- 6.49
  #j_coupling_mat[5,4] <- -11.72
  #j_coupling_mat[2,1] <- -11.72
  #j_coupling_mat[3,1] <- 6.49
  #j_coupling_mat[3,2] <- 4.43
  #j_coupling_mat[4,3] <- 4.43
  #j_coupling_mat[5,3] <- 6.49
  #j_coupling_mat[5,4] <- -11.72
  
  # this ordering numbers gives the best agreement with the text book
  # eg:
  # get_mol_paras("glyc") %>% sim_mol(ft = 300e6, N = 1024*8) %>%
  # lb(1, 1) %>% plot(xlim = c(3.9, 3.4))
  
  j_coupling_mat[2,1] <- -11.72
  j_coupling_mat[3,1] <- 4.43
  j_coupling_mat[3,2] <- 6.49
  j_coupling_mat[4,3] <- 6.49
  j_coupling_mat[5,3] <- 4.43
  j_coupling_mat[5,4] <- -11.72
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "In Vivo NMR Spectroscopy: Principles and Techniques,
              Robin A. de Graaf"
  
  paras <- list(spin_groups = list(spin_group_a), name = "Glyc",
                source = source, full_name = "Glycerol")
  class(paras) <- "mol_parameters"
  paras
}

get_ins_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Ins",
                source = source, full_name = "myo-Inositol")
  class(paras) <- "mol_parameters"
  paras
}

get_ins_rt_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5
  nucleus <- rep("1H", 6)
  chem_shift <- c(3.5217, 4.0538, 3.5217, 3.608, 3.265, 3.608)
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153. Modified by MW for
              room temperature phantom scans."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Ins", source = source,
                full_name = "myo-Inositol")
  class(paras) <- "mol_parameters"
  paras
}

get_lip09_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("Lip09", 0.89, "1H", 3, 0.14 * ft / 1e6, 1,
                             "Lipid at 0.9 ppm")
  paras$source <- "LCModel manual."
  paras
}

get_lip13a_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("Lip13a", 1.28, "1H", 2, 0.15 * ft / 1e6, 1,
                             "Lipid at 1.3 ppm - component a")
  paras$source <- "LCModel manual."
  paras
}

get_lip13b_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("Lip13b", 1.28, "1H", 2, 0.089 * ft / 1e6, 1,
                             "Lipid at 1.3 ppm - component b")
  paras$source <- "LCModel manual."
  paras
}

get_lip20_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("Lip20", 2.04, "1H", 1.33, 0.15 * ft / 1e6, 1,
                             "Lipid at 2.0 ppm")
  paras_b <- get_uncoupled_mol("Lip20", 2.25, "1H", 0.67, 0.15 * ft / 1e6, 1)
  paras_c <- get_uncoupled_mol("Lip20", 2.8, "1H", 0.87, 0.2 * ft / 1e6, 1)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  paras$spin_groups[[3]] <- paras_c$spin_groups[[1]]
  paras$source <- "LCModel manual."
  paras
}

get_mm09_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("MM09", 0.91, "1H", 3, 0.14 * ft / 1e6, 1,
                             "Macromolecule at 0.9 ppm")
  paras$source <- "LCModel manual."
  paras
}

get_mm12_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("MM12", 1.21, "1H", 2, 0.15 * ft / 1e6, 1,
                             "Macromolecule at 1.2 ppm")
  paras$source <- "LCModel manual."
  paras
}

get_mm14_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("MM14", 1.43, "1H", 2, 0.17 * ft / 1e6, 1,
                             "Macromolecule at 1.4 ppm")
  paras$source <- "LCModel manual."
  paras
}

get_mm17_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("MM17", 1.67, "1H", 2, 0.15 * ft / 1e6, 1,
                             "Macromolecule at 1.7 ppm")
  paras$source <- "LCModel manual."
  paras
}

get_mm20_paras <- function(ft, ...) {
  paras <- get_uncoupled_mol("MM20", 2.08, "1H", 1.33, 0.15 * ft / 1e6, 1,
                             "Macromolecule at 2.0 ppm")
  paras_b <- get_uncoupled_mol("MM20", 2.25, "1H", 0.33, 0.2 * ft / 1e6, 1)
  paras_c <- get_uncoupled_mol("MM20", 1.95, "1H", 0.33, 0.15 * ft / 1e6, 1)
  paras_d <- get_uncoupled_mol("MM20", 3.0, "1H", 0.4, 0.2 * ft / 1e6, 1)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  paras$spin_groups[[3]] <- paras_c$spin_groups[[1]]
  paras$spin_groups[[4]] <- paras_d$spin_groups[[1]]
  paras$source <- "LCModel manual."
  paras
}

get_mm_3t_paras <- function(ft, ...) {
  paras   <- get_uncoupled_mol("MM", 0.90, "1H", 0.72, 21.20 / 128 * ft / 1e6,
                               1, "Macromolecule signal at 3 Tesla")
  paras_b <- get_uncoupled_mol("MM", 1.21, "1H", 0.28, 19.16 / 128 * ft / 1e6,
                               1)
  paras_c <- get_uncoupled_mol("MM", 1.38, "1H", 0.38, 15.90 / 128 * ft / 1e6,
                               1)
  paras_d <- get_uncoupled_mol("MM", 1.63, "1H", 0.05,  7.50 / 128 * ft / 1e6,
                               1)
  paras_e <- get_uncoupled_mol("MM", 2.01, "1H", 0.45, 29.03 / 128 * ft / 1e6,
                               1)
  paras_f <- get_uncoupled_mol("MM", 2.09, "1H", 0.36, 20.53 / 128 * ft / 1e6,
                               1)
  paras_g <- get_uncoupled_mol("MM", 2.25, "1H", 0.36, 17.89 / 128 * ft / 1e6,
                               1)
  paras_h <- get_uncoupled_mol("MM", 2.61, "1H", 0.04,  5.30 / 128 * ft / 1e6,
                               1)
  paras_i <- get_uncoupled_mol("MM", 2.96, "1H", 0.20, 14.02 / 128 * ft / 1e6,
                               1)
  paras_j <- get_uncoupled_mol("MM", 3.11, "1H", 0.11, 17.89 / 128 * ft / 1e6,
                               1)
  paras_k <- get_uncoupled_mol("MM", 3.67, "1H", 0.64, 33.52 / 128 * ft / 1e6,
                               1)
  paras_l <- get_uncoupled_mol("MM", 3.80, "1H", 0.07, 11.85 / 128 * ft / 1e6,
                               1)
  paras_m <- get_uncoupled_mol("MM", 3.96, "1H", 1.00, 37.48 / 128 * ft / 1e6,
                               1)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  paras$spin_groups[[3]] <- paras_c$spin_groups[[1]]
  paras$spin_groups[[4]] <- paras_d$spin_groups[[1]]
  paras$spin_groups[[5]] <- paras_e$spin_groups[[1]]
  paras$spin_groups[[6]] <- paras_f$spin_groups[[1]]
  paras$spin_groups[[7]] <- paras_g$spin_groups[[1]]
  paras$spin_groups[[8]] <- paras_h$spin_groups[[1]]
  paras$spin_groups[[9]] <- paras_i$spin_groups[[1]]
  paras$spin_groups[[10]] <- paras_j$spin_groups[[1]]
  paras$spin_groups[[11]] <- paras_k$spin_groups[[1]]
  paras$spin_groups[[12]] <- paras_l$spin_groups[[1]]
  paras$spin_groups[[13]] <- paras_m$spin_groups[[1]]
  paras$source <- "Birch et al Magn Reson Med. 2017 Jan; 77(1): 34-43."
  paras
}

get_naa_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153. Exlcuding the resonance
              at 7.82 ppm."
  
  paras <- list(spin_groups = list(spin_group_a,spin_group_b), name = "NAA",
                source = source, full_name = "N-acetylaspartate")
  class(paras) <- "mol_parameters"
  paras
}

get_naa_rt_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5
  nucleus <- c("1H")
  chem_shift <- c(2.008)
  j_coupling_mat <- matrix(0, 1, 1)
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 3,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H")
  chem_shift <- c(4.3817, 2.681, 2.4845)
  j_coupling_mat <- matrix(0, 3, 3)
  j_coupling_mat[2,1] <- 3.861
  j_coupling_mat[3,1] <- 9.821
  j_coupling_mat[3,2] <- -15.592
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153. Modified by MW for
              room temperature phantom scans. Exlcuding the resonance at 7.82
              ppm."
  
  paras <- list(spin_groups = list(spin_group_a,spin_group_b), name = "NAA",
                source = source, full_name = "N-acetylaspartate")
  class(paras) <- "mol_parameters"
  paras
}

get_naa2_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5
  nucleus <- c("1H")
  chem_shift <- c(2.008)
  j_coupling_mat <- matrix(0, 1, 1)
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 3,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H", "1H")
  chem_shift <- c(4.3817, 2.6727, 2.4863, 7.8205)
  j_coupling_mat <- matrix(0, 4, 4)
  j_coupling_mat[2,1] <- 3.861
  j_coupling_mat[3,1] <- 9.821
  j_coupling_mat[3,2] <- -15.592
  j_coupling_mat[4,1] <- 6.4
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153. Including the
              resonance at 7.82 ppm."
  
  paras <- list(spin_groups = list(spin_group_a,spin_group_b), name = "NAA",
                source = source, full_name = "N-acetylaspartate")
  class(paras) <- "mol_parameters"
  paras
}

get_naag_ch3_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5
  paras <- get_uncoupled_mol("NAAG", 2.042, "1H", 3, lw, lg,
                             "N-acetylaspartylglutamate")
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153. Note, only the acetyl
              moiety is simualted."
  
  paras$source <- source
  paras
}

get_naag_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5
  nucleus <- c("1H")
  chem_shift <- c(2.042)
  j_coupling_mat <- matrix(0, 1, 1)
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 3,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H")
  chem_shift <- c(4.607, 2.721, 2.519)
  j_coupling_mat <- matrix(0, 3, 3)
  j_coupling_mat[2,1] <- 4.412
  j_coupling_mat[3,1] <- 9.515
  j_coupling_mat[3,2] <- -15.910
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  nucleus <- c("1H", "1H", "1H", "1H", "1H", "1H")
  chem_shift <- c(4.128, 2.049, 1.881, 2.180, 2.190, 7.950)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- 4.61
  j_coupling_mat[3,1] <- 8.42
  j_coupling_mat[6,1] <- 7.46
  j_coupling_mat[3,2] <- -14.28
  j_coupling_mat[4,2] <- 10.56
  j_coupling_mat[5,2] <- 6.09
  j_coupling_mat[4,3] <- 4.9
  j_coupling_mat[5,3] <- 11.11
  j_coupling_mat[5,4] <- -15.28
  spin_group_c <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Chemical shifts are from the Govindaraju paper, J-couplings are
              from : Characterisation of the 1H and 13C NMR spectra of N-
              acetylaspartylglutamate and its detection in urine from patients
              with Canavan disease, J Pharm Biomed Anal. 2003 Mar
              10;31(3):455-63."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b, spin_group_c),
                name = "NAAG", source = source,
                full_name = "N-acetylaspartylglutamate")
  class(paras) <- "mol_parameters"
  paras
}

get_pch_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for
              brain metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b), name = "PCh",
                source = source, full_name = "Phosphocholine")
  
  class(paras) <- "mol_parameters"
  paras
}

get_pcr_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("PCr", 3.029, "1H", 3, lw, lg,
                             full_name = "Phosphocreatine")
  paras_b <- get_uncoupled_mol("PCr", 3.930, "1H", 2, lw, lg)
  paras$spin_groups[[2]] <- paras_b$spin_groups[[1]]
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras$source <- source
  paras
}

get_peth_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- c("1H", "1H", "1H", "1H", "31P", "14N")
  chem_shift <- c(3.9765, 3.9765, 3.216, 3.216, 0, 0)
  j_coupling_mat <- matrix(0, 6, 6)
  j_coupling_mat[2,1] <- -14.56
  j_coupling_mat[3,1] <- 3.182
  j_coupling_mat[4,1] <- 6.716
  j_coupling_mat[5,1] <- 7.288
  j_coupling_mat[6,1] <- 0.464
  j_coupling_mat[3,2] <- 7.204
  j_coupling_mat[4,2] <- 2.98
  j_coupling_mat[5,2] <- 7.088
  j_coupling_mat[6,2] <- 0.588
  j_coupling_mat[4,3] <- -14.71
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for
              brain metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "PEth",
                source = source, full_name = "Phosphoethanolamine")
  class(paras) <- "mol_parameters"
  paras
}

get_ser_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 3)
  chem_shift <- c(3.8347, 3.9379, 3.9764)
  j_coupling_mat <- matrix(0, 3, 3)
  j_coupling_mat[2,1] <- 5.979
  j_coupling_mat[3,1] <- 3.561
  j_coupling_mat[3,2] <- -12.254
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for
              brain metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Ser", source = source,
                full_name = "Serine")
  class(paras) <- "mol_parameters"
  paras
}

get_sins_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("sIns", 3.34, "1H", 6, lw, lg, "scyllo-Inositol")
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras$source <- source
  paras
}

get_h2o_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("H2O", 4.65, "1H", 2, lw, lg, "water")
  
  source <- "Singlet at 4.65 ppm."
  
  paras$source <- source
  paras
}

get_msm_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 1
  paras <- get_uncoupled_mol("MSM", 3.142, "1H", 6, lw, lg,
                             "Methylsulfonylmethane")
  
  source <- "Kaiser LG, Russell D, Maschmeyer T, Redfern RL, Inglis BA.
    Methylsulfonylmethane (MSM): A chemical shift reference for 1 H MRS of human
    brain. Magn Reson Med. 2020 Apr;83(4):1157-1167. doi: 10.1002/mrm.27997.
    Epub 2019 Sep 30. PMID: 31566256."
  
  paras$source <- source
  paras
}

get_ace_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("Ace", 2.222, "1H", 3, lw, lg,
                             "Acetone")
  
  source <- "In Vivo NMR Spectroscopy Principles and Techniques, Robin A. de
             Graaf, Third Edition, 2019, Wiley."
  
  paras$source <- source
  paras
}

get_pyr_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("Pyr", 2.3580, "1H", 3, lw, lg, "Pyruvate")
  
  source <- "In Vivo NMR Spectroscopy Principles and Techniques, Robin A. de
             Graaf, Third Edition, 2019, Wiley."
  
  paras$source <- source
  paras
}

get_suc_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("Suc", 2.3920, "1H", 4, lw, lg, "Succinate")
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras$source <- source
  paras
}

get_tau_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for
               brain metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Tau", source = source,
                full_name = "Taurine")
  class(paras) <- "mol_parameters"
  paras
}

get_thr_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 3)
  chem_shift <- c(3.578, 4.246, 1.316)
  j_coupling_mat <- matrix(0, 3, 3)
  j_coupling_mat[2,1] <- 4.917
  j_coupling_mat[3,2] <- 6.35
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "In Vivo NMR Spectroscopy: Principles and Techniques,
             Robin A. de Graaf"
  
  paras <- list(spin_groups = list(spin_group_a), name = "Thr", source = source,
                full_name = "Threonine")
  class(paras) <- "mol_parameters"
  paras
}

get_lac_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
              metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Lac", source = source,
                full_name = "Lactate")
  class(paras) <- "mol_parameters"
  paras
}

get_lac_rt_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5 
  nucleus <- rep("1H", 4)
  chem_shift <- c(4.0974, 1.3142, 1.3142, 1.3142)
  j_coupling_mat <- matrix(0, 4, 4)
  
  j_coupling_mat[2,1] <- 6.833
  j_coupling_mat[3,1] <- 6.833
  j_coupling_mat[4,1] <- 6.833
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
             metabolites. NMR Biomed. 2000; 13:129-153. Modified by MW for
             room temperature phantom scans."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Lac", source = source,
                full_name = "Lactate")
  class(paras) <- "mol_parameters"
  paras
}

get_glu_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 5)
  chem_shift <- c(3.7433, 2.0375, 2.1200, 2.3378, 2.352)
  j_coupling_mat <- matrix(0, 5, 5)
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
             metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Glu", source = source,
                full_name = "Glutamate")
  class(paras) <- "mol_parameters"
  paras
}

get_glu_rt_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 0.5 
  nucleus <- rep("1H", 5)
  chem_shift <- c(3.75, 2.0475, 2.1200, 2.3378, 2.352)
  j_coupling_mat <- matrix(0, 5, 5)
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
             metabolites. NMR Biomed. 2000; 13:129-153. Modified by MW for
             room temperature phantom scans."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Glu", source = source,
                full_name = "Glutamate")
  class(paras) <- "mol_parameters"
  paras
}

get_a_glc_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
             metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "aGlc",
                source = source, full_name = "alpha-D-glucose")
  class(paras) <- "mol_parameters"
  paras
}

get_b_glc_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 7)
  chem_shift <- c(4.630, 3.230, 3.473, 3.387, 3.450, 3.882, 3.707)
  j_coupling_mat <- matrix(0, 7, 7)
  j_coupling_mat[2,1] <- 8.0
  j_coupling_mat[3,2] <- 9.1
  j_coupling_mat[4,3] <- 9.4
  j_coupling_mat[5,4] <- 8.9
  j_coupling_mat[6,5] <- 1.6
  j_coupling_mat[7,5] <- 5.4
  j_coupling_mat[7,6] <- -12.3
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
             metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a), name = "bGlc",
                source = source, full_name = "beta-D-glucose")
  class(paras) <- "mol_parameters"
  paras
}

get_glc_paras <- function(lw = NULL, lg = 0, ...) {
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
                       j_coupling_mat = j_coupling_mat, scale_factor = 0.36,
                       lw = lw, lg = lg)
  
  nucleus <- rep("1H", 7)
  chem_shift <- c(4.630, 3.230, 3.473, 3.387, 3.450, 3.882, 3.707)
  j_coupling_mat <- matrix(0, 7, 7)
  j_coupling_mat[2,1] <- 8.0
  j_coupling_mat[3,2] <- 9.1
  j_coupling_mat[4,3] <- 9.4
  j_coupling_mat[5,4] <- 8.9
  j_coupling_mat[6,5] <- 1.6
  j_coupling_mat[7,5] <- 5.4
  j_coupling_mat[7,6] <- -12.3
  
  spin_group_b <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 0.64,
                       lw = lw, lg = lg)
  
  source <- "Proton NMR chemical shifts and coupling constants for brain
             metabolites. NMR Biomed. 2000; 13:129-153. This is a combination
             of alpha-glc (36%) and beta-glc (64%)."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b), name = "Glc",
                source = source, full_name = "Glucose")
  class(paras) <- "mol_parameters"
  paras
}

get_gpc_paras <- function(lw = NULL, lg = 0, ...) {
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
  
  source <- "Corrigendum: Proton NMR chemical shifts and coupling constants for
             brain metabolites. NMR Biomed. 2000; 13:129-153."
  
  paras <- list(spin_groups = list(spin_group_a, spin_group_b, spin_group_c),
                name = "GPC", source = source,
                full_name = "Glycerophosphocholine")
  class(paras) <- "mol_parameters"
  paras
}

get_2hg_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 5)
  chem_shift <- c(4.022, 1.825, 1.977, 2.221, 2.272)
  j_coupling_mat <- matrix(0, 5, 5)
  j_coupling_mat[2,1] <- 7
  j_coupling_mat[3,1] <- 4.1
  j_coupling_mat[3,2] <- -14
  j_coupling_mat[4,2] <- 5.3
  j_coupling_mat[5,2] <- 10.6
  j_coupling_mat[4,3] <- 10.4
  j_coupling_mat[5,3] <- 6
  j_coupling_mat[5,4] <- -15
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "2-hydroxyglutarate detection by magnetic resonance spectroscopy
             in subjects with IDH-mutated gliomas. Nat Med. 2012 Jan
             26;18(4):624-9. Numbers are given in the supplementary
             materials."
  
  paras <- list(spin_groups = list(spin_group_a), name = "2HG", source = source,
                full_name = "2-hydroxyglutarate")
  class(paras) <- "mol_parameters"
  paras
}

get_cit_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("1H", 2)
  chem_shift <- c(2.6735, 2.5265)
  j_coupling_mat <- matrix(0, 2, 2)
  j_coupling_mat[2,1] <- -16.1
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "Optimizing PRESS localized citrate detection at 3 Tesla. Magn Reson
             Med. 2005 Jul;54(1):51-8 and van der Graaf M, Heerschap A. Effect
             of Cation Binding on the Proton Chemical Shifts and the Spin-Spin
             Coupling Constant of Citrate. J Magn Reson B. 1996
             Jul;112(1):58-62. Central frequency assumed to be 2.6 ppm, delta
             chemical shift of 0.147 ppm and j-coupling value of 16.1 Hz."
  
  paras <- list(spin_groups = list(spin_group_a), name = "Cit", source = source,
                full_name = "Citrate")
  class(paras) <- "mol_parameters"
  paras
}

get_atp_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  nucleus <- rep("31P", 3)
  chem_shift <- c(-7.56, -16.18, -2.53)
  j_coupling_mat <- matrix(0, 3, 3)
  j_coupling_mat[2,1] <- 16.3
  j_coupling_mat[3,2] <- 15.4
  
  spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
                       j_coupling_mat = j_coupling_mat, scale_factor = 1,
                       lw = lw, lg = lg)
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras <- list(spin_groups = list(spin_group_a), name = "ATP", source = source,
                full_name = "Adenosine triphosphate")
  
  class(paras) <- "mol_parameters"
  paras
}

# get_atp_31p_paras <- function(lw = NULL, lg = 0, ...) {
#   if (is.null(lw)) lw = 2
#   nucleus <- rep("31P", 3)
#   chem_shift <- c(-7.616, -16.26, -2.6)
#   j_coupling_mat <- matrix(0, 3, 3)
#   j_coupling_mat[2,1] <- 16.3
#   j_coupling_mat[3,2] <- 16.1
#   
#   spin_group_a <- list(nucleus = nucleus, chem_shift = chem_shift, 
#                        j_coupling_mat = j_coupling_mat, scale_factor = 1,
#                        lw = lw, lg = lg)
#   
#   source <- "Unknown"
#   
#   paras <- list(spin_groups = list(spin_group_a), name = "ATP", source = source,
#                 full_name = "Adenosine triphosphate")
#   
#   class(paras) <- "mol_parameters"
#   paras
# }

get_gpc_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("GPC", 2.94, "31P", 1, lw, lg,
                             "Glycerophosphocholine")
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras$source <- source
  paras
}

get_gpe_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("GPE", 3.49, "31P", 1, lw, lg,
                             "Glycerol phosphorylethanolamine")
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras$source <- source
  paras
}

get_nadh_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("NADH", -8.13, "31P", 1, lw, lg,
                             "Nicotinamide adenine dinucleotide, reduced")
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras$source <- source
  paras
}

get_nadp_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("NAD+", -8.31, "31P", 1, lw, lg,
                             "Nicotinamide adenine dinucleotide, oxidized")
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras$source <- source
  paras
}

get_pch_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("PCh", 6.23, "31P", 1, lw, lg, "Phosphorylcholine")
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras$source <- source
  paras
}

get_pcr_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("PCr", 0, "31P", 1, lw, lg, "Phosphocreatine")
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras$source <- source
  paras
}

get_pe_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("PE", 6.77, "31P", 1, lw, lg,
                             "Phosphorylethanolamine ")
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras$source <- source
  paras
}

get_pi_31p_paras <- function(lw = NULL, lg = 0, ...) {
  if (is.null(lw)) lw = 2
  paras <- get_uncoupled_mol("PI", 4.84, "31P", 1, lw, lg,
                             "Inorganic phosphate")
  
  source <- "NMR Biomed. 2015 June; 28(6): 633-641."
  
  paras$source <- source
  paras
}