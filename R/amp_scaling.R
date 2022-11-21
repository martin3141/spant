#' Apply water reference scaling to a fitting results object to yield metabolite 
#' quantities in millimolar (mM) units (mol / kg of tissue water).
#'
#' Details of this method can be found in "Use of tissue water as a
#' concentration reference for proton spectroscopic imaging" by Gasparovic et al
#' MRM 2006 55(6):1219-26. 1.5 Tesla relaxation assumptions are taken from this
#' paper. For 3 Tesla data, relaxation assumptions are taken from "NMR 
#' relaxation times in the human brain at 3.0 tesla" by Wansapura et al J Magn
#' Reson Imaging 1999 9(4):531-8. 
#' 
#' @param fit_result result object generated from fitting.
#' @param ref_data water reference MRS data object.
#' @param p_vols a numeric vector of partial volumes expressed as percentages.
#' For example, a voxel containing 100% white matter tissue would use : 
#' p_vols = c(WM = 100, GM = 0, CSF = 0).
#' @param te the MRS TE in seconds.
#' @param tr the MRS TR in seconds.
#' @param ... additional arguments to get_td_amp function.
#' @return A \code{fit_result} object with a rescaled results table.
#' @export
scale_amp_molal_pvc <- function(fit_result, ref_data, p_vols, te, tr, ...){
  
  if (!identical(dim(fit_result$data$data)[2:6], dim(ref_data$data)[2:6])) {
    stop("Mismatch between fit result and reference data dimensions.")
  }
  
  # check if res_tab_unscaled exists, and if not create it
  if (is.null(fit_result$res_tab_unscaled)) {
    fit_result$res_tab_unscaled <- fit_result$res_tab
  } else {
    fit_result$res_tab <- fit_result$res_tab_unscaled
  }
  
  B0 <- round(fit_result$data$ft / 42.58e6, 1)
  corr_factor <- get_corr_factor(te, tr, B0, p_vols[["GM"]], p_vols[["WM"]],
                                 p_vols[["CSF"]])
  
  amp_cols <- fit_result$amp_cols
  
  w_amp <- as.numeric(get_td_amp(ref_data, ...))
  
  if (length(w_amp) != nrow(fit_result$res_tab)) {
    stop("Mismatch between fit result and reference data.")
  }
  
  fit_result$res_tab$w_amp <- w_amp
  
  fit_result$res_tab$GM_vol    <- p_vols[["GM"]]
  fit_result$res_tab$WM_vol    <- p_vols[["WM"]]
  fit_result$res_tab$CSF_vol   <- p_vols[["CSF"]]
  if ("Other" %in% names(p_vols)) {
    fit_result$res_tab$Other_vol <- p_vols[["Other"]]
  }
  fit_result$res_tab$GM_frac   <- p_vols[["GM"]] / 
                                 (p_vols[["GM"]] + p_vols[["WM"]])
  
  fit_result$res_tab_unscaled$GM_vol    <- p_vols[["GM"]]
  fit_result$res_tab_unscaled$WM_vol    <- p_vols[["WM"]]
  fit_result$res_tab_unscaled$CSF_vol   <- p_vols[["CSF"]]
  if ("Other" %in% names(p_vols)) {
    fit_result$res_tab_unscaled$Other_vol <- p_vols[["Other"]]
  }
  fit_result$res_tab_unscaled$GM_frac   <- p_vols[["GM"]] / 
                                          (p_vols[["GM"]] + p_vols[["WM"]])
  
  # append tables with %GM, %WM, %CSF and %Other
  pvc_cols <- 6:(5 + amp_cols * 2)
  fit_result$res_tab[, pvc_cols] <- fit_result$res_tab[, pvc_cols] *
                                    corr_factor / w_amp
  
  return(fit_result)
}
  
#' Apply water reference scaling to a fitting results object to yield metabolite 
#' quantities in millimolar (mM) units (mol / Litre of tissue).
#' 
#' See the LCModel manual section on water-scaling for details on the
#' assumptions and relevant references. Use this type of concentration scaling
#' to compare fit results with LCModel and TARQUIN defaults. Otherwise
#' scale_amp_molal_pvc is generally the preferred method.
#' 
#' @param fit_result a result object generated from fitting.
#' @param ref_data water reference MRS data object.
#' @param w_att water attenuation factor (default = 0.7). Assumes water T2 of
#' 80ms and a TE = 30 ms. exp(-30ms / 80ms) ~ 0.7.
#' @param w_conc assumed water concentration (default = 35880). Default value
#' corresponds to typical white matter. Set to 43300 for gray matter, and 55556 
#' for phantom measurements.
#' @param ... additional arguments to get_td_amp function.
#' @return a \code{fit_result} object with a rescaled results table.
#' @export
scale_amp_molar <- function(fit_result, ref_data, w_att = 0.7, w_conc = 35880,
                            ...) {
  
  if (!identical(dim(fit_result$data$data)[2:6], dim(ref_data$data)[2:6])) {
    stop("Mismatch between fit result and reference data dimensions.")
  }
  
  # check if res_tab_unscaled exists, and if not create it
  if (is.null(fit_result$res_tab_unscaled)) {
    fit_result$res_tab_unscaled <- fit_result$res_tab
  } else {
    fit_result$res_tab <- fit_result$res_tab_unscaled
  }
  
  w_amp <- as.numeric(get_td_amp(ref_data, ...))
  
  if (length(w_amp) != nrow(fit_result$res_tab)) {
    stop("Mismatch between fit result and reference data.")
  }
  
  fit_result$res_tab$w_amp <- w_amp
  
  amp_cols <- fit_result$amp_cols
  ws_cols <- 6:(5 + amp_cols * 2)
  
  fit_result$res_tab[, ws_cols] <- (fit_result$res_tab[, ws_cols] * w_att *
                                    w_conc / w_amp)
    
  fit_result
}

#' Scale metabolite amplitudes as a ratio to the unsuppressed water amplitude.
#' @param fit_result a result object generated from fitting.
#' @param ref_data a water reference MRS data object.
#' @param ... additional arguments to get_td_amp function.
#' @return a \code{fit_result} object with a rescaled results table.
#' @export
scale_amp_water_ratio <- function(fit_result, ref_data, ...) {
  
  if (!identical(dim(fit_result$data$data)[2:6], dim(ref_data$data)[2:6])) {
    stop("Mismatch between fit result and reference data dimensions.")
  }
  
  # check if res_tab_unscaled exists, and if not create it
  if (is.null(fit_result$res_tab_unscaled)) {
    fit_result$res_tab_unscaled <- fit_result$res_tab
  } else {
    fit_result$res_tab <- fit_result$res_tab_unscaled
  }
  
  w_amp <- as.numeric(get_td_amp(ref_data, ...))
  
  if (length(w_amp) != nrow(fit_result$res_tab)) {
    stop("Mismatch between fit result and reference data.")
  }
  
  fit_result$res_tab$w_amp <- w_amp
  
  amp_cols <- fit_result$amp_cols
  ws_cols <- 6:(5 + amp_cols * 2)
  
  fit_result$res_tab[, ws_cols] <- fit_result$res_tab[, ws_cols] / w_amp
  
  fit_result
}

#' Scale fitted amplitudes to a ratio of signal amplitude.
#' @param fit_result a result object generated from fitting.
#' @param name the signal name to use as a denominator (usually, "tCr" or 
#' "tNAA").
#' @return a \code{fit_result} object with a rescaled results table.
#' @export
scale_amp_ratio <- function(fit_result, name) {
  
  # check if res_tab_unscaled exists, and if not create it
  if (is.null(fit_result$res_tab_unscaled)) {
    fit_result$res_tab_unscaled <- fit_result$res_tab
  } else {
    fit_result$res_tab <- fit_result$res_tab_unscaled
  }
  
  ratio_amp <- as.numeric(fit_result$res_tab[[name]])
  
  amp_cols <- fit_result$amp_cols
  ws_cols <- 6:(5 + amp_cols * 2)
  
  fit_result$res_tab[, ws_cols] <- fit_result$res_tab[, ws_cols] / ratio_amp
  
  fit_result
}

#' Scale fitted amplitudes to a ratio of signal amplitude.
#' @param fit_result a result object generated from fitting.
#' @param value the number use as a denominator.
#' @return a \code{fit_result} object with a rescaled results table.
#' @export
scale_amp_ratio_value <- function(fit_result, value) {
  
  # check if res_tab_unscaled exists, and if not create it
  if (is.null(fit_result$res_tab_unscaled)) {
    fit_result$res_tab_unscaled <- fit_result$res_tab
  } else {
    fit_result$res_tab <- fit_result$res_tab_unscaled
  }
  
  ratio_amp <- value
  
  amp_cols <- fit_result$amp_cols
  ws_cols <- 6:(5 + amp_cols * 2)
  
  fit_result$res_tab[, ws_cols] <- fit_result$res_tab[, ws_cols] / ratio_amp
  
  fit_result
}

get_corr_factor <- function(te, tr, B0, gm_vol, wm_vol, csf_vol) {
  # Correction factor calcualted according to the method of Gasparovic et al (MRM 55:1219-1226 2006)
  # NOTE - gives concs as Mol/kg of water NOT Mol/liter of tissue like default LCM/TQN analysis.
  if ((B0 == 3.0) | (B0 == 2.9)) {
    # Wanasapura values given in Harris paper
    t1_gm    <- 1.331
    t2_gm    <- 0.110
    t1_wm    <- 0.832
    t2_wm    <- 0.0792
    t1_csf   <- 3.817
    t2_csf   <- 0.503
    t1_metab <- 1.15
    t2_metab <- 0.3
  } else if (B0 == 1.5) {
    # values from Gasparovic 2006 MRM paper
    t1_gm    <- 1.304
    t2_gm    <- 0.093
    t1_wm    <- 0.660
    t2_wm    <- 0.073
    t1_csf   <- 2.93
    t2_csf   <- 0.23
    t1_metab <- 1.15
    t2_metab <- 0.3
  } else {
    stop("Error. Relaxation values not available for this field strength.")
  }
  
  # MR-visable water densities
  gm_vis  <- 0.78
  wm_vis  <- 0.65
  csf_vis <- 0.97
  
  # molal concentration (moles/gram) of MR-visible water
  water_conc <- 55510.0
  
  # fractions of water attributable to GM, WM and CSF
  f_gm  <- gm_vol * gm_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  f_wm  <- wm_vol * wm_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  f_csf <- csf_vol * csf_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  #This might give the result in Mol/kg?
  #f_gm  <- gm_vol  * gm_vis   / ( gm_vol + wm_vol + csf_vol )
  #f_wm  <- wm_vol  * wm_vis   / ( gm_vol + wm_vol + csf_vol )
  #f_csf <- csf_vol * csf_vis  / ( gm_vol + wm_vol + csf_vol )
  
  # Relaxtion attenuation factors
  R_h2o_gm  <- exp(-te / t2_gm) * (1.0 - exp(-tr / t1_gm))
  R_h2o_wm  <- exp(-te / t2_wm) * (1.0 - exp(-tr / t1_wm))
  R_h2o_csf <- exp(-te / t2_csf) * (1.0 - exp(-tr / t1_csf))
  R_metab   <- exp(-te / t2_metab) * (1.0 - exp(-tr / t1_metab))
  
  corr_factor <- ((f_gm * R_h2o_gm + f_wm * R_h2o_wm + f_csf * R_h2o_csf) / 
                 ((1 - f_csf) * R_metab)) * water_conc
  
  return(corr_factor)
}

#' Convert default LCM/TARQUIN concentration scaling to molal units with partial 
#' volume correction.
#' @param fit_result a \code{fit_result} object to apply partial volume 
#' correction.
#' @param p_vols a numeric vector of partial volumes expressed as percentages.
#' For example, a voxel containing 100% white matter tissue would use : 
#' p_vols = c(WM = 100, GM = 0, CSF = 0).
#' @param te the MRS TE.
#' @param tr the MRS TR.
#' @return a \code{fit_result} object with a rescaled results table.
#' @export
apply_pvc <- function(fit_result, p_vols, te, tr){
  
  # check if res_tab_unscaled exists, and if not create it
  if (is.null(fit_result$res_tab_unscaled)) {
    fit_result$res_tab_unscaled <- fit_result$res_tab
  } else {
    fit_result$res_tab <- fit_result$res_tab_unscaled
  }
  
  B0 <- round(fit_result$data$ft / 42.58e6,1)
  corr_factor <- get_corr_factor(te, tr, B0, p_vols[["GM"]], p_vols[["WM"]],
                                 p_vols[["CSF"]])
  
  amp_cols <- fit_result$amp_cols
  default_factor <- 35880 * 0.7
  fit_result$res_tab$GM_vol <- p_vols[["GM"]]
  fit_result$res_tab$WM_vol <- p_vols[["WM"]]
  fit_result$res_tab$CSF_vol <- p_vols[["CSF"]]
  
  if ("Other" %in% names(p_vols)) {
    fit_result$res_tab$Other_vol <- p_vols[["Other"]]
  }
  
  # append tables with %GM, %WM, %CSF and %Other
  pvc_cols <- 6:(5 + amp_cols * 2)
  fit_result$res_tab[, pvc_cols] <- fit_result$res_tab[, pvc_cols] /
                                    default_factor * corr_factor
  return(fit_result)
}
