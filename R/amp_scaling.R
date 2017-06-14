#' Apply partial volume correction to a fitting result object.
#' @param fit_result result object generated from fitting.
#' @param ref_data water reference MRS data object.
#' @param p_vols a numeric vector of partial volumes.
#' @param te the MRS TE.
#' @param tr the MRS TR.
#' @return A \code{fit_result} object with a res_tab_molal_pvc data table added.
#' @export
scale_amp_molal_pvc <- function(fit_result, ref_data, p_vols, te, tr){
  B0 = round(fit_result$data$ft / 42.58e6,1)
  corr_factor <- get_corr_factor(te, tr, B0, p_vols[["GM"]], p_vols[["WM"]],
                                 p_vols[["CSF"]])
  
  amp_cols = fit_result$amp_cols
  
  w_amp <- as.numeric(get_td_amp(ref_data))
  fit_result$res_tab$w_amp = w_amp
  
  fit_result$res_tab$GM_vol = p_vols[["GM"]]
  fit_result$res_tab$WM_vol = p_vols[["WM"]]
  fit_result$res_tab$CSF_vol = p_vols[["CSF"]]
  fit_result$res_tab$Other_vol = p_vols[["Other"]]
  fit_result$res_tab_molal_pvc <- fit_result$res_tab
  
  # append tables with %GM, %WM, %CSF and %Other
  pvc_cols <- 6:(5 + amp_cols * 2)
  fit_result$res_tab_molal_pvc[, pvc_cols] <- fit_result$res_tab_molal_pvc[, pvc_cols] *
                                              corr_factor / w_amp
  return(fit_result)
}
  
#' Apply water reference scaling to a fitting results object to yield metabolite 
#' quantities in millimolar (mM) units (mol/litre).
#' @param fit_result A result object generated from fitting.
#' @param ref_data Water reference MRS data object.
#' @param w_att Water attenuation factor (default = 0.7).
#' @param w_conc Assumed water concentration (default = 35880).
#' @return A \code{fit_result} object with a res_tab_molar data table added.
#' @export
scale_amp_molar <- function(fit_result, ref_data, w_att = 0.7, 
                                 w_conc = 35880) {
  
  w_amp <- as.numeric(get_td_amp(ref_data))
  
  fit_result$res_tab$w_amp = w_amp
  
  amp_cols = fit_result$amp_cols
  ws_cols <- 6:(5 + amp_cols * 2)
  
  fit_result$res_tab_molar <- fit_result$res_tab
  fit_result$res_tab_molar[, ws_cols] <- (fit_result$res_tab_molar[, ws_cols] * w_att * 
                                   w_conc / w_amp)
  
  fit_result
}

#' Scale metabolite amplitudes as a ratio to the unsupressed water amplitude.
#' @param fit_result a result object generated from fitting.
#' @param ref_data a water reference MRS data object.
#' @return a \code{fit_result} object with a res_tab_water_ratio data table added.
#' @export
scale_amp_water_ratio <- function(fit_result, ref_data) {
  
  w_amp <- as.numeric(get_td_amp(ref_data))
  
  fit_result$res_tab$w_amp = w_amp
  
  amp_cols = fit_result$amp_cols
  ws_cols <- 6:(5 + amp_cols * 2)
  
  fit_result$res_tab_water_ratio <- fit_result$res_tab
  fit_result$res_tab_water_ratio[, ws_cols] <- (fit_result$res_tab_water_ratio[, ws_cols] / w_amp)
  
  fit_result
}

#' Scale fitted amplitudes to a ratio of signal amplitude.
#' @param fit_result a result object generated from fitting.
#' @param name the signal name to use as a denominator (usually, "TCr" or 
#' "TNAA").
#' @return a \code{fit_result} object with a res_tab_ratio data table added.
#' @export
scale_amp_ratio <- function(fit_result, name) {
  
  ratio_amp <- as.numeric(fit_result$res_tab[name])
  
  amp_cols = fit_result$amp_cols
  ws_cols <- 6:(5 + amp_cols * 2)
  
  fit_result$res_tab_ratio <- fit_result$res_tab
  fit_result$res_tab_ratio[, ws_cols] <- (fit_result$res_tab_ratio[, ws_cols] / ratio_amp)
  
  fit_result
}


get_corr_factor <- function(te, tr, B0, gm_vol, wm_vol, csf_vol) {
  # Correction factor calcualted according to the method of Gasparovic et al (MRM 55:1219-1226 2006)
  # NOTE - gives concs as Mol/kg of water NOT Mol/liter of tissue like default LCM/TQN analysis.
  if (B0 == 3.0) {
    # Wanasapura values given in Harris paper
    t1_gm    = 1.331
    t2_gm    = 0.110
    t1_wm    = 0.832
    t2_wm    = 0.0792
    t1_csf   = 3.817
    t2_csf   = 0.503
    t1_metab = 1.15
    t2_metab = 0.3
  } else if (B0 == 1.5) {
    # values from Gasparovic 2006 MRM paper
    t1_gm    = 1.304
    t2_gm    = 0.093
    t1_wm    = 0.660
    t2_wm    = 0.073
    t1_csf   = 2.93
    t2_csf   = 0.23
    t1_metab = 1.15
    t2_metab = 0.3
  } else {
    stop("Error. Relaxation values not available for this field strength.")
  }
  
  # MR-visable water densities
  gm_vis  = 0.78
  wm_vis  = 0.65
  csf_vis = 0.97
  
  # molal concentration (moles/gram) of MR-visible water
  water_conc = 55510.0
  
  # fractions of water attributable to GM, WM and CSF
  f_gm  = gm_vol * gm_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  f_wm  = wm_vol * wm_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  f_csf = csf_vol * csf_vis / 
          (gm_vol * gm_vis + wm_vol * wm_vis + csf_vol * csf_vis)
  
  #This might give the result in Mol/kg?
  #f_gm  = gm_vol  * gm_vis   / ( gm_vol + wm_vol + csf_vol )
  #f_wm  = wm_vol  * wm_vis   / ( gm_vol + wm_vol + csf_vol )
  #f_csf = csf_vol * csf_vis  / ( gm_vol + wm_vol + csf_vol )
  
  # Relaxtion attenuation factors
  R_h2o_gm  = exp(-te / t2_gm) * (1.0 - exp(-tr / t1_gm))
  R_h2o_wm  = exp(-te / t2_wm) * (1.0 - exp(-tr / t1_wm))
  R_h2o_csf = exp(-te / t2_csf) * (1.0 - exp(-tr / t1_csf))
  R_metab   = exp(-te / t2_metab) * (1.0 - exp(-tr / t1_metab))
  
  corr_factor = ((f_gm * R_h2o_gm + f_wm * R_h2o_wm + f_csf * R_h2o_csf) / 
                ((1 - f_csf) * R_metab)) * water_conc
  
  return(corr_factor)
}

#' Convert default LCM/TARQUIN concentration scaling to molal units with partial 
#' volume correction.
#' @param fit_result a \code{fit_result} object to apply partial volume 
#' correction.
#' @param p_vols a numeric vector of partial volumes.
#' @param te the MRS TE.
#' @param tr the MRS TR.
#' @return a \code{fit_result} object with an added results table: 
#' "res_tab_molal_pvc".
#' @export
apply_pvc <- function(fit_result, p_vols, te, tr){
  #te = result$data$te
  B0 = round(fit_result$data$ft / 42.58e6,1)
  corr_factor <- get_corr_factor(te, tr, B0, p_vols[["GM"]], p_vols[["WM"]],
                                 p_vols[["CSF"]])
  
  amp_cols = fit_result$amp_cols
  default_factor = 35880 * 0.7
  fit_result$res_tab$GM_vol = p_vols[["GM"]]
  fit_result$res_tab$WM_vol = p_vols[["WM"]]
  fit_result$res_tab$CSF_vol = p_vols[["CSF"]]
  fit_result$res_tab$Other_vol = p_vols[["Other"]]
  fit_result$res_tab_molal_pvc <- fit_result$res_tab
  
  # append tables with %GM, %WM, %CSF and %Other
  pvc_cols <- 6:(5 + amp_cols * 2)
  fit_result$res_tab_molal_pvc[, pvc_cols] <- fit_result$res_tab_molal_pvc[, pvc_cols] /
                                    default_factor * corr_factor
  return(fit_result)
}
