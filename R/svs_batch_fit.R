#' Standard SVS 1H brain analysis pipeline.
#' @param metab filepath or mrs_data object containing MRS metabolite data.
#' @param w_ref filepath or mrs_data object containing MRS water reference data.
#' @param output_dir directory path to output fitting results.
#' @param basis precompiled basis set object to use for analysis.
#' @param p_vols a numeric vector of partial volumes expressed as percentages.
#' Defaults to 100% white matter. A voxel containing 100% gray matter tissue
#' would use : p_vols = c(WM = 0, GM = 100, CSF = 0).
#' @param append_basis names of extra signals to add to the default basis. Eg 
#' append_basis = c("peth", "cit"). Cannot be used with precompiled basis sets.
#' @param remove_basis names of signals to remove from the basis. Cannot be used
#' with precompiled basis sets.
#' @param dfp_corr perform dynamic frequency and phase correction using the RATS
#' method.
#' @param omit_bad_dynamics detect and remove bad dynamics.
#' @param te metabolite mrs data echo time in seconds. If not supplied this will
#' be guessed from the metab data file.
#' @param tr metabolite mrs data repetition time in seconds. If not supplied
#' this will be guessed from the metab data file.
#' @param output_ratio optional string to specify a metabolite ratio to output.
#' Defaults to "tCr" and multiple metabolites may be specified for multiple
#' outputs. Set as NULL to omit.
#' @param ecc option to perform water reference based eddy current correction,
#' defaults to FALSE.
#' @param abfit_opts options to pass to ABfit.
#' @param verbose output potentially useful information.
#' @export
svs_1h_brain_analysis_dev <- function(metab, w_ref = NULL, output_dir = NULL,
                                      basis = NULL, p_vols = NULL,
                                      append_basis = NULL, remove_basis = NULL,
                                      dfp_corr = TRUE, omit_bad_dynamics = TRUE,
                                      te = NULL, tr = NULL,
                                      output_ratio = "tCr", ecc = FALSE,
                                      abfit_opts = NULL, verbose = FALSE) {
  
  # TODO
  # Auto sequence detection and override option
  # Realistic PRESS sim for B0 > 2.9T
  
  if (!is.null(basis) & !is.null(append_basis)) {
    stop("basis and append_basis options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(basis) & !is.null(remove_basis)) {
    stop("basis and remove_basis options cannot both be set. Use one or the other.")
  }
  
  # read the data file if not already an mrs_data object
  if (class(metab)[[1]] == "mrs_data") {
    if (is.null(output_dir)) output_dir <- paste0("mrs_res_",
                                                  format(Sys.time(),
                                                         "%Y-%M-%d_%H%M%S"))
  } else if (dir.exists(metab)) {
    if (is.null(output_dir)) {
      output_dir <- sub("\\.", "_", basename(metab))
      output_dir <- paste0(output_dir, "_results")
    }
    metab <- read_ima_dyn_dir(metab) 
  } else {
    if (is.null(output_dir)) {
      output_dir <- sub("\\.", "_", basename(metab))
      output_dir <- paste0(output_dir, "_results")
    }
    metab <- read_mrs(metab)
  }
  
  if (verbose) cat(paste0("Output directory : ", output_dir, "\n"))
  
  # read the ref data file if not already an mrs_data object
  if (is.def(w_ref) & (class(w_ref)[[1]] != "mrs_data")) {
    if (dir.exists(w_ref)) {
      w_ref <- read_ima_dyn_dir(w_ref) 
    } else {
      w_ref <- read_mrs(w_ref)
    }
  }
  
  # check for GE style data
  if (identical(class(metab), c("list", "mrs_data"))) {
    x     <- metab
    metab <- x$metab
    w_ref <- x$ref
  }
  
  # create the output dir if it doesn't exist
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  if (is.null(tr)) tr <- tr(metab)
      
  # check we have what we need for concentration scaling
  if (is.null(te)) {
    te <- te(metab)
    if (is.null(te)) stop("Unable to determine echo-time from the data file, please pass as an argument.")
  }
  
  # combine coils if needed
  if (Ncoils(metab) > 1) {
    coil_comb_res <- comb_coils(metab, w_ref)
    if (is.null(w_ref)) {
      metab <- coil_comb_res
    } else {
      metab <- coil_comb_res$metab
      w_ref <- coil_comb_res$ref
    }
  }
  
  metab_pre_dfp_corr <- metab
  
  # dynamic frequency and phase correction
  if (dfp_corr & (Ndyns(metab) > 1)) {
    metab <- rats(metab, zero_freq_shift_t0 = TRUE, xlim = c(4, 1.8))
    metab_post_dfp_corr <- metab
  }
  
  metab <- mean_dyns(metab)
 
  if (!is.null(w_ref)) w_ref <- mean_dyns(w_ref)
  
  # eddy current correction
  if (ecc & (!is.null(w_ref))) metab <- ecc(metab, w_ref)
  
  # simulate a basis if needed
  if (is.null(basis)) {
    
    mol_list_chars <- c("m_cr_ch2", "ala", "asp", "cr", "gaba", "glc", "gln",
                        "gsh", "glu", "gpc", "ins", "lac", "lip09", "lip13a",
                        "lip13b", "lip20", "mm09", "mm12", "mm14", "mm17",
                        "mm20", "naa", "naag", "pch", "pcr", "sins", "tau")
    
    if (!is.null(append_basis)) mol_list_chars <- c(mol_list_chars,
                                                    append_basis)
    
    if (!is.null(remove_basis)) {
      inds <- which(mol_list_chars == remove_basis)
      mol_list_chars <- mol_list_chars[-inds]
    }
    
    mol_list <- get_mol_paras(mol_list_chars, ft = metab$ft)
    
    basis <- sim_basis(mol_list, acq_paras = metab, pul_seq = seq_press_ideal)
    
    if (verbose) print(basis)
  }
  
  # fitting
  fit_res <- fit_mrs(metab, basis = basis, opts = abfit_opts)
    
  phase_offset <- fit_res$res_tab$phase
  shift_offset <- fit_res$res_tab$shift
  
  # check for poor dynamics
  if (omit_bad_dynamics & (Ndyns(metab_pre_dfp_corr) > 1)) {
    if (dfp_corr) {
      dyn_data <- shift(phase(metab_post_dfp_corr, phase_offset), shift_offset)
    } else {
      dyn_data <- shift(phase(metab_pre_dfp_corr, phase_offset), shift_offset)
    }
    
    dyn_data_proc <- bc_poly(crop_spec(zf(lb(dyn_data, 2)), c(3.5, 1.8)), 2)
    
    peak_height <- spec_op(dyn_data_proc, operator = "max")
    
    peak_height <- peak_height / max(peak_height) * 100
    
    bad_shots      <- peak_height < 75
    bad_shots_n    <- sum(bad_shots)
    bad_shots_perc <- bad_shots_n / length(peak_height) * 100
    
    grDevices::png(file.path(output_dir, "drift_plot_peak_height.png"),
                   res = 2 * 72, height = 2 * 480, width = 2 * 480)
    graphics::image(dyn_data_proc)
    grDevices::dev.off()
    
    grDevices::pdf(file.path(output_dir, "dynamic_peak_height.pdf"))
    graphics::plot(peak_height, type = "l", ylim = c(0, 100),
                   ylab = "Max peak height (%)", xlab = "Dynamic")
    graphics::abline(h = 75, lty = 2)
    grDevices::dev.off()
    
    if (bad_shots_n > 0) {
      subset <- which(!bad_shots)
      
      cat(paste0(bad_shots_n, " bad shots (", round(bad_shots_perc), 
                 "%) detected.\n"))
      
      # remove bad shots and refit
      metab_pre_dfp_corr  <- get_dyns(metab_pre_dfp_corr, subset)
      
      if (dfp_corr) {
        metab_post_dfp_corr <- get_dyns(metab_post_dfp_corr, subset)
        metab <- metab_post_dfp_corr
      } else {
        metab <- metab_pre_dfp_corr
      }
      
      metab <- mean_dyns(metab)
  
      # eddy current correction
      if (ecc & (!is.null(w_ref))) metab <- ecc(metab, w_ref)
      
      # fitting
      if (verbose) cat("Refitting without bad shots.\n")
      fit_res <- fit_mrs(metab, basis = basis, opts = abfit_opts)
      
      phase_offset <- fit_res$res_tab$phase
      shift_offset <- fit_res$res_tab$shift
    } else {
      if (verbose) cat("No bad shots detected.\n")
    }
  }
  
  grDevices::pdf(file.path(output_dir, "fit_plot.pdf"))
  plot(fit_res)
  grDevices::dev.off()
  
  if (Ndyns(metab_pre_dfp_corr) > 1) {
    
    grDevices::png(file.path(output_dir, "drift_plot.png"), res = 2 * 72,
                   height = 2 * 480, width = 2 * 480)
    graphics::image(lb(phase(metab_pre_dfp_corr, phase_offset), 2),
                    xlim = c(4, 0.5))
    grDevices::dev.off()
    
    if (dfp_corr) {
      grDevices::png(file.path(output_dir, "drift_plot_dfp_corr.png"),
                     res = 2 * 72, height = 2 * 480, width = 2 * 480)
      graphics::image(lb(phase(metab_post_dfp_corr, phase_offset), 2),
                      xlim = c(4, 0.5))
      grDevices::dev.off()
    }
  }
  
  # output unscaled results
  utils::write.csv(fit_res$res_tab, file.path(output_dir,
                                              "fit_res_unscaled.csv"))
 
  # output ratio results if requested 
  if (!is.null(output_ratio)) {
    for (output_ratio_element in output_ratio) {
      fit_res_rat <- scale_amp_ratio(fit_res, output_ratio_element)
      
      file_out <- file.path(output_dir, paste0("fit_res_", output_ratio_element,
                                               "_ratio.csv"))
      
      utils::write.csv(fit_res_rat$res_tab, file_out)
    }
  }
  
  if (!is.null(w_ref)) {
    # assume 100% white matter if not told otherwise
    if (is.null(p_vols)) p_vols <- c(WM = 100, GM = 0, CSF = 0)
    fit_res_molal <- scale_amp_molal_pvc(fit_res, w_ref, p_vols, te, tr)
    file_out <- file.path(output_dir, "fit_res_molal_conc.csv")
    utils::write.csv(fit_res_molal$res_tab, file_out)
  }
  
  return(fit_res)
}

#' #' Standard SVS 1H brain analysis pipeline.
#' #' @param metab filepath or mrs_data object containing MRS metabolite data.
#' #' @param basis basis set object to use for analysis.
#' #' @param w_ref filepath or mrs_data object containing MRS water reference data.
#' #' @param mri_seg filepath or nifti object containing segmented MRI data.
#' #' @param mri filepath or nifti object containing anatomical MRI data.
#' #' @param output_dir directory path to output fitting results.
#' #' @param extra data.frame with one row containing additional information to be
#' #' attached to the fit results table.
#' #' @param decimate option to decimate the input data by a factor of two. The 
#' #' default value of NULL does not perform decimation unless the spectral width
#' #' is greater than 20 PPM.
#' #' @param rats_corr option to perform rats correction, defaults to TRUE.
#' #' @param ecc option to perform water reference based eddy current correction,
#' #' defaults to FALSE.
#' #' @param comb_dyns option to combine dynamic scans, defaults to TRUE.
#' #' @param hsvd_filt option to apply hsvd water removal, defaults to FALSE.
#' #' @param scale_amps option to scale metabolite amplitude estimates, defaults to
#' #' TRUE.
#' #' @param te metabolite mrs data echo time in seconds.
#' #' @param tr metabolite mrs data repetition time in seconds.
#' #' @param preproc_only only perform the preprocessing steps and omit fitting.
#' #' The preprocessed metabolite data will be returned in this case.
#' #' @return a fit_result or mrs_data object depending on the preproc_only option.
#' #' @param method analysis method to use, see fit_mrs help.
#' #' @param opts options to pass to the analysis method.
#' #' @export
#' svs_1h_brain_analysis <- function(metab, basis = NULL, w_ref = NULL,
#'                                   mri_seg = NULL, mri = NULL, output_dir = NULL,
#'                                   extra = NULL, decimate = NULL,
#'                                   rats_corr = TRUE, ecc = FALSE,
#'                                   comb_dyns = TRUE, hsvd_filt = FALSE,
#'                                   scale_amps = TRUE, te = NULL, tr = NULL,
#'                                   preproc_only = FALSE, method = "ABFIT",
#'                                   opts = NULL) {
#'   
#'   if (!preproc_only & is.null(basis)) stop("basis argument not set")
#'   
#'   # read the data file if not already an mrs_data object
#'   if (class(metab)[[1]] != "mrs_data") metab <- read_mrs(metab)
#'   
#'   # TODO check for combined metab + water ref. data, eg GE
#'   
#'   # remove fs_path types as they cause problems with comb_fit_list_fit_tables
#'   extra[] = lapply(extra, as.character)
#'   
#'   # create the output dir if it doesn't exist
#'   if (is.def(output_dir)) {
#'     if(!dir.exists(output_dir)) dir.create(output_dir)
#'   }
#'   
#'   # read the ref data file if not already an mrs_data object
#'   if (is.def(w_ref) & (class(w_ref)[[1]] != "mrs_data")) {
#'     w_ref <- read_mrs(w_ref)
#'   }
#'   
#'   if (is.def(mri) & (!("niftiImage" %in% class(mri)))) {
#'     mri <- readNifti(mri)
#'   }
#'   
#'   if (is.def(mri_seg) & (!("niftiImage" %in% class(mri_seg)))) {
#'     mri_seg <- readNifti(mri_seg)
#'   }
#'   
#'   if (is.def(mri)) RNifti::orientation(mri) <- "RAS"
#'   
#'   if (is.def(mri_seg)) RNifti::orientation(mri_seg) <- "RAS"
#'   
#'   # combine coils if needed
#'   if (Ncoils(metab) > 1) {
#'     coil_comb_res <- comb_coils(metab, w_ref)
#'     if (is.null(w_ref)) {
#'       metab <- coil_comb_res
#'     } else {
#'       metab <- coil_comb_res$metab
#'       w_ref <- coil_comb_res$ref
#'     }
#'   }
#'   
#'   # if the spectral width exceeds 20 PPM, then it should be ok
#'   # to decimate the signal to improve analysis speed
#'   if (is.null(decimate) & ((fs(metab) / metab$ft * 1e6) > 20)) {
#'     decimate <- TRUE
#'   } else {
#'     decimate <- FALSE
#'   }
#'   
#'   if (decimate) {
#'     metab <- decimate_mrs_fd(metab)
#'     if (!is.null(w_ref)) w_ref <- decimate_mrs_fd(w_ref)
#'   }
#'   
#'   # rats
#'   if (rats_corr & (Ndyns(metab)> 1)) metab <- rats(metab)$corrected
#'   
#'   # TODO plot of shifts?
#'   
#'   # combine dynamic scans
#'   if (comb_dyns) metab <- mean_dyns(metab)
#'   
#'   # eddy current correction
#'   if (ecc & (!is.null(w_ref))) metab <- ecc(metab, w_ref)
#'   
#'   # HSVD residual water removal
#'   if (hsvd_filt) metab <- hsvd_filt(metab)
#'   
#'   # end here if we're only interested in the preprocessed data
#'   if (preproc_only) return(metab)
#'   
#'   fit_res <- fit_mrs(metab, basis = basis, extra = extra, method = method,
#'                      opts = opts)
#'   
#'   # plot the fit result and output csv
#'   if (is.def(output_dir)) {
#'     grDevices::pdf(file.path(output_dir, "fit_plot.pdf"))
#'     plot(fit_res)
#'     grDevices::dev.off()
#'     utils::write.csv(fit_res$res_tab, file.path(output_dir, "fit_res.csv"))
#'   }
#'   
#'   # plot the voxel location on the mri
#'   if (is.def(mri) & is.def(output_dir)) {
#'     # generate the svs voi in the segmented image space
#'     voi <- get_svs_voi(metab, mri)
#'     plot_voi_overlay(mri, voi, file.path(output_dir, "vox_plot.png"))
#'   }
#'   
#'   # the calculate water reference frequency, linewidth and water suppression
#'   # efficiency
#'   if (is.def(w_ref)) {
#'     # phase correct based on the phase of the first data point
#'     w_ref      <- fp_phase_correct(w_ref)
#'     w_info_re  <- peak_info(w_ref, xlim = c(6, 4), mode = "real")
#'     w_info_mod <- peak_info(w_ref, xlim = c(6, 4), mode = "mod")
#'     
#'     fit_res$res_tab$water_freq_mod_ppm  <- w_info_mod$freq_ppm
#'     fit_res$res_tab$water_freq_real_ppm <- w_info_re$freq_ppm
#'     fit_res$res_tab$water_fwhm_mod_ppm  <- w_info_mod$fwhm_ppm
#'     fit_res$res_tab$water_fwhm_real_ppm <- w_info_re$fwhm_ppm
#'     
#'     # search for highest peak in the water suppressed data 0.1 ppm either side
#'     # of the water frequency
#'     water_freq <- w_info_mod$freq_ppm
#'     resid_water_peak_mod <- peak_info(metab, xlim = c(water_freq + 0.1,
#'                                              water_freq - 0.1), mode = "mod")
#'     fit_res$res_tab$water_sup_efficiency_percent <- resid_water_peak_mod$height /
#'                                                   w_info_mod$height * 100
#'   }
#'   
#'   if (scale_amps) {
#'     if (is.def(w_ref) & is.def(mri_seg)) {
#'       
#'       if (is.null(te)) stop("te not given, amplitude scaling failed")
#'         
#'       if (is.null(tr)) stop("tr not given, amplitude scaling failed")
#'       
#'       # generate the svs voi in the segmented image space
#'       voi_seg <- get_svs_voi(metab, mri_seg)
#'       
#'       # calculate partial volumes
#'       seg <- get_voi_seg(voi_seg, mri_seg)
#'       # do pvc
#'       fit_res <- scale_amp_molal_pvc(fit_res, w_ref, seg, te, tr)
#'       
#'       if (is.def(output_dir)) {
#'         plot_voi_overlay_seg(mri_seg, voi_seg, file.path(output_dir,
#'                                                      "vox_seg_plot.png"))
#'       }
#'       
#'     } else if (is.def(w_ref) & !is.def(mri_seg)) {
#'       # do straight w scaling default LCM style
#'       fit_res <- scale_amp_molar(fit_res, w_ref)
#'     } else {
#'       # scale to tCr
#'       fit_res <- scale_amp_ratio(fit_res, "tCr")
#'     }
#'   }
#'   
#'   return(fit_res)
#' }
#' 
#' # TODO
#' # what if metab fname results in metab + ref file, eg GE data?
#' # options to add, fit_opts, fit_method, format (GE, Siemens etc)
#' 
#' #' Batch interface to the standard SVS 1H brain analysis pipeline.
#' #' @param metab_list list of file paths or mrs_data objects containing MRS 
#' #' metabolite data.
#' #' @param w_ref_list list of file paths or mrs_data objects containing MRS 
#' #' water reference data.
#' #' @param mri_seg_list list of file paths or nifti objects containing segmented
#' #' MRI data.
#' #' @param mri_list list of file paths or nifti objects containing anatomical
#' #' MRI data.
#' #' @param output_dir_list list of directory paths to output fitting results.
#' #' @param extra a data frame with the same number of rows as metab_list,
#' #' containing additional information to be attached to the fit results table.
#' #' @param ... additional options to be passed to the svs_1h_brain_analysis
#' #' function.
#' #' @return a list of fit_result objects.
#' #' @export
#' svs_1h_brain_batch_analysis <- function(metab_list, w_ref_list = NULL,
#'                                         mri_seg_list = NULL, mri_list = NULL,
#'                                         output_dir_list = NULL, extra = NULL, 
#'                                         ...) {
#'   
#'   # check input is sensible
#'   metab_n <- length(metab_list)
#'   
#'   if (is.def(w_ref_list) & length(w_ref_list) != metab_n) {
#'     stop("Incorrect number of w_ref_list items.")
#'   }
#'   
#'   if (is.def(mri_list) & length(mri_list) != metab_n) {
#'     stop("Incorrect number of mri_list items.")
#'   }
#'   
#'   if (is.def(mri_seg_list) & length(mri_seg_list) != metab_n) {
#'     stop("Incorrect number of mri_seg_list items.")
#'   }
#'   
#'   if (is.def(output_dir_list) & length(output_dir_list) != metab_n) {
#'     stop("Incorrect number of output_dir_list items.")
#'   }
#'   
#'   if (is.def(extra) & nrow(extra) != metab_n) {
#'     stop("Incorrect number of rows in extra.")
#'   } else {
#'     # split into a list suitable for mapply
#'     extra <- split(extra, seq(nrow(extra)))
#'   }
#'   
#'   # check file paths exist TODO write a function to return bad paths
#'   
#'   if (is.null(w_ref_list)) w_ref_list <- vector("list", metab_n)
#'   
#'   if (is.null(mri_seg_list)) mri_seg_list <- vector("list", metab_n)
#'   
#'   if (is.null(mri_list)) mri_list <- vector("list", metab_n)
#'   
#'   if (is.null(output_dir_list)) output_dir_list <- vector("list", metab_n)
#'   
#'   if (is.null(extra)) extra <- vector("list", metab_n)
#'   
#'   fit_list <- mapply(svs_1h_brain_analysis, metab = metab_list,
#'                      w_ref = w_ref_list, mri_seg = mri_seg_list, mri = mri_list,
#'                      output_dir = output_dir_list, extra = extra,
#'                      MoreArgs = list(...), SIMPLIFY = FALSE, USE.NAMES = FALSE)
#'   
#'   return(fit_list)
#' }
