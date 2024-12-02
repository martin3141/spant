#' Standard SVS 1H brain analysis pipeline.
#' 
#' Note this function is still under development and liable to changes.
#' 
#' @param metab path or mrs_data object containing MRS metabolite data.
#' @param w_ref path or mrs_data object containing MRS water reference data.
#' @param output_dir directory path to output fitting results.
#' @param external_basis precompiled basis set object to use for analysis.
#' @param p_vols a numeric vector of partial volumes expressed as percentages.
#' Defaults to 100% white matter. A voxel containing 100% gray matter tissue
#' would use : p_vols = c(WM = 0, GM = 100, CSF = 0).
#' @param format Override automatic data format detection. See format argument
#' in [read_mrs()] for permitted values.
#' @param pul_seq Pulse sequence to use for basis simulation. Can be one of the
#' following values : "press", "press_ideal", "press_shaped", "steam" or
#' "slaser". If "press" then "press_ideal" will be assumed unless the magnetic
#' field is stronger that 2.8 Tesla, "press_shaped" will be assumed for 2.9 
#' Tesla and above. 
#' @param TE metabolite mrs data echo time in seconds. If not supplied this will
#' be guessed from the metab data file.
#' @param TR metabolite mrs data repetition time in seconds. If not supplied
#' this will be guessed from the metab data file.
#' @param TE1 PRESS or sLASER sequence timing parameter in seconds.
#' @param TE2 PRESS or sLASER sequence timing parameter in seconds.
#' @param TE3 sLASER sequence timing parameter in seconds.
#' @param TM STEAM mixing time parameter in seconds.
#' @param append_basis names of extra signals to add to the default basis. Eg 
#' append_basis = c("peth", "cit"). Cannot be used with precompiled basis sets.
#' @param remove_basis grep expression to match names of signals to remove from
#' the basis. For example: use "*" to remove all signals, "^mm|^lip" to remove
#' all macromolecular and lipid signals, "^lac" to remove lactate. This operation
#' is performed before signals are added with append_basis. Cannot be used with
#' precompiled basis sets.
#' @param dfp_corr perform dynamic frequency and phase correction using the RATS
#' method.
#' @param output_ratio optional string to specify a metabolite ratio to output.
#' Defaults to "tCr" and multiple metabolites may be specified for multiple
#' outputs. Set as NULL to omit.
#' @param ecc option to perform water reference based eddy current correction,
#' defaults to FALSE.
#' @param fit_opts options to pass to ABfit.
#' @param fit_subset specify a subset of dynamics to analyse, for example
#' 1:16 would only fit the first 16 dynamic scans.
#' @param legacy_ws perform and output legacy water scaling compatible with
#' default LCModel and TARQUIN behaviour. See w_att and w_conc arguments to 
#' change the default assumptions. Default value is FALSE.
#' @param w_att water attenuation factor (default = 0.7) for legacy water 
#' scaling. Assumes water T2 of 80ms and a TE = 30 ms. exp(-30ms / 80ms) ~ 0.7.
#' @param w_conc assumed water concentration (default = 35880) for legacy water 
#' scaling. Default value corresponds to typical white matter. Set to 43300 for
#' gray matter, and 55556 for phantom measurements.
#' @param use_basis_cache Pre-cache basis sets to reduce analysis speed. Can be
#' one of the following : "auto", "all" or "none". The default value of "auto" 
#' will only use the cache for 3T PRESS - which generally requires more detailed
#' simulation due to high CSD.
#' @param verbose output potentially useful information.
#' @examples
#' metab <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
#'                      package = "spant")
#' w_ref <- system.file("extdata", "philips_spar_sdat_W.SDAT",
#'                      package = "spant")
#' \dontrun{
#' fit_result <- svs_1h_brain_analysis(metab, w_ref, "fit_res_dir")
#' }
#' @export
fit_svs <- function(metab, w_ref = NULL, output_dir = NULL,
                    external_basis = NULL, p_vols = NULL, format = NULL,
                    pul_seq = NULL, TE = NULL, TR = NULL, TE1 = NULL,
                    TE2 = NULL, TE3 = NULL, TM = NULL, append_basis = NULL,
                    remove_basis = NULL, dfp_corr = TRUE, output_ratio = "tCr",
                    ecc = FALSE, fit_opts = NULL, fit_subset = NULL, 
                    legacy_ws = FALSE, w_att = 0.7, w_conc = 35880,
                    use_basis_cache = "auto", verbose = FALSE) {
  
  argg <- c(as.list(environment()))
  
  if (!is.null(external_basis) & !is.null(append_basis)) {
    stop("external_basis and append_basis options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(external_basis) & !is.null(remove_basis)) {
    stop("external_basis and remove_basis options cannot both be set. Use one or the other.")
  }
  
  # read the data file if not already an mrs_data object
  if (class(metab)[[1]] == "mrs_data") {
    if (is.null(output_dir)) {
      output_dir <- paste0("mrs_res_", format(Sys.time(), "%Y-%M-%d_%H%M%S"))
    }
  } else {
    metab_path <- metab
    metab      <- read_mrs(metab, format = format)
    if (is.null(output_dir)) {
      output_dir <- sub("\\.", "_", basename(metab_path))
      output_dir <- paste0(output_dir, "_results")
    }
  }
  
  if (verbose) cat(paste0("Output directory : ", output_dir, "\n"))
  
  # read the ref data file if not already an mrs_data object
  if (is.def(w_ref) & (class(w_ref)[[1]] != "mrs_data")) {
    w_ref <- read_mrs(w_ref, format = format)
  }
  
  # check for GE style data with metabolite and water reference data
  # contained in the same file
  if (identical(class(metab), c("list", "mrs_data"))) {
    x     <- metab
    metab <- x$metab
    if (is.null(w_ref)) {
      if (verbose) cat("Using water reference data within metab file.\n")
      w_ref <- x$ref
    } else {
      if (verbose) cat("Overriding water reference data within metab file.\n")
    } 
  }
  
  # by this point we know if we have water reference data available
  if (is.null(w_ref)) {
    if (verbose) cat("Water reference is not available.\n")
    w_ref_available = FALSE
  } else {
    if (verbose) cat("Water reference data is available.\n")
    w_ref_available = TRUE
  }
  
  # create the output dir if it doesn't exist
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  # try to get TE and TR parameters from the data if not passed in
  if (is.null(TR)) TR <- tr(metab)
  if (is.null(TE)) TE <- te(metab)
  
  # check we have what's needed for standard water concentration scaling
  if (w_ref_available) {
    if (is.null(TR)) stop("Please provide seqeuence TR argument for water concentration scaling.")
    if (is.null(TE)) stop("Please provide seqeuence TE argument for water concentration scaling.")
  }
  
  # combine coils if needed
  if (Ncoils(metab) > 1) {
    coil_comb_res <- comb_coils_svs_gls(metab, w_ref) # may not want to use
    if (is.null(w_ref)) {                             # w_ref in some cases?
      metab <- coil_comb_res
    } else {
      metab <- coil_comb_res$metab
      w_ref <- coil_comb_res$ref
    }
  }
  
  # extract a subset of dynamic scans if specified
  if (!is.null(fit_subset)) {
    metab <- get_dyns(metab, fit_subset) 
  }
  
  # dynamic frequency and phase correction
  metab_pre_dfp_corr <- metab
  if (dfp_corr & (Ndyns(metab) > 1)) {
    metab <- rats(metab, zero_freq_shift_t0 = TRUE, xlim = c(4, 1.8))
    metab_post_dfp_corr <- metab
  }
  
  # take the mean of the metabolite data
  metab <- mean_dyns(metab)
 
  # take the mean of the water reference data
  if (w_ref_available) w_ref <- mean_dyns(w_ref)
  
  # calculate the water suppression efficiency
  # the ratio of the residual water peak height
  # relative to the height of the unsuppressed water signal
  
  if (w_ref_available) {
    w_ref_peak   <- peak_info(w_ref, xlim = c(7, 3), mode = "mod")
    w_ref_height <- w_ref_peak$height[1]
    w_ref_freq   <- w_ref_peak$freq_ppm[1]
    x_range      <- c(w_ref_freq - 0.2, w_ref_freq + 0.2)
    metab_water_height <- spec_op(metab, xlim = x_range, operator = "max",
                                  mode = "mod")[1]
    ws_efficiency <- metab_water_height / w_ref_height * 100
  }
  
  # eddy current correction
  if (ecc & w_ref_available) metab <- ecc(metab, w_ref)
  
  # simulate a basis if needed
  if (is.null(external_basis)) {
    
    if (is.null(TE)) stop("Could not determine the sequence echo time. Please provide the TE argument.")
    
    # list of "standard" signals to include in the basis set
    mol_list_chars <- c("m_cr_ch2", "ala", "asp", "cr", "gaba", "glc", "gln",
                        "gsh", "glu", "gpc", "ins", "lac", "lip09", "lip13a",
                        "lip13b", "lip20", "mm09", "mm12", "mm14", "mm17",
                        "mm20", "naa", "naag", "pch", "pcr", "sins", "tau")
    
    # option to remove signals
    if (!is.null(remove_basis)) {
        inds <- grep(remove_basis, mol_list_chars)
        if (length(inds) == 0) stop("No signals matching remove_basis found.")
        mol_list_chars <- mol_list_chars[-inds]
    }
    
    # option to append signals
    if (!is.null(append_basis)) mol_list_chars <- c(mol_list_chars,
                                                    append_basis)
    
    # probably set remove_basis to * and forgot to use append_basis
    if (is.null(mol_list_chars)) stop("No basis signals named for simulation.")
    
    # get the parameters
    mol_list <- get_mol_paras(mol_list_chars, ft = metab$ft)
    
    # check parameters are consistent and infer any missing values
    sim_paras <- check_sim_paras(pul_seq, metab, TE1, TE2, TE3, TE, TM)
    
    # determine if basis caching should be used
    if (use_basis_cache == "always") {
      use_basis_cache = TRUE
    } else if (use_basis_cache == "never") {
      use_basis_cache = FALSE
    } else if (use_basis_cache == "auto") {
      if (sim_paras$pul_seq == "press_shaped") {
        use_basis_cache = TRUE
      } else {
        use_basis_cache = FALSE
      }
    } else {
      stop("incorrect value for use_basis_cache, should be: 'auto', 'never' or 'always'") 
    }
    
    if (sim_paras$pul_seq == "press_ideal") {
      if (verbose) cat("Simulating ideal PRESS sequence.\n")
      basis <- sim_basis(mol_list, acq_paras = metab,
                         pul_seq = seq_press_ideal, TE1 = sim_paras$TE1,
                         TE2 = sim_paras$TE2, use_basis_cache = use_basis_cache,
                         verbose = verbose)
    } else if (sim_paras$pul_seq == "press_shaped") {
      if (verbose) cat("Simulating shaped PRESS sequence.\n")
      pulse_file <- system.file("extdata", "press_refocus.pta",
                                package = "spant")
      basis <- sim_basis(mol_list, acq_paras = metab,
                         pul_seq = seq_press_2d_shaped, TE1 = sim_paras$TE1,
                         TE2 = sim_paras$TE2, use_basis_cache = use_basis_cache,
                         verbose = verbose, pulse_file = pulse_file,
                         pulse_dur = 5e-3, pulse_file_format = "pta")
    } else if (sim_paras$pul_seq == "steam") {
      if (verbose) cat("Simulating STEAM sequence.\n")
      basis <- sim_basis(mol_list, acq_paras = metab,
                         pul_seq = seq_steam_ideal_cof, TE = sim_paras$TE,
                         TM = sim_paras$TM, use_basis_cache = use_basis_cache,
                         verbose = verbose)
    } else if (sim_paras$pul_seq == "slaser") {
      if (verbose) cat("Simulating sLASER sequence.\n")
      basis <- sim_basis(mol_list, acq_paras = metab,
                         pul_seq = seq_slaser_ideal, TE1 = sim_paras$TE1,
                         TE2 = sim_paras$TE2, TE3 = sim_paras$TE3,
                         use_basis_cache = use_basis_cache, verbose = verbose)
    }
    if (verbose) print(basis)
  } else {
    basis <- external_basis
  }
  
  if (is.null(fit_opts)) fit_opts <- abfit_reg_opts()
  
  # fitting
  fit_res <- fit_mrs(metab, basis = basis, opts = fit_opts)
    
  phase_offset <- fit_res$res_tab$phase
  shift_offset <- fit_res$res_tab$shift
  
  if (w_ref_available) fit_res$res_tab$ws_eff <- ws_efficiency
  
  # check for poor dynamics - needs work, so not used at the moment!
  # omit_bad_dynamics <- FALSE # not an option yet
  # if (omit_bad_dynamics & (Ndyns(metab_pre_dfp_corr) > 1)) {
  #   if (dfp_corr) {
  #     dyn_data <- shift(phase(metab_post_dfp_corr, phase_offset), shift_offset)
  #   } else {
  #     dyn_data <- shift(phase(metab_pre_dfp_corr, phase_offset), shift_offset)
  #   }
  #   
  #   dyn_data_proc <- bc_poly(crop_spec(zf(lb(dyn_data, 2)), c(3.5, 1.8)), 2)
  #   
  #   peak_height <- spec_op(dyn_data_proc, operator = "max")
  #   
  #   peak_height <- peak_height / max(peak_height) * 100
  #   
  #   bad_shots      <- peak_height < 75
  #   bad_shots_n    <- sum(bad_shots)
  #   bad_shots_perc <- bad_shots_n / length(peak_height) * 100
  #   
  #   grDevices::png(file.path(output_dir, "drift_plot_peak_height.png"),
  #                  res = 2 * 72, height = 2 * 480, width = 2 * 480)
  #   graphics::image(dyn_data_proc)
  #   grDevices::dev.off()
  #   
  #   grDevices::pdf(file.path(output_dir, "dynamic_peak_height.pdf"))
  #   graphics::plot(peak_height, type = "l", ylim = c(0, 100),
  #                  ylab = "Max peak height (%)", xlab = "Dynamic")
  #   graphics::abline(h = 75, lty = 2)
  #   grDevices::dev.off()
  #   
  #   if (bad_shots_n > 0) {
  #     subset <- which(!bad_shots)
  #     
  #     cat(paste0(bad_shots_n, " bad shots (", round(bad_shots_perc), 
  #                "%) detected.\n"))
  #     
  #     # remove bad shots and refit
  #     metab_pre_dfp_corr  <- get_dyns(metab_pre_dfp_corr, subset)
  #     
  #     if (dfp_corr) {
  #       metab_post_dfp_corr <- get_dyns(metab_post_dfp_corr, subset)
  #       metab <- metab_post_dfp_corr
  #     } else {
  #       metab <- metab_pre_dfp_corr
  #     }
  #     
  #     metab <- mean_dyns(metab)
  # 
  #     # eddy current correction
  #     if (ecc & (!is.null(w_ref))) metab <- ecc(metab, w_ref)
  #     
  #     # fitting
  #     if (verbose) cat("Refitting without bad shots.\n")
  #     fit_res <- fit_mrs(metab, basis = basis, opts = abfit_opts)
  #     
  #     phase_offset <- fit_res$res_tab$phase
  #     shift_offset <- fit_res$res_tab$shift
  #   } else {
  #     if (verbose) cat("No bad shots detected.\n")
  #   }
  # }
  
  # keep unscaled results
  res_tab_unscaled <- fit_res$res_tab
 
  # output ratio results if requested 
  if (!is.null(output_ratio)) {
    for (output_ratio_element in output_ratio) {
      fit_res_rat <- scale_amp_ratio(fit_res, output_ratio_element)
      
      res_tab_ratio <- fit_res_rat$res_tab
      file_out <- file.path(output_dir, paste0("fit_res_", output_ratio_element,
                                               "_ratio.csv"))
      
      utils::write.csv(res_tab_ratio, file_out)
    }
  } else {
    res_tab_ratio <- NULL  
  }
  
  if (w_ref_available) {
    # assume 100% white matter unless told otherwise
    if (is.null(p_vols)) p_vols <- c(WM = 100, GM = 0, CSF = 0)
    fit_res_molal <- scale_amp_molal_pvc(fit_res, w_ref, p_vols, TE, TR)
    res_tab_molal <- fit_res_molal$res_tab
    file_out <- file.path(output_dir, "fit_res_molal_conc.csv")
    utils::write.csv(res_tab_molal, file_out)
    if (legacy_ws) {
      fit_res_legacy <- scale_amp_legacy(fit_res, w_ref, w_att, w_conc)
      res_tab_legacy <- fit_res_legacy$res_tab
      file_out <- file.path(output_dir, "fit_res_legacy_conc.csv")
      utils::write.csv(res_tab_legacy, file_out)
    } else {
      res_tab_legacy <- NULL
    }
  } else {
    res_tab_legacy <- NULL 
    res_tab_molal  <- NULL 
  }
  
  # add water amplitude and PVC info to the unscaled output
  res_tab_unscaled <- cbind(res_tab_unscaled, w_amp = res_tab_molal$w_amp,
                            GM_vol = res_tab_molal$GM_vol,
                            WM_vol = res_tab_molal$WM_vol,
                            CSF_vol = res_tab_molal$CSF_vol,
                            GM_frac = res_tab_molal$GM_frac)
  
  utils::write.csv(res_tab_unscaled, file.path(output_dir,
                                               "fit_res_unscaled.csv"))
  
  # prepare dynamic data for plotting
  if (Ndyns(metab_pre_dfp_corr) > 1) {
    # phase according to the fit results
    dyn_data_uncorr <- phase(metab_pre_dfp_corr, fit_res$res_tab$phase)
    # correct chem. shift scale according to the fit results
    dyn_data_uncorr <- shift(dyn_data_uncorr, fit_res$res_tab$shift,
                             units = "ppm")
    # add 2 Hz LB
    dyn_data_uncorr <- lb(dyn_data_uncorr, 2)
    
    if (dfp_corr) {
      # phase according to the fit results
      dyn_data_corr <- phase(metab_post_dfp_corr, fit_res$res_tab$phase)
      # correct chem. shift scale according to the fit results
      dyn_data_corr <- shift(dyn_data_corr, fit_res$res_tab$shift,
                             units = "ppm")
      # add 2 Hz LB
      dyn_data_corr <- lb(dyn_data_corr, 2)
    } else {
      dyn_data_corr <- NULL
    }
  } else {
    dyn_data_uncorr <- NULL
    dyn_data_corr   <- NULL
  }
  
  # data needed to produce the output html report
  results <- list(fit_res = fit_res, argg = argg,
                  w_ref_available = w_ref_available,
                  w_ref = w_ref,
                  res_tab_unscaled = res_tab_unscaled,
                  res_tab_ratio = res_tab_ratio,
                  res_tab_legacy = res_tab_legacy,
                  res_tab_molal = res_tab_molal,
                  dyn_data_uncorr = dyn_data_uncorr,
                  dyn_data_corr = dyn_data_corr) 
  
  rmd_file <- system.file("rmd", "svs_report.Rmd", package = "spant")
  
  rmd_out_f <- file.path(tools::file_path_as_absolute(output_dir), "report")
  
  rmarkdown::render(rmd_file, params = results, output_file = rmd_out_f,
                    quiet = !verbose)
  
  return(fit_res)
}

check_sim_paras <- function(pul_seq, metab, TE1, TE2, TE3, TE, TM) {
  
  # try to guess the pulse sequence from the MRS data if not specified
  if (is.null(pul_seq)) {
    if (!is.null(metab$meta$PulseSequenceType)){
      pul_seq <- metab$meta$PulseSequenceType
    } else {
      warning(paste0("Could not determine the pulse sequence, so assuming ",
                     "PRESS. Provide the pul_seq argument to stop this ",
                     "warning."))
      pul_seq <- "press"
    }
  }
  
  # force lowercase
  pul_seq <- tolower(pul_seq)
  
  if (pul_seq == "press") {
    B0 <- round(metab$ft / 42.58e6, 1)
    if (B0 > 2.8) {
      pul_seq <- "press_shaped"
    } else {
      pul_seq <- "press_ideal"
    }
  }
  
  if (pul_seq == "press_ideal") {
    if (is.null(TE1)) {
      TE1 <- 0.0126
      # warning("TE1 assumed to be 0.0126s. Provide the TE1 argument to stop this warning.")
    }
    if (is.null(TE2)) TE2 <- TE - TE1
    if ((TE1 + TE2) != TE) warning("TE, TE1 and TE2 do not match.")
    return(list(pul_seq = pul_seq, TE1 = TE1, TE2 = TE2))
  } else if (pul_seq == "press_shaped") {
    if (is.null(TE1)) {
      TE1 <- 0.0126
      # warning("TE1 assumed to be 0.0126s.")
    }
    if (is.null(TE2)) TE2 <- TE - TE1
    if ((TE1 + TE2) != TE) warning("TE, TE1 and TE2 do not match.")
    return(list(pul_seq = pul_seq, TE1 = TE1, TE2 = TE2))
  } else if (pul_seq == "steam") {
    if (is.null(TM)) {
      TM <- 0.01
      warning("TM assumed to be 0.01s.")
    }
    return(list(pul_seq = pul_seq, TE = TE, TM = TM))
  } else if (pul_seq == "slaser") {
    if (is.null(TE1)) {
      if (!is.null(metab$meta$TE1)) {
        TE1 <- metab$meta$TE1
      } else {
        TE1 <- 0.0008
        warning("TE1 assumed to be 0.0008s.")
      }
    }
    if (is.null(TE2)) {
      if (!is.null(metab$meta$TE2)) {
        TE2 <- metab$meta$TE2
      } else {
        TE2 <- 0.0011
        warning("TE2 assumed to be 0.0011s.")
      }
    }
    if (is.null(TE3)) {
      if (!is.null(metab$meta$TE3)) {
        TE3 <- metab$meta$TE3
      } else {
        TE3 <- 0.0009
        warning("TE3 assumed to be 0.0009s.")
      }
    }
    if ((TE1 + TE2 + TE3) != TE) warning("TE, TE1, TE2 and TE3 do not match.")
    return(list(pul_seq = pul_seq, TE1 = TE1, TE2 = TE2, TE3 = TE3))
  } else {
    stop(paste0("pul_seq not supported : ", pul_seq))
  }
}
