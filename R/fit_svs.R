#' Standard SVS 1H brain analysis pipeline.
#' 
#' Note this function is still under development and liable to changes.
#' 
#' @param metab path or mrs_data object containing MRS metabolite data.
#' @param w_ref path or mrs_data object containing MRS water reference data.
#' @param output_dir directory path to output fitting results.
#' @param basis precompiled basis set object to use for analysis.
#' @param p_vols a numeric vector of partial volumes expressed as percentages.
#' Defaults to 100% white matter. A voxel containing 100% gray matter tissue
#' would use : p_vols = c(WM = 0, GM = 100, CSF = 0).
#' @param format Override automatic data format detection. See format argument
#' in [read_mrs()] for permitted values.
#' @param pul_seq Pulse sequence to use for basis simulation. Can be one of the
#' following values : "press_ideal", "press_shaped", "steam" or "slaser".
#' @param TE1 PRESS or sLASER sequence timing parameter in seconds.
#' @param TE2 PRESS or sLASER sequence timing parameter in seconds.
#' @param TE3 sLASER sequence timing parameter in seconds.
#' @param TM STEAM mixing time parameter in seconds.
#' @param append_basis names of extra signals to add to the default basis. Eg 
#' append_basis = c("peth", "cit"). Cannot be used with precompiled basis sets.
#' @param remove_basis names of signals to remove from the basis. Cannot be used
#' with precompiled basis sets.
#' @param dfp_corr perform dynamic frequency and phase correction using the RATS
#' method.
#' @param TE metabolite mrs data echo time in seconds. If not supplied this will
#' be guessed from the metab data file.
#' @param TR metabolite mrs data repetition time in seconds. If not supplied
#' this will be guessed from the metab data file.
#' @param output_ratio optional string to specify a metabolite ratio to output.
#' Defaults to "tCr" and multiple metabolites may be specified for multiple
#' outputs. Set as NULL to omit.
#' @param ecc option to perform water reference based eddy current correction,
#' defaults to FALSE.
#' @param fit_opts options to pass to ABfit.
#' @param legacy_ws perform and output legacy water scaling compatible with
#' default LCModel and TARQUIN behaviour. See w_att and w_conc arguments to 
#' change the default assumptions. Default value is FALSE.
#' @param w_att water attenuation factor (default = 0.7) for legacy water 
#' scaling. Assumes water T2 of 80ms and a TE = 30 ms. exp(-30ms / 80ms) ~ 0.7.
#' @param w_conc assumed water concentration (default = 35880) for legacy water 
#' scaling. Default value corresponds to typical white matter. Set to 43300 for
#' gray matter, and 55556 for phantom measurements.
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
fit_svs <- function(metab, w_ref = NULL, output_dir = NULL, basis = NULL,
                    p_vols = NULL, format = NULL, pul_seq = NULL, TE1 = NULL,
                    TE2 = NULL, TE3 = NULL, TM = NULL, append_basis = NULL,
                    remove_basis = NULL, dfp_corr = TRUE, TE = NULL, TR = NULL,
                    output_ratio = "tCr", ecc = FALSE, fit_opts = NULL,
                    legacy_ws = FALSE, w_att = 0.7, w_conc = 35880,
                    verbose = FALSE) {
  
  argg <- c(as.list(environment()))
  
  # TODO
  # Implement and test format, pul_seq, TE1, TE2, TE3 and TM arguments.
  # Check reading Siemens dynamic data is sensible.
  # Realistic PRESS sim for B0 > 2.9T.
  # Add an option to select a subset of dynamics.
  
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
  # } else if (dir.exists(metab)) {
    # ima dyn directory
  #   if (is.null(output_dir)) {
  #    output_dir <- sub("\\.", "_", basename(metab))
  #    output_dir <- paste0(output_dir, "_results")
      # output_dir <- file.path(normalizePath(dirname(metab)), output_dir)
   # }
    #metab <- read_ima_dyn_dir(metab) 
  } else {
    metab <- read_mrs(metab, format = format)
    if (is.null(output_dir)) {
      output_dir <- sub("\\.", "_", basename(metab))
      output_dir <- paste0(output_dir, "_results")
    }
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
  
  # try get TE and TR from the data if not passed in
  if (is.null(TR)) TR <- tr(metab)
  if (is.null(TE)) TE <- te(metab)
  
  # check we have what's needed for standard water concentration scaling
  if (!is.null(w_ref)) {
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
    
    TE1 <- 0.0126
    TE2 <- te - TE1
    
    warning(paste0("Basis not specified, so assuming PRESS sequence with ",
                   "ideal pulses, ", "TE1 = ", TE1, "s, TE2 = ", TE2, "s."))
    
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
    
    basis <- sim_basis(mol_list, acq_paras = metab, pul_seq = seq_press_ideal,
                       TE1 = TE1, TE2 = TE2)
    
    if (verbose) print(basis)
  }
  
  if (is.null(fit_opts)) fit_opts <- abfit_reg_opts()
  
  # fitting
  fit_res <- fit_mrs(metab, basis = basis, opts = fit_opts)
    
  phase_offset <- fit_res$res_tab$phase
  shift_offset <- fit_res$res_tab$shift
  
  # check for poor dynamics
  omit_bad_dynamics <- FALSE # not an option yet
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
  
  # grDevices::pdf(file.path(output_dir, "fit_plot.pdf"))
  # plot(fit_res)
  # grDevices::dev.off()
  
  # if (Ndyns(metab_pre_dfp_corr) > 1) {
  #   
  #   grDevices::png(file.path(output_dir, "drift_plot.png"), res = 2 * 72,
  #                  height = 2 * 480, width = 2 * 480)
  #   graphics::image(lb(phase(metab_pre_dfp_corr, phase_offset), 2),
  #                   xlim = c(4, 0.5))
  #   grDevices::dev.off()
  #   
  #   if (dfp_corr) {
  #     grDevices::png(file.path(output_dir, "drift_plot_dfp_corr.png"),
  #                    res = 2 * 72, height = 2 * 480, width = 2 * 480)
  #     graphics::image(lb(phase(metab_post_dfp_corr, phase_offset), 2),
  #                     xlim = c(4, 0.5))
  #     grDevices::dev.off()
  #   }
  # }
  
  # output unscaled results
  res_tab_unscaled <- fit_res$res_tab
  utils::write.csv(res_tab_unscaled, file.path(output_dir,
                                               "fit_res_unscaled.csv"))
 
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
  
  if (!is.null(w_ref)) {
    # assume 100% white matter if not told otherwise
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
  
  results <- list(fit_res = fit_res, argg = argg,
                  w_ref_available = w_ref_available, w_ref = w_ref,
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