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
                    p_vols = NULL, append_basis = NULL, remove_basis = NULL,
                    dfp_corr = FALSE, omit_bad_dynamics = FALSE, te = NULL,
                    tr = NULL, output_ratio = "tCr", ecc = FALSE,
                    abfit_opts = NULL, verbose = FALSE) {
  
  # TODO
  # Auto sequence detection.
  # Realistic PRESS sim for B0 > 2.9T.
  
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
    coil_comb_res <- comb_coils_svs_gls(metab, w_ref)
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