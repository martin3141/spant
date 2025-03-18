#' Standard SVS 1H brain analysis pipeline.
#' 
#' Note this function is still under development and liable to changes.
#' 
#' @param input path or mrs_data object containing MRS data.
#' @param w_ref path or mrs_data object containing MRS water reference data.
#' @param output_dir directory path to output fitting results.
#' @param mri filepath or nifti object containing anatomical MRI data.
#' @param mri_seg filepath or nifti object containing segmented MRI data.
#' @param external_basis precompiled basis set object to use for analysis.
#' @param append_external_basis append the external basis with the internally
#' generated one. Useful for adding experimentally acquired baseline signals to
#' internally simulated basis sets. Defaults to FALSE - meaning only signals 
#' from the external basis will be used in analysis.
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
#' @param pre_align perform simple frequency alignment to known reference peaks.
#' @param dfp_corr perform dynamic frequency and phase correction using the RATS
#' method.
#' @param output_ratio optional string to specify a metabolite ratio to output.
#' Defaults to "tCr". Multiple metabolites may be specified for multiple
#' outputs. Set to NA to omit.
#' @param ecc option to perform water reference based eddy current correction,
#' defaults to FALSE.
#' @param hsvd_width set the width of the HSVD filter in Hz. Note the applied
#' width is between -width and +width Hz, with 0 Hz being defined at the centre
#' of the spectral width. Default is disabled (set to NULL), 30 Hz is a
#' reasonable value.
#' @param decimate option on decimate the data by a factor of 2 before analysis.
#' Defaults to FALSE.
#' @param trunc_fid_pts number of points to truncate the input data by in the
#' time-domain. E.g. setting to 1024 will ensure data with more time-domain
#' points will be truncated to a length of 1024. Defaults to NULL, where
#' truncation is not performed.
#' @param fit_method can be "ABFIT-REG" or "LCMODEL. Defaults to "ABFIT-REG".
#' @param fit_opts options to pass to the fitting method.
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
#' @param summary_measures output an additional table with a subset of
#' metabolite levels, eg c("tNAA", "tNAA/tCr", "tNAA/tCho", "Lac/tNAA").
#' @param dyn_av_block_size perform temporal averaging with the specified block
#' size. Defaults to NULL, eg average across all dynamic scans.
#' @param dyn_av_scheme a numeric vector of sequential integers (starting at 1),
#' with the same length as the number of dynamic scans in the metabolite data.
#' For example: c(1, 1, 2, 1, 1, 3, 1, 1).
#' @param dyn_av_scheme_file a file path containing a single column of
#' sequential integers (starting at 1) with the same length as the number of
#' dynamic scans in the metabolite data. File may be formatted as .xlsx, .xls,
#' text or csv format.
#' @param lcm_bin_path set the path to LCModel binary.
#' @param plot_ppm_xlim plotting ppm axis limits in the html results.
#' results.
#' @param extra_output write extra output files for generating custom plots.
#' Defaults to FALSE.
#' @param verbose output potentially useful information.
#' @examples
#' metab <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
#'                      package = "spant")
#' w_ref <- system.file("extdata", "philips_spar_sdat_W.SDAT",
#'                      package = "spant")
#' out_dir <- file.path("~", "fit_svs_result")
#' \dontrun{
#' fit_result <- fit_svs(metab, w_ref, out_dir)
#' }
#' @export
fit_svs <- function(input, w_ref = NULL, output_dir = NULL, mri = NULL,
                    mri_seg = NULL, external_basis = NULL,
                    append_external_basis = FALSE, p_vols = NULL,
                    format = NULL, pul_seq = NULL, TE = NULL, TR = NULL,
                    TE1 = NULL, TE2 = NULL, TE3 = NULL, TM = NULL,
                    append_basis = NULL, remove_basis = NULL, pre_align = TRUE,
                    dfp_corr = TRUE, output_ratio = NULL, ecc = FALSE,
                    hsvd_width = NULL, decimate = FALSE, trunc_fid_pts = NULL,
                    fit_method = NULL, fit_opts = NULL, fit_subset = NULL,
                    legacy_ws = FALSE, w_att = 0.7, w_conc = 35880,
                    use_basis_cache = "auto", summary_measures = NULL,
                    dyn_av_block_size = NULL, dyn_av_scheme = NULL,
                    dyn_av_scheme_file = NULL, lcm_bin_path = NULL,
                    plot_ppm_xlim = NULL, extra_output = FALSE,
                    verbose = FALSE) {
  
  argg  <- c(as.list(environment()))
  
  metab <- input
  
  if (!is.null(dyn_av_scheme) & !is.null(dyn_av_scheme_file)) {
    print(dyn_av_scheme)
    print(dyn_av_scheme_file)
    stop("dyn_av_scheme and dyn_av_scheme_file options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(external_basis) & !is.null(append_basis) & !append_external_basis) {
    stop("external_basis and append_basis options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(external_basis) & !is.null(remove_basis) & !append_external_basis) {
    stop("external_basis and remove_basis options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(dyn_av_block_size) & !is.null(dyn_av_scheme)) {
    stop("dyn_av_block_size and dyn_av_scheme options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(dyn_av_block_size) & !is.null(dyn_av_scheme_file)) {
    stop("dyn_av_block_size and dyn_av_scheme_file options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(mri_seg) & !is.null(p_vols)) {
    warning("mri_seg and pvols options have both been set. Only p_vols will be used for partial volume correction calculation.")
  }
  
  # read the data file if not already an mrs_data object
  if (class(metab)[[1]] == "mrs_data") {
    if (is.null(output_dir)) {
      output_dir <- paste0("mrs_res_", format(Sys.time(), "%Y-%M-%d_%H%M%S"))
    }
  } else {
    metab_path <- metab
    if (verbose) cat(paste0("Reading MRS input data : ", metab,"\n"))
    metab      <- read_mrs(metab, format = format)
    if (is.null(output_dir)) {
      output_dir <- gsub("\\.", "_", basename(metab_path))
      output_dir <- gsub("#", "_", output_dir)
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
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  } else {
    warning(paste0("Output directory already exists : ", output_dir))
  }
  
  # check the mri data if specified 
  if (is.def(mri) & (!("niftiImage" %in% class(mri)))) {
    mri <- readNifti(mri)
  }
  
  # reorientate mri
  if (is.def(mri)) RNifti::orientation(mri) <- "RAS"

  # check the mri_seg data if specified 
  if (is.def(mri_seg) & (!("niftiImage" %in% class(mri_seg)))) {
    mri_seg <- readNifti(mri_seg)
  }
  
  # reorientate mri_seg
  if (is.def(mri_seg)) RNifti::orientation(mri_seg) <- "RAS"
  
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
  
  # decimate if specified
  if (decimate) {
    if (verbose) cat("Decimating data.\n")
    metab <- decimate_mrs_fd(metab)
  }
  
  # extract a subset of dynamic scans if specified
  if (!is.null(fit_subset)) {
    metab <- get_dyns(metab, fit_subset) 
  }
  
  metab_pre_dfp_corr <- metab
  
  # pre-alignment
  if (pre_align) {
    metab <- align(metab, c(2.01, 3.03, 3.22), max_shift = 40)
    if (Ndyns(metab) > 1) metab_post_dfp_corr <- metab
  }
  
  # rats correction
  if (dfp_corr & (Ndyns(metab) > 1)) {
    metab <- rats(metab, zero_freq_shift_t0 = TRUE, xlim = c(4, 1.8))
    metab_post_dfp_corr <- metab
  }
  
  if (!exists("metab_post_dfp_corr")) metab_post_dfp_corr <- NULL
  
  # read the dynamic averaging scheme from a file if specified
  if (!is.null(dyn_av_scheme_file)) {
    file_ext <- tools::file_ext(dyn_av_scheme_file)
    if (file_ext == "xls" | file_ext == "xlsx") {
      scheme_tab <- suppressMessages(readxl::read_excel(dyn_av_scheme_file,
                                     col_names = FALSE, col_types = "numeric"))
      dyn_av_scheme <- as.integer(scheme_tab[[1]])
    } else {
      dyn_av_scheme <- as.integer(utils::read.csv(dyn_av_scheme_file,
                                                  header = FALSE)[[1]])
    }
  }
  
  # take the mean of the metabolite data
  if (!is.null(dyn_av_block_size)) {
    metab <- mean_dyn_blocks(metab, dyn_av_block_size)
  } else if (!is.null(dyn_av_scheme)) {
    if (length(dyn_av_scheme) != Ndyns(metab)) {
      stop(paste0("dyn_av_scheme is the wrong length. Currently : ",
                  length(dyn_av_scheme),", should be : ", Ndyns(metab)))
    }
    dyn_av_scheme <- as.integer(dyn_av_scheme)
    max_dyn       <- max(dyn_av_scheme)
    metab_list    <- vector("list", length = max_dyn)
    for (n in 1:max_dyn) {
      subset <- which(dyn_av_scheme == n) 
      metab_list[[n]] <- mean_dyns(get_dyns(metab, subset))
    }
    metab <- append_dyns(metab_list)
  } else {
    metab <- mean_dyns(metab)
  }
 
  # take the mean of the water reference data
  # if (w_ref_available) w_ref <- mean_dyns(w_ref)
  
  # extract the first dynamic of the water reference data
  # better for Dinesh sLASER where the water peak can shift before and after
  # the metabolite data collection
  if (w_ref_available) w_ref <- get_dyns(w_ref, 1)
  
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
  if (is.null(external_basis) | append_external_basis) {
    
    if (is.null(TE)) stop("Could not determine the sequence echo time. Please provide the TE argument.")
    
    # list of "standard" signals to include in the basis set
    mol_list_chars <- c("m_cr_ch2", "ala", "asp", "cr", "gaba", "glc", "gln",
                        "gsh", "glu", "gpc", "ins", "lac", "lip09", "lip13a",
                        "lip13b", "lip20", "mm09", "mm12", "mm14", "mm17",
                        "mm20", "naa", "naag", "pch", "pcr", "sins", "tau")
    
    # option to remove signals
    if (!is.null(remove_basis)) {
        inds <- grep(remove_basis, mol_list_chars)
        if (length(inds) == 0) {
          print(mol_list_chars)
          stop("No signals (as listed above) matching remove_basis found.")
        }
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
      # round B0 to 5 s.f. for effective basis caching
      metab$ft <- signif(metab$ft, 5)
      # regen mol_list with updated B0
      mol_list <- get_mol_paras(mol_list_chars, ft = metab$ft)
      basis <- sim_basis(mol_list, acq_paras = metab,
                         pul_seq = seq_press_2d_shaped, TE1 = sim_paras$TE1,
                         TE2 = sim_paras$TE2, use_basis_cache = use_basis_cache,
                         verbose = verbose, pulse_file = pulse_file,
                         pulse_dur = 5e-3, pulse_file_format = "pta",
                         auto_scale = TRUE)
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
  }
  
  if (!is.null(external_basis)) {
    
    if (!inherits(external_basis, "basis_set")) {
      if (inherits(external_basis, "character")) {
        # TODO - could make this a bit more flexible, eg read a directory of
        # NIfTI MRS files if a dir name is passed, assume an MRS data file if
        # extension isn't .BASIS or .basis...
        external_basis <- read_basis(external_basis)
      } else {
        stop("Unrecognised exernal_basis object.")
      }
    }
    
    if (append_external_basis) {
      external_basis <- resample_basis(external_basis, metab)
      basis <- append_basis(basis, external_basis)
    } else {
      basis <- resample_basis(external_basis, metab)
    }
  }
  
  # check the fit_method is sane
  if (!is.null(fit_method)) {
    fit_method <- toupper(fit_method)
    allowed <- c("ABFIT-REG", "LCMODEL")
    if (!(fit_method %in% allowed)) {
      print(allowed)
      stop("Error, incorrect fit method, must be one of the above.")
    }
  }
  
  if (is.null(fit_method)) {
    fit_method <- "ABFIT"
    if (is.null(fit_opts)) fit_opts <- abfit_reg_opts()
  } else {
    if (fit_method == "ABFIT-REG") {
      fit_method <- "ABFIT"
      if (is.null(fit_opts)) fit_opts <- abfit_reg_opts()
    }
  }
  
  # ask LCModel not to simulate any additional signals by defaults
  if (fit_method == "LCMODEL" & is.null(fit_opts)) fit_opts <- c("NSIMUL=0")
  
  if (is.null(output_ratio)) {
    if (fit_method == "LCMODEL") {
      output_ratio <- "Cr.PCr"
    } else {
      output_ratio <- "tCr"
    }
  }
  
  # output_ratio of NA means we only want unscaled values
  if (anyNA(output_ratio)) output_ratio <- NULL
  
  # filter residual water
  if (!is.null(hsvd_width)) {
    if (verbose) cat("Applying HSVD filter.\n")
    metab <- hsvd_filt(metab, xlim = c(-hsvd_width, hsvd_width))
  }
  
  # truncate the FID if option is set
  if (!is.null(trunc_fid_pts)) {
    if (verbose) cat("Truncating FID.\n")
    metab <- crop_td_pts(metab, end = trunc_fid_pts)
    basis_mrs <- crop_td_pts(basis2mrs_data(basis), end = trunc_fid_pts)
    basis <- mrs_data2basis(basis_mrs, names = basis$names)
  }
  
  # set path to the LCModel binary
  if (!is.null(lcm_bin_path)) set_lcm_cmd(lcm_bin_path)
  
  # water referencing is performed internally by LCModel
  if (fit_method == "LCMODEL" & w_ref_available) {
    w_ref_fit <- w_ref 
  } else {
    w_ref_fit <- NULL
  }
  
  # fitting
  if (verbose) cat("Starting fitting.\n")
  fit_res <- fit_mrs(metab = metab, basis = basis, method = fit_method,
                     w_ref = w_ref_fit, opts = fit_opts)
  if (verbose) cat("Fitting complete.\n")
    
  phase_offset <- fit_res$res_tab$phase
  shift_offset <- fit_res$res_tab$shift
  
  if (w_ref_available) fit_res$res_tab$ws_eff <- ws_efficiency
  
  # keep unscaled results
  res_tab_unscaled <- fit_res$res_tab
  
  # assume 100% white matter unless told otherwise
  if (is.null(p_vols) & is.null(mri_seg)) {
    p_vols <- c(WM = 100, GM = 0, CSF = 0)
  }
  
  if (is.null(p_vols) & !is.null(mri_seg)) {
    # generate the svs voi in the segmented image space
    voi_seg <- get_svs_voi(metab, mri_seg)   
    
    # calculate partial volumes
    p_vols <- get_voi_seg(voi_seg, mri_seg)
  }
  
  # add an "Other" component to p_vols if missing (to keep things consistent)
  if (!is.null(p_vols)) if (!("Other" %in% names(p_vols))) p_vols["Other"] <- 0
 
  # output ratio results if requested 
  if (!is.null(output_ratio)) {
    for (output_ratio_element in output_ratio) {
      fit_res_rat <- scale_amp_ratio(fit_res, output_ratio_element,
                                     use_mean_value = TRUE)
      
      fit_res_rat$res_tab <- append_p_vols(fit_res_rat$res_tab, p_vols)
      
      res_tab_ratio <- fit_res_rat$res_tab
      file_out <- file.path(output_dir, paste0("fit_res_", output_ratio_element,
                                               "_ratio.csv"))
      
      utils::write.csv(res_tab_ratio, file_out)
    }
  } else {
    res_tab_ratio <- NULL  
  }
  
  # perform water reference amplitude scaling
  if (w_ref_available) {
    if (fit_method != "LCMODEL") {
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
      fit_res_molal <- scale_amp_molar2molal_pvc(fit_res, p_vols, TE, TR)
      res_tab_molal <- fit_res_molal$res_tab
      file_out <- file.path(output_dir, "fit_res_molal_conc.csv")
      utils::write.csv(res_tab_molal, file_out)
      res_tab_legacy <- NULL 
    }
  } else {
    res_tab_legacy <- NULL 
    res_tab_molal  <- NULL 
  }
  
  # add PVC info to the unscaled output
  res_tab_unscaled <- append_p_vols(res_tab_unscaled, p_vols)
  
  if (fit_method != "LCMODEL") {
    utils::write.csv(res_tab_unscaled, file.path(output_dir,
                                                 "fit_res_unscaled.csv"))
  }
  
  # prepare dynamic data for plotting
  if (Ndyns(metab_pre_dfp_corr) > 1) {
    # phase according to the fit results
    dyn_data_uncorr <- phase(metab_pre_dfp_corr, mean(fit_res$res_tab$phase))
    # correct chem. shift scale according to the fit results
    dyn_data_uncorr <- shift(dyn_data_uncorr, mean(fit_res$res_tab$shift),
                             units = "ppm")
    # add 2 Hz LB
    dyn_data_uncorr <- lb(dyn_data_uncorr, 2)
    
    if (!is.null(metab_post_dfp_corr)) {
      # phase according to the fit results
      dyn_data_corr <- phase(metab_post_dfp_corr, mean(fit_res$res_tab$phase))
      # correct chem. shift scale according to the fit results
      dyn_data_corr <- shift(dyn_data_corr, mean(fit_res$res_tab$shift),
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
  
  # generate a summary table
  if (is.null(summary_measures)) {
    summary_tab <- NULL
  } else {
    # extract the values from the fit results
    if (w_ref_available) {
      summary_tab <- parse_summary(summary_measures, fit_res_molal, " (mM)") 
    } else {
      summary_tab <- parse_summary(summary_measures, fit_res, " (a.u.)") 
    }
    summary_tab$values <- format(summary_tab$values, digits = 3)
  }
  
  # data needed to produce the output html report
  results <- list(fit_res = fit_res, argg = argg,
                  w_ref_available = w_ref_available,
                  w_ref = w_ref,
                  output_ratio = output_ratio,
                  res_tab_unscaled = res_tab_unscaled,
                  res_tab_ratio = res_tab_ratio,
                  res_tab_legacy = res_tab_legacy,
                  res_tab_molal = res_tab_molal,
                  dyn_data_uncorr = dyn_data_uncorr,
                  dyn_data_corr = dyn_data_corr,
                  summary_tab = summary_tab,
                  plot_ppm_xlim = plot_ppm_xlim,
                  mri = mri,
                  mri_seg = mri_seg,
                  p_vols = p_vols)
  
  if (Ndyns(metab) == 1) {
    rmd_file <- system.file("rmd", "svs_report.Rmd", package = "spant")
  } else {
    rmd_file <- system.file("rmd", "dyn_svs_report.Rmd", package = "spant")
  }
  
  rmd_out_f <- file.path(tools::file_path_as_absolute(output_dir), "report")
  
  if (verbose) cat("Generating html report.\n")
  rmarkdown::render(rmd_file, params = results, output_file = rmd_out_f,
                    quiet = !verbose)
  
  if (extra_output) {
    if (verbose) cat("Writing extra output files.\n")
    saveRDS(results, file = file.path(output_dir, "fit_res_data.rds"))
    utils::write.csv(results$fit_res$fits[[1]],
                     file = file.path(output_dir, "fit_plot_data.csv"),
                     row.names = FALSE)
  }
  
  if (verbose) cat("fit_svs finished.\n")
  
  return(fit_res)
}

check_sim_paras <- function(pul_seq, metab, TE1, TE2, TE3, TE, TM,
                            press_TE1_guess = 0.0126) {
  
  # try to guess the pulse sequence from the MRS data if not specified
  if (is.null(pul_seq)) {
    if (!is.null(metab$meta$PulseSequenceType)){
      pul_seq <- metab$meta$PulseSequenceType
    } else {
      warning(paste0("Could not determine the pulse sequence, so assuming\n",
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
      TE1 <- press_TE1_guess
      # warning("TE1 assumed to be 0.0126s. Provide the TE1 argument to stop this warning.")
    }
    if (is.null(TE2)) TE2 <- TE - TE1
    
    if (!isTRUE(all.equal(TE1 + TE2, TE))) {
      warning("TE, TE1 and TE2 do not match.")
    }
    
    return(list(pul_seq = pul_seq, TE1 = TE1, TE2 = TE2))
  } else if (pul_seq == "press_shaped") {
    if (is.null(TE1)) {
      TE1 <- press_TE1_guess
      # warning("TE1 assumed to be 0.0126s.")
    }
    if (is.null(TE2)) TE2 <- TE - TE1
    
    if (!isTRUE(all.equal(TE1 + TE2, TE))) {
      warning("TE, TE1 and TE2 do not match.")
    }
      
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
    if (!isTRUE(all.equal(TE1 + TE2 + TE3, TE))) {
      warning("TE, TE1, TE2 and TE3 do not match.")
    }
    return(list(pul_seq = pul_seq, TE1 = TE1, TE2 = TE2, TE3 = TE3))
  } else {
    stop(paste0("pul_seq not supported : ", pul_seq))
  }
}

parse_summary <- function(measures, fit_res, units) {
  res_names <- colnames(fit_res$res_tab)
  val_num   <- length(measures)
  values    <- rep(NA, val_num)
  
  for (n in 1:val_num) {
    measure <- gsub(" ", "", measures[n])

    if (length(grep("/", measure)) > 0) {
      # looks like a ratio
      numerator   <- strsplit(measure, "/")[[1]][1]
      denominator <- strsplit(measure, "/")[[1]][2]
      values[n]   <- as.numeric(fit_res$res_tab[numerator] / 
                                  fit_res$res_tab[denominator])
      measures[n] <- measure
    } else {
      # looks like a single value
      values[n]   <- as.numeric(fit_res$res_tab[measure])
      measures[n] <- paste0(measure, units)
    }
  }
  return(data.frame(measures = measures, values = values)) 
}

#' GUI interface for the standard SVS 1H brain analysis pipeline, this is a 
#' work in progress, and not ready for serious use.
fit_svs_gui <-function() {
  
  run_fit <- function() {
    pb <- tcltk::tkProgressBar("running analysis...")
    tcltk::setTkProgressBar(pb, 0)
    
    print(tcltk::tkget(wsup_path))
    
    fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT",
                         package = "spant")
    svs <- read_mrs(fname)
    basis <- sim_basis_1h_brain_press(svs)
    fit_result <- fit_mrs(svs, basis)
    tcltk::setTkProgressBar(pb, 100)
    close(pb)
    response <- tcltk::tk_messageBox("yesno", "Analysis completed, run another?")
    if (response == "no") tcltk::tkdestroy(tt) 
  }
  
  wsup_file_chooser <- function() {
    wsup_path_str <- tcltk::tk_choose.files()
    
    if (!identical(wsup_path_str, character())) {
      tcltk::tkconfigure(wsup_path, textvariable = tcltk::tclVar(wsup_path_str))
      
      # change output dir if not set
      if (identical(as.character(tcltk::tkget(output_path)), character())) {
        tcltk::tkconfigure(output_path,
                           textvariable = tcltk::tclVar(dirname(wsup_path_str)))
      }
    }
  }
  
  wsup_dir_chooser <- function() {
    wsup_path_str <- tcltk::tk_choose.dir()
    
    if (!identical(wsup_path_str, character())) {
      tcltk::tkconfigure(wsup_path, textvariable = tcltk::tclVar(wsup_path_str))
      
      # change output dir if not set
      if (identical(as.character(tcltk::tkget(output_path)), character())) {
        tcltk::tkconfigure(output_path,
                           textvariable = tcltk::tclVar(wsup_path_str))
      }
    }
  }
  
  wref_file_chooser <- function() {
    wref_path_str <- tcltk::tk_choose.files()
    
    if (!identical(wref_path_str, character())) {
      tcltk::tkconfigure(wref_path, textvariable = tcltk::tclVar(wref_path_str))
    }
  }
  
  wref_dir_chooser <- function() {
    wref_path_str <- tcltk::tk_choose.dir()
    
    if (!identical(wref_path_str, character())) {
      tcltk::tkconfigure(wref_path, textvariable = tcltk::tclVar(wref_path_str))
    }
  }
  
  output_dir_chooser <- function() {
    output_path_str <- tcltk::tk_choose.dir()
    
    if (!identical(output_path_str, character())) {
      tcltk::tkconfigure(output_path,
                         textvariable = tcltk::tclVar(output_path_str))
    }
  }
  
  clear <- function() {
    tcltk::tkconfigure(wsup_path,   textvariable = tcltk::tclVar(""))
    tcltk::tkconfigure(wref_path,   textvariable = tcltk::tclVar(""))
    tcltk::tkconfigure(output_path, textvariable = tcltk::tclVar(""))
  }
  
  tt <- tcltk::tktoplevel()
  tcltk::tktitle(tt) <- "spant GUI"
  
  heading <- tcltk::tklabel(tt, text = "spant SVS MRS analysis")
 
  wsup_lab          <- tcltk::tklabel(tt, text = "Water supressed data")
  wsup_file_button  <- tcltk::tkbutton(tt, text = "choose file",
                                       command = wsup_file_chooser)
  
  wsup_dir_button   <- tcltk::tkbutton(tt, text = "choose dir",
                                       command = wsup_dir_chooser)
  
  wsup_path         <- tcltk::tkentry(tt, width = 60)
  
  wref_lab          <- tcltk::tklabel(tt, text = "Water reference data")
  wref_file_button  <- tcltk::tkbutton(tt, text = "choose file",
                                       command = wref_file_chooser)
  
  wref_dir_button   <- tcltk::tkbutton(tt, text = "choose dir",
                                       command = wref_dir_chooser)
  
  wref_path         <- tcltk::tkentry(tt, width = 60)
  
  output_lab        <- tcltk::tklabel(tt, text = "Output directory")
  output_path       <- tcltk::tkentry(tt, width = 60)
  
  output_dir_button <- tcltk::tkbutton(tt, text = "choose dir",
                                       command = output_dir_chooser)
  
  run_button <- tcltk::tkbutton(tt, text = "run analysis", width = 20,
                                height = 2, command = run_fit)
  
  exit_button <- tcltk::tkbutton(tt, text = "exit",
                                 command = function() tcltk::tkdestroy(tt))
  
  clear_button <- tcltk::tkbutton(tt, text = "clear paths", command = clear)
  
  dummy_lab    <- tcltk::tklabel(tt, text = "")
  
  tcltk::tkgrid(heading, columnspan = 4, pady = 10)
  tcltk::tkgrid(wsup_lab, wsup_path, wsup_file_button, wsup_dir_button)
  tcltk::tkgrid(wref_lab, wref_path, wref_file_button, wref_dir_button)
  tcltk::tkgrid(output_lab, output_path, dummy_lab, output_dir_button)
  tcltk::tkgrid(run_button, columnspan = 3, pady = 10)
  tcltk::tkgrid(clear_button, exit_button, columnspan = 3, pady = 20)
  
  tcltk::tkgrid.configure(wsup_lab, wref_lab, output_lab, sticky = "e")
  tcltk::tkgrid.configure(clear_button, exit_button, sticky = "e")
  
}

