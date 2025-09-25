#' Edited SVS 1H brain analysis pipeline.
#' 
#' Note this function is still under development and liable to changes.
#' 
#' @param input path or mrs_data object containing MRS data.
#' @param w_ref path or mrs_data object containing MRS water reference data.
#' @param output_dir directory path to output fitting results.
#' @param mri filepath or nifti object containing anatomical MRI data.
#' @param mri_seg filepath or nifti object containing segmented MRI data.
#' @param external_basis precompiled basis set object to use for analysis.
#' @param p_vols a numeric vector of partial volumes expressed as percentages.
#' Defaults to 100% white matter. A voxel containing 100% gray matter tissue
#' would use : p_vols = c(WM = 0, GM = 100, CSF = 0).
#' @param format Override automatic data format detection. See format argument
#' in [read_mrs()] for permitted values.
#' @param editing_type can be one of : "gaba_1.9" or "gsh_4.54". Defaults to 
#' "gaba_1.9".
#' @param editing_scheme describes the dynamic data ordering. Can be one of:
#' 'on-off-blocks', 'on-off-interleaved', 'off-on-blocks' or
#' 'off-on-interleaved'.
#' @param invert_edit_on set to TRUE to invert the edit-on sub-spectra.
#' @param invert_edit_off set to TRUE to invert the edit-off sub-spectra.
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
#' @param append_basis_ed_off names of extra signals to add to the default
#' basis. Eg append_basis_ed_off = c("peth", "cit"). Cannot be used with
#' precompiled basis sets.
#' @param remove_basis_ed_off grep expression to match names of signals to
#' remove from the basis. For example: use "*" to remove all signals, "^mm|^lip"
#' to remove all macromolecular and lipid signals, "^lac" to remove lactate.
#' This operation is performed before signals are added with
#' append_basis_ed_off. Cannot be used with precompiled basis sets.
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
#' @param fit_opts_edited options to pass to the fitting method for the 
#' edited spectrum.
#' @param fit_opts_ed_off options to pass to the fitting method for the
#' edit-off spectrum.
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
fit_svs_edited <- function(input, w_ref = NULL, output_dir = NULL, mri = NULL,
                           mri_seg = NULL, external_basis = NULL, p_vols = NULL,
                           format = NULL, editing_type = "gaba_1.9",
                           editing_scheme = NULL,
                           invert_edit_on = NULL, invert_edit_off = NULL,
                           pul_seq = NULL, TE = NULL, TR = NULL,
                           TE1 = NULL, TE2 = NULL, TE3 = NULL, TM = NULL,
                           append_basis_ed_off = NULL,
                           remove_basis_ed_off = NULL,
                           pre_align = TRUE, dfp_corr = TRUE,
                           output_ratio = NULL, ecc = FALSE,
                           hsvd_width = NULL, decimate = FALSE,
                           trunc_fid_pts = NULL, fit_opts_edited = NULL,
                           fit_opts_ed_off = NULL,
                           fit_subset = NULL, legacy_ws = FALSE, w_att = 0.7,
                           w_conc = 35880, use_basis_cache = "auto",
                           summary_measures = NULL, dyn_av_block_size = NULL,
                           dyn_av_scheme = NULL, dyn_av_scheme_file = NULL,
                           plot_ppm_xlim = NULL, extra_output = FALSE,
                           verbose = FALSE) {
  
  warning("fit_svs_exited is under active development and liable to significant changes.")
  
  argg  <- c(as.list(environment()))
  
  metab <- input
  
  if (!is.null(dyn_av_scheme) & !is.null(dyn_av_scheme_file)) {
    print(dyn_av_scheme)
    print(dyn_av_scheme_file)
    stop("dyn_av_scheme and dyn_av_scheme_file options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(external_basis) & !is.null(append_basis_ed_off)) {
    stop("external_basis and append_basis_ed_off options cannot both be set. Use one or the other.")
  }
  
  if (!is.null(external_basis) & !is.null(remove_basis_ed_off)) {
    stop("external_basis and remove_basis_ed_off options cannot both be set. Use one or the other.")
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
  
  # determine edit-on / edit-off scans, likely to be vendor and sequence
  # dependant
  if (is.null(editing_scheme))  editing_scheme  <- "off-on-interleaved"
  if (is.null(invert_edit_on))  invert_edit_on  <- FALSE
  if (is.null(invert_edit_off)) invert_edit_off <- FALSE
  
  if (editing_scheme == "off-on-interleaved") {
    ed_off <- get_odd_dyns(metab)
    ed_on  <- get_even_dyns(metab)
  } else if (editing_scheme == "on-off-interleaved") {
    ed_on  <- get_odd_dyns(metab)
    ed_off <- get_even_dyns(metab)
  } else if (editing_scheme == "off-on-blocks") {
    ed_off <- get_fh_dyns(metab)
    ed_on  <- get_sh_dyns(metab)
  } else if (editing_scheme == "on-off-blocks") {
    ed_on  <- get_fh_dyns(metab)
    ed_off <- get_sh_dyns(metab)
  } else {
    print(editing_scheme)
    stop("Incorrect editing_scheme string.")
  }
  
  editing_types <- c("gaba_1.9", "gsh_4.54")
  if (!(editing_type %in% editing_types)) {
    print(editing_types)
    stop("editing_type not recognised, should be one of the above.")
  }
  
  if (invert_edit_off) ed_off <- phase(ed_off, 180)
  if (invert_edit_on)  ed_on  <- phase(ed_on,  180)
  
  # reorganise dynamics into two blocks
  metab <- append_dyns(ed_on, ed_off)
  
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
    if (w_ref_available) w_ref <- decimate_mrs_fd(w_ref)
  }
  
  ed_on  <- get_fh_dyns(metab)
  ed_off <- get_sh_dyns(metab)
  
  # extract a subset of dynamic scans if specified
  if (!is.null(fit_subset)) ed_off <- get_dyns(ed_off, fit_subset) 
  if (!is.null(fit_subset)) ed_on  <- get_dyns(ed_on,  fit_subset) 
  
  ed_off_pre_dfp_corr <- ed_off
  ed_on_pre_dfp_corr  <- ed_on
  
  # pre-alignment
  if (pre_align) {
    ed_off <- align(ed_off, c(2.01, 3.03, 3.22), max_shift = 40)
    if (Ndyns(ed_off) > 1) ed_off_post_dfp_corr <- ed_off
    
    if (editing_type == "gaba_1.9") {
      ed_on <- align(ed_on, c(3.03, 3.22), max_shift = 40)
    } else {
      ed_on <- align(ed_on, c(2.01, 3.03, 3.22), max_shift = 40)
    }
    
    if (Ndyns(ed_on) > 1) ed_on_post_dfp_corr <- ed_on
  }
  
  # rats correction
  if (dfp_corr & (Ndyns(ed_off) > 1)) {
    ed_off <- rats(ed_off, zero_freq_shift_t0 = TRUE, xlim = c(4, 1.8))
    ed_off_post_dfp_corr <- ed_off
    ed_on  <- rats(ed_on, zero_freq_shift_t0 = TRUE, xlim = c(4, 1.8))
    ed_on_post_dfp_corr <- ed_on
  }
  
  if (!exists("ed_off_post_dfp_corr")) ed_off_post_dfp_corr <- NULL
  if (!exists("ed_on_post_dfp_corr"))  ed_on_post_dfp_corr <- NULL
  
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
    ed_off <- mean_dyn_blocks(ed_off, dyn_av_block_size)
    ed_on  <- mean_dyn_blocks(ed_on,  dyn_av_block_size)
  } else if (!is.null(dyn_av_scheme)) {
    if (length(dyn_av_scheme) != Ndyns(ed_off)) {
      stop(paste0("dyn_av_scheme is the wrong length. Currently : ",
                  length(dyn_av_scheme),", should be : ", Ndyns(ed_off)))
    }
    dyn_av_scheme <- as.integer(dyn_av_scheme)
    max_dyn       <- max(dyn_av_scheme)
    ed_off_list   <- vector("list", length = max_dyn)
    ed_on_list    <- vector("list", length = max_dyn)
    for (n in 1:max_dyn) {
      subset <- which(dyn_av_scheme == n) 
      ed_off_list[[n]] <- mean_dyns(get_dyns(ed_off, subset))
      ed_on_list[[n]]  <- mean_dyns(get_dyns(ed_on,  subset))
    }
    ed_off <- append_dyns(ed_off_list)
    ed_on  <- append_dyns(ed_on_list)
  } else {
    ed_off <- mean_dyns(ed_off)
    ed_on  <- mean_dyns(ed_on)
  }
 
  # take the mean of the water reference data
  # if (w_ref_available) w_ref <- mean_dyns(w_ref)
  
  # align dynamic water reference scans to the first dynamic using rats and take
  # the mean
  if (w_ref_available) {
    if (Ndyns(w_ref) > 1) {
      w_ref_ref <- get_dyns(w_ref, 1)
      w_ref <- rats(w_ref, xlim = c(5.3, 4), ref = w_ref_ref)
      w_ref <- mean_dyns(w_ref)
    }
  }
  
  # extract the first dynamic of the water reference data
  # better for Dinesh sLASER where the water peak can shift before and after
  # the metabolite data collection
  # if (w_ref_available) w_ref <- get_dyns(w_ref, 1)
  
  # calculate the water suppression efficiency
  # the ratio of the residual water peak height
  # relative to the height of the unsuppressed water signal
  if (w_ref_available) {
    w_ref_peak   <- peak_info(w_ref, xlim = c(7, 3), mode = "mod")
    w_ref_height <- w_ref_peak$height[1]
    w_ref_freq   <- w_ref_peak$freq_ppm[1]
    x_range      <- c(w_ref_freq - 0.2, w_ref_freq + 0.2)
    ed_off_water_height <- spec_op(ed_off, xlim = x_range, operator = "max",
                                   mode = "mod")[1]
    ws_efficiency <- ed_off_water_height / w_ref_height * 100
  }
  
  # eddy current correction
  if (ecc & w_ref_available) {
    ed_off <- ecc(ed_off, w_ref)
    ed_on  <- ecc(ed_on,  w_ref)
  }
  
  # simulate a basis if needed
  if (is.null(external_basis)) {
    
    if (is.null(TE)) stop("Could not determine the sequence echo time. Please provide the TE argument.")
    
    # list of "standard" signals to include in the basis set
    mol_list_chars <- c("m_cr_ch2", "ala", "asp", "cr", "gaba", "glc", "gln",
                        "gsh", "glu", "gpc", "ins", "lac", "lip09", "lip13a",
                        "lip13b", "lip20", "mm09", "mm12", "mm14", "mm17",
                        "mm20", "naa", "naag", "pch", "pcr", "sins", "tau")
    
    # option to remove signals
    if (!is.null(remove_basis_ed_off)) {
        inds <- grep(remove_basis_ed_off, mol_list_chars)
        if (length(inds) == 0) {
          print(mol_list_chars)
          stop("No signals (as listed above) matching remove_basis_ed_off found.")
        }
        mol_list_chars <- mol_list_chars[-inds]
    }
    
    # option to append signals
    if (!is.null(append_basis_ed_off)) mol_list_chars <- c(mol_list_chars,
                                                           append_basis_ed_off)
    
    # probably set remove_basis_ed_off to * and forgot to use append_basis_ed_off
    if (is.null(mol_list_chars)) stop("No basis signals named for simulation.")
    
    # get the parameters
    mol_list <- get_mol_paras(mol_list_chars, ft = ed_off$ft)
    
    # check parameters are consistent and infer any missing values
    sim_paras <- check_sim_paras(pul_seq, ed_off, TE1, TE2, TE3, TE, TM)
    
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
      basis <- sim_basis(mol_list, acq_paras = ed_off,
                         pul_seq = seq_press_ideal, TE1 = sim_paras$TE1,
                         TE2 = sim_paras$TE2, use_basis_cache = use_basis_cache,
                         verbose = verbose)
    } else if (sim_paras$pul_seq == "press_shaped") {
      if (verbose) cat("Simulating shaped PRESS sequence.\n")
      pulse_file <- system.file("extdata", "press_refocus.pta",
                                package = "spant")
      # round B0 to 5 s.f. for effective basis caching
      ed_off$ft <- signif(ed_off$ft, 5)
      # regen mol_list with updated B0
      mol_list <- get_mol_paras(mol_list_chars, ft = ed_off$ft)
      basis <- sim_basis(mol_list, acq_paras = ed_off,
                         pul_seq = seq_press_2d_shaped, TE1 = sim_paras$TE1,
                         TE2 = sim_paras$TE2, use_basis_cache = use_basis_cache,
                         verbose = verbose, pulse_file = pulse_file,
                         pulse_dur = 5e-3, pulse_file_format = "pta",
                         auto_scale = TRUE)
    } else if (sim_paras$pul_seq == "steam") {
      if (verbose) cat("Simulating STEAM sequence.\n")
      basis <- sim_basis(mol_list, acq_paras = ed_off,
                         pul_seq = seq_steam_ideal_cof, TE = sim_paras$TE,
                         TM = sim_paras$TM, use_basis_cache = use_basis_cache,
                         verbose = verbose)
    } else if (sim_paras$pul_seq == "slaser") {
      if (verbose) cat("Simulating sLASER sequence.\n")
      basis <- sim_basis(mol_list, acq_paras = ed_off,
                         pul_seq = seq_slaser_ideal, TE1 = sim_paras$TE1,
                         TE2 = sim_paras$TE2, TE3 = sim_paras$TE3,
                         use_basis_cache = use_basis_cache, verbose = verbose)
    }
    if (verbose) print(basis)
  } else {
    basis <- external_basis
  }
  
  # if (is.null(fit_opts_edited)) fit_opts_edited <- abfit_reg_opts()
  
  if (is.null(fit_opts_edited)) {
    if (editing_type == "gsh_4.54") {
      fit_opts_edited <- abfit_reg_opts(auto_bl_flex = FALSE, bl_ed_pppm = 3,
                                        ppm_left = 3.5, ppm_right = 1.8,
                                        pre_align = FALSE)
    } else {
      fit_opts_edited <- abfit_reg_opts(auto_bl_flex = FALSE, bl_ed_pppm = 3,
                                        pre_align = FALSE)
    }
  }
  
  if (is.null(fit_opts_ed_off)) fit_opts_ed_off <- abfit_reg_opts()
  
  if (is.null(output_ratio)) output_ratio <- "tCr"
  
  # output_ratio of NA means we only want unscaled values
  if (anyNA(output_ratio)) output_ratio <- NULL
  
  if (editing_type == "gaba_1.9") {
    # align ed_on and ed_off based on the residual water signal
    ed_on  <- rats(ed_on, mean_dyns(ed_off), xlim = c(4.8, 4.5))
  } else if (editing_type == "gsh_4.54") {
    # align ed_on and ed_off based on the residual tNAA
    ed_on  <- rats(ed_on, mean_dyns(ed_off), xlim = c(1.9, 2.1))
  }
  
  # take the mean rather than just straight subtraction
  edited <- (ed_on - ed_off) / 2
  
  # filter residual water
  if (!is.null(hsvd_width)) {
    if (verbose) cat("Applying HSVD filter.\n")
    ed_off <- hsvd_filt(ed_off, xlim = c(-hsvd_width, hsvd_width))
    ed_on  <- hsvd_filt(ed_on,  xlim = c(-hsvd_width, hsvd_width))
    edited <- hsvd_filt(edited, xlim = c(-hsvd_width, hsvd_width))
  }
  
  # truncate the FID if option is set
  if (!is.null(trunc_fid_pts)) {
    if (verbose) cat("Truncating FID.\n")
    ed_off <- crop_td_pts(ed_off, end = trunc_fid_pts)
    ed_on  <- crop_td_pts(ed_on,  end = trunc_fid_pts)
    edited <- crop_td_pts(edited, end = trunc_fid_pts)
    if (w_ref_available) w_ref <- crop_td_pts(w_ref, end = trunc_fid_pts)
    basis_mrs <- crop_td_pts(basis2mrs_data(basis), end = trunc_fid_pts)
    basis <- mrs_data2basis(basis_mrs, names = basis$names)
  }
  
  # edit-off fitting
  if (verbose) cat("Starting edit-off fitting.\n")
  fit_res <- fit_mrs(metab = ed_off, basis = basis, opts = fit_opts_ed_off)
  if (verbose) cat("Edit-off fitting complete.\n")
  
  if (editing_type == "gaba_1.9") {
    # edited fitting
    # NAA is 1.5 (rather than 3) because it is zero in the edited data due to
    # the GABA editing pulse
    mol_list <- list(get_uncoupled_mol("MM09",   0.92, "1H",   1,   12, 1),
                     get_uncoupled_mol("NAA",    2.01, "1H",  -1.5,  3, 0),
                     get_uncoupled_mol("Glx_A",  2.31, "1H",   1,    3, 0),
                     get_uncoupled_mol("Glx_B",  2.40, "1H",   1,    3, 0),
                     # 2 Gaus model
                     get_uncoupled_mol("GABA_A", 2.95, "1H",   1.5,   12, 1),
                     get_uncoupled_mol("GABA_B", 3.04, "1H",   1.5,   12, 1),
                     # 1 Gaus model
                     # get_uncoupled_mol("GABA",   3.00, "1H",   3,    22, 1),
                     get_uncoupled_mol("Glx_C",  3.72, "1H",   1,   2.5, 0),
                     get_uncoupled_mol("Glx_D",  3.8,  "1H",   1,   2.5, 0))
  } else if (editing_type == "gsh_4.54") {
    ang <- 1i / 180 * pi
    damp_adj <- 0.7
    # mol_list <- list(get_uncoupled_mol("GSH",  2.960, "1H", 1, 9, 0),
    mol_list <- list(get_uncoupled_mol("GSH",  2.960, "1H", 1, 14, 1),
                     get_uncoupled_mol("P237", 2.372, "1H", exp(23 * ang),
                                       4.0 - damp_adj, 0),
                     get_uncoupled_mol("P246", 2.455, "1H", exp(-5 * ang),
                                       6.7 - damp_adj, 0),
                     get_uncoupled_mol("P257", 2.570, "1H", exp(-68 * ang),
                                       4.8 - damp_adj, 0),
                     get_uncoupled_mol("P265", 2.649, "1H", exp(76 * ang),
                                       5.5 - damp_adj, 0),
                     # get_uncoupled_mol("P275", 2.746, "1H", exp(-12 * ang),
                     #                   5.8 - damp_adj, 0))
                     get_uncoupled_mol("P275", 2.746, "1H", exp(-5 * ang),
                                       5.8 - damp_adj, 0))
  }
  
  basis_ed <- sim_basis(mol_list, acq_paras = edited)
  
  if (verbose) cat("Starting edited fitting.\n")
  fit_res_ed <- fit_mrs(metab = edited, basis = basis_ed,
                        opts = fit_opts_edited)
  if (verbose) cat("Edited fitting complete.\n")
    
  phase_offset <- fit_res$res_tab$phase
  shift_offset <- fit_res$res_tab$shift
  
  if (w_ref_available) fit_res$res_tab$ws_eff <- ws_efficiency
  
  # keep unscaled results
  res_tab_unscaled    <- fit_res$res_tab
  res_tab_ed_unscaled <- fit_res_ed$res_tab
  
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
      
      value <- mean(as.numeric(fit_res$res_tab[[output_ratio_element]]))
      fit_res_rat <- scale_amp_ratio_value(fit_res, value)
      
      fit_res_rat$res_tab <- append_p_vols(fit_res_rat$res_tab, p_vols)
      
      res_tab_ratio <- fit_res_rat$res_tab
      file_out <- file.path(output_dir, paste0("fit_res_edit_off_",
                                        output_ratio_element, "_ratio.csv"))
      utils::write.csv(res_tab_ratio, file_out, row.names = FALSE)
      
      # edited ratio
      fit_res_ed_rat <- scale_amp_ratio_value(fit_res_ed, value)
      
      fit_res_ed_rat$res_tab <- append_p_vols(fit_res_ed_rat$res_tab, p_vols)
      # fit_res_ed_rat$res_tab <- append_mpress_gaba(fit_res_ed_rat$res_tab)
      
      res_tab_ed_ratio <- fit_res_ed_rat$res_tab
      file_out <- file.path(output_dir, paste0("fit_res_edited_",
                                        output_ratio_element, "_ratio.csv"))
      utils::write.csv(res_tab_ed_ratio, file_out, row.names = FALSE)
    }
  } else {
    res_tab_ratio    <- NULL
    res_tab_ed_ratio <- NULL
  }
  
  # perform water reference amplitude scaling
  if (w_ref_available) {
    
    # edit off
    fit_res_molal <- scale_amp_molal_pvc(fit_res, w_ref, p_vols, TE, TR)
    res_tab_molal <- fit_res_molal$res_tab
    file_out <- file.path(output_dir, "fit_res_edit_off_molal_conc.csv")
    utils::write.csv(res_tab_molal, file_out, row.names = FALSE)
    
    # edited
    fit_res_ed_molal <- scale_amp_molal_pvc(fit_res_ed, w_ref, p_vols, TE, TR)
    # fit_res_ed_molal$res_tab <- append_mpress_gaba(fit_res_ed_molal$res_tab)
    res_tab_ed_molal <- fit_res_ed_molal$res_tab
    file_out <- file.path(output_dir, "fit_res_edited_molal_conc.csv")
    utils::write.csv(res_tab_ed_molal, file_out, row.names = FALSE)
    
    if (legacy_ws) {
      # edit off
      fit_res_legacy <- scale_amp_legacy(fit_res, w_ref, w_att, w_conc)
      res_tab_legacy <- fit_res_legacy$res_tab
      file_out <- file.path(output_dir, "fit_res_edit_off_legacy_conc.csv")
      utils::write.csv(res_tab_legacy, file_out, row.names = FALSE)
      
      # edited
      fit_res_ed_legacy <- scale_amp_legacy(fit_res_ed, w_ref, w_att, w_conc)
      # fit_res_ed_legacy$res_tab <- append_mpress_gaba(fit_res_ed_legacy$res_tab)
      res_tab_ed_legacy <- fit_res_ed_legacy$res_tab
      file_out <- file.path(output_dir, "fit_res_edited_legacy_conc.csv")
      utils::write.csv(res_tab_ed_legacy, file_out, row.names = FALSE)
    } else {
      res_tab_legacy <- NULL
      res_tab_ed_legacy <- NULL
    }
  } else {
    res_tab_legacy <- NULL 
    res_tab_molal  <- NULL 
    res_tab_ed_legacy <- NULL 
    res_tab_ed_molal  <- NULL 
  }
  
  # add PVC info to the unscaled output and write to csv
  res_tab_unscaled <- append_p_vols(res_tab_unscaled, p_vols)
  utils::write.csv(res_tab_unscaled, file.path(output_dir,
                                               "fit_res_edit_off_unscaled.csv"),
                   row.names = FALSE)
  res_tab_ed_unscaled <- append_p_vols(res_tab_ed_unscaled, p_vols)
  # res_tab_ed_unscaled <- append_mpress_gaba(res_tab_ed_unscaled)
  utils::write.csv(res_tab_ed_unscaled, file.path(output_dir,
                                               "fit_res_edited_unscaled.csv"),
                   row.names = FALSE)
  
  # prepare dynamic data for plotting
  if (Ndyns(ed_off_pre_dfp_corr) > 1) {
    # phase according to the fit results
    dyn_data_uncorr_ed_off <- phase(ed_off_pre_dfp_corr,
                                    mean(fit_res$res_tab$phase))
    # correct chem. shift scale according to the fit results
    dyn_data_uncorr_ed_off <- shift(dyn_data_uncorr_ed_off,
                                    mean(fit_res$res_tab$shift), units = "ppm")
    # add 2 Hz LB
    dyn_data_uncorr_ed_off <- lb(dyn_data_uncorr_ed_off, 2)
    
    if (!is.null(ed_off_post_dfp_corr)) {
      # phase according to the fit results
      dyn_data_corr_ed_off <- phase(ed_off_post_dfp_corr,
                                    mean(fit_res$res_tab$phase))
      # correct chem. shift scale according to the fit results
      dyn_data_corr_ed_off <- shift(dyn_data_corr_ed_off,
                                    mean(fit_res$res_tab$shift),
                                    units = "ppm")
      # add 2 Hz LB
      dyn_data_corr_ed_off <- lb(dyn_data_corr_ed_off, 2)
    } else {
      dyn_data_corr_ed_off <- NULL
    }
  } else {
    dyn_data_uncorr_ed_off <- NULL
    dyn_data_corr_ed_off   <- NULL
  }
  if (Ndyns(ed_on_pre_dfp_corr) > 1) {
    # phase according to the fit results
    dyn_data_uncorr_ed_on <- phase(ed_on_pre_dfp_corr,
                                   mean(fit_res$res_tab$phase))
    # correct chem. shift scale according to the fit results
    dyn_data_uncorr_ed_on <- shift(dyn_data_uncorr_ed_on,
                                   mean(fit_res$res_tab$shift), units = "ppm")
    # add 2 Hz LB
    dyn_data_uncorr_ed_on <- lb(dyn_data_uncorr_ed_on, 2)
    
    if (!is.null(ed_on_post_dfp_corr)) {
      # phase according to the fit results
      dyn_data_corr_ed_on <- phase(ed_on_post_dfp_corr,
                                    mean(fit_res$res_tab$phase))
      # correct chem. shift scale according to the fit results
      dyn_data_corr_ed_on <- shift(dyn_data_corr_ed_on,
                                   mean(fit_res$res_tab$shift),
                                  units = "ppm")
      # add 2 Hz LB
      dyn_data_corr_ed_on <- lb(dyn_data_corr_ed_on, 2)
    } else {
      dyn_data_corr_ed_on <- NULL
    }
  } else {
    dyn_data_uncorr_ed_on <- NULL
    dyn_data_corr_ed_on   <- NULL
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
  results <- list(fit_res = fit_res,
                  fit_res_ed = fit_res_ed,
                  argg = argg,
                  w_ref_available = w_ref_available,
                  w_ref = w_ref,
                  ed_on = ed_on,
                  output_ratio = output_ratio,
                  res_tab_unscaled = res_tab_unscaled,
                  res_tab_ratio = res_tab_ratio,
                  res_tab_legacy = res_tab_legacy,
                  res_tab_molal = res_tab_molal,
                  res_tab_ed_unscaled = res_tab_ed_unscaled,
                  res_tab_ed_ratio = res_tab_ed_ratio,
                  res_tab_ed_legacy = res_tab_ed_legacy,
                  res_tab_ed_molal = res_tab_ed_molal,
                  dyn_data_uncorr_ed_off = dyn_data_uncorr_ed_off,
                  dyn_data_corr_ed_off = dyn_data_corr_ed_off,
                  dyn_data_uncorr_ed_on = dyn_data_uncorr_ed_on,
                  dyn_data_corr_ed_on = dyn_data_corr_ed_on,
                  summary_tab = summary_tab,
                  plot_ppm_xlim = plot_ppm_xlim,
                  mri = mri,
                  mri_seg = mri_seg,
                  p_vols = p_vols,
                  editing_type = editing_type)
  
  rmd_file <- system.file("rmd", "svs_edited_report.Rmd", package = "spant")
  
  rmd_out_f <- file.path(tools::file_path_as_absolute(output_dir), "report")
  
  if (verbose) cat("Generating html report.\n")
  rmarkdown::render(rmd_file, params = results, output_file = rmd_out_f,
                    quiet = !verbose)
  
  saveRDS(results, file = file.path(output_dir,
                                    "spant_fit_svs_edited_data.rds"))
  
  if (extra_output) {
    if (verbose) cat("Writing extra output files.\n")
    
    warning("extra_output doesn't do anything at the moment...")
    
    # utils::write.csv(results$fit_res$fits[[1]],
    #                  file = file.path(output_dir, "fit_plot_data_edit_off.csv"),
    #                  row.names = FALSE)
    # utils::write.csv(results$fit_res_ed$fits[[1]],
    #                  file = file.path(output_dir, "fit_plot_data_edited.csv"),
    #                  row.names = FALSE)
  }
  
  if (verbose) cat("fit_svs_edited finished.\n")
  
  return(list(fit_res_ed, fit_res))
}


append_mpress_gaba <- function(res_tab) {
  res_tab["GABA"] <- res_tab["GABA_A"] + res_tab["GABA_B"]
  res_tab["Glx"]  <- res_tab["Glx_C"]  + res_tab["Glx_D"]
  return(res_tab)
}