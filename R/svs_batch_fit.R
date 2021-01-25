#' @export
svs_1h_analysis <- function(metab, basis, w_ref = NULL, mri_seg = NULL,
                            decimate = FALSE, rats_corr = TRUE, ecc = FALSE,
                            comb_dyns = TRUE, hsvd_filt = FALSE,
                            scale_amps = TRUE) {
  
  # combine coils if needed
  if (Ncoils(metab) > 1) {
    coil_comb_res <- comb_coils(metab, w_ref)
    if (is.null(w_ref)) {
      metab <- coil_comb_res
    } else {
      metab <- coil_comb_res$metab
      w_ref   <- coil_comb_res$ref
    }
  }
  
  # decimate x2 if requested
  if (decimate) {
    metab <- decimate_mrs_fd(metab)
    if (!is.null(w_ref)) w_ref <- decimate_mrs_fd(w_ref)
  }
  
  # rats
  if (rats_corr) metab <- rats(metab)$corrected
  
  # TODO plot of shifts?
  
  # eddy current correction
  if (ecc & (!is.null(w_ref))) metab <- ecc(metab, w_ref)
  
  # combine dynamic scans
  if (comb_dyns) metab <- mean_dyns(metab)
  
  # HSVD residual water removal
  if (hsvd_filt) metab <- hsvd_filt(metab)
  
  # TODO fitting
  fit_res <- fit_mrs(metab, basis = basis)
  
  # if (scale_amps) {
  #   if (is.null(w_ref)) {
  #     
  #   } else {
  #     if (is.null())
  #     
  #   }
  # }
  
  # if water file but no seg file then scale using
  
  # otherwise scale to tCr
  
  # TODO optional pvc
  
  # TODO optional amp scaling
  
  return(fit_res)
}

# TODO
# what if metab fname results in metab + ref file, eg GE data?
# options to add, fit_opts, fit_method, format (GE, Siemens etc)
# allow lists of mrs_data objects, or file paths as input
# add extra option for id vars

#' @export
svs_1h_batch_analysis <- function(metab_paths, w_ref_paths = NULL,
                                  mri_paths = NULL, mri_seg_paths = NULL,
                                  output_dirs = NULL, ...) {
  
  # check input it sensible
  metab_n <- length(metab_paths)
  
  if (!is.null(w_ref_paths) & length(w_ref_paths) != metab_n) {
    stop("Incorrect number of ref items.")
  }
  
  if (!is.null(mri_paths) & length(mri_paths) != metab_n) {
    stop("Incorrect number of mri items.")
  }
  
  if (!is.null(mri_seg_paths) & length(mri_seg_paths) != metab_n) {
    stop("Incorrect number of mri_seg items.")
  }
  
  if (!is.null(output_dirs) & length(output_dirs) != metab_n) {
    stop("Incorrect number of output_dir items.")
  }
  
  # check file paths exist TODO write a function to return bad paths
  # metab_paths[which(!file.exists(metab_paths))]
  
  # convert file paths to a list of mrs_data objects
  metab_list <- lapply(metab_paths, read_mrs)
  
  if (is.null(w_ref_paths)) {
    w_ref_list <- vector("list", metab_n)
  } else {
    w_ref_list <- lapply(w_ref_paths, read_mrs)
  }
  
  fit_list <- mapply(svs_1h_analysis, metab = metab_list, ref = w_ref_list,
                     mri_seg = mri_seg_list, MoreArgs = ..., SIMPLIFY = FALSE)
  
  return(fit_list)
}
