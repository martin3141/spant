#' @export
svs_1h_fit <- function(metab, ref = NULL, decimate = FALSE, rats_corr = TRUE,
                       ecc = FALSE, comb_dyns = TRUE, hsvd_filt = FALSE) {
  
  # combine coils if needed
  if (Ncoils(metab) > 1) {
    coil_comb_res <- comb_coils(metab, ref)
    if (is.null(ref)) {
      metab <- coil_comb_res
    } else {
      metab <- coil_comb_res$metab
      ref   <- coil_comb_res$ref
    }
  }
  
  # decimate x2 if requested
  if (decimate) {
    metab <- decimate_mrs_fd(metab)
    if (!is.null(ref)) ref <- decimate_mrs_fd(ref)
  }
  
  # rats
  if (rats_corr) metab <- rats(metab)$corrected
  
  # TODO plot of shifts?
  
  # eddy current correction
  if (ecc & (!is.null(ref))) metab <- ecc(metab, ref)
  
  # combine dynamic scans
  if (comb_dyns) metab <- mean_dyns(metab)
  
  # HSVD residual water removal
  if (hsvd_filt) metab <- hsvd_filt(metab)
  
  if (is.null(ref)) {
    return(metab)
  } else {
    return(list(metab = metab, ref = ref))
  }
}

# TODO
# what if metab fname results in metab + ref file, eg GE data?
# options to add, fit_opts, fit_method, format (GE, Siemens etc)
# allow lists of mrs_data objects, or file paths as input
# add extra option for id vars

#' @export
svs_1h_batch_fit <- function(metab_paths, ref_paths = NULL, mri_paths = NULL,
                             mri_seg_paths = NULL, output_dirs = NULL, ...) {
  
  # check input it sensible
  metab_n <- length(metab_paths)
  
  if (!is.null(ref_paths) & length(ref_paths) != metab_n) {
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
  
  if (is.null(ref_paths)) {
    ref_list <- vector("list", metab_n)
  } else {
    ref_list <- lapply(ref_paths, read_mrs)
  }
  
  fit_list <- mapply(svs_1h_fit, metab = metab_list, ref = ref_list,
                     MoreArgs = ..., SIMPLIFY = FALSE)
  
  return(fit_list)
}
