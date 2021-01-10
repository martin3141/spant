#' @export
svs_1h_preproc <- function(metab, ref = NULL, decimate = FALSE,
                           rats_corr = TRUE, ecc = FALSE, comb_dyns = TRUE,
                           hsvd_filt = FALSE) {
  
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
svs_1h_batch_fit <- function(metab, ref = NULL, basis,
                             preproc_fn = svs_1h_preproc, ...) {
  
  if (!is.null(ref_fname) & (length(metab_fname) != length(ref_fname))) {
    stop("Number of metabolite and reference files do not match.")
  }
  
  
}