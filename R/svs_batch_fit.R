# TODO add option to explicity specify the type of metabolite scaling to perform
# include "auto" version which guesses based on the input

#' @export
svs_1h_analysis <- function(metab, basis = NULL, w_ref = NULL, mri_seg = NULL,
                            mri = NULL, output_dir = NULL, extra = NULL,
                            decimate = NULL, rats_corr = TRUE, ecc = FALSE,
                            comb_dyns = TRUE, hsvd_filt = FALSE,
                            scale_amps = TRUE, te = NULL, tr = NULL,
                            preproc_only = FALSE) {
  
  if (!preproc_only & is.null(basis)) stop("basis argument not set")
  
  # read the data file if not already an mrs_data object
  if (class(metab)[[1]] != "mrs_data") metab <- read_mrs(metab)
  
  # remove fs_path types as they cause problems with comb_fit_list_fit_tables
  extra[] = lapply(extra, as.character)
  
  # read the ref data file if not already an mrs_data object
  if (is.def(w_ref) & (class(w_ref)[[1]] != "mrs_data")) {
    w_ref <- read_mrs(w_ref)
  }
  
  if (is.def(mri) & (!("niftiImage" %in% class(mri)))) mri <- readNifti(mri)
  
  if (is.def(mri_seg) & (!("niftiImage" %in% class(mri_seg)))) {
    mri_seg <- readNifti(mri_seg)
  }
  
  # TODO reorient images?
  
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
  
  # if the spectral width exceeds 20 PPM, then it should be ok
  # to decimate the signal to improve analysis speed
  if (is.null(decimate) & ((fs(metab) / matab$ft * 1e6) > 20)) {
    decimate <- TRUE
  } else {
    decimate <- FALSE
  }
  
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
  
  # end here if we're only interested in the preprocessed data
  if (preproc_only) return(metab)
  
  # TODO fitting options
  fit_res <- fit_mrs(metab, basis = basis, extra = extra)
  
  if (scale_amps) {
    if (is.def(w_ref) & is.def(mri_seg)) {
      
      if (is.null(te)) stop("te not given, amplitude scaling failed")
        
      if (is.null(tr)) stop("tr not given, amplitude scaling failed")
      
      # generate the svs voi in the segmented image space
      voi <- get_svs_voi(metab, mri_seg)
      # calculate partial volumes
      seg <- get_voi_seg(voi, mri_seg)
      # do pvc
      fit_res <- scale_amp_molal_pvc(fit_res, w_ref, seg, te, tr)
    } else if (is.def(w_ref) & !is.def(mri_seg)) {
      # do straight w scaling default LCM style
      fit_res <- scale_amp_molar(fit_res, w_ref)
    } else {
      # scale to tCr
      fit_res <- scale_amp_ratio(fit_res, "tCr")
    }
  }
  return(fit_res)
}

# TODO
# what if metab fname results in metab + ref file, eg GE data?
# options to add, fit_opts, fit_method, format (GE, Siemens etc)
# allow lists of mrs_data objects, or file paths as input
# add extra option for id vars

#' @export
svs_1h_batch_analysis <- function(metab_list, w_ref_list = NULL,
                                  mri_list = NULL, mri_seg_list = NULL,
                                  output_dir = NULL, extra = NULL, ...) {
  
  # check input it sensible
  metab_n <- length(metab_list)
  
  if (is.def(w_ref_list) & length(w_ref_list) != metab_n) {
    stop("Incorrect number of w_ref_list items.")
  }
  
  if (is.def(mri_list) & length(mri_list) != metab_n) {
    stop("Incorrect number of mri_list items.")
  }
  
  if (is.def(mri_seg_list) & length(mri_seg_list) != metab_n) {
    stop("Incorrect number of mri_seg_list items.")
  }
  
  if (is.def(output_dir) & length(output_dir) != metab_n) {
    stop("Incorrect number of output_dir_list items.")
  }
  
  if (is.def(extra) & length(extra) != metab_n) {
    stop("Incorrect number of rows in extra.")
  }
  
  # check file paths exist TODO write a function to return bad paths
  
  if (is.null(w_ref_list)) w_ref_list <- vector("list", metab_n)
  
  if (is.null(mri_seg_list)) mri_seg_list <- vector("list", metab_n)
  
  if (is.null(mri_list)) mri_list <- vector("list", metab_n)
  
  if (is.null(output_dir)) output_dir_list <- vector("list", metab_n)
  
  fit_list <- mapply(svs_1h_analysis, metab = metab_list, w_ref = w_ref_list,
                     mri_seg = mri_seg_list, mri = mri_list,
                     output_dir = output_dir_list, extra = extra,
                     MoreArgs = list(...), SIMPLIFY = FALSE, USE.NAMES = FALSE)
  
  return(fit_list)
}
