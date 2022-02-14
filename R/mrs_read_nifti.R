read_mrs_nifti <- function(fname, extra) {
  
  fname_low <- tolower(fname)
  
  # check the file extension is sensible
  if (!stringr::str_ends(fname_low, ".nii.gz") &
      !stringr::str_ends(fname_low, ".nii")) {
    stop("filename argument must end in .nii.gz or .nii")
  }
  
  # check file exists
  if (!file.exists(fname)) {
    cat(fname)
    stop("File not found.")
  } 
  
  # read the nifti file
  nii_data <- RNifti::readNifti(fname)
  
  # get the dimensions
  pixdim <- nii_data$pixdim
  
  # read array values 
  data <- nii_data[]
  
  # add any missing dimensions
  if (length(dim(data)) < 6) {
    zero_dims <- rep(1, 6 - length(dim(data)))
    dim(data) <- c(dim(data), zero_dims)
  }
  
  # reorder the dimensions
  # NIFTI default dimension ordering is X, Y, Z, FID, coil, dynamic, indirect
  # spant default dimension ordering is (dummy,) X, Y, Z, dynamic, coil, FID 
  data <- aperm(data, c(1, 2, 3, 6, 5, 4)) 
  
  # add a dummy dimension
  dim(data) <- c(1, dim(data))
  
  ext_char <- RNifti::extension(readNifti(fname), 44, "character")
  
  if (is.null(ext_char)) stop("NIfTI extension header for MRS not found.")
  
  # read the json file
  json_data <- jsonlite::fromJSON(ext_char)
  
  # TODO
  # if ("dim_5" %in% names(json_data)) stop("NIfTI MRS non-default dimensions are not currently supported")
  # if ("dim_6" %in% names(json_data)) stop("NIfTI MRS non-default dimensions are not currently supported")
  # if ("dim_7" %in% names(json_data)) stop("NIfTI MRS non-default dimensions are not currently supported")
  
  # read voxel dimensions, dwell time and time between dynamic scans
  res <- c(NA, pixdim[2], pixdim[3], pixdim[4], pixdim[6], NA, pixdim[5])
  
  # affine and position information
  xform_mat <- RNifti::xform(nii_data)
  col_vec <- xform_mat[1:3, 1] / sum(xform_mat[1:3, 1] ^ 2) ^ 0.5 * c(-1, -1, 1)
  row_vec <- xform_mat[1:3, 2] / sum(xform_mat[1:3, 2] ^ 2) ^ 0.5 * c(-1, -1, 1)
  sli_vec <- crossprod_3d(row_vec, col_vec)
  pos_vec <- xform_mat[1:3, 4] * c(-1, -1, 1)
  affine  <- xform_mat
  
  attributes(affine) <- list(dim = dim(affine))
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  # if (is.null(json_data$EchoTime)) {
  #   te <- json_data$EchoTime
  # } else {
  #   te <- json_data$EchoTime / 1e3
  # }
  
  ft <- json_data$SpectrometerFrequency * 1e6
  
  # read the nucleus
  nuc <- json_data$ResonantNucleus
  
  # TODO get ref from a lookup table of defaults depending on "nuc" when not in
  # the json sidecar
  ref <- def_ref()
  
  # get all metadata
  meta <- json_data
  # remove any data that is explicitly part of the mrs_data structure
  # TODO add ref when we decide what it is called
  # meta$TxOff <- NULL
  meta$TransmitterFrequency <- NULL
  meta$ResonantNucleus <- NULL
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = freq_domain, affine = affine,
                       meta = meta, extra = extra)
  
  return(mrs_data)
}