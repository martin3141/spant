read_mrs_nifti <- function(fname) {
  
  # check the file extension is sensible
  if (stringr::str_sub(fname, -7) != ".nii.gz") {
    stop("filename argument must end in .nii.gz")
  }
  
  # get the file name of the json sidecar file
  fname_json <- stringr::str_c(stringr::str_sub(fname, 1, -7), "json")
  
  # check both files exist
  if (!file.exists(fname)) {
    cat(fname)
    stop("File not found.")
  } else if (!file.exists(fname_json)) {
    cat(fname_json)
    stop("json file not found.")
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
  data <- aperm(data, c(1, 2, 3, 5, 6, 4)) 
  
  # add a dummy dimension
  dim(data) <- c(1, dim(data))
  
  # read the json file
  json_data <- jsonlite::fromJSON(fname_json)
  
  if ("dim_5" %in% json_data) stop("NIFTI MRS non-default dimensions are not currently supported")
  if ("dim_6" %in% json_data) stop("NIFTI MRS non-default dimensions are not currently supported")
  if ("dim_7" %in% json_data) stop("NIFTI MRS non-default dimensions are not currently supported")
  
  # read voxel dimensions, dwell time and time between dynamic scans
  res <- c(NA, pixdim[2], pixdim[3], pixdim[4], pixdim[6], NA, pixdim[5])
  
  # affine and position information
  xform_mat <- RNifti::xform(nii_data)
  col_vec <- xform_mat[1:3, 1] / sum(xform_mat[1:3, 1] ^ 2) ^ 0.5 * c(-1, -1, 1)
  row_vec <- xform_mat[1:3, 2] / sum(xform_mat[1:3, 2] ^ 2) ^ 0.5 * c(-1, -1, 1)
  sli_vec <- crossprod_3d(col_vec, row_vec)
  pos_vec <- xform_mat[1:3, 4] * c(-1, -1, 1)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  te  <- json_data$EchoTime / 1e3
  ft  <- json_data$TransmitterFrequency * 1e6
  
  if (exists("SpectralWidth", where = json_data)) {
    # check json and nifti header values for the sampling frequency are consistent
    if (abs((1 / pixdim[5]) - json_data$SpectralWidth) > 0.01) {
      stop("Sampling frequencies in the json sidecar and NIFTI header differ by greater than 0.01 Hz")
    }
    res[7] <- 1 / json_data$SpectralWidth # prefer the json value if available
                                          # due to higher precision than the
                                          # NIFTI header (double vs float)
  }
 
  # TODO read from the file 
  nuc <- def_nuc()
  
  # TODO get ref from a lookup table of defaults when not in the json sidecar
  ref <- def_ref()
  
  mrs_data <- list(ft = ft, data = data, resolution = res,
                   te = te, ref = ref, nuc = nuc, row_vec = row_vec,
                   col_vec = col_vec, sli_vec = sli_vec, pos_vec = pos_vec,
                   freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}