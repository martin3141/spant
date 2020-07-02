read_mrs_nifti <- function(fname) {
  
  if (stringr::str_sub(fname, -7) != ".nii.gz") {
    stop("filename argument must end in .nii.gz")
  }
  
  # get fname of the json sidecar file
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
  data <- aperm(data, c(1, 2, 3, 5, 6, 4)) # TODO check coils and dyns are the
                                           # right way round
  
  # add a dummy dimension
  dim(data) <- c(1, dim(data))
  
  # read the json file
  json_data <- jsonlite::fromJSON(fname_json)
  
  # TODO add voxel dims
  res <- c(NA, NA, NA, NA, 1, NA, pixdim[5])
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  te  <- json_data$EchoTime / 1e3
  ft  <- json_data$TransmitterFrequency * 1e6
  
  # check json and nifti header values for the sampling frequency are consistent
  if (abs((1 / pixdim[5]) - json_data$SpectralWidth) > 0.01) {
    stop("Sampling frequencies in the json sidecar and NIFTI header differ by greater than 0.01 Hz")
  }
  
  # TODO get from a lookup table of defaults when not in the json sidecar
  ref <- def_acq_paras()$ref
  
  mrs_data <- list(ft = ft, data = data, resolution = res,
                   te = te, ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}