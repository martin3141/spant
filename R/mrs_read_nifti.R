read_mrs_nifti <- function(fname) {
  
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
  
  res <- c(NA, NA, NA, NA, 1, NA, pixdim[5])
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  te <- json_data$EchoTime / 1e3
  
  proton_gr <- 42.5774785182e6 # TODO other nuclei
  ft <- json_data$MagneticFieldStrength * proton_gr
  ref <- def_acq_paras()$ref
  
  mrs_data <- list(ft = ft, data = data, resolution = res,
                   te = te, ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}

# this is slow and not used, but kept as a reference
vaxf2numeric <- function(raw) {
  sign  <- rawShift(raw[2] & as.raw(0x80), -7)
  sign  <- readBin(sign, "integer", size = 1, signed = F)
  expon <- readBin(rawShift(raw[2] & as.raw(0x7f), 1), "integer", size = 1,
                   signed = F)
  
  expon <- expon + readBin(rawShift(raw[1] & as.raw(0x80), -7), "integer",
                           size = 1, signed = F)
  
  frac  <- bitwShiftL(readBin(raw[1] & as.raw(0x7f), "integer", size = 1,
                              signed = F), 16)
  
  frac  <- frac + bitwShiftL(readBin(raw[4], "integer", size = 1,
                             signed = F), 8)
  
  frac  <- frac + readBin(raw[3], "integer", size = 1, signed = F)
  
  if (0 < expon) {
    val <- ((-1) ^ sign) * (0.5 + (frac / 16777216)) * (2 ^ (expon - 128))
  } else if ((expon == 0) & (sign == 0)) {
    val <- 0
  } else {
    val <- 0
    warning("Unusual VAX number found, corrupted file?")
  }
  val
}

# this is slow and not used, but kept as a reference
read_sdat_slow <- function(fname) {
  fbytes <- file.size(fname)
  Npts <- fbytes / 4
  raw <- readBin(fname, "raw", fbytes)
  vec <- rep(NA, Npts)
  for (n in 1:Npts) {
    fpnt <- (n - 1) * 4 + 1
    vec[n] <- vaxf2numeric(raw[fpnt:(fpnt + 4)])
  }
  vec[seq(1, Npts, 2)] - vec[seq(2, Npts, 2)] * 1i
}

read_sdat <- function(fname) {
  fbytes <- file.size(fname)
  Npts <- fbytes / 4
  raw <- readBin(fname,"raw",fbytes)
  # reorder bytes
  raw <- raw[c(rbind(seq(3, fbytes, 4), seq(4, fbytes, 4), 
                     seq(1, fbytes, 4), seq(2, fbytes, 4)))]
  
  vec <- readBin(raw, "double", size = 4, endian = "little", n = Npts) / 4
  vec[seq(1, Npts, 2)] - vec[seq(2, Npts, 2)] * 1i
}