read_dicom <- function(fname, verbose) {
  
  # read full file as raw 
  fsize <- file.info(fname)$size
  fraw  <- readBin(fname, "raw", n = as.integer(fsize), endian = "little")
  
  # read the SOP class and manufacturer tags
  res <- dicom_reader(fraw, tags = list(sop_class_uid = "0008,0016",
                                        manuf         = "0008,0070"))
  
  sop_class_uid <- rawToChar(res$sop_class_uid)
  manuf         <- rawToChar(res$manuf)
  
  if (grepl("SIEMENS", manuf)) {
    if (sop_class_uid == "1.3.12.2.1107.5.9.1") {
      # SiemensPrivateCSA Non-ImageStorage - AKA ima format
      return(read_ima(fraw))
    } else if (sop_class_uid == "1.2.840.10008.5.1.4.1.1.4.2") { 
      # MR Spectroscopy Storage
      return(read_siemens_dicom(fraw))
    } else {
      stop(paste0("Unsupported SOP class UID : ", sop_class_uid,
                  ". This doesn't look like MRS data."))
    }
  } else if (grepl("Philips", manuf)) {
    if (sop_class_uid == "1.3.46.670589.11.0.0.12.1") {
      return(read_philips_priv_dicom(fraw))
    } else if (sop_class_uid == "1.2.840.10008.5.1.4.1.1.4.2") {
      return(read_philips_dicom(fraw))
    } else {
      stop(paste0("Unsupported SOP class UID : ", sop_class_uid,
                  ". This doesn't look like MRS data."))
    }
  } else {
    stop(paste0("DICOM MRS manufacturer format: \"", trimws(manuf),
                "\" is not currently supported."))
  }
}

read_siemens_dicom <- function(fraw) {
  
  # list of tags to pull from the dicom file
  tags <- list(data    = "5600,0020",
               fs      = "0018,9052",
               ft      = "0018,9098",
               ipp     = "0020,0032",
               iop     = "0020,0037",
               pixsp   = "0028,0030",
               slice_t = "0018,0050",
               rows    = "0028,0010",
               cols    = "0028,0011",
               slices  = "0018,9159",
               Npt     = "0028,9002",
               te      = "0018,9082")

  # read em out
  dcm_res  <- dicom_reader(fraw, tags)
  
  # sort out the right datatypes
  fs  <- readBin(dcm_res$fs, "double")
  ft  <- readBin(dcm_res$ft, "double") * 1e6
  ipp <- as.numeric(strsplit(rawToChar(dcm_res$ipp), "\\\\")[[1]])
  iop <- as.numeric(strsplit(rawToChar(dcm_res$iop), "\\\\")[[1]])
  pixsp <- as.numeric(strsplit(rawToChar(dcm_res$pixsp), "\\\\")[[1]])
  slice_t <- as.numeric(rawToChar(dcm_res$slice_t))
  rows <- readBin(dcm_res$rows, "integer", size = 2, signed = FALSE)
  cols <- readBin(dcm_res$cols, "integer", size = 2, signed = FALSE)
  slices <- readBin(dcm_res$slices, "integer", size = 2, signed = FALSE)
  N <- readBin(dcm_res$Npt, "integer")
  te <- readBin(dcm_res$te, "double") / 1e3
  
  row_ori <- iop[1:3]
  col_ori <- iop[4:6]
  sli_vec <- crossprod_3d(row_ori, col_ori)
  pos_vec <- ipp
  
  fids <- rows * cols * slices
  
  if (N * 2 * fids * 4 != length(dcm_res$data)) {
    stop("Unexpected number of data points.")
  }
  
  # read the fid points 
  raw_vec <- readBin(dcm_res$data, "double", length(dcm_res$data) / 4, size = 4,
                     endian = "little")
  
  data <- raw_vec[c(TRUE, FALSE)] - 1i * raw_vec[c(FALSE, TRUE)]
  
  dim(data) <- c(N, rows, cols, slices, 1, 1, 1)
  data <- aperm(data, c(7, 2, 3, 4, 5, 6, 1))
  
  res <- c(NA, pixsp[2], pixsp[1], slice_t, 1, NA, 1 / fs)
  ref <- def_ref()
  
  # TODO determine from the data
  nuc <- def_nuc()
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  pos_vec_affine <- pos_vec
  
  affine <- cbind(c(row_ori * res[2], 0),
                  c(col_ori * res[3], 0),
                  c(sli_vec * res[4], 0),
                  c(pos_vec_affine, 1))
  affine[1:2,] <- -affine[1:2,]
  
  meta <- list(EchoTime = te)
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = freq_domain, affine = affine,
                       meta = meta)
  
  return(mrs_data)
}

read_philips_dicom <- function(fraw) {
  
  # list of tags to pull from the dicom file
  tags <- list(data    = "5600,0020",
               fs      = "0018,9052",
               ft      = "0018,9098",
               ipp     = "0020,0032",
               iop     = "0020,0037",
               pixsp   = "0028,0030",
               slice_t = "0018,0050",
               rows    = "0028,0010",
               cols    = "0028,0011",
               slices  = "0018,9159",
               Nframes = "0028,0008",
               Npt     = "0028,9002",
               te      = "0018,9082")

  # read em out
  dcm_res  <- dicom_reader(fraw, tags)
  
  # sort out the right datatypes
  fs  <- readBin(dcm_res$fs, "double")
  ft  <- readBin(dcm_res$ft, "double") * 1e6
  ipp <- as.numeric(strsplit(rawToChar(dcm_res$ipp), "\\\\")[[1]])
  iop <- as.numeric(strsplit(rawToChar(dcm_res$iop), "\\\\")[[1]])
  pixsp <- as.numeric(strsplit(rawToChar(dcm_res$pixsp), "\\\\")[[1]])
  slice_t <- as.numeric(rawToChar(dcm_res$slice_t))
  rows <- readBin(dcm_res$rows, "integer", size = 2, signed = FALSE)
  cols <- readBin(dcm_res$cols, "integer", size = 2, signed = FALSE)
  slices <- readBin(dcm_res$slices, "integer", size = 2, signed = FALSE)
  dyns <- as.numeric(rawToChar(dcm_res$Nframes))
  N <- readBin(dcm_res$Npt, "integer")
  te <- readBin(dcm_res$te, "double") / 1e3
  
  row_ori <- iop[1:3]
  col_ori <- iop[4:6]
  sli_vec <- crossprod_3d(row_ori, col_ori)
  pos_vec <- ipp
  
  fids <- rows * cols * slices * dyns
  
  if (N * 2 * fids * 4 != length(dcm_res$data)) {
    warning("Unexpected number of data points.")
  }
  
  # read the fid points 
  raw_vec <- readBin(dcm_res$data, "double", length(dcm_res$data) / 4, size = 4,
                     endian = "little")
  
  data <- raw_vec[c(TRUE, FALSE)] - 1i * raw_vec[c(FALSE, TRUE)]
  
  dim(data) <- c(N, rows, cols, slices, dyns, 1, 1)
  data <- aperm(data, c(7, 2, 3, 4, 5, 6, 1))
  
  res <- c(NA, pixsp[2], pixsp[1], slice_t, 1, NA, 1 / fs)
  ref <- def_ref()
  
  # TODO determine from the data
  nuc <- def_nuc()
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  pos_vec_affine <- pos_vec
  
  affine <- cbind(c(row_ori * res[2], 0),
                  c(col_ori * res[3], 0),
                  c(sli_vec * res[4], 0),
                  c(pos_vec_affine, 1))
  affine[1:2,] <- -affine[1:2,]
  
  meta <- list(EchoTime = te)
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = freq_domain, affine = affine,
                       meta = meta)
  
  return(mrs_data)
}
  

read_philips_priv_dicom <- function(fraw) {
  
  # list of tags to pull from the dicom file
  tags <- list(data    = "2005,1270",
               fs      = "2005,1030",
               ft      = "2001,1083",
               te      = "2005,1310")
  
  dcm_res  <- dicom_reader(fraw, tags)
  
  fs  <- readBin(dcm_res$fs, "double", size = 4, n = 2)[1]
  ft  <- as.numeric(rawToChar(dcm_res$ft)) * 1e6
  te  <- readBin(dcm_res$te, "double", size = 4) / 1e3
  
  raw_vec <- readBin(dcm_res$data, "double", length(dcm_res$data) / 4, size = 4)
  data    <- raw_vec[c(TRUE, FALSE)] - 1i * raw_vec[c(FALSE, TRUE)]
  
  N <- length(data)
  
  dim(data) <- c(N, 1, 1, 1, 1, 1, 1)
  data <- aperm(data, c(7, 2, 3, 4, 5, 6, 1))
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  ref <- def_ref()
  
  # TODO determine from the data
  nuc <- def_nuc()
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  meta <- list(EchoTime = te)
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = freq_domain, affine = NULL,
                       meta = meta)
  return(mrs_data)
}