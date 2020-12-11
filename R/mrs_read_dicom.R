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
      stop("Incorrect SOP class UID, this doesn't look like MRS data.")
    }
  } else {
    stop("DICOM MRS manufacturer format not currently supported.") 
  }
}

read_siemens_dicom <- function(fraw) {
  
  # list of tags to pull from the dicom file
  tags <- list(data = "5600,0020",
               fs = "0018,9052",
               ft = "0018,9098",
               ipp = "0020,0032",
               iop = "0020,0037",
               pixsp = "0028,0030",
               slice_t = "0018,0050",
               rows = "0028,0010",
               cols = "0028,0011",
               slices = "0018,9159",
               Npt = "0028,9002",
               te = "0018,9082")

  # read em out
  res  <- dicom_reader(fraw, tags)
  
  # sort out the right datatypes
  fs  <- readBin(res$fs, "double")
  ft  <- readBin(res$ft, "double") * 1e6
  ipp <- as.numeric(strsplit(rawToChar(res$ipp), "\\\\")[[1]])
  iop <- as.numeric(strsplit(rawToChar(res$iop), "\\\\")[[1]])
  pixsp <- as.numeric(strsplit(rawToChar(res$pixsp), "\\\\")[[1]])
  slice_t <- as.numeric(rawToChar(res$slice_t))
  rows <- readBin(res$rows, "integer", size = 2, signed = FALSE)
  cols <- readBin(res$cols, "integer", size = 2, signed = FALSE)
  slices <- readBin(res$slices, "integer", size = 2, signed = FALSE)
  N <- readBin(res$Npt, "integer")
  te <- readBin(res$te, "double") / 1e3
  
  row_ori <- iop[1:3]
  col_ori <- iop[4:6]
  sli_vec <- crossprod_3d(row_ori, col_ori)
  pos_vec <- ipp
  
  fids <- rows * cols * slices
  
  if (N * 2 * fids * 4 != length(res$data)) {
    stop("Unexpected number of data points.")
  }
  
  # read the fid points 
  raw_vec <- readBin(res$data, what = "double", n = N * 2 * fids, size = 4,
                     endian = "little")
  
  data <- raw_vec[c(TRUE, FALSE)] - 1i * raw_vec[c(FALSE, TRUE)]
  
  dim(data) <- c(N, rows, cols, slices, 1, 1, 1)
  data <- aperm(data, c(7, 2, 3, 4, 5, 6, 1))
  
  res <- c(NA, pixsp[1], pixsp[2], slice_t, 1, NA, 1 / fs)
  ref <- def_ref()
  
  # TODO determine from the data
  nuc <- def_nuc()
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  pos_vec_affine <- pos_vec + row_ori * res[3] / 2 + col_ori * res[2] / 2
  
  affine <- cbind(c(row_ori * res[3], 0),
                  c(col_ori * res[2], 0),
                  c(sli_vec * res[4], 0),
                  c(pos_vec_affine, 1))
  affine[1:2,] <- -affine[1:2,]
  
  mrs_data <- list(ft = ft, data = data, resolution = res, te = te, ref = ref, 
                   nuc = nuc, row_vec = row_ori, col_vec = col_ori,
                   sli_vec = sli_vec, pos_vec = pos_vec,
                   freq_domain = freq_domain, affine = affine)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}