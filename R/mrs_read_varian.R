read_varian <- function(fname, extra) {
  
  # open the file
  con <- file(fname, "rb")
  
  # data file header
  nblocks   <- readBin(con, "int", size = 4, endian = "big")
  ntraces   <- readBin(con, "int", size = 4, endian = "big")
  np        <- readBin(con, "int", size = 4, endian = "big")
  ebytes    <- readBin(con, "int", size = 4, endian = "big")
  tbytes    <- readBin(con, "int", size = 4, endian = "big")
  bbytes    <- readBin(con, "int", size = 4, endian = "big")
  vers_id   <- readBin(con, "int", size = 2, endian = "big")
  status    <- readBin(con, "int", size = 2, endian = "big")
  nbheaders <- readBin(con, "int", size = 4, endian = "big")
  
  # block header
  scale     <- readBin(con, "int", size = 2, endian = "big")
  bstatus   <- readBin(con, "int", size = 2, endian = "big")
  index     <- readBin(con, "int", size = 2, endian = "big")
  mode      <- readBin(con, "int", size = 2, endian = "big")
  ctcount   <- readBin(con, "int", size = 4, endian = "big")
  lpval     <- readBin(con, "double", size = 4, endian = "big")
  rpval     <- readBin(con, "double", size = 4, endian = "big")
  lvl       <- readBin(con, "double", size = 4, endian = "big")
  tlt       <- readBin(con, "double", size = 4, endian = "big")
  
  status_bits <- intToBits(status)
  
  # warning - this bit is untested
  if (status_bits[12]) {
    raw_vec <- readBin(con, "double", size = 4, n = np / 2, endian = "big")
  } else if (status_bits[13]) {
    raw_vec <- readBin(con, "int", size = 4, n = np / 2, endian = "big")
  } else {
    raw_vec <- readBin(con, "int", size = 2, n = np / 2, endian = "big")
  }
  
  data <- raw_vec[c(TRUE, FALSE)] + 1i * raw_vec[c(FALSE, TRUE)]
  
  # close the fid file
  close(con)
  
  # read the procpar file
  fnamepp <- file.path(dirname(fname), "procpar")
  
  pp_txt <- readLines(fnamepp)
  
  # find the spectral width
  fs_line <- which(substr(pp_txt, 1, 3) == "sw ")
  fs <- as.numeric(strsplit(pp_txt[fs_line + 1], " ")[[1]][2])
  
  # find the transmitter frequency
  ft_line <- which(substr(pp_txt, 1, 5) == "sfrq ")
  ft <- as.numeric(strsplit(pp_txt[ft_line + 1], " ")[[1]][2]) * 1e6
  
  dim(data) <- c(1, 1, 1, 1, 1, 1, length(data))
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  
  ref <- def_ref()
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  te <- NA
  nuc <- def_nuc()
  
  meta <- list(Manufacturer = "Varian")
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = freq_domain, affine = NULL,
                       meta = meta, extra = extra)
  
  return(mrs_data)
}