# TODO test MEGA-PRESS and CSI
read_pfile <- function(fname) {
  # check the file size
  fbytes <- file.size(fname)
  
  hdr <- read_pfile_header(fname)
  con <- file(fname, "rb")
  seek(con, hdr$off_data)
  endian <- "little"
  
  # calculate number of data points from file size
  Npts <- (fbytes - hdr$off_data) / 4
  raw_pts <- readBin(con, "int", n = Npts, size = 4, endian = endian)
  close(con)
  
  # calculate coil elements
  coils <- 0
  for (n in seq(1, 8, 2)) {
    if ((hdr$rcv[n] != 0) || (hdr$rcv[n + 1] != 0)) {
      coils <- coils + 1 + hdr$rcv[n + 1] - hdr$rcv[n]
    }
  }
  
  if (coils == 0) coils <- 1
  
  expt_pts <- coils * (hdr$nframes * hdr$nechoes + 1) * hdr$frame_size * 2
  
  if (expt_pts != Npts) warning("Unexpected number of data points.")
  
  data <- raw_pts[c(TRUE, FALSE)] + 1i * raw_pts[c(FALSE, TRUE)]
  
  data <- array(data, dim = c(hdr$frame_size, hdr$nechoes * hdr$nframes + 1, coils, 1, 1, 1, 1))
  data <- aperm(data, c(7,6,5,4,2,3,1))
  
  # remove the empty frame at the start of each coil
  data <- data[,,,,-1,,,drop = FALSE]
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / hdr$spec_width)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  ref <- def_acq_paras()$ref
  
  mrs_data <- list(ft = hdr$ps_mps_freq / 10, data = data, resolution = res,
                   te = hdr$te, ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  
  if (hdr$rhuser19 > 0) {
    ref_mrs <- get_dyns(mrs_data, 1:hdr$rhuser)
    metab_mrs <- get_dyns(mrs_data, (1 + hdr$rhuser):dyns(mrs_data))
  } else {
    ref_mrs <- NA
    metab_mrs <- mrs_data
  }
  
  list(metab = metab_mrs, ref = ref_mrs)
}
  
read_pfile_header <- function(fname) {
  endian <- "little"
  vars <- get_pfile_vars()
  con <- file(fname, "rb")
  vars$hdr_rev <- readBin(con, "numeric", size = 4, endian = endian)
  
  loc <- get_pfile_dict(vars$hdr_rev)
  
  seek(con, loc$off_data)
  vars$off_data <- readBin(con, "int", size = 4, endian = endian)
  
  seek(con, loc$nechoes)
  vars$nechoes <- readBin(con, "int", size = 2, endian = endian)
  
  seek(con, loc$nframes)
  vars$nframes <- readBin(con, "int", size = 2, endian = endian)
  
  seek(con, loc$frame_size)
  vars$frame_size <- readBin(con, "int", size = 2, signed = FALSE, 
                             endian = endian)
  
  seek(con, loc$rcv)
  vars$rcv <- readBin(con, "int", n = 8, size = 2, endian = endian)
  
  seek(con, loc$rhuser19)
  vars$rhuser19 <- readBin(con, "numeric", size = 4, endian = endian)
  
  seek(con, loc$spec_width)
  vars$spec_width <- readBin(con, "numeric", size = 4, endian = endian)
  
  seek(con, loc$csi_dims)
  vars$csi_dims <- readBin(con, "int", size = 2, endian = endian)
  
  seek(con, loc$xcsi)
  vars$xcsi <- readBin(con, "int", size = 2, endian = endian)
  
  seek(con, loc$ycsi)
  vars$ycsi <- readBin(con, "int", size = 2, endian = endian)
  
  seek(con, loc$zcsi)
  vars$zcsi <- readBin(con, "int", size = 2, endian = endian)
  
  seek(con, loc$ps_mps_freq)
  # read as int
  ps_mps_freq_bits <- intToBits(readBin(con, "int", size = 4, endian = endian))
  # convert to uint
  vars$ps_mps_freq <- sum(2^.subset(0:31, as.logical(ps_mps_freq_bits)))
  
  seek(con, loc$te)
  vars$te <- readBin(con, "int", size = 4, endian = endian)
  
  close(con)
  
  vars
}

get_pfile_vars <- function() {
  vars <- vector(mode = "list", length = 14)
  names(vars) <- c("hdr_rev", "off_data", "nechoes", "nframes", "frame_size", 
                   "rcv", "rhuser19", "spec_width", "csi_dims", "xcsi", "ycsi",
                   "zcsi", "ps_mps_freq", "te")
  vars
}

get_pfile_dict <- function(hdr_rev) {
  loc <- get_pfile_vars()
  
  if (floor(hdr_rev) > 25) {
    loc$hdr_rev <- 0
    loc$off_data <- 4
    loc$nechoes <- 146
    loc$nframes <- 150
    loc$frame_size <- 156
    loc$rcv <- 264
    loc$rhuser19 <- 356
    loc$spec_width <- 432
    loc$csi_dims <- 436
    loc$xcsi <- 438
    loc$ycsi <- 440
    loc$zcsi <- 442
    loc$ps_mps_freq <- 488
    loc$te <- 1148
  } else if ((floor(hdr_rev) > 11) && (floor(hdr_rev) < 25)) {
    loc$hdr_rev <- 0
    loc$off_data <- 1468
    loc$nechoes <- 70
    loc$nframes <- 74
    loc$frame_size <- 80
    loc$rcv <- 200
    loc$rhuser19 <- 292
    loc$spec_width <- 368
    loc$csi_dims <- 372
    loc$xcsi <- 374
    loc$ycsi <- 376
    loc$zcsi <- 378
    loc$ps_mps_freq <- 424
    loc$te <- 1212
  } else {
    stop(paste("Error, pfile version not supported :", hdr_rev))
  }
  loc
}
