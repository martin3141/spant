# TODO test MEGA-PRESS and CSI
# when nechoes = 2 we're probally dealing with MEGA-PRESS data, where edited
# pairs are not interleaved but occupy the first and second half of the data

read_pfile <- function(fname, n_ref_scans = NULL, verbose, extra) {
  # check the file size
  fbytes <- file.size(fname)
  
  hdr <- read_pfile_header(fname)
  
  # override w ref scans detected in the header 
  if (!is.null(n_ref_scans)) hdr$rhuser19 <- n_ref_scans
  
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
  
  # each echo starts with a frame of zeros followed by w_ref frames followed by
  # ws frames
  expt_pts <- coils * (hdr$nframes * hdr$nechoes + hdr$nechoes) * hdr$frame_size * 2
  
  if (expt_pts != Npts) {
    warning("Unexpected number of data points.")
    cat(paste("Expecting   :", Npts, "points based on file size.\n"))
    cat(paste("Expecting   :", expt_pts, "points based on header information.\n"))
    cat(paste("Coils       :", coils, "\n"))
    cat(paste("nframes     :", hdr$nframes, "\n"))
    cat(paste("nechoes     :", hdr$nechoes), "\n")
    cat(paste("frame_size  :", hdr$frame_size), "\n")
    cat(paste("w_frames    :", hdr$rhuser19), "\n")
    cat(paste("Header rev. :", hdr$hdr_rev), "\n")
  }
  
  if (verbose) {
    cat(paste("Expecting   :", Npts, "points based on file size.\n"))
    cat(paste("Expecting   :", expt_pts, "points based on header information.\n"))
    cat(paste("Coils       :", coils, "\n"))
    cat(paste("nframes     :", hdr$nframes, "\n"))
    cat(paste("nechoes     :", hdr$nechoes), "\n")
    cat(paste("frame_size  :", hdr$frame_size), "\n")
    cat(paste("w_frames    :", hdr$rhuser19), "\n")
    cat(paste("Header rev. :", hdr$hdr_rev), "\n")
    cat(paste("Seq. name   :", hdr$seq_name), "\n")
    cat(paste("CSI dims    :", hdr$csi_dims), "\n")
    cat(paste("xcsi        :", hdr$xcsi), "\n")
    cat(paste("ycsi        :", hdr$ycsi), "\n")
    cat(paste("zcsi        :", hdr$zcsi), "\n\n")
  }
  
  data <- raw_pts[c(TRUE, FALSE)] + 1i * raw_pts[c(FALSE, TRUE)]
  
  if (hdr$csi_dims == 0) {
    # looks like SVS
    dyns <- hdr$nechoes * hdr$nframes + hdr$nechoes
    data <- array(data, dim = c(hdr$frame_size, dyns, coils, 1, 1, 1, 1))
    data <- aperm(data, c(7, 6, 5, 4, 2, 3, 1))
    
    # remove the empty frame at the start of each coil for each echo
    rem <- seq(from = 1, to = dyns, by = dyns / hdr$nechoes)
    data <- data[,,,,-rem,,,drop = FALSE]
  } else {
    # looks like MRSI
    
    # the + 1 is an empty FID at the start of each coil 
    data <- array(data, dim = c(hdr$frame_size, hdr$xcsi * hdr$ycsi + 1, coils,
                                1, 1, 1, 1))
    
    data <- aperm(data, c(7, 6, 5, 4, 2, 3, 1))
    
    # remove the empty FID
    data <- data[,,,,-1,,,drop = FALSE]
    
    # reshape
    new_dim <- dim(data)
    new_dim[2] <- hdr$xcsi
    new_dim[3] <- hdr$ycsi
    new_dim[5] <- 1
    
    dim(data) <- new_dim
    
    # phase cycling correction?
    k_sp_corr_x <- array(1, dim(data))
    k_sp_corr_x[, c(T, F),,,,,] <- -1
    data <- data * k_sp_corr_x
}
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / hdr$spec_width)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  ref <- def_ref()
  nuc <- def_nuc()
  
  
  meta <- list(EchoTime = hdr$te,
               RepetitionTime = hdr$tr,
               SequenceName = hdr$seq_name,
               Manufacturer = "GE")
  
  if (toupper(hdr$seq_name) == "PROBE-P") {
    meta <- append(meta, list(PulseSequenceType = "press"))
  }
  
  if (toupper(hdr$seq_name) == "PROBE-S") {
    meta <- append(meta, list(PulseSequenceType = "steam"))
  }
  
  mrs_data <- mrs_data(data = data, ft = hdr$ps_mps_freq / 10, resolution = res,
                       ref = ref, nuc = nuc, freq_domain = freq_domain,
                       affine = NULL, meta = meta, extra = extra)
  
  if (hdr$rhuser19 > 0) {
    # split water and metab data for each echo
    wref_inds <- rep(FALSE, Ndyns(mrs_data) / hdr$nechoes)
    wref_inds[1:hdr$rhuser] <- TRUE
    wref_inds <- rep(wref_inds, hdr$nechoes)
    
    ref_mrs <- get_dyns(mrs_data, which(wref_inds))
    metab_mrs <- get_dyns(mrs_data, which(!wref_inds))
    
    if (hdr$csi_dims > 0) {
      ref_mrs   <- img2kspace_xy(ref_mrs)
      metab_mrs <- img2kspace_xy(metab_mrs)
    }
  } else {
    ref_mrs <- NA
    metab_mrs <- mrs_data
    
    if (hdr$csi_dims > 0) {
      metab_mrs <- img2kspace_xy(metab_mrs)
    }
  }
  
  out <- list(metab = metab_mrs, ref = ref_mrs)
  class(out) <- c("list", "mrs_data")
  return(out)
}
  
read_pfile_header <- function(fname) {
  endian <- "little"
  vars <- get_pfile_vars()
  con <- file(fname, "rb")
  vars$hdr_rev <- readBin(con, "numeric", size = 4, endian = endian)
  
  loc <- get_pfile_dict(vars$hdr_rev, con)
  
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
  vars$te <- readBin(con, "int", size = 4, endian = endian) / 1e6
  
  seek(con, loc$tr)
  vars$tr <- readBin(con, "int", size = 4, endian = endian) / 1e6
  
  seek(con, loc$seq_name)
  vars$seq_name <- readBin(con, "character", size = 33)
  
  close(con)
  
  vars
}

get_pfile_vars <- function() {
  vars <- vector(mode = "list", length = 16)
  names(vars) <- c("hdr_rev", "off_data", "nechoes", "nframes", "frame_size", 
                   "rcv", "rhuser19", "spec_width", "csi_dims", "xcsi", "ycsi",
                   "zcsi", "ps_mps_freq", "te", "tr", "seq_name")
  vars
}

# following are useful resources for GE p-files
# https://github.com/SIVICLab/sivic/blob/master/libs/src/svkGEPFileReader.cc
# https://github.com/chenkonturek/MRS_MRI_libs/blob/master/MRS_lib/io/mrs_readGEpfile.m
# TARQUIN source is an ok reference for older versions.
# Ralph Noeske also suggests GERead.m from GANNET as a good place to look for
# newer versions.
# Other option is to use a hex editor, and search for known values
# https://hexed.it/

get_pfile_dict <- function(hdr_rev, con) {
  loc <- get_pfile_vars()
  
  if (floor(hdr_rev) == 11L) { # version 12
    warning("pfile version is untested :", hdr_rev)
    loc$hdr_rev     <- 0
    loc$off_data    <- 4
    loc$nechoes     <- 146
    loc$nframes     <- 150
    loc$frame_size  <- 156
    loc$rcv         <- 264
    loc$rhuser19    <- 356
    loc$spec_width  <- 432
    loc$csi_dims    <- 436
    loc$xcsi        <- 438
    loc$ycsi        <- 440
    loc$zcsi        <- 442
    loc$ps_mps_freq <- 488
    loc$te          <- 1148
    loc$tr          <- 199236
    loc$seq_name    <- 199812
  } else if (floor(hdr_rev) == 14L) {
    loc$hdr_rev     <- 0
    loc$off_data    <- 1468
    loc$nechoes     <- 70
    loc$nframes     <- 74
    loc$frame_size  <- 80
    loc$rcv         <- 200
    loc$rhuser19    <- 292
    loc$spec_width  <- 368
    loc$csi_dims    <- 372
    loc$xcsi        <- 374
    loc$ycsi        <- 376
    loc$zcsi        <- 378
    loc$ps_mps_freq <- 424
    loc$te          <- 144580
    loc$tr          <- 144572
    loc$seq_name    <- 145132
  } else if (floor(hdr_rev) == 21L) {
    loc$hdr_rev     <- 0
    loc$off_data    <- 1468
    loc$nechoes     <- 70
    loc$nframes     <- 74
    loc$frame_size  <- 80
    loc$rcv         <- 200
    loc$rhuser19    <- 292
    loc$spec_width  <- 368
    loc$csi_dims    <- 372
    loc$xcsi        <- 374
    loc$ycsi        <- 376
    loc$zcsi        <- 378
    loc$ps_mps_freq <- 424
    loc$te          <- 1212
    loc$tr          <- 148944
    loc$seq_name    <- 149520
  } else if (floor(hdr_rev) == 24L) {
    loc$hdr_rev     <- 0
    loc$off_data    <- 1468
    loc$nechoes     <- 70
    loc$nframes     <- 74
    loc$frame_size  <- 80
    loc$rcv         <- 200
    loc$rhuser19    <- 292
    loc$spec_width  <- 368
    loc$csi_dims    <- 372
    loc$xcsi        <- 374
    loc$ycsi        <- 376
    loc$zcsi        <- 378
    loc$ps_mps_freq <- 424
    # loc$te          <- 1212
    loc$te          <- 148404
    loc$tr          <- 148396
    loc$seq_name    <- 148972
  } else if (floor(hdr_rev) == 26L) {
    loc$hdr_rev     <- 0
    loc$off_data    <- 4
    loc$nechoes     <- 146
    loc$nframes     <- 150
    loc$frame_size  <- 156
    loc$rcv         <- 264
    loc$rhuser19    <- 356
    loc$spec_width  <- 432
    loc$csi_dims    <- 436
    loc$xcsi        <- 438
    loc$ycsi        <- 440
    loc$zcsi        <- 442
    loc$ps_mps_freq <- 488
    loc$te          <- 1148
    loc$tr          <- 199236
    loc$seq_name    <- 199812
  } else if (floor(hdr_rev) == 28L) {
    # some files from this version are corrupt due to a bug in the GE software
    # (relating to auto-coil selection "Air Touch")
    # don't spend ages trying to fix files that are inherently borked
    loc$hdr_rev     <- 0
    loc$off_data    <- 4
    loc$nechoes     <- 146
    loc$nframes     <- 150
    loc$frame_size  <- 156
    loc$rcv         <- 264
    loc$rhuser19    <- 356
    loc$spec_width  <- 432
    loc$csi_dims    <- 436
    loc$xcsi        <- 438
    loc$ycsi        <- 440
    loc$zcsi        <- 442
    loc$ps_mps_freq <- 488
    loc$te          <- 1148
    loc$tr          <- 207428
    loc$seq_name    <- 208004
  } else {
    warning("pfile version is untested :", hdr_rev)
    loc$hdr_rev     <- 0
    loc$off_data    <- 4
    loc$nechoes     <- 146
    loc$nframes     <- 150
    loc$frame_size  <- 156
    loc$rcv         <- 264
    loc$rhuser19    <- 356
    loc$spec_width  <- 432
    loc$csi_dims    <- 436
    loc$xcsi        <- 438
    loc$ycsi        <- 440
    loc$zcsi        <- 442
    loc$ps_mps_freq <- 488
    loc$te          <- 1148
    loc$tr          <- 207428
    loc$seq_name    <- 208004
  }
  return(loc)
}
