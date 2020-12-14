read_spar_sdat <- function(fname) {
  # generate matching SPAR and SDAT files
  ext <- stringr::str_sub(fname, -5)
  name <- stringr::str_sub(fname, 1, -6)
  
  if ( ext == ".SPAR" ) {
    spar <- fname
    sdat <- paste0(name, ".SDAT")
  } else if ( ext == ".SDAT" ) {
    sdat <- fname
    spar <- paste0(name, ".SPAR")
  } else if ( ext == ".spar" ) {
    spar <- fname
    sdat <- paste0(name, ".sdat")
  } else if ( ext == ".sdat" ) {
    sdat <- fname
    spar <- paste0(name, ".spar")
  } else {
    stop("Incorrect file extension.")
  }
  
  # check both files exist
  if (!file.exists(spar)) {
    cat(spar)
    stop("SPAR file not found.")
  } else if (!file.exists(sdat)) {
    cat(sdat)
    stop("SDAT file not found.")
  }
  
  paras <- utils::read.delim(spar, sep = ":", comment.char = "!",
                             header = FALSE, strip.white = TRUE,
                             stringsAsFactors = FALSE)
                    
  #N <- as.integer(paras$V2[which(paras$V1 == "samples")])
  N <- as.numeric(paras$V2[which(paras$V1 == "dim1_pnts")])
  dyns <- as.integer(paras$V2[which(paras$V1 == "rows")])
  ft <- as.numeric(paras$V2[which(paras$V1 == "synthesizer_frequency")])
  fs <- as.numeric(paras$V2[which(paras$V1 == "sample_frequency")])
  te <- as.numeric(paras$V2[which(paras$V1 == "echo_time")]) * 1e-3
  ap_oc <- as.numeric(paras$V2[which(paras$V1 == "ap_off_center")])
  lr_oc <- as.numeric(paras$V2[which(paras$V1 == "lr_off_center")])
  cc_oc <- as.numeric(paras$V2[which(paras$V1 == "cc_off_center")])
  ap_an <- as.numeric(paras$V2[which(paras$V1 == "ap_angulation")])
  lr_an <- as.numeric(paras$V2[which(paras$V1 == "lr_angulation")])
  cc_an <- as.numeric(paras$V2[which(paras$V1 == "cc_angulation")])
  ap_size <- as.numeric(paras$V2[which(paras$V1 == "ap_size")])
  lr_size <- as.numeric(paras$V2[which(paras$V1 == "lr_size")])
  cc_size <- as.numeric(paras$V2[which(paras$V1 == "cc_size")])
  sli_thick <- as.numeric(paras$V2[which(paras$V1 == "slice_thickness")])
  pe_fov <- as.numeric(paras$V2[which(paras$V1 == "phase_encoding_fov")])
  cols <- as.numeric(paras$V2[which(paras$V1 == "dim2_pnts")])
  rows <- as.numeric(paras$V2[which(paras$V1 == "dim3_pnts")])
  slices <- as.numeric(paras$V2[which(paras$V1 == "nr_of_slices_for_multislice")])
  #cols <- as.numeric(paras$V2[which(paras$V1 == "SUN_dim2_pnts")])
  #rows <- as.numeric(paras$V2[which(paras$V1 == "SUN_dim3_pnts")])
  
  # May be useful...
  # slices <- as.numeric(paras$V2[which(paras$V1 == "nr_of_slices_for_multislice")])
  # avs <- as.integer(paras$V2[which(paras$V1 == "averages")])
  # dim1_pts <- as.numeric(paras$V2[which(paras$V1 == "dim1_pnts")])
  # dim1_pts <- as.numeric(paras$V2[which(paras$V1 == "SUN_dim1_pnts")])
  # nuc <- as.numeric(paras$V2[which(paras$V1 == "nucleus")])
  
  # the following can be true for non localised acquisitions 
  if (length(ap_an) == 0) ap_an <- 0
  if (length(lr_an) == 0) lr_an <- 0
  if (length(cc_an) == 0) cc_an <- 0
  if (length(ap_size) == 0) ap_size <- 0
  if (length(lr_size) == 0) lr_size <- 0
  if (length(cc_size) == 0) cc_size <- 0
  if (length(ap_oc) == 0) ap_oc <- 0
  if (length(lr_oc) == 0) lr_oc <- 0
  if (length(cc_oc) == 0) cc_oc <- 0
  
  true_row   <- c(1,0,0)
  true_col   <- c(0,1,0)
  true_slice <- c(0,0,1)
  
  row_ori <- rotate_vec(true_row, true_slice, cc_an * pi / 180)
  row_ori <- rotate_vec(row_ori, true_col, ap_an * pi / 180)
  row_ori <- rotate_vec(row_ori, true_row, lr_an * pi / 180)
  
  col_ori <- rotate_vec(true_col, true_slice, cc_an * pi / 180)
  col_ori <- rotate_vec(col_ori, true_col, ap_an * pi / 180)
  col_ori <- rotate_vec(col_ori, true_row, lr_an * pi / 180)
  
  sli_vec <- crossprod_3d(row_ori, col_ori)
  
  pos_vec <- c(lr_oc, ap_oc, cc_oc)
  
  data_vec <- read_sdat(sdat)
  
  # SVS or MRSI?
  if ((rows == 1) & (cols == 1)) {
    row_dim   <- ap_size
    col_dim   <- lr_size
    slice_dim <- cc_size
  } else {
    dyns <- 1
    row_dim   <- pe_fov / cols
    col_dim   <- pe_fov / cols
    slice_dim <- sli_thick
    pos_vec <- (pos_vec - col_ori * row_dim * 0.5 * (rows - 1) - 
                row_ori * col_dim * 0.5 * (cols - 1))
  }
  
  #data <- array(data_vec,dim = c(1, cols, rows, slices, N, 1, dyns)) 
  data <- array(data_vec,dim = c(N, cols, rows, slices, dyns, 1, 1)) 
  data <- aperm(data,c(6, 2, 3, 4, 5, 7, 1))
  
  res <- c(NA, row_dim, col_dim, slice_dim, 1, NA, 1 / fs)
  ref <- def_ref()
  
  # TODO get from the data file
  nuc <- def_nuc()
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, te = te,
                       ref = ref, nuc = nuc, freq_domain = freq_domain,
                       affine = NULL, meta = NULL)
  
  return(mrs_data)
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