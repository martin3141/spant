#' Read MRS data from a file.
#' @param fname The filename of the dpt format MRS data.
#' @param format A string describing the data format. May be one of the 
#' following : "spar_sdat", "dpt".
#' @return An MRS data object.
#' @examples
#' fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
#' mrs_data <- read_mrs(fname, format = "spar_sdat")
#' print(mrs_data)
#' @export
read_mrs <- function(fname, format) {
  if (format == "spar_sdat") {
    return(read_spar_sdat(fname))
  } else if (format == "dpt") {
    return(read_mrs_dpt(fname))
  } else {
    stop("Unrecognised file format.")
  }
}

#' Read MRS data stored in dangerplot (dpt) v3 format.
#' @param fname The filename of the dpt format MRS data.
#' @return An MRS data object.
#' @examples
#' \dontrun{
#' mrs_data <- read_mrs_dpt(system.file("extdata","svs.dpt",package="spant"))
#' }
read_mrs_dpt <- function(fname) {
  header <- utils::read.table(fname, nrows = 15, as.is = TRUE)
  
  # Check dpt version number
  dpt_ver <- header$V2[1]
  if (dpt_ver != "3.0") {
    stop("Error, dangerplot version is not supported (!=3.0).")
  }
  
  N <- as.integer(header$V2[2])
  fs <- as.double(header$V2[3])
  ft <- as.double(header$V2[4])
  phi0 <- as.double(header$V2[5])
  phi1 <- as.double(header$V2[6])
  ref <- as.double(header$V2[7])
  te <- as.double(header$V2[8])
  rows <- as.integer(header$V2[9])
  cols <- as.integer(header$V2[10])
  slices <- as.integer(header$V2[11])
  pix_sp <- header$V2[12]
  if (pix_sp == "Unknown") {
    row_dim <- NA
    col_dim <- NA
  } else {
    row_dim <- as.double(strsplit(pix_sp, "\\\\")[[1]][1])
    col_dim <- as.double(strsplit(pix_sp, "\\\\")[[1]][2])
  }
  slice_dim_str <- header$V2[13]
  if (slice_dim_str == "Unknown") {
    slice_dim <- NA
  } else {
    slice_dim <- as.double(slice_dim_str)
  }
  
  if (header$V2[14] == "Unknown") {
    IOP <- NA
  } else {
    IOP <- as.double(strsplit(header$V2[14], "\\\\")[[1]])
  }
  
  if (header$V2[15] == "Unknown") {
    IPP <- NA
  } else {
    IPP <- as.double(strsplit(header$V2[15], "\\\\")[[1]])
  }
  
  if (!is.na(IOP[1])) {
    row_vec <- IOP[1:3]
    col_vec <- IOP[4:6]
  } else {
    row_vec = NA  
    col_vec = NA  
  }
  pos_vec <- IPP
  
  # read the data points  
  raw_data <- utils::read.table(fname, skip = 16, as.is = TRUE)
  raw_data_cplx <- raw_data$V1 + raw_data$V2 * 1i
  # construct the data array
  data_arr <- as.array(raw_data_cplx)
  
  # TODO - special case for Philips fMRS
  # (ws,w), x, y, z, t, coil, spec
  dim(data_arr) <- c(1, N, rows, cols, slices, 1, 1)
  data_arr = aperm(data_arr,c(1, 4, 3, 5, 6, 7, 2))
  
  if (dim(data_arr)[2] > 1 && dim(data_arr)[3] == 1) {
    warning("Data is 1D, assuming dynamic MRS format.")
    data_arr = aperm(data_arr,c(1, 5, 3, 4, 2, 6, 7))
  }
  
  # resolution information
  # x, y, z, t, coil, spec
  res <- c(NA, row_dim, col_dim, slice_dim, 1, NA, 1 / fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)
  
  mrs_data <- list(ft = ft, data = data_arr, resolution = res, te = te,
                   ref = ref, row_vec = row_vec, col_vec = col_vec,
                   pos_vec = pos_vec, freq_domain = freq_domain)
  class(mrs_data) <- "mrs_data"
  mrs_data
}

#' Read MRS data using the TARQUIN software package.
#' @param fname The filename containing the MRS data.
#' @param fname_ref A second filename containing reference MRS data.
#' @param format The format of the MRS data. Can be one of the following:
#' siemens, philips, ge, dcm, dpt, rda, lcm, varian, bruker, jmrui_txt.
#' @param id An optional ID string.
#' @param group An optional group string.
#' @return MRS data object.
#' @examples
#' fname <- system.file("extdata","philips_spar_sdat_WS.SDAT",package="spant")
#' \dontrun{
#' mrs_data <- read_mrs_tqn(fname, format="philips")
#' }
#' @export
read_mrs_tqn <- function(fname, fname_ref = NA, format, id = NA, group = NA) {
  # check the input file exists
  if (!file.exists(fname)) {
    print(fname)
    stop("Error, above input file does not exist.")    
  }
  
  # specify some temp file names
  ws_fname <- tempfile()
  ws_fname <- gsub(' ', '" "', ws_fname) # this is for spaces on windows
  w_fname <- tempfile()
  w_fname <- gsub(' ', '" "', w_fname) # this is for spaces on windows
  fname <- gsub(' ', '" "', fname) # this is for spaces on windows
  cmd = paste(getOption("spant.tqn_cmd"), "--input", fname, "--format", format,
                        "--write_raw_v3", ws_fname, "--write_raw_w_v3",
                        w_fname, "--rw_only", "true","--dyn_av","none", 
                        "--dyn_av_w", "none") #,"2>&1")
  
  if (!is.na(fname_ref)) {
    if (!file.exists(fname_ref)) {
      print(fname_ref)
      stop("Error, above input file does not exist.")    
    }
    cmd = paste(cmd, "--input_w", fname_ref)
  }
  
  #cmd = as.character(cat(cmd))
  #print(class(cmd))
  #print(cmd)
  res = system(cmd, intern = TRUE)
  
  if (!file.exists(ws_fname)) {
    print(res)
    print(cmd)
    stop("Error loading data with above TARQUIN command.")
  }
  
  main <- read_mrs_dpt(ws_fname)
  
  if (is.na(id)) {
    id = fname
  }
  
  main$fname = fname
  main$fname_ref = fname_ref
  main$id = id
  main$group = group
  
  if (file.exists(w_fname)) {
    ref <- read_mrs_dpt(w_fname)
    main$data <- combine_metab_ref(main, ref)
    #main$data <- abind::abind(main$data, ref$data, along=1)
  }
  
  return(main)
}

#' Write MRS data object to file in dangerplot (dpt) v2 format.
#' @param fname The filename of the output dpt format MRS data.
#' @param mrs_data Object to be written to file.
#' @examples
#' \dontrun{
#' mrs_data <- write_mrs_dpt_v2("my_mrs_data.dpt", my_mrs_data)
#' }
#' @export
write_mrs_dpt_v2 <- function(fname, mrs_data) {
  sig <- mrs_data$data[1, 1, 1, 1, 1, 1,]
  N <- length(sig)
  fs <- 1 / mrs_data$resolution[7]
  ft <-  mrs_data$ft
  ref <- mrs_data$ref
  te <- mrs_data$te
  sink(fname)
  cat("Dangerplot_version\t2.0\n")
  cat(paste("Number_of_points\t", N, "\n", sep = ""))
  cat(paste("Sampling_frequency\t", fs, "\n", sep = ""))
  cat(paste("Transmitter_frequency\t", ft, "\n", sep = ""))
  cat("Phi0\t0.0\n")
  cat("Phi1\t0.0\n")
  cat(paste("PPM_reference\t", ref, "\n", sep = ""))
  cat(paste("Echo_time\t", te, "\n", sep = ""))
  cat("Real_FID\tImag_FID\n")
  for (n in 1:N) {
    cat(paste(format(Re(sig[n]), scientific = TRUE), "\t", format(Im(sig[n]),
              scientific = TRUE), '\n', sep = ""))
  }
  sink()
}

# stolen from interweb
write.mat <- function(mat, codes, sep = "", ...) {
  s <- do.call(sprintf, unname(c(paste(codes, collapse = ""),
                                 as.data.frame(mat))))
  if (length(list(...)) > 0) cat(s, sep = sep, ...) else s
}

write_mrs_lcm_raw <- function(fname, mrs_data) {
  sig <- mrs_data$data[1, 1, 1, 1, 1, 1,]
  N <- length(sig)
  sink(fname)
  cat(" $NMID\n")
  cat(" ID='Simulated Data', FMTDAT='(2E15.6)'\n")
  cat(" VOLUME=1\n")
  cat(" TRAMP=1\n")
  cat(" $END\n")
  for (n in 1:N) {
    cat(" ")
    cat(noquote(formatC(c(Re(sig[n]), Im(sig[n])), width = 14, format = "E",
                          digits = 6)))
    cat("\n")
  }
  sink()
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
                    
  N <- as.integer(paras$V2[which(paras$V1 == "samples")])
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
  cols <- as.numeric(paras$V2[which(paras$V1 == "SUN_dim2_pnts")])
  rows <- as.numeric(paras$V2[which(paras$V1 == "SUN_dim3_pnts")])
  
  # May be useful...
  # slices <- as.numeric(paras$V2[which(paras$V1 == "nr_of_slices_for_multislice")])
  # avs <- as.integer(paras$V2[which(paras$V1 == "averages")])
  # dim1_pts <- as.numeric(paras$V2[which(paras$V1 == "dim1_pnts")])
  # dim1_pts <- as.numeric(paras$V2[which(paras$V1 == "SUN_dim1_pnts")])
  # nuc <- as.numeric(paras$V2[which(paras$V1 == "nucleus")])
  
  true_row   <- c(1,0,0)
  true_col   <- c(0,1,0)
  true_slice <- c(0,0,1)
  
  row_ori <- rotate_vec(true_row, true_slice, cc_an * pi / 180)
  row_ori <- rotate_vec(row_ori, true_col, ap_an * pi / 180)
  row_ori <- rotate_vec(row_ori, true_row, lr_an * pi / 180)
  
  col_ori <- rotate_vec(true_col, true_slice, cc_an * pi / 180)
  col_ori <- rotate_vec(col_ori, true_col, ap_an * pi / 180)
  col_ori <- rotate_vec(col_ori, true_row, lr_an * pi / 180)
  
  pos_vec <- c(lr_oc, ap_oc, cc_oc)
  
  data_vec <- read_sdat(sdat)
  
  # TODO update dims for MRSI
  data <- array(data_vec,dim = c(1, 1, 1, 1, N, 1, dyns)) 
  data = aperm(data,c(1, 2, 3, 4, 7, 6, 5))
  
  # SVS or MRSI?
  if ((rows == 1) & (cols == 1)) {
    row_dim   <- lr_size
    col_dim   <- ap_size
    slice_dim <- cc_size
  } else {
    row_dim   <- pe_fov / cols
    col_dim   <- pe_fov / cols
    slice_dim <- sli_thick
    pos_vec <- (pos_vec - col_ori * row_dim * 0.5 * (rows - 1) - 
                row_ori * col_dim * 0.5 * (cols - 1))
  }
  
  res <- c(NA, row_dim, col_dim, slice_dim, 1, NA, 1 / fs)
  ref <- get_def_acq_paras()$ref
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  mrs_data <- list(ft = ft, data = data, resolution = res, te = te, ref = ref, 
                   row_vec = row_ori, col_vec = col_ori, pos_vec = pos_vec, 
                   freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}

#' @export
read_list_data <- function(fname, ft = 127e6, fs = 2000, ref = 4.65) {
  # generate matching data and list files
  ext <- stringr::str_sub(fname, -5)
  name <- stringr::str_sub(fname, 1, -6)
  
  if ( ext == ".list" ) {
    list <- fname
    data <- paste0(name, ".data")
  } else if ( ext == ".data" ) {
    data <- fname
    list <- paste0(name, ".list")
  } else {
    stop("Incorrect file extension.")
  }
  
  # check both files exist
  if (!file.exists(list)) {
    cat(list)
    stop("list file not found.")
  } else if (!file.exists(data)) {
    cat(data)
    stop("data file not found.")
  }
  
  # read list file as text
  txt <- as.array(readLines(list))
  
  N_txt <- ".    0    0    0  F-resolution"
  N_ind <- which(apply(txt, 1, startsWith, N_txt))
  N <- as.numeric(strsplit(txt[N_ind], ":")[[1]][2])
  
  data_ind_start_txt <- "# === START OF DATA VECTOR INDEX"
  data_ind_start <- which(apply(txt, 1, startsWith, data_ind_start_txt))
  data_ind_end_txt <- "# === END OF DATA VECTOR INDEX"
  data_ind_end <- which(apply(txt, 1, startsWith, data_ind_end_txt))
  data_ind_tab <- read.table(text = txt[(data_ind_start + 3):(data_ind_end - 1)])
  col_names <- strsplit(txt[data_ind_start + 2], "\\s+")[[1]][2:22]
  colnames(data_ind_tab) <- col_names
  
  fid_num <- nrow(data_ind_tab)
  chans <- max(data_ind_tab$chan) + 1
  ref_inds <- which(data_ind_tab$typ == "STD" & data_ind_tab$mix == 1)
  metab_inds <- which(data_ind_tab$typ == "STD" & data_ind_tab$mix == 0)
  
  ref_N <- length(ref_inds) 
  ref_start <- (ref_inds[1] - 1) * N + 1
  ref_end   <- ref_inds[ref_N] * N
  
  metab_N <- length(metab_inds) 
  metab_start <- (metab_inds[1] - 1) * N + 1
  metab_end   <- metab_inds[metab_N] * N
  
  raw_vec <- readBin(data, what = "double", n = 2 * N * (fid_num), size = 4,
                     endian = "little")
  
  cplx_vec <- raw_vec[c(TRUE, FALSE)] - 1i * raw_vec[c(FALSE, TRUE)]
  
  ref_data <- cplx_vec[ref_start:ref_end]
  dim(ref_data) <- c(N, chans, ref_N/chans, 1, 1, 1, 1)
  ref_data <- aperm(ref_data, c(7,6,5,4,3,2,1))
  
  metab_data <- cplx_vec[metab_start:metab_end]
  dim(metab_data) <- c(N, chans, metab_N/chans, 1, 1, 1, 1)
  metab_data <- aperm(metab_data, c(7,6,5,4,3,2,1))
  
  #res <- rep(NA, 7)
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)
  
  metab_mrs <- list(ft = ft, data = metab_data, resolution = res, te = NA,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  class(metab_mrs) <- "mrs_data"
  
  
  ref_mrs <- list(ft = ft, data = ref_data, resolution = res, te = NA,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  class(ref_mrs) <- "mrs_data"
  
  list(metab_mrs = metab_mrs, ref_mrs = ref_mrs)
}

#' @export
read_rda <- function(fname) {
  con = file(fname, "r")
  n = 1
  while (TRUE) {
    line = readLines(con, n = 1)
    if (startsWith(line, ">>> End of header <<<")) {
      data_pos <- seek(con)
      break
    }
    n = n + 1
  }
  
  # go back to the start
  seek(con, 0)
  txt <- read.delim(con, sep = ":", nrows = (n - 2), header = FALSE,
                    strip.white = TRUE, stringsAsFactors = FALSE,
                    comment.char = ">")
  close(con)
  
  # open in binary mode
  con <- file(fname, "rb")
  # skip the text bit
  seek(con, data_pos, "start", rw = "rb")
  raw_vec <- readBin(con, what = "double", n = 1024*2, size = 8, endian = "little")
  close(con)
  
}
