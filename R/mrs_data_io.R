#' Read MRS data from a file.
#' @param fname The filename of the dpt format MRS data.
#' @param format A string describing the data format. May be one of the 
#' following : "spar_sdat", "rda", "pfile", "list_data", "paravis", "dpt".
#' @param ft Transmitter frequency in Hz (required for list_data format).
#' @param fs Sampling frequency in Hz (required for list_data format).
#' @param ref Reference value for ppm scale (required for list_data format).
#' @return An MRS data object.
#' @examples
#' fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
#' mrs_data <- read_mrs(fname, format = "spar_sdat")
#' print(mrs_data)
#' @export
read_mrs <- function(fname, format, ft = NULL, fs = NULL, ref = NULL) {
  if (format == "spar_sdat") {
    return(read_spar_sdat(fname))
  } else if (format == "rda") {
    return(read_rda(fname))
  } else if (format == "pfile") {
    return(read_pfile(fname))
  } else if (format == "list_data") {
    if (is.null(ft)) stop("Please specify ft parameter for list_data format")
    if (is.null(fs)) stop("Please specify fs parameter for list_data format")
    if (is.null(ref)) stop("Please specify ref parameter for list_data format")
    return(read_list_data(fname, ft, fs, ref))
  } else if (format == "dpt") {
    return(read_mrs_dpt(fname))
  } else if (format == "paravis") {
    return(read_paravis_raw(fname))
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
  data = aperm(data,c(6, 2, 3, 4, 5, 7, 1))
  
  res <- c(NA, row_dim, col_dim, slice_dim, 1, NA, 1 / fs)
  ref <- def_acq_paras()$ref
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  mrs_data <- list(ft = ft, data = data, resolution = res, te = te, ref = ref, 
                   row_vec = row_ori, col_vec = col_ori, pos_vec = pos_vec, 
                   freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}

read_list_data <- function(fname, ft, fs, ref) {
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
  data_ind_tab <- utils::read.table(text = txt[(data_ind_start + 3):(data_ind_end - 1)])
  col_names <- strsplit(txt[data_ind_start + 2], "\\s+")[[1]][2:22]
  colnames(data_ind_tab) <- col_names
  
  fid_num <- nrow(data_ind_tab)
  chans <- max(data_ind_tab$chan) + 1
  ref_inds <- which(data_ind_tab$typ == "STD" & data_ind_tab$mix == 1)
  metab_inds <- which(data_ind_tab$typ == "STD" & data_ind_tab$mix == 0)
  noise_inds <- which(data_ind_tab$typ == "NOI" & data_ind_tab$mix == 0)
  
  ref_N <- length(ref_inds) 
  ref_start <- (ref_inds[1] - 1) * N + 1
  ref_end   <- ref_inds[ref_N] * N
  
  metab_N <- length(metab_inds) 
  metab_start <- (metab_inds[1] - 1) * N + 1
  metab_end   <- metab_inds[metab_N] * N
  
  noise_N <- length(noise_inds) 
  noise_start <- (noise_inds[1] - 1) * N + 1
  noise_end   <- noise_inds[noise_N] * N
  
  raw_vec <- readBin(data, what = "double", n = 2 * N * (fid_num), size = 4,
                     endian = "little")
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)
  
  cplx_vec <- raw_vec[c(TRUE, FALSE)] - 1i * raw_vec[c(FALSE, TRUE)]
  
  if (is.na(ref_start)) {
    ref_mrs <- NA
  } else {
    ref_data <- cplx_vec[ref_start:ref_end]
    dim(ref_data) <- c(N, chans, ref_N/chans, 1, 1, 1, 1)
    ref_data <- aperm(ref_data, c(7,6,5,4,3,2,1))
    
    ref_mrs <- list(ft = ft, data = ref_data, resolution = res, te = NA,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
    class(ref_mrs) <- "mrs_data"
  }
  
  metab_data <- cplx_vec[metab_start:metab_end]
  dim(metab_data) <- c(N, chans, metab_N/chans, 1, 1, 1, 1)
  metab_data <- aperm(metab_data, c(7,6,5,4,3,2,1))
  
  noise_data <- cplx_vec[noise_start:noise_end]
  dim(noise_data) <- c(N, chans, noise_N/chans, 1, 1, 1, 1)
  noise_data <- aperm(noise_data, c(7,6,5,4,3,2,1))
  
  metab_mrs <- list(ft = ft, data = metab_data, resolution = res, te = NA,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  class(metab_mrs) <- "mrs_data"
  
  
  
  
  noise_mrs <- list(ft = ft, data = noise_data, resolution = res, te = NA,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  class(noise_mrs) <- "mrs_data"
  
  list(metab = metab_mrs, ref = ref_mrs, noise = noise_mrs)
}

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
  txt <- utils::read.delim(con, sep = ":", nrows = (n - 2), header = FALSE,
                    strip.white = TRUE, stringsAsFactors = FALSE,
                    comment.char = ">")
  close(con)
  
  N <- as.integer(txt$V2[which(txt$V1 == "VectorSize")])
  fs <- 1e6 / as.numeric(txt$V2[which(txt$V1 == "DwellTime")])
  ft <- 1e6 * as.numeric(txt$V2[which(txt$V1 == "MRFrequency")])
  te <- as.numeric(txt$V2[which(txt$V1 == "TE")]) / 1e3
  #avgs <- as.integer(txt$V2[which(txt$V1 == "NumberOfAverages")]) 
  rows <- as.numeric(txt$V2[which(txt$V1 == "CSIMatrixSize[0]")])
  cols <- as.numeric(txt$V2[which(txt$V1 == "CSIMatrixSize[1]")])
  slices <- as.numeric(txt$V2[which(txt$V1 == "CSIMatrixSize[2]")])
  
  fids <- rows * cols * slices
  
  # open in binary mode
  con <- file(fname, "rb")
  # skip the text bit
  seek(con, data_pos, "start", rw = "rb")
  raw_vec <- readBin(con, what = "double", n = N * 2 * fids, size = 8,
                     endian = "little")
  close(con)
  
  data <- raw_vec[c(TRUE, FALSE)] + 1i * raw_vec[c(FALSE, TRUE)]
  
  dim(data) <- c(N, rows, cols, slices, 1, 1, 1)
  data <- aperm(data, c(7,2,3,4,5,6,1))
  
  res <- c(NA, rows, cols, slices, 1, NA, 1 / fs)
  ref <- def_acq_paras()$ref
  row_ori = NA
  col_ori = NA
  pos_vec = NA
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  mrs_data <- list(ft = ft, data = data, resolution = res, te = te, ref = ref, 
                   row_vec = row_ori, col_vec = col_ori, pos_vec = pos_vec, 
                   freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}

read_paravis_raw <- function(fname) {
  # find the method file in the same directory
  method_fname <- file.path(dirname(fname), "method")
  
  if (!file.exists(method_fname)) {
    cat(method_fname)
    stop("method file not found.")
  }
  
  # read paramters
  lines <- utils::read.delim(method_fname, sep = "=", header = FALSE, 
                      stringsAsFactors = FALSE)
  
  reps <- as.integer(get_para_val(lines, "##$PVM_NRepetitions"))
  avgs <- as.integer(get_para_val(lines, "##$PVM_NAverages"))
  dynamics <- reps * avgs
  N <- as.integer(get_para_val(lines, "##$PVM_DigNp"))
  fs <- as.double(get_para_val(lines, "##$PVM_DigSw"))
  shift <- as.integer(get_para_val(lines, "##$PVM_DigShift"))
  coils <- as.integer(get_para_val(lines, "##$PVM_EncNReceivers"))
  ft_str <- lines$V1[1 + which(lines$V1 == "##$PVM_FrqRef")]
  ft <- as.double(strsplit(ft_str, " ")[[1]][1]) * 1e6
  te <- as.double(get_para_val(lines, "##$PVM_EchoTime")) / 1e3
  
  expected_Npts <- dynamics * N * 2 * coils
    
  # read the raw data file 
  fbytes <- file.size(fname)
  Npts <- fbytes / 4
  
  if (Npts != expected_Npts) warning("Unexpected number of data points.")
  
  raw_vec <- readBin(fname, "int", size = 4, n = Npts)
  data <- raw_vec[c(TRUE, FALSE)] - 1i * raw_vec[c(FALSE, TRUE)]
  
  dim(data) <- c(N, coils, dynamics, 1, 1, 1, 1)
  data <- aperm(data, c(7,6,5,4,3,2,1))
  
  # move dig. filter guff to end of the FID (the Bruker way of doing things?)
  #data <- abind::abind(data[,,,,,,(shift + 1):N,drop = FALSE], 
  #                     data[,,,,,,1:shift,drop = FALSE], along = 7)
  
  filt_pts <- data[,,,,,,shift:1,drop = FALSE]
  second_part <- data[,,,,,,(shift + 1):N, drop = FALSE]
  data <- second_part 
  data[,,,,,,1:shift] <- data[,,,,,,1:shift, drop = FALSE] - filt_pts
  data <- abind::abind(data[,,,,,,,drop = FALSE], 
                       array(0, dim = dim(filt_pts)), along = 7)
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  ref <- def_acq_paras()$ref
  
  mrs_data <- list(ft = ft, data = data, resolution = res, te = te,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}

get_para_val <- function(lines, name_str) {
  lines$V2[which(lines$V1 == name_str)]
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
  vars$ps_mps_freq <- sum(2^.subset(0:32, as.logical(ps_mps_freq_bits)))
  
  seek(con, loc$te)
  vars$te <- readBin(con, "int", size = 4, endian = endian)
  
  close(con)
  
  vars
}