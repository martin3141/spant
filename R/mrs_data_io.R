#' Read MRS data from a file.
#' @param fname filename of the dpt format MRS data.
#' @param format string describing the data format. May be one of the 
#' following : "spar_sdat", "rda", "ima", "twix", "pfile", "list_data",
#' "paravis", "dpt".
#' @param ft transmitter frequency in Hz (required for list_data format).
#' @param fs sampling frequency in Hz (required for list_data format).
#' @param ref reference value for ppm scale (required for list_data format).
#' @param n_ref_scans override the number of water reference scans detected in
#' the file header (GE p-file only).
#' @return MRS data object.
#' @examples
#' fname <- system.file("extdata", "philips_spar_sdat_WS.SDAT", package = "spant")
#' mrs_data <- read_mrs(fname, format = "spar_sdat")
#' print(mrs_data)
#' @export
read_mrs <- function(fname, format, ft = NULL, fs = NULL, ref = NULL,
                     n_ref_scans = NULL) {
  
  if (format == "spar_sdat") {
    return(read_spar_sdat(fname))
  } else if (format == "rda") {
    return(read_rda(fname))
  } else if (format == "ima") {
    return(read_ima(fname))
  } else if (format == "twix") {
    return(read_twix(fname))
  } else if (format == "pfile") {
    return(read_pfile(fname, n_ref_scans))
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
#' @param fname filename of the dpt format MRS data.
#' @return MRS data object.
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
#' @param fname the filename containing the MRS data.
#' @param fname_ref a second filename containing reference MRS data.
#' @param format format of the MRS data. Can be one of the following:
#' siemens, philips, ge, dcm, dpt, rda, lcm, varian, bruker, jmrui_txt.
#' @param id optional ID string.
#' @param group optional group string.
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
    #main$data <- comb_metab_ref(main, ref)
    #main$data <- abind::abind(main$data, ref$data, along=1)
    return(list(metab = main, ref = ref))
  } else {
    return(main)
  }
}

#' Write MRS data object to file in dangerplot (dpt) v2 format.
#' @param fname the filename of the output dpt format MRS data.
#' @param mrs_data object to be written to file.
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

#' Write MRS data object to file in a RAW format compatible with LCModel.
#' @param fname the filename of the output RAW format MRS data.
#' @param mrs_data object to be written to file.
#' @param id text string to identify the data.
#' @examples
#' \dontrun{
#' mrs_data <- write_mrs_lcm_raw("my_mrs_data.RAW", my_mrs_data)
#' }
#' @export
write_mrs_lcm_raw <- function(fname, mrs_data, id = NA) {
  sig <- mrs_data$data[1, 1, 1, 1, 1, 1,]
  N <- length(sig)
  sink(fname)
  cat(" $NMID\n")
  if (is.na(id)) id <- fname
  cat(paste(" ID='", id, "', FMTDAT='(2E15.6)'\n", sep = ""))
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