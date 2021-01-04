read_ima <- function(fraw, verbose = FALSE, extra) {
  
  tags <- list(ascii_hdr = "0029,1120", spec_data = "7FE1,1010")
  res  <- dicom_reader(fraw, tags)
  
  vars <- read_siemens_txt_hdr(res$ascii_hdr, "vd")
  
  # calculate expected size of data in bytes - assuming complex 4byte floats
  data_size <- vars$x_pts * vars$y_pts * vars$z_pts * vars$N * 4 * 2
  
  raw_pts <- readBin(res$spec_data, what = "double", n = data_size, size = 4L)
  
  # make complex
  data <- raw_pts[c(TRUE, FALSE)] + 1i * raw_pts[c(FALSE, TRUE)]
  
  data <- array(data, dim = c(vars$N, 1, 1, vars$z_pts, vars$y_pts, vars$x_pts, 
                              1))
  
  data <- aperm(data, c(7,5,6,4,3,2,1))
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  # get the resolution and geom info
  paras <- calc_siemens_paras(vars, TRUE)
  
  meta = list(EchoTime = vars$te)
  
  mrs_data <- mrs_data(data = data, ft = vars$ft, resolution = paras$res,
                       ref = paras$ref, nuc = paras$nuc,
                       freq_domain = freq_domain, affine = paras$affine,
                       meta = meta, extra = extra)
  
  return(mrs_data)
}

#' Read a directory containing Siemens MRS IMA files and combine along the coil
#' dimension. Note that the coil ID is inferred from the sorted file name and
#' should be checked when consistency is required between two directories.
#' @param dir data directory path.
#' @param extra an optional data frame to provide additional variables for use
#' in subsequent analysis steps, eg id or grouping variables.
#' @return mrs_data object.
#' @export
read_ima_coil_dir <- function(dir, extra = NULL) {
  files <- list.files(dir, full.names = TRUE)
  #warning("coil ordering is based on file name only.")
  files <- sort(files)
  mrs_list <- lapply(files, read_mrs, format = "ima", verbose = TRUE,
                     extra = extra)
  mrs_data <- append_coils(mrs_list)
  return(mrs_data)
}

#' Read a directory containing Siemens MRS IMA files and combine along the
#' dynamic dimension. Note that the coil ID is inferred from the sorted file
#' name and should be checked when consistency is required.
#' @param dir data directory path.
#' @param extra an optional data frame to provide additional variables for use
#' in subsequent analysis steps, eg id or grouping variables.
#' @return mrs_data object.
#' @export
read_ima_dyn_dir <- function(dir, extra = NULL) {
  files <- list.files(dir, full.names = TRUE)
  #warning("coil ordering is based on file name only.")
  files <- sort(files)
  mrs_list <- lapply(files, read_mrs, format = "ima", verbose = TRUE,
                     extra = extra)
  mrs_data <- append_dyns(mrs_list)
  return(mrs_data)
}