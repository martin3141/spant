read_ima <- function(fname) {
  vars <- read_siemens_txt_hdr(fname, "vd")
  
  # calculate expected size of data in bytes - assuming complex 4byte floats
  data_size <- vars$x_pts * vars$y_pts * vars$z_pts * vars$N * 4 * 2
  
  dcm_info <- oro.dicom::readDICOM(fname, pixelData = FALSE)$hdr
  dcm_info <- as.data.frame(dcm_info)
  padding <- 0
  if (dcm_info[nrow(dcm_info),3] == "DataSetTrailingPadding") {
    padding <- as.numeric(dcm_info[nrow(dcm_info),5]) + 12
  }
  
  con <- file(fname, 'rb')
  # assume data points are at the end of the file 
  seek(con, -data_size - padding, origin = "end")
  #print(seek(con, where = NA))
  raw_pts <- readBin(con, "numeric", size = 4L, n = (vars$x_pts * vars$y_pts * 
                                                     vars$z_pts * vars$N * 2),
                     endian = "little")
  close(con)
  
  # make complex
  data <- raw_pts[c(TRUE, FALSE)] + 1i * raw_pts[c(FALSE, TRUE)]
  
  data <- array(data, dim = c(vars$N, 1, 1, vars$z_pts, vars$y_pts, vars$x_pts, 
                              1))
  
  data <- aperm(data, c(7,6,5,4,3,2,1))
   
  res <- c(NA, vars$x_dim / vars$x_pts, vars$y_dim / vars$y_pts,
           vars$z_dim / vars$z_pts, 1, NA, 1 / vars$fs * 2)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  ref <- def_acq_paras()$ref
  
  mrs_data <- list(ft = vars$ft, data = data, resolution = res,
                   te = vars$te, ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}