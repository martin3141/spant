read_ima <- function(fname) {
  vars <- read_siemens_txt_hdr(fname, "vd")
  
  # calculate expected size of data in bytes - assuming complex 4byte floats
  data_size <- vars$N * 4 * 2
  
  con <- file(fname, 'rb')
  # assume data points are at the end of the file 
  seek(con, -data_size, origin = "end")
  raw_pts <- readBin(con, "numeric", size = 4L, n = (vars$N * 2), endian = "little")
  close(con)
  
  # make complex
  data <- raw_pts[c(TRUE, FALSE)] + 1i * raw_pts[c(FALSE, TRUE)]
  
  data <- array(data, dim = c(vars$N, 1, 1, 1, 1, 1, 1))
  data <- aperm(data, c(7,6,5,4,3,2,1))
   
  res <- c(NA, NA, NA, NA, 1, NA, 1 / vars$fs * 2)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  ref <- def_acq_paras()$ref
  
  mrs_data <- list(ft = vars$ft, data = data, resolution = res,
                   te = vars$te, ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}