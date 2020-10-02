read_lcm_raw <- function(fname, ft, fs, ref) {
  in_nmid <- FALSE 
  con <- file(fname, "rb")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (endsWith(line, "$NMID")) {
      in_nmid <- TRUE
    } else if (endsWith(line, "$END") && in_nmid) {
      x <- utils::read.fortran(con, "2F15.0")
      break
    }
  }
  close(con)
  
  N <- nrow(x)
  data <- x$V1 + x$V2 * 1i
  dim(data) <- c(1, 1, 1, 1, 1, 1, N)
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  te <- NA
  nuc <- def_nuc()
  
  mrs_data <- list(ft = ft, data = data, resolution = res,
                   te = NA, ref = ref, nuc = nuc, row_vec = NA, col_vec = NA,
                   sli_vec = NA, pos_vec = NA, freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}