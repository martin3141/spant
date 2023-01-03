read_lcm_raw <- function(fname, ft, fs, ref, extra) {
  in_nmid <- FALSE 
  con <- file(fname, "r")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (grepl("$NMID", line, fixed = TRUE)) {       # in the header
      in_nmid <- TRUE
    } else if (endsWith(line, "$END") && in_nmid) { # in the numerical part
      x <- utils::read.table(con, fill = TRUE)
      break
    }
  }
  close(con)
  
  data <- as.vector(t(as.matrix(x)))
  data <- stats::na.omit(data)
  N <- length(data) / 2
  data <-      data[seq(1, 2 * N, 2)] +
          1i * data[seq(2, 2 * N, 2)]
  
  dim(data) <- c(1, 1, 1, 1, 1, 1, N)
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  nuc <- def_nuc()
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = freq_domain, affine = NULL,
                       meta = NULL, extra = extra)
  
  return(mrs_data)
}