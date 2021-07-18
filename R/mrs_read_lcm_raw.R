read_lcm_raw <- function(fname, ft, fs, ref, extra) {
  in_nmid <- FALSE 
  con <- file(fname, "rb")
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (endsWith(line, "$NMID")) {
      in_nmid <- TRUE
    } else if (endsWith(line, "$END") && in_nmid) {
      
      # ARC 2021-07-15: auto-detect column layout {{{
      # read a single line, then rewind
      fp <- seek(con, origin='current')
      l1 <- readLines(con, n=1, warn=FALSE);
      fpn <- seek(con, origin='start', where=fp);
  
      tokens <- strsplit(trimws(l1),"[[:space:]]+")[[1]];
      cols <- length(tokens);
      width <- ceiling(nchar(l1)/cols);
      fmt <- sprintf("%dF%d.0", cols, width);
      # }}}
          
      x <- utils::read.fortran(con, fmt);
      break
    }
  }
  close(con)

  data <- as.vector(t(as.matrix(x)))
  N <- length(data)/2
  data <- data[seq(1, 2 * N, 2)] +
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
