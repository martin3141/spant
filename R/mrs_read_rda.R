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