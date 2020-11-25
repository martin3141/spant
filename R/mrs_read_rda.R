read_rda <- function(fname) {
  
  # find out where the raw data points start by looking for
  # ">>> End of header <<<"
  con <- file(fname, "rb")
  n <- 1
  while (TRUE) {
    line <- scan(con, "character", nlines = 1, sep = "\n", quiet = TRUE)
    if (startsWith(line, ">>> End of header <<<")) {
      data_pos <- seek(con)
      break
    }
    n = n + 1
  }
  close(con)
  
  # go back to the start and read in the data parameters
  con <- file(fname, "r")
  txt <- utils::read.delim(con, sep = ":", nrows = (n - 2), header = FALSE,
                    strip.white = TRUE, stringsAsFactors = FALSE,
                    comment.char = ">")
  close(con)
  
  row_ori <- rep(NA, 3)
  col_ori <- rep(NA, 3)
  pos_vec <- rep(NA, 3)
  
  N <- as.integer(txt$V2[which(txt$V1 == "VectorSize")])
  fs <- 1e6 / as.numeric(txt$V2[which(txt$V1 == "DwellTime")])
  ft <- 1e6 * as.numeric(txt$V2[which(txt$V1 == "MRFrequency")])
  te <- as.numeric(txt$V2[which(txt$V1 == "TE")]) / 1e3
  rows <- as.numeric(txt$V2[which(txt$V1 == "CSIMatrixSize[0]")])
  cols <- as.numeric(txt$V2[which(txt$V1 == "CSIMatrixSize[1]")])
  slices <- as.numeric(txt$V2[which(txt$V1 == "CSIMatrixSize[2]")])
  
  pos_vec[1] <- as.numeric(txt$V2[which(txt$V1 == "PositionVector[0]")])
  pos_vec[2] <- as.numeric(txt$V2[which(txt$V1 == "PositionVector[1]")])
  pos_vec[3] <- as.numeric(txt$V2[which(txt$V1 == "PositionVector[2]")])
  
  row_ori[1] <- as.numeric(txt$V2[which(txt$V1 == "RowVector[0]")])
  row_ori[2] <- as.numeric(txt$V2[which(txt$V1 == "RowVector[1]")])
  row_ori[3] <- as.numeric(txt$V2[which(txt$V1 == "RowVector[2]")])
  
  col_ori[1] <- as.numeric(txt$V2[which(txt$V1 == "ColumnVector[0]")])
  col_ori[2] <- as.numeric(txt$V2[which(txt$V1 == "ColumnVector[1]")])
  col_ori[3] <- as.numeric(txt$V2[which(txt$V1 == "ColumnVector[2]")])
  
  col_vox_dim <- as.numeric(txt$V2[which(txt$V1 == "PixelSpacingCol")])
  row_vox_dim <- as.numeric(txt$V2[which(txt$V1 == "PixelSpacingRow")])
  slice_vox_dim <- as.numeric(txt$V2[which(txt$V1 == "PixelSpacing3D")])
  
  pos_vec <- pos_vec + row_ori * row_vox_dim / 2 + col_ori * col_vox_dim / 2
  sli_vec <- crossprod_3d(col_ori, row_ori)
  
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
  data <- aperm(data, c(7, 2, 3, 4, 5, 6, 1))
  
  res <- c(NA, row_vox_dim, col_vox_dim, slice_vox_dim, 1, NA, 1 / fs)
  ref <- def_ref()
  
  # TODO determine from the data
  nuc <- def_nuc()
  
  # freq domain vector
  freq_domain <- rep(FALSE, 7)
  
  mrs_data <- list(ft = ft, data = data, resolution = res, te = te, ref = ref, 
                   nuc = nuc, row_vec = row_ori, col_vec = col_ori,
                   sli_vec = sli_vec, pos_vec = pos_vec,
                   freq_domain = freq_domain)
  
  class(mrs_data) <- "mrs_data"
  mrs_data
}