read_paravis_raw <- function(fname) {
  # find the method file in the same directory
  method_fname <- file.path(dirname(fname), "method")
  
  if (!file.exists(method_fname)) {
    cat(method_fname)
    stop("method file not found.")
  }
  
  # read parameters
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