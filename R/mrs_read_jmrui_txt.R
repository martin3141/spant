read_mrs_jmrui_txt <- function(fname, extra) {
 
  # read the full file into memory
  txt_lines <- readLines(fname)
  
  # find the start of the data
  n_data_start <- which(txt_lines == "Signal and FFT") + 3
  n            <- length(txt_lines)
  
  # read some basic parameters
  header <- utils::read.delim(textConnection(txt_lines[1:(n_data_start - 4)]),
                              sep = ":")
  
  header <- as.data.frame(t(header))
  
  fs <- 1 / as.numeric(header$SamplingInterval) * 1e3
  ft <- as.numeric(header$TransmitterFrequency)
  
  mrs_data_table <- utils::read.delim(textConnection(txt_lines[n_data_start:n]),
                                      header = FALSE)
  
  data <- mrs_data_table$V1 - 1i * mrs_data_table$V2
  
  N <- length(data)
  
  dim(data) <- c(1, 1, 1, 1, 1, 1, N)
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)

  nuc <- def_nuc()
  
  ref <- def_ref()
  
  mrs_data <- mrs_data(data = data, ft = ft, resolution = res, ref = ref,
                       nuc = nuc, freq_domain = freq_domain, affine = NULL,
                       meta = NULL, extra = extra)
  
  return(mrs_data)
}
