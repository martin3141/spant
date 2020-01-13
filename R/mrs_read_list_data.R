read_list_data <- function(fname, ft, fs, ref) {
  # generate matching data and list files
  ext <- stringr::str_sub(fname, -5)
  name <- stringr::str_sub(fname, 1, -6)
  
  if ( ext == ".list" ) {
    list <- fname
    data <- paste0(name, ".data")
  } else if ( ext == ".data" ) {
    data <- fname
    list <- paste0(name, ".list")
  } else {
    stop("Incorrect file extension.")
  }
  
  # check both files exist
  if (!file.exists(list)) {
    cat(list)
    stop("list file not found.")
  } else if (!file.exists(data)) {
    cat(data)
    stop("data file not found.")
  }
  
  # read list file as text
  txt <- as.array(readLines(list))
  
  N_txt <- ".    0    0    0  F-resolution"
  N_ind <- which(apply(txt, 1, startsWith, N_txt))
  N <- as.numeric(strsplit(txt[N_ind], ":")[[1]][2])
  
  data_ind_start_txt <- "# === START OF DATA VECTOR INDEX"
  data_ind_start <- which(apply(txt, 1, startsWith, data_ind_start_txt))
  data_ind_end_txt <- "# === END OF DATA VECTOR INDEX"
  data_ind_end <- which(apply(txt, 1, startsWith, data_ind_end_txt))
  data_ind_tab <- utils::read.table(text = txt[(data_ind_start + 3):(data_ind_end - 1)])
  col_names <- strsplit(txt[data_ind_start + 2], "\\s+")[[1]][2:22]
  colnames(data_ind_tab) <- col_names
  
  fid_num <- nrow(data_ind_tab)
  chans <- max(data_ind_tab$chan) + 1
  ref_inds <- which(data_ind_tab$typ == "STD" & data_ind_tab$mix == 1)
  metab_inds <- which(data_ind_tab$typ == "STD" & data_ind_tab$mix == 0)
  noise_inds <- which(data_ind_tab$typ == "NOI" & data_ind_tab$mix == 0)
  
  ref_N <- length(ref_inds) 
  ref_start <- (ref_inds[1] - 1) * N + 1
  ref_end   <- ref_inds[ref_N] * N
  
  metab_N <- length(metab_inds) 
  metab_start <- (metab_inds[1] - 1) * N + 1
  metab_end   <- metab_inds[metab_N] * N
  
  noise_N <- length(noise_inds) 
  noise_start <- (noise_inds[1] - 1) * N + 1
  noise_end   <- noise_inds[noise_N] * N
  
  raw_vec <- readBin(data, what = "double", n = 2 * N * (fid_num), size = 4,
                     endian = "little")
  
  res <- c(NA, NA, NA, NA, 1, NA, 1 / fs)
  
  # freq domain vector vector
  freq_domain <- rep(FALSE, 7)
  
  cplx_vec <- raw_vec[c(TRUE, FALSE)] - 1i * raw_vec[c(FALSE, TRUE)]
  
  if (is.na(ref_start)) {
    ref_mrs <- NA
  } else {
    ref_data <- cplx_vec[ref_start:ref_end]
    dim(ref_data) <- c(N, chans, ref_N/chans, 1, 1, 1, 1)
    ref_data <- aperm(ref_data, c(7,6,5,4,3,2,1))
    
    ref_mrs <- list(ft = ft, data = ref_data, resolution = res, te = NA,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
    class(ref_mrs) <- "mrs_data"
  }
  
  metab_data <- cplx_vec[metab_start:metab_end]
  dim(metab_data) <- c(N, chans, metab_N/chans, 1, 1, 1, 1)
  metab_data <- aperm(metab_data, c(7,6,5,4,3,2,1))
  
  noise_data <- cplx_vec[noise_start:noise_end]
  dim(noise_data) <- c(N, chans, noise_N/chans, 1, 1, 1, 1)
  noise_data <- aperm(noise_data, c(7,6,5,4,3,2,1))
  
  metab_mrs <- list(ft = ft, data = metab_data, resolution = res, te = NA,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  class(metab_mrs) <- "mrs_data"
  
  noise_mrs <- list(ft = ft, data = noise_data, resolution = res, te = NA,
                   ref = ref, row_vec = NA, col_vec = NA,
                   pos_vec = NA, freq_domain = freq_domain)
  class(noise_mrs) <- "mrs_data"
  
  list(metab = metab_mrs, ref = ref_mrs, noise = noise_mrs)
}