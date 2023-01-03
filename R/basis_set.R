#' Generate a basis file using TARQUIN.
#' @param basis_file filename of the basis file to be generated.
#' @param metab_data MRS data object to match the generated basis parameters.
#' @param opts list of options to pass to TARQUIN.
#' @examples
#' \dontrun{
#' write_basis_tqn('test.basis',mrs_data,c("--echo","0.04"))
#' }
#' @export
write_basis_tqn <- function(basis_file, metab_data, opts = NULL) {
  # write the metab data to a temp dpt file
  metab_file <- tempfile()
  write_mrs_dpt_v2(metab_file, metab_data)
  
  # run TARQUIN
  cmd = paste(getOption("spant.tqn_cmd"), "--input",metab_file, "--format dpt",
              "--no_fit true", "--water_width 0", "--auto_ref false",
              "--auto_phase false", "--output_basis_lcm", basis_file,
              paste(opts,collapse = " "))
  
  res = system(cmd, intern = TRUE)
  
  if (!file.exists(basis_file)) {
    print(res)
    print(cmd)
    stop("Error loading data with above TARQUIN command.")
  }
}

#' Simulate a basis file using TARQUIN.
#' @param fs sampling frequency
#' @param ft transmitter frequency
#' @param N number of data points
#' @param ref chemical shift reference
#' @param opts list of options to pass to TARQUIN.
#' @examples
#' \dontrun{
#' write_basis_tqn('test.basis',mrs_data,c("--echo","0.04"))
#' }
#' @export
sim_basis_tqn <- function(fs = def_fs(), ft = def_ft(), N = def_N(),
                          ref = def_ref(), opts = NULL) {
  
  metab_data <- sim_zero(fs = fs, ft = ft, N = N, ref = ref)
  
  # write the metab data to a temp dpt file
  metab_file <- tempfile()
  write_mrs_dpt_v2(metab_file, metab_data)
  
  # specify a basis file 
  basis_file <- tempfile()
  
  # run TARQUIN
  cmd = paste(getOption("spant.tqn_cmd"), "--input", metab_file, "--format dpt",
              "--no_fit true", "--water_width 0", "--auto_ref false",
              "--auto_phase false", "--output_basis_lcm", basis_file,
              paste(opts,collapse = " "))
  
  res = system(cmd, intern = TRUE)
  
  if (!file.exists(basis_file)) {
    print(res)
    print(cmd)
    stop("Error loading data with above TARQUIN command.")
  }
  
  read_basis(basis_file, ref = ref)
}

#' @export
print.basis_set <- function(x, ...) {
  cat("Basis set parameters\n")
  cat("-------------------------------\n")
  cat(paste(c("Trans. freq (MHz)       : ", round(x$ft * 1e-6, 1), "\n")), sep = "")
  cat(paste(c("Data points             : ", dim(x$data)[1], "\n")), sep = "")
  cat(paste(c("Sampling frequency (Hz) : ", x$fs, "\n")), sep = "")
  cat(paste(c("Elements                : ", dim(x$data)[2], "\n\n")), sep = "")
  cat(paste(c("Names\n")), sep = "")
  cat("-------------------------------\n")
  cat(x$names, sep = ",", fill = 31)
}

#' @export
plot.basis_set <- function(x, n = 1, ...) {
  mrs_basis <- basis2mrs_data(x)
  graphics::plot(mrs_basis, dyn = n, ...)
  #matplot(Re(x$data),type="l",lty=1,col=1)
}

#' @export
stackplot.basis_set <- function(x, ...) {
  mrs_basis <- basis2mrs_data(x)
  stackplot(mrs_basis, ...)
}

#' Read a basis file in LCModel .basis format.
#' @param basis_file path to basis file.
#' @param ref assumed ppm reference value.
#' @param sort_basis sort the basis set based on signal names.
#' @return basis object.
#' @export
read_basis <- function(basis_file, ref = def_ref(), sort_basis = TRUE) {
  
  # open the file
  con  <- file(basis_file, open = "r")
  names <- vector()
  data <- vector()
  
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    # reading main header
    if (startsWith(line, " NDATAB = ")) {
      N <- as.integer(strsplit(trimws(line), "\\s+")[[1]][3])
    } else if (startsWith(line, " HZPPPM = ")) {
      bas_ft <- strsplit(trimws(line), "\\s+")[[1]][3]
      bas_ft <- as.double(gsub(",", "", bas_ft)) * 1e6
    } else if (startsWith(line, " BADELT = ")) {
      bas_fs <- strsplit(trimws(line), "\\s+")[[1]][3]
      bas_fs <- 1 / as.double(gsub(",", "", bas_fs))
    } else if (startsWith(line, " FMTBAS = ")) {
      format <- (strsplit(trimws(line), " ")[[1]][3])
      format <- gsub(",", "", format)
      format <- gsub("'", "", format)
      format <- gsub("\\(", "", format)
      format <- gsub("\\)", "", format)
    } else if (endsWith(line, "$BASIS")) { # we're at the end of the main header
      cols <- as.numeric(substr(format, 1, 1))
      data_lines <- ceiling(2 * N / cols)
      while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
        if (startsWith(line, " ID = ")) {
          # reading metabolite header
          id <- (strsplit(trimws(line), " ")[[1]][3])
          id <- gsub(",", "", id)
          id <- gsub("'", "", id)
          names <- c(names, id)
        } else if (endsWith(line, "$END")) {
          # read data points
          x <- utils::read.table(con, nrows = data_lines, fill = TRUE)
          data_pts <- as.vector(t(as.matrix(x)))
          data_pts <- stats::na.omit(data_pts)
          data_pts <-      data_pts[seq(1, 2 * N, 2)] +
                      1i * data_pts[seq(2, 2 * N, 2)]
          
          # ifftshift
          data_pts <- pracma::ifftshift(data_pts)
          data     <- cbind(data, data_pts) # slow?
          break
        }
      }
    }
  }
  close(con)
  
  basis_set <- list(data = data, N = N, fs = bas_fs, ft = bas_ft, 
                    names = names, ref = ref)
  
  class(basis_set) <- "basis_set"
  
  if (sort_basis) basis_set <- sort_basis(basis_set)
  
  dimnames(basis_set$data) <- NULL
  
  return(basis_set)
}

#' Write a basis object to an LCModel .basis formatted file.
#' @param basis basis object to be exported.
#' @param basis_file path to basis file to be generated.
#' @param fwhmba parameter used by LCModel.
#' @export
write_basis <- function(basis, basis_file, fwhmba = 0.1) {
  mrs_data <- basis2mrs_data(basis)
  N <- Npts(mrs_data)
  
  sink(basis_file)
  cat(" $SEQPAR\n", sep = "")
  cat(" FWHMBA =  ", fwhmba, " ,\n", sep = "")
  cat(" HZPPPM =  ",mrs_data$ft * 1e-6,",\n", sep = "")
  cat(" ECHOT =  ",mrs_data$meta$EchoTime * 1000,",\n", sep = "")
  cat(" SEQ = 'PRESS' $END\n", sep = "")
  cat(" $BASIS1\n", sep = "")
  cat(" IDBASI = 'SPANT',\n", sep = "")
  cat(" FMTBAS = '(6E13.5)',\n", sep = "")
  cat(" BADELT =  ", 1 / fs(mrs_data), ",\n", sep = "")
  cat(" NDATAB = ", Npts(mrs_data), " $END\n",sep = "")
  
  for (n in 1:Ndyns(mrs_data)) {
    name = basis$names[n]
    cat(" $NMUSED\n")
    cat(" FILERAW = '", name,"',\n", sep = "")
    cat(" METABO_CONTAM = ' ',\n")
    cat(" METABO_SINGLET = ' ',\n")
    cat(" AUTOPH = F,\n")
    cat(" AUTOSC = F,\n")
    cat(" DO_CONTAM = F,\n")
    cat(" FLATEN = F,\n")
    cat(" CONCSC = -1.,\n")
    cat(" DEGZER =  0.,\n")
    cat(" DEGPAP =  0.,\n")
    cat(" DEGPPM =  0.,\n")
    cat(" HWDPHA =  0.100000001,\n")
    cat(" HWDSCA =  0.100000001,\n")
    cat(" PPMAPP =  0.5 -0.5,\n")
    cat(" PPMAPP_CONTAM =  0.  0.,\n")
    cat(" PPMBAS =  0.00999999978,\n")
    cat(" PPMPK =  0.00999999978,\n")
    cat(" PPMPK_CONTAM =  0.,\n")
    cat(" PPMPHA =  0.,\n")
    cat(" PPMSCA =  8.43999958,\n")
    cat(" PPM_SPLIT = -999.,\n")
    cat(" RINTEG =  0.,\n")
    cat(" SDPNTS =  1.,\n")
    cat(" XTRASH =  0. $END\n")
    cat(" $BASIS\n")
    cat(" ID = '", name, "',\n", sep = "")
    cat(" METABO = '", name, "',\n", sep = "")
    cat(" CONC = 1.0,\n")
    cat(" TRAMP = 1.0,\n")
    cat(" VOLUME = 1.0,\n")
    cat(" ISHIFT = 0 $END\n")
    sig <- pracma::fftshift(mrs_data$data[1, 1, 1, 1, n, 1,])
    for (n in seq(1, N, 3)) {
      cat(" ")
      cat(noquote(formatC(c(Re(sig[n]), Im(sig[n])), width = 12,format = "E",
                            digits = 5)))
      
      if ((n + 1) < N) {
        cat(" ")
        cat(noquote(formatC(c(Re(sig[n + 1]), Im(sig[n + 1])), width = 12, 
                              format = "E", digits = 5)))
      }
      if ((n + 2) < N) {
        cat(" ")
        cat(noquote(formatC(c(Re(sig[n + 2]), Im(sig[n + 2])), width = 12, 
                              format = "E", digits = 5)))
      }
      cat("\n")
    }
  }
  sink()
}
  
#' Convert a basis object to an mrs_data object - where basis signals are spread
#' across the dynamic dimension.
#' @param basis basis set object.
#' @param sum_elements return the sum of basis elements (logical)
#' @param amps a vector of scaling factors to apply to each basis element.
#' @param shifts a vector of frequency shifts (in ppm) to apply to each basis
#' element.
#' @return an mrs_data object with basis signals spread across the dynamic 
#' dimension or summed.
#' @export
basis2mrs_data <- function(basis, sum_elements = FALSE, amps = NULL,
                           shifts = NULL) {
  
  res <- mat2mrs_data(t(basis$data), fs = basis$fs, ft = basis$ft,
                      ref = basis$ref, fd = TRUE)
  
  n_sigs <- Ndyns(res)
  
  # scale basis elements
  if (!is.null(amps)) {
    if (n_sigs != length(amps)) {
      stop(paste("Error, length of amps does not match the number of basis elements :", dim(basis$data)[2]))
    }
    for (n in 1:n_sigs) {
      res$data[,,,,n ,,] <- res$data[,,,,n ,,] * amps[n]
    }
  }
  
  # shift basis elements
  if (!is.null(shifts)) {
    if (n_sigs != length(shifts)) {
      stop(paste("Error, length of amps does not match the number of basis elements :", dim(basis$data)[2]))
    }
    res <- shift(res, shifts)
  }
  
  if (sum_elements) res <- sum_dyns(res)
  
  res
}

#' Convert an mrs_data object to basis object - where basis signals are spread
#' across the dynamic dimension in the MRS data.
#' @param mrs_data mrs_data object with basis signals spread across the dynamic dimension.
#' @param names list of names corresponding to basis signals.
#' @return basis set object.
#' @export
mrs_data2basis <- function(mrs_data, names) {
  # transform to FD
  if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
  }
  
  data <- matrix(mrs_data$data[1, 1, 1, 1,, 1,], 
            nrow = Npts(mrs_data), ncol = Ndyns(mrs_data), byrow = TRUE)
  
  basis_set <- list(data = data, N = Npts(mrs_data), fs = fs(mrs_data), 
                    ft = mrs_data$ft, names = names, ref = mrs_data$ref)
  
  class(basis_set) <- "basis_set"
  basis_set
}

#' Sort the basis-set elements alphabetically.
#' @param basis input basis.
#' @return sorted basis.
#' @export
sort_basis <- function(basis) {
  names_sorted <- sort(basis$names, index.return = TRUE)
  basis$names <- names_sorted$x
  basis$data <- basis$data[, names_sorted$ix]
  basis
}

#' Combine a pair of basis set objects.
#' @param basis_a first basis.
#' @param basis_b second basis.
#' @return combined basis set object.
#' @export
append_basis <- function(basis_a, basis_b) {
  basis_a_mrs <- basis2mrs_data(basis_a)
  basis_b_mrs <- basis2mrs_data(basis_b)
  basis_out_mrs <- append_dyns(basis_a_mrs, basis_b_mrs)
  basis_out <- mrs_data2basis(basis_out_mrs, c(basis_a$names, basis_b$names))
  return(basis_out)
}

#' Apply frequency shifts to basis set signals.
#' @param basis the basis to apply the shift to.
#' @param shifts a vector of frequency shifts to apply in ppm units. Must be the
#' same length as there are basis elements, or one value to be applied to all
#' elements.
#' @return modified basis set object.
#' @export
shift_basis <- function(basis, shifts) {
  if (length(shifts) != ncol(basis$data)) {
    if (length(shifts == 1)) {
      rep(shifts, ncol(basis$data))    
    } else {
      stop("the length of shifts does not match the number of basis signals")
    }
  }
  basis_mrs_data <- basis2mrs_data(basis)
  basis_mrs_data <- shift(basis_mrs_data, shifts)
  mrs_data2basis(basis_mrs_data, basis$names)
}

#' Resample a basis-set to match a mrs_data acquisition.
#' @param basis the basis to be resampled.
#' @param mrs_data the mrs_data to match the number of data points and sampling
#' frequency.
#' @param ref_freq_match apply a frequency shift to the basis to match the
#' reference frequency (usually 4.65 or 4.68) of the mrs_data.
#' @return resampled basis set object.
#' @export
resample_basis <- function(basis, mrs_data, ref_freq_match = TRUE) {
  N  <- basis$N
  fs <- basis$fs
  t_orig <- seq(from = 0, to = (N - 1) / fs, by = 1 / fs)
  t_new  <- seconds(mrs_data)
  
  # create resampled basis 
  basis_resamp    <- basis
  basis_resamp$N  <- Npts(mrs_data)
  basis_resamp$fs <- fs(mrs_data)
  basis_resamp$data <- matrix(nrow = basis_resamp$N, ncol = length(basis$names))
  
  # inverse FT back to the time-domain
  basis$data <- ift_shift_mat(basis$data)
  
  # loop through basis elements and resample
  for (n in 1:length(basis$names)) {
    
    resamp_ele <- stats::spline(t_orig, Re(basis$data[,n]), xout = t_new)$y +
                  1i * stats::spline(t_orig, Im(basis$data[,n]), xout = t_new)$y
    basis_resamp$data[,n] <- resamp_ele
  }
  
  # back to freq. domain
  basis_resamp$data <- ft_shift_mat(basis_resamp$data)
  
  if (ref_freq_match) {
    basis_resamp <- shift_basis(basis_resamp, basis$ref - mrs_data$ref)
    basis_resamp$ref <- mrs_data$ref
  }
  
  return(basis_resamp)
}

#' Return a subset of the input basis.
#' @param basis input basis.
#' @param names basis set elements to keep in the returned object. 
#' @param invert set to true to return all basis elements except those given in
#' the names argument.
#' @return a subset of the input basis.
#' @export
get_basis_subset <- function(basis, names, invert = FALSE) {
  matches <- rep(FALSE, length(basis$names))
  for (n in 1:length(names)) {
    match <- basis$names == names[n]
    if (any(match)) {
      matches <- matches | match
    } else {
      stop(names[n], " not found")
    }
  }
  
  if (invert) matches <- !matches
  
  basis$data  <- basis$data[, matches, drop = FALSE]
  basis$names <- basis$names[matches]
  return(basis)
}

#' Scale a basis-set to be consistent with spant assumptions for water scaling.
#' 
#' For correct water scaling, spant assumes the time-domain amplitude (t = 0)
#' for a single proton is 0.5. Internally simulated basis-sets will be correctly
#' scaled, however imported basis-sets should be assumed to be un-scaled and 
#' this function should be used. Note that the singlet specified should only
#' contain one resonance, and that any additional signals (eg TSP or residual 
#' water) will result in incorrect scaling. Therefore, only simulated basis sets
#' are appropriate for use with this function.
#' 
#' @param basis basis set to be scaled.
#' @param name the name of the singlet to be used as a scaling reference.
#' @param protons the number of MRS visible protons contributing to the singlet
#' resonance.
#' @return a scaled basis.
#' @export
scale_basis_from_singlet <- function(basis, name, protons) {
  basis_singlet <- get_basis_subset(basis, name)
  mrs_singlet   <- basis2mrs_data(basis_singlet)
  sc_factor     <- get_td_amp(mrs_singlet, nstart = 2) / protons * 2
  basis$data    <- basis$data / as.numeric(sc_factor)
  return(basis)
}

#' Make a basis-set object from a directory containing LCModel formatted RAW
#' files.
#' @param dir_path path to the directory containing LCModel RAW files. One file
#' per signal.
#' @param ft transmitter frequency in Hz.
#' @param fs sampling frequency in Hz.
#' @param ref reference value for ppm scale.
#' @return a basis-set object.
#' @export
make_basis_from_raw <- function(dir_path, ft, fs, ref) {
  files <- list.files(dir_path)
  raw_files <- grep("\\.RAW$", files, value = TRUE, ignore.case = TRUE)
  raw_path <- file.path(dir_path, raw_files)
  mrs_data_list <- read_mrs(raw_path, ft = ft, fs = fs, ref = ref,
                            format = "lcm_raw")
  mrs_data_dynamic <- append_dyns(mrs_data_list)
  names <- tools::file_path_sans_ext(raw_files)
  basis <- mrs_data2basis(mrs_data_dynamic, names = names)
  return(basis)
}
