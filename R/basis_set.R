#' Generate a basis file using TARQUIN.
#' @param basis_file The filename of the basis file to be generated.
#' @param metab_data The MRS data object to match the generated basis parameters.
#' @param opts List of options to pass to TARQUIN.
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
#' @param opts List of options to pass to TARQUIN.
#' @examples
#' \dontrun{
#' write_basis_tqn('test.basis',mrs_data,c("--echo","0.04"))
#' }
#' @export
sim_basis_tqn <- function(fs = 2000, ft = 128e6, N = 1024, ref = 4.65, 
                          opts = NULL) {
  
  metab_data <- sim_zeros(fs = fs, ft = ft, N = N, ref = ref)
  
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
  cat(paste(c("Trans. freq (MHz)       : ", x$ft * 1e-6, "\n")), sep = "")
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
#' @param basis_file Path to basis file.
#' @param ref Assumed ppm reference value.
#' @return Basis object.
#' @export
read_basis <- function(basis_file, ref = 4.65) {
  con  <- file(basis_file, open = "r")
  names <- vector()
  data <- vector()
  
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (startsWith(line, " NDATAB = ")) {
      N <- as.integer(strsplit(trimws(line), "\\s+")[[1]][3])
      data_lines <- ceiling(2*N/6) # assume 6 cols
    } else if (startsWith(line, " HZPPPM = ")) {
      bas_ft <- strsplit(trimws(line), "\\s+")[[1]][3]
      bas_ft <- as.double(gsub(",", "", bas_ft))*1e6
    } else if (startsWith(line, " BADELT = ")) {
      bas_fs <- strsplit(trimws(line), "\\s+")[[1]][3]
      bas_fs <- 1/as.double(gsub(",", "",bas_fs))
    } else if (endsWith(line, "$BASIS")) {
      while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
        if (startsWith(line, " ID = ")) {
          id <- (strsplit(trimws(line), " ")[[1]][3])
          id <- gsub(",", "",id)
          id <- gsub("'", "",id)
          names <- c(names, id)
        } else if (endsWith(line, "$END")) {
          x <- utils::read.fortran(con, "6F13.0", n = data_lines)
          data_pts <- as.vector(t(as.matrix(x)))
          data_pts <- data_pts[seq(1, 2 * N, 2)] + 1i*data_pts[seq(2, 2 * N, 2)]
          data_pts <- pracma::ifftshift(data_pts)
          data <- cbind(data,data_pts)
          break
        }
      }
    }
  }
  close(con)
  
  basis_set <- list(data = data, N = N, fs = bas_fs, ft = bas_ft, 
                    names = names, ref = ref)
  
  class(basis_set) <- "basis_set"
  sort_basis(basis_set)
}

#' Write a basis object to an LCModel .basis formatted file.
#' @param basis Basis object to be exported.
#' @param basis_file Path to basis file to be generated.
#' @export
write_basis <- function(basis, basis_file) {
  mrs_data <- basis2mrs_data(basis)
  N <- N(mrs_data)
  
  sink(basis_file)
  cat(" $SEQPAR\n", sep = "")
  cat(" FWHMBA =  ", 0.00621699, " ,\n", sep = "")
  cat(" HZPPPM =  ",mrs_data$ft * 1e-6,",\n", sep = "")
  cat(" ECHOT =  ",mrs_data$te * 1000,",\n", sep = "")
  cat(" SEQ = 'PRESS' $END\n", sep = "")
  cat(" $BASIS1\n", sep = "")
  cat(" IDBASI = 'SPANT',\n", sep = "")
  cat(" FMTBAS = '(6E13.5)',\n", sep = "")
  cat(" BADELT =  ", 1 / fs(mrs_data), ",\n", sep = "")
  cat(" NDATAB = ", N(mrs_data), " $END\n",sep = "")
  
  for (n in 1:dyns(mrs_data)) {
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
  
basis2mrs_data <- function(basis) {
  mat2mrs_data(basis$data, fs = basis$fs, ft = basis$ft, ref = basis$ref, 
               fd = TRUE)
}

mrs_data2basis <- function(mrs_data, names) {
  # transform to FD
  if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
  }
  
  data <- t(as.matrix(mrs_data$data[1, 1, 1, 1,, 1,]))
  basis_set <- list(data = data, N = N(mrs_data), fs = fs(mrs_data), 
                    ft = mrs_data$ft, names = names, ref = mrs_data$ref)
  
  class(basis_set) <- "basis_set"
  basis_set
}

sort_basis <- function(basis) {
  names_sorted <- sort(basis$names, index.return = TRUE)
  basis$names <- names_sorted$x
  basis$data <- basis$data[, names_sorted$ix]
  basis
}