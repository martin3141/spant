#' spant: spectroscopy analysis tools.
#'
#' spant provides a set of tools for reading, visualising and processing 
#' Magnetic Resonance Spectroscopy (MRS) data.
#'
#' To learn more about spant, start with the vignettes:
#' `browseVignettes(package = "spant")`
#' 
#' For a full list of functions:
#' `help(package=spant,help_type="html")`
#' 
"_PACKAGE"

# next chunk is for Rcpp and RcppEigen
#' @useDynLib spant
#' @import RcppEigen
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

#' Example MEGA-PRESS data with significant B0 drift.
"spant_mpress_drift"

.onLoad <- function(libname, pkgname) {
  # set default options unless already set
  # by the user
  op <- options()
  op.spant <- list(
    spant.def_fs  = 2000,
    spant.def_ft  = 127.8e6,
    spant.def_N   = 1024,
    spant.def_ref = 4.65,
    spant.tqn_cmd = "tarquin",
    spant.lcm_cmd = file.path(Sys.getenv("HOME"), ".lcmodel", "bin", "lcmodel") 
  )
  toset <- !(names(op.spant) %in% names(op))
  if (any(toset)) options(op.spant[toset])
  
  invisible()
}

#' Check the TARQUIN binary can be run
#' @export
check_tqn <- function() {
  result <- tryCatch({
    sys_res <- suppressWarnings(system(getOption("spant.tqn_cmd"),
                                       intern = TRUE, ignore.stderr = TRUE))
  }, error = function(e) {
    return(NA)
  })
  
  if (!is.na(result[1])) {
    sys_res <- suppressWarnings(system(getOption("spant.tqn_cmd"),
                                       intern = TRUE, ignore.stderr = TRUE))
    
    tqn_ver <- strsplit(sys_res[2],"\\s+")[[1]][6]
    cat(paste("TARQUIN version",tqn_ver ,"was found successfully."))
  } else {
    stop("TARQUIN software is not functioning with the following command setting:\n", 
          getOption("spant.tqn_cmd"), "\nTry changing the path with the 'set_tqn_cmd' function.")
  }
}

#' Check LCModel can be run
#' @export
check_lcm <- function() {
  if (file.exists(getOption("spant.lcm_cmd"))) {
    cat("LCModel program sucessfully found.") 
  } else {
    stop("LCModel program not found in the following location:\n", 
          getOption("spant.lcm_cmd"),
         "\nif in a non-standard location try changing with the 'set_lcm_cmd' function.")
  }
}


#' Set the command to run the TARQUIN command-line program.
#' @param cmd Path to binary.
#' @export
set_tqn_cmd <- function(cmd) {
  options(spant.tqn_cmd = cmd)
}

#' Set the command to run the LCModel command-line program.
#' @param cmd Path to binary.
#' @export
set_lcm_cmd <- function(cmd) {
  options(spant.lcm_cmd = cmd)
}

#' Return (and optionally modify using the input arguments) a list of the 
#' default acquisition parameters.
#' @param ft specify the transmitter frequency in Hz.
#' @param fs specity the sampling frequency in Hz.
#' @param N specify the number of data points in the spectral dimension.
#' @param ref specify the reference value for ppm scale.
#' @return A list containing the following elements:
#' * ft Transmitter frequency in Hz.
#' * fs Sampling frequency in Hz.
#' * N Number of data points in the spectral dimension.
#' * ref Reference value for ppm scale.
#' @export
def_acq_paras <- function(ft  = getOption("spant.def_ft"),
                          fs  = getOption("spant.def_fs"),
                          N   = getOption("spant.def_N"),
                          ref = getOption("spant.def_ref")) {
  list(ft = ft, fs = fs, N = N, ref = ref)
}

#' Return the default reference value for ppm scale.
#' @return Reference value for ppm scale.
#' @export
def_ref <- function() {
  options()$spant.def_ref
}

#' Return the default sampling frequency in Hz.
#' @return Sampling frequency in Hz.
#' @export
def_fs <- function() {
  options()$spant.def_fs
}

#' Return the default transmitter frequency in Hz.
#' @return Transmitter frequency in Hz.
#' @export
def_ft <- function() {
  options()$spant.def_ft
}

#' Return the default number of data points in the spectral dimension.
#' @return Number of data points in the spectral dimension.
#' @export
def_N <- function() {
  options()$spant.def_N
}

#' Set the default acquisition parameters.
#' @param ft Transmitter frequency in Hz.
#' @param fs Sampling frequency in Hz.
#' @param N Number of data points in the spectral dimension.
#' @param ref Reference value for ppm scale.
#' @export
set_def_acq_paras <- function(ft  = getOption("spant.def_ft"),
                              fs  = getOption("spant.def_fs"),
                              N   = getOption("spant.def_N"),
                              ref = getOption("spant.def_ref")) {
  
  options(spant.def_ft  = ft)
  options(spant.def_fs  = fs)
  options(spant.def_N   = N)
  options(spant.def_ref = ref)
}

