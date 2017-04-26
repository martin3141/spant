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
    spant.lcm_cmd = "lcmodel" 
  )
  toset <- !(names(op.spant) %in% names(op))
  if (any(toset)) options(op.spant[toset])
  
  invisible()
}

#' Set the path to the TARQUIN command-line program.
#' @param path Path to binary.
#' @export
set_tqn_path <- function(path) {
  options(spant.tqn_cmd = path)
}

#' Set the path to the LCModel command-line program.
#' @param path Path to binary.
#' @export
set_lcm_path <- function(path) {
  options(spant.lcm_cmd = path)
}

#' Return a list of the default acquisition parameters.
#' @return A list containing the following elements:
#' * ft Trasmitter frequency in Hz.
#' * fs Sampling frequency in Hz.
#' * N Number of data points in the spectral dimension.
#' * ref Reference value for ppm scale.
#' @export
get_def_acq_paras <- function() {
  op <- options()
  list(ft = op$spant.def_ft, fs = op$spant.def_fs, N = op$spant.def_N,
      ref = op$spant.def_ref)
}

#' Set the default acquisition parameters.
#' @param ft Trasmitter frequency in Hz.
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

