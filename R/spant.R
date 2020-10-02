#' spant: spectroscopy analysis tools.
#'
#' spant provides a set of tools for reading, visualising and processing 
#' Magnetic Resonance Spectroscopy (MRS) data.
#' 
#' To get started with spant, take a look at the introduction vignette:
#' 
#' `vignette("spant-intro", package="spant")`
#' 
#' Full list of vignettes:
#' 
#' `browseVignettes(package = "spant")`
#' 
#' Full list of functions:
#' 
#' `help(package = spant, help_type = "html")`
#' 
#' An online version of the documentation is available from:
#' 
#' \url{https://martin3141.github.io/spant/}
#' 
"_PACKAGE"

#' @useDynLib spant

.onLoad <- function(libname, pkgname) {
  # set default options unless already set
  # by the user
  op <- options()
  op.spant <- list(
    spant.def_fs  = 2000,
    spant.def_ft  = 127.8e6,
    spant.def_N   = 1024,
    spant.def_ref = 4.65,
    spant.def_nuc = "1H",
    spant.tqn_cmd = "tarquin",
    spant.lcm_cmd = file.path(Sys.getenv("HOME"), ".lcmodel", "bin", "lcmodel"), 
    spant.precomp_verbose = TRUE,
    spant.precomp_mode = "default"
  )
  toset <- !(names(op.spant) %in% names(op))
  if (any(toset)) options(op.spant[toset])
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("spant", utils::packageVersion("spant")))
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
#' @param cmd path to binary.
#' @export
set_tqn_cmd <- function(cmd) {
  options(spant.tqn_cmd = cmd)
}

#' Set the command to run the LCModel command-line program.
#' @param cmd path to binary.
#' @export
set_lcm_cmd <- function(cmd) {
  options(spant.lcm_cmd = cmd)
}

#' Return (and optionally modify using the input arguments) a list of the 
#' default acquisition parameters.
#' @param ft specify the transmitter frequency in Hz.
#' @param fs specify the sampling frequency in Hz.
#' @param N specify the number of data points in the spectral dimension.
#' @param ref specify the reference value for ppm scale.
#' @param nuc specify the resonant nucleus.
#' @return A list containing the following elements:
#' * ft transmitter frequency in Hz.
#' * fs sampling frequency in Hz.
#' * N number of data points in the spectral dimension.
#' * ref reference value for ppm scale.
#' * nuc resonant nucleus.
#' @export
def_acq_paras <- function(ft  = getOption("spant.def_ft"),
                          fs  = getOption("spant.def_fs"),
                          N   = getOption("spant.def_N"),
                          ref = getOption("spant.def_ref"),
                          nuc = getOption("spant.def_nuc")) {
  list(ft = ft, fs = fs, N = N, ref = ref, nuc = nuc)
}

#' Return the default reference value for ppm scale.
#' @return reference value for ppm scale.
#' @export
def_ref <- function() {
  options()$spant.def_ref
}

#' Return the default sampling frequency in Hz.
#' @return sampling frequency in Hz.
#' @export
def_fs <- function() {
  options()$spant.def_fs
}

#' Return the default transmitter frequency in Hz.
#' @return transmitter frequency in Hz.
#' @export
def_ft <- function() {
  options()$spant.def_ft
}

#' Return the default number of data points in the spectral dimension.
#' @return number of data points in the spectral dimension.
#' @export
def_N <- function() {
  options()$spant.def_N
}

#' Return the default nucleus.
#' @return number of data points in the spectral dimension.
#' @export
def_nuc <- function() {
  options()$spant.def_nuc
}

#' Set the default acquisition parameters.
#' @param ft transmitter frequency in Hz.
#' @param fs sampling frequency in Hz.
#' @param N number of data points in the spectral dimension.
#' @param ref reference value for ppm scale.
#' @export
set_def_acq_paras <- function(ft  = getOption("spant.def_ft"),
                              fs  = getOption("spant.def_fs"),
                              N   = getOption("spant.def_N"),
                              ref = getOption("spant.def_ref"),
                              nuc = getOption("spant.nuc")) {
  
  options(spant.def_ft  = ft)
  options(spant.def_fs  = fs)
  options(spant.def_N   = N)
  options(spant.def_ref = ref)
  options(spant.nuc     = nuc)
}