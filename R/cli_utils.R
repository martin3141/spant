#' Install the spant command-line interface scripts to a system path.
#' 
#' @param path optional path to install the scripts. Defaults to : 
#' "/usr/local/bin".
#' 
#' @export
install_cli <- function(path = NULL) {
  
  if (is.null(path)) {
    path <- "/usr/local/bin"
  }
  
  fit_svs_path <- file.path(path, "spant_fit_svs")
    
  package_fit_svs_path <- system.file('cli_scripts', 'spant_fit_svs',
                                      package = 'spant')
  
  # copy script from package directory
  file.copy(package_fit_svs_path, fit_svs_path, overwrite = TRUE)
  
  # make executable
  Sys.chmod(fit_svs_path, mode = "755")
  
  cat("spant_fit_svs sucessfully installed to : ", path, "\n", sep = "")
}

# https://www.r-bloggers.com/2015/06/identifying-the-os-from-r/
get_os <- function() {
  sysinf <- Sys.info()
  if (!is.null(sysinf)) {
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "macos"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "macos"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
