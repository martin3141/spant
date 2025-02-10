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
}