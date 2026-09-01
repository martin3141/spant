#' Install the spant command-line interface scripts to a system path.
#'
#' This should be run following each new install of spant to ensure consistency.
#' Typical command line usage : sudo Rscript -e "spant::install_cli()"
#' Note sudo resets the HOME environment variable by default, which can cause
#' "there is no package called 'spant'" errors if spant is only installed in
#' your personal library. Use "sudo -E Rscript -e \"spant::install_cli()\""
#' to preserve HOME and avoid this issue.
#'
#' @param path optional path to install the scripts. Defaults to :
#' "/usr/local/bin".
#' @export
install_cli <- function(path = NULL) {

  if (is.null(path)) path <- "/usr/local/bin"

  cli_scripts_dir <- system.file('cli_scripts', package = 'spant')
  cli_scripts <- list.files(cli_scripts_dir)

  ver_line <- paste0("ver <- \"", utils::packageVersion("spant"), "\"")

  for (script_name in cli_scripts) {
    script_path <- file.path(path, script_name)

    package_script_path <- file.path(cli_scripts_dir, script_name)

    # embed the spant version number into the cli script
    script_txt <- readLines(package_script_path)
    script_txt <- c("#!/usr/bin/env Rscript --vanilla", "", ver_line, "",
                    script_txt)

    # write to filesystem
    write(script_txt, file = script_path)

    # make executable
    Sys.chmod(script_path, mode = "755")

    cat(script_name, " successfully installed to : ", path, "\n", sep = "")
  }
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
