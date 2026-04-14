#' Return the spant resources directory.
#' 
#' Usually used to store additional software, datasets and templates. An attempt
#' will be made to create the folder if not found. The default resources
#' directory is "spant_resources" located in the users HOME directory.
#' @return full path to the spant resources directory.
#' @export
get_spant_resources_dir <- function() {
  
 deps_dir <- file.path(path.expand("~"), "spant_resources")
 
 # create if doesn't exist
 if (!dir.exists(deps_dir)) dir.create(deps_dir)
 
 return(deps_dir)
}

#' Install ANTs / ANTsX from : https://github.com/ANTsX/ANTs/releases
#' @param platform see the releases page for supported platforms. Common
#' platforms include : "ubuntu-24.04-X64-gcc" or "macos-14-ARM64-clang".
#' @param version ANTs version to install, defaults to "2.6.5".
#' @export
install_ants <- function(platform, version = "2.6.5") {
  
  spant_resources <- get_spant_resources_dir()
 
  dl_file <- file.path(spant_resources, paste0("ants-",version, ".zip"))
  
  url <- paste0("https://github.com/ANTsX/ANTs/releases/download/v",
                 version, "/ants-", version, "-", platform, ".zip")
  utils::download.file(url, destfile = dl_file, mode = "wb")
  
  # use unzip rather than internal to preserve execute permissions
  utils::unzip(dl_file, exdir = spant_resources, unzip = "unzip")
  
  # delete zip file
  file.remove(dl_file)
}

#' Install the Oasis brain template from : 
#' https://doi.org/10.6084/m9.figshare.915436.v2
#' @export
install_oasis_template <- function() {
  
  spant_resources <- get_spant_resources_dir()
 
  dl_file <- file.path(spant_resources, "Oasis.zip")
  url <- "https://ndownloader.figshare.com/files/3133832"
  utils::download.file(url, destfile = dl_file, mode = "wb")
  
  utils::unzip(dl_file, exdir = spant_resources)
  
  # delete zip file
  file.remove(dl_file)
}

#' Set the ANTs installation directory location.
#' @param dir ANTs installation directory.
#' @export
set_ants_dir <- function(dir) {
  dir <- path.expand(dir)
  
  if (!dir.exists(dir)) stop("ANTs directory does not exist.")
  
  eg_util <- file.path(dir, "bin", "antsBrainExtraction.sh")
  
  if (!file.exists(eg_util)) stop("antsBrainExtraction.sh not found")
  
  options(spant.ants_dir = dir)
}

#' Return the ANTs installation directory, or throw an error if not found.
#' 
#' Will check and return the "spant.ants_dir" option set by set_ants_dir. If
#' not set, will search the spant_resources directory for ANTs and return the most 
#' recent version if multiple are found.
#' 
#' @return ANTs installation directory.
#' @export
get_ants_dir <- function() {
  if (is.null(getOption("spant.ants_dir"))) {
    spant_resources <- get_spant_resources_dir()
    ants_dirs <- Sys.glob(file.path(spant_resources, "ants-?.?.?"))
    
    if ((length(ants_dirs)) == 0) {
      stop("ANTs not found, try using set_ants_dir() or install_ants() to rectify.")
    } else if ((length(ants_dirs)) == 1) {
      return(ants_dirs[1])
    } else {
      return(sort(ants_dirs, decreasing = TRUE)[1])
    }
  } else {
    return(getOption("spant.ants_dir")) 
  }
}