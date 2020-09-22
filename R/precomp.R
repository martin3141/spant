#' Save function results to file and load on subsequent calls to avoid repeat 
#' computation.
#' @param file file name to write the results.
#' @param fun function to run.
#' @param ... arguments to be passed to fun.
#' @export
precomp <- function(file, fun, ...) {
  
  # check if the filename has an extension, and add ".rds" if not 
  bnf <- basename(file)
  if (!grepl("\\.", bnf)) file <- paste0(file, ".rds")
  
  precomp_global_mode    <- getOption("spant.precomp_mode")
  precomp_global_verbose <- getOption("spant.precomp_verbose")
  
  # act based on the mode
  if (precomp_global_mode == "default") {
    if (file.exists(file)) {
      if (precomp_global_verbose) cat("reading precomputed result:", file)
      result <- readRDS(file)
    } else {
      result <- do.call(fun, list(...)) 
      if (precomp_global_verbose) cat("writing precomputed result:", file)
      saveRDS(result, file)
    }
  } else if (precomp_global_mode == "overwrite") {
    if (file.exists(file)) {
      result <- do.call(fun, list(...)) 
      if (precomp_global_verbose) cat("overwriting precomputed result:", file)
      saveRDS(result, file)
    } else {
      result <- do.call(fun, list(...)) 
      if (precomp_global_verbose) cat("writing precomputed result:", file)
      saveRDS(result, file)
    }
  } else if (precomp_global_mode == "clean") {
    if (file.exists(file)) {
      if (precomp_global_verbose) cat("reading and deleting precomputed result:", file)
      result <- readRDS(file)
      file.remove(file)
    } else {
      stop(paste0("precomputed file not found: ", file))
    }
  } else if (precomp_global_mode == "disabled") {
    result <- do.call(fun, list(...)) 
  } else {
    stop("I don't belong here.")
  }
      
  return(result)
}

#' Set the precompute mode.
#' @param mode can be one of: "default", "overwrite", "clean" or "disabled".
#' @export
set_precomp_mode <- function(mode = NA) {
  
  if (!((mode == "default") | (mode == "overwrite") | (mode == "clean") |
        (mode == "disabled"))) {
    
    stop("function argument must be one of \"default\", \"overwrite\", \"clean\" or \"disabled\"")
  }
  
  options(spant.precomp_mode = mode)
}

#' Set the verbosity of the precompute function.
#' @param verbose can be TRUE or FALSE.
#' @export
set_precomp_verbose <- function(verbose = NA) {
  
  if (!is.logical(verbose)) stop("function argument must be logical (TRUE or FALSE)")
  
  if (is.na(verbose)) stop("function argument must be logical (TRUE or FALSE)")
  
  options(spant.precomp_verbose = verbose)
}