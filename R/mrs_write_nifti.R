#' Write MRS data object to file in NIFTI format.
#' @param fname the filename of the output NIFTI MRS data.
#' @param mrs_data object to be written to file.
#' @export
write_mrs_nifti <- function(fname, mrs_data) {
  if (class(mrs_data) != "mrs_data") stop("data object is not mrs_data format")
  
  if (stringr::str_sub(fname, -7) != ".nii.gz") {
    stop("filename argument must end in .nii.gz")
  }
  
  # convert to nii
  data_points <- mrs_data$data
  
  # drop the first dummy dimension
  data_points <- abind::adrop(data_points, 1)
  
  # reorder the dimensions to X, Y, Z, FID, coil, dynamics, indirect
  data_points <- aperm(data_points, c(1, 2, 3, 6, 5, 4))
  
  # add a 7th dimension
  dim(data_points) <- c(dim(data_points), 1)
  
  # convert to nii
  mrs_nii <- RNifti::asNifti(data_points)
  
  # get the geometry information
  #affine  <- get_mrs_affine(mrs_data, 1.5, 1.5, 1.5) # old version
  affine  <- mrs_data$affine
  
  # voxel dimensions
  mrs_pixdim <- mrs_data$resolution[2:4]
  dwell_time <- mrs_data$resolution[7]
  dyn_interval <- mrs_data$resolution[5]
  mrs_nii$pixdim <- c(-1, mrs_pixdim, dwell_time, dyn_interval, 0, -1)
  
  # set the sform
  mrs_nii <- RNifti::`sform<-`(mrs_nii, structure(affine, code = 2L))
  mrs_nii$intent_name <- "mrs_v0_2"
  
  # set the nucleus to a default value if not specified in mrs_data
  if (!exists("nuc", where = mrs_data)) mrs_data$nuc <- def_nuc()
  
  # create the R list to be exported as json
  json_list <- list(TransmitterFrequency = mrs_data$ft / 1e6,
                    ResonantNucleus = mrs_data$nuc,
                    SpectralWidth = 1 / dwell_time,
                    EchoTime = jsonlite::unbox(mrs_data$te * 1e3)) 
  
  RNifti::extension(mrs_nii, 44) <- jsonlite::toJSON(json_list, digits = NA)
  
  # write nifti to disk
  RNifti::writeNifti(mrs_nii, fname, version = 2)
}