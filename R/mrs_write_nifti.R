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
  
  # reorder the dimensions to x, y, z, t, dynamics, coil 
  data_points <- aperm(data_points, c(1, 2, 3, 6, 4, 5))
  
  # add a 7th dimension
  dim(data_points) <- c(dim(data_points), 1)
  
  # convert to nii
  mrs_nii <- RNifti::asNifti(data_points)
  
  # get the geometry information
  affine  <- get_mrs_affine(mrs_data, 1.5, 1.5, 1.5)
  
  # voxel dimensions
  mrs_pixdim <- mrs_data$resolution[2:4]
  dwell_time <- mrs_data$resolution[7]
  mrs_nii$pixdim <- c(0, mrs_pixdim, dwell_time, 0, 0, 0)
  
  # set the qform
  mrs_nii <- RNifti::`qform<-`(mrs_nii, structure(affine, code = 2L))
  mrs_nii$qform_code <- 1
  
  # write nifti to disk
  RNifti::writeNifti(mrs_nii, fname)
  
  # get fname of the json sidecar file
  fname_json <- stringr::str_c(stringr::str_sub(fname, 1, -7), "json")
  
  # create the R list to be exported as json
  proton_gr <- 42.5774785182e6 # TODO other nuclei
  Bzero_t   <- mrs_data$ft / proton_gr
  json_list <- list(MagneticFieldStrength = Bzero_t,
                    EchoTime = mrs_data$te * 1e3) 
  
  export_json <- jsonlite::toJSON(json_list, pretty = TRUE, auto_unbox = TRUE)
  
  # write json to disk
  write(export_json, fname_json)
}