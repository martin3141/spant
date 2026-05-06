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
#' platforms include : "macos-14-ARM64-clang", "almalinux9-X64-gcc",
#' "centos7-X64-gcc" or "ubuntu-24.04-X64-gcc". Note, It may be necessary to
#' increase the timeout with "options(timeout = 1000)" on slower connections.
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

#' Install FaceOff from : https://github.com/srikash/FaceOff, requires ANTs to
#' be installed to work.
#' @export
install_faceoff <- function() {
  
  spant_resources <- get_spant_resources_dir()
 
  dl_file <- file.path(spant_resources, "FaceOff-2.0.zip")
  url <- "https://github.com/srikash/FaceOff/archive/refs/tags/2.0.zip"
  utils::download.file(url, destfile = dl_file, mode = "wb")
  
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

#' Segment T1 weighted MRI data using the ANTs binaries and write to file.
#' @param mri_path path to the volumetric T1 data.
#' @param out_dir optional output directory. Defaults to the same directory
#' as mri_path if not specified.
#' @export
segment_t1_ants <- function(mri_path, out_dir = NULL) {
  
  if (is.null(out_dir)) {
    dir_path <- dirname(mri_path)
  } else {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    dir_path <- out_dir
  }
  
  resources_dir <- get_spant_resources_dir()
  
  # check ANTs is available
  ants_dir <- get_ants_dir()
  
  # check Oasis template is available
  oasis_dir <- file.path(resources_dir, "MICCAI2012-Multi-Atlas-Challenge-Data")
  
  if (!dir.exists(oasis_dir)) {
    stop("Oasis brain template directory not found. Try 'install_oasis_template()' to install.")
  }
  
  temp_path <- tempfile()
  
  # create an empty file to keep the directory alive
  file.create(paste0(temp_path, "_temp.txt"))
  
  brain_extraction_template <- file.path(oasis_dir, "T_template0.nii.gz")
  brain_extraction_probability_mask <- file.path(oasis_dir,
                            "T_template0_BrainCerebellumProbabilityMask.nii.gz")
  brain_extraction_registration_mask <- file.path(oasis_dir,
                           "T_template0_BrainCerebellumRegistrationMask.nii.gz")

  args <- paste0("-d 3 -k 1 -c 3x1x2x3 -a ",
                 mri_path," -e ",
                 brain_extraction_template, " -m ",
                 brain_extraction_probability_mask, " -f ",
                 brain_extraction_registration_mask, " -o ",
                 temp_path)
  
  env <- paste0("PATH=", ants_dir, "/bin:/usr/bin:/bin")

  system2(command = "antsBrainExtraction.sh", args = args, env = env)
  
  brain <- paste0(temp_path, "BrainExtractionBrain.nii.gz")
  file.copy(from = brain, to = file.path(dir_path, "t1_brain.nii.gz"))
  
  csf <- paste0(temp_path, "BrainExtractionCSF.nii.gz")
  # file.copy(from = csf, to = file.path(dir_path, "t1_csf.nii.gz"))
  
  wm  <- paste0(temp_path, "BrainExtractionWM.nii.gz")
  # file.copy(from = wm, to = file.path(dir_path, "t1_wm.nii.gz"))
  
  gm  <- paste0(temp_path, "BrainExtractionGM.nii.gz")
  # file.copy(from = gm, to = file.path(dir_path, "t1_gm.nii.gz"))
  
  # seg  <- paste0(temp_path, "BrainExtractionSegment.nii.gz")
  # file.copy(from = seg, to = file.path(dir_path, "t1_wm_gm.nii.gz"))
  
  # ImageMath 3 Summed_MRI.nii.gz + MRI_1.nii.gz MRI_2.nii.gz
  
  # combine CSF, WM and GM into a single volume
  seg  <- file.path(dir_path, "t1_seg.nii.gz")
  args <- paste0("3 ", seg, " + ", wm, " ", gm)
  system2(command = "ImageMath", args = args, env = env)
  
  args <- paste0("3 ", seg, " + ", csf, " ", seg)
  system2(command = "ImageMath", args = args, env = env)
}

#' Segment T1 weighted MRI data using the rpyANTs interface to ANTs and write to
#' file.
#' @param mri_path path to the volumetric T1 data.
#' @param out_dir optional output directory. Defaults to the same directory
#' as mri_path if not specified.
#' @export
segment_t1_rpyants <- function(mri_path, out_dir = NULL) {
  
  if (is.null(out_dir)) {
    dir_path <- dirname(mri_path)
  } else {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    dir_path <- out_dir
  }
  
  # antspynet brain extraction step followed by recommended AtroposN4 method
  # from :
  # https://github.com/ntustison/antsAtroposN4Example/blob/master/antsAtroposN4Command.R
  
  t1      <- rpyANTs::ants$image_read(mri_path)
  t1_mask <- rpyANTs::antspynet_brain_extraction(t1, modality = "t1")
  
  outer_iters <- 5
  weight_mask <- NULL
  
  for (n in 1:outer_iters) {
    
    t1_n4 <- rpyANTs::ants$n4_bias_field_correction(t1, t1_mask,
                                                    weight_mask = weight_mask)
    atropos_seg <- rpyANTs::ants$atropos(a = t1_n4, x = t1_mask,
                                         m = '[0.2,1x1x1]')
    
    if (n != outer_iters) {
      # only use gm and wm probabilities for weight mask
      weight_mask <- atropos_seg$probabilityimages[[1]] *
        (1 - atropos_seg$probabilityimages[[0]]) *
        (1 - atropos_seg$probabilityimages[[2]]) +
        atropos_seg$probabilityimages[[2]] *
        (1 - atropos_seg$probabilityimages[[0]]) *
        (1 - atropos_seg$probabilityimages[[1]])
    }
  }
  
  t1_nii               <- readNifti(mri_path)
  atropos_n4_seg_out   <- t1_nii
  atropos_n4_seg_out[] <- atropos_seg$segmentation[]
  writeNifti(atropos_n4_seg_out, file.path(dir_path, "t1_seg.nii.gz"))
}

# #' Segment T1 weighted MRI data using the ANTsR interface to ANTs and write to
# #' file.
# #' @param mri_path path to the volumetric T1 data.
# #' @param out_dir optional output directory. Defaults to the same directory
# #' as mri_path if not specified.
# segment_t1_antsr <- function(mri_path, out_dir = NULL) {
#   
#   if (is.null(out_dir)) {
#     dir_path <- dirname(mri_path)
#   } else {
#     dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#     dir_path <- out_dir
#   }
#   
#   # brain extraction step followed by recommended AtroposN4 method from :
#   # https://github.com/ntustison/antsAtroposN4Example/blob/master/antsAtroposN4Command.R
#   
#   t1         <- ANTsR::antsImageRead(mri_path)
#   t1_n4_init <- ANTsR::abpN4(t1)
#   
#   tem     <- ANTsR::antsImageRead("~/spant_templates/oasis/T_template0.nii.gz")
#   temmask <- ANTsR::antsImageRead("~/spant_templates/oasis/T_template0_BrainCerebellumProbabilityMask.nii.gz")
#   
#   brain   <- ANTsR::abpBrainExtraction(img = t1, tem = tem,
#                                        temmask = temmask, regtype = "SyN")
#   
#   outer_iters <- 5
#   weight_mask <- NULL
#   
#   for (n in 1:outer_iters) {
#     
#     t1_n4 <- ANTsR::n4BiasFieldCorrection(t1, mask = brain$bmask,
#                                           weightMask = weight_mask)
#     atropos_seg <- ANTsR::atropos(a = t1_n4, x = brain$bmask, m = '[0.2,1x1x1]')
#     
#     if (n != outer_iters) {
#       # only use gm and wm probabilities for weight mask
#       weight_mask <- atropos_seg$probabilityimages[[2]] *
#         (1 - atropos_seg$probabilityimages[[1]]) *
#         (1 - atropos_seg$probabilityimages[[3]]) +
#         atropos_seg$probabilityimages[[3]] *
#         (1 - atropos_seg$probabilityimages[[1]]) *
#         (1 - atropos_seg$probabilityimages[[2]])
#     }
#   }
#   
#   t1_nii               <- readNifti(mri_path)
#   atropos_n4_seg_out   <- t1_nii
#   atropos_n4_seg_out[] <- atropos_seg$segmentation[]
#   writeNifti(atropos_n4_seg_out, file.path(dir_path, "t1_seg.nii.gz"))
# }

#' Deface a T1 weighted head scan using the FaceOff method described in :
#' https://github.com/srikash/FaceOff. Requires ANTs to be installed.
#' @param mri_path path to the volumetric T1 data.
#' @param out_dir optional output directory. Defaults to the same directory
#' as mri_path if not specified.
#' @export
faceoff <- function(mri_path, out_dir = NULL) {
  
  if (is.null(out_dir)) {
    dir_path <- dirname(mri_path)
  } else {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    dir_path <- out_dir
  }
  
  resources_dir <- get_spant_resources_dir()
  
  # check ANTs is available
  ants_dir <- get_ants_dir()
  
  # check FaceOff is available
  faceoff_dir <- file.path(resources_dir, "FaceOff-2.0")
  
  if (!dir.exists(faceoff_dir)) {
    stop("FaceOff directory not found. Try 'install_faceoff()' to install.")
  }
  
  command <- "FaceOff"
  args    <- paste0("-i ", mri_path)
  
  antspath_env <- paste0("ANTSPATH=", ants_dir, "/bin")
  path_env     <- paste0("PATH=", faceoff_dir, ":", ants_dir,
                         "/bin:/usr/bin:/bin")
  
  env <- c(antspath_env, path_env, "ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1")
  
  system2(command = command, args = args, env = env)
  
  file_base <- tools::file_path_sans_ext(mri_path, compression = TRUE)
  
  defaced_out     <- paste0(file_base, "_defaced.nii.gz")
  deface_mask_out <- paste0(file_base, "_defaceMask.nii.gz")
  out_name        <- file.path(dir_path, "mri_deface.nii.gz")
  
  file.copy(defaced_out, out_name)
  file.remove(deface_mask_out)
}