#' Generate a SVS acquisition volume from an \code{mrs_data} object.
#' @param mrs_data MRS data.
#' @param target_mri optional image data to match the intended volume space.
#' @return volume data as a nifti object.
#' @export
get_svs_voi <- function(mrs_data, target_mri) {
  affine <- get_mrs_affine(mrs_data)
  raw_data <- array(1, c(mrs_data$resolution[2:4]))
  voi <- RNifti::retrieveNifti(raw_data)
  voi <- RNifti::`sform<-`(voi, structure(affine, code = 2L))
  
  if (missing(target_mri)) {
    warning("Target MRI data has not been specified.")  
  } else {
    voi <- resample_voi(voi, target_mri)
  }
  voi  
}

#' Generate a MRSI voxel from an \code{mrs_data} object.
#' @param mrs_data MRS data.
#' @param target_mri optional image data to match the intended volume space.
#' @param x_pos x voxel coordinate.
#' @param y_pos y voxel coordinate.
#' @param z_pos z voxel coordinate.
#' @return volume data as a nifti object.
#' @export
get_mrsi_voxel <- function(mrs_data, target_mri, x_pos, y_pos, z_pos) {
  affine <- get_mrs_affine(mrs_data, x_pos, y_pos, z_pos)
  raw_data <- array(1, c(mrs_data$resolution[2:4]))
  voi <- RNifti::retrieveNifti(raw_data)
  voi <- RNifti::`sform<-`(voi, structure(affine, code = 2L))
  
  if (missing(target_mri)) {
    warning("Target MRI data has not been specified.")  
  } else {
    voi <- resample_voi(voi, target_mri)
  }
  voi  
}

# Generate a MRSI voxel PSF from an \code{mrs_data} object.
# @param mrs_data MRS data.
# @param target_mri optional image data to match the intended volume space.
# @param x_pos x voxel coordinate.
# @param y_pos y voxel coordinate.
# @param z_pos z voxel coordinate.
# @return volume data as a nifti object.
# get_mrsi_voxel_xy_psf_old <- function(mrs_data, target_mri, x_pos, y_pos, z_pos) {
#   affine <- get_mrs_affine(mrs_data, x_pos - 0.5, y_pos - 0.5, z_pos)
#   nom_vox_res <- mrs_data$resolution[2:4]
#   nom_vox_res[1:2] <- nom_vox_res[1:2] * 2 # double the size in x-y direction
# 
#   psf_mat <- pracma::repmat(signal::hamming(nom_vox_res[1]), nom_vox_res[1], 1)
#   psf_mat <- psf_mat * t(psf_mat)
#   psf_vec <- rep(psf_mat, nom_vox_res[3])
#   vox_psf <- array(psf_vec, nom_vox_res)
#   
#   voi <- RNifti::retrieveNifti(vox_psf)
#   voi <- RNifti::`sform<-`(voi, structure(affine, code = 2L))
#   
#   if (missing(target_mri)) {
#     warning("Target MRI data has not been specified.")  
#   } else {
#     voi <- resample_voi(voi, target_mri)
#   }
#   voi  
# }

#' Generate a MRSI voxel PSF from an \code{mrs_data} object.
#' @param mrs_data MRS data.
#' @param target_mri optional image data to match the intended volume space.
#' @param x_pos x voxel coordinate.
#' @param y_pos y voxel coordinate.
#' @param z_pos z voxel coordinate.
#' @return volume data as a nifti object.
#' @export
get_mrsi_voxel_xy_psf <- function(mrs_data, target_mri, x_pos, y_pos, z_pos) {
  FOV <- mrs_data$resolution[2] * Nx(mrs_data)
  
  psf_spatial_factor <- 8 # create a psf with spatial extent 8 times the nominal
                          # voxel size
  
  aff_off <- psf_spatial_factor / 2 - 0.5 # affine position offset
  affine <- get_mrs_affine(mrs_data, x_pos - aff_off, y_pos - aff_off, z_pos)
  nom_vox_res <- mrs_data$resolution[2:4]
  nom_vox_res[1:2] <- nom_vox_res[1:2] * psf_spatial_factor

  start_pt <- FOV / 2 - nom_vox_res[1] / 2 + 1
  subset <- seq(start_pt, by = 1, length.out = nom_vox_res[1])
  
  psf_mat <- Re(get_2d_psf(FOV = FOV, mat_size = Nx(mrs_data)))[subset, subset]
  
  psf_vec <- rep(psf_mat, nom_vox_res[3])
  vox_psf <- array(psf_vec, nom_vox_res)
  
  voi <- RNifti::retrieveNifti(vox_psf)
  voi <- RNifti::`sform<-`(voi, structure(affine, code = 2L))
  
  if (missing(target_mri)) {
    warning("Target MRI data has not been specified.")  
  } else {
    voi <- resample_voi(voi, target_mri)
  }
  voi  
}

#' Generate a MRSI VOI from an \code{mrs_data} object.
#' @param mrs_data MRS data.
#' @param target_mri optional image data to match the intended volume space.
#' @param map optional voi intensity map.
#' @param ker kernel to rescale the map data to the target_mri. Default value is
#' mmand::boxKernel(), use mmand::mnKernel() for a smoothed map.
#' @return volume data as a nifti object.
#' @export
get_mrsi_voi <- function(mrs_data, target_mri = NULL, map = NULL,
                         ker = mmand::boxKernel()) {
  #affine <- get_mrs_affine(mrs_data)
  affine <- get_mrs_affine(mrs_data)
  
  rows   <- dim(mrs_data$data)[2]
  cols   <- dim(mrs_data$data)[3]
  slices <- dim(mrs_data$data)[4]
  
  voi_dim <- c(rows * mrs_data$resolution[2],
               cols * mrs_data$resolution[3],
               slices * mrs_data$resolution[4])
  
  voi_dim <- round(voi_dim)
  
  if (is.null(map)) {
    raw_data <- array(1, voi_dim)
  } else {
    # assume 2D xy map for now
    #resamp_map <- mmand::rescale(drop(map), 10, mmand::mnKernel())
    
    # deal with Infs
    if (any(is.infinite(map))) map[is.infinite(map)] <- NA
    
    # this is to keep dimensions consistent
    rounded_res <- voi_dim[1] / rows
    
    if (any(is.na(map)) && (ker$name != "box")) {
      resamp_map <- interpolate_nas(drop(map), rounded_res, 2)
    } else {
      resamp_map <- mmand::rescale(drop(map), rounded_res, ker)
    }
    
    if (dim(resamp_map)[1] != voi_dim[1]) {
      stop("get_mrsi_voi - inconsistent matrix dimensions found following interpolation.")
    }
    
    if (dim(resamp_map)[2] != voi_dim[2]) {
      stop("get_mrsi_voi - inconsistent matrix dimensions found following interpolation.")
    }
      
    raw_data <- array(resamp_map, voi_dim) 
  }
  
  voi <- RNifti::retrieveNifti(raw_data)
  voi <- RNifti::`sform<-`(voi, structure(affine, code = 2L))
  
  if (is.null(target_mri)) {
    RNifti::orientation(voi) <- "RAS"
    voi <- voi[]
  } else {
    voi <- resample_voi(voi, target_mri)
  }
  voi
}

#' Resample a VOI to match a target image space using nearest-neighbour
#' interpolation.
#' @param voi volume data as a nifti object.
#' @param mri image data as a nifti object.
#' @return volume data as a nifti object.
#' @export
resample_voi <- function(voi, mri) {
  RNiftyReg::niftyreg.linear(voi, mri, nLevels = 0, interpolation = 0,
                             init = diag(4))$image
}

#' Resample an image to match a target image space.
#' @param source image data as a nifti object.
#' @param target image data as a nifti object.
#' @param interp interpolation parameter, see nifyreg.linear definition.
#' @return resampled image data as a nifti object.
#' @export
resample_img <- function(source, target, interp = 3L) {
  res <- RNiftyReg::niftyreg.linear(source, target, nLevels = 0,
                                    init = diag(4),
                                    interpolation = interp)$image
  return(res)
}

#' Plot a volume as an image overlay.
#' @export
#' @param mri image data as a nifti object or path to data file.
#' @param voi volume data as a nifti object or path to data file.
#' @param export_path optional path to save the image in png format.
#' @param zlim underlay intensity limits.
#' @param ... additional arguments to the ortho3 function.
plot_voi_overlay <- function(mri, voi, export_path = NULL, zlim = NULL, ...) {
  
  if ("character" %in% class(mri)) mri <- RNifti::readNifti(mri)
  
  if ("character" %in% class(voi)) {
    mrs_data <- read_mrs(voi)
    voi <- get_svs_voi(mrs_data, mri)
  } else if (class(voi)[1] == "mrs_data") {
    voi <- get_svs_voi(voi, mri)
  }
  
  # check the image orientation etc is the same
  check_geom(voi, mri)
 
  # get the centre of gravity coords
  plot_col <- grDevices::heat.colors(10)
  
  if (is.null(zlim)) zlim <- stats::quantile(mri, probs = c(0, 0.999))
  
  if (!is.null(export_path)) grDevices::png(export_path)
  ortho3(mri, voi, col_ol = plot_col, zlim = zlim, zlim_ol = c(0.99, 2),
         alpha = 0.4, colourbar = FALSE, ...)
  if (!is.null(export_path)) grDevices::dev.off()
}

#' Plot a volume as an overlay on a segmented brain volume.
#' @param mri_seg segmented brain volume as a nifti object.
#' @param voi volume data as a nifti object.
#' @param export_path optional path to save the image in png format.
#' @param ... additional arguments to the ortho3 function.
#' @export
plot_voi_overlay_seg <- function(mri_seg, voi, export_path = NULL, ...) {
  
  if ("character" %in% class(mri_seg)) mri_seg <- RNifti::readNifti(mri_seg)
  
  if ("character" %in% class(voi)) {
    mrs_data <- read_mrs(voi)
    voi <- get_svs_voi(mrs_data, mri_seg)
  } else if (class(voi)[1] == "mrs_data") {
    voi <- get_svs_voi(voi, mri_seg)
  }
  
  # check the image orientation etc is the same
  check_geom(voi, mri_seg)
  
  pvs <- get_voi_seg(voi, mri_seg)
  table <- paste("WM\t\t=  ", sprintf("%.1f", pvs[["WM"]]), "%\nGM\t\t=  ", 
                 sprintf("%.1f", pvs[["GM"]]), "%\nCSF\t=  ", 
                 sprintf("%.1f", pvs[["CSF"]]), "%\nOther\t=  ", 
                 sprintf("%.1f", pvs[["Other"]]),'%', sep = "")
  
  plot_col <- add_alpha(grDevices::heat.colors(10), 0.4)
  
  if (!is.null(export_path)) grDevices::png(export_path)
  
  ortho3(mri_seg, voi, col_ol = plot_col, zlim_ol = c(0.99, 2), alpha = 0.4,
         colourbar = FALSE, zlim = range(mri_seg), ...)
  
  graphics::par(xpd = NA)
  graphics::text(x = 0.55, y = 0.22, labels = c(table), col = "white", pos = 4, 
                 cex = 1.2)
  
  if (!is.null(export_path)) grDevices::dev.off()
  
  return(pvs)
}

#' Return the white matter, gray matter and CSF composition of a volume.
#' @param voi volume data as a nifti object.
#' @param mri_seg segmented brain volume as a nifti object.
#' @return a vector of partial volumes expressed as percentages.
#' @export
get_voi_seg <- function(voi, mri_seg) {
  # check the image orientation etc is the same
  check_geom(voi, mri_seg)
  vals <- mri_seg[voi == 1]
  pvs <- summary(factor(vals, levels = c(0, 1, 2, 3), 
         labels = c("Other", "CSF", "GM", "WM"))) / sum(voi) * 100
  pvs <- as.data.frame(t(pvs))
  return(pvs)
}

#' Return the white matter, gray matter and CSF composition of a volume.
#' @param psf volume data as a nifti object.
#' @param mri_seg segmented brain volume as a nifti object.
#' @return a vector of partial volumes expressed as percentages.
#' @export
get_voi_seg_psf <- function(psf, mri_seg) {
  # check the image orientation etc is the same
  check_geom(psf, mri_seg)
  
  # check the segmentation data has only integer values
  if (!all(mri_seg == as.integer(mri_seg))) {
    stop("Segmentation file has non-integer values.") 
  }
  
  mask <- (abs(psf) > 0)
  vals <- mri_seg[mask]
  
  
  other <- sum(as.numeric(vals == 0) * psf[mask])
  csf   <- sum(as.numeric(vals == 1) * psf[mask])
  gm    <- sum(as.numeric(vals == 2) * psf[mask])
  wm    <- sum(as.numeric(vals == 3) * psf[mask])
  vec <- c(other, csf, gm, wm)
  vec <- 100 * vec / sum(vec)  # normalise as a %
  names(vec) <- c("Other", "CSF", "GM", "WM")
  vec <- as.data.frame(t(vec))
  return(vec)
}

#' Convert SPM style segmentation files to a single categorical image where
#' the numerical values map as: 0) Other, 1) CSF, 2) GM and 3) WM.
#' @param fname any of the segmentation files (eg c1_MY_T1.nii).
#' @return nifti object.
#' @export
spm_pve2categorical <- function(fname) {
  # check the file name makes sense
  if (substr(basename(fname), 1, 1) != "c") stop("Error, filename does not begin with 'c'.")
   
  if (!(substr(basename(fname), 2, 2) %in% c("1","2","3","4","5","6"))) {
    stop("Error, filename does not follow expected pattern.")
  }
  
  # check all files are present and correct
  for (n in 1:6) {
    base <- basename(fname)
    substr(base,2, 2) <- as.character(n)
    fname <- file.path(dirname(fname), base)
    if (!file.exists(fname)) stop(paste("Error, file does not exist :",fname))
  }
  
  cat("Reading segmentation images...\n")
  # read the first file 
  substr(base,2, 2) <- "1"
  fname <- file.path(dirname(fname), base)
  x <- RNifti::readNifti(fname)
  
  # new structure
  new_dim <- c(dim(x),6)
  combined <- array(rep(NA, prod(new_dim)), dim = new_dim)
  combined[,,,1] <- x[]
  
  for (n in 2:6) {
    substr(base,2,2) <- as.character(n)
    fname <- file.path(dirname(fname), base)
    x <- RNifti::readNifti(fname)
    combined[,,,n] <- x[]
  }
  
  cat("Combining segmentation images...\n")
  # convert probabilties to catagorical
  # combine c4, c5, c6 to improve speed
  combined[,,,4] <- combined[,,,4] + combined[,,,5] + combined[,,,6]
  combined <- combined[,,,1:4]
  cat <- apply(combined, c(1,2,3), which.max) # THIS IS SLOOW
  
  # convert labels to FSL fast style output
  
  # SPM convention
  # c1 - GM; c2 - WM; c3 - CSF; c4,c5,c6 - OTHER 
  
  # FSL convention
  # 0 - OTHER; 1 - CSF; 2 - GM; 3 - WM 
  
  cat("Remap to FSL convention...\n")
  cat_fsl <- array(rep(NA, prod(new_dim[1:3])), dim = new_dim[1:3])
  cat_fsl <- ifelse(cat == 1, 2, cat_fsl)
  cat_fsl <- ifelse(cat == 2, 3, cat_fsl)
  cat_fsl <- ifelse(cat == 3, 1, cat_fsl)
  cat_fsl <- ifelse(cat == 4, 0, cat_fsl)
  cat("Done\n")
  
  x[] <- cat_fsl
  x
}

#' Generate an affine for nifti generation.
#' @param mrs_data input data.
#' @param x_pos x_position coordinate.
#' @param y_pos y_position coordinate.
#' @param z_pos z_position coordinate.
#' @return affine matrix.
get_mrs_affine <- function(mrs_data, x_pos = 1, y_pos = 1, z_pos = 1) {
  affine <- mrs_data$affine
  
  affine[1:3, 1] <- l2_norm_vec(affine[1:3, 1])
  affine[1:3, 2] <- l2_norm_vec(affine[1:3, 2])
  affine[1:3, 3] <- l2_norm_vec(affine[1:3, 3])
  affine[1:3, 4] <- affine[1:3, 4] - mrs_data$affine[1:3, 1] / 2 -
                                     mrs_data$affine[1:3, 2] / 2 -
                                     mrs_data$affine[1:3, 3] / 2
  
  affine[1:3, 4] <- affine[1:3, 4] + (x_pos - 1) * mrs_data$affine[1:3, 1] +
                                     (y_pos - 1) * mrs_data$affine[1:3, 2] +
                                     (z_pos - 1) * mrs_data$affine[1:3, 3]
  
  return(affine)
}

# check two nifti images are in the same space
check_geom <- function(a, b) {
  eq_xform <- max(Mod(RNifti::xform(a) - RNifti::xform(b))) < 1e-5
  eq_dim <- identical(dim(a), dim(b))
  eq_pix_dim <- identical(RNifti::pixdim(a), RNifti::pixdim(b))
  
  if ( !eq_xform | !eq_dim | !eq_pix_dim ) {
    print(RNifti::xform(a))
    print(RNifti::xform(b))
    print(dim(a))
    print(dim(b))
    print(RNifti::pixdim(a))
    print(RNifti::pixdim(b))
    stop("Inconsistant image geometry.")
  }
} 

#' Calculate the centre of gravity for an image containing 0 and 1's.
#' @param voi nifti object.
#' @return triplet of x,y,z coordinates.
#' @export
get_voi_cog <- function(voi) {
  as.integer(colMeans(which(voi == 1, arr.ind = TRUE)))
}

#' Flip the x data dimension order of a nifti image. This corresponds to
#' flipping MRI data in the left-right direction, assuming the data in save in
#' neurological format (can check with fslorient program).
#' @param x nifti object to be processed.
#' @return nifti object with reversed x data direction.
#' @export
nifti_flip_lr <- function(x) {
  lr_dim <- dim(x)[1]
  x[] <- x[lr_dim:1,,]
  return(x)
}

#' Reslice a nifti object to match the orientation of mrs data.
#' @param mri nifti object to be resliced.
#' @param mrs mrs_data object for the target orientation.
#' @param interp interpolation parameter, see nifyreg.linear definition.
#' @return resliced imaging data.
#' @export
reslice_to_mrs <- function(mri, mrs, interp = 3L) {
  mrs_affine <- RNifti::xform(get_mrsi_voi(mrs))
  dummy <- mri
  new_affine <- RNifti::xform(mri)
  new_affine[1:3, 1:3] <- mrs_affine[1:3, 1:3]
  # find the central point in the mri
  centre <- new_affine[1:3, 4] + RNifti::xform(mri)[1:3, 1:3] %*% 
                                 (dim(mri)[1:3] * RNifti::pixdim(mri)[1:3]) / 2
  
  new_pos <- centre - new_affine[1:3, 1:3] %*% (dim(mri)[1:3] *
                                                RNifti::pixdim(mri)[1:3]) / 2
  
  new_affine[1:3, 4] <- new_pos
  #RNifti::sform(dummy) <- new_affine
  dummy <- RNifti::`sform<-`(dummy, structure(new_affine, code = 2L))
  resample_img(mri, dummy, interp)
}

#' Calculate the partial volume estimates for each voxel in a 2D MRSI dataset.
#' 
#' Localisation is assumed to be perfect in the z direction and determined by
#' the ker input in the x-y direction.
#' 
#' @param mrs_data 2D MRSI data with multiple voxels in the x-y dimension.
#' @param mri_seg MRI data with values corresponding to the segmentation class. 
#' Must be 1mm isotropic resolution.
#' @param ker MRSI PSF kernel in the x-y direction compatible with the mmand 
#' package, eg: mmand::shapeKernel(c(10, 10), type = "box").
#' @return a data frame of partial volume estimates and individual segmentation
#' maps.
#' @export
get_mrsi2d_seg <- function(mrs_data, mri_seg, ker) {
  # reformat seg data into the same space as the full MRSI voi
  voi <- get_mrsi_voi(mrs_data)
  mri_seg_crop <- resample_img(mri_seg, voi, interp = 0L)
  
  # check mri_seg_crop has only integer values
  if (!all(mri_seg_crop == as.integer(mri_seg_crop))) {
    stop("Segmentation file has non-integer values.") 
  }
  
  mri_seg_0 <- rowMeans(mri_seg_crop == 0, dims = 2)
  mri_seg_1 <- rowMeans(mri_seg_crop == 1, dims = 2)
  mri_seg_2 <- rowMeans(mri_seg_crop == 2, dims = 2)
  mri_seg_3 <- rowMeans(mri_seg_crop == 3, dims = 2)
  
  mri_seg_sum <- mri_seg_0 + mri_seg_1 + mri_seg_2 + mri_seg_3
  
  mri_seg_0 <- mri_seg_0 / mri_seg_sum
  mri_seg_1 <- mri_seg_1 / mri_seg_sum
  mri_seg_2 <- mri_seg_2 / mri_seg_sum
  mri_seg_3 <- mri_seg_3 / mri_seg_sum
  
  mri_seg_0_blur <- mmand::morph(mri_seg_0, ker, operator = "*", merge = "sum")
  mri_seg_1_blur <- mmand::morph(mri_seg_1, ker, operator = "*", merge = "sum")
  mri_seg_2_blur <- mmand::morph(mri_seg_2, ker, operator = "*", merge = "sum")
  mri_seg_3_blur <- mmand::morph(mri_seg_3, ker, operator = "*", merge = "sum")
  
  mri_seg_blur_sum <- mri_seg_0_blur + mri_seg_1_blur + mri_seg_2_blur +
                      mri_seg_3_blur
  
  mri_seg_0_blur <- mri_seg_0_blur / mri_seg_blur_sum
  mri_seg_1_blur <- mri_seg_1_blur / mri_seg_blur_sum
  mri_seg_2_blur <- mri_seg_2_blur / mri_seg_blur_sum
  mri_seg_3_blur <- mri_seg_3_blur / mri_seg_blur_sum
  
  x_res <- mrs_data$resolution[2]
  y_res <- mrs_data$resolution[3]
  
  x_seq <- seq(from = x_res / 2, by = x_res, length.out = Nx(mrs_data))
  y_seq <- seq(from = y_res / 2, by = y_res, length.out = Ny(mrs_data))
    
  mri_seg_0_final <- pracma::fliplr(pracma::flipud(mri_seg_0_blur[x_seq, y_seq]))
  mri_seg_1_final <- pracma::fliplr(pracma::flipud(mri_seg_1_blur[x_seq, y_seq]))
  mri_seg_2_final <- pracma::fliplr(pracma::flipud(mri_seg_2_blur[x_seq, y_seq]))
  mri_seg_3_final <- pracma::fliplr(pracma::flipud(mri_seg_3_blur[x_seq, y_seq]))
  
  # gray matter fraction
  GMF <- 100 * as.numeric(mri_seg_2_final) /
                     (as.numeric(mri_seg_2_final) + as.numeric(mri_seg_3_final))
  
  # gray matter + white matter
  GM_WM <- 100 * (as.numeric(mri_seg_2_final) + as.numeric(mri_seg_3_final))
  
  seg_table <- data.frame(Other = 100 * as.numeric(mri_seg_0_final),
                      CSF   = 100 * as.numeric(mri_seg_1_final),
                      GM    = 100 * as.numeric(mri_seg_2_final),
                      WM    = 100 * as.numeric(mri_seg_3_final),
                      GMF   = GMF,
                      GM_WM = GM_WM)
  
  map_dim <- c(1, Nx(mrs_data), Ny(mrs_data), 1, 1, 1)
  
  gm_map    <- array(seg_table$GM,    map_dim)
  wm_map    <- array(seg_table$WM,    map_dim)
  csf_map   <- array(seg_table$CSF,   map_dim)
  other_map <- array(seg_table$Other, map_dim)
  gmf_map   <- array(seg_table$GMF,   map_dim)
  gm_wm_map <- gm_map + wm_map
  
  res <- list(seg_table = seg_table, gm_map = gm_map, wm_map = wm_map,
              csf_map = csf_map, other_map = other_map, gmf_map = gmf_map,
              gm_wm_map = gm_wm_map) 
  
  return(res)
}