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
#' @return volume data as a nifti object.
#' @export
get_mrsi_voi <- function(mrs_data, target_mri) {
  affine <- get_mrs_affine(mrs_data)
  
  rows   <- dim(mrs_data$data)[2]
  cols   <- dim(mrs_data$data)[3]
  slices <- dim(mrs_data$data)[4]
  
  voi_dim <- c(rows * mrs_data$resolution[2],
               cols * mrs_data$resolution[3],
               slices * mrs_data$resolution[4])
  
  raw_data <- array(1, voi_dim)
  voi <- RNifti::retrieveNifti(raw_data)
  voi <- RNifti::`sform<-`(voi, structure(affine, code = 2L))
  
  if (missing(target_mri)) {
    RNifti::orientation(voi) <- "RAS"
    voi <- voi[]
  } else {
    voi <- resample_voi(voi, target_mri)
  }
  voi
}

#' Resample a VOI to match a target image space using nearest-neighbor
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
#' @return resampled image data as a nifti object.
#' @export
resample_img <- function(source, target) {
  res <- RNiftyReg::niftyreg.linear(source, target, nLevels = 0,
                                    init = diag(4))$image
  return(res)
}

#' Plot a volume as an image overlay.
#' @export
#' @param voi volume data as a nifti object.
#' @param mri image data as a nifti object.
#' @param flip_lr flip the image in the left-right direction.
plot_voi_overlay <- function(voi, mri, flip_lr = TRUE) {
  # check the image orientation etc is the same
  check_geom(voi, mri)
 
  # swap L/R direction 
  if (flip_lr) {
    voi <- nifti_flip_lr(voi)
    mri <- nifti_flip_lr(mri)
  }
  
  # get the centre of gravity coords
  vox_inds <- get_voi_cog(voi)
  plot_col <- add_alpha(grDevices::heat.colors(10), 0.4)
  mri_oro <- neurobase::robust_window(oro.nifti::nifti(mri))
  neurobase::ortho2(mri_oro, oro.nifti::nifti(voi), xyz = vox_inds, 
                    col.y = plot_col, zlim.y = c(1, 2))
}

#' Plot a volume as an overlay on a segmented brain volume.
#' @param voi volume data as a nifti object.
#' @param mri_seg segmented brain volume as a nifti object.
#' @param flip_lr flip the image in the left-right direction.
#' @export
plot_voi_overlay_seg <- function(voi, mri_seg, flip_lr = TRUE) {
  # check the image orientation etc is the same
  check_geom(voi, mri_seg)
  
  # swap L/R direction 
  if (flip_lr) {
    voi     <- nifti_flip_lr(voi)
    mri_seg <- nifti_flip_lr(mri_seg)
  }
  
  # get the centre of gravity coords
  vox_inds <- get_voi_cog(voi)
  
  pvs <- get_voi_seg(voi, mri_seg)
  table <- paste("WM\t\t=  ", sprintf("%.1f", pvs[["WM"]]), "%\nGM\t\t=  ", 
                 sprintf("%.1f", pvs[["GM"]]), "%\nCSF\t=  ", 
                 sprintf("%.1f", pvs[["CSF"]]), "%\nOther\t=  ", 
                 sprintf("%.1f", pvs[["Other"]]),'%', sep = "")
  
  plot_col <- add_alpha(grDevices::heat.colors(10), 0.4)
  neurobase::ortho2(oro.nifti::nifti(mri_seg), oro.nifti::nifti(voi),
                    xyz = vox_inds, col.y = plot_col, zlim.y = c(1, 2))
  
  graphics::par(xpd = NA)
  graphics::text(x = 0.55, y = 0.12, labels = c(table), col = "white", pos = 4, 
                 cex = 1.5)
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
  mask <- (psf > 0)
  vals <- mri_seg[mask]
  other <- sum(as.numeric(vals == 0) * psf[mask])
  csf   <- sum(as.numeric(vals == 1) * psf[mask])
  gm    <- sum(as.numeric(vals == 2) * psf[mask])
  wm    <- sum(as.numeric(vals == 3) * psf[mask])
  vec <- c(other, csf, gm, wm)
  vec <- 100 * vec / sum(vec)  # normalise as a %
  names(vec) <- c("Other", "CSF", "GM", "WM")
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

# generate an sform affine for nifti generation
get_mrs_affine <- function(mrs_data, x_pos = 1, y_pos = 1, z_pos = 1) {
  # l1 norm
  col_vec <- mrs_data$col_vec/sqrt(sum(mrs_data$col_vec ^ 2))
  row_vec <- mrs_data$row_vec/sqrt(sum(mrs_data$row_vec ^ 2))
  col_vec[1:2] <- col_vec[1:2] * -1
  row_vec[1:2] <- row_vec[1:2] * -1
  slice_vec <- crossprod_3d(row_vec, col_vec)
  slice_vec <- slice_vec / sqrt(sum(slice_vec ^ 2))
  pos_vec <- mrs_data$pos_vec
  pos_vec[1:2] <- pos_vec[1:2] * -1
  affine <- diag(4)
  affine[1:3, 1] <- col_vec
  affine[1:3, 2] <- row_vec
  affine[1:3, 3] <- slice_vec
  rows <- dim(mrs_data$data)[2]
  cols <- dim(mrs_data$data)[3]
  slices <- dim(mrs_data$data)[4]
  
  affine[1:3, 4] <- pos_vec -
                    (mrs_data$resolution[2] * (-(x_pos - 1) + 0.5)) * row_vec -
                    (mrs_data$resolution[3] * (-(y_pos - 1) + 0.5)) * col_vec -
                    (mrs_data$resolution[4] * (-(z_pos - 1) + 0.5)) * slice_vec
  return(affine)
}

# check two nifti images are in the same space
check_geom <- function(a, b) {
  eq_xform <- max(Mod(RNifti::xform(a) - RNifti::xform(b))) < 1e-6
  eq_dim <- identical(dim(a), dim(b))
  eq_pix_dim <- identical(RNifti::pixdim(a), RNifti::pixdim(b))
  
  if ( !eq_xform | !eq_dim | !eq_pix_dim ) {
    stop("Inconsistant image geometry.")
  }
} 

# VOI centre of gravity
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
#' @return resliced imaging data.
#' @export
reslice_to_mrs <- function(mri, mrs) {
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
  resample_img(mri, dummy)
}