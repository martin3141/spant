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

#' Resample a VOI to match a target image space.
#' @param voi volume data as a nifti object.
#' @param mri image data as a nifti object.
#' @return volume data as a nifti object.
#' @export
resample_voi <- function(voi, mri) {
  reg_res <- RNiftyReg::niftyreg.linear(voi, mri, nLevels = 0, 
                                        interpolation = 0, init = diag(4))$image
}

#' Plot a volume as an image overlay.
#' @export
#' @param voi volume data as a nifti object.
#' @param mri image data as a nifti object.
plot_voi_overlay <- function(voi, mri) {
  # check the image orientation etc is the same
  check_geom(voi, mri)
  
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
#' @export
plot_voi_overlay_seg <- function(voi, mri_seg) {
  # check the image orientation etc is the same
  check_geom(voi, mri_seg)
  
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

# generate an sform affine for nifti generation
get_mrs_affine <- function(mrs_data) {
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
  affine[1:3, 4] <- pos_vec - mrs_data$resolution[2] / 2 * col_vec -
                    mrs_data$resolution[3] / 2 * row_vec -
                    mrs_data$resolution[4] / 2 * slice_vec 
  
  return(affine)
}

# check two nifti images are in the same space
check_geom <- function(a, b) {
  eq_xform <- identical(RNifti::xform(a), RNifti::xform(b))
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