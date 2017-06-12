#' @export
get_svs_voi <- function(mrs_data, target_mri = NA) {
  affine <- get_mrs_affine(mrs_data)
  raw_data <- array(1, c(mrs_data$resolution[2:4]))
  voi <- RNifti::retrieveNifti(raw_data)
  voi <- RNifti::`sform<-`(voi, structure(affine, code = 2L))
  
  if (is.na(target_mri)) {
    warning("Target MRI data has not been specified.")  
  } else {
    voi <- resample_voi(voi, target_mri)
  }
  
  voi  
}

#' @export
resample_voi <- function(voi, mri) {
  reg_res <- RNiftyReg::niftyreg.linear(voi, mri, nLevels = 0, 
                                        interpolation = 0, init = diag(4))$image
}

#' @export
plot_svs_overlay <- function(svs_voi, mri_data) {
  # TODO check the image orientation etc is the same
  # get the centre of gravity coords
  vox_inds <- get_voi_cog(svs_voi)
  plot_col <- add_alpha(grDevices::heat.colors(10), 0.4)
  neurobase::ortho2(oro.nifti::nifti(mri_data), oro.nifti::nifti(svs_voi),
                    xyz = vox_inds, col.y = plot_col, zlim.y = c(1, 2))
}

#' @export
plot_svs_overlay_seg <- function(svs_voi, mri_data) {
  # TODO check the image orientation etc is the same
  # get the centre of gravity coords
  vox_inds <- get_voi_cog(svs_voi)
  
  pvs <- get_voi_seg(svs_voi, mri_data)
  table <- paste("WM\t\t=  ", sprintf("%.1f", pvs[["WM"]]), "%\nGM\t\t=  ", 
                 sprintf("%.1f", pvs[["GM"]]), "%\nCSF\t=  ", 
                 sprintf("%.1f", pvs[["CSF"]]), "%\nOther\t=  ", 
                 sprintf("%.1f", pvs[["Other"]]),'%', sep = "")
  
  plot_col <- add_alpha(grDevices::heat.colors(10), 0.4)
  neurobase::ortho2(oro.nifti::nifti(mri_data), oro.nifti::nifti(svs_voi),
                    xyz = vox_inds, col.y = plot_col, zlim.y = c(1, 2))
  
  graphics::par(xpd = NA)
  graphics::text(x = 0.55, y = 0.12, labels = c(table), col = "white", pos = 4, 
                 cex = 1.5)
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

#' @export
get_vox_pvcs <- function(mrs_data, seg_mri) {
  mrs_nii <- get_svs_voi(mrs_data)
  res_mrs <- resample_voi(mrs_nii, seg_mri)
  pvs <- get_voi_seg(res_mrs, seg_mri)
  return(pvs)  
}

#' @export
apply_pvc <- function(result, pvcs, tr){
  te = result$data$te
  B0 = round(result$data$ft / 42.58e6,1)
  corr_factor <- get_corr_factor(te, tr, B0, pvcs[["GM"]], pvcs[["WM"]],
                                 pvcs[["CSF"]])
  
  amp_cols = result$amp_cols
  default_factor = 35880 * 0.7
  result$results$GM_vol = pvcs[["GM"]]
  result$results$WM_vol = pvcs[["WM"]]
  result$results$CSF_vol = pvcs[["CSF"]]
  result$results$Other_vol = pvcs[["Other"]]
  result$results_pvc <- result$results
  pvc_cols <- 6:(5 + amp_cols * 2)
  result$results_pvc[, pvc_cols] <- result$results_pvc[, pvc_cols] /
                                    default_factor * corr_factor
  
  # append tables with %GM, %WM, %CSF and %Other
  
  return(result)
}

get_voi_cog <- function(voi) {
  as.integer(colMeans(which(voi == 1, arr.ind = TRUE)))
}

get_voi_seg <- function(voi, seg_data) {
  vals <- seg_data[voi == 1]
  pvs <- summary(factor(vals, levels = c(0, 1, 2, 3), 
        labels = c("Other", "CSF", "GM", "WM"))) / sum(voi) * 100
  return(pvs)
}