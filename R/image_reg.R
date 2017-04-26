
get_svs_voi <- function(mrs_data, mrsi_data = NA) {
  affine <- get_mrs_affine(mrs_data)
  raw_data <- array(1, c(mrs_data$resolution[2:4]))
  nifti_data <- RNifti::retrieveNifti(raw_data)
  #RNifti::sform(nifti_data) <- structure(affine, code = 2L)
  RNifti::`sform<-`(nifti_data,structure(affine, code = 2L))
  #RNifti::pixunits(nifti_data) <- c("mm")
  affine[1:3, 4] = c(50, 50, 50)
  return(nifti_data)
}

resample_svi_voi <- function(svs_voi, target) {
  reg_res <- RNiftyReg::niftyreg.linear(svs_voi, target, nLevels = 0, 
                                        interpolation = 0, init = diag(4))$image
}

# stolen from the interweb
add_alpha <- function(col, alpha = 1){
  if (missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, grDevices::col2rgb) / 255, 2, 
        function(x) 
          grDevices::rgb(x[1], x[2], x[3], alpha = alpha))  
}

plot_vox_overlay <- function(mrs, mri) {
  # TODO test if the mrs needs resampling before doing it
  res_mrs <- resample_svi_voi(mrs, mri)
  vox_inds <- get_vox_cog(res_mrs)
  plot_col <- add_alpha(grDevices::heat.colors(10), 0.4)
  neurobase::ortho2(oro.nifti::nifti(mri), oro.nifti::nifti(res_mrs),
                    xyz = vox_inds, col.y = plot_col, zlim.y = c(1, 2))
}

plot_vox_overlay_seg <- function(mrs, mri) {
  # TODO test if the mrs needs resampling before doing it
  res_mrs <- resample_svi_voi(mrs, mri)
  vox_inds <- get_vox_cog(res_mrs)
  pvs <- get_vox_seg(res_mrs, mri)
  table <- paste("WM\t\t=  ", sprintf("%.1f", pvs[["WM"]]), "%\nGM\t\t=  ", 
                 sprintf("%.1f", pvs[["GM"]]), "%\nCSF\t=  ", 
                 sprintf("%.1f", pvs[["CSF"]]), "%\nOther\t=  ", 
                 sprintf("%.1f", pvs[["Other"]]),'%', sep = "")
  
  plot_col <- add_alpha(grDevices::heat.colors(10), 0.4)
  neurobase::ortho2(oro.nifti::nifti(mri), oro.nifti::nifti(res_mrs),
                    xyz = vox_inds, col.y = plot_col, zlim.y = c(1, 2))
  
  graphics::par(xpd = NA)
  graphics::text(x = 0.55, y = 0.12, labels = c(table), col = "white", pos = 4, 
                 cex = 1.5)
  x = 1
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

get_vox_pvcs <- function(mrs_data, seg_mri) {
  mrs_nii <- get_svs_voi(mrs_data)
  res_mrs <- resample_svi_voi(mrs_nii, seg_mri)
  pvs <- get_vox_seg(res_mrs, seg_mri)
  return(pvs)  
}

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

get_vox_cog <- function(vox_data) {
  as.integer(colMeans(which(vox_data == 1, arr.ind = TRUE)))
}

get_vox_seg <- function(vox_data, seg_data) {
  vox_num = sum(vox_data)
  vals <- seg_data[vox_data == 1]
  pvs <- summary(factor(vals, levels = c(0, 1, 2, 3), 
        labels = c("Other", "CSF", "GM", "WM"))) / sum(vox_data) * 100
  return(pvs)
}

# Compute the vector cross product between x and y, and return the components
# indexed by i. Stolen from: 
# http://stackoverflow.com/questions/15162741/what-is-rs-crossproduct-function
crossprod_3d <- function(x, y, i = 1:3) {
  # Project inputs into 3D, since the cross product only makes sense in 3D.
  To3D <- function(x) utils::head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)
  
  # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
  # so on).  Index3D() lets us do that using R's convention of 1-based (rather
  # than 0-based) arrays.
  Index3D <- function(i) (i - 1) %% 3 + 1
  
  # The i'th component of the cross product is:
  # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat the indices cyclically.
  return(x[Index3D(i + 1)] * y[Index3D(i + 2)] -
            x[Index3D(i + 2)] * y[Index3D(i + 1)])
}