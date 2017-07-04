mrs_ica <- function(mrs_data, n, xlim = NULL, scale = "ppm", mode = "real",
                    ref = 1, slice = 1, mask_cutoff = 20,
                    mask_map = NULL, ...) {
  
  # transform to FD
  if (!is_fd(mrs_data)) {
      mrs_data <- td2fd(mrs_data)
  }
  
  if (scale == "ppm") {
    x_scale <- ppm(mrs_data)
  } else if (scale == "hz") {
    x_scale <- hz(mrs_data)
  } else if (scale == "points") {
    x_scale <- pts(mrs_data)
  }
  
  if (is.null(xlim)) {
    xlim <- c(x_scale[1], x_scale[N(mrs_data)])
  }
  
  subset <- get_seg_ind(x_scale, xlim[1], xlim[2])
  
  old_dim <- dim(mrs_data)
  new_dim <- c(1 * old_dim[2] * old_dim[3] * length(slice) * old_dim[5] *
               old_dim[6], length(subset)) 
  
  data_reshape <- array(mrs_data$data[ref,,, slice,,, subset], dim = new_dim)
  
  if (mode == "real") {
    data_reshape <- Re(data_reshape)
  } else if (mode == "imag") {
    data_reshape <- Im(data_reshape)
  } else if (mode == "abs") {
    data_reshape <- Mod(data_reshape)
  }
  
  if (is.null(mask_map)) {
    res <- fastICA::fastICA(data_reshape, n, ...) 
    mask_map_mask = NULL
  } else {
    mask_map_mask <- mask_map > max(mask_map[ref,,, slice,,]) * 
                                mask_cutoff / 100
    
    ica_mask <- array(mask_map_mask[ref,,, slice,,],dim = new_dim[1])
    res <- fastICA::fastICA(data_reshape[ica_mask,], n, ...)  
    res$Srecon <- array(NA, dim = c(dim(data_reshape)[1], n))
    res$Srecon[ica_mask,] <- res$S
  }
  
  # flip the spectra that have negative data points
  neg_inds <- which(apply(res$A, 1, max) + apply(res$A, 1, min) < 0)
  res$A[neg_inds,] <- -res$A[neg_inds,]
  # invert coeffs to match (S matrix)
  res$S[,neg_inds] <- -res$S[,neg_inds]
  
  list(res = res, x_scale = x_scale[subset], mask_map = mask_map_mask)
}