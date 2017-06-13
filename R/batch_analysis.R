
plot_batch_analysis_svs <- function(results_list, output_f) {
  
  grDevices::pdf(output_f)
  for (result in results_list) {
    graphics::plot(result$fits[[1]], xlim = c(4,0.5), plt_title = TRUE, 
                   main = result$data$id)
  }
  grDevices::dev.off()
}

# batch_analysis_svs <- function(metab_files, ref_files = NULL, format, ids,
#                                output_dir = NULL,  mri_files = NULL,
#                                mri_seg_files = NA, tr = NA,
#                                method = "TARQUIN", basis_file = NULL,
#                                opts = NULL) {
#   
#   # TODO check input files exist
#   # TODO check all vectors are the same length
#   # TODO check ids are unique
#   
#   # create the output directory
#   if (!is.null(output_dir)) {
#     dir.create(file.path(output_dir))
#   }
#   
#   n_results <- length(metab_files)
#   
#   results_list <- list()
#   
#   for (n in 1:n_results) {
#     cat("Running analysis ", n, " of ", n_results, "\n")
#     mrs_data <- read_mrs_tqn(metab_files[n], ref_files[n], format = format,
#                              id = ids[n])
#     
#     results <- fit_mrs(mrs_data, method = method,
#                        basis = basis_file, opts = opts)
#     
#     # do pv correction
#     if (!is.null(mri_seg_files[1])) {
#       mrs_nii <- get_svs_voi(mrs_data)
#       seg_nii <- RNifti::readNifti(mri_seg_files[n])
#       voxel_pvcs <- get_vox_pvcs(mrs_data, seg_nii)
#       # TODO make the tr an option
#       results <- apply_pvc(results, voxel_pvcs, tr = 2.0) 
#     }
#     
#     results_list <- c(results_list, list(results))
#     
#     # generate some output files
#     if (!is.null(output_dir)) {
#       dir.create(file.path(output_dir, ids[n]))
#       grDevices::pdf(file.path(output_dir, ids[n], "fit.pdf"))
#       graphics::plot(results)
#       grDevices::dev.off()
#       
#       output_csv(results, file.path(output_dir, ids[n], "results.csv"))
#       
#       if (!is.null(mri_files[1])) {
#         mrs_nii <- get_svs_voi(mrs_data)
#         mri_nii <- RNifti::readNifti(mri_files[n])
#         grDevices::png(file.path(output_dir, ids[n], "vox_overlay.png"))
#         plot_vox_overlay(mrs_nii, mri_nii)      
#         grDevices::dev.off()
#       }
#       
#       if (!is.null(mri_seg_files[1])) {
#         grDevices::png(file.path(output_dir,ids[n], "vox_overlay_seg.png"))
#         plot_vox_overlay_seg(mrs_nii, seg_nii)      
#         grDevices::dev.off()
#         
#         output_csv(results, file.path(output_dir, ids[n], "results_pvc.csv"),
#                    pvc = TRUE)
#       }
#     }
#   }
#   return(results_list)
# }

combine_results <- function(results_list, pvc = FALSE) {
  n_results <- length(results_list)
  out_table <- data.frame()
  
  for (n in 1:n_results) {
    if (pvc == TRUE) {
      temp_res = results_list[[n]]$results_pvc
    } else {
      temp_res = results_list[[n]]$results
    }
    temp_res$id <- results_list[[n]]$data$id
    temp_res$fname <- results_list[[n]]$data$fname
    out_table <- rbind(out_table, temp_res)
  }
  
  return(out_table)
}
