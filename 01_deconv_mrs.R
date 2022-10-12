library(spant)    # mrs processing
library(ggplot2)  # nice plots
library(readr)    # csv export
library(stringr)  # string handling
library(ragg)     # image export

# remove old variables
rm(list = ls())

# get a list of available subjects
subj_vec <- Sys.glob("../sub-??") |> str_sub(-2)
subj_n   <- length(subj_vec)

# we know there are 3 runs per subject
subj_vec <- rep(subj_vec, each = 3)
run_vec  <- rep(c("1", "2", "3"), subj_n)

# loop over all data
for (n in 1:length(subj_vec)) {
  subj <- subj_vec[n]
  run  <- run_vec[n]
  
  # generate paths
  mrs_path <- paste0("../sub-", subj, "/mrs/sub-", subj, "_run-", run,
                     "_mrsi.nii.gz")
  
  t1_path  <- paste0("../sub-", subj, "/anat/sub-", subj, "_T1w.nii.gz")
  
  seg_path <- paste0("../sub-", subj, "/anat/T1_brain_seg.nii.gz")
  
  mrs_data <- mrs_path |> read_mrs()
  t1_data  <- t1_path  |> readNifti()
  seg_data <- seg_path |> readNifti()
  
  # create an output folder
  output_folder <- paste0("deconv_results/", "subj-" ,subj, "_run-", run)
  
  dir.create(output_folder, recursive = TRUE)
  
  # extract the central 8x8 matrix of spectra
  mrs_data_cropped <- mrs_data |> crop_xy(8, 8)
  
  # crop the corner voxels as they are close to the diagonal saturation bands
  mask <- matrix(FALSE, 8, 8)
  mask[1, 1] <- mask[1, 8] <- mask[8, 1] <- mask[8, 8] <- TRUE
  mrs_data_cropped <- mrs_data_cropped |> mask_xy_mat(mask)
  
  ref <- get_voxel(mrs_data_cropped, x_pos = 5, y_pos = 5) |> align(2.01)
  
  mrs_data_rats <- spant::rats(mrs_data_cropped, ref = ref, xlim = c(4, 0.5))
  
  # plot the grid and save as an image
  grid_path <- file.path(output_folder, "rats_grid.tiff")
  agg_tiff(grid_path, width = 1000, height = 1025, scaling = 1.5)
  mrs_data_rats |> zf() |> gridplot(xlim = c(4, 0.5))
  dev.off()
  
  mrs_data_pre_proc <- mrs_data_rats |> zf() |> crop_spec(c(4, 0.2)) |>
                       bc_als(lambda = 1e3)
  
  # mrs_data_pre_proc <- mrs_data_cropped |> zf() |> crop_spec(c(4, 0.2)) |>
  #                      bc_als(lambda = 1e3)
  
  scale_res <- mrs_data_pre_proc |> scale_spec(ret_scale_factor = TRUE)
  
  # apply scaling factor to uncropped spectra
  mrs_data_scaled <- scale_mrs_amp(mrs_data_rats, scale_res$scale_factor)
  # mrs_data_scaled <- scale_mrs_amp(mrs_data_cropped, scale_res$scale_factor)
  
  mrs_data_pre_proc <- scale_res$mrs_data
  
  # stackplot(mrs_data_pre_proc, xlim = c(4.2, 0.5))
  
  # get the approximate MRSI excitation region
  mrsi_mask     <- get_mrsi_voi(mrs_data |> crop_xy(10, 10), seg_data)
  
  # mask the segmented data by the excitation region
  seg_data_crop <- seg_data * mrsi_mask
  
  # make a plot
  
  seg_img_path <- file.path(output_folder, "segmented_t1.tiff")
  agg_tiff(seg_img_path, width = 1000, height = 1025, scaling = 1.5)
  ortho3(t1_data, seg_data_crop, colourbar = FALSE, zlim_ol = c(0, 3),
         crosshairs = FALSE, orient_lab = FALSE)
  dev.off()
  
  # calculate the PSF for convolution with the segmentation data
  mat_size <- Nx(mrs_data_pre_proc)
  psf_ker  <- Re(get_2d_psf(FOV = mrs_data_pre_proc$resolution[2] * mat_size,
                            mat_size = mat_size))
  
  seg_res <- get_mrsi2d_seg(mrs_data_pre_proc, seg_data, psf_ker)
  
  # reslice MRI to match MRSI
  t1_data_resliced <- reslice_to_mrs(t1_data, mrs_data_pre_proc)
  
  gmf_map <- seg_res$gmf_map
  
  glu_tcr_map  <- int_spec(mrs_data_pre_proc, xlim = c(2.3, 2.4)) /
    int_spec(mrs_data_pre_proc, xlim = c(3.1, 2.9))
  
  tcho_tcr_map <- int_spec(mrs_data_pre_proc, xlim = c(3.1, 3.2)) /
    int_spec(mrs_data_pre_proc, xlim = c(3.1, 2.9))
  
  tnaa_tcr_map <- int_spec(mrs_data_pre_proc, xlim = c(1.9, 2.1)) /
    int_spec(mrs_data_pre_proc, xlim = c(3.1, 2.9))
  
  # plot(gmf_map, glu_tcr_map)
  # plot(gmf_map, tcho_tcr_map)
  # plot(gmf_map, tnaa_tcr_map)
  
  map_resamp <- get_mrsi_voi(mrs_data_pre_proc, t1_data_resliced, glu_tcr_map)
  map_resamp <- get_mrsi_voi(mrs_data_pre_proc, t1_data_resliced, gmf_map)
  
  # ortho3(t1_data_resliced, map_resamp)
  
  # sort by GMF and stackplot
  order <- sort(seg_res$seg_table$GMF, index.return = TRUE)$ix
  # mrs_data_pre_proc |> collapse_to_dyns() |> get_subset(dyn_set = order) |>
  #                      stackplot(y_offset = 10, x_offset = 0)
  
  spec <- spec_decomp(mrs_data_pre_proc, seg_res$seg_table$WM,
                      seg_res$seg_table$GM)
  
  spec_full <- spec_decomp(mrs_data_scaled, seg_res$seg_table$WM,
                           seg_res$seg_table$GM)
  
  write_mrs(spec_full$wm, file.path(output_folder, "wm.nii.gz"), force = TRUE)
  write_mrs(spec_full$gm, file.path(output_folder, "gm.nii.gz"), force = TRUE)
  
}