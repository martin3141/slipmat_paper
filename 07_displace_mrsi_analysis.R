library(spant)      # mrs processing
library(ggplot2)    # nice plots
library(ggfortify)  # nice plots
library(stringr)    # string stuff
library(cowplot)    # multipart figs
library(ragg)       # nice image export

# remove old variables
rm(list = ls())

paths <- Sys.glob("./displace_results_0mm/subj-??_run-?/??.nii.gz") |> sort()

subj   <- str_extract(paths, "subj.{3}")    |> str_sub(6, 8)
run    <- str_extract(paths, "run.{2}")     |> str_sub(5, 6)
tissue <- str_extract(paths, ".{2}.nii.gz") |> str_sub(1, 2)

id_table <- data.frame(subj, run, tissue)

all_mrs <- read_mrs(paths)

mean_spec <- all_mrs |> mean_mrs_list()

rats_corr <- all_mrs |> rats(ref = mean_spec, ret_corr_only = TRUE)

# find the broadest lw
lws <- rats_corr |> bc_constant(xlim = c(0, -1)) |> append_dyns() |> 
       peak_info() |> (\(x) x$fwhm_ppm)() |> max()

# round-up
max_lw <- ceiling(lws * 1000) / 1000

# Second phase of SLIPMAT
proc_spec <-
  rats_corr                     |> 
  bc_constant(xlim = c(0, -1))  |>
  set_lw(max_lw)                |>
  zf()                          |>
  crop_spec(c(4, 0.2))          |>
  bc_als(lambda = 1e3)          |>
  scale_spec()

proc_mat <- proc_spec |> mrs_data2mat() |> Re()

gm_mat <- proc_mat[c(TRUE,  FALSE), ]
wm_mat <- proc_mat[c(FALSE, TRUE),  ]

# combine spectral pairs
gm_wm_concat <- cbind(gm_mat, wm_mat)

pca_res <- prcomp(gm_wm_concat, scale. = FALSE)

vec_n <- proc_mat |> ncol()

pc1_gm <- pca_res$rotation[1:vec_n, 1] |> 
  vec2mrs_data(fs = fs(proc_spec[[1]]), ft = proc_spec[[1]]$ft,
               ref = proc_spec[[1]]$ref, fd = TRUE)

pc1_wm <- pca_res$rotation[(vec_n + 1):(vec_n * 2), 1] |> 
  vec2mrs_data(fs = fs(proc_spec[[1]]), ft = proc_spec[[1]]$ft,
               ref = proc_spec[[1]]$ref, fd = TRUE)

pc2_gm <- pca_res$rotation[1:vec_n, 2] |>
  vec2mrs_data(fs = fs(proc_spec[[1]]), ft = proc_spec[[1]]$ft,
               ref = proc_spec[[1]]$ref, fd = TRUE)

pc2_wm <- pca_res$rotation[(vec_n + 1):(vec_n * 2), 2] |>
  vec2mrs_data(fs = fs(proc_spec[[1]]), ft = proc_spec[[1]]$ft,
               ref = proc_spec[[1]]$ref, fd = TRUE)

p1 <- autoplot(pca_res, data = id_table[c(TRUE, FALSE),], colour = 'subj',
               x = 1, y = 2)

p2 <- ~plot(pc1_gm, xlim = c(4, 0.2), bl_lty = 1)
p3 <- ~plot(pc1_wm, xlim = c(4, 0.2), bl_lty = 1)
p4 <- ~plot(pc2_gm, xlim = c(4, 0.2), bl_lty = 1)
p5 <- ~plot(pc2_wm, xlim = c(4, 0.2), bl_lty = 1)

# agg_tiff("fig_5_XX.tiff", width = 1800, height = 800, scaling = 2)
# pc_loadings <- plot_grid(p2, p3, p4, p5, labels = c('B', 'C', 'D', 'E'),
#                          label_size = 14) |> print()
# plot_grid(p1, pc_loadings, labels = c('A'), label_size = 14) |> print()
# dev.off()


get_pca_displaced <- function(paths, pca_res) {
  paths <- paths |> sort()
  
  subj   <- str_extract(paths, "subj.{3}")    |> str_sub(6, 8)
  run    <- str_extract(paths, "run.{2}")     |> str_sub(5, 6)
  tissue <- str_extract(paths, ".{2}.nii.gz") |> str_sub(1, 2)
  
  id_table <- data.frame(subj, run, tissue)
  
  all_mrs <- read_mrs(paths)
  
  mean_spec <- all_mrs |> mean_mrs_list()
  
  rats_corr <- all_mrs |> rats(ref = mean_spec, ret_corr_only = TRUE)
  
  # Second phase of SLIPMAT
  proc_spec <-
    rats_corr                           |> 
    bc_constant(xlim = c(0, -1))        |>
    set_lw(max_lw, mask_narrow = FALSE) |>
    zf()                                |>
    crop_spec(c(4, 0.2))                |>
    bc_als(lambda = 1e3)                |>
    scale_spec()
  
  proc_mat <- proc_spec |> mrs_data2mat() |> Re()
  
  gm_mat <- proc_mat[c(TRUE,  FALSE), ]
  wm_mat <- proc_mat[c(FALSE, TRUE),  ]
  
  # combine spectral pairs
  gm_wm_concat <- cbind(gm_mat, wm_mat)
  
  pca_res_disp <- pca_res
  
  pca_res_disp$x <- predict(pca_res, gm_wm_concat)
  
  return(pca_res_disp)
}

paths <- Sys.glob("./displace_results_0mm/subj-??_run-?/??.nii.gz")
pca_res_disp_0 <- get_pca_displaced(paths, pca_res)

paths <- Sys.glob("./displace_results_2mm/subj-??_run-?/??.nii.gz")
pca_res_disp_2 <- get_pca_displaced(paths, pca_res)

paths <- Sys.glob("./displace_results_5mm/subj-??_run-?/??.nii.gz")
pca_res_disp_5 <- get_pca_displaced(paths, pca_res)

paths <- Sys.glob("./displace_results_8mm/subj-??_run-?/??.nii.gz")
pca_res_disp_8 <- get_pca_displaced(paths, pca_res)

paths <- Sys.glob("./displace_results_10mm/subj-??_run-?/??.nii.gz")
pca_res_disp_10 <- get_pca_displaced(paths, pca_res)

paths <- Sys.glob("./displace_results_15mm/subj-??_run-?/??.nii.gz")
pca_res_disp_15 <- get_pca_displaced(paths, pca_res)

p1 <- autoplot(pca_res_disp_0, data = id_table[c(TRUE, FALSE),],
               colour = 'subj', variance_percentage = FALSE,
               main = "0mm displacement") + theme(legend.position = "none")

p2 <- autoplot(pca_res_disp_2, data = id_table[c(TRUE, FALSE),],
               colour = 'subj', variance_percentage = FALSE,
               main = "2mm displacement") + theme(legend.position = "none")

p3 <- autoplot(pca_res_disp_5, data = id_table[c(TRUE, FALSE),],
               colour = 'subj', variance_percentage = FALSE,
               main = "5mm displacement")

p4 <- autoplot(pca_res_disp_8, data = id_table[c(TRUE, FALSE),],
               colour = 'subj', variance_percentage = FALSE,
               main = "8mm displacement") + theme(legend.position = "none")

p5 <- autoplot(pca_res_disp_10, data = id_table[c(TRUE, FALSE),],
               colour = 'subj', variance_percentage = FALSE,
               main = "10mm displacement") + theme(legend.position = "none")

p6 <- autoplot(pca_res_disp_15, data = id_table[c(TRUE, FALSE),],
               colour = 'subj', variance_percentage = FALSE,
               main = "15mm displacement")

agg_tiff("fig_7.tiff", width = 1400, height = 900, scaling = 2)
plot_grid(p1, p2, p3, p4, p5, p6, rel_widths = c(1, 1, 1.25))
dev.off()
