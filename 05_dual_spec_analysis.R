library(spant)      # mrs processing
library(ggplot2)    # nice plots
library(ggfortify)  # nice plots
library(stringr)    # string stuff
library(cowplot)    # multipart figs
library(ragg)       # nice image export

# remove old variables
rm(list = ls())

paths <- Sys.glob("./deconv_results/subj-??_run-?/??.nii.gz") |> sort()

subj   <- str_extract(paths, "subj.{3}")    |> str_sub(6, 8)
run    <- str_extract(paths, "run.{2}")     |> str_sub(5, 6)
tissue <- str_extract(paths, ".{2}.nii.gz") |> str_sub(1, 2)

id_table <- data.frame(subj, run, tissue)

all_mrs <- read_mrs(paths)

mean_spec <- all_mrs |> mean_mrs_list()

rats_corr <- all_mrs |> rats(ref = mean_spec, ret_corr_only = TRUE)

# Second phase of SLIPMAT
proc_spec <-
  rats_corr                     |> 
  bc_constant(xlim = c(0, -1))  |>
  set_lw(0.06)                  |>
  zf()                          |>
  crop_spec(c(4, 0.2))          |>
  bc_als(lambda = 1e3)          |>
  scale_spec()

proc_mat <- proc_spec |> mrs_data2mat() |> Re()

gm_mat <- proc_mat[c(TRUE,  FALSE), ]
wm_mat <- proc_mat[c(FALSE, TRUE),  ]

# combine spectral pairs
gm_wm_concat <- cbind(gm_mat, wm_mat)

# plot example feature vector

agg_tiff("feature_vec.tiff", width = 1000, height = 800, scaling = 2.5)
par(mar = c(4.5, 4.5, 2, 2))
gm_wm_concat[1,] |> plot(type = "l", ylab = "Intensity (au)", bty = "l")
dev.off()

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

agg_tiff("fig_5.tiff", width = 1800, height = 800, scaling = 2)
pc_loadings <- plot_grid(p2, p3, p4, p5, labels = c('B', 'C', 'D', 'E'),
                         label_size = 14) |> print()
plot_grid(p1, pc_loadings, labels = c('A'), label_size = 14) |> print()
dev.off()

F_vals_wm <- rep(NA, ncol(wm_mat))
p_vals_wm <- rep(NA, ncol(wm_mat))
for (n in 1:ncol(wm_mat)) {
  F_vals_wm[n] <- summary(aov(wm_mat[,n] ~ id_table[c(TRUE, FALSE),]$subj))[[1]][1,4]
  p_vals_wm[n] <- summary(aov(wm_mat[,n] ~ id_table[c(TRUE, FALSE),]$subj))[[1]][1,5]
}

F_vals_wm_spec <- F_vals_wm |> vec2mrs_data(fs = fs(proc_spec[[1]]),
                                            ft = proc_spec[[1]]$ft,
                                            ref = proc_spec[[1]]$ref,
                                            fd = TRUE)

F_vals_gm <- rep(NA, ncol(gm_mat))
p_vals_gm <- rep(NA, ncol(gm_mat))
for (n in 1:ncol(gm_mat)) {
  F_vals_gm[n] <- summary(aov(gm_mat[,n] ~ id_table[c(TRUE, FALSE),]$subj))[[1]][1,4]
  p_vals_gm[n] <- summary(aov(gm_mat[,n] ~ id_table[c(TRUE, FALSE),]$subj))[[1]][1,5]
}

F_vals_gm_spec <- F_vals_gm |> vec2mrs_data(fs = fs(proc_spec[[1]]),
                                            ft = proc_spec[[1]]$ft,
                                            ref = proc_spec[[1]]$ref,
                                            fd = TRUE)

p1 <- ~plot(F_vals_wm_spec, y_scale = TRUE, yaxis_lab = "F-statistic",
            xlim = c(4, 0.2))

p2 <- ~plot(F_vals_gm_spec, y_scale = TRUE, yaxis_lab = "F-statistic",
            xlim = c(4, 0.2))

agg_tiff("fig_6.tiff", width = 1400, height = 800, scaling = 2)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 14) |> print()
dev.off()

# fs <- fs(F_vals_gm_spec)
# N <- Npts(F_vals_gm_spec)
# ref <- F_vals_gm_spec$ref
# ft <- F_vals_gm_spec$ft
# metab <- get_mol_paras("naag") |> sim_mol(pul_seq = seq_slaser_ideal, fs = fs,
#                                        N = N, ref = ref, ft = ft, xlim = c(4, 0.2))
# stackplot(append_dyns(metab, F_vals_gm_spec), xlim = c(4, 0.2), y_offset = 20)

