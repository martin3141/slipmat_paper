library(spant)      # mrs processing
library(ggplot2)    # nice plots
library(ggfortify)  # nice plots
library(stringr)    # string stuff
library(cowplot)    # multipart figs
library(ragg)       # nice image export

# remove old variables
rm(list = ls())

paths <- Sys.glob("./deconv_results_3_comp/subj-??_run-?/*.nii.gz") |> sort()

subj   <- str_extract(paths, "subj.{3}")    |> str_sub(6, 8)
run    <- str_extract(paths, "run.{2}")     |> str_sub(5, 6)
tissue <- paths |> basename() |> str_remove(".nii.gz")

tissue <- factor(tissue, levels = c("gm", "cwm", "pwm"))
subj   <- factor(subj)

id_table <- data.frame(subj, run, tissue)

all_mrs <- read_mrs(paths)

mean_spec <- all_mrs |> mean_mrs_list()

rats_corr <- all_mrs |> rats(ref = mean_spec)

# find the broadest lw
lws <- rats_corr |> bc_constant(xlim = c(0, -1)) |> append_dyns() |> 
       peak_info() |> (\(x) x$fwhm_ppm)() |> max()

# round-up
max_lw <- ceiling(lws * 1000) / 1000

proc_spec <-
  rats_corr                        |> 
  bc_constant(xlim = c(0, -1))     |>
  set_lw(max_lw)                   |> 
  zf()                             |>
  crop_spec(c(4, 0.2))             |>
  bc_als(lambda = 1e3)             |>
  scale_spec()

proc_mat <- proc_spec |> mrs_data2mat() |> Re()

pca_res <- prcomp(proc_mat, scale. = FALSE)

pca_res$x[,1] <- -pca_res$x[,1]
pca_res$rotation[,1] <- -pca_res$rotation[,1]

# autoplot(pca_res, data = id_table, colour = 'subj', shape = 'tissue',
#          x = 1, y = 2)

p1 <- autoplot(pca_res, data = id_table, colour = 'subj', shape = 'tissue',
               x = 1, y = 2) + guides(shape = guide_legend(order = 1),
                                      col = guide_legend(order = 2))

pc1 <- pca_res$rotation[,1] |> vec2mrs_data(fs = fs(proc_spec[[1]]),
                                            ft = proc_spec[[1]]$ft,
                                            ref = proc_spec[[1]]$ref,
                                            fd = TRUE)

pc2 <- pca_res$rotation[,2] |> vec2mrs_data(fs = fs(proc_spec[[1]]),
                                            ft = proc_spec[[1]]$ft,
                                            ref = proc_spec[[1]]$ref,
                                            fd = TRUE)

p2 <- ~plot(pc2, y_scale = TRUE, xlim = c(4, 0.2), bl_lty = 1)

agg_tiff("fig_8.tiff", width = 1200, height = 550, scaling = 2)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 14) |> print()
dev.off()

# autoplot(pca_res, data = id_table, colour = 'subj', shape = 'run')
