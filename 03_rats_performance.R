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

id_table_rats <- data.frame(subj, run, tissue)

all_mrs_rats <- read_mrs(paths)

# all_mrs_rats |> append_dyns() |> peak_info() |> (\(x) x$fwhm_ppm)() |> max()

# all_mrs_rats |> bc_constant(c(0, -1)) |> append_dyns() |> peak_info() |> 
#                 (\(x) x$fwhm_ppm)() |> max()

# basic bl correction and scaling
all_mrs_rats <- all_mrs_rats |>
                bc_constant(c(0, -1)) |> 
                scale_spec(xlim = c(1.8, 2.2), operator = "max")

paths <- Sys.glob("./deconv_results_no_rats/subj-??_run-?/??.nii.gz") |> sort()

subj   <- str_extract(paths, "subj.{3}")    |> str_sub(6, 8)
run    <- str_extract(paths, "run.{2}")     |> str_sub(5, 6)
tissue <- str_extract(paths, ".{2}.nii.gz") |> str_sub(1, 2)

id_table_no_rats <- data.frame(subj, run, tissue)

all_mrs_no_rats <- read_mrs(paths)

# basic bl correction and scaling
all_mrs_no_rats <- all_mrs_no_rats |>
                   bc_constant(c(0, -1)) |> 
                   scale_spec(xlim = c(1.8, 2.2), operator = "max")

all_mrs_rats_gm <- all_mrs_rats[c(TRUE, FALSE)]
all_mrs_rats_wm <- all_mrs_rats[c(FALSE, TRUE)]

all_mrs_no_rats_gm <- all_mrs_no_rats[c(TRUE, FALSE)]
all_mrs_no_rats_wm <- all_mrs_no_rats[c(FALSE, TRUE)]

# measure linewidths
rats_gm_peak_info <- all_mrs_rats_gm |> append_dyns() |> peak_info()
rats_wm_peak_info <- all_mrs_rats_wm |> append_dyns() |> peak_info()
no_rats_gm_peak_info <- all_mrs_no_rats_gm |> append_dyns() |> peak_info()
no_rats_wm_peak_info <- all_mrs_no_rats_wm |> append_dyns() |> peak_info()

# measure SNR
rats_gm_snr <- all_mrs_rats_gm |> append_dyns() |> calc_spec_snr()
rats_wm_snr <- all_mrs_rats_wm |> append_dyns() |> calc_spec_snr()
no_rats_gm_snr <- all_mrs_no_rats_gm |> append_dyns() |> calc_spec_snr()
no_rats_wm_snr <- all_mrs_no_rats_wm |> append_dyns() |> calc_spec_snr()

# plot typical improvement for gm/wm rats vs no rats
eg_spec <- 22

p1 <- function() {
  plot(all_mrs_no_rats_gm[[eg_spec]] |> zf(), xlim = c(4, 0.2),
       restore_def_par = FALSE)
  snr <- no_rats_gm_snr[eg_spec]
  lw  <- no_rats_gm_peak_info$fwhm_ppm[eg_spec]
  string <- paste0("LW = ", round(lw, 3), " ppm\nSNR = ", round(snr))
  text(1.8, 0.6, string, adj = 0)
  text(1.35, 0.22, "Lac")
}

p2 <- function() {
  plot(all_mrs_no_rats_wm[[eg_spec]] |> zf(), xlim = c(4, 0.2),
       restore_def_par = FALSE)
  snr <- no_rats_wm_snr[eg_spec]
  lw  <- no_rats_wm_peak_info$fwhm_ppm[eg_spec]
  string <- paste0("LW = ", round(lw, 3), " ppm\nSNR = ", round(snr))
  text(1.8, 0.6, string, adj = 0)
}

p3 <- function() {
  plot(all_mrs_rats_gm[[eg_spec]] |> zf(), xlim = c(4, 0.2),
       restore_def_par = FALSE)
  snr <- rats_gm_snr[eg_spec]
  lw  <- rats_gm_peak_info$fwhm_ppm[eg_spec]
  string <- paste0("LW = ", round(lw, 3), " ppm\nSNR = ", round(snr))
  text(1.8, 0.6, string, adj = 0)
  text(1.35, 0.20, "Lac")
  text(2.45, 0.45, expression(symbol("\255")*Glu))
  text(3.6, 0.46, expression(symbol("\255")*Ins))
  text(3.3, 0.62, expression(symbol("\257")*tCho))
}

p4 <- function() {
  plot(all_mrs_rats_wm[[eg_spec]] |> zf(), xlim = c(4, 0.2),
       restore_def_par = FALSE)
  snr <- rats_wm_snr[eg_spec]
  lw  <- rats_wm_peak_info$fwhm_ppm[eg_spec]
  string <- paste0("LW = ", round(lw, 3), " ppm\nSNR = ", round(snr))
  text(1.8, 0.6, string, adj = 0)
  text(2.45, 0.32, expression(symbol("\257")*Glu))
  text(3.6, 0.40, expression(symbol("\257")*Ins))
  text(3.3, 0.72, expression(symbol("\255")*tCho))
  text(2.35, 0.50, expression(symbol("\255")*NAAG))
}

agg_tiff("fig_3.tiff", width = 1200, height = 700, scaling = 2)
plot_grid(p1, p2, p3, p4, labels = c('A', 'B', 'C', 'D'), label_size = 14) |>
  print()
dev.off()

rats_snr <- c(as.numeric(rats_gm_snr), as.numeric(rats_wm_snr)) 

mean_rats_snr <- rats_snr |> mean() |> round()

no_rats_snr <- c(as.numeric(no_rats_gm_snr), as.numeric(no_rats_wm_snr))

mean_no_rats_snr <- no_rats_snr |> mean() |> round()

rats_lw <- c(as.numeric(rats_gm_peak_info$fwhm_ppm),
             as.numeric(rats_wm_peak_info$fwhm_ppm))

mean_rats_lw <- rats_lw |> mean() |> round(3)

no_rats_lw <- c(as.numeric(no_rats_gm_peak_info$fwhm_ppm),
                as.numeric(no_rats_wm_peak_info$fwhm_ppm))

mean_no_rats_lw <- no_rats_lw |> mean(na.rm = TRUE) |> round(3)

paste0("no rats mean SNR : ", mean_no_rats_snr) |> cat()
paste0("\nrats mean SNR    : ", mean_rats_snr) |> cat()
paste0("\nno rats mean LW  : ", mean_no_rats_lw, " ppm") |> cat()
paste0("\nrats mean LW     : ", mean_rats_lw, " ppm") |> cat()

t.test(rats_snr, no_rats_snr, paired = TRUE)

t.test(rats_lw, no_rats_lw, paired = TRUE)
