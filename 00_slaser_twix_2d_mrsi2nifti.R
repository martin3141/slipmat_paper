library(spant)

# remove old variables
rm(list = ls())

file_paths     <- Sys.glob("../sub-??/mrs/*.dat")
file_paths_out <- gsub(".dat", ".nii.gz", file_paths)

for (n in 1:length(file_paths)) {
  
  paste0("Processing ", n, " of ", length(file_paths), " : ", file_paths[n], 
         "\n") |> cat()
  
  twix_raw <- read_mrs(file_paths[n], verbose = TRUE)
  
  twix_proc <- twix_raw |> decimate_mrs_fd() |> recon_twix_2d_mrsi() |> 
               comb_coils()
  
  write_mrs(twix_proc, file_paths_out[n], force = TRUE)
}
