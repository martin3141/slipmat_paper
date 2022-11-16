# SLIPMAT paper

R code to generate results from:

https://www.biorxiv.org/content/10.1101/2022.11.15.516599v1

Data available from:

https://doi.org/10.5281/zenodo.7189140

Place these scripts in a folder at the root level of the BIDS data folder, eg same level as dataset_description.json.

Note that FSL BET followed by FSL FAST must be run on each sub-??_T1w.nii.gz file for these scripts to work. Following commands were used:

`bet sub-??_T1w.nii.gz T1_brain -B`

`fast T1_brain -B`
