# SLIPMAT paper

R code to generate results from:

https://www.biorxiv.org/content/10.1101/2022.11.15.516599v1

Data available from:

https://doi.org/10.5281/zenodo.7189139

Place these scripts in a folder at the root level of the BIDS data folder, eg same level as dataset_description.json.

Note that FSL BET followed by FSL FAST must be run on each sub-??_T1w.nii.gz file for these scripts to work. Following commands were used for each T1 dataset:

`bet sub-??_T1w.nii.gz T1_brain -B`

`fast T1_brain`

## Paper citation

Olivia Vella, Andrew P. Bagshaw, Martin Wilson,
SLIPMAT: A pipeline for extracting tissue-specific spectral profiles from 1H MR spectroscopic imaging data,
NeuroImage,
Volume 277,
2023,
120235,
ISSN 1053-8119,
https://doi.org/10.1016/j.neuroimage.2023.120235
