# foodwebs_vs_land_use
R Scripts to reproduce results of "Signatures of land use intensity on european tetrapod food-web architectures"

First download the dataset here: https://zenodo.org/record/5831144

To reproduce the data preprocessing from the raw data, run script **Preproc_and_make_nets.R**, which loads the Rdata file **raw_data**.

To reproduce manuscript results / Figures, run script **analysis_pipeline.R**, which loads the files **preprocessed_data** and **TrophicNetworksList**.
