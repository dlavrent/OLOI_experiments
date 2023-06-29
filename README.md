# OLOI_experiments/
> Matthew A. Churgin, Danylo Lavrentovich
>
> de Bivort Lab

Experimental data and analysis scripts for *Neural correlates of individual odor preference in **Drosophila***.

This GitHub repository contains raw behavioral and 2-photon calcium data, along with MATLAB scripts that generate the figures. Mappings between MATLAB scripts and figure panels are listed at the above link. A complete repository, with large raw microscopy images too large for sharing via GitHub, is available by contacting the authors.

### Directory structure

A brief walkthrough of the directories in this repository:

**utilities/** contains helpful auxiliary functions and .mat files used throughout the repository. If you download this repository, please set your local path in `SET_DIRECTORY_PATH.m` and run it to generate `analysis_dir_path.mat`, which is loaded by all MATLAB analysis scripts to keep file paths consistent

**ORNvsPN_analysis_ALLDATA/**, **ORN_analysis_oct_vs_air_choice/**, **ORN_analysis/**, **PN_analysis_oct_vs_air_choice/**, **PN_analysis/** store calcium imaging data and perform analyses. Main analysis script that generates most of the figures is `ORNvsPN_analysis_ALLDATA/plotRandomSubsetPNandORN.m`

**IHC/** contains Brp-Short immunohistochemistry analysis

**rho_signal_inference/** contains a Python script that runs simulations to infer the extent of latent neural predictor - latent behavior correlation based on the calcium/Brp-Short linear models


