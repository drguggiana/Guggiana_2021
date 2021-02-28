# Repository for the article ["Distributed chromatic processing at the interface between retina and brain in the larval zebrafish", Guggiana et al. 2021, Curr. Bio.](https://doi.org/10.1016/j.cub.2021.01.088)

## Requirements
- MATLAB (the code was written using version R2018b but will most likely work in newer versions)
- External MATLAB packages
  - uipickfiles

## General instructions
Please follow this readme to be able to reproduce the results in the figures of the publication.

In brief, the steps are as follows

1. Clone/download this repository into a target folder and create the folder structure
2. Download the preprocessed ROIs from the [corresponding repository](http://dx.doi.org/10.17632/szj869h34m.1) in Mendeley Data
3. Run the clustering analysis
4. Find the desired analyses from the list below and run the indicated code

## 1) Download the repository
Click on the clone/download button on the right and choose a target folder to download the code to. Afterwards, edit the script Paths.m with the folder where the repository was cloned to (assign it as a string value to the variable root_path). Run the script so the basic folder structure is created. The datasets present are:

- p8_gc6s: Tectum data for Figure 4
- p8_SynG6s: AF10 data for Figure 4
- p17b_gc6s: RA and Tectum data for Figures 1-3
- p17b_syngc6s: RGC and AF10 data for Figures 1-3
- p17b_syngc6s_conv_tau_0.2: delayed RGC and AF10 data for Figures S1G and S2H
- p17b_h2b6s: RA and Tectum data for Figures S1F and S2F
- p17bdownsample_gc6s: spatially downsampled RA and Tectum data for Figures S1H and S2G

## 2) Download the data
Go to the Mendeley Data repository for this publication and download the preprocessed ROIs to the folder Stage2.

## 3) Cluster the data
Run the script Stage3_Cluster. Select all the datasets. This will automatically cluster the ROIs as in our publication and generate the files that will be used for most of the subsequent analyses.

## 4) Find the desired analysis and follow the instructions
For all scripts, when running them, a pop-up window will prompt selection of the desired datasets. Select the target dataset and then run the cells with the analysis of interest. Listed below are the figure panels corresponding to each script.

### Stage4_tracePlotting
1. 1F, 4C
2. 2D top
3. 1E
4. 2C
5. S2A-B

### Stage5_plotUMAP
1. 3C, S3C-D
2. 4F, S4D

### Stage6_plotPCAtraj
1. 3F
2. 3G
3. S3E-F S4C

### Stage8_Model
1. S2C

### Stage9_features
1. 1H
2. 1I

### Stage11_correlation
1. 3D, S4A
2. S3A
3. S3B
4. S4B
5. 3E
6. S2D-E
