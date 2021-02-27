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
Click on the clone/download button on the right and choose a target folder to download the code to. Afterwards, edit the script Paths.m with the folder where the repository was cloned to (assign it as a string value to the variable root_path). Run the script so the basic folder structure is created.

## 2) Download the data
Go to the Mendeley Data repository for this publication and download the preprocessed ROIs to the folder Stage2.

## 3) Cluster the data
Run the script Stage3_Cluster. Select all the datasets. This will automatically cluster the ROIs as in our publication and generate the files that will be used for most of the subsequent analyses.

## 4) Find the desired analysis and follow the instructions
