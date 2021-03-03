# Repository for the article ["Distributed chromatic processing at the interface between retina and brain in the larval zebrafish", Guggiana et al. 2021, Curr. Bio.](https://doi.org/10.1016/j.cub.2021.01.088)

## Requirements
- MATLAB (the code was written using version R2018b but will most likely work in newer versions)
- External MATLAB packages (place in Guggiana_2021/Subscripts and Utilities/external_packages , as the code was built to run isolated from the rest of the MATLAB path)
  - [SPASM](https://www.jstatsoft.org/article/view/v084i10)
  - [UMAP](https://de.mathworks.com/matlabcentral/fileexchange/71902-uniform-manifold-approximation-and-projection-umap)

## General instructions
Please follow this readme to be able to reproduce the results in the figures of the publication.

In brief, the steps are as follows

1. Clone/download this repository into a target folder and create the folder structure
2. Download the preprocessed files from the associated Mendeley Data Repository
3. Run the clustering analysis
4. Find the desired analyses from the list below and run the indicated code

## 1) Download/clone the repository
First make sure all external requirements listed above are installed. Click on the clone/download button on the right and choose a target folder to download the code to (alternatively use git to clone the repository). Afterwards, edit the script Paths.m with the folder where the repository was cloned to (assign it as a string value to the variable root_path).

## 2) Download the data
Go to the [Mendeley Data repository](http://dx.doi.org/10.17632/szj869h34m.1) for this publication, download the zip file and uncompress the entire Analysis folder into the main code folder (Guggiana_2021, should be next to Subscripts and Utilities)

As a reference, the ROI datasets present (found in the rois folder) are:

- p8_gc6s: Tectum data for Figure 4
- p8_SynG6s: AF10 data for Figure 4
- p17b_gc6s: RA and Tectum data for Figures 1-3
- p17b_syngc6s: RGC and AF10 data for Figures 1-3
- p17b_syngc6s_conv_tau_0.2: delayed RGC and AF10 data for Figures S1G and S2H
- p17b_h2b6s: RA and Tectum data for Figures S1F and S2F
- p17bdownsample_gc6s: spatially downsampled RA and Tectum data for Figures S1H and S2G

## 3) Cluster the data
Run the script Stage3_Cluster. Select all the datasets. This will automatically cluster the ROIs as in our publication and generate the files that will be used for most of the subsequent analyses. (this step can take a while so heads up)

## 4) Find the desired analysis and follow the instructions
For all scripts, start by running the first cell (cell 0). A pop-up window will prompt selection of the desired datasets. Select the target dataset and then run the cells with the analysis of interest (for example, for Supplementary Figure 1F, one requires the p17b_gc6s and the p17b_h2b6s datasets). Listed below are the figure panels corresponding to each script. Also note that Stage10_registration and Stage13_convolution require special downloads, listed in their instructions.

### Script1_tracePlotting
1. 1F, 4C, S1F-H middle top
2. 2D top
3. 1E
4. 2C
5. S2A-B

### Script2_plotUMAP
1. 3C, S3C-D
2. 4F, S4D

### Script3_plotPCAtraj
1. 3F
2. 3G
3. S3E-F S4C

### Script4b_Classify
- Provided solely as reference, as the required classifiers are provided with the dataset.

### Script4_plotClassifier
1. S2J, S2K
2. 2E
3. 3H
4. 3I
5. 3J, 4G, S1F-H

### Script5_Model
1. S2C

### Script6_features
1. 1H
2. 1I

### Script7_registration
- Run cells 0 AND 1

1. 1G, 4D
2. 2D_right
3. 1J, S1D
4. 1K, 4E


### Script8_correlation
1. 3D, S4A
2. S3A
3. S3B
4. S4B
5. 3E
6. S2D-E

### Script9_convolution
- Download the data from [Zhou et al. 2020](https://datadryad.org/stash/dataset/doi:10.5061/dryad.7sqv9s4pm)
- Place the files in Guggiana_2021/Analysis/reference

1. S2I
2. S1C

### Script10_mixROIs
1. S2F-H
2. S1F-H middle bottom
