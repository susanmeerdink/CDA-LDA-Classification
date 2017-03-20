# CDA-LDA-Classification

Project related to CDA and LDA classification (cross platform)

### Dependencies
---

Classification Part 1 requires ENVI's IMSL license, which was on Raider's 32-bit ENVI + IDL
(currently on ariel will only run 8 mins - will be moved to Raider)

Need to run IDL code on IDL 64 bit without ENVI

### Code Order
---

#### Building spectral libraries:

  Use ENVI and VIPER Tools 2.0 to build library (see steps below)

1. `fillingOutENVISpectralLibraryMetadata.py`
2. `strip_redundant_spectra.pro`
3. `sample_spectral_lib.pro`

#### Developing CDA coefficients:
  * `Develop_CDA_coeff.m` (this calls three functions listed below)
    * `ReadSpecLib_wfullmeta.m`
    * `CDA_manova.m`
    * `writeCDAvars.m`

#### Classification:
  * `create_cda_image_wo_control_file.pro`
  * `image_classify_w_lda_wo_control_file_Part1.pro`
  * `image_classify_w_lda_wo_control_file_Part2.pro`

#### Accuracy Assessment
  * Use ENVI to Extract ROI information from classified image
    * `pullingOutROIInfo.py`
    * `fillingOutENVISpectralLibraryMetadata.py`
    * `pullingOutROIval.py`
    * `calcclassificationresults.py`
    * `LoopFilesPullClassificationStats.py`
    * `LoopFilesPullCDAVars.py`


### DESCRIPTION OF CODE
---

##### Building spectral libraries:

1. Load image that spectra will be selected from in ENVI + IDL.
2. Load Vector of polygons -> Vector -> Open Vector File
3. Export vector data into ROI
     1. File > Export Layers to ROI
     2. Select Data File to Associate with new ROIs
     3. Export EVF Layers to ROI -> Convert each record of an EVF layer to new ROI -> Polygon ID
4. Open ROI Tool (to check that ROIs have transferred)
     1. Image Display -> Tools -> Region of Interest -> ROI Tool
5. Open VIPER Tools 2.0 (Spectral -> VIPER Tools)
6. Create/Manage Spectral Library -> Build Library from ROIs
     * Note you must have image with ROIs pulled up in a display
7. Create Metadata for spectral library
8. Use python code `fillingOutENVISpectralLibraryMetadata.py` to fill out metadata 
9. Use `strip_redundant_spectra.pro` to remove duplicate spectra 
10. Use `sample_spectral_lib.pro` with a control file to create training and validation libraries from an input spectral library using user-defined thresholds for sampling from each dominant species randomly. Two libraries with associated metadata are output, one training data set and one validation data set.

#### Developing CDA Coefficients:

1. `Develop_CDA_coeff.m`: this calls three functions listed below
  1. `ReadSpecLib_wfullmeta.m`: this function reads a spectral library (created in ENVI VIPER tools) so that it can be worked with to develop CDA coefficients.
  2. `CDA_manova.m`: This function creates the CDA coefficients and does accuracy assessment
  3. `writeCDAvars.m`: This function creates a .CSV file with the CDA coefficients

#### Classification:

1. `create_cda_lib_wo_control_file.pro`: This step creates the CDA training spectral library
2. `create_cda_image_wo_control_file.pro`: This program takes a set of CDA coefficients and a reflectance image and creates an output image of the CDA variables.
3. `image_classify_w_lda_wo_control_file_Part1.pro`: This step reads in the CDA training spectral library and runs LDA classification and outputs library to CSV
4. `image_classify_w_lda_wo_control_file_Part2.pro`: This step classifies the CDA image created previously and classifies the image using the LDA classification library created in step above.

#### Accuracy Assessment

1. Use ENVI to Extract ROI information from classified image
2. `pullingOutROIInfo.py`: Reads through ENVI ASCII file (output from ROI extraction) and puts into .CSV for easy reading
3. Use python code `fillingOutENVISpectralLibraryMetadata.py` to fill out metadata 
4. `pullingOutROIval.py`: Reads through .CSV of ROI information and pulls out rows that are to be used in validation
5. `calcclassificationresults.py`: This code reads in ROI formated CSV (formatted using `pullingOutROIVal.py`) and calculates classification results including producers, users, overall, and kappa accuracy results are output to a CSV
6. OPTIONAL: If you've run multiple CDA-LDA classifications and want to compare model results use `LoopFilesPullClassificationStats.py` which reads through .CSVs and pulls out the producers, users, overall, and kappa accuracy results are output to a CSV
7. OPTIONAL: If you've run multiple CDA-LDA classifications and want to compare CDA coefficients use `LoopFilesPullCDAVars.py` which reads through .CSVs and pulls out CDA variables.


