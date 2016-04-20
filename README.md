# CDA-LDA-Classification
Project related to CDA and LDA classification (cross platform)

**Dependencies**
image_classify_w_lda.pro requires ENVI's IMSL license, which was on Raider's 32-bit ENVI + IDL
(currently on ariel will only run 8 mins - will be moved to Raider)
Need to run IDL code on IDL 64 bit without ENVI

**Code Order**
In MATLAB:
Develop_CDA_coeff.m (this calls three functions)
  ReadSpecLib_wfullmeta.m
  CDA_manova.m
  writeCDAvars.m
In IDL:
create_cda_image_wo_control_file.pro
image_classify_w_lda_wo_control_file_Part1.pro
image_classify_w_lda_wo_control_file_Part2.pro

**INPUTS**
This project requires a spectral library built in ENVI using VIPER Tools 2.0.

**DESCRIPTION OF CODE**
Develop_CDA_coeff.m: this calls three functions listed below
ReadSpecLib_wfullmeta.m: this function reads a spectral library (created in ENVI VIPER tools) so that it can be worked with to develop CDA coefficients.
CDA_manova.m: This function creates the CDA coefficients and does accuracy assessment
writeCDAvars.m: This function creates a .csv file with the CDA coefficients
create_cda_lib_wo_control_file.pro: This step creates the CDA training spectral library
create_cda_image_wo_control_file.pro: This program takes a set of CDA coefficients and a reflectance image and creates an output image of the CDA variables.
image_classify_w_lda_wo_control_file_Part1.pro: This step reads in the CDA training spectral library and runs LDA classification and outputs library to csv
image_classify_w_lda_wo_control_file_Part2.pro This step classifies the CDA image created previously and classifies the image using the LDA classification library created in step above.
