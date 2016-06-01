%% Developing CDA Coefficients to Use in LDA classification
%Susan Meerdink
%11/21/2015

%% Import the training & validation spectral libraries into Matlab using ReadSpecLib_wfullmeta.m
%INFO:
% this function reads a spectral library (created in ENVI VIPER tools) so
% that it can be worked with to develop CDA coefficients.

%INPUTS:
%inlibfilebase: a string or string variable name of the spectral library file; must include path name.
directory = 'R:\users\susan.meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS\'; %Set directory
%directory = 'R:\users\susan.meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS & MASTER\'; %Set directory
%directory = 'R:\users\susan.meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\Combined\'; %Set directory

filename = '20130411_Spectral_Library_AVIRIS_sorted'; %Set filename
%filename = '20130606_Spectral_Library_AVIRIS_sorted'; %Set filename
%filename = '20131125_Spectral_Library_AVIRIS_sorted'; %Set filename
%filename = '20140416_Spectral_Library_AVIRIS_sorted'; %Set filename
%filename = '20140606_Spectral_Library_AVIRIS_sorted'; %Set filename
%filename = '20140829_Spectral_Library_AVIRIS_sorted'; %Set filename
%filename = '20130411_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
%filename = '20130606_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
%filename = '20131125_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
%filename = '20140416_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
%filename = '20140606_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
%filename = '20140829_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
%filename = '2013&2014_Spectral_Library_AVIRIS'; %Set filename
%filename = '2013&2014_Spectral_Library_AVIRIS&MASTER'; %Set filename

inlibfilebaseTrain = strcat(directory,filename,'_train'); %set path for training library
inlibfilebaseVal = strcat(directory,filename,'_valid'); %set path for validation library

% OUTPUTS:
% _Train: all variables associated with the training library
% _Val: all variables associated with the validation librar
% outlib_goodbands: a named variable into which the spectral library data is saved.
% all_wl: vector of All wavelengths
% good_wl: Only the  good wavelengths corresponding to the output library
% bbl: Bad Band List, used to select good wavelengths
% metadata: contains the metadata for the spectral library
% n_cols: Number of columns
% metadata_fields: The fields that are inthe metadata 
% nspec: The total number of spectra in the library

% FUNCTION:
% read in training library
[outlib_goodbands_Train,all_wl_Train,good_wl_Train,bbl_Train,metadata_Train,n_cols_Train,metadata_fields_Train,nspec_Train] = ReadSpecLib_wfullmeta(inlibfilebaseTrain);
% read in validation library
[outlib_goodbands_Val,all_wl_Val,good_wl_Val,bbl_Val,metadata_Val,n_cols_Val,metadata_fields_Val,nspec_Val] = ReadSpecLib_wfullmeta(inlibfilebaseVal);

%% Run CDA classification using CDA_manova.m
% INFO:
% This step creates the CDA coefficients and does accuracy assessment

% INPUTS:
% trainlib: Your training library read in using ReadSpecLib_wfullmeta.m
% train_group: pull out the group from the metadata (species, etc)
% validlib: Your validation library read in using ReadSpecLib_wfullmeta.m
% valid_group: pull out the group from the metadata (species, etc)
trainlib = outlib_goodbands_Train;
train_group = table2cell(metadata_Train(:,10)); %Pull out species for grouping - it is column 9 in metadata table
validlib = outlib_goodbands_Val;
valid_group = table2cell(metadata_Val(:,10));%Pull out species for grouping - it is column 9 in metadata table

% OUTPUTS:
% accstats: accuracy stats (overall, cappa, error matrix, producers/users, valid class(what it was classified as), valid group (truth))
% cdastats: comes out of manova1 function, eigenvectors, eigen values, all cda outputs
% canon_vars_Train: cda verison of your spectral library, what you are going to use to classify images, key output of CDA

% FUNCTION:
[accstats,cdastats,canon_vars_Train] = CDA_manova(trainlib,train_group,validlib,valid_group);

%% Run writeCDAvars.m
% This step creates a .csv file with the CDA coefficients
% doesn't return anything in matlab
% all it does is put coefficients in csv so that it is ready for IDL

% INPUTS:
% cdastats: This is the output from the CDA_manova.m function
% outfilebase: This is the location of the output spectral library
% classes: This is the number of classes (or groups) from your spectral
% library. The function will calculate the number, just give list of class
% names

outfilebase = strcat(directory,filename,'_CDAvars');
classes = cdastats.gnames; %Set classes equal to the number of group names

writeCDAvars(cdastats,outfilebase,classes)

%% Run Linear Discriminant Analysis Classification
%%%%%%%%AFTER IDL STUFF - NOT DONE EITHER %%%%%%%%%%%%%%%%%%%
% % INPUTS  
% % directory = file directory pathname
% directory = 'E:\Meerdink\Dropbox\Code\CDA_Work_With_Keely\Data\';
% % fileTrainLib = file name of the training library
% fileTrainLib = 'SB_alldates_20092011_train';%NO .sli
% % fileImage = file name of image to be classified
% fileImage ='f130411t01p00r08rfl_hpc18_mask';
% % className = class column name
% className = 'Dominant';
% 
% image_classify_w_lda(directory,fileTrainLib,fileImage,className)