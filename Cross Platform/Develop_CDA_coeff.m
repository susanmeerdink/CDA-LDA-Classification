%% Developing CDA Coefficients to Use in LDA classification
%Susan Meerdink
%11/21/2015

%% Import the training & validation spectral libraries into Matlab using ReadSpecLib_wfullmeta.m
%INFO:
% this function reads a spectral library (created in ENVI VIPER tools) so
% that it can be worked with to develop CDA coefficients.

%INPUTS:
%inlibfilebase: a string or string variable name of the spectral library file; must include path name.
directory = 'F:\Classification-Products\1 - Spectral Library\Combined Single Date\'; %Set directory

filename = '140416_spectral_library'; %Set filename

inlibfilebaseTrainSpec = strcat(directory,filename,'_calibration_spectra'); %set path for training library
inlibfilebaseTrainMeta = strcat(directory,filename,'_calibration_metadata'); %set path for training library
inlibfilebaseValSpec = strcat(directory,filename,'_validation_spectra'); %set path for validation library
inlibfilebaseValMeta = strcat(directory,filename,'_validation_metadata'); %set path for validation library

% read in training library
metadata_Train = readtable(strcat(inlibfilebaseTrainMeta,'.csv'));
spectra_Train = readtable(strcat(inlibfilebaseTrainSpec,'.csv'),'ReadVariableNames',0);
outlib_goodbands_Train = cell2mat(table2cell(spectra_Train(2:end,6:229)));

% read in validation library
metadata_Val = readtable(strcat(inlibfilebaseValMeta,'.csv'));
spectra_Val = readtable(strcat(inlibfilebaseValSpec,'.csv'),'ReadVariableNames',0);
outlib_goodbands_Val = cell2mat(table2cell(spectra_Val(2:end,6:229)));

%Other Variables:
bbl = [0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0];
goodbandind = find(bbl);
%% Run CDA classification using CDA_manova.m
% INFO:
% This step creates the CDA coefficients and does accuracy assessment

% INPUTS:
% trainlib: Your training library read in using ReadSpecLib_wfullmeta.m
% train_group: pull out the group from the metadata (species, etc)
% validlib: Your validation library read in using ReadSpecLib_wfullmeta.m
% valid_group: pull out the group from the metadata (species, etc)
trainlib = outlib_goodbands_Train(:,goodbandind);
train_group = cell2mat(table2cell(metadata_Train(:,15))); %Pull out species for grouping 

validlib = outlib_goodbands_Val(:,goodbandind);
valid_group = cell2mat(table2cell(metadata_Val(:,15)));%Pull out species for grouping 

% OUTPUTS:
% accstats: accuracy stats (overall, cappa, error matrix, producers/users, valid class(what it was classified as), valid group (truth))
% cdastats: comes out of manova1 function, eigenvectors, eigen values, all cda outputs
% canon_vars_Train: cda verison of your spectral library, what you are going to use to classify images, key output of CDA

% FUNCTION:
[accstats,cdastats,canon_vars_Train] = CDA_manova(trainlib,train_group,validlib,valid_group);

% PLOT:
input = validlib * cdastats.eigenvec(:,1:23); 
valid_class = classify(input,canon_vars_Train,train_group);
scatter(valid_group, valid_class)

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

%% Load ENVI Image
imageFile = 'I:\AVIRIS\FL03\7 - Other Files\FL03_f131125t01p00r14rfl_hpc18_bottomhalf_BIL';
imageData = enviread(imageFile);
disp('Done Reading in ENVI Image')
%% Run Linear Discriminant Analysis Classification
outputImage = zeros(length(imageData.y),length(imageData.x));
for r = 1:length(imageData.y) %loop through rows/lines
    temp = single(reshape(imageData.z(r,:,:),[length(imageData.x),224])); %grab row of data
    temp = temp(:,goodbandind);
    tempVars = temp* cdastats.eigenvec(:,1:18);
    valid_class=classify(tempVars,canon_vars_Train,train_group);
    outputImage(r,:) = valid_class';
end
disp('Done Classifying Image')
%% Write ENVI file
outputLoc = 'I:\AVIRIS\FL03\7 - Other Files\FL03_f131125t01p00r14_CDALDA_classification';
headerLoc = 'I:\AVIRIS\FL03\7 - Other Files\FL03_f131125t01p00r14_CDALDA_classification.hdr';
info = enviinfo(outputImage);
%info.data_type = 'ENVI Classification';
%info.class_names = unique(metadata_Train(:,10));
%info.map_info = imageData.info.map_info;
enviwrite(outputImage,info,outputLoc,headerLoc);