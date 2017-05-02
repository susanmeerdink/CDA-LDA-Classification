%% Developing CDA Coefficients to Use in LDA classification
% Using this code to trouble shoot CDA coefficient development
% Susan Meerdink
% 4/26/17

%% Import the training & validation spectral libraries into Matlab using ReadSpecLib_wfullmeta.m

%INPUTS:
%inlibfilebase: a string or string variable name of the spectral library file; must include path name.
directory = 'F:\Classification-Products\1 - Spectral Library\Combined Single Date\'; %Set directory

filename = '140416_spectral_library'; %Set filename

inlibfilebaseTrainSpec = strcat(directory,filename,'_calibration_spectra'); %set path for training library
inlibfilebaseTrainMeta = strcat(directory,filename,'_calibration_metadata'); %set path for training library
inlibfilebaseValSpec = strcat(directory,filename,'_validation_spectra'); %set path for validation library
inlibfilebaseValMeta = strcat(directory,filename,'_validation_metadata'); %set path for validation library
inlibfilebaseSpec = strcat(directory,filename,'_spectra'); %set path for library (all samples)
inlibfilebaseMeta = strcat(directory,filename,'_metadata'); %set path for library (all samples)

% read in all samples library
metadata_All = readtable(strcat(inlibfilebaseMeta,'.csv'));
spectra_All = readtable(strcat(inlibfilebaseSpec,'.csv'),'ReadVariableNames',0);
outlib_goodbands_All = cell2mat(table2cell(spectra_All(2:end,6:229)));

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

%% Plot Spectra
uniqueGroup = unique(cell2mat(table2cell(metadata_All(:, 15))));
for g = 1: size(uniqueGroup,1)
    
    index = find(cell2mat(table2cell(metadata_All(:,15))) == uniqueGroup(g));
    spectra = cell2mat(table2cell(spectra_All(index,6:end)));
    meanSpec = mean(spectra);
    minSpec = min(spectra);
    maxSpec = max(spectra);
    figure()
    hold on
    plot(wavelength, meanSpec)
    plot(wavelength, minSpec)
    plot(wavelength, maxSpec)
    hold off
end

%% Plot Polygons
uniquePoly = unique(table2cell(metadata_All(:, 3)));
for p = 285: 338%size(uniquePoly,1)
    figure()
    index = find(strcmp(uniquePoly(p), table2cell(metadata_All(:,3))));
    spectra = cell2mat(table2cell(spectra_All(index,6:end)));
    plot(wavelength, spectra)
    t = title(strcat(uniquePoly(p),' (n= ', num2str(size(index,1)), ')'));
    set(t,'Interpreter','none');
end

%% Plot Dominant with mean from polygon
uniqueGroup = unique(table2cell(metadata_All(:, 15)));
domNum = cell(size(uniqueGroup,1),3);
domNum(:,1) = uniqueGroup;
for g = 1: size(uniqueGroup,1)
    numPixel = 0;
    indexGroup = find(strcmp(uniqueGroup(g), table2cell(metadata_All(:,15))));
    meta = metadata_All(indexGroup, :);
    spectra = cell2mat(table2cell(spectra_All(indexGroup, 6:end))); 
    uniquePoly = unique(table2cell(meta(:, 3)));
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    for p = 1: size(uniquePoly,1)
        numPixel = numPixel + size(index,1);
        index = find(strcmp(uniquePoly(p), table2cell(meta(:,3))));
        meanSpec = mean(spectra(index,1:end));  
        plot(wavelength,meanSpec)
    end
    title(uniqueGroup(g))
    legend(uniquePoly, 'Location','BestOutside')
    hold off
    
    domNum(g,2) = num2cell(size(uniquePoly,1));
    domNum(g,3) = num2cell(numPixel);
end

%% Run CDA classification using CDA_manova.m
% INFO:
% This step creates the CDA coefficients and does accuracy assessment

% INPUTS:
% trainlib: Your training library read in using ReadSpecLib_wfullmeta.m
% train_group: pull out the group from the metadata (species, etc)
% validlib: Your validation library read in using ReadSpecLib_wfullmeta.m
% valid_group: pull out the group from the metadata (species, etc)
trainlib = outlib_goodbands_Train(:,goodbandind);
train_group = cell2mat(table2cell(metadata_Train(:,16))); %Pull out species for grouping 

validlib = outlib_goodbands_Val(:,goodbandind);
valid_group = cell2mat(table2cell(metadata_Val(:,16)));%Pull out species for grouping 

% OUTPUTS:
% accstats: accuracy stats (overall, cappa, error matrix, producers/users, valid class(what it was classified as), valid group (truth))
% cdastats: comes out of manova1 function, eigenvectors, eigen values, all cda outputs
% canon_vars_Train: cda verison of your spectral library, what you are going to use to classify images, key output of CDA

% FUNCTION:
[accstats,cdastats,canon_vars_Train] = CDA_manova(trainlib,train_group,validlib,valid_group);

% PLOT:
input = validlib * cdastats.eigenvec(:,1:22); 
valid_class = classify(input,canon_vars_Train,train_group);
scatter(valid_group, valid_class)
input = trainlib * cdastats.eigenvec(:,1:22); 
train_class = classify(input,canon_vars_Train,train_group);
scatter(train_group,train_class)
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