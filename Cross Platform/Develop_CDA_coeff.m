%% Developing CDA Coefficients to Use in LDA classification
%Susan Meerdink
%11/21/2015

%% Import the training & validation spectral libraries into Matlab using ReadSpecLib_wfullmeta.m
%INFO:
% this function reads a spectral library (created in ENVI VIPER tools) so
% that it can be worked with to develop CDA coefficients.

%INPUTS:
%inlibfilebase: a string or string variable name of the spectral library file; must include path name.
directory = 'I:\Classification-Products\Combined\1 - Spectral Library\'; %Set directory

filename = 'f130411_FL02&FL03_AVIRIS_spectral_library'; %Set filename

inlibfilebaseTrainSpec = strcat(directory,filename,'_spectra_calibration_horz'); %set path for training library
inlibfilebaseTrainMeta = strcat(directory,filename,'_metadata_calibration'); %set path for training library
inlibfilebaseValSpec = strcat(directory,filename,'_spectra_validation'); %set path for validation library
inlibfilebaseValMeta = strcat(directory,filename,'_metadata_validation'); %set path for validation library

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
%[outlib_goodbands_Train,all_wl_Train,good_wl_Train,bbl_Train,metadata_Train,n_cols_Train,metadata_fields_Train,nspec_Train] = ReadSpecLib_wfullmeta(inlibfilebaseTrain);
metadata_Train = readtable(strcat(inlibfilebaseTrainMeta,'.csv'));
spectra_Train = readtable(strcat(inlibfilebaseTrainSpec,'.csv'));
outlib_goodbands_Train = cell2mat(table2cell(spectra_Train(:,2:225)));

%open header file and read in all data
% inlibfile_hdr=strcat(inlibfilebaseTrain,'.hdr');
% inlibfile_hdr_info=fileread(inlibfile_hdr);
% %locate & store # spectra and # bands
% nspecline_search='[^\n]*lines[^\n]*';
% nspecline=regexp(inlibfile_hdr_info,nspecline_search,'match');
% nspecline_parse=regexp(nspecline{1},'\ ','split');
% speccells=length(nspecline_parse{1});
% nspectxt=nspecline_parse{speccells};
% nspec=str2num(nspectxt);
% 
% nbandsline_search='[^\n]*samples[^\n]*';
% nbandsline=regexp(inlibfile_hdr_info,nbandsline_search,'match');
% nbandsline_parse=regexp(nbandsline{1},'\ ','split');
% nbandscells=length(nbandsline_parse);
% nbandstxt=nbandsline_parse{nbandscells};
% nbands=str2num(nbandstxt);
% 
% %open spectral library file & read in
% inlibfile=strcat(inlibfilebaseTrain,'.sli');
% inlibfileID=fopen(inlibfile);
% fulllib=fread(inlibfileID,[nbands,nspec],'double=>float32');
% fclose(inlibfileID);
% fulllib=fulllib'; %now it's 1 row per spectrum and 1 col per band
% outlib_goodbands_Train=fulllib;

% read in validation library
%[outlib_goodbands_Val,all_wl_Val,good_wl_Val,bbl_Val,metadata_Val,n_cols_Val,metadata_fields_Val,nspec_Val] = ReadSpecLib_wfullmeta(inlibfilebaseVal);
metadata_Val = readtable(strcat(inlibfilebaseValMeta,'.csv'));
spectra_Val = readtable(strcat(inlibfilebaseValSpec,'.csv'));
outlib_goodbands_Val = cell2mat(table2cell(spectra_Val(:,2:225)));

% %open header file and read in all data
% inlibfile_hdr=strcat(inlibfilebaseVal,'.hdr');
% inlibfile_hdr_info=fileread(inlibfile_hdr);
% %locate & store # spectra and # bands
% nspecline_search='[^\n]*lines[^\n]*';
% nspecline=regexp(inlibfile_hdr_info,nspecline_search,'match');
% nspecline_parse=regexp(nspecline{1},'\ ','split');
% speccells=length(nspecline_parse{1});
% nspectxt=nspecline_parse{speccells};
% nspec=str2num(nspectxt);
% 
% nbandsline_search='[^\n]*samples[^\n]*';
% nbandsline=regexp(inlibfile_hdr_info,nbandsline_search,'match');
% nbandsline_parse=regexp(nbandsline{1},'\ ','split');
% nbandscells=length(nbandsline_parse);
% nbandstxt=nbandsline_parse{nbandscells};
% nbands=str2num(nbandstxt);
% 
% %open spectral library file & read in
% inlibfile=strcat(inlibfilebaseVal,'.sli');
% inlibfileID=fopen(inlibfile);
% fulllib=fread(inlibfileID,[nbands,nspec],'double=>float32');
% fclose(inlibfileID);
% fulllib=fulllib'; %now it's 1 row per spectrum and 1 col per band
% outlib_goodbands_Val=fulllib;
% 
% % split data into chunks delimited by {} signs
% remain = inlibfile_hdr_info;
% c=1;
% while true
%    [str, remain] = strtok(remain, '{}');
%    if isempty(str),  break;  end
%    str_store{c}=str;
%    c=c+1;
% end
% n_cells=size(str_store,2);
% 
% %locate & store bad bands info
% for i=1:n_cells
%     check=strfind(str_store{i},'bbl');
%     if check > 0
%         bbl_col=i+1;
%     end
% end
% s2 = regexp(str_store{bbl_col}, '\,', 'split');
% for c=1:nbands
%     bbl(c)=str2num(s2{c});
% end

%Other Variables:
bbl = {0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0};
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
train_group = cell2mat(table2cell(metadata_Train(:,12))); %Pull out species for grouping - it is column 9 in metadata table
%train_group = table2cell(metadata_Train(:,12)); %Pull out species for grouping - it is column 9 in metadata table
validlib = outlib_goodbands_Val(:,goodbandind);
valid_group = cell2mat(table2cell(metadata_Val(:,12)));%Pull out species for grouping - it is column 9 in metadata table
%valid_group = table2cell(metadata_Val(:,12));%Pull out species for grouping - it is column 9 in metadata table


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

%% Load ENVI Image
imageFile = 'I:\AVIRIS\FL03\7 - Other Files\FL03_f130411t01p00r07rfl_hpc18_bottomhalf_BIL';
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

%% Write ENVI file
outputLoc = 'I:\AVIRIS\FL03\7 - Other Files\FL03_f130411t01p00r07_CDALDA_classification';
headerLoc = 'I:\AVIRIS\FL03\7 - Other Files\FL03_f130411t01p00r07_CDALDA_classification.hdr';
info = enviinfo(outputImage);
%info.data_type = 'ENVI Classification';
%info.class_names = unique(metadata_Train(:,10));
%info.map_info = imageData.info.map_info;
enviwrite(outputImage,info,outputLoc,headerLoc);

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