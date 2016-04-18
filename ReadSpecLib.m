%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes an input spectral library file name including path %%
% directory and a column name for the classes to be used in the          %%
% corresponding library metadata .csv file and reads in the spectral     %%
% library data including the wavelengths, bad bands list and metadata.   %%
% It applies the bad bands list and returns an array of the spectral     %%
% data [nspec, ngoodbands], a variable called 'textdata' which is the    %%
% metadata, a vector of group number (by class), size [nspec], and a     %%
% cell array vector of the unique class names (index matches the group   %%
% number).                                                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%
%  INPUTS    %
%%%%%%%%%%%%%%
%inlibfilebase: 
% a string or string variable name of the spectral library
% file; must include path name.
%
%class_col_name:
% a string or string variable name of the column name to use for grouping
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%  OUTPUTS   %
%%%%%%%%%%%%%%
%outlib_goodbands: 
% a named variable into which the spectral library data
% is saved.
%
%group_ind:
% a vector of group indices for each spectrum.
%
%ind2classlist:
% a cell array of class names; the indices correspond to the values in
% group_ind.
%
%good_wl:
% a vector of the wavelengths corresponding to the output library
%%%%%%%%%%%%%%


function [outlib_goodbands,group_ind,ind2classlist,good_wl] = ...
            ReadSpecLib(inlibfilebase,class_col_name)
    
%% Read in key info from header (nspec, nbands, bbl, wavelengths, etc.)

%open header file and read in all data
inlibfile_hdr=horzcat(inlibfilebase,'.hdr');
inlibfile_hdr_info=fileread(inlibfile_hdr);


%locate & store # spectra and # bands
nspecline_search='[^\n]*lines[^\n]*';
nspecline=regexp(inlibfile_hdr_info,nspecline_search,'match');
nspecline_parse=regexp(nspecline{1},'\ ','split');
speccells=length(nspecline_parse{1});
nspectxt=nspecline_parse{speccells};
nspec=str2num(nspectxt);
nbandsline_search='[^\n]*samples[^\n]*';
nbandsline=regexp(inlibfile_hdr_info,nbandsline_search,'match');
nbandsline_parse=regexp(nbandsline{1},'\ ','split');
nbandscells=length(nbandsline_parse);
nbandstxt=nbandsline_parse{nbandscells};
nbands=str2num(nbandstxt);


% split data into chunks delimited by {} signs
remain = inlibfile_hdr_info;
c=1;
while true
   [str, remain] = strtok(remain, '{}');
   if isempty(str),  break;  end
   str_store{c}=str;
   c=c+1;
end
n_cells=size(str_store,2);

%locate & store bad bands info
for i=1:n_cells
    check=strfind(str_store{i},'bbl');
    if check > 0
        bbl_col=i+1;
    end
end
s2 = regexp(str_store{bbl_col}, '\,', 'split');
for c=1:nbands
    bbl(c)=str2num(s2{c});
end
goodbandind=find(bbl);

%locate & store wavelength info
for i=1:n_cells
    check=strfind(str_store{i},'wavelength');
    if check > 0
        wl_col=i+1;
    end
end
s3 = regexp(str_store{wl_col}, '\,', 'split');
for c=1:nbands
    wl(c)=str2num(s3{c});
end
good_wl=wl(goodbandind);

%% Read in Spectral Library & Metadata

%open spectral library file & read in
inlibfileID=fopen(horzcat(inlibfilebase,'.sli'));
fulllib=fread(inlibfileID,[nbands,nspec],'float32=>float32');
fclose(inlibfileID);
fulllib=fulllib'; %now it's 1 row per spectrum and 1 col per band
outlib_goodbands=fulllib(:,goodbandind);


% Import the metadata 
inlibfile_metadata=horzcat(inlibfilebase,'.csv');
inlib_metadata_all=importdata(inlibfile_metadata);
inlib_metadata=regexp(inlib_metadata_all, '\,', 'split');
n_cols=size(inlib_metadata{1});
n_cols=n_cols(2);
classcol_name=class_col_name;
curr_match=strfind(inlib_metadata{1},classcol_name);
for c=1:n_cols
    if curr_match{c} > 0
        class_col=c;
    end
end

%determine classes and create group index
for x=2:nspec+1
    spec_class_list(x-1)=inlib_metadata{x}(class_col);
end
group_ind=grp2idx(spec_class_list);
ind2classlist=unique(spec_class_list);

return;