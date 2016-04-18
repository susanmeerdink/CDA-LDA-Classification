%This script should be run after CDA_manova to write the canonical
%coefficients to a .csv (ngoodbands rows by n_vars columns) file.

function writeCDAvars(cdastats,outfilebase,classes)

% Find number of unique classes and calculate the number of function to
% use. (ngroups-1).
n_vars=length(classes)-1;

% Extract the function variables from cdastats and write them to .csv
canoncoeffs=cdastats.eigenvec(:,1:n_vars);
coeff_file=horzcat(outfilebase,'_CDAcoeffs.csv');
csvwrite(coeff_file,canoncoeffs);

% %Create an array with the canonical training variables and an added column
% %for their group #
% trainvars_file=horzcat(outfilebase,'_CDAtrainvars.csv');
% 
% csvwrite(trainvars_file,CDAtrainvars);


%must just output a file with training vars and group membership (ngroups-1
%+ a group column)
%another file with CDA vals
return
