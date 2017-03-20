pro image_classify_w_lda_wo_control_file_Part1

;; This program classifies an image using a training spectral library.
;; It runs using a control file with the inputs specified below. The program uses IDL's discriminant anlaysis (run on Raider, ENVIIDL 4.7 32 bit)
;; function to determine LDA coefficients. The image is read in line by line, and the LDA coefficients are applied.
;; These values yield an LDA score for each class. Each pixel in the image is assigned to the class with the highest 
;; LDA score. The program outputs the classfication image.

; Last updated Phil Dennison to fix minor errors based on older code 15 Mar 2013.
; Last updated Keely Roth to add distance measure for computing class confidences UNDERWAY

;This code no longer uses control file, use hard input below
;Control File Contains
; Line 1 = file directory pathname
; Line 2 = file name of the training library
; Line 3 = file name of image to be classified
; Line 4 = class column name

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Begin Main Program:
close, /ALL
;cfile=envi_pickfile(title='Select a Control file', filter='*.ctl')

;declare inputs
;;;;;;;;;;;;;;;Declare Directory for Training Library;;;;;;;;;;;;;;;;;;;;;;;
pnameTrain = strarr(1)
pnameTrain = 'I:\Classification-Products\FL03\2 - CDA Variables\'; %Set directory
;pnameTrain = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS & MASTER\'; %Set directory
;pnameTrain = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\Combined\'; %Set directory

;;;;;;;;;;;;;;;Declare file name of the training library;;;;;;;;;;;;;;;;;;;;;;;
trainlib_file = strarr(1)
trainlib_file =  'f140829_AVIRIS_spectral_library'; %Set filename
trainlib_file = trainlib_file + '_calibration_CDA.sli' ;Variable for CDA coefficients name

;;;;;;;NEED FOR Part 2 NOT Part 1;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;Declare Directory for Image;;;;;;;;;;;;;;;;;;;;;;;
;pnameImage = strarr(1)
;pnameImage = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Images\Mosaic\'; %Set directory
;
;;;;;;;;;;;;;;;;Declare file name of image to be classified;;;;;;;;;;;;;;;;;;;;;;;
;image_file=strarr(1)
;image_file = '20130411_SBFrontRange_Mosaic'; %Set filename
;;image_file = '20130606_SBFrontRange_Mosaic'; %Set filename
;;image_file = '20131125_SBFrontRange_Mosaic'; %Set filename
;;image_file = '20140416_SBFrontRange_Mosaic'; %Set filename
;;image_file = '20140606_SBFrontRange_Mosaic'; %Set filename
;;image_file = '20140829_SBFrontRange_Mosaic'; %Set filename;
;;image_file = '20130411_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;;image_file = '20130606_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;;image_file = '20131125_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;;image_file = '20140416_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;;image_file = '20140606_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;;image_file = '20140829_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = image_file + '_CDAvars'


;read in control file information
;openr,lun,cfile,/GET_LUN
;readf,lun,pname,trainlib_file,image_file,key_col_name
key_col_name = strarr(1)
key_col_name = 'Dominant' ;Declare class column Name
trainlib_file = pnameTrain + trainlib_file
;image_file = pnameImage + image_file
;CDAcoeff_file=pname+CDAcoeff_file


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; CREATE CLASSIFICATION VARIABLES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;read in training spectral library
envi_open_file,trainlib_file,r_fid=train_fid
nspec=intarr(1)
nb=intarr(1)
envi_file_query,train_fid,spec_names=specnames,nl=nspec,ns=nb,bbl=bbl,wl=wl,dims=tldims
;read in training library
trainlib=fltarr(nb,nspec)

trainlib=envi_get_data(fid=train_fid, pos=0, dims=tldims)
;check for bad bands list, and if there is one,subset spectra to just good bands

if size(bbl,/N_ELEMENTS) ne -1 then begin ;if there is a bad band list then...
  trainlib_goodbands = trainlib[where(bbl eq 1,count),*] 
endif else begin ;if there isn't a bad band list
  bbl = make_array(nb,/BYTE,value=1)
endelse
    

;open training library metadata & read in lines
;root=strsplit(trainlib_file,".",/extract) ;Removed this line of code because there is a period in all users names on Ophelia
root = STRMID(trainlib_file,0, STRLEN(trainlib_file)-4)
metadata_file=root[0]+'.csv'
metadata_lines=strarr(nspec)
metadata_header=strarr(1)
get_lun, u2
openr, u2, metadata_file
readf, u2, metadata_header
readf, u2, metadata_lines
close, u2
free_lun, u2

;find number of metadata columns & key column index
col_headings = strsplit(metadata_header, ',', /extract)
ncols_met=n_elements(col_headings)
for x=0,ncols_met-1 do col_headings(x)=strtrim(col_headings(x),2)
key_col_search=strmatch(col_headings,key_col_name,/FOLD_CASE)
key_col=where(key_col_search eq 1)

;read in spectra class from metadata
specclass=strarr(nspec)
for i = 0L, nspec-1 do begin
  temp = strsplit(metadata_lines(i), ',', /extract)
  specclass(i) = temp(key_col)
  specclass(i)=strtrim(specclass(i),2)
  print,specclass(i)
endfor

;Create a numeric index for the class & determine number of classes
classlist=specclass[uniq(specclass)]          ;string vector of unique classes in the training data
;classlist = ['ADFA', 'ARCA-SALE', 'ARGL', 'BAPI', 'BRNI', 'CECU', 'CEME', 'CESP' ,'CISP', 'ERFA', 'EUSP', 'IRGR', 'MAGF' ,'MARSH' ,'PEAM' ,'PISA' ,'PLRA', 'QUAG', 'QUDO', 'ROCK' ,'SOIL', 'UMCA' ,'AGRES', 'URBAN']
classlist = ['ADFA','AGRES','ARCA-SALE','ARGL','BAPI','BRNI','CECU','CEME','CESP','CISP','ERFA','MAGF','PISA','PLRA','QUDO','ROCK','SOIL','UMCA']
n_groups=n_elements(classlist)
classlist_ind=indgen(n_groups)     ;numeric index for each class
;specclass_ind = fltarr(nspec,1)
specclass_ind = make_array(nspec,1,VALUE = 0,/INTEGER)

;Issues with assigning the group value above so doing code differently
for e = 0,size(specclass,/N_ELEMENTS)-1 do begin ;Loop through the list of species row by row
  for c = 0,n_groups-1 do begin ;Loop through the class list to find a match
    if strmatch(specclass[e],classlist[c]) EQ 1 then begin ;If the strings match
      specclass_ind[e] = classlist_ind[c]+1 ;Then assign that class to the index
      break
    endif 
  endfor
endfor

trainspec = fltarr(nspec,nb+1) ; create an array that will contain both the training spectra and a last column that contains the class number
trainspec(*,0:nb-1) = transpose(trainlib) ; put in the training spectra
;train_canon_vars = fltarr(nspec,nb+1)
;train_canon_vars(*,0:nb-1) = transpose(trainlib)

;print, size(trainspec,/DIMENSIONS)
trainspec(*,nb) = specclass_ind  ; put in the class index
;train_canon_vars(*,nb) = specclass_ind
trainspec_wclass = trainspec   ; assign to Keely's variable she uses below
;trainspec_wclass = train_canon_vars

;Run LDA to get classification functions
print, 'starting lda'
IMSL_DISCR_ANALYSIS,$ ;(http://www.exelisvis.com/docs/IMSL_DISCR_ANALYSIS.html)
  trainspec_wclass,$ ;Two-dimensional array of size n_rows by n_variables + 1 containing the data in n_rows and n_variables = number of variables to be used in the discrimination. 
  (n_groups),$ ;Number of groups in the data.(n_groups-1)
  coefficients = LDAcoeffs,$ ;Named variable into which a two-dimensional array of size n_groups by (n_variables + 1) containing the linear discriminant coefficients is stored.
  class_member = tclass,$ ;Named variable into which an one-dimensional integer array of length n_rows containing the group to which the observation was classified is stored.
  class_table = terrmat, $ ;Named variable into which a two-dimensional array of size n_groups by n_groups containing the classification table is stored. 
  method = 3 ;Method of discrimination. Linear, pooled covariance computed, reclassification classification method
print,'discriminant analysis run'
correct = where(tclass eq specclass_ind)
overall_acc = 100*n_elements(correct)/n_elements(specclass_ind)
OA = string(overall_acc)
OA = strtrim(OA,2)

;Write LDAcoeffs to csv
WHILE (((N = STRPOS(trainlib_file, '.sli'))) NE -1) DO STRPUT, trainlib_file, '_LDA.csv', N
outfile = trainlib_file
WRITE_CSV,outfile,LDAcoeffs
end