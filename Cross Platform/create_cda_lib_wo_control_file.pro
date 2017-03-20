pro create_cda_lib_wo_control_file


;; This program classifies an image using a EMing spectral library and a set of CDA coefficients output from Matlab.
;; It runs using a control file with the inputs specified below. The program applies the CDA coefficients to the 
;; EMing library spectral data to calculate canonical variables. It then inputs these to IDL's discriminant anlaysis
;; function to determine LDA coefficients. The image is read in line by line, and the CDA & LDA coefficients are applied.
;; These values yield an LDA score for each class. Each pixel in the image is assigned to the class with the highest 
;; LDA score. The program outputs the classfication image.

;This file no longer uses control file, input variables below
;Control File Contains
; Line 1 = file directory pathname
; Line 2 = file name of the endmember library
; Line 3 = CDA canon coeffs file name
; Line 4 = class column name


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Begin Main Program:
close, /ALL
;cfile=envi_pickfile(title='Select a Control file', filter='*.ctl')

;declare inputs
;;;;;;;;;;;;;;;Declare Directory;;;;;;;;;;;;;;;;;;;;;;;
pname=strarr(1)
;pname = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS\'; %Set directory
;pname = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS & MASTER\'; %Set directory
pname = 'I:\Classification-Products\FL03\2 - CDA Variables\'; %Set directory

;;;;;;;;;;;;;;;Declare filename of endmember library;;;;;;;;;;;;;;;;;;;;;;;
EMlib_file = strarr(1)
EMlib_file =  'f140829_AVIRIS_spectral_library'; %Set filename

EMlib_file = EMlib_file + '_calibration.sli' ;Variable for CDA coefficients name

;;;;;;;;;;;;;;;Declare filename of CDA Coefficients File;;;;;;;;;;;;;;;;;;;;;;;
CDAcoeff_file = strarr(1)
CDAcoeff_file =  'f140829_AVIRIS_spectral_library'; %Set filename
CDAcoeff_file = CDAcoeff_file + '_CDAvars_CDAcoeffs.csv' ;Variable for CDA coefficients name

;;;;;;;;;;;;;;;Declare class column name;;;;;;;;;;;;;;;;;;;;;;;
key_col_name = strarr(1)
key_col_name = 'Dominant' ;Same for all libraries (doesn't change)

;read in control file information
;openr,lun,cfile,/GET_LUN
;readf,lun,pname,EMlib_file,CDAcoeff_file,key_col_name
EMlib_file = pname + EMlib_file
;root=strsplit(EMlib_file,'.',/extract) ;Removed this line of code because there is a period in all users names on Ophelia
EMlib_filebase = STRMID(EMlib_file,0, STRLEN(EMlib_file)-4)
;EMlib_filebase=root[0];Removed this line of code because there is a period in all users names on Ophelia
coeff_name_parts = strsplit(CDAcoeff_file,'_',/extract)
coeff_name = coeff_name_parts[1]
CDAcoeff_file = pname+CDAcoeff_file

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; CREATE CLASSIFICATION VARIABLES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;read in EMing spectral library
envi_open_file,EMlib_file,r_fid=EM_fid
nspec = intarr(1)
nb=intarr(1)
envi_file_query,EM_fid,spec_names=specnames,nl=nspec,ns=nb,bbl=bbl,wl=wl,dims=tldims,interleave=interleave,data_type=data_type,file_type=file_type
;read in EMing library
EMlib = fltarr(nb,nspec)

EMlib=envi_get_data(fid=EM_fid, pos=0, dims=tldims)
;subset spectra to just good bands
EMlib_goodbands=EMlib[where(bbl eq 1,count),*]

;open EMing library metadata & read in lines
;root=strsplit(EMlib_file,".",/extract)
metadata_file=EMlib_filebase+'.csv'
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
endfor

;Create a numeric index for the class & determine number of classes
classlist=specclass[uniq(specclass)]          ;string vector of unique classes in the EMing data
n_groups = n_elements(classlist)
;n_groups = 18 ;Just for AVIRIS with 2013 & 2014
classlist_ind=indgen(n_groups)     ;numeric index for each class
specclass_ind=intarr(nspec)
for x=0,n_groups-1 do begin
  classmatch=where(strmatch(specclass,classlist[x],/FOLD_CASE))
  specclass_ind[classmatch]=classlist_ind[x]+1  ;vector with numeric class index for each spectrum in EMing lib
endfor

;open CDAcoeff file and read in array (nvars cols by ngoodbands rows)
CDAcoeffs = fltarr(n_groups-1,n_elements(where(bbl eq 1)))
;CDAcoeffs=fltarr(n_groups-1,n_elements(where(bbl eq 1))) ;Original
get_lun, u1
openr, u1, CDAcoeff_file
readf, u1, CDAcoeffs
close, u1
free_lun, u1


;Multiply the CDA coeffs through the EMing library to yield the EMing canon. variables
;EMlib is nbands (m) x nspec (n)
;CDAcoeffs is ng-1(p) x nbands(m)
;EM_canon_vars (pxn) = CDAcoeffs (pxm)* EMlib (mxn)
;CDA_EMs=#canonvars x #spec
;
;CDAcoeffs(p x m)*EMlib(m x n)=EM_canon_vars(p x n)
;EM_canon_vars is ng-1(p) x nspec(n)
EM_canon_vars=CDAcoeffs # EMlib_goodbands
nvars=n_elements(EM_canon_vars[*,1])

;;Format the variables for LDA & add class membership #
;EM_canon_vars=transpose(EM_canon_vars) ;should have a column for each spectrum and ng-1 rows
;EM_canon_vars_wclass=fltarr(nspec,n_groups) ;create a new column for the classID #
;EM_canon_vars_wclass[*,n_groups-1]=specclass_ind   ;add classID#
;EM_canon_vars_wclass[*,0:n_groups-2]=EM_canon_vars  ;put in EM_canon_vars
;
;;Determine if the number of observations per class limits the number of discriminant variables to use
;
;;Run LDA to get classification functions
;IMSL_DISCR_ANALYSIS,EM_canon_vars_wclass,n_groups,coefficients=LDAcoeffs,class_member=tclass,class_table=terrmat, method=3
;print,'discriminant analysis run'
;correct=where(tclass eq specclass_ind)
;overall_acc=100*n_elements(correct)/n_elements(specclass_ind)
;OA=string(overall_acc)
;OA=strtrim(OA,2)
;
;;;;;;;;;;;;;;;;;;;;  IMAGE CLASSIFICATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;;open image & get basic info
;envi_open_file,image_file,r_fid=img_fid
;ENVI_FILE_QUERY, img_fid,BBL=img_bbl,DATA_TYPE=data_type, DIMS=img_dims,INTERLEAVE=interleave,NB=imgnb,NL=imgnl,NS=imgns,WL=imgwl 
;map_info=envi_get_map_info(fid=img_fid)
;
;;create an output image array
;class_image=intarr(imgns,imgnl)
;
;;process image based on interleave (0=BSQ,1=BIL,2=BIP)
;If interleave eq 0 then begin
;    ;use envi_get_data with multiple calls
;    curr_line=fltarr(imgns,imgnb)
;    for l=0,imgnl-1 do begin      ;for each line in the image
;      for b=0,imgnb-1 do begin    ;for each band in the line
;        curr_line[*,b]=envi_get_data(dims=[-1L,0,imgns-1,l,l],fid=img_fid,pos=b)
;      endfor
;      ;apply bbl
;      curr_line_data=curr_line[*,where(bbl eq 1)]   ;n samples x m bands
;      ;multiply CDA coeffs through to get CDA variables
;      trans_CDAcoeffs=transpose(CDAcoeffs)   ;now m bands x p functions
;      curr_line_CDAvars=curr_line_data # trans_CDAcoeffs   ;n samples x p variables
;      ;multiply LDA coeffs through get a classification score for each class
;      for p=0,imgns-1 do begin        ;for each pixel
;        LDAscore=fltarr(n_groups)     ;vector of LDA scores for each class
;        LDAvals=fltarr(1,n_groups)    ;LDA values for each variable + a constant
;        curr_pix_CDAvars=fltarr(n_groups-1)  
;        curr_pix_CDAvars=curr_line_CDAvars[p,*] ;current pixel's CDA variables
;        for g=0,n_groups-1 do begin
;          LDAvals[1:n_groups-1]=curr_pix_CDAvars*LDAcoeffs[g,1:n_groups-1]
;          LDAvals[0]=LDAcoeffs[g,0]
;          LDAscore[g]=total(LDAvals)
;        endfor
;        ;determine class membership based on highest LDA score
;        maxscore=max(LDAscore)
;        pix_class_ind=where(LDAscore eq maxscore)
;        class_image[p,l]=pix_class_ind+1   ;add one because 0=unclassified  
;      endfor
;    endfor
;endif else begin
;    ;use envi_get_slice for line by line
;    curr_line=fltarr(imgns,imgnb)
;    goodbands=where(bbl eq 1)  ;apply bbl before reading data in
;    for l=0,imgnl-1 do begin      ;for each line in the image
;      curr_line=envi_get_slice(/BIL,line=l,fid=img_fid,pos=goodbands,xs=0,xe=imgns-1)  ;n samples x m bands
;      ;multiply CDA coeffs through to get CDA variables
;      trans_CDAcoeffs=transpose(CDAcoeffs)   ;now m bands x p functions
;      curr_line_CDAvars=curr_line # trans_CDAcoeffs   ;n samples x p variables
;      ;multiply LDA coeffs through get a classification score for each class
;      for p=0,imgns-1 do begin        ;for each pixel
;        LDAscore=fltarr(n_groups)     ;vector of LDA scores for each class
;        LDAvals=fltarr(n_groups)    ;LDA values for each variable + a constant
;        curr_pix_CDAvars=fltarr(n_groups-1)  
;        curr_pix_CDAvars=curr_line_CDAvars[p,*] ;current pixel's CDA variables
;        for g=0,n_groups-1 do begin
;          LDAvals[1:n_groups-1]=curr_pix_CDAvars*LDAcoeffs[g,1:n_groups-1]
;          LDAvals[0]=LDAcoeffs[g,0]
;          LDAscore[g]=total(LDAvals)
;        endfor
;        ;determine class membership based on highest LDA score
;        maxscore=max(LDAscore)
;        pix_class_ind=where(LDAscore eq maxscore)
;        class_image[p,l]=pix_class_ind+1   ;add one because 0=unclassified  
;      endfor
;    endfor
;endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; WRITE OUTPUT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;create a new classification image & notes
;outfile=EMlib_filebase+'_'+'.sli'
;notes='CDA/LDA classification results using '+EMlib_file+' as EMing data. EMing accuracy was '+OA+'%.'
;
;;write output image array to new class image
;classlist_out=strarr(n_groups+1)
;classlist_out[1:n_groups]=classlist
;classlist_out[0]='unclassified'
;file_type=envi_file_type('ENVI Classification')
;lookup=intarr(3,n_groups+1)
;lookup[*,0]=[0,0,0]
;  for i=1,n_groups do begin 
;    envi_get_rgb_triplets, i+2, r, g, b 
;    lookup[*,i] = [r,g,b] 
;  endfor

 outlib = EMlib_filebase + '_CDA.sli'
  openw,1,outlib
  writeu,1,EM_canon_vars
  close,1
  new_head_name = EMlib_filebase+ '_CDA.hdr'
  file_type=envi_file_type('ENVI Spectral Library')
  envi_setup_head, fname=new_head_name,ns=nvars, nl=nspec, nb=1, spec_names=specnames, data_type=data_type,$
   interleave=interleave, file_type=file_type, /write
      
;;create header with all info including EMing accuracy
;outfile_hdr=outfile+'.hdr'
;envi_setup_head,class_names=classlist_out,data_type=2,bnames=['classification results'],file_type=file_type,fname=outfile_hdr,$$
;   interleave=0,map_info=map_info,nb=1,nl=imgnl,ns=imgns,num_classes=n_groups+1,offset=0,descrip=notes,lookup=lookup,/WRITE

print,'Library Done'


end