pro image_classify_w_lda_wo_control_file_Part2

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
pnameTrain = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS\'; %Set directory
;pnameTrain = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS & MASTER\'; %Set directory
;pnameTrain = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\Combined\'; %Set directory

;;;;;;;;;;;;;;;Declare file name of the training library;;;;;;;;;;;;;;;;;;;;;;;
trainlib_file = strarr(1)
trainlib_file =  '20130411_Spectral_Library_AVIRIS_sorted'; %Set filename
;trainlib_file =  '20130606_Spectral_Library_AVIRIS_sorted'; %Set filename
;trainlib_file =  '20131125_Spectral_Library_AVIRIS_sorted'; %Set filename
;trainlib_file =  '20140416_Spectral_Library_AVIRIS_sorted'; %Set filename
;trainlib_file =  '20140606_Spectral_Library_AVIRIS_sorted'; %Set filename
;trainlib_file =  '20140829_Spectral_Library_AVIRIS_sorted'; %Set filename
;trainlib_file =  '20130411_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;trainlib_file =  '20130606_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;trainlib_file =  '20131125_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;trainlib_file =  '20140416_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;trainlib_file =  '20140606_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;trainlib_file =  '20140829_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;trainlib_file =  '2013&2014_Spectral_Library_AVIRIS'; %Set filename
;trainlib_file =  '2013&2014_Spectral_Library_AVIRIS&MASTER'; %Set filename
trainlib_file = trainlib_file + '_train_Spectral.sli' ;Variable for CDA coefficients name

;;;;;;;;;;;;;;;Declare Directory for Image;;;;;;;;;;;;;;;;;;;;;;;
pnameImage = strarr(1)
pnameImage = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Images\Mosaic\'; %Set directory

;;;;;;;;;;;;;;;Declare file name of image to be classified;;;;;;;;;;;;;;;;;;;;;;;
image_file=strarr(1)
image_file = '20130411_SBFrontRange_Mosaic'; %Set filename
;image_file = '20130606_SBFrontRange_Mosaic'; %Set filename
;image_file = '20131125_SBFrontRange_Mosaic'; %Set filename
;image_file = '20140416_SBFrontRange_Mosaic'; %Set filename
;image_file = '20140606_SBFrontRange_Mosaic'; %Set filename
;image_file = '20140829_SBFrontRange_Mosaic'; %Set filename;
;image_file = '20130411_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20130606_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20131125_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20140416_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20140606_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20140829_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
image_file = image_file + '_CDAvars'


trainlib_file = pnameTrain + trainlib_file
image_file = pnameImage + image_file
;CDAcoeff_file=pname+CDAcoeff_file

;open CDAcoeff file and read in array (nvars cols by ngoodbands rows)
LDAcoeffs = fltarr(n_groups,n_groups) ;change back nb to n_groups-1 for AVIRIS
get_lun, u1
openr, u1, LDAcoeff_file
readf, u1, LDAcoeffs
close, u1
free_lun, u1

;;;;;;;;;;;;;;;;;;;  IMAGE CLASSIFICATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;open image & get basic info
envi_open_file,image_file,r_fid=img_fid
ENVI_FILE_QUERY, img_fid,BBL = img_bbl,DATA_TYPE = data_type, DIMS=img_dims,INTERLEAVE=interleave,NB=imgnb,NL=imgnl,NS=imgns,WL=imgwl 
img_map_info=envi_get_map_info(fid=img_fid)

;create an output image array of all zeros
class_image = make_array(imgns,imgnl,/INTEGER,value=0)

;create an output LDA scores array
lda_score_image = make_array(imgns,imgnl,n_groups,/FLOAT,value=0)

;process image based on interleave (0=BSQ,1=BIL,2=BIP)
If interleave eq 0 then begin
    ;use envi_get_data with multiple calls
    curr_line = fltarr(imgns,imgnb)
    for l=0,imgnl-1 do begin      ;for each line in the image
      for b=0,imgnb-1 do begin    ;for each band in the line
        curr_line[*,b]= envi_get_data(dims=[-1L,0,imgns-1,l,l],fid=img_fid,pos=b)
      endfor
      ;apply bbl
      curr_line_data=curr_line[*,where(bbl eq 1)]   ;n samples x m bands
      ngoodbands = n_elements(curr_line_data[0,*])
;      ;multiply CDA coeffs through to get CDA variables
;      trans_CDAcoeffs=transpose(CDAcoeffs)   ;now m bands x p functions
;      curr_line_CDAvars=curr_line_data # trans_CDAcoeffs   ;n samples x p variables
;      
;     ;ignore background pixels in the line
      goodpix=where(total(curr_line_data,2)ne 0)
      if total(goodpix) eq -1 then class_image[*,l]=0 else begin
      ngoodpix=n_elements(goodpix)
      ;multiply LDA coeffs through get a classification score for each class
      for p=0,ngoodpix-1 do begin        ;for each pixel
        LDAscore=fltarr(n_groups)     ;vector of LDA scores for each class
        LDAvals=fltarr(1,ngoodbands+1)    ;LDA values for each variable + a constant
;        curr_pix_CDAvars=fltarr(n_groups-1)  
;        curr_pix_CDAvars=curr_line_CDAvars[p,*] ;current pixel's CDA variables
         curr_pix_spectrum=curr_line_data[goodpix[p],*]
        for g=0,n_groups-1 do begin
          LDAvals[1:ngoodbands]=curr_pix_spectrum*LDAcoeffs[g,1:ngoodbands]
          LDAvals[0]=LDAcoeffs[g,0]
          LDAscore[g]=total(LDAvals)
        endfor
        ;determine class membership based on highest LDA score
        maxscore=max(LDAscore)
        pix_class_ind=where(LDAscore eq maxscore)
		lda_score_image[goodpix[p],l,*] = LDAscore ;record the LDA score for each class
        class_image[goodpix[p],l]=pix_class_ind+1   ;add one because 0=unclassified  
      endfor
      endelse
    endfor
endif else begin
    ;use envi_get_slice for line by line
    curr_line=fltarr(imgns,imgnb)
    goodbands=where(bbl eq 1)  ;apply bbl before reading data in
    ngoodbands=n_elements(goodbands)
    for l=0,imgnl-1 do begin      ;for each line in the image
      curr_line=envi_get_slice(/BIL,line=l,fid=img_fid,pos=goodbands,xs=0,xe=imgns-1)  ;n samples x m bands
;     ;ignore background pixels in the line
      goodpix=where(total(curr_line_data,2)ne 0)
      if total(goodpix) eq -1 then class_image[*,l]=0 else begin
      ngoodpix=n_elements(goodpix)
      ;      ;multiply CDA coeffs through to get CDA variables
;      trans_CDAcoeffs=transpose(CDAcoeffs)   ;now m bands x p functions
;      curr_line_CDAvars=curr_line # trans_CDAcoeffs   ;n samples x p variables
      ;multiply LDA coeffs through get a classification score for each class
      for p=0,ngoodpix-1 do begin        ;for each pixel
        LDAscore=fltarr(n_groups)     ;vector of LDA scores for each class
        LDAvals=fltarr(n_groups)    ;LDA values for each variable + a constant
        curr_pix_spectrum=fltarr(ngoodbands) 
        curr_pix_spectrum=curr_line[goodpix[p],*] ;current pixel's CDA variables
        for g=0,n_groups-1 do begin
          LDAvals[1:ngoodbands]=curr_pix_spectrum*LDAcoeffs[g,1:ngoodbands]
          LDAvals[0]=LDAcoeffs[g,0]
          LDAscore[g]=total(LDAvals)
        endfor
        ;determine class membership based on highest LDA score
        maxscore=max(LDAscore)
        pix_class_ind=where(LDAscore eq maxscore)
        class_image[p,l]=pix_class_ind+1   ;add one because 0=unclassified  
		lda_score_image[p,l,*] = LDAscore ;record the LDA score for each class
      endfor
      endelse
    endfor
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; WRITE OUTPUT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;create a new classification image & notes
outfile = image_file+'_class'
notes='CDA/LDA classification results using '+trainlib_file+' as training data. Training accuracy was '+OA+'%.'

;write output image array to new class image
classlist_out=strarr(n_groups+1)
classlist_out[1:n_groups]=classlist
classlist_out[0]='unclassified'
file_type=envi_file_type('ENVI Classification')
lookup=intarr(3,n_groups+1)
lookup[*,0]=[0,0,0]
  for i=1,n_groups do begin 
    envi_get_rgb_triplets, i+2, r, g, b 
    lookup[*,i] = [r,g,b] 
  endfor

envi_write_envi_file, class_image,class_names=classlist_out,out_name=outfile,data_type=2, $
      file_type=file_type,descrip=notes, nb=1,nl=imgnl,ns=imgns,num_classes=n_groups+1,lookup=lookup
      
;create header with all info including training accuracy
outfile_hdr=outfile+'.hdr'
envi_setup_head,class_names=classlist_out,data_type=2,bnames=['classification results'],file_type=file_type,fname=outfile_hdr,$$
   interleave=0,map_info=img_map_info,nb=1,nl=imgnl,ns=imgns,num_classes=n_groups+1,offset=0,descrip=notes,lookup=lookup,/WRITE

;create a new scores image
outscores=image_file+'_scores'

; write outscores image
file_type1 = envi_file_type('ENVI Standard')
envi_write_envi_file, lda_score_image, bnames = classlist, nb=n_groups, nl=imgnl, ns=imgns, out_name=outscores
outscores_hdr=outscores+'.hdr'
envi_setup_head,data_type=4,bnames=classlist,file_type=file_type1,fname=outscores_hdr,$$
   interleave=0,map_info=img_map_info,nb=n_groups,nl=imgnl,ns=imgns,offset=0,/WRITE


   
print,'Image classified'


end