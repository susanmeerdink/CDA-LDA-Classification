pro create_cda_image_wo_control_file
;Updated Susan Meerdink 3/15/2016

; This program takes a set of CDA coefficients and a reflectance image and creates an output image of the CDA variables.

;This program no longer takes control file, but hard codes variables (see below)
;Control File Contains
; Line 1 = file directory pathname
; Line 2 = file name of the image
; Line 3 = CDA canon coeffs file name

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;Begin Main Program:
close, /ALL
;cfile=envi_pickfile(title='Select a Control file', filter='*.ctl')

;;;;;;;;;;;;;;;Declare Directory for Images;;;;;;;;;;;;;;;;;;;;;;;
pname = 'I:\AVIRIS\FL03\6 - Spectral Correction Files\'; %Set directory

;;;;;;;;;;;;;;;Declare Directory for CDA Coefficients;;;;;;;;;;;;;;;;;;;;;;;
pnameCDAFile = 'I:\Classification-Products\FL03\2 - CDA Variables\'; %Set directory
;pnameCDAFile = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS & MASTER\'; %Set directory

;;;;;;;;;;;;;;;Declare FileName;;;;;;;;;;;;;;;;;;;;;;;
image_file = strarr(1) ;Variable for image name
image_file = 'FL03_f140829t01p00r12_refl_hpc18_BIL_bottomhalf'; %Set filename

;Other inputs
CDAcoeff_file =  'f140829_AVIRIS_spectral_library'; %Set filename
CDAcoeff_file = CDAcoeff_file + '_CDAvars_CDAcoeffs.csv' ;Variable for CDA coefficients name

;read in control file information
CDAcoeff_file = pnameCDAFile + CDAcoeff_file ;Set cda coefficient file name
image_file = pname + image_file ;Set image filename

;open CDAcoeff file and read in array (nvars cols by ngoodbands rows)
CDAcoeffs_struct = read_csv(CDAcoeff_file)
nvars = n_tags(CDAcoeffs_struct)
ngoodbands = n_elements(CDAcoeffs_struct.field01)
CDAcoeffs=fltarr(nvars,ngoodbands)
for t=0,nvars-1 do CDAcoeffs[t,*]=CDAcoeffs_struct.(t)
;
;;;;;;;;;;;;;;;;;;;;  IMAGE CLASSIFICATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;open image & get basic info
envi_open_file,image_file,r_fid=img_fid
;ENVI_FILE_QUERY, img_fid,BBL=img_bbl,DATA_TYPE=data_type, DIMS=img_dims,INTERLEAVE=interleave,NB=imgnb,NL=imgnl,NS=imgns,WL=imgwl, file_type=file_type ;Original
ENVI_FILE_QUERY, img_fid,BBL = img_bbl, DATA_TYPE=data_type, DIMS=img_dims,INTERLEAVE=interleave,NB=imgnb,NL=imgnl,NS=imgns,WL=imgwl, file_type=file_type
map_info=envi_get_map_info(fid=img_fid)
outfile = image_file+'_CDAvars'

;img_bbl= make_array(191,1,VALUE = 1) ; For Combined AVIRIS + MASTER
;img_bbl= make_array(181,1,VALUE = 1) ;For AVIRIS
;create an output image array
;CDA_image=fltarr(imgns,imgnl,nvars)
;create an output line by line array
CDA_line=fltarr(imgns,nvars)
openw,1,outfile

;process image based on interleave (0=BSQ,1=BIL,2=BIP)
If interleave eq 0 then begin
  
    print,'Running CDA1'
    ;use envi_get_data with multiple calls
    curr_line=fltarr(imgns,imgnb)
    for l=0,imgnl-1 do begin      ;for each line in the image
      for b=0,imgnb-1 do begin    ;for each band in the line
        curr_line[*,b]=envi_get_data(dims=[-1L,0,imgns-1,l,l],fid=img_fid,pos=b)
      endfor
      ;apply bbl
      curr_line_data=curr_line[*,where(img_bbl eq 1)]   ;n samples x m bands
      ;ignore background pixels in the line
      curr_line_CDAvars=make_array(imgns,n_elements(CDAcoeffs[*,0]),/FLOAT,value=0)
      ;print,size(curr_line_data,/N_ELEMENTS)
      ;goodpix=where(total(curr_line,2)ne 0)
      ;goodpix=where(total(curr_line_data,1)ne 0)
      goodpix=where(total(curr_line_data,2)ne 0) ; ORIGINAL
      ;if total(goodpix) eq -1 then CDA_image[*,l,*]=curr_line_CDAvars else begin
      ;print,goodpix
      if total(goodpix) eq -1 then CDA_line=curr_line_CDAvars else begin
      ;multiply CDA coeffs through to get CDA variables
        trans_CDAcoeffs=transpose(CDAcoeffs)   ;now m bands x p functions
        ;print, size(trans_CDAcoeffs)
        curr_line_CDAvars[goodpix,*] = curr_line_data[goodpix,*] # trans_CDAcoeffs   ;n samples x p variables
        ;CDA_image[*,l,*]=curr_line_CDAvars     ;n samples x p variables
        ;CDA_line=curr_line_CDAvars
        ;writeu,1,CDA_line 
      endelse
      CDA_line=curr_line_CDAvars
      writeu,1,CDA_line
    endfor
;    endfor
endif else begin
    ;use envi_get_slice for line by line
    print,'Running CDA2'
    curr_line=fltarr(imgns,imgnb)
    goodbands=where(img_bbl eq 1)  ;apply bbl before reading data in
 for l=0,imgnl-1 do begin      ;for each line in the image
      curr_line=envi_get_slice(/BIL,line=l,fid=img_fid,pos=goodbands,xs=0,xe=imgns-1)  ;n samples x m bands
      ;ignore background pixels in the line
      curr_line_CDAvars=make_array(imgns,n_elements(CDAcoeffs[*,0]),/FLOAT,value=0)
      goodpix=where(total(curr_line,2)ne 0)
      ;if total(goodpix) eq -1 then CDA_image[*,l,*]=curr_line_CDAvars else begin
      if total(goodpix) eq -1 then CDA_line=curr_line_CDAvars else begin
        ;multiply CDA coeffs through to get CDA variables
        trans_CDAcoeffs=transpose(CDAcoeffs)   ;now m bands x p functions
        curr_line_CDAvars[goodpix,*]=curr_line[goodpix,*] # trans_CDAcoeffs   ;n samples x p variables
        ;CDA_image[*,l,*]=curr_line_CDAvars   ;add one because 0=unclassified  
      endelse
      CDA_line=curr_line_CDAvars
      writeu,1,CDA_line
 endfor
;    endfor
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; WRITE OUTPUT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


close,1 

;;create header with all info including validing accuracy
outfile_hdr=outfile+'.hdr'
file_type = ENVI_FILE_TYPE('ENVI Standard')
envi_setup_head,data_type=4,file_type=file_type,fname=outfile_hdr,$
   interleave=1,map_info=map_info,nb=nvars,nl=imgnl,ns=imgns,offset=0,/WRITE

print,'Image Done'


end