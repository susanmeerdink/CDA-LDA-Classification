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
pname = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Images\Mosaic\'; %Set directory

;;;;;;;;;;;;;;;Declare Directory for CDA Coefficients;;;;;;;;;;;;;;;;;;;;;;;
pnameCDAFile = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS\'; %Set directory
;pnameCDAFile = 'H:\users\meerdink\Dropbox\AAG_2016_Research\Spectral Libraries\AVIRIS & MASTER\'; %Set directory

;;;;;;;;;;;;;;;Declare FileName;;;;;;;;;;;;;;;;;;;;;;;
image_file = strarr(1) ;Variable for image name
;image_file = '20130411_SBFrontRange_Mosaic'; %Set filename
;image_file = '20130606_SBFrontRange_Mosaic'; %Set filename
;image_file = '20131125_SBFrontRange_Mosaic'; %Set filename
;image_file = '20140416_SBFrontRange_Mosaic'; %Set filename
;image_file = '20140606_SBFrontRange_Mosaic'; %Set filename
image_file = '20140829_SBFrontRange_Mosaic'; %Set filename
;image_file = '20130411_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20130606_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20131125_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20140416_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20140606_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename
;image_file = '20140829_SBFrontRange_Mosaic_AVIRIS&MASTER'; %Set filename


;Other inputs
;CDAcoeff_file =  '20130411_Spectral_Library_AVIRIS_sorted'; %Set filename
;CDAcoeff_file =  '20130606_Spectral_Library_AVIRIS_sorted'; %Set filename
;CDAcoeff_file =  '20131125_Spectral_Library_AVIRIS_sorted'; %Set filename
;CDAcoeff_file =  '20140416_Spectral_Library_AVIRIS_sorted'; %Set filename
;CDAcoeff_file =  '20140606_Spectral_Library_AVIRIS_sorted'; %Set filename
CDAcoeff_file =  '20140829_Spectral_Library_AVIRIS_sorted'; %Set filename
;CDAcoeff_file =  '20130411_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;CDAcoeff_file =  '20130606_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;CDAcoeff_file =  '20131125_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;CDAcoeff_file =  '20140416_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;CDAcoeff_file =  '20140606_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;CDAcoeff_file =  '20140829_Spectral_Library_AVIRIS&MASTER_sorted'; %Set filename
;CDAcoeff_file =  '2013&2014_Spectral_Library_AVIRIS'; %Set filename
;CDAcoeff_file =  '2013&2014_Spectral_Library_AVIRIS&MASTER'; %Set filename
CDAcoeff_file = CDAcoeff_file + '_CDAvars_CDAcoeffs.csv' ;Variable for CDA coefficients name

;read in control file information
;openr,lun,cfile,/GET_LUN
;readf,lun,pname,image_file,CDAcoeff_file
CDAcoeff_file = pnameCDAFile + CDAcoeff_file ;Set cda coefficient file name
image_file = pname + image_file ;Set image filename


;open CDAcoeff file and read in array (nvars cols by ngoodbands rows)
CDAcoeffs_struct = read_csv(CDAcoeff_file)
nvars = n_tags(CDAcoeffs_struct)
ngoodbands = n_elements(CDAcoeffs_struct.field01)
CDAcoeffs=fltarr(nvars,ngoodbands)
for t=0,nvars-1 do CDAcoeffs[t,*]=CDAcoeffs_struct.(t)


;CDAcoeffs=fltarr(n_groups-1,n_elements(where(bbl eq 1)))
;get_lun, u1
;openr, u1, CDAcoeff_file
;readf, u1, CDAcoeffs
;close, u1
;free_lun, u1


;Multiply the CDA coeffs through the validing library to yield the validing canon. variables
;validlib is nbands (m) x nspec (n)
;CDAcoeffs is ng-1(p) x nbands(m)
;valid_canon_vars (pxn) = CDAcoeffs (pxm)* validlib (mxn)
;CDA_valids=#canonvars x #spec
;
;CDAcoeffs(p x m)*validlib(m x n)=valid_canon_vars(p x n)
;;valid_canon_vars is ng-1(p) x nspec(n)
;valid_canon_vars=CDAcoeffs # validlib_goodbands
;nvars=n_elements(valid_canon_vars[*,1])

;;Format the variables for LDA & add class mvalidbership #
;valid_canon_vars=transpose(valid_canon_vars) ;should have a column for each spectrum and ng-1 rows
;valid_canon_vars_wclass=fltarr(nspec,n_groups) ;create a new column for the classID #
;valid_canon_vars_wclass[*,n_groups-1]=specclass_ind   ;add classID#
;valid_canon_vars_wclass[*,0:n_groups-2]=valid_canon_vars  ;put in valid_canon_vars
;
;;Determine if the number of observations per class limits the number of discriminant variables to use
;
;;Run LDA to get classification functions
;IMSL_DISCR_ANALYSIS,valid_canon_vars_wclass,n_groups,coefficients=LDAcoeffs,class_mvalidber=tclass,class_table=terrmat, method=3
;print,'discriminant analysis run'
;correct=where(tclass eq specclass_ind)
;overall_acc=100*n_elvalidents(correct)/n_elvalidents(specclass_ind)
;OA=string(overall_acc)
;OA=strtrim(OA,2)
;
;;;;;;;;;;;;;;;;;;;;  IMAGE CLASSIFICATION ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;open image & get basic info
envi_open_file,image_file,r_fid=img_fid
;ENVI_FILE_QUERY, img_fid,BBL=img_bbl,DATA_TYPE=data_type, DIMS=img_dims,INTERLEAVE=interleave,NB=imgnb,NL=imgnl,NS=imgns,WL=imgwl, file_type=file_type ;Original
ENVI_FILE_QUERY, img_fid,DATA_TYPE=data_type, DIMS=img_dims,INTERLEAVE=interleave,NB=imgnb,NL=imgnl,NS=imgns,WL=imgwl, file_type=file_type
map_info=envi_get_map_info(fid=img_fid)
outfile = image_file+'_CDAvars'

;img_bbl= make_array(191,1,VALUE = 1) ; For Combined AVIRIS + MASTER
img_bbl= make_array(186,1,VALUE = 1) ;For AVIRIS
;create an output image array
;CDA_image=fltarr(imgns,imgnl,nvars)
;create an output line by line array
CDA_line=fltarr(imgns,nvars)
openw,1,outfile

;process image based on interleave (0=BSQ,1=BIL,2=BIP)
If interleave eq 0 then begin
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
        CDA_line=curr_line_CDAvars
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
;        ;determine class mvalidbership based on highest LDA score
;        maxscore=max(LDAscore)
;        pix_class_ind=where(LDAscore eq maxscore)
        ;add one because 0=unclassified 
        writeu,1,CDA_line 
      endelse
    endfor
;    endfor
endif else begin
    ;use envi_get_slice for line by line
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
        CDA_line=curr_line_CDAvars
        writeu,1,CDA_line 
      endelse
 endfor
;    endfor
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; WRITE OUTPUT ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;create a new classification image & notes
;outfile=validlib_filebase+'_CDA.sli'
;notes='CDA/LDA classification results using '+validlib_file+' as validing data. validing accuracy was '+OA+'%.'
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

; outlib=validlib_filebase+'_CDA.sli'
;  openw,1,outlib
;  writeu,1,valid_canon_vars
;  close,1
;  new_head_name = validlib_filebase+'_CDA.hdr'
;  file_type=envi_file_type('ENVI Spectral Library')
;  envi_setup_head, fname=new_head_name,ns=nvars, nl=nspec, nb=1, spec_names=specnames, data_type=data_type,$
;   bbl=bbl,interleave=interleave, wl=wl, offset=offset, file_type=file_type, /write
      

;envi_write_envi_file, CDA_image,out_name=outfile,data_type=4, $
;      file_type=file_type, nb=nvars,nl=imgnl,ns=imgns   

close,1 

;;create header with all info including validing accuracy
outfile_hdr=outfile+'.hdr'
file_type = ENVI_FILE_TYPE('ENVI Standard')
envi_setup_head,data_type=4,file_type=file_type,fname=outfile_hdr,$
   interleave=1,map_info=map_info,nb=nvars,nl=imgnl,ns=imgns,offset=0,/WRITE

print,'Image Done'


end