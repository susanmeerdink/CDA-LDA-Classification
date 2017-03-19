pro sample_spectral_lib_v7

  COMPILE_OPT IDL2
  ;MUST SET DATA SORT TO BE BY DOMINANT(1) then POLYGON(2)

  ;Use this code with a control file to create training and validation libraries from an input spectral library
  ;using user-defined thresholds for sampling from each dominant species randomly. Two libraries with associated
  ;metadata are output, one training data set and one validation data set.
  ;
  ;Created 1/27/10 by Keely Roth and last updated 12/1/2015
  ;This code was altered so that a certain percent of polygons were set aside at beginning
  ;
  ;Updated Kenneth Dudley changed write_csv, to printf to prevent "" addition breaking compatibility with VIPER tools.
  ; Cleaned up some code to make program less verbose, changed print outputs (to console) to make more sense to th user,
  ; ADDED internal sorting, the input sort order of the library no longer matters, and will not be modified by the program
  ;   train/valid library will maintain the same sort order as the input library.
  ;**************************************************************************************
  ;control file holds the following:
  ;Line 1 is the directory pathname
  ;Line 2 is the input spectral library name (no file extension attached)
  ;Line 3 is the class (column heading) to stratify sampling by
  ;Line 4 is the proportion (0-1) of polygons to be held aside for independent validation (20%)
  ;Line 5 is the absolute (integer) sampling limit for a polygon
  ;Line 6 is the proportional (0-1) sampling limit for a polygon

  ;Begin Program:
  close, /ALL
  cfile=envi_pickfile(title='Select a Control file', filter='*.txt')

  ;;;;;;;;read in control file information;;;;;;;;;;;;;;;;;;;;;;
  openr,1,cfile
  pname=strarr(1)
  readf,1,pname
  in_lib=strarr(1)
  readf,1,in_lib
  in_lib = pname + in_lib
  in_lib_metadata = in_lib + '.csv'
  in_lib_spectra = in_lib + '.sli'
  train_lib=strarr(1);Create string array to hold training library name
  key_class=strarr(1);Create string array to class to stratify sampling by
  readf,1,key_class ;read in class
  ind_limit = fltarr(1);Create float array to hold proportional (0-1) of polygons to be held aside for independent validation
  readf,1,ind_limit ;read in proportion for indenpendent validation
  num_limit=intarr(1);Create integer array to hold absolute (integer) sampling limit for a polygon
  readf,1,num_limit ;read in absolute limit
  prop_limit=fltarr(1);Create float array to hold proportional (0-1) sampling limit for a polygon
  readf,1,prop_limit ;read in proportional limit
  out_summ_file=strarr(1)

  ;;;;;;;;;create training library variables/arrays;;;;;;;;;;;;;;;;;;
  train_lib=in_lib+'_train'
  train_lib_metadata=train_lib+'.csv'
  train_lib_spectra=train_lib+'.sli'
  valid_lib=strarr(1)
  ;;;;;;;;;;create validation libraray variables/arrays;;;;;;;;;;;;;;;
  valid_lib=in_lib+'_valid'
  valid_lib_metadata=valid_lib+'.csv'
  valid_lib_spectra=valid_lib+'.sli'
  ;;;;;;;;;;create independent validation libraray variables/arrays;;;;;;;;;;;;;;;
  independent_lib = in_lib+'_independent'
  independent_lib_metadata = independent_lib+'.csv'
  independent_lib_spectra = independent_lib+'.sli'
  ;;;;;;;;;create output summary report file name;;;;;;;;;;;;;;;;;;;;;;;;;
  out_summ_file=in_lib+'_sample_summary.txt'
  
  ;;;;;;;;;;; open the spectral library and read in data;;;;;;;;;;;;;;;;;;;;;;;;;;
  envi_open_file, in_lib_spectra, r_fid=sli_fid
  envi_file_query, sli_fid, dims=dims, nb=nb, ns=ns, nl=nl, spec_names=spec_names, $
    interleave=interleave,byte_order=byte_order, bbl=bbl, wl=wl, data_type=data_type, $
    offset=offset, file_type=file_type
  speclib = envi_get_data(fid=sli_fid, dims=dims, pos=[0])

  ;;;;;;;;;;; open the metadata file and read in columns ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Read in header separately, since read_csv will not find the header if next line does not contain numerics
  header_csv = read_csv(in_lib_metadata, NUM_RECORDS=0, RECORD_START=0) 
  col_headings = strarr(n_tags(header_csv))
  col_headings = strtrim(col_headings,2)
  for i=0, n_tags(header_csv)-1 do col_headings[i] = header_csv.(i)

  ;Read in rest of metadata
  in_sli_metadata=read_csv(in_lib_metadata,count=nspec, RECORD_START=1)
  n_cols=n_elements(col_headings)
  ; add an array for train/valid denotation
  voidvals=strarr(nspec)
  voidvals[*]='o'
  void=create_struct('t_v',voidvals)
  in_sli_metadata_updated=create_struct(in_sli_metadata,void);The variable will hold the new updated spectral library
  
  ;************************************************************************************

  ;;;;;Find dominants, and associated polygons (in the proper order) for any library, sorted or unsorted.;;;;;;;;;;;;;;;;;;;
  dominant_index=where(strmatch(col_headings,key_class,/FOLD_CASE) eq 1)
  in_sli_metadata_updated.(dominant_index) = strtrim(strupcase(in_sli_metadata_updated.(dominant_index)),2)
  dominants = in_sli_metadata_updated.(dominant_index)[UNIQ(in_sli_metadata_updated.(dominant_index), SORT(in_sli_metadata_updated.(dominant_index)))]
  nf1uniq = n_elements(dominants)
  polygon_index=where(strmatch(col_headings,'Polygon_ID',/FOLD_CASE) eq 1)
  ; generate output order
  polygons_in_lib = !NULL
  for i=0,nf1uniq-1 do begin
    index = where(in_sli_metadata_updated.(dominant_index) eq dominants[i])
    field2sub = in_sli_metadata_updated.(polygon_index)[index]
    uniq_poly = field2sub[uniq(field2sub, sort(field2sub))]
    polygons_in_lib = [polygons_in_lib,uniq_poly]
  endfor
  
  ;;;;;;;;;;;;;;store the unique dominants present in library;;;;;;;;;;;;;;;;;;;;;
  n_doms = n_elements(dominants)
  dominant_summary={name:strarr(n_doms),n_polys:intarr(n_doms),total_pix:intarr(n_doms),train_pix:intarr(n_doms),perc_sample:fltarr(n_doms)}
  ;Dominant summary is a structure that holds multiple variables including:
  ;name = (string) holds the dominant species name
  ;n_polys = (integer) holds the number of polygons for a dominant species
  ;total_pix = (integer) hold the total number of pixels of the dominant species
  ;train_pix = (integer) holds the number of pixels used in training for the dominant species
  ;perc_sample = (float) holds the percentage of pixels used in training for thedominant species
  dominant_summary.name = dominants ;Set the dominant species name

  ;;;store the unique polygons present in the library;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  total_polys=n_elements(polygons_in_lib)
  polygon_summary={name:strarr(total_polys),dominant:strarr(total_polys),n_pix:intarr(total_polys),$
    limit:intarr(total_polys),flag:bytarr(total_polys),train_pix:intarr(total_polys),perc_sample:fltarr(total_polys),independent:strarr(total_polys)}
  ;Polygon summary is a structure that holds multiple variables including:
  ;name = (string) holds the polygons name
  ;dominant = (string) holds the dominant species name of this polygon
  ;n_pix = (integer) holds the number of pixels of the polygon
  ;limit = (integer) 
  ;flag = (byte)
  ;train_pix = (integer) holds the number of pixels in the polygon used for training
  ;perc_sample = (float) holds the percentage of pixels in the polygon used for training
  ;independent = (string) yes means the polygon was held aside for independent validation while no means it was split into training and validation
  polygon_summary.name = polygons_in_lib ;Set the name of the polygons
  
  ;;;summarize the metadata for dominants and polygons;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  poly_offset=0L ; Set counter for writing to differing polygon indices
  ;for each dominant
  print, FORMAT='(A-12,A-10)', 'Class', 'Polygons'
  for i=0L,n_doms-1 do begin
    ;assign dominant name
    curr_dom = dominant_summary.name[i]
    ;determine # polygons in the dominant
    polys = in_sli_metadata_updated.(polygon_index)[where(in_sli_metadata_updated.(dominant_index) eq curr_dom)]
    polys_in_dom = polys[UNIQ(polys, SORT(polys))]
    n_polys = n_elements(polys_in_dom)
    dominant_summary.n_polys[i] = n_polys
    ;begin a pixel counter for this dominant
    pix_count=0L ; Cumulative count of number of pixels for a dominant class, sum of between polygons
    ;for each polygon
    for j=0L,n_polys-1 do begin
      curr_poly = polys_in_dom[j]
      ;determine pixels present and create structure to hold
      pixel_names = in_sli_metadata_updated.(0)[where(in_sli_metadata_updated.(polygon_index) eq curr_poly)]
      n_pix_c = n_elements(pixel_names)
      pix_count = pix_count+n_pix_c ; This gives the overall pixel count for class
      pp_index = poly_offset+j ; index of storage for number of pixels for each polygon within each dominant class
      polygon_summary.n_pix[pp_index] = n_pix_c ; Store number of pixels in a polygon for a specific dominant class
      polygon_summary.dominant[pp_index] = curr_dom ; Store the associated dominant class with the above number of pixels
      ;read each pixel name and spectral library index into the structure
    endfor
    poly_offset = pp_index+1L
    print, FORMAT='(A-12,I-8)', curr_dom, n_polys
    dominant_summary.total_pix[i]=pix_count
  endfor
  
  ;;;Calculate polygon sampling limits;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  for i=0,total_polys-1 do begin
    if (polygon_summary.n_pix[i] ge (1/prop_limit)*num_limit) then begin
      polygon_summary.limit[i]= num_limit
    endif else begin
      polygon_summary.limit[i] = floor(polygon_summary.n_pix[i]*prop_limit)
      if polygon_summary.limit[i] eq 0 then polygon_summary.flag[i]=1 ; no sampling?
    endelse
  endfor
  
  
  for i=0L, n_doms-1 do begin ;Loop through the dominant species
    ;;;Determine which polyons will be used for independent validation and which ones will be split into training/validation;;;;;;;;;;;;;;;;
    ;generating M random numbers, check for duplicates, and generate M-n more, accumulating until you have enough
    M = round((ind_limit*dominant_summary.n_polys[i])) ;Number of Polygons to be pulled out for independent validation
    len = dominant_summary.n_polys[i] ;Total number of polygons for a dominant species
    n = M

    count = 1
    WHILE n GT 0 DO BEGIN
      inds = make_array(n,1)
      inds = round(RandomU(seed, n)*len)
      u = inds[Uniq(inds)];following expression returns a copy of the sorted array with duplicate adjacent elements removed
      n = M - N_Elements(u)
      inds = u
      count = count+1
      if count GE 10000 then break
    ENDWHILE

    random_ind_index = make_array(len,1,/BYTE,VALUE = 1);create a byte array of 1s that is the length of polygons
    for x = 0,(dominant_summary.n_polys[i]-1) do begin ;loop through number of polygons
      for z = 0,(size(inds,/N_ELEMENTS)-1) do begin ;loop through random indices
        if x EQ inds[z] then begin ;If the polygon index matches the random index set it to zero
          random_ind_index[x] = 0
        endif
      endfor
    endfor
    
    spectra_indices = where(in_sli_metadata_updated.(dominant_index) eq dominant_summary.name[i]) ;Get the indices that belong to the dominant specie
    train_count = 0L
    closed_flag_array = bytarr(dominant_summary.n_polys[i]) ;creates a byte array of the number of polygons for the dominant specie
    closed_flag_array += byte(1)
    random_indices = lonarr(dominant_summary.total_pix[i]) ;create a longword array of all the pixels of a dominant specie
    r_ind_count = 0L
    
    ;loop through polygons
    for p = 0, (dominant_summary.n_polys[i]-1) do begin ;loop through polygons in dominant specie
      
      if random_ind_index[p] EQ 0 then begin ;If the random index is set to zero pull polygon out for independent validation
        print, dominant_summary.name[i]
        polygon_summary.independent[p] = 1 ;set independent validation flag to yes
        print,polygon_summary.(7)[p]
        polygon_summary.train_pix[p] = 0 ;Set training pixels to 0
        for a = 0, (polygon_summary.n_pix[i]-1) do begin ;Loop through all pixels of polygon 
          in_sli_metadata_updated.(n_cols)[spectra_indices[a]]='i' ;update the pixel metadata to independent validation
        endfor
      endif else begin ;If the random index is set to 1 split polygon into training and validation
        polygon_summary.independent[p] = 0 ;set independent validation flag to no

        ;Set random index to pull out training data
        new_random: random_multiplier=randomu(S,1)
        random_index = round((polygon_summary.n_pix[i]-1)*random_multiplier)
        repeat_indices = where(random_indices eq random_index[0])
        test_num = total(repeat_indices)
        if test_num ge 0 then begin
          goto,new_random
        endif else begin
          random_indices[r_ind_count] = random_index
        endelse
        r_ind_count++

        sampled_poly_name = in_sli_metadata_updated.(polygon_index)[spectra_indices[random_index]] ;Find the sample polygon name
        sampled_poly = where(polygon_summary.name eq sampled_poly_name[0])
        if polygon_summary.flag[sampled_poly] ne 1 then begin
          in_sli_metadata_updated.(n_cols)[spectra_indices[random_index]]='t';update the pixel metadata to training
          polygon_summary.train_pix[sampled_poly]++ ;increase the count that is keeping track how many pixels are used in training
          train_count++
          if polygon_summary.train_pix[sampled_poly] eq polygon_summary.limit[sampled_poly] then polygon_summary.flag[sampled_poly]=1
        endif
      endelse  
      
    endfor ;End of loop for polygonsr
    dominant_summary.train_pix[i] = train_count
    dominant_summary.perc_sample[i] = float(dominant_summary.train_pix[i])/float(dominant_summary.total_pix[i])*100 
  endfor ;End of loop for dominant species
  
  for j=0,total_polys-1 do polygon_summary.perc_sample[j]=float(polygon_summary.train_pix[j])/float(polygon_summary.n_pix[j])*100
  
  ;;;set all non-training pixels to validation;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  non_sampled_spectra = where(in_sli_metadata_updated.(n_cols) ne 't' AND in_sli_metadata_updated.(n_cols) ne 'i')
  in_sli_metadata_updated.(n_cols)[non_sampled_spectra]='v'
  train_spectra_indices = where(in_sli_metadata_updated.(n_cols) eq 't',n_train_spectra)
  valid_spectra_indices = where(in_sli_metadata_updated.(n_cols) eq 'v',n_valid_spectra)
  independent_spectra_indices = where(in_sli_metadata_updated.(n_cols) eq 'i',n_independent_spectra)

  ;split metadata into two arrays for output
  library_metadata_train=strarr(n_tags(in_sli_metadata_updated),n_train_spectra)
  library_metadata_valid=strarr(n_tags(in_sli_metadata_updated),n_valid_spectra)
  library_metadata_independent=strarr(n_tags(in_sli_metadata_updated),n_independent_spectra)
  for i=0,n_cols do begin
    library_metadata_train[i,*]=in_sli_metadata_updated.(i)[train_spectra_indices]
    library_metadata_valid[i,*]=in_sli_metadata_updated.(i)[valid_spectra_indices]
    library_metadata_independent[i,*]=in_sli_metadata_updated.(i)[independent_spectra_indices]
  endfor
  new_col_headings=strarr(n_cols+1)
  new_col_headings[0:n_cols-1]=col_headings
  new_col_headings[n_cols]='train_valid_independent'
  
  ;;;;;;;;;;;;;;;;;;write output files;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;Training Library
  ;spectral libraries + header files + metadata files
  get_lun, a
  openw, a, train_lib_spectra
  writeu, a, speclib[*,train_spectra_indices]
  close, a
  free_lun, a

  train_head_name = train_lib + '.hdr'
  envi_setup_head, fname=train_head_name,ns=ns, nl=n_train_spectra, nb=nb, spec_names=library_metadata_train[0,*], data_type=data_type,$
    bbl=bbl,interleave=interleave, wl=wl, offset=offset, file_type=file_type, /write

  library_metadata_train = strtrim(library_metadata_train,2)
  openw, u1, train_lib_metadata, /get_lun
  printf, u1, strjoin(new_col_headings,',',/single)
  for i=0L,n_train_spectra-1 do printf,u1, strjoin(library_metadata_train[*,i],',',/single)
  close,u1
  free_lun,u1

  ;Validation Library
  get_lun, b
  openw, b, valid_lib_spectra
  writeu, b, speclib[*,valid_spectra_indices]
  close, b
  free_lun, b

  valid_head_name = valid_lib + '.hdr'
  envi_setup_head, fname=valid_head_name,ns=ns, nl=n_valid_spectra, nb=nb, spec_names=library_metadata_valid[0,*], data_type=data_type,$
    bbl=bbl,interleave=interleave, wl=wl, offset=offset, file_type=file_type, /write

  library_metadata_valid = strtrim(library_metadata_valid,2)
  openw, u1, valid_lib_metadata, /get_lun
  printf, u1, strjoin(new_col_headings,',',/single)
  for i=0L,n_valid_spectra-1 do printf,u1, strjoin(library_metadata_valid[*,i],',',/single)
  close,u1
  free_lun,u1

  ;Validation Library
  get_lun, b
  openw, b, independent_lib_spectra
  writeu, b, speclib[*,independent_spectra_indices]
  close, b
  free_lun, b

  independent_head_name = independent_lib + '.hdr'
  envi_setup_head, fname=independent_head_name,ns=ns, nl=n_independent_spectra, nb=nb, spec_names=library_metadata_valid[0,*], data_type=data_type,$
    bbl=bbl,interleave=interleave, wl=wl, offset=offset, file_type=file_type, /write

  library_metadata_independent = strtrim(library_metadata_independent,2)
  openw, u1, independent_lib_metadata, /get_lun
  printf, u1, strjoin(new_col_headings,',',/single)
  for i=0L,n_independent_spectra-1 do printf,u1, strjoin(library_metadata_independent[*,i],',',/single)
  close,u1
  free_lun,u1

  ;update full library
  ;write_csv, in_lib_metadata, in_sli_metadata_updated, header=new_col_headings
  ncol = n_tags(in_sli_metadata_updated)
  meta_array=strarr(ncol,nl)
  for i=0, ncol-1 do begin ; convert metadata file to a string array
    meta_array[i,*] = in_sli_metadata_updated.(i)
  endfor
  meta_array = strtrim(meta_array,2) ; Remove extra spaces carried over from metadata (if they exist)

  ;Alternative to write_csv
  openw, u1, in_lib_metadata, /get_lun
  printf, u1, strjoin(new_col_headings,',',/single)
  for i=0L,nl-1 do printf,u1, strjoin(meta_array[*,i],',',/single)
  close,u1
  free_lun,u1

  ;sampling summary report
  get_lun, c
  openw,c,out_summ_file
  printf,c,'Class Summary'
  printf,c,tag_names(dominant_summary)
  for i=0,n_doms-1 do begin
    printf,c,dominant_summary.(0)[i],'  ',dominant_summary.(1)[i],'  ',dominant_summary.(2)[i],'  ',dominant_summary.(3)[i],'  ',dominant_summary.(4)[i]
  endfor
  printf,c,'-------------------------------'
  printf,c,'Polygon Summary'
  printf,c,tag_names(polygon_summary)
  for i=0,total_polys-1 do begin
    
    printf,c,polygon_summary.(0)[i],' ',polygon_summary.(1)[i],'  ',polygon_summary.(2)[i],'  ',polygon_summary.(3)[i],'  ',polygon_summary.(4)[i],'  ',polygon_summary.(5)[i],'  ',polygon_summary.(6)[i],'  ',polygon_summary.(7)[i]
  endfor
  close, c
  free_lun, c

  ; load spectral libraries into ENVI memory
  envi_open_file, train_lib_spectra, r_fid=t_fid
  envi_file_query, t_fid, dims=dims, nb=nb, ns=ns, nl=nl, spec_names=spec_names, $
    interleave=interleave,byte_order=byte_order, bbl=bbl, wl=wl, data_type=data_type, $
    offset=offset, file_type=file_type

  envi_open_file, valid_lib_spectra, r_fid=v_fid
  envi_file_query, v_fid, dims=dims, nb=nb, ns=ns, nl=nl, spec_names=spec_names, $
    interleave=interleave,byte_order=byte_order, bbl=bbl, wl=wl, data_type=data_type, $
    offset=offset, file_type=file_type

  print, 'Libraries created.'
end