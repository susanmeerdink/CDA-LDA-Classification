pro strip_redundant_spectra_v6
  COMPILE_OPT IDL2
  
  ;Last update by Keely Roth on 8/12/10
  ;Update by Kenneth Dudley on 2/22/14 Add support for very large spectral libraries
  ;Update by Kenneth Dudley on 4/19/14 to change csv metadata read in to use read_csv which prevents quotation "" issues
  ; Streamlined calculation of brightness values and concatenating the values to metadata array
  ; Removed scaling, which was not being used at all
  ; This program changes the method from comapring neighbors, to comparing each spectra with all other spectra "below" it.
  ;   This method will find multiple duplicates for any given spectra. Even if they are not adjacent. Spectra are sorted by the
  ;   total of each spectra, such that if the total is not equal, we prevent checking the remaining spectra, this is much faster.
  ; Sorting is reverted to maintain the basic sort scheme of the original input dataset.
  ; *************************************************************************************
  ; INPUT VARIABLES
  close, /all ; close open file units
  
  ; spectral library name
  sli_file_name =envi_pickfile(title='Select a Spectral Library', filter='*.sli')
  root = STRMID(sli_file_name,0, STRLEN(sli_file_name)-4)
  out_lib = root[0]+'_nodupev6.sli'
  
  ;**************************************************************************************
  
  ; open the spectral library and read in data
  envi_open_file, sli_file_name, r_fid=sli_fid
  envi_file_query, sli_fid, dims=dims, nb=nb, ns=ns, nl=nl, spec_names=spec_names, $
    interleave=interleave,byte_order=byte_order, bbl=bbl, wl=wl, data_type=data_type, $
    offset=offset, file_type=file_type
  speclib = envi_get_data(fid=sli_fid, dims=dims, pos=[0])
  
  total_libs = (total(speclib,1))
  sort_spec_lib = sort(total_libs) ; Sort by totals of each spectral library
  revert_sort = sort(sort_spec_lib) ; Use the mapping vector to re-sort data back into original positions
  checksort = total_libs[sort_spec_lib]
  speclib_find = speclib[*,[sort_spec_lib]]
  
  ;Begin testing spectra for duplicates
  counter = 1L
  tracker = bytarr(nl)
  totest = fltarr(ns)
  for i=0L, nl-1L do begin
    totest[0] = speclib_find[*,i] ; Spectra to compare to remaining spectra
    for j=counter, nl-1L do begin
      if checksort[i] EQ checksort[j] then begin ; Check if totals are equal, when they are not we can break this loop
        result = array_equal(totest,speclib_find[*,j], /NO_TYPECONV) ; Compare spectra to all spectra below it
        if result EQ 1 then begin
          tracker[i]=byte(1) ; track location of duplicate spectra
          break ; if a duplicate is found stop looking in this loop
        endif
      endif else break
    endfor
    counter++
  endfor
  
  tracker = tracker[revert_sort] ; Reverts the sort to re-index based on original input locations
  are_there_duplicates = where(tracker EQ 1, /NULL, COMPLEMENT=non_duplicate_loc) ; finds duplicate and non-duplicate index locations
  new_nl = n_elements(non_duplicate_loc)
  if are_there_duplicates EQ !NULL then begin
    print, "There are no duplicate spectra in this library, quitting..."
    return
  endif
  
  ;Read in metadata file
  meta_data_file_name = root[0] + '.csv'
  header_csv = read_csv(meta_data_file_name, NUM_RECORDS=0, RECORD_START=0) ; Read in header separately, since read_csv will not find the header if next line does not contain numerics
  header = strarr(n_tags(header_csv))
  header = strtrim(header,2)
  for i=0, n_tags(header_csv)-1 do header[i] = header_csv.(i)
  
  metadata_csvread = read_csv(meta_data_file_name, RECORD_START=1) ; Start reading at line 1, after line 0 (the header).
  ncol = n_tags(metadata_csvread)
  meta_array=strarr(ncol,nl)
  for i=0, ncol-1 do meta_array[i,*] = metadata_csvread.(i) ; convert metadata file to a string array
  meta_array = strtrim(meta_array,2) ; Remove extra spaces carried over from metadata (if they exist)
  no_dupe_meta = meta_array[*,non_duplicate_loc] ; subset metadata array to only those elements which are not duplicates
  
  ;Name for the new output library minus duplicates
  root = strsplit(out_lib,'.', /extract)
  new_meta_data_file_name = root[0] + '.csv'
  
  ;  write_csv, new_meta_data_file_name, no_dupe_meta, HEADER=header ; write new output csv file
  
  ;Alternative to write_csv
  openw, u1, new_meta_data_file_name, /get_lun
  printf, u1, strjoin(header,',',/single)
  for i=0L,new_nl-1 do printf,u1, strjoin(no_dupe_meta[*,i],',',/single)
  close,u1
  free_lun,u1
  
  new_head_name = root[0] + '.hdr'
  new_spec_names = spec_names[non_duplicate_loc]
  new_spec_lib = speclib[*,non_duplicate_loc]
  
  openw, 1, out_lib
  writeu, 1, new_spec_lib
  close, 1
  
  envi_setup_head, fname=new_head_name,ns=ns, nl=new_nl, nb=nb, spec_names=new_spec_names, data_type=data_type,$
    bbl=bbl,interleave=interleave, wl=wl, offset=offset, file_type=file_type, /write
    
  print, "Input #: ", strtrim(nl,2)
  print, "Output #: ", strtrim(new_nl,2)
  print, "# Removed: ", strtrim(nl-new_nl,2)
end