pro compare_parts_to_mosaic $
   , old_maps = old_maps

; DIRECTORIES

  root_imaging_dir = '../'

  if n_elements(version) eq 0 then $
     version = '3'
  
  if version eq '3' then begin
     vstring = 'v3'
     release_dir = root_imaging_dir+'release/'+vstring+'/'
  endif else begin
     print, "Version not recognized. Returning."
     return
  endelse

  process_dir = release_dir+'process/'

  old_dir = '../release/v3/delivery/broad_maps/'

; GALAXIES
  
; ... look up the list of galaxies

  readcol, 'ms_file_key.txt', comment='#', format='A,X,A' $
           , ms_file_gal, ms_file_array, /silent
  gals = ms_file_gal[sort(ms_file_gal)]
  gals = gals[uniq(gals)]
  n_gals = n_elements(gals)
  
; ... look up the cases with nonstandard directory names

  readcol, 'dir_key.txt', comment='#', format='A,A' $
           , dir_key_gal, dir_key_dir
  dir_for_gal = gals
  for ii = 0, n_elements(dir_key_gal)-1 do begin
     ind = where(dir_key_gal[ii] eq gals, ct)
     if ct eq 0 then continue
     dir_for_gal[ind] = (dir_key_dir[ii])
  endfor

  fullgal_list = [dir_for_gal, gals]
  fullgal_list = fullgal_list[sort(fullgal_list)]
  fullgal_list = fullgal_list[uniq(fullgal_list)]
  n_fullgals = n_elements(fullgal_list)

  gal_for_fullgals = fullgal_list
  for ii = 0, n_elements(dir_key_gal)-1 do begin
     ind = where(dir_key_gal[ii] eq fullgal_list, ct)
     if ct eq 0 then continue
     gal_for_fullgals[ind] = (dir_key_dir[ii])
  endfor

  multipart_gals = dir_for_gal[uniq(dir_for_gal, sort(dir_for_gal))]
  multipart_gals = multipart_gals[sort(multipart_gals)]

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER MULTI-PART GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(old_maps) then begin
     this_dir = old_dir
  endif else begin
     this_dir = process_dir
  endelse

  n_mpg = n_elements(multipart_gals)
  for ii = 0, n_mpg-1 do begin
     this_gal = (multipart_gals[ii])[0]
     this_ind = where(dir_key_dir eq this_gal, this_ct)
     if this_ct eq 0 then continue
     parts = dir_key_gal[this_ind]
     n_parts = n_elements(parts)
     print, this_gal, " : ", parts

     infile = this_dir+this_gal+'_7m+tp_co21_broad_mom0.fits'
     einfile = this_dir+this_gal+'_7m+tp_co21_broad_emom0.fits'
     if file_test(infile) eq 0 then begin
        print, infile, ' not found.'
        continue
     endif

     mosaic = readfits(infile, mos_hdr)
     emosaic = readfits(einfile, emos_hdr)
     target_res = sxpar(mos_hdr, 'BMAJ')*3600.

     for jj = 0, n_parts-1 do begin
        this_part = parts[jj]
        infile = process_dir+this_part+'_7m+tp_co21_broad_mom0.fits'
        einfile = process_dir+this_part+'_7m+tp_co21_broad_emom0.fits'
        if file_test(infile) eq 0 then begin
           print, infile, ' not found.'
           continue
        endif
        part = readfits(infile, part_hdr)
        epart = readfits(einfile, epart_hdr)

        if sxpar(part_hdr,'BMAJ')*1.02*3600. lt target_res then begin
           conv_with_gauss $
              , data=part, hdr=part_hdr $
              , target_beam=target_res*[1,1,0] $
              , out_data=smooth_part, out_hdr=smooth_part_hdr $
              , worked=worked
        endif else begin
           print,""
           print,"CONVOLUTION DID NOT WORK. COPYING."
           print,""
           smooth_part = part
           smooth_part_hdr = part_hdr
        endelse
        hastrom, smooth_part, smooth_part_hdr, mos_hdr $
                 , interp=2, cubic=-0.5 $
                 , missing=!values.f_nan
        hastrom, epart, epart_hdr, mos_hdr $
                 , interp=2, cubic=-0.5 $
                 , missing=!values.f_nan

        overlap = where(finite(smooth_part) and (smooth_part ne 0) $
                        and finite(mosaic) and (mosaic ne 0), overlap_ct)
        if overlap_ct eq 0 then continue
        
        if n_elements(x) eq 0 then begin
           g = replicate(this_gal, overlap_ct)
           x = smooth_part[overlap]
           ex = epart[overlap]
           y = mosaic[overlap]
           ey = emosaic[overlap]
        endif else begin
           g = [g, replicate(this_gal, overlap_ct)]
           x = [x,smooth_part[overlap]]
           ex = [ex,epart[overlap]]
           y = [y, mosaic[overlap]]
           ey = [ey, emosaic[overlap]]
        endelse
        
        plot, x, y, ps=3 $
              , xrange=[1,1d3], yrange=[1,1d3], /xlo, /ylo
     endfor

  endfor

  ind = where(y gt 1 and x gt 1, ct)  
  print, median((y/x)[ind])
  print, mad((y/x)[ind])
  fasthist, alog10((y/x)[ind])

  ind = where(y gt 5 and x gt 5, ct)  
  print, median((y/x)[ind])
  print, mad((y/x)[ind])
  fasthist, alog10((y/x)[ind])

  n_mpg = n_elements(multipart_gals)
  for ii = 0, n_mpg-1 do begin
     this_gal = (multipart_gals[ii])[0]
     ind = where(g eq this_gal and x gt 1 and y gt 1, ct)
     if ct eq 0 then continue
     print, this_gal, median((y/x)[ind]), mad(alog10((y/x)[ind]))
  endfor

  stop

end
