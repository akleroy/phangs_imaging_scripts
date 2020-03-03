pro inspect_residuals $
   , pause=pause $
   , start = start_num

;+
;
; Scripts to look at the residuals.
;
;-

; DIRECTORIES

  root_imaging_dir = '../'

  if n_elements(version) eq 0 then $
     version = '2'
  
  if version eq '2' then begin
     vstring = 'v2'
     release_dir = root_imaging_dir+'release/'+vstring+'/'
  endif else begin
     print, "Version not recognized. Returning."
     return
  endelse

  if n_elements(start_num) eq 0 then $
     start_num = 0 

; GALAXIES
  
; ... look up the list of galaxies

  readcol, 'ms_file_key.txt', comment='#', format='A,X,A' $
           , ms_file_gal, ms_file_array
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

; ARRAYS
  
  ;array_list = ['7m', '12m', '12m+7m']
  array_list = ['7m']
  n_array = n_elements(array_list)

; PRODUCTS
  
  product_list = ['co21']
  n_product = n_elements(product_list)

; SKIP

  only = []
  skip = []

; HAS A CENTRAL BRIGHT SOURCE

  bright_center = $
     ['ngc1097_1' $
      , 'ngc1097_2' $
      , 'ngc1300_1' $
      , 'ngc1300_2' $
      , 'ngc1317' $
      , 'ngc1365' $
      , 'ngc1433_1' $
      , 'ngc1433_2' $
      , 'ngc1512' $
      , 'ngc1566_1' $
      , 'ngc1566_2' $
      , 'ngc1637' $
      , 'ngc1672' $
      , 'ngc2566' $
      , 'ngc2903_1' $
      , 'ngc2903_2' $
      , 'ngc2903_3' $
      , 'ngc2997_1' $
      , 'ngc2997_2' $
      , 'ngc2997_3' $
      , 'ngc3351' $
      , 'ngc3507' $
      , 'ngc3626' $    
      , 'ngc3627_1' $
      , 'ngc3627_2' $
      , 'ngc4293' $
      , 'ngc4303' $
      , 'ngc4321_1' $
      , 'ngc4321_2' $
      , 'ngc4457' $
      , 'ngc4535' $
      , 'ngc4536_1' $
      , 'ngc4536_2' $
      , 'ngc4548' $
      , 'ngc4569' $      
      , 'ngc4579' $
      , 'ngc4826' $
      , 'ngc4941' $
      , 'ngc4951' $
      , 'ngc5128' $
      , 'ngc5248_1' $
      , 'ngc5248_2' $      
      , 'ngc5643_1' $
      , 'ngc5643_2' $
      , 'ngc6300' $
      , 'ngc7496' $
     ]


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = start_num, n_gals-1 do begin

     if n_elements(skip) gt 0 then begin
        if total(gals[ii] eq skip) gt 0 then begin
           print, "Skipping "+gals[ii]
           continue
        endif
     endif

     if n_elements(only) gt 0 then begin
        if total(gals[ii] eq only) eq 0 then begin
           print, "Skipping "+gals[ii]
           continue
        endif
     endif
     
     print, ""
     print, "Making clean mask for "+gals[ii]+' (Galaxy '+str(ii)+')'
     print, ""

     dir = release_dir+'process/'
     
     if n_elements(only) gt 0 then $
        if total(only eq gals[ii]) eq 0 then continue

     fname = '../'+dir_for_gal[ii]+'/'+gals[ii]+'_7m_co21_residual.fits'
     if file_test(fname) eq 0 then begin
        print, "No file found for "+gals[ii]
        print, "Continuing."
        continue
     endif
     
     cube = readfits(fname, hdr)

     mask = readfits('../'+dir_for_gal[ii]+'/'+gals[ii]+'_7m_co21_multiscale_mask.fits', mask_hdr)

     !p.multi=[0,3,2]

     rms = mad(cube)
     disp_fac = 5.0

     loadct, 33
     disp, max(cube, dim=3, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(mask,3,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(mask,3,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube, dim=2, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(mask,2,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(mask,2,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube, dim=1, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(mask,1,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(mask,1,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube*(mask eq 0), dim=3, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(mask,3,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(mask,3,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube*(mask eq 0), dim=2, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(mask,2,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(mask,2,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube*(mask eq 0), dim=1, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(mask,1,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(mask,1,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     if keyword_set(pause) then begin
        print, "Clean mask for "+gals[ii]+". Key to continue."
        test = get_kbrd(1)
     endif

  endfor


end
