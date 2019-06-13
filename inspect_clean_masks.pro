pro inspect_clean_masks $
   , version=version $
   , nopause=nopause $
   , inspect=do_inspect $
   , start = start_num $
   , array = array $
   , skip = skip $
   , only = only

;+
;
; Create clean masks for use with imaging.
;
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DIRECTORIES AND SETUP
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  root_imaging_dir = '../'

; ... look up the version to based the masks on

  if n_elements(version) eq 0 then $
     version = '3'
  
  if version eq '2' then begin
     vstring = 'v2'
     release_dir = root_imaging_dir+'release/'+vstring+'/'
  endif else if version eq '3' then begin
     vstring = 'v3'
     release_dir = root_imaging_dir+'release/'+vstring+'/'
  endif else begin
     print, "Version not recognized. Returning."
     return
  endelse

  if n_elements(start_num) eq 0 then $
     start_num = 0 

; ... look up the list of galaxies

  readcol, 'ms_file_key.txt', comment='#', format='A,X,A' $
           , ms_file_gal, ms_file_array, /silent
  gals = ms_file_gal[sort(ms_file_gal)]
  gals = gals[uniq(gals)]
  n_gals = n_elements(gals)
  
; ... look up the cases with nonstandard directory names

  readcol, 'dir_key.txt', comment='#', format='A,A' $
           , dir_key_gal, dir_key_dir, /silent
  dir_for_gal = gals
  for ii = 0, n_elements(dir_key_gal)-1 do begin
     ind = where(dir_key_gal[ii] eq gals, ct)
     if ct eq 0 then continue
     dir_for_gal[ind] = (dir_key_dir[ii])
  endfor

; ... define array to use (default to 7m)

  if n_elements(array) eq 0 then $
     array = '7m'

; ... define products
  
  product_list = ['co21']
  n_product = n_elements(product_list)

; ... lists of galaxies to skip or focus on

  if n_elements(only) eq 0 then $
     only = []
  if n_elements(skip) eq 0 then $
     skip = []

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = start_num, n_gals-1 do begin

     dir_ind = where(dir_key_gal eq gals[ii], dir_ct)
     is_multipart = dir_ct ge 1
     if is_multipart then begin
        this_gal = dir_key_dir[dir_ind]
     endif else begin
        this_gal = gals[ii]
     endelse

     if n_elements(skip) gt 0 then begin
        if total(this_gal eq skip) gt 0 or $
           total(gals[ii] eq skip) gt 0 then begin
           ;print, "Skipping "+gals[ii]
           continue
        endif
     endif

     if n_elements(only) gt 0 then begin
        if (total(gals[ii] eq only) eq 0) and $
           (total(this_gal eq only) eq 0) $
        then begin
           ;print, "Skipping "+gals[ii]
           continue
        endif
     endif
     
     print, ""
     print, "Inspecting clean mask for "+gals[ii]+' (Galaxy '+str(ii)+')'
     print, ""

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; VERIFY THE CLEAN MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     clean_mask_fname = '../clean_masks/'+this_gal+'_co21_clean_mask.fits'

     if file_test(clean_mask_fname) eq 0 then begin
        print, "Missing clean mask: ", clean_mask_fname
        print, "... continuing to next target."
        continue
     endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE DATA CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     dir = release_dir+'process/'

     cube_fname = dir+this_gal+'_'+array+'+tp_co21_flat_round_k.fits'
     
     if file_test(cube_fname) eq 0 then begin
        cube_fname = dir+this_gal+'_'+array+'_co21_flat_round_k.fits'
     endif

     if file_test(cube_fname) eq 0 then begin
        print, "Missing data cube: ", cube_fname
        print, "... continuing to next target."
        continue
     endif     
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ AND ALIGN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     mask = readfits(clean_mask_fname, mask_hdr)
     
     cube = readfits(cube_fname, cube_hdr)

     sz_mask = size(mask)
     sz_cube = size(cube)

     if sz_mask[3] ne sz_cube[3] then begin
        print, "Spectral axes disagree. Stopping (for now)."
        print, "Implement regridding here if needed."
        stop
     endif

     if sz_mask[1] ne sz_cube[1] or $
        sz_mask[2] ne sz_cube[2] then begin       
        cube_hastrom $
           , data = mask $
           , hdr_in = mask_hdr $
           , outcube = aligned_mask $
           , target_hdr = cube_hdr $   
           , pinterp=0 $
           , pmissing=0 $
           , operation='POS'
     endif else begin
        aligned_mask = mask
     endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHOW
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
     
     !p.multi=[0,3,2]

     disp_cube = cube
     rms = mad(disp_cube)
     disp_fac = 5.0

     loadct, 33
     disp, max(disp_cube, dim=3, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(aligned_mask,3,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(aligned_mask,3,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(disp_cube, dim=2, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(aligned_mask,2,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(aligned_mask,2,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(disp_cube, dim=1, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(aligned_mask,1,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(aligned_mask,1,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(disp_cube*(aligned_mask eq 0), dim=3, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(aligned_mask,3,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(aligned_mask,3,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(disp_cube*(aligned_mask eq 0), dim=2, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(aligned_mask,2,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(aligned_mask,2,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(disp_cube*(aligned_mask eq 0), dim=1, /nan), /xs, /ys, max=rms*disp_fac
     contour, total(aligned_mask,1,/nan) gt 0, /overplot, lev=[1], thick=5, color=cgcolor('black')
     contour, total(aligned_mask,1,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')
     
     if keyword_set(nopause) eq 0 then begin
        print, "Clean mask for "+gals[ii]+". Key to continue."
        test = get_kbrd(1)
     endif

  endfor


end
