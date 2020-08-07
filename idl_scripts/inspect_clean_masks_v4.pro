pro inspect_clean_masks_v4 $
   , version=version $
   , nopause=nopause $
   , existence=existence $
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

; ... look up the version to based the masks on

  imroot_dir = '../../phangs-alma/imaging/'
  pproot_dir = '../../phangs-alma/postprocess/'
  cleanmask_dir = '../../phangs-alma/cleanmasks/'     

  if n_elements(start_num) eq 0 then $
     start_num = 0 


; ... look up the list of galaxies

  readcol, '../phangsalma_keys/ms_file_key.txt', comment='#', format='A,X,A' $
           , ms_file_gal, ms_file_array, /silent

  gals = ms_file_gal[sort(ms_file_gal)]
  gals = gals[uniq(gals)]
  n_gals = n_elements(gals)
  
; ... look up the cases with nonstandard directory names

  readcol, '../phangsalma_keys/dir_key.txt', comment='#', format='A,A' $
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
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; VERIFY THE CLEAN MASK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     clean_mask_fname = cleanmask_dir+this_gal+'_co21_clean_mask.fits'

     if file_test(clean_mask_fname) eq 0 then begin
        print, "Missing clean mask: ", clean_mask_fname
        print, "... continuing to next target."
        continue
     endif

     if keyword_set(existence) then begin
        continue
     endif

     print, ""
     print, "Inspecting clean mask for "+gals[ii]+' (Galaxy '+str(ii)+')'
     print, ""

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ THE DATA CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

     dir = pproot_dir

     cube_fname = $
        dir+this_gal+'/'+ $
        this_gal+'_'+array+'+tp_co21_pbcorr_trimmed_k.fits'
     
     if file_test(cube_fname) eq 0 then begin
        cube_fname = $
           dir+this_gal+'/'+ $
           this_gal+'_'+array+'_co21_pbcorr_trimmed_k.fits'
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

     cube_hastrom $
        , data = mask $
        , hdr_in = mask_hdr $
        , outcube = aligned_mask $
        , target_hdr = cube_hdr $   
        , pinterp=0 $
        , pmissing=0 $
        , operation='BOT'

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
