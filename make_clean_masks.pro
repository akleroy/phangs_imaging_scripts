pro make_clean_masks $
   , pause=pause

;+
;
; Scripts to build the imaged data into data cubes suitable for
; analysis.
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
  
  array_list = ['7m', '12m', '12m+7m']
  n_array = n_elements(array_list)

; PRODUCTS
  
  product_list = ['co21']
  n_product = n_elements(product_list)

; SKIP

  only = ['ngc1097','ngc4579','ngc1097_1','ngc1097_2']
  skip = ['ngc1365']

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, n_gals-1 do begin

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

     message, "Making clean mask for "+gals[ii], /info

     dir = release_dir+'process/'
     
     if n_elements(only) gt 0 then $
        if total(only eq gals[ii]) eq 0 then continue

;     fname = dir+gals[ii]+'_7m+tp_co21_flat_round_k.fits'
;     if file_test(fname) eq 0 then $

     fname = dir+gals[ii]+'_7m_co21_flat_round_k.fits'
     if file_test(fname) eq 0 then begin
        print, "No file found for "+gals[ii]
        print, "Continuing."
        continue
     endif
        
     cube = readfits(fname, hdr)
     
     hires_name = file_search(dir+gals[ii]+'_co21_12m+7m+tp_flat_round_k.fits', count=hires_ct)
     if hires_ct eq 1 then begin
        hires = readfits(hires_name, hires_hdr)
        rms_hires = mad(hires)
        ppbeam = calc_pixperbeam(hdr=hires_hdr)
        make_cprops_mask $
           , indata=hires $
           , inrms = rms_hires $
           , lo_thresh = 2 $
           , hi_thresh = 5. $
           , hi_nchan = 2 $
           , min_area = 1.0*ppbeam $
           , outmask=hires_mask
        conv_with_gauss $
           , data=hires_mask*1.0 $
           , hdr=hires_hdr $
           , target = [1,1,0.] * 20.0 $
           , out_data=out_hires_mask $
           , /perbeam        
        hires_mask = out_hires_mask ge 0.5          
        cube_hastrom $
           , data=hires_mask $
           , hdr=hires_hdr $
           , target=hdr $
           , outcube=out_hires_mask
        hires_mask = out_hires_mask gt 0.5
     endif

     make_noise_cube $
        , cube_in = cube $
        , out_cube = rms_cube $
        , /twod_only $
        , show=0B $
        , /iterate

     ppbeam = calc_pixperbeam(hdr=hdr)

     fac = 1.0
     if gals[ii] eq 'ngc1365' then $
        fac = 0.25
     if gals[ii] eq 'ngc1566' then $
        fac = 0.25
     if gals[ii] eq 'ngc1672' then $
        fac = 0.25

     make_cprops_mask $
        , indata=cube $
        , inrms = rms_cube $
        , lo_thresh = 2 $
        , hi_thresh = 5. $
        , hi_nchan = 2 $
        , min_area = fac*ppbeam $
        , outmask=mask

     conv_with_gauss $
        , data=mask*1.0 $
        , hdr=hdr $
        , target = [1,1,0.] * 20.0 $
        , out_data=out_mask $
        , /perbeam        
     mask = out_mask ge 0.5

     conv_with_gauss $
        , data=cube $
        , hdr=hdr $
        , target = [1,1,0.] * 30.0 $
        , out_data= tp_cube $
        , out_hdr = tp_hdr
     
     make_noise_cube $
        , cube_in = tp_cube $
        , out_cube = rms_tp_cube $
        , /twod_only $
        , show=0B $
        , /iterate

     tp_thresh = 10.
     if gals[ii] eq 'ic5332' then $
        tp_thresh = 7.
     
     tp_ppbeam = calc_pixperbeam(hdr=tp_hdr)
     make_cprops_mask $
        , indata=tp_cube $
        , inrms = rms_tp_cube $
        , lo_thresh = 3 $
        , hi_thresh = tp_thresh $
        , hi_nchan = 2 $
        , min_area = 3.0*tp_ppbeam $
        , min_pix = 6.0*tp_ppbeam $
        , outmask=tp_mask

     conv_with_gauss $
        , data=tp_mask*1.0 $
        , hdr=tp_hdr $
        , target = [1,1,0.] * 33.0 $
        , out_data=out_tp_mask $
        , /perbeam        
     tp_mask = out_tp_mask ge 1.0

     if hires_ct eq 1 then $
        mask = (hires_mask + mask + tp_mask) ge 1 $
     else $
        mask = (mask + tp_mask) ge 1

     mask = grow_mask(mask, iters=3, /z_only)
     if gals[ii] eq 'ngc1365' then $
        mask = grow_mask(mask, iters=5, /z_only)
     if gals[ii] eq 'ngc1672' then $
        mask = grow_mask(mask, iters=10, /z_only)
     if gals[ii] eq 'ngc3627' then $
        mask = grow_mask(mask, iters=5, /z_only)
     mask = grow_mask(mask, iters=5, /xy_only)

     !p.multi=[0,3,2]

     loadct, 33
     disp, max(cube, dim=3, /nan), /xs, /ys, max=0.1
     contour, total(mask,3,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube, dim=2, /nan), /xs, /ys, max=0.1
     contour, total(mask,2,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube, dim=1, /nan), /xs, /ys, max=0.1
     contour, total(mask,1,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube*(mask eq 0), dim=3, /nan), /xs, /ys
     contour, total(mask,3,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube*(mask eq 0), dim=2, /nan), /xs, /ys
     contour, total(mask,2,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     loadct, 33
     disp, max(cube*(mask eq 0), dim=1, /nan), /xs, /ys
     contour, total(mask,1,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

     sxaddpar, hdr, 'BUNIT', 'MASK'
     writefits $
        , dir+gals[ii]+'_co21_clean_mask.fits' $
        , float(mask*1.0), hdr
     
     if keyword_set(pause) then begin
        print, "Clean mask for "+gals[ii]+". Key to continue."
        test = get_kbrd(1)
     endif

  endfor


end
