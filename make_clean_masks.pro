pro make_clean_masks $
   , pause=pause $
   , inspect=do_inspect $
   , start = start_num

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
  
  array_list = ['7m', '12m', '12m+7m']
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
      , 'ngc1300' $
      , 'ngc1317' $
      , 'ngc1365' $
      , 'ngc1433' $
      , 'ngc1512' $
      , 'ngc1566' $
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
      , 'ngc3627north' $
      , 'ngc3627south' $
      , 'ngc4293' $
      , 'ngc4303' $
      , 'ngc4321north' $
      , 'ngc4321south' $
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

  empty = $
     ['ngc3239' $
     ]

  high_seed = $
     ['ngc3627north' $
      , 'ngc3627south' $
      , 'ngc4321north' $
      , 'ngc4321south' $
      , 'ngc4579']

  low_seed = $
     ['ic5332' $
      , 'ngc1317' $
      , 'ngc3239' $
      , 'ngc5042' $
      , 'ngc5068north' $
      , 'ngc5068south' $
      , 'ngc7456' $
     ]

  low_thresh = $
     ['ngc1317' $
      , 'ngc5042' $
      , 'ngc5068north' $
      , 'ngc5068south' $
      , 'ngc7456' $
     ]

  expand_in_space = $
     ['ngc1317' $
      , 'ngc5042' $
      , 'ngc5068north' $
      , 'ngc5068south' $
      , 'ngc7456']

  threed_noise = $
     ['ngc3627north' $
      , 'ngc3627south' $
      , 'ngc4321north' $
      , 'ngc4321south' $
      , 'ngc4579']

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
     
     clean_mask_fname = '../clean_masks/'+gals[ii]+'_co21_clean_mask.fits'

     if n_elements(only) gt 0 then $
        if total(only eq gals[ii]) eq 0 then continue

     fname = dir+gals[ii]+'_7m_co21_flat_round.fits'
     if file_test(fname) eq 0 then begin
        print, "No file found for "+gals[ii]
        print, "Continuing."
        continue
     endif
     
     cube = readfits(fname, hdr)

     if keyword_set(do_inspect) eq 0 then begin 

        ;if total(gals[ii] eq threed_noise) gt 0 then begin
        ;   print, "Three-d noise."
        ;   twod_only = 0B
        ;endif else begin
        ;   twod_only = 1B
        ;endelse

        make_noise_cube $
           , cube_in = cube $
           , out_cube = rms_cube $
           , twod_only=twod_only $
           , show=0B $
           , /iterate

        ppbeam = calc_pixperbeam(hdr=hdr)

        fac = 1.0

        hi_thresh = 5.0

        if total(gals[ii] eq high_seed) gt 0 then begin
           print, "High seed"
           hi_thresh = 10.0
        endif

        if total(gals[ii] eq low_seed) gt 0 then begin
           print, "Low seed"
           hi_thresh = 5.0
           if gals[ii] eq 'ic5332' then hi_thresh = 7.0
           if gals[ii] eq 'ngc5042' then hi_thresh = 4.0
           if gals[ii] eq 'ngc7456' then hi_thresh = 3.5
        endif

        lo_thresh = 2        

        if total(gals[ii] eq low_thresh) gt 0 then begin
           print, "Low threshold"
           lo_thresh = 1.5
        endif

        make_cprops_mask $
           , indata=cube $
           , inrms = rms_cube $
           , lo_thresh = lo_thresh $
           , hi_thresh = hi_thresh $
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

        ; LOW RESOLUTION VERSION

        conv_with_gauss $
           , data=cube $
           , hdr=hdr $
           , target = [1,1,0.] * 30.0 $
           , out_data= lowres_cube $
           , out_hdr = lowres_hdr
        
        make_noise_cube $
           , cube_in = lowres_cube $
           , out_cube = rms_lowres_cube $
           , twod_only=twod_only $
           , show=0B $
           , /iterate

        hi_thresh_lowres = 10.0

        if total(gals[ii] eq high_seed) gt 0 then begin
           print, "High seed"
           hi_thresh_lowres = 20.0
        endif

        if total(gals[ii] eq low_seed) gt 0 then begin
           print, "Low seed"
           hi_thresh_lowres = 5.0
           if gals[ii] eq 'ic5332' then hi_thresh_lowres = 7.0
           if gals[ii] eq 'ngc5042' then hi_thresh_lowres = 4.0
           if gals[ii] eq 'ngc7456' then hi_thresh_lowres = 3.5
        endif

        lo_thresh_lowres = 3.0
        if total(gals[ii] eq low_thresh) gt 0 then begin
           print, "Low threshold"
           lo_thresh_lowres = 2.0
           if gals[ii] eq 'ngc5042' then lo_thresh_lowres = 1.5
           if gals[ii] eq 'ngc7456' then lo_thresh_lowres = 1.5
        endif        

        lowres_ppbeam = calc_pixperbeam(hdr=lowres_hdr)
        make_cprops_mask $
           , indata=lowres_cube $
           , inrms = rms_lowres_cube $
           , lo_thresh = lo_thresh_lowres $
           , hi_thresh = hi_thresh_lowres $
           , hi_nchan = 2 $
           , min_area = 3.0*lowres_ppbeam $
           , min_pix = 6.0*lowres_ppbeam $
           , outmask=lowres_mask

        conv_with_gauss $
           , data=lowres_mask*1.0 $
           , hdr=lowres_hdr $
           , target = [1,1,0.] * 33.0 $
           , out_data=out_lowres_mask $
           , /perbeam        
        if total(gals[ii] eq expand_in_space) gt 0 then begin
           lowres_mask = out_lowres_mask ge 0.5
        endif else begin
           lowres_mask = out_lowres_mask ge 1.0
        endelse

        mask = (mask + lowres_mask) ge 1

        if gals[ii] eq 'ngc7456' then begin
           mask[*,*,0:50] = 0B
        endif

        mask = grow_mask(mask, iters=5, /z_only)
        if total(gals[ii] eq bright_center) gt 0 then begin
           mask = grow_mask(mask, iters=5, /z_only)
        endif

        if total(gals[ii] eq empty) gt 0 then begin
           print, "Cube looks empty - setting mask to unity."
           mask = finite(mask)*1B
        endif

     endif else begin
        
        if file_test(clean_mask_fname) eq 0 then begin
           print, "Inspection mode but mask not found : ", clean_mask_fname
        endif
        mask = readfits(clean_mask_fname, hdr)

     endelse

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

     sxaddpar, hdr, 'BUNIT', 'MASK'
     writefits $
        , dir+gals[ii]+'_co21_clean_mask.fits' $
        , float(mask*1.0), hdr
     
     if keyword_set(pause) then begin
        print, "Clean mask for "+gals[ii]+". Key to continue."
        test = get_kbrd(1)
     endif

     writefits, clean_mask_fname $
                , mask, hdr

  endfor


end
