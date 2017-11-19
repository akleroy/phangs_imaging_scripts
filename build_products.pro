pro build_products $
   , version=version $
   , only=only $
   , skip=skip $
   , just_array=just_array $
   , noise=do_noise $
   , mask=do_masks $
   , collapse=do_collapse $
   , target_res=target_res

;+
;
; Scripts to build the data cubes into maps and products for release.
;
;-

; THRESHOLDS FOR MOMENT CREATION

  mom1_thresh = 2.5d
  mom0_thresh = 3.0d

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

; ARRAYS
  
  array_list = ['7m', '12m', '12m+7m']
  n_array = n_elements(array_list)

  fullarray_list = ['7m', '7m+tp', '12m', '12m+tp', '12m+7m', '12m+7m+tp']
  n_fullarray = n_elements(fullarray_list)

; PRODUCTS
  
  product_list = ['co21']
  n_product = n_elements(product_list)

; RESOLUTIONS

  if n_elements(target_res) eq 0 then begin
     target_res = [-1, 45, 60, 80, 100, 120, 500, 750, 1000]
  endif
  n_res = n_elements(target_res)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ESTIMATE THE NOISE FOR EACH CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_noise) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'ESTIMATE NOISE', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     ext_to_process = $
        ['flat_round_k' $
         , 'pbcorr_round_k']
     n_ext = n_elements(ext_to_process)
     
     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        message, "Estimating noise for for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_process[ll]

                 for zz = 0, n_res -1 do begin
                    
                    if target_res[zz] eq -1 then begin
                       res_str = ''
                       print, "Native resolution."
                    endif else begin
                       res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                       print, "Resolution "+res_str
                    endelse

                    in_file = release_dir+'process/'+ $
                              this_gal+'_'+this_array+'_'+ $
                              this_product+'_'+this_ext+res_str+'.fits'
                    
                    test = file_search(in_file, count=found)
                    if found eq 0 then begin
                       message, 'File '+in_file+' not found.', /info
                       continue
                    endif                 
                    
                    cube = readfits(in_file, hdr)
                    
                    make_noise_cube $
                       , cube_in = cube $
                       , out_cube = rms_cube $
                       , box = 5 $
                       , /twod_only $
                       , /show $
                       , /iterate
                    
                    out_file = $
                       release_dir+'process/'+ $
                       this_gal+'_'+this_array+'_'+ $
                       this_product+'_'+'noise_'+ $
                       this_ext+res_str+'.fits'
                    writefits, out_file, rms_cube, cube_hdr                    

                 endfor
                 
              endfor

           endfor

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_masks) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'BUILD MASKS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        message, "Estimating noise for for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; AT EACH RESOLUTION MAKE A MASK HOLDING BRIGHT SIGNAL
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

              for zz = 0, n_res -1 do begin
                 
                 if target_res[zz] eq -1 then begin
                    res_str = ''
                    print, "Native resolution."
                 endif else begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                    print, "Resolution "+res_str
                 endelse
                 
                 noise_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_'+'noise_flat_round_k'+ $
                    res_str+'.fits'

                 test = file_search(noise_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+noise_file+' not found.', /info
                    continue
                 endif                 
                 rms_cube = readfits(noise_file, noise_hdr)
                 
                 cube_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_flat_round_k'+ $
                    res_str+'.fits'
                 
                 if found eq 0 then begin
                    message, 'File '+cube_file+' not found.', /info
                    continue
                 endif                 
                 cube = readfits(cube_file, cube_hdr)
                 
                 ppbeam = calc_pixperbeam(hdr=cube_hdr)         
                 
                 make_cprops_mask $
                    , indata=cube $
                    , inrms = rms_cube $
                    , lo_thresh = 2 $
                    , hi_thresh = 3.5 $
                    , hi_nchan = 3 $
                    , min_area = ppbeam $
                    , outmask=mask

                 !p.multi=[0,2,1]
                 viridis
                 disp, max(cube, dim=3, /nan)
                 contour, max(mask, dim=3), lev=[1], /overplot
                 
                 disp, max(cube, dim=2, /nan)
                 contour, max(mask, dim=2), lev=[1], /overplot
                 
                 mask_hdr = cube_hdr
                 sxaddpar, mask_hdr, 'BUNIT', 'MASK'

                 out_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_signalmask'+ $
                    res_str+'.fits'

                 writefits, out_file $
                            , mask, mask_hdr
                 
              endfor

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; HYBRIDIZE THE MASKS, COMBINING BRIGHT SIGNAL AND BROAD LOW RES
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

              lores_fname = $
                 release_dir+'process/'+ $
                 this_gal+'_'+this_array+'_'+ $
                 this_product+ $
                 '_signalmask_500pc.fits'

              test = file_search(lores_fname, count=found)
              if found eq 0 then begin
                 message, 'File '+lores_fname+' not found.', /info
                 continue
              endif                 
              mask_lores = readfits(lores_fname)

              for zz = 0, n_res do begin
                 
                 if zz lt n_res then begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                    print, "Resolution "+res_str
                 endif else begin
                    res_str = ''
                    print, "Native resolution"
                 endelse
                 
                 hires_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_signalmask'+ $
                    res_str+'.fits'
                 
                 test = file_search(hires_fname, count=found)
                 if found eq 0 then begin
                    message, 'File '+hires_fname+' not found.', /info
                    continue
                 endif                 
                 mask_hires = readfits(hires_fname, mask_hdr)
                 
                 hybrid = mask_hires or mask_lores

                 hybrid_fname= $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_hybridmask'+ $
                    res_str+'.fits'
                 
                 writefits, hybrid_fname $
                            , hybrid, mask_hdr                 

              endfor

           endfor

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COLLAPSE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_collapse) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COLLAPSE INTO MOMENTS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     dir = release_dir+'process/'
     
     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        message, "Estimating noise for for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; AT EACH RESOLUTION MAKE A MASK HOLDING BRIGHT SIGNAL
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

              for zz = 0, n_res -1 do begin
                 
                 if target_res[zz] eq -1 then begin
                    res_str = ''
                    print, "Native resolution."
                 endif else begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                    print, "Resolution "+res_str
                 endelse
                 
                 noise_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_'+'noise_flat_round_k'+ $
                    res_str+'.fits'

                 test = file_search(noise_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+noise_file+' not found.', /info
                    continue
                 endif                 
                 rms_cube = readfits(noise_file, noise_hdr)
                 
                 cube_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_pbcorr_round_k'+ $
                    res_str+'.fits'
                 
                 test = file_search(cube_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+cube_file+' not found.', /info
                    continue
                 endif                 
                 cube = readfits(cube_file, cube_hdr)
                 
                 mask_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_hybridmask'+ $
                    res_str+'.fits'

                 test = file_search(mask_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+mask_file+' not found.', /info
                    continue
                 endif                 
                 mask = readfits(mask_file, mask_hdr)
                 
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; CALCULATE MOMENTS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 collapse_cube $
                    , cube=cube $
                    , hdr=cube_hdr $
                    , mask=mask $
                    , noise=rms_cube $
                    , mom0 = mom0 $
                    , e_mom0 = e_mom0 $
                    , mom1 = mom1 $
                    , e_mom1 = e_mom1 $
                    , mom2 = mom2 $
                    , e_mom2 = e_mom2 $
                    , ew = ew $
                    , e_ew = e_ew $
                    , var = var $
                    , e_var = e_var $
                    , tpeak = tpeak

                 blank_ind = where(total(finite(cube),3) eq 0)

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PEAK TEMPERATURE IN ONE AND FIVE CHANNELS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 tpeak_mask = mask
                 sz = size(mask)
                 for pp = 0, sz[3]-1 do $
                    tpeak_mask[*,*,pp] = total(mask[*,*,pp]) ge 1

                 tpeak = max(cube*tpeak_mask, dim=3, /nan)
                 tpeak[blank_ind] = !values.f_nan
                 tpeak_hdr = twod_head(cube_hdr)
                 sxaddpar, tpeak_hdr, 'BUNIT', 'K'
                 
                 tpeak_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_tpeak'+ $
                    res_str+'.fits'

                 writefits, tpeak_fname, tpeak, tpeak_hdr
                 
                 tpeak_12p5 = max(smooth(cube,[1,1,5],/nan,/edge_wrap)*tpeak_mask, dim=3, /nan)
                 tpeak_12p5[blank_ind] = !values.f_nan
                 tpeak_hdr = twod_head(cube_hdr)
                 sxaddpar, tpeak_hdr, 'BUNIT', 'K'

                 tpeak12p5_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_tpeak'+ $
                    res_str+'.fits'

                 writefits, tpeak12p5_fname, tpeak_12p5, tpeak_hdr

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; MOMENT 0
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 mom0_hdr = twod_head(cube_hdr)
                 sxaddpar, mom0_hdr, 'BUNIT', 'K*KM/S'
                 mom0[blank_ind] = !values.f_nan
                 
                 mom0_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_mom0'+ $
                    res_str+'.fits'

                 writefits, mom0_fname $
                            , mom0, mom0_hdr
                 
                 emom0_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_emom0'+ $
                    res_str+'.fits'

                 writefits, emom0_fname $
                            , e_mom0, mom0_hdr
                 
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; MOMENT 1
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 blank_mom1 = where((mom0 le e_mom0*mom0_thresh) $
                                    or (e_mom1 gt mom1_thresh), mom1_ct)
                 if mom1_ct gt 0 then begin
                    mom1[blank_mom1] = !values.f_nan
                    e_mom1[blank_mom1] = !values.f_nan
                 endif

                 mom1_hdr = twod_head(cube_hdr)
                 sxaddpar, mom1_hdr, 'BUNIT', 'KM/S'
                 mom1[blank_ind] = !values.f_nan

                 mom1_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_mom1'+ $
                    res_str+'.fits'

                 writefits, mom1_fname $
                            , mom1, mom1_hdr
                 
                 emom1_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_emom1'+ $
                    res_str+'.fits'
                 
                 writefits, emom1_fname $
                            , e_mom1, mom1_hdr

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; DISPLAY
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 !p.multi = [0, 2, 2]
                 loadct, 33
                 disp, tpeak, /sq
                 disp, tpeak_12p5, /sq
                 disp, mom0, /sq
                 disp, mom1, /sq

              endfor

           endfor
           
        endfor

     endfor

  endif

end
