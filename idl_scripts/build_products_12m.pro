pro build_products_12m $
   , version=version $
   , only=only $
   , skip=skip $
   , start_with=start_with $
   , stop_at=stop_at $
   , just_array=just_array $
   , clear_convolve=do_clear_convolve $
   , convolve=do_conv_to_res $
   , prelim=do_prelim $
   , noise=do_noise $
   , mask=do_masks $
   , hybrid=do_hybrid $
   , collapse=do_collapse $
   , highlevel=do_highlevel $
   , target_res=target_res

;+
;
; Scripts to build the data cubes into maps and products for release.
;
;-

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

; SAMPLE TABLE
  sample_version_string = '1.2'
  sample_tab_file = 'phangs_sample_table_v1p2.fits'  
  sample_tab = mrdfits(sample_tab_file, 1, h)

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

; ARRAYS
  
  array_list = ['12m', '12m+7m']
  n_array = n_elements(array_list)

  fullarray_list = ['12m', '12m+tp', '12m+7m', '12m+7m+tp']
  n_fullarray = n_elements(fullarray_list)

; PRODUCTS
  
  product_list = ['co21']
  n_product = n_elements(product_list)

; RESOLUTIONS

  if n_elements(target_res) eq 0 then begin
     target_res = [-1, 60, 90, 120, 150, 500]
; target_res = [-1, 45, 60, 80, 100, 120, 500, 750, 1000, 1250, 1500, 2000]
  endif
  n_res = n_elements(target_res)

; PRELIMINARY RESOLUTION FOR MASKING
  prelim_res = 500

; LOW RESOLUTION STRING FOR USE IN BROAD MASKS
  lowres_string = '500pc'

; TUNING PARAMETERS FOR THE VELOCITY FIELD

  vfield_reject_thresh = 30.
  mom0_thresh_for_mom1 = 2.0

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; IF REQUESTED, WIPE PREVIOUS VERSIONS OF THE CONVOLUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_clear_convolve) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'REMOVING PREVIOUS VERSIONS OF THE CONVOLUTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     

     dir = release_dir+'process/'
     tol = 0.1

     ;s = gal_data(gal_for_fullgals)

     ext_to_process = $
        ['flat_round_k' $
         , 'pbcorr_round_k']
     n_ext = n_elements(ext_to_process)
     
     first = 0B
     last = 0B

     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        if n_elements(start_with) gt 0 then begin
           if this_gal eq start_with then first = 1B
           if first eq 0 then continue
        endif

        if n_elements(stop_at) gt 0 then begin
           if this_gal eq stop_at then last = 1B
           if last eq 1B then continue
        endif

        message, "Wiping convolved cubes for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue
           
           flist = file_search(release_dir+'process/'+ $
                               this_gal+'_'+this_array+'_*pc.fits' $
                               , count=fct)
           
           if fct eq 0 then continue

           for kk = 0, fct-1 do begin

              command = 'rm -rf '+flist[kk]
              print, command
              spawn, command

           endfor

        endfor

     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'DONE CLEARING PREVIOUS CONVOLUTIONS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE TO SPECIFIC RESOLUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_conv_to_res) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLVING TO FIXED PHYSICAL RESOLUTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     
     
     dir = release_dir+'process/'
     tol = 0.1

     s = gal_data(gal_for_fullgals)

     ext_to_process = $
        ['flat_round_k' $
         , 'pbcorr_round_k']
     n_ext = n_elements(ext_to_process)

     first = 0B
     last = 0B
     
     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        if n_elements(start_with) gt 0 then begin
           if this_gal eq start_with then first = 1B
           if first eq 0 then continue
        endif

        if n_elements(stop_at) gt 0 then begin
           if this_gal eq stop_at then last = 1B
           if last eq 1B then continue
        endif

        sample_ind = where(s[ii].pgc eq sample_tab.pgc, sample_ct)
        if sample_ct eq 0 then begin
           print, "Galaxy not found in sample file. Stopping."
        endif
        adopted_distance = sample_tab[sample_ind].dist

        message, "Convolving the cubes for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_process[ll]
                 
                 in_file = release_dir+'process/'+ $
                           this_gal+'_'+this_array+'_'+ $
                           this_product+'_'+this_ext+'.fits'
                 
                 test = file_search(in_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+in_file+' not found.', /info
                    continue
                 endif                 

                 cube = readfits(in_file, hdr)

                 sxaddpar, hdr, 'SAMPVER', sample_version_string
                 sxaddpar, hdr, 'SAMPTAB', sample_tab_file
                 sxaddpar, hdr, 'DIST', adopted_distance, 'MPC / USED IN CONVOLUTION'

                 current_res_pc = adopted_distance*!dtor*sxpar(hdr, 'BMAJ')*1d6
                 
                 for zz = 0, n_res -1 do begin
                    
                    res_str = strcompress(str(target_res[zz]),/rem)+'pc'
                    out_file = release_dir+'process/'+ $
                               this_gal+'_'+this_array+'_'+ $
                               this_product+'_'+this_ext+'_'+res_str+'.fits'
                    target_res_as = target_res[zz]/(adopted_distance*1d6)/!dtor*3600.d
                    
                    if current_res_pc gt (1.0+tol)*target_res[zz] then begin
                       print, strupcase(this_gal)+": Resolution too coarse. Skipping."
                       continue
                    endif
                    
                    if abs(current_res_pc - target_res[zz])/target_res[zz] lt tol then begin
                       print, strupcase(this_gal)+": I will call ", current_res_pc, " ", target_res[zz]
                       writefits, out_file, cube, hdr
                    endif else begin
                       print, strupcase(this_gal)+": I will convolve ", current_res_pc $
                              , " to ", target_res[zz]
                       conv_with_gauss $
                          , data=cube $
                          , hdr=hdr $
                          , target_beam=target_res_as*[1,1,0] $
                          , out_file=out_file
                    endelse
                    
                 endfor
                 
              endfor
              
           endfor
           
        endfor

     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLVING TO FIXED PHYSICAL RESOLUTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GENERATE A LOW RESOLUTION MASK FROM THE FLAT CUBE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_prelim) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'PRELIMINARY MASK', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     ext_to_process = $
        ['flat_round_k']
     n_ext = n_elements(ext_to_process)
     
     first = 0B
     last = 0B
     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        if n_elements(start_with) gt 0 then begin
           if this_gal eq start_with then first = 1B
           if first eq 0 then continue
        endif

        if n_elements(stop_at) gt 0 then begin
           if this_gal eq stop_at then last = 1B
           if last eq 1B then continue
        endif

        message, "A low res initial mask for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_process[ll]

                 res_str = '_'+strcompress(str(prelim_res),/rem)+'pc'
                 print, "Resolution "+res_str

                 in_file = release_dir+'process/'+ $
                           this_gal+'_'+this_array+'_'+ $
                           this_product+'_'+this_ext+res_str+'.fits'
                 
                 test = file_search(in_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+in_file+' not found.', /info
                    continue
                 endif                 
                 
                 cube = readfits(in_file, cube_hdr)
                 
                 make_noise_cube $
                    , cube_in = cube $
                    , out_cube = rms $
                    , /zero_only $
                    , /show $
                    , /iterate
                 
                 ppbeam = calc_pixperbeam(hdr=cube_hdr)         
                 
                 make_cprops_mask $
                    , indata=cube $
                    , inrms = rms $
                    , lo_thresh = 2 $
                    , hi_thresh = 3.5 $
                    , hi_nchan = 3 $
                    , min_area = ppbeam $
                    , outmask=prelim_mask

                 out_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    'prelim_mask.fits'
                 writefits, out_file, prelim_mask, cube_hdr
                 
              endfor
              
           endfor

        endfor

     endfor

  endif

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

     first = 0B
     last = 0B

     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        if n_elements(start_with) gt 0 then begin
           if this_gal eq start_with then first = 1B
           if first eq 0 then continue
        endif

        if n_elements(stop_at) gt 0 then begin
           if this_gal eq stop_at then last = 1B
           if last eq 1B then continue
        endif

        message, "Estimating noise for for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           prelim_mask_file = $
              release_dir+'process/'+ $
              this_gal+'_'+this_array+'_'+ $
              'prelim_mask.fits'
           
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

                    if file_test(prelim_mask_file) then begin
                       prelim_mask = readfits(prelim_mask_file)
                    endif else begin
                       prelim_mask = finite(cube) eq 0 
                    endelse

                    make_noise_cube $
                       , cube_in = cube $
                       , out_cube = rms_cube $
                       , mask_in = prelim_mask $
                       , box=11, spec_box=5 $
                       , /show $
                       , /iterate

                                ;, twod_only, box=5 $
                                ;, /oned_only, spec_box=2 $
                    
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
; BUILD MASKS HOLDING BRIGHT SIGNAL AT EACH RESOLUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_masks) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'BUILD MASKS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     first = 0B
     last = 0B
     
     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        if n_elements(start_with) gt 0 then begin
           if this_gal eq start_with then first = 1B
           if first eq 0 then continue
        endif

        if n_elements(stop_at) gt 0 then begin
           if this_gal eq stop_at then last = 1B
           if last eq 1B then continue
        endif

        message, "Building simple masks for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

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
                 
                 test = file_search(cube_file, count=found)
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
                    , lo_nchan = 2 $
                    , hi_thresh = 3.5 $
                    , hi_nchan = 3 $
                    , min_area = ppbeam $
                    , outmask=mask                    

                 make_cprops_mask $
                    , indata=-1.*cube $
                    , inrms = rms_cube $
                    , lo_thresh = 2 $
                    , lo_nchan = 2 $
                    , hi_thresh = 3.5 $
                    , hi_nchan = 3 $
                    , min_area = ppbeam $
                    , outmask=negmask

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

                 out_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_negmask'+ $
                    res_str+'.fits'

                 writefits, out_file $
                            , negmask, mask_hdr
                 
              endfor

           endfor
           
        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; HYBRIDIZE THE MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_hybrid) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'HYBRIDIZE MASKS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     first = 0B
     last = 0B

     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        if n_elements(start_with) gt 0 then begin
           if this_gal eq start_with then first = 1B
           if first eq 0 then continue
        endif

        if n_elements(stop_at) gt 0 then begin
           if this_gal eq stop_at then last = 1B
           if last eq 1B then continue
        endif

        message, "Creating hybrid (cross scale) masks for for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              lores_fname = $
                 release_dir+'process/'+ $
                 this_gal+'_'+this_array+'_'+ $
                 this_product+ $
                 '_signalmask_'+lowres_string+'.fits'
              
              test = file_search(lores_fname, count=found)
              if found eq 0 then begin
                 message, 'File '+lores_fname+' not found.', /info
                 continue
              endif                 
              mask_lores = readfits(lores_fname)

              for zz = 0, n_res -1 do begin
                 
                 if target_res[zz] eq -1 then begin
                    res_str = ''
                    print, "Native resolution."
                 endif else begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                    print, "Resolution "+res_str
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
; COLLAPSE INTO A SIMPLE SET OF MOMENT MAPS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_collapse) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COLLAPSE INTO MOMENTS USING SIMPLE TECHNIQUES', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     dir = release_dir+'process/'

     first = 0B
     last = 0B

     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        if n_elements(start_with) gt 0 then begin
           if this_gal eq start_with then first = 1B
           if first eq 0 then continue
        endif

        if n_elements(stop_at) gt 0 then begin
           if this_gal eq stop_at then last = 1B
           if last eq 1B then continue
        endif

        message, "Making moments for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              for zz = 0, n_res -1 do begin
                 
                 if target_res[zz] eq -1 then begin
                    res_str = ''
                    print, "Native resolution."
                 endif else begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                    print, "Resolution "+res_str
                 endelse
                 
; Read the noise. Note that this is constructed in the "flat" cube to
; make masking easier, so we need to correct for that. We will do this
; inelegantly by dividing the flat and pbcorr cube at the appropriate
; resolution to construct a correction factor.

; REPLACE THIS WITH THE PRIMARY BEAM CORRECTED VERSION WHEN AVAILABLE

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
                 
; Read the cube itself.

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

                 cube_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_pbcorr_round_k'+ $
                    res_str+'.fits'

; Now read the mask appropriate for this resolution. For this step, we
; use only the signalmask and construct simple moments.
                 
                 mask_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_signalmask'+ $
                    res_str+'.fits'

                 test = file_search(mask_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+mask_file+' not found.', /info
                    continue
                 endif                 
                 mask = readfits(mask_file, mask_hdr)
                 
; Blank regions of the moment map with no coverage. This is different
; than zeros, which indicate a region covered but not in the mask.

                 blank_ind = where(total(finite(cube),3) eq 0)

; Calculate moments.

; For the moment 0, moment 1, and the equivalent width, grow the mask
; by two channels (5 km/s) along the velocity axis.

                 mom0_mask = $
                    grow_mask(mask, iters=2, /z_only)
                 
                 collapse_cube $
                    , cube=cube $
                    , hdr=cube_hdr $
                    , mask=mom0_mask $
                    , noise=rms_cube $
                    , mom0 = mom0 $
                    , e_mom0 = e_mom0 $
                    , mom1 = mom1 $
                    , e_mom1 = e_mom1 $
                    , ew = ew $
                    , e_ew = e_ew $
                    , tpeak = tpeak                                  

                 mom0_hdr = twod_head(cube_hdr)
                 sxaddpar, mom0_hdr, 'BUNIT', 'K*KM/S'
                 mom0[blank_ind] = !values.f_nan
                 
                 mom0_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_mom0'+ $
                    res_str+'.fits'

                 writefits, mom0_fname, mom0, mom0_hdr
                 
                 emom0_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_emom0'+ $
                    res_str+'.fits'

                 writefits, emom0_fname, e_mom0, mom0_hdr

                 ew_hdr = twod_head(cube_hdr)
                 sxaddpar, ew_hdr, 'BUNIT', 'KM/S'
                 ew[blank_ind] = !values.f_nan

                 ew_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_ew'+ $
                    res_str+'.fits'

                 writefits, ew_fname, ew, ew_hdr
                 
                 eew_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_eew'+ $
                    res_str+'.fits'
                 
                 writefits, eew_fname, e_ew, ew_hdr

; For the moment 2 map use the signal-only version of the cube (i.e.,
; no extension along the velocity axis). Also produce peak-to-edge
; measurements that can be used to correct the line width.

                 collapse_cube $
                    , cube=cube $
                    , hdr=cube_hdr $
                    , mask=mask $
                    , noise=rms_cube $
                    , mom1 = mom1 $
                    , e_mom1 = e_mom1 $
                    , mom2 = mom2 $
                    , e_mom2 = e_mom2 $
                    , tpeak = tpeak $                  
                    , tmin = tmin $
                    , vquad = vquad $
                    , e_vquad = e_vquad $
                    , vpeak = vpeak

                 peak_to_edge = tpeak/tmin

                 mom2_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_mom2'+ $
                    res_str+'.fits'
                 
                 mom2_hdr = twod_head(cube_hdr)
                 sxaddpar, mom2_hdr, 'BUNIT', 'KM/S'
                 mom2[blank_ind] = !values.f_nan

                 writefits, mom2_fname, mom2, mom2_hdr

                 p2e_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_p2e'+ $
                    res_str+'.fits'
                 
                 p2e_hdr = twod_head(cube_hdr)
                 sxaddpar, p2e_hdr, 'BUNIT', 'NONE'
                 peak_to_edge[blank_ind] = !values.f_nan

                 writefits, p2e_fname, peak_to_edge, p2e_hdr

                 mom1_hdr = twod_head(cube_hdr)
                 sxaddpar, mom1_hdr, 'BUNIT', 'KM/S'
                 mom1[blank_ind] = !values.f_nan

                 mom1_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_mom1'+ $
                    res_str+'.fits'

                 writefits, mom1_fname, mom1, mom1_hdr
                 
                 emom1_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_emom1'+ $
                    res_str+'.fits'
                 
                 writefits, emom1_fname, e_mom1, mom1_hdr

                 vpeak_hdr = twod_head(cube_hdr)
                 sxaddpar, vpeak_hdr, 'BUNIT', 'KM/S'
                 vpeak[blank_ind] = !values.f_nan

                 vpeak_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_vpeak'+ $
                    res_str+'.fits'
                 
                 writefits, vpeak_fname, vpeak, vpeak_hdr                 

                 vquad_hdr = twod_head(cube_hdr)
                 sxaddpar, vquad_hdr, 'BUNIT', 'KM/S'
                 vquad[blank_ind] = !values.f_nan

                 vquad_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_vquad'+ $
                    res_str+'.fits'
                 
                 writefits, vquad_fname, vquad, vquad_hdr                 

                 evquad_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_evquad'+ $
                    res_str+'.fits'
                 
                 writefits, evquad_fname, e_vquad, vquad_hdr

; The mask used for the peak temperature is right now the full width
; of the mask across the galaxy. This can be improved at least a bit
; while still staying general. One improvement would be to use the
; maximum window size anywhere in the galaxy. A better map would come
; from the very low resolution masks, but this is less totally
; agnostic.

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
                    this_product+'_strict_tpeak'+ $
                    res_str+'.fits'

                 writefits, tpeak_fname, tpeak, tpeak_hdr
                 
                 tpeak_12p5 = max(smooth(cube,[1,1,5],/nan,/edge_wrap)*tpeak_mask, dim=3, /nan)
                 tpeak_12p5[blank_ind] = !values.f_nan
                 tpeak_hdr = twod_head(cube_hdr)
                 sxaddpar, tpeak_hdr, 'BUNIT', 'K'

                 tpeak12p5_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_tpeak12p5kms'+ $
                    res_str+'.fits'

                 writefits, tpeak12p5_fname, tpeak_12p5, tpeak_hdr
                                  
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; DISPLAY
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 !p.multi = [0, 3, 2]
                 loadct, 33
                 disp, tpeak_12p5, /sq
                 disp, mom0, /sq
                 disp, ew, /sq
                 disp, mom1, /sq
                 disp, vpeak, /sq
                 disp, vquad, /sq
                 ;disp, mom2, /sq

              endfor

           endfor
           
        endfor

     endfor

  endif  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE MAPS USING MORE SOPHISTICATED MASKING TECHNIQUES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_highlevel) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COLLAPSE INTO FANCIER MOMENTS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     dir = release_dir+'process/'

     first = 0B
     last = 0B

     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        if n_elements(start_with) gt 0 then begin
           if this_gal eq start_with then first = 1B
           if first eq 0 then continue
        endif

        if n_elements(stop_at) gt 0 then begin
           if this_gal eq stop_at then last = 1B
           if last eq 1B then continue
        endif

        message, "Making high level products for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

; Loop over resolutions

              for zz = 0, n_res -1 do begin
                 
                 if target_res[zz] eq -1 then begin
                    res_str = ''
                    print, "Native resolution."
                 endif else begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                    print, "Resolution "+res_str
                 endelse
                 
; Read the noise. Note that this is constructed in the "flat" cube to
; make masking easier, so we need to correct for that. We will do this
; inelegantly by dividing the flat and pbcorr cube at the appropriate
; resolution to construct a correction factor.

; REPLACE THIS WIT THE PRIMARY BEAM CORRECTED VERSION WHEN AVAILABLE
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
                 
; Read the cube itself.

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

                 cube_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_pbcorr_round_k'+ $
                    res_str+'.fits'

; Now read the mask appropriate for this resolution. We want to use
; the "hybrid" masks, which combine the low low resolution detections
; with the signal mask at this resolution. If a hybrid mask is not
; present, just skip this combination. The strict moments will still
; exist.
                 
                 mask_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_hybridmask'+ $
                    res_str+'.fits'

                 if file_test(mask_file) eq 0 then begin
                    message, 'File '+mask_file+' not found. Skipping this step.', /info
                    continue
                 endif                 
                 mask = readfits(mask_file, mask_hdr)
                              
; Blank regions of the moment map with no coverage. This is different
; than zeros, which indicate a region covered but not in the mask.

                 blank_ind = where(total(finite(cube),3) eq 0)
    
; Calculate moments.

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

; Moment 0
                 
                 mom0_hdr = twod_head(cube_hdr)
                 sxaddpar, mom0_hdr, 'BUNIT', 'K*KM/S'
                 mom0[blank_ind] = !values.f_nan
                 
                 mom0_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_broad_mom0'+ $
                    res_str+'.fits'
                 
                 writefits, mom0_fname, mom0, mom0_hdr
                 
                 emom0_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_broad_emom0'+ $
                    res_str+'.fits'

                 writefits, emom0_fname, e_mom0, mom0_hdr
                 
; Moment 1

; The recipe here gets a little complicated. We always take the
; "strict" moment1 if we have one. But then we also calculate a
; moment1 from the broader mask. We only keep this, however, if it
; agrees with a prior on the velocity field within some tolerance. For
; now, this is the low resolution velocity field. We also require some
; nominal signal to noise in the moment0 to believe that we have a
; reasonable moment1.
;

; Read the "strict" velocity field and the low resolution prior.

                 vfield_strict_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_mom1'+ $
                    res_str+'.fits'                    

                 e_vfield_strict_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_emom1'+ $
                    res_str+'.fits'                    

                 if file_test(vfield_strict_fname) eq 0 then begin
                    message, 'No strict velocity field information found ... '+ $
                             vfield_strict_fname+' ... skipping.', /info
                    continue
                 endif

                 vfield_strict = readfits(vfield_strict_fname, prior_hdr)
                 e_vfield_strict = readfits(e_vfield_strict_fname, prior_hdr)

                 vfield_prior_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_strict_mom1'+ $
                    '_'+lowres_string+'.fits'

                 if file_test(vfield_prior_fname) eq 0 then begin
                    message, 'No prior velocity field information found ... '+ $
                             vfield_prior_fname+' ... skipping.', /info
                    continue
                 endif

                 vfield_prior = readfits(vfield_prior_fname, prior_hdr)

; Now build the combined moment1 map

                 has_strict_vfield = $
                    finite(vfield_strict)
                 
                 mom1_reject_ind = $
                    where((has_strict_vfield eq 0) and $
                          ((abs(mom1 - vfield_prior) gt vfield_reject_thresh) or $
                           (mom0 le e_mom0*mom0_thresh_for_mom1))$
                          , mom1_reject_ct)
                 if mom1_reject_ct gt 0 then begin
                    mom1[mom1_reject_ind] = !values.f_nan
                    e_mom1[mom1_reject_ind] = !values.f_nan
                 endif

                 strict_ind = where(has_strict_vfield, strict_ct)
                 if strict_ct gt 0 then begin
                    mom1[strict_ind] = vfield_strict[strict_ind]
                    e_mom1[strict_ind] = e_vfield_strict[strict_ind]
                 endif

                 mom1_hdr = twod_head(cube_hdr)
                 sxaddpar, mom1_hdr, 'BUNIT', 'KM/S'

                 mom1_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_broad_mom1'+ $
                    res_str+'.fits'

                 writefits, mom1_fname $
                            , mom1, mom1_hdr
                 
                 emom1_fname = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+ $
                    this_product+'_broad_emom1'+ $
                    res_str+'.fits'
                 
                 writefits, emom1_fname $
                            , e_mom1, mom1_hdr

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; DISPLAY
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 !p.multi = [0, 2, 1]
                 loadct, 33
                 disp, mom0, /sq
                 disp, mom1, /sq
                 !p.multi = [0, 1, 1]

              endfor

           endfor
           
        endfor

     endfor

  endif

end
