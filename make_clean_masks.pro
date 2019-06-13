pro make_clean_masks $
   , nopause=nopause $
   , inspect=do_inspect $
   , start = start_num $
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

; ... define arrays
  
  array_list = ['7m', '12m', '12m+7m']
  n_array = n_elements(array_list)

; ... define products
  
  product_list = ['co21']
  n_product = n_elements(product_list)

; ... lists of galaxies to skip or focus on

  if n_elements(only) eq 0 then $
     only = []
  if n_elements(skip) eq 0 then $
     skip = []

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; GALAXY-SPECIFIC TUNING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; ... these galaxies have central bright sources

  has_bright_center = $
     ['circinus' $
      , 'ngc0253' $
      , 'ngc1097' $
      , 'ngc1300' $
      , 'ngc1317' $
      , 'ngc1365' $
      , 'ngc1433' $
      , 'ngc1512' $
      , 'ngc1566' $
      , 'ngc1637' $
      , 'ngc1672' $
      , 'ngc2566' $
      , 'ngc2903' $
      , 'ngc2997' $
      , 'ngc3351' $
      , 'ngc3507' $
      , 'ngc3626' $    
      , 'ngc3627' $
      , 'ngc4293' $
      , 'ngc4303' $
      , 'ngc4321' $
      , 'ngc4457' $
      , 'ngc4535' $
      , 'ngc4536' $
      , 'ngc4548' $
      , 'ngc4569' $
      , 'ngc4579' $
      , 'ngc4826' $
      , 'ngc4941' $
      , 'ngc4951' $
      , 'ngc5128' $
      , 'ngc5248' $
      , 'ngc5643' $
      , 'ngc6300' $
      , 'ngc7496' $
     ]

  is_super_faint = $
     ['ngc0300' $
     ]

  empty = $
     ['ngc3239' $
     ]

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
        if total(gals[ii] eq skip) gt 0 then begin
           print, "Skipping "+gals[ii]
           continue
        endif
     endif

     if n_elements(only) gt 0 then begin
        if (total(gals[ii] eq only) eq 0) and $
           (total(this_gal eq only) eq 0) $
        then begin
           print, "Skipping "+gals[ii]
           continue
        endif
     endif
     
     print, ""
     print, "Making clean mask for "+gals[ii]+' (Galaxy '+str(ii)+')'
     print, ""

     dir = release_dir+'process/'
     print, dir, version, vstring

     clean_mask_fname = '../clean_masks/'+this_gal+'_co21_clean_mask.fits'

;if n_elements(only) gt 0 then $
;   if total(only eq gals[ii]) eq 0 then continue

     cube_fname = dir+this_gal+'_7m+tp_co21_flat_round_k.fits'
     if file_test(cube_fname) eq 0 then begin
        cube_fname = dir+this_gal+'_7m_co21_flat_round_k.fits'
     endif
     if file_test(cube_fname) eq 0 then begin
        cube_fname = dir+this_gal+'_12m_co21_flat_round_k.fits'
     endif
     if file_test(cube_fname) eq 0 then begin
        print, ""
        print, "WARNING!!!"
        print, ""
        print, "No cube found for "+this_gal
        print, "... looking for "+cube_fname
        print, "Continuing."
        continue
     endif
     cube = readfits(cube_fname, cube_hdr)

     conv_with_gauss $
        , data=cube*1.0 $
        , hdr=cube_hdr $
        , target = [1,1,0.] * 33.0 $
        , out_data=lowres_cube $
        , out_hdr=lowres_hdr $
        , /perbeam        
     lowres_cube = smooth(lowres_cube, [1,1,20], /nan)
     
     lowres_ppbeam = calc_pixperbeam(hdr=lowres_hdr)
     rms_lowres_cube = mad(lowres_cube)

     if total(this_gal[0] eq is_super_faint) ge 1 then begin
        lo_thresh_lowres = 3
        hi_thresh_lowres = 5
     endif else begin
        lo_thresh_lowres = 3
        hi_thresh_lowres = 10
     endelse

     make_cprops_mask $
        , indata=lowres_cube $
        , inrms = rms_lowres_cube $
        , lo_thresh = lo_thresh_lowres $
        , hi_thresh = hi_thresh_lowres $
        , hi_nchan = 2 $
        , min_area = 2.0*lowres_ppbeam $
        , min_pix = 4.0*lowres_ppbeam $
        , outmask=mask ;lowres_mask
     
     mask = grow_mask(mask, /z_only, iter=4)

     if this_gal[0] eq 'ngc0300' then begin
        mask[*,*,0:25] = 0B
        mask[*,*,150:*] = 0B
        mask = grow_mask(mask, /xy_only, iter=20)
     endif

     if total(this_gal[0] eq has_bright_center) ge 1 then begin

        s = gal_data(this_gal)
        make_axes, lowres_hdr, ri=ri, di=di
        rad = sphdist(ri, di, s.ra_deg, s.dec_deg, /deg)
        center_mask = rad le 20./3600.
        sz = size(lowres_cube)
        for kk = 0, sz[3]-1 do begin
           mask[*,*,kk] = mask[*,*,kk] or center_mask
        endfor

     endif

     ;conv_with_gauss $
     ;   , data=lowres_mask*1.0 $
     ;   , hdr=lowres_hdr $
     ;   , target = [1,1,0.] * 33.0 $
     ;   , out_data=out_lowres_mask $
     ;   , /perbeam        
     ;if total(gals[ii] eq expand_in_space) gt 0 then begin
     ;   lowres_mask = out_lowres_mask ge 0.5
     ;endif else begin
     ;   lowres_mask = out_lowres_mask ge 1.0
     ;endelse

     writefits, clean_mask_fname, mask, lowres_hdr

     inspect_clean_masks $
        , nopause=nopause $
        , only=this_gal $
        , array='7m'

     ;inspect_clean_masks $
     ;   , nopause=nopause $
     ;   , only=this_gal $
     ;   , array='12m+7m'

  endfor


end
