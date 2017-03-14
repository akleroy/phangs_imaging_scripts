pro build_release_v0p5 $
   , just=just $
   , reset=reset $
   , copy=copy $
   , rescale=rescale $
   , process=process $
   , pbcorr=pbcorr $
   , roundbeam=roundbeam $
   , convolve=do_convolve $
   , smooth=do_smooth $
   , noise=do_noise $
   , mask=do_masks $
   , collapse=do_collapse $
   , compile=do_compile

;+
;
; Scripts to build the imaged data into a data release.
;
;-

; DIRECTORIES
  version = '0.5'
  vstring = 'v0p5'
  root_imaging_dir = '../'
  release_dir = root_imaging_dir+'release/'+vstring+'/'

; GALAXIES
  dirs = $
     ['ic5332' $
      , 'ngc0628' $
      , 'ngc1087' $
      , 'ngc1300' $
      , 'ngc1385' $
      , 'ngc1433' $
      , 'ngc1672' $
      , 'ngc2835' $
      , 'ngc3351' $
      , 'ngc3627' $
      , 'ngc4254' $
      , 'ngc4303' $
      , 'ngc4321' $
      , 'ngc4535' $
      , 'ngc5068' $
      , 'ngc6744' $
      , 'ngc6744' $
     ]
  gals = $
     ['ic5332' $
      , 'ngc0628' $
      , 'ngc1087' $
      , 'ngc1300' $
      , 'ngc1385' $
      , 'ngc1433' $
      , 'ngc1672' $
      , 'ngc2835' $
      , 'ngc3351' $
      , 'ngc3627' $
      , 'ngc4254' $
      , 'ngc4303' $
      , 'ngc4321' $
      , 'ngc4535' $
      , 'ngc5068' $
      , 'ngc6744north' $
      , 'ngc6744south' $
     ]

  n_gals = n_elements(gals)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RESET THE DIRECTORY STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(reset) then begin

     spawn, 'rm -rf '+release_dir+'raw/'
     spawn, 'mkdir '+release_dir+'raw/'
     
     spawn, 'rm -rf '+release_dir+'process/'
     spawn, 'mkdir '+release_dir+'process/'
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(copy) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COPYING DATA', /info     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     message, 'I am copying the raw data to build a release.', /info
     message, '... deleting and remaking the raw data directory', /info

     ext_to_copy = $
        ['_co21_pb.fits' $
         , '_co21_residual.fits' $
         , '_co21.fits' $
         , '_co21_dirty.fits']
     n_ext = n_elements(ext_to_copy)

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if just ne gals[ii] then continue

        message, '... copying data for '+gals[ii], /info        
        for jj = 0, n_ext-1 do begin
           spawn, 'cp '+root_imaging_dir+dirs[ii]+'/'+ $
                  gals[ii]+ext_to_copy[jj]+' '+ $
                  release_dir+'raw/.'
        endfor
     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED COPYING DATA', /info     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRIMARY BEAM CORRECTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(pbcorr) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'PRIMARY BEAM CORRECTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     print, "!!!ALERT!!! I am applying the mosaic bug patch for 4.7.1"
     print, "!!!ALERT!!! This will be WRONG for imaging in another version"

     pb_limit = 0.5

     for ii = 0, n_gals-1 do begin

        message, "Applying primary beam correction for "+gals[ii], /info

        if n_elements(just) gt 0 then $
           if just ne gals[ii] then continue

;       Read our data

        message, '... primary beam correction for '+gals[ii], /info        
        
        clean_cube = readfits(release_dir+'raw/'+$
                              gals[ii]+'_co21.fits', clean_hdr)

        resid_cube = readfits(release_dir+'raw/'+$
                              gals[ii]+'_co21_residual.fits', resid_hdr)

        pb_cube =  readfits(release_dir+'raw/'+$
                            gals[ii]+'_co21_pb.fits', pb_hdr)
        
;       Blank everything blow the pb_limit
        
;        blank_ind = where(pb_cube eq sxpar(pb_hdr, 'BLANK'),
;        blank_ct)
        blank_ind = where(pb_cube eq pb_cube[0,0,0], blank_ct)
        if blank_ct gt 0 then begin
           clean_cube[blank_ind] = !values.f_nan
           resid_cube[blank_ind] = !values.f_nan
           pb_cube[blank_ind] = !values.f_nan
        endif

        loadct, 0
        maxpb = max(pb_cube,dim=3,/nan)
        disp, maxpb, min=0., max=1.0, /sq
        contour, finite(maxpb), lev=[1], /overplot, color=cgcolor('red')

        ind = where(pb_cube lt pb_limit)        
        clean_cube[ind] = !values.f_nan

;       Correct for the mosaic scaling bug

        model_cube = clean_cube - resid_cube
        clean_cube = model_cube*pb_cube + resid_cube

        sxdelpar, clean_hdr, 'BLANK'
        sxdelpar, clean_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_flat.fits', clean_cube, clean_hdr

        clean_cube = model_cube + resid_cube/pb_cube

        sxdelpar, clean_hdr, 'BLANK'
        sxdelpar, clean_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_pbcorr.fits', clean_cube, clean_hdr

        sxdelpar, resid_hdr, 'BLANK'
        sxdelpar, resid_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_resid.fits', resid_cube, resid_hdr

     endfor     

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ROUND BEAM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(roundbeam) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLUTION TO ROUND BEAM', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if just ne gals[ii] then continue
        
        beam_table = mrdfits(release_dir+'raw/'+gals[ii]+'_co21.fits',1)
        hdr = headfits(release_dir+'raw/'+gals[ii]+'_co21.fits')

        bmaj_max = max(beam_table.bmaj, /nan)
        pix = abs(sxpar(hdr,'CDELT1'))
        target_bmaj = sqrt(max(beam_table.bmaj,/nan)^2+(pix*2.*3600.)^2)

        ;for kk = 0, n_elements(beam_table)-1 do begin
        ;   dummy = $
        ;      calc_conv_beam(start=[beam_table[kk].bmaj, beam_table[kk].bmin, beam_table[kk].bpa] $
        ;                     , target = [1,1,0.]*target_bmaj)
        ;   print, "Plane "+str(kk)+" "+str(dummy[0])+" "+str(dummy[1])+" "+str(dummy[2])
        ;endfor

        ;stop

        ext_to_convolve = $
           ['_co21_flat' $
            , '_co21_pbcorr' $
            , '_co21_resid']
        n_ext = n_elements(ext_to_convolve)
        
        for jj = 0, n_ext-1 do begin

           cube = float(readfits(release_dir+'process/'+gals[ii]+ext_to_convolve[jj]+'.fits', hdr))
           print, "convolving "+gals[ii]+ext_to_convolve[jj]

           sz = size(cube)
           for kk = 0, sz[3]-1 do begin
              counter, kk, sz[3], 'Plane '
              
              copy_hdr = twod_head(hdr)
              plane = cube[*,*,kk]
              conv_with_gauss $
                 , data=plane $
                 , hdr=copy_hdr $
                 , start_beam= $
                 [beam_table[kk].bmaj, beam_table[kk].bmin, beam_table[kk].bpa] $
                 , target = $
                 [1,1,0.]*target_bmaj $
                 , out_data=out_plane $
                 , /quiet
              cube[*,*,kk] = out_plane
              
           endfor

           sxaddpar, hdr, 'BPA', 0.0
           sxaddpar, hdr, 'BMAJ', target_bmaj/3600.
           sxaddpar, hdr, 'BMIN', target_bmaj/3600.
           sxdelpar, hdr, 'CASAMBM'
           writefits, release_dir+'process/'+gals[ii]+ext_to_convolve[jj]+'_round.fits', cube, hdr

        endfor

     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE TO 3 TIMES BEAM RESOLUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_smooth) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLUTION TO LOWER RESOLUTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if just ne gals[ii] then continue
                
        ext_to_convolve = $
           ['_co21_flat_round' $
            , '_co21_pbcorr_round' $
            , '_co21_resid_round']
        n_ext = n_elements(ext_to_convolve)
        
        for jj = 0, n_ext-1 do begin

           dir = release_dir+'process/'
           in_cube = dir+gals[ii]+ext_to_convolve[jj]+".fits"
           out_cube = dir+gals[ii]+ext_to_convolve[jj]+"_smoothed.fits"
           print, "convolving "+in_cube

           hdr = headfits(in_cube)
           target_bmaj = 3.*sxpar(hdr,'BMAJ')*3600.
           print, "Convolving to "+str(target_bmaj)

           conv_with_gauss $
              , data=in_cube $
              , target = [1.,1.,0.]*target_bmaj $
              , out_file=out_cube

        endfor

     endfor
     
  endif
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ESTIMATE NOISE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_noise) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'ESTIMATE NOISE', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if just ne gals[ii] then continue
        
        gal = gals[ii]

        ext = $
           ['_co21_clip_round' $
            , '_co21_clip_round_smoothed' $
           ]
        n_ext = n_elements(ext)

        for jj = 0, n_ext-1 do begin

           cube = readfits('../release/v0p4/process/'+gal+ext[jj]+'.fits', cube_hdr)

           make_noise_cube $
              , cube_in = cube $
              , out_cube = rms_cube $
              , box = 21 $
              , /twod_only $
              , /show $
              , /iterate
           
           writefits, '../release/v0p4/process/'+gal+ext[jj]+'_noise.fits', rms_cube, cube_hdr

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
     
     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if just ne gals[ii] then continue
        
        gal = gals[ii]
     
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
; MAKE MASKS FOR BOTH RESOLUTIONS
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   
        ext = $
           ['_co21_clip_round' $
            , '_co21_clip_round_smoothed' $
           ]
        n_ext = n_elements(ext)
        
        for jj = 0, n_ext-1 do begin

           cube = readfits('../release/v0p4/process/'+gal+ext[jj]+'.fits', cube_hdr)
           rms_cube = readfits('../release/v0p4/process/'+gal+ext[jj]+'_noise.fits', noise_hdr)

           ppbeam = calc_pixperbeam(hdr=cube_hdr)         

           make_cprops_mask $
              , indata=cube $
              , inrms = rms_cube $
              , lo_thresh = 2 $
              , hi_thresh = 3 $
              , hi_nchan = 3 $
              , min_area = ppbeam $
              , outmask=mask
           
           !p.multi=[0,2,1]
           loadct, 33
           disp, max(cube, dim=3, /nan)
           contour, max(mask, dim=3), lev=[1], /overplot

           disp, max(cube, dim=2, /nan)
           contour, max(mask, dim=2), lev=[1], /overplot
           
           mask_hdr = cube_hdr
           sxaddpar, mask_hdr, 'BUNIT', 'MASK'
           writefits, '../release/v0p4/process/'+gal+ext[jj]+'_mask.fits', mask, mask_hdr

           if jj eq 0 then $
              join_mask = mask $
           else $
              join_mask = join_mask or mask

        endfor   

        mask_hdr = cube_hdr
        sxaddpar, mask_hdr, 'BUNIT', 'MASK'
        writefits, '../release/v0p4/process/'+gal+'_joint_mask.fits', join_mask, mask_hdr
        
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
; SMOOTH THE MASK TO MAKE A CLEAN MASK
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

        conv_with_gauss $
           , data = join_mask*1.0 $
           , hdr=mask_hdr $
           ;, start_beam = [0,0,0] $
           , target = [15.,15.,0.] $
           , out_data= smooth_mask $
           , /perbeam
        smooth_mask = smooth_mask gt 1.0

        !p.multi=[0,2,1]
        loadct, 33
        disp, max(cube, dim=3, /nan)
        contour, max(smooth_mask, dim=3), lev=[1], /overplot
        
        disp, max(cube, dim=2, /nan)
        contour, max(smooth_mask, dim=2), lev=[1], /overplot

;        print, "This is the candidate clean mask. Hit key to continue."
;        ch = get_kbrd(1)

        writefits, '../release/v0p4/process/'+gal+'_large_mask.fits', smooth_mask, mask_hdr

     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COLLAPSE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_collapse) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COLLAPSE INTO MOMENTS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if just ne gals[ii] then continue
        
        gal = gals[ii]

        pbcor = readfits('../release/v0p4/process/'+gal+'_co21_pbcorr_round.fits', hdr)
        mask = readfits('../release/v0p4/process/'+gal+'_co21_clip_round_smoothed_mask.fits', mask_hdr)

        jtok = calc_jtok(hdr=hdr)
        pbcor *= jtok
        sxaddpar, hdr, 'BUNIT', 'K'
        sxaddpar, hdr, 'JTOK', jtok
        writefits, '../release/v0p4/process/'+gal+'_co21_correct.fits', pbcor, hdr

        collapse_cube $
           , cube=pbcor $
           , hdr=hdr $
           , mask=mask $
           , mom0 = mom0 $
           , e_mom0 = e_mom0 $
           , mom1 = mom1 $
           , e_mom1 = e_mom1

        tpeak_mask = mask
        sz = size(mask)
        for kk = 0, sz[3]-1 do $
           tpeak_mask[*,*,kk] = total(mask[*,*,kk]) ge 1

        tpeak = max(pbcor*tpeak_mask, dim=3, /nan)
        tpeak_hdr = twod_head(hdr)
        sxaddpar, tpeak_hdr, 'BUNIT', 'K'

        writefits, '../release/v0p4/process/'+ $
                   gal+'_co21_tpeak.fits', tpeak, tpeak_hdr

        tpeak_12p5 = max(smooth(pbcor,[1,1,5],/nan,/edge_wrap)*tpeak_mask, dim=3, /nan)
        tpeak_hdr = twod_head(hdr)
        sxaddpar, tpeak_hdr, 'BUNIT', 'K'
        writefits, '../release/v0p4/process/'+ $
                   gal+'_co21_tpeak_12p5kms.fits', tpeak_12p5, tpeak_hdr
        
        mom0_hdr = twod_head(hdr)
        sxaddpar, mom0_hdr, 'BUNIT', 'K*KM/S'
        writefits, '../release/v0p4/process/'+ $
                   gal+'_co21_mom0.fits', mom0, mom0_hdr
        writefits, '../release/v0p4/process/'+ $
                   gal+'_co21_emom0.fits', e_mom0, mom0_hdr
        
        mom1_hdr = twod_head(hdr)
        sxaddpar, mom1_hdr, 'BUNIT', 'KM/S'
        writefits, '../release/v0p4/process/'+ $
                   gal+'_co21_mom1.fits', mom1, mom1_hdr
        writefits, '../release/v0p4/process/'+ $
                   gal+'_co21_emom1.fits', e_mom1, mom1_hdr
                
        !p.multi = [0, 2, 2]
        loadct, 33
        disp, tpeak, /sq
        disp, tpeak_12p5, /sq
        disp, mom0, /sq
        disp, mom1, /sq

     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SHUFFLE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMPILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_compile) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COMPILE A RELEASE', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
 
     dir = '../release/v0p4/'

     ext_in = $
        ['_co21_correct.fits' $
         ,'_co21_clip_round_smoothed_mask.fits' $
         ,'_co21_clip_round_mask.fits' $
         ,'_large_mask.fits' $
         ,'_co21_mom0.fits' $
         ,'_co21_emom0.fits' $
         ,'_co21_mom1.fits' $
         ,'_co21_emom1.fits' $
         ,'_co21_tpeak.fits' $
         ,'_co21_tpeak_12p5kms.fits' $
        ]

     ext_out = $
        ['_co21.fits' $
         , '_co21_smoothedmask.fits' $
         , '_co21_brightmask.fits' $
         , '_co21_widemask.fits' $
         ,'_co21_mom0.fits' $
         ,'_co21_emom0.fits' $
         ,'_co21_mom1.fits' $
         ,'_co21_emom1.fits' $
         ,'_co21_tpeak.fits' $
         ,'_co21_tpeak_12p5kms.fits' $
        ]

     n_ext = n_elements(ext_in)

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if just ne gals[ii] then continue
        
        gal = gals[ii]

        print, 'Copying '+gal

        for jj = 0, n_ext-1 do begin

           spawn, 'rm -rf '+dir+'delivery/'+gal+ext_out[jj]

           spawn, 'cp '+dir+'process/'+gal+ext_in[jj]+' '+dir+'delivery/'+gal+ext_out[jj]

        endfor

     endfor

  endif

end
