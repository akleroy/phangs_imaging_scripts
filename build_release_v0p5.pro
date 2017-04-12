pro build_release_v0p5 $
   , just=just $
   , reset=do_reset $
   , copy=do_copy $
   , pbcorr=do_pbcorr $
   , roundbeam=do_roundbeam $
   , stage=do_stage_feather $
   , feather=do_copy_feather $
   , convolve=do_convolve $
   , smooth=do_smooth $
   , noise=do_noise $
   , mask=do_masks $
   , collapse=do_collapse $
   , shuffle=do_shuffle $
   , convtores=do_conv_to_res $
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
      , 'ngc1512' $
      , 'ngc1566' $
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
      , 'ngc1512' $
      , 'ngc1566' $
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
  array = $
     ['7M' $
      , '12M+7M' $
      , '7M' $
      , '7M' $
      , '7M' $
      , '7M' $
      , '7M' $
      , '7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
      , '12M+7M' $
     ]

  n_gals = n_elements(gals)

  mom1_thresh = 2.5d
  mom0_thresh = 3.0d

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RESET THE DIRECTORY STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_reset) then begin

     spawn, 'rm -rf '+release_dir+'raw/'
     spawn, 'mkdir '+release_dir+'raw/'
     
     spawn, 'rm -rf '+release_dir+'process/'
     spawn, 'mkdir '+release_dir+'process/'

     spawn, 'rm -rf '+release_dir+'feather/'
     spawn, 'mkdir '+release_dir+'feather/'

     spawn, 'rm -rf '+release_dir+'delivery/'
     spawn, 'mkdir '+release_dir+'delivery/'
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_copy) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COPYING DATA', /info     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     message, 'I am copying the raw data to build a release.', /info

     ext_to_copy = $
        ['_co21_pb.fits' $
         , '_co21_residual.fits' $
         , '_co21.fits' $
         , '_co21_dirty.fits']
     n_ext = n_elements(ext_to_copy)

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

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

  if keyword_set(do_pbcorr) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'PRIMARY BEAM CORRECTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     print, "!!!ALERT!!! I am applying the mosaic bug patch for 4.7.1"
     print, "!!!ALERT!!! This will be WRONG for imaging in another version"

     pb_limit = 0.5

     for ii = 0, n_gals-1 do begin

        message, "Applying primary beam correction for "+gals[ii], /info

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

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
        sxaddpar, clean_hdr, 'ARRAY', array[ii]
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_flat.fits', clean_cube, clean_hdr

        clean_cube = model_cube + resid_cube/pb_cube

        sxdelpar, clean_hdr, 'BLANK'
        sxdelpar, clean_hdr, 'CASAMBM'
        sxaddpar, clean_hdr, 'ARRAY', array[ii]
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_pbcorr.fits', clean_cube, clean_hdr

        sxdelpar, resid_hdr, 'BLANK'
        sxdelpar, resid_hdr, 'CASAMBM'
        sxaddpar, clean_hdr, 'ARRAY', array[ii]
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_resid.fits', resid_cube, resid_hdr

     endfor     

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ROUND BEAM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_roundbeam) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLUTION TO ROUND BEAM', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue
        
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
                 , /perbeam $
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
; COPY TO STAGE THE FEATHERING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_stage_feather) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'STAGING FEATHER CALL', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     spawn, 'rm -rf '+release_dir+'feather/'
     spawn, 'mkdir '+release_dir+'feather/'     

     in_dir = release_dir+'process/'
     out_dir = release_dir+'feather/'

     readcol $
        , 'singledish_key.txt' $
        , format='A,A', comment='#' $
        , sd_gal, sd_fname

     for ii = 0, n_gals-1 do begin

        message, 'Copying data to feather '+gals[ii], /info

        spawn, 'cp '+in_dir+gals[ii]+'_co21_pbcorr_round.fits '+$
               out_dir+'.'        

        sd_ind = where(sd_gal eq gals[ii], sd_ct)
        if sd_ct eq 0 then begin
           message, 'I did not find a single dish key entry for '+gals[ii], /info
        endif else begin
           spawn, 'cp '+sd_fname[sd_ind]+' '+out_dir+gals[ii]+'_tp.fits'
        endelse

     endfor

     spawn, 'cp feather_script.py '+out_dir+'/.'

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FEATHERED AND SINGLE DISH DATA INTO THE DIRECTORY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_copy_feather) then begin

     readcol $
        , 'singledish_key.txt' $
        , format='A,A', comment='#' $
        , sd_gal, sd_fname

     in_dir = release_dir+'feather/'
     out_dir = release_dir+'process/'

     for ii = 0, n_gals-1 do begin

        message, 'Copying feathered and single dish data for '+gals[ii], /info

        spawn, 'rm -rf '+$
               out_dir+gals[ii]+'_co21_feather_pbcorr.fits'
        spawn, 'cp '+in_dir+gals[ii]+'_co21_feathered.fits '+$
               out_dir+gals[ii]+'_co21_feather_pbcorr.fits'

        sd_ind = where(sd_gal eq gals[ii], sd_ct)
        if sd_ct eq 0 then begin
           message, 'I did not find a single dish key entry for '+gals[ii], /info
        endif else begin
           spawn, 'rm -rf '+$
                  out_dir+gals[ii]+'_tp.fits'
           spawn, 'cp '+sd_fname[sd_ind]+' '+out_dir+gals[ii]+'_tp.fits'
        endelse

        if sd_ct gt 0 then begin

           message, 'Cleaning feathered data for '+gals[ii], /info

           cube = readfits(out_dir+gals[ii]+'_co21_feather_pbcorr.fits', hdr)

           template = readfits(out_dir+gals[ii]+'_co21_flat_round.fits', temp_hdr)
           blank_ind = where(finite(template) eq 0 or $
                             abs(cube - sxpar(hdr,'BLANK')) lt 1d-6, blank_ct)
           if blank_ct gt 0 then $
              cube[blank_ind] = !values.f_nan
           sxaddpar, hdr, 'ARRAY', '12M+7M+TP'

           writefits, out_dir+gals[ii]+'_co21_feather_pbcorr.fits', cube, hdr

           cube = readfits(out_dir+gals[ii]+'_tp.fits', hdr)

           sxaddpar, hdr, 'ARRAY', 'TP'
           blank_ind = where(abs(cube - sxpar(hdr,'BLANK')) lt 1d-6, blank_ct)
           if blank_ct gt 0 then $
              cube[blank_ind] = !values.f_nan
           
           writefits, out_dir+gals[ii]+'_tp.fits', cube, hdr

           jtok = calc_jtok(hdr=hdr)
           sxaddpar, hdr, 'JTOK', jtok
           cube *= jtok
           sxaddpar, hdr, 'BUNIT', 'K'

           writefits, out_dir+gals[ii]+'_tp_k.fits', cube, hdr

        endif

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
           if total(just eq gals[ii]) eq 0 then continue
        
        ext_to_convolve = $
           ['_co21_flat_round' $
            , '_co21_pbcorr_round' $
            , '_co21_resid_round' $
            , '_co21_feather_pbcorr']
        n_ext = n_elements(ext_to_convolve)
        
        for jj = 0, n_ext-1 do begin

           dir = release_dir+'process/'
           in_cube = dir+gals[ii]+ext_to_convolve[jj]+".fits"
           out_cube = dir+gals[ii]+ext_to_convolve[jj]+"_smoothed.fits"
           print, "convolving "+in_cube

           fname = file_search(in_cube, count=file_ct)
           if file_ct eq 0 then begin
              message, 'No file '+in_cube, /info
              continue
           endif

           hdr = headfits(in_cube)
           target_bmaj = 3.*sxpar(hdr,'BMAJ')*3600.
           print, "Convolving to "+str(target_bmaj)

           conv_with_gauss $
              , data=in_cube $
              , target = [1.,1.,0.]*target_bmaj $
              , out_file=out_cube $
              , /perbeam

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

     dir = release_dir+'process/'

     for ii = 0, n_gals-1 do begin

        message, 'Estimating noise for '+gals[ii], /info

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue
        
        gal = gals[ii]

        ext = $
           ['_co21_flat_round' $
            , '_co21_flat_round_smoothed' $
            , '_co21_feather_pbcorr' $
            , '_co21_feather_pbcorr_smoothed' $
           ]
        n_ext = n_elements(ext)

        for jj = 0, n_ext-1 do begin

           fname = file_search(dir+gal+ext[jj]+'.fits', count=fct)

           if fct eq 0 then begin
              message, "No file found for "+fname, /info
              continue
           endif

           cube = readfits(fname, cube_hdr)
           
           make_noise_cube $
              , cube_in = cube $
              , out_cube = rms_cube $
              , box = 21 $
              , /twod_only $
              , /show $
              , /iterate
           
           writefits, dir+gal+ext[jj]+'_noise.fits', rms_cube, cube_hdr

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

     dir = release_dir+'process/'
     
     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue
        
        gal = gals[ii]
        
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
; MAKE MASKS FOR BOTH RESOLUTIONS
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

        ext = $
           ['_co21_flat_round' $
            , '_co21_flat_round_smoothed' $
            , '_co21_feather_pbcorr' $
            , '_co21_feather_pbcorr_smoothed' $
           ]
        n_ext = n_elements(ext)
        
        for jj = 0, n_ext-1 do begin
           
           fname = file_search(dir+gal+ext[jj]+'.fits', count=fct)

           if fct eq 0 then begin
              message, "No file found for "+fname, /info
              continue
           endif

           cube = readfits(dir+gal+ext[jj]+'.fits', cube_hdr)
           rms_cube = readfits(dir+gal+ext[jj]+'_noise.fits', noise_hdr)

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
           writefits, dir+gal+ext[jj]+'_mask.fits', mask, mask_hdr

        endfor   

; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
; COMBINE THE MASKS
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
        
        first = 1B
        for jj = 0, n_ext-1 do begin
           
           fname = file_search(dir+gal+ext[jj]+'.fits', count=fct)
           
           if fct eq 0 then begin
              message, "No file found for "+fname, /info
              continue
           endif
           
           this_mask = readfits(dir+gal+ext[jj]+'_mask.fits', mask_hdr)
           
           if first then begin
              joint_mask = this_mask
              first = 0B
           endif else begin
              joint_mask = (joint_mask + this_mask) ge 1
              writefits, dir+gal+'_joint_mask.fits', joint_mask, mask_hdr 
           endelse
        endfor
        
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
; SMOOTH THE JOINT MASK TO MAKE A CLEAN MASK
; -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

        target_hdr = headfits(dir+gal + '_co21_flat_round_smoothed.fits')
        target_bmaj = sxpar(target_hdr, 'BMAJ')*3600.*3.

        conv_with_gauss $
           , data = joint_mask*1.0 $
           , hdr=mask_hdr $
           , target = [1.,1.,0.]*target_bmaj $
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
        
        writefits, dir+gal+'_large_mask.fits', smooth_mask, mask_hdr
        
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
     
     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue
        
        gal = gals[ii]

        fname = file_search(dir+gal+'_co21_feather_pbcorr.fits', count=fct)        
        if fct eq 0 or $
           gal eq 'ngc1087' or $
           gal eq 'ngc4254' or $
           gal eq 'ngc5068' $
        then begin
           pbcor = readfits(dir+gal+'_co21_pbcorr_round.fits', hdr)
        endif else begin
           pbcor = readfits(fname, hdr)
        endelse

        mask = readfits(dir+gal+'_joint_mask.fits', mask_hdr)

        jtok = calc_jtok(hdr=hdr)
        pbcor *= jtok
        sxaddpar, hdr, 'BUNIT', 'K'
        sxaddpar, hdr, 'JTOK', jtok
        writefits, dir+gal+'_co21_correct.fits', pbcor, hdr

        collapse_cube $
           , cube=pbcor $
           , hdr=hdr $
           , mask=mask $
           , mom0 = mom0 $
           , e_mom0 = e_mom0 $
           , mom1 = mom1 $
           , e_mom1 = e_mom1

        blank_mom1 = where((mom0 le e_mom0*mom0_thresh) $
                           or (e_mom1 gt mom1_thresh), mom1_ct)
        if mom1_ct gt 0 then begin
           mom1[blank_mom1] = !values.f_nan
           e_mom1[blank_mom1] = !values.f_nan
        endif

        tpeak_mask = mask
        sz = size(mask)
        for kk = 0, sz[3]-1 do $
           tpeak_mask[*,*,kk] = total(mask[*,*,kk]) ge 1

        fin_map = total(finite(pbcor)) gt 0
        blank_ind = where(fin_map eq 0)

        tpeak = max(pbcor*tpeak_mask, dim=3, /nan)
        tpeak[blank_ind] = !values.f_nan
        tpeak_hdr = twod_head(hdr)
        sxaddpar, tpeak_hdr, 'BUNIT', 'K'

        writefits, dir+ $
                   gal+'_co21_tpeak.fits', tpeak, tpeak_hdr

        tpeak_12p5 = max(smooth(pbcor,[1,1,5],/nan,/edge_wrap)*tpeak_mask, dim=3, /nan)
        tpeak_12p5[blank_ind] = !values.f_nan
        tpeak_hdr = twod_head(hdr)
        sxaddpar, tpeak_hdr, 'BUNIT', 'K'
        writefits, dir+ $
                   gal+'_co21_tpeak_12p5kms.fits', tpeak_12p5, tpeak_hdr
        
        mom0_hdr = twod_head(hdr)
        sxaddpar, mom0_hdr, 'BUNIT', 'K*KM/S'
        mom0[blank_ind] = !values.f_nan
        writefits, dir+$
                   gal+'_co21_mom0.fits', mom0, mom0_hdr
        writefits, dir+$
                   gal+'_co21_emom0.fits', e_mom0, mom0_hdr
        
        mom1_hdr = twod_head(hdr)
        sxaddpar, mom1_hdr, 'BUNIT', 'KM/S'
        mom1[blank_ind] = !values.f_nan
        writefits, dir+$
                   gal+'_co21_mom1.fits', mom1, mom1_hdr
        writefits, dir+$
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
; SHUFFLE AND CARRY OUT A WEIGHTED SHUFFLE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_shuffle) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'SHUFFLING CUBES', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     dir = release_dir+'process/'

     n_ext = n_elements(ext_in)

     for ii = 0, n_gals-1 do begin

        gal = gals[ii]

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        mom1 = $
           readfits(dir+gal+'_co21_mom1.fits', mom1_hdr)
        cube = $
           readfits(dir+gal+'_co21_correct.fits', cube_hdr)
        make_axes, cube_hdr, vaxis=vaxis, /vonly

        deltav = abs(vaxis[1]-vaxis[0])
        new_vaxis = findgen(201)*deltav*0.5
        new_vaxis -= mean(new_vaxis)
        
        shuffled_cube = $
           shuffle( $
           spec=cube $
           , vaxis=vaxis $
           , zero=mom1*1d3 $
           , target_vaxis=new_vaxis)

        shuffle_hdr = cube_hdr
        sxaddpar, shuffle_hdr, 'NAXIS3', n_elements(new_vaxis)
        sxaddpar, shuffle_hdr, 'CDELT3', new_vaxis[1]-new_vaxis[0]
        sxaddpar, shuffle_hdr, 'CRVAL3', new_vaxis[0]
        sxaddpar, shuffle_hdr, 'CRPIX3', 1
        
        writefits, dir+gal+'_co21_shuffle.fits' $
                   , shuffled_cube, shuffle_hdr  

        disp, max(shuffled_cube,dim=2,/nan)

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE TO SPECIFIC RESOLUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_conv_to_res) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLVING TO PHYSICAL RESOLUTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     

     dir = '../release/v0p5/process/'

     target_res = [45, 50, 60, 80, 100, 120, 150, 250, 500]
     n_res = n_elements(target_res)
     tol = 0.1

     for ii = 0, n_gals-1 do begin

        gal = gals[ii]

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        cube = readfits(dir+gal+'_co21_correct.fits', cube_hdr)
        s = gal_data(dirs[ii])
        sxaddpar, cube_hdr, 'DIST', s.dist_mpc, 'MPC'
        
        current_res = s.dist_mpc*!dtor*sxpar(cube_hdr, 'BMAJ')*1d6
        for jj = 0, n_res -1 do begin
           res_str = strcompress(str(target_res[jj]),/rem)+'pc'
           out_name = dir+strlowcase(gal)+'_co21_correct_'+res_str+'.fits'
           target_res_as = target_res[jj]/(s.dist_mpc*1d6)/!dtor*3600.d
           if current_res gt (1.0+tol)*target_res[jj] then begin
              print, "Resolution too coarse. Skipping."
              continue
           endif
           if abs(current_res - target_res[jj])/target_res[jj] lt tol then begin
              print, "I will call ", current_res, " ", target_res[jj]
              writefits, out_name, cube, cube_hdr           
           endif else begin
              print, "I will convolve ", current_res, " to ", target_res[jj]
              conv_with_gauss $
                 , data=cube $
                 , hdr=cube_hdr $
                 , target_beam=target_res_as*[1,1,0] $
                 , out_file=out_name              
           endelse
           
        endfor    

     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMPILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_compile) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COMPILE A RELEASE', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     dir = '../release/v0p5/'

     ext_in = $
        ['_co21_correct.fits' $
         ,'_co21_resid_round.fits' $
         ,'_co21_flat_round_smoothed_mask.fits' $
         ,'_co21_flat_round_mask.fits' $
         ,'_large_mask.fits' $
         ,'_co21_shuffle.fits' $
         ,'_co21_mom0.fits' $
         ,'_co21_emom0.fits' $
         ,'_co21_mom1.fits' $
         ,'_co21_emom1.fits' $
         ,'_co21_tpeak.fits' $
         ,'_co21_tpeak_12p5kms.fits' $
         ,'_tp_k.fits' $
        ]

     ext_out = $
        ['_co21.fits' $
         ,'_co21_resid.fits' $
         ,'_co21_smoothedmask.fits' $
         ,'_co21_brightmask.fits' $
         ,'_co21_widemask.fits' $
         ,'_co21_shuffled.fits' $
         ,'_co21_mom0.fits' $
         ,'_co21_emom0.fits' $
         ,'_co21_mom1.fits' $
         ,'_co21_emom1.fits' $
         ,'_co21_tpeak.fits' $
         ,'_co21_tpeak_12p5kms.fits' $
         ,'_co21_tp.fits' $
        ]

     n_ext = n_elements(ext_in)

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue
        
        gal = gals[ii]

        print, 'Copying '+gal

        for jj = 0, n_ext-1 do begin

           spawn, 'rm -rf '+dir+'delivery/'+gal+ext_out[jj]

           spawn, 'cp '+dir+'process/'+gal+ext_in[jj]+' '+dir+'delivery/'+gal+ext_out[jj]

        endfor

     endfor

  endif

end
