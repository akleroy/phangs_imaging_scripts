pro build_release_v0p6 $
   , just=just $
   , only_7m=only_7m $
   , reset=do_reset $
   , copy=do_copy $
   , pbcorr=do_pbcorr $
   , roundbeam=do_roundbeam $
   , taper_tp=do_taper_tp $
   , stage=do_stage_feather $
   , feather=do_copy_feather $
   , special=do_special $
   , sanitize=do_sanitize $
   , correct=do_correct $
   , heracles=do_heracles_smooth $
   , test=do_test $
   , smooth=do_smooth $
   , noise=do_noise $
   , mask=do_masks $
   , collapse=do_collapse $
   , shuffle=do_shuffle $
   , convolve=do_conv_to_res $
   , compile=do_compile

;+
;
; Scripts to build the imaged data into a data release.
;
;-

; DIRECTORIES
  version = '0.6'
  vstring = 'v0p6'
  root_imaging_dir = '../'
  release_dir = root_imaging_dir+'release/'+vstring+'/'

; GALAXIES
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
      , 'ngc6744' $
     ]

  has_12m = $
     ['ic5332' $
      , 'ngc0628' $
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
     ]

  two_part = $
     ['ngc3627' $
      , 'ngc4254' $
      , 'ngc4321' $
      , 'ngc5068' $
      , 'ngc6744' $
     ]

  n_gals = n_elements(gals)

  mom1_thresh = 2.5d
  mom0_thresh = 3.0d

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; RESET THE DIRECTORY STRUCTURE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_reset) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'RESETTING RELEASE', /info     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

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
        ['_pb.fits' $
         , '_residual.fits' $
         , '.fits' $
         , '_dirty.fits']

     n_ext = n_elements(ext_to_copy)

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        message, '... copying data for '+gals[ii], /info        
        
        for kk = 0, 2 do begin
           
           if kk eq 0 then begin
              array = '_7m'
           endif

           if kk eq 1 then begin
              array = '_12m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           if kk eq 2 then begin
              array = '_12m+7m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           for jj = 0, n_ext-1 do begin
              
              if total(gals[ii] eq two_part) eq 0 then begin                 
                 galname = [gals[ii]]
              endif else begin
                 galname = [gals[ii]+'north', gals[ii]+'south']
              endelse
              
              for zz = 0, n_elements(galname)-1 do begin                     
                 spawn, 'cp '+root_imaging_dir+gals[ii]+'/'+ $
                        galname[zz]+'_co21'+array+ext_to_copy[jj]+' '+ $
                        release_dir+'raw/.'
              endfor

           endfor

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

     for ii = 0, n_gals-1 do begin

        message, "Applying primary beam correction for "+gals[ii], /info

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue
        
        for kk = 0, 2 do begin
           
           if kk eq 0 then begin
              array = '_7m'
           endif

           if kk eq 1 then begin
              array = '_12m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           if kk eq 2 then begin
              array = '_12m+7m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif
           
           pb_limit = 0.25

           if total(gals[ii] eq two_part) eq 0 then begin                 
              galname = [gals[ii]]
           endif else begin
              galname = [gals[ii]+'north', gals[ii]+'south']
           endelse

           for jj = 0, n_elements(galname)-1 do begin

              flat_cube = readfits(release_dir+'raw/'+$
                                   galname[jj]+'_co21'+array+'.fits', flat_hdr)

              pb_cube =  readfits(release_dir+'raw/'+$
                                  galname[jj]+'_co21'+array+'_pb.fits', pb_hdr)

              blank_ind = where(pb_cube lt pb_limit)
              flat_cube[blank_ind] = !values.f_nan

              writefits, release_dir+'process/'+$
                         galname[jj]+'_co21'+array+'_flat.fits' $
                         , flat_cube, flat_hdr
              
              pbcorr_cube = flat_cube/pb_cube
              pbcorr_hdr = flat_hdr
              writefits, release_dir+'process/'+$
                         galname[jj]+'_co21'+array+'_pbcorr.fits' $
                         , pbcorr_cube, pbcorr_hdr
              
           endfor

        endfor     

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
        
        for kk = 0, 2 do begin
           
           if kk eq 0 then begin
              array = '_7m'
           endif

           if kk eq 1 then begin
              array = '_12m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           if kk eq 2 then begin
              array = '_12m+7m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif
           
           if total(gals[ii] eq two_part) eq 0 then begin                 
              galname = [gals[ii]]
           endif else begin
              galname = [gals[ii]+'north', gals[ii]+'south']
           endelse

           for jj = 0, n_elements(galname)-1 do begin
              
              hdr = headfits(release_dir+'process/'+$
                             galname[jj]+'_co21'+array+'_flat.fits')
              
              bmaj = sxpar(hdr, 'BMAJ')*3600.
              pix = abs(sxpar(hdr,'CDELT1'))
              target_bmaj = sqrt(bmaj^2+(pix*2.*3600.)^2)
              
              ext_to_convolve = $
                 ['_flat' $
                  , '_pbcorr']
              n_ext = n_elements(ext_to_convolve)
              
              for zz = 0, n_ext-1 do begin

                 fname = $
                    release_dir+'process/'+galname[jj]+ $
                    '_co21'+array+ext_to_convolve[zz]+ $
                    '.fits'
                 cube = float(readfits(fname, hdr))
                 print, "convolving "+fname

                 sz = size(cube)
                 conv_with_gauss $
                    , data=cube $
                    , hdr=hdr $
                    , target =[1,1,0.]*target_bmaj $
                    , out_data=out_cube $
                    , out_hdr=out_hdr $
                    , /perbeam

                 outfile = $
                    release_dir+'process/'+galname[jj]+ $
                    '_co21'+array+ext_to_convolve[zz]+ $
                    '_round.fits'
                 writefits, outfile, out_cube, out_hdr

              endfor

           endfor

        endfor
        
     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY AND TAPER THE TOTAL POWER WITH THE PRIMARY BEAM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_taper_tp) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COPYING AND PROCESSING TOTAL POWER DATA', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     

     readcol $
        , 'singledish_key.txt' $
        , format='A,A', comment='#' $
        , sd_gal, sd_fname

     out_dir = release_dir+'process/'

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        if total(gals[ii] eq two_part) eq 0 then begin                 
           galname = [gals[ii]]
        endif else begin
           galname = [gals[ii]+'north', gals[ii]+'south']
        endelse

        for jj = 0, n_elements(galname)-1 do begin           

           sd_ind = where(sd_gal eq galname[jj], sd_ct)
           if sd_ct eq 0 then begin
              message, 'I did not find a single dish key entry for '+galname[jj], /info
              continue
           endif else begin
              spawn, 'rm -rf '+$
                     out_dir+galname[jj]+'_tp.fits'
              spawn, 'cp '+sd_fname[sd_ind]+' '+out_dir+galname[jj]+'_tp.fits'
           endelse

           cube = readfits(out_dir+galname[jj]+'_tp.fits', hdr)
           
           sxaddpar, hdr, 'ARRAY', 'TP'
           blank_ind = where(abs(cube - sxpar(hdr,'BLANK')) lt 1d-6, blank_ct)
           if blank_ct gt 0 then $
              cube[blank_ind] = !values.f_nan
           
           writefits, out_dir+galname[jj]+'_tp.fits', cube, hdr
           
           jtok = calc_jtok(hdr=hdr)
           hdr_k = hdr
           sxaddpar, hdr_k, 'JTOK', jtok
           cube_k = cube*jtok
           sxaddpar, hdr_k, 'BUNIT', 'K'

           writefits, out_dir+galname[jj]+'_tp_k.fits', cube_k, hdr_k

           for kk = 0, 1 do begin
              
              if kk eq 0 then begin
                 array = '_7m'
              endif

              if kk eq 1 then begin
                 if keyword_set(only_7m) then $
                    continue
                 if total(gals[ii] eq has_12m) eq 0B then $
                    continue
                 array = '_12m+7m'
              endif

              pb_cube =  readfits(release_dir+'raw/'+$
                                  galname[jj]+'_co21'+array+'_pb.fits', pb_hdr)
              
              cube_hastrom $
                 , data = cube $
                 , hdr_in = hdr $
                 , outcube = aligned_tp $
                 , outhdr = new_hdr_tp $
                 , target_hdr = pb_hdr
              
              writefits $
                 , out_dir+galname[jj]+'_tp_aligned'+array+'.fits' $
                 , aligned_tp, new_hdr_tp
              
              tapered_tp = pb_cube * aligned_tp

              writefits $
                 , out_dir+galname[jj]+'_tp_tapered'+array+'.fits' $
                 , tapered_tp, new_hdr_tp
              
           endfor

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

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        message, 'Copying data to feather '+gals[ii], /info

        if total(gals[ii] eq two_part) eq 0 then begin                 
           galname = [gals[ii]]
        endif else begin
           galname = [gals[ii]+'north', gals[ii]+'south']
        endelse

        for jj = 0, n_elements(galname)-1 do begin           
           
           for kk = 0, 1 do begin
              
              if kk eq 0 then begin
                 array = '_7m'
                 spawn, 'cp '+in_dir+galname[jj]+'_co21_7m_flat_round.fits '+$
                        out_dir+'.'                 
                 spawn, 'cp '+in_dir+galname[jj]+'_tp_tapered_7m.fits '+$
                        out_dir+'.'
              endif

              if kk eq 1 then begin
                 array = '_12m+7m'
                 if keyword_set(only_7m) then $
                    continue
                 if total(gals[ii] eq has_12m) eq 0B then $
                    continue
                 spawn, 'cp '+in_dir+galname[jj]+'_co21_12m+7m_flat_round.fits '+$
                        out_dir+'.'                 
                 spawn, 'cp '+in_dir+galname[jj]+'_tp_tapered_12m+7m.fits '+$
                        out_dir+'.'
              endif

           endfor

        endfor

     endfor

     spawn, 'cp feather_script_7m.py '+out_dir+'/.'
     spawn, 'cp feather_script_12m.py '+out_dir+'/.'

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FEATHERED AND SINGLE DISH DATA INTO THE DIRECTORY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_copy_feather) then begin

     in_dir = release_dir+'feather/'
     out_dir = release_dir+'process/'

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        message, 'Copying data to feather '+gals[ii], /info

        if total(gals[ii] eq two_part) eq 0 then begin                 
           galname = [gals[ii]]
        endif else begin
           galname = [gals[ii]+'north', gals[ii]+'south']
        endelse

        for jj = 0, n_elements(galname)-1 do begin           
           
           for kk = 0, 1 do begin
              
              if kk eq 0 then begin
                 array = '_7m'
                 cube = $
                    readfits(in_dir+galname[jj]+ $
                             '_co21'+array+'_feathered.fits', hdr)
                 sxaddpar, hdr, 'ARRAY', '7M+TP'
                 writefits, $
                    out_dir+galname[jj]+'_co21_7m+tp_flat_round_taper.fits' $
                    , cube, hdr
                 pb_cube = $
                    readfits(release_dir+'raw/'+$
                             galname[jj]+'_co21'+array+'_pb.fits', pb_hdr)
                 cube = cube/pb_cube
                 writefits, $
                    out_dir+galname[jj]+'_co21_7m+tp_pbcorr_round_taper.fits' $
                    , cube, hdr
              endif

              if kk eq 1 then begin
                 array = '_12m+7m'
                 cube = $
                    readfits(in_dir+galname[jj]+ $
                             '_co21'+array+'_feathered.fits', hdr)
                 sxaddpar, hdr, 'ARRAY', '12M+7M+TP'
                 writefits, $
                    out_dir+galname[jj]+'_co21_12m+7m+tp_flat_round_taper.fits' $
                    , cube, hdr
                 pb_cube = $
                    readfits(release_dir+'raw/'+$
                             galname[jj]+'_co21_12m+7m_pb.fits', pb_hdr)
                 cube = cube/pb_cube
                 writefits, $
                    out_dir+galname[jj]+'_co21_12m+7m+tp_pbcorr_round_taper.fits' $
                    , cube, hdr
              endif

           endfor

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEAL WITH SOME SPECIAL CASES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_special) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'DEALING WITH SPECIAL CASES', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     message, '', /info
     message, 'Merging NGC 6744 into a single cube.', /info
     message, '', /info
     
     dir = release_dir+'process/'
     
     ext = $
        ['_co21_flat_round' $
         , '_co21_pbcorr_round' $
         , '_co21_resid_round' $
         , '_co21_feather_pbcorr' $
        ]
     n_ext = n_elements(ext)

     for jj = 0, n_elements(ext)-1 do begin
        cube_north = readfits(dir+'ngc6744north'+ext[jj]+'.fits', hdr_north)
        if cube_north[0] eq -1 then continue
        cube_south = readfits(dir+'ngc6744south'+ext[jj]+'.fits', hdr_south)
        if cube_south[0] eq -1 then continue
        
        if sxpar(hdr_north,'BMAJ') lt sxpar(hdr_south,'BMAJ') then begin
           print, sxpar(hdr_north,'BMAJ'), sxpar(hdr_south,'BMAJ')
           conv_with_gauss $
              , data=cube_north $
              , hdr=hdr_north $
              , target_beam=sxpar(hdr_south, 'BMAJ')*[1,1,0] $
              , out_data=new_cube_north $
              , out_hdr = new_hdr_north
           cube_north = new_cube_north
           hdr_north = new_hdr_north
        endif else begin
           print, sxpar(hdr_north,'BMAJ'), sxpar(hdr_south,'BMAJ')
           conv_with_gauss $
              , data=cube_south $
              , hdr=hdr_south $
              , target_beam=sxpar(hdr_north, 'BMAJ')*3600.*[1,1,0] $
              , out_data=new_cube_south $
              , out_hdr = new_hdr_south
           cube_south = new_cube_south
           hdr_south = new_hdr_south
        endelse

;       Combine them
        new_hdr = hdr_north
        sxaddpar, new_hdr, 'CTYPE1', 'RA---SIN'
        sxaddpar, new_hdr, 'CRVAL1', 287.442083
        sxaddpar, new_hdr, 'NAXIS1', 2500
        sxaddpar, new_hdr, 'CRPIX1', long(floor(2500/2))

        sxaddpar, new_hdr, 'CTYPE2', 'DEC--SIN'  
        sxaddpar, new_hdr, 'CRVAL2', -63.857528
        sxaddpar, new_hdr, 'NAXIS2', 3000
        sxaddpar, new_hdr, 'CRPIX2', long(floor(3000/2))

        cube_hastrom $
           , data = cube_south $
           , hdr_in = hdr_south $
           , outcube = new_cube_south $
           , outhdr = new_hdr_south $
           , target_hdr = new_hdr

        cube_hastrom $
           , data = cube_north $
           , hdr_in = hdr_north $
           , outcube = new_cube_north $
           , outhdr = new_hdr_north $
           , target_hdr = new_hdr
        
        new_cube = new_cube_south
        new_cube[*,1500:*,*] = new_cube_north[*,1500:*,*]
        writefits, dir+'ngc6744'+ext[jj]+'.fits', new_cube, new_hdr

     endfor

;    Total power          
     cube_north = readfits(dir+'ngc6744north_tp_k.fits', hdr_north)     
     cube_south = readfits(dir+'ngc6744south_tp_k.fits', hdr_south)

;    Combine them
     new_hdr = hdr_north
     sxaddpar, new_hdr, 'CTYPE1', 'RA---SIN'
     sxaddpar, new_hdr, 'CRVAL1', 287.442083
     sxaddpar, new_hdr, 'NAXIS1', 130
     sxaddpar, new_hdr, 'CRPIX1', floor(150/2)
     
     sxaddpar, new_hdr, 'CTYPE2', 'DEC--SIN'  
     sxaddpar, new_hdr, 'CRVAL2', -63.857528
     sxaddpar, new_hdr, 'NAXIS2', 180
     sxaddpar, new_hdr, 'CRPIX2', floor(200/2)
     
     cube_hastrom $
        , data = cube_south $
        , hdr_in = hdr_south $
        , outcube = new_cube_south $
        , outhdr = new_hdr_south $
        , target_hdr = new_hdr

     cube_hastrom $
        , data = cube_north $
        , hdr_in = hdr_north $
        , outcube = new_cube_north $
        , outhdr = new_hdr_north $
        , target_hdr = new_hdr
     
     new_cube = new_cube_south
     new_cube[*,90:*,*] = new_cube_north[*,90:*,*]
     writefits, dir+'ngc6744_tp_k.fits', new_cube, new_hdr

;    Free up memory
     cube_north = -1
     cube_south = -1
     new_cube_north = -1
     new_cube_south = -1
     new_cube = -1

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SANITIZE CUBES - TRIM EMPTY SPACE, ETC.
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_sanitize) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'SANITIZING CUBES FOR FINAL PRODUCT CONSTRUCTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_final_gals-1 do begin
        
        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        gal = final_gals[ii]

        ext_to_trim = $
           ['_co21_flat_round' $
            , '_co21_pbcorr_round' $
            , '_co21_resid_round' $
            , '_co21_feather_pbcorr']
        n_ext = n_elements(ext_to_trim)
        
        for jj = 0, n_ext-1 do begin

           dir = release_dir+'process/'
           in_cube = dir+gal+ext_to_trim[jj]+".fits"
           out_file = dir+gal+ext_to_trim[jj]+"_trimmed.fits"
           print, "Trimming "+in_cube

           cube_trim $
              , data = in_cube $
              , outfile = out_file

           test = headfits(out_file)
           pix_per_beam = abs(sxpar(test, 'BMAJ') /  sxpar(test, 'CDELT1'))
           print, "I calculated "+str(pix_per_beam)+" pixels per beam."
           if pix_per_beam gt 6 then begin
              print, "... I will rebin."
              cube_hrebin $
                 , data = out_file $
                 , outfile = out_file $
                 , factor = 2
           endif                      

        endfor

        message, "I will now align the total power to the interferometer", /info
        
        tp = readfits(dir+gal+'_tp_k.fits', tp_hdr)
        target_hdr = headfits(out_file)
        cube_hastrom $
           , data = tp $
           , hdr_in = tp_hdr $
           , outcube = new_tp $
           , outhdr = new_tp_hdr $
           , target_hdr = target_hdr
        writefits, dir+gal+'_tp_k_align.fits', new_tp, new_tp_hdr

     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVERT TO KELVIN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_correct) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVERTING TO KELVIN', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     

     dir = release_dir+'process/'

     for ii = 0, n_final_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        
        gal = final_gals[ii]

        ext_to_convert = $
           ['_co21_flat_round' $
            , '_co21_pbcorr_round' $
            , '_co21_resid_round' $
            , '_co21_feather_pbcorr']
        n_ext = n_elements(ext_to_convert)
        
        for jj = 0, n_ext-1 do begin        

           cube = readfits(dir+gal+ext_to_convert[jj]+'_trimmed.fits', hdr)           
           jtok = calc_jtok(hdr=hdr)

           cube *= jtok
           sxaddpar, hdr, 'BUNIT', 'K'
           sxaddpar, hdr, 'JTOK', jtok
           
           writefits, dir+gal+ext_to_convert[jj]+'_correct.fits', cube, hdr
           
        endfor
        
     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE TO SPECIFIC RESOLUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_conv_to_res) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLVING TO PHYSICAL RESOLUTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     

     dir = release_dir+'process/'
     target_res = [60, 80, 120, 500, 750]
     n_res = n_elements(target_res)
     tol = 0.1

     for ii = 0, n_final_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        
        gal = final_gals[ii]

        ext_to_convolve = $
           ['_co21_flat_round' $
            , '_co21_pbcorr_round' $
            , '_co21_resid_round' $
            , '_co21_feather_pbcorr']
        n_ext = n_elements(ext_to_convolve)
        
        for jj = 0, n_ext-1 do begin        

           cube = readfits(dir+gal+ext_to_convolve[jj]+'_correct.fits', hdr)           
           s = gal_data(dirs[ii])
           sxaddpar, hdr, 'DIST', s.dist_mpc, 'MPC / USED IN CONVOLUTION'           
           current_res = s.dist_mpc*!dtor*sxpar(hdr, 'BMAJ')*1d6

           for kk = 0, n_res -1 do begin

              res_str = strcompress(str(target_res[kk]),/rem)+'pc'
              out_name = dir+strlowcase(gal)+ext_to_convolve[jj]+'_'+res_str+'.fits'
              target_res_as = target_res[kk]/(s.dist_mpc*1d6)/!dtor*3600.d

              if current_res gt (1.0+tol)*target_res[jj] then begin
                 print, strupcase(gal)+": Resolution too coarse. Skipping."
                 continue
              endif

              if abs(current_res - target_res[kk])/target_res[kk] lt tol then begin
                 print, strupcase(gal)+": I will call ", current_res, " ", target_res[kk]
                 writefits, out_name, cube, hdr           
              endif else begin                 
                 print, strupcase(gal)+": I will convolve ", current_res, " to ", target_res[kk]
                 conv_with_gauss $
                    , data=cube $
                    , hdr=hdr $
                    , target_beam=target_res_as*[1,1,0] $
                    , out_file=out_name              
              endelse
              
           endfor

        endfor

     endfor
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE TO THE HERACLES RESOLTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_heracles_smooth) then begin

     dir = release_dir+'process/'

     for ii = 0, n_final_gals-1 do begin
        
        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        
        gal = final_gals[ii]

        ext_to_convolve = $
           ['_co21_feather_pbcorr']
        n_ext = n_elements(ext_to_convolve)
        
        for jj = 0, n_ext-1 do begin        
           
           cube = readfits(dir+gal+ext_to_convolve[jj]+'_correct.fits', hdr)           
           s = gal_data(dirs[ii])
           current_res = sxpar(hdr,'BMAJ')*3600.

           n_res = 1
           for kk = 0, n_res -1 do begin

              out_name = dir+strlowcase(gal)+ext_to_convolve[jj]+'_heracles.fits'
              target_res_as = 13.3
              tol = 0.1

              if current_res gt (1.0+tol)*target_res_as then begin
                 print, strupcase(gal)+": Resolution too coarse. Skipping."
                 continue
              endif

              if abs(current_res - target_res_as)/target_res_as lt tol then begin
                 print, strupcase(gal)+": I will call ", current_res, " ", target_res_as
                 writefits, out_name, cube, hdr           
              endif else begin                 
                 print, strupcase(gal)+": I will convolve ", current_res, " to ", target_res_as
                 conv_with_gauss $
                    , data=cube $
                    , hdr=hdr $
                    , target_beam=target_res_as*[1,1,0] $
                    , out_file=out_name              
              endelse
              
           endfor
           
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
     
     for ii = 0, n_final_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        gal = final_gals[ii]

        ext_to_convolve = $
           ['_co21_flat_round_trimmed' $
            , '_co21_pbcorr_round_trimmed' $
            , '_co21_resid_round_trimmed' $
            , '_co21_feather_pbcorr_trimmed']
        n_ext = n_elements(ext_to_convolve)
        
        for jj = 0, n_ext-1 do begin

           dir = release_dir+'process/'
           in_cube = dir+gal+ext_to_convolve[jj]+".fits"
           out_cube = dir+gal+ext_to_convolve[jj]+"_smoothed.fits"
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

     for ii = 0, n_final_gals-1 do begin

        message, 'Estimating noise for '+final_gals[ii], /info

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        
        gal = final_gals[ii]

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
     
     for ii = 0, n_final_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        
        gal = final_gals[ii]
        
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
     
     for ii = 0, n_final_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        
        gal = final_gals[ii]

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

     for ii = 0, n_final_gals-1 do begin

        gal = final_gals[ii]

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue

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

     for ii = 0, n_final_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue
        
        gal = final_gals[ii]

        print, 'Copying '+gal

        for jj = 0, n_ext-1 do begin

           spawn, 'rm -rf '+dir+'delivery/'+gal+ext_out[jj]

           spawn, 'cp '+dir+'process/'+gal+ext_in[jj]+' '+dir+'delivery/'+gal+ext_out[jj]

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; INFORMATIONAL
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_test) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'TEST', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     

     dir = '../release/v0p5/process/'

     for ii = 0, n_final_gals-1 do begin

        gal = final_gals[ii]

        if n_elements(just) gt 0 then $
           if total(just eq final_gals[ii]) eq 0 then continue

        h = headfits(dir+gal+'_co21_pbcorr_round_trimmed.fits')
        s = gal_data(dirs[ii])
        current_res = s.dist_mpc*!dtor*sxpar(h, 'BMAJ')*1d6
        print, gal, ' starts at resolution ', current_res

     endfor
     
  endif

end
