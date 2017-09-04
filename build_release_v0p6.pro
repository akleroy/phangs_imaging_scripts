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
   , merge=do_merge $
   , sanitize=do_sanitize $
   , clean_mask=do_clean_mask $
   , convolve=do_conv_to_res $
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
      , 'ngc1365' $
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
;     , 'ngc1365' $
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
;     , 'ngc5128' $
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

  target_res = [45, 60, 80, 100, 120, 500, 750, 1000]
  n_res = n_elements(target_res)

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

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue
        
        message, "Applying primary beam correction for "+gals[ii], /info

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

        message, 'Copying data from feathered '+gals[ii], /info

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

                 template = readfits(out_dir+galname[jj]+ $
                                     '_co21_7m_flat_round.fits', temp_hdr)
                 blank_ind = where(finite(template) eq 0 or $
                                   abs(cube - sxpar(hdr,'BLANK')) lt 1d-6, blank_ct)
                 if blank_ct gt 0 then $
                    cube[blank_ind] = !values.f_nan

                 sxaddpar, hdr, 'ARRAY', '7M+TP'
                 writefits, $
                    out_dir+galname[jj]+'_co21_7m+tp_flat_round.fits' $
                    , cube, hdr
                 pb_cube = $
                    readfits(release_dir+'raw/'+$
                             galname[jj]+'_co21'+array+'_pb.fits', pb_hdr)
                 cube = cube/pb_cube
                 writefits, $
                    out_dir+galname[jj]+'_co21_7m+tp_pbcorr_round.fits' $
                    , cube, hdr
              endif

              if kk eq 1 then begin
                 if keyword_set(only_7m) then $
                    continue
                 if total(gals[ii] eq has_12m) eq 0B then $
                    continue
                 array = '_12m+7m'
                 cube = $
                    readfits(in_dir+galname[jj]+ $
                             '_co21'+array+'_feathered.fits', hdr)

                 template = readfits(out_dir+galname[jj]+ $
                                     '_co21_12m+7m_flat_round.fits', temp_hdr)
                 blank_ind = where(finite(template) eq 0 or $
                                   abs(cube - sxpar(hdr,'BLANK')) lt 1d-7, blank_ct)
                 if blank_ct gt 0 then $
                    cube[blank_ind] = !values.f_nan

                 sxaddpar, hdr, 'ARRAY', '12M+7M+TP'
                 writefits, $
                    out_dir+galname[jj]+'_co21_12m+7m+tp_flat_round.fits' $
                    , cube, hdr
                 pb_cube = $
                    readfits(release_dir+'raw/'+$
                             galname[jj]+'_co21_12m+7m_pb.fits', pb_hdr)
                 cube = cube/pb_cube
                 writefits, $
                    out_dir+galname[jj]+'_co21_12m+7m+tp_pbcorr_round.fits' $
                    , cube, hdr
              endif

           endfor

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MERGE MULTI-FIELD CUBES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_merge) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'MERGING MULTI-PART MOSAICS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     message, '', /info
     message, 'Merging two part cubes into a single cube.', /info
     message, '', /info
     
     dir = release_dir+'process/'
     
     readcol $
        , 'twopart_fields.txt', format='A,F,F,F,F,I' $
        , merge_name, merge_ra, merge_dec, merge_dra, merge_ddec, merge_copy_tp $
        , comment='#'

     in_dir = release_dir+'process/'
     out_dir = release_dir+'process/'

     for ii = 0, n_gals-1 do begin

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        if total(gals[ii] eq two_part) eq 0 then begin                 
           continue
        endif else begin
           galname = [gals[ii]+'north', gals[ii]+'south']
        endelse
        message, 'Merging '+gals[ii], /info

        message, 'Merging multiple parts for '+gals[ii], /info

        for kk = 0, 4 do begin
           
           if kk eq 0 then begin
              array = '_7m'
           endif

           if kk eq 1 then begin
              array = '_7m+tp'
           endif
           
           if kk eq 2 then begin
              array = '_12m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           if kk eq 3 then begin
              array = '_12m+7m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           if kk eq 4 then begin
              array = '_12m+7m+tp'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; LOOP OVER FILES OF INTEREST
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
           
           ext = $
              ['_flat_round' $
               , '_pbcorr_round' $
              ]
           n_ext = n_elements(ext)           
           
           for jj = 0, n_elements(ext)-1 do begin

              if n_elements(galname) ne 2 then begin
                 message, "I expect two parts but did not find two parts", /info
              endif

              part1 = $
                 readfits(in_dir+galname[0]+ $
                          '_co21'+array+ext[jj]+'.fits', part1_hdr)

              part2 = $
                 readfits(in_dir+galname[1]+ $
                          '_co21'+array+ext[jj]+'.fits', part2_hdr)
              
              if array eq '_7m' or array eq '_7m+tp' then begin
                 pb_cube1 =  readfits(release_dir+'raw/'+$
                                      galname[0]+'_co21_7m_pb.fits', pb1_hdr)
                 pb_cube2 =  readfits(release_dir+'raw/'+$
                                      galname[1]+'_co21_7m_pb.fits', pb2_hdr)
              endif else if array eq '_12m' then begin
                 pb_cube1 =  readfits(release_dir+'raw/'+$
                                      galname[0]+'_co21_12m_pb.fits', pb1_hdr)
                 pb_cube2 =  readfits(release_dir+'raw/'+$
                                      galname[1]+'_co21_12m_pb.fits', pb2_hdr)
              endif else begin
                 pb_cube1 =  readfits(release_dir+'raw/'+$
                                      galname[0]+'_co21_12m+7m_pb.fits', pb1_hdr)
                 pb_cube2 =  readfits(release_dir+'raw/'+$
                                      galname[1]+'_co21_12m+7m_pb.fits', pb2_hdr)
              endelse

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; MATCH THE BEAMS BETWEEN THE TWO PARTS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

              beam1 = sxpar(part1_hdr, 'BMAJ')
              beam2 = sxpar(part2_hdr, 'BMAJ')              

              pix = abs(sxpar(part2_hdr,'CDELT1'))
              target_bmaj = sqrt((beam1 > beam2)^2+pix^2)*3600.
              
              conv_with_gauss $
                 , data=part1 $
                 , hdr=part1_hdr $
                 , target =[1,1,0.]*target_bmaj $
                 , out_data=new_part1 $
                 , out_hdr=new_part1_hdr $
                 , /perbeam $
                 , worked=worked
              if worked eq 0 then begin
                 message, 'Problem with convolution in merging.', /info
                 stop
              endif 

              nan_ind = where(finite(part1) eq 0, nan_ct)
              if nan_ct gt 0 then new_part1[nan_ind] = !values.f_nan
              part1 = new_part1
              part1_hdr = new_part1_hdr

              conv_with_gauss $
                 , data=part2 $
                 , hdr=part2_hdr $
                 , target =[1,1,0.]*target_bmaj $
                 , out_data=new_part2 $
                 , out_hdr=new_part2_hdr $
                 , /perbeam $
                 , worked=worked
              if worked eq 0 then begin
                 message, 'Problem with convolution in merging.', /info
                 stop
              endif 
              
              nan_ind = where(finite(part2) eq 0, nan_ct)
              if nan_ct gt 0 then new_part2[nan_ind] = !values.f_nan
              part2 = new_part2
              part2_hdr = new_part2_hdr

              target_hdr = part1_hdr

              print, "Beam 1, Beam 2, Target Beam: ", beam1*3600., beam2*3600, target_bmaj
              print, "Beam in target header: ", sxpar(target_hdr, 'bmaj')*3600.

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; GENERATE A NEW TARGET HEADER AND ALIGN TO THE TARGET
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

              merge_ind = where(gals[ii] eq strcompress(merge_name,/rem) $
                                , merge_ct)
              if merge_ct eq 0 then begin
                 message, 'Failed to find a match in the merge table. Stopping.', /info
                 stop
              endif
              cdelt = abs(sxpar(target_hdr,'CDELT1'))

              npix_ra = abs(ceil(merge_dra[merge_ind]/3600. / cdelt))
              crpix_ra = npix_ra/2.0

              npix_dec = abs(ceil(merge_ddec[merge_ind]/3600. / cdelt))
              crpix_dec = npix_dec/2.0

              sxaddpar, target_hdr, 'CTYPE1', 'RA---SIN'
              sxdelpar, target_hdr, 'CRVAL1'
              sxaddpar, target_hdr, 'CRVAL1', double((merge_ra[merge_ind]*1.0)[0])
              sxdelpar, target_hdr, 'NAXIS1'              
              sxaddpar, target_hdr, 'NAXIS1', long(npix_ra[0]), after='NAXIS'
              sxdelpar, target_hdr, 'CRPIX1'
              sxaddpar, target_hdr, 'CRPIX1', crpix_ra[0]*1.0

              sxaddpar, target_hdr, 'CTYPE2', 'DEC--SIN'  
              sxdelpar, target_hdr, 'CRVAL2'
              sxaddpar, target_hdr, 'CRVAL2', double((merge_dec[merge_ind]*1.0)[0])
              sxdelpar, target_hdr, 'NAXIS2'
              sxaddpar, target_hdr, 'NAXIS2', long(npix_dec[0]), after='NAXIS1'
              sxdelpar, target_hdr, 'CRPIX2'
              sxaddpar, target_hdr, 'CRPIX2', crpix_dec[0]*1.0
              
              cube_hastrom $
                 , data = part1 $
                 , hdr_in = part1_hdr $
                 , target_hdr = target_hdr $
                 , outcube = new_part1 $
                 , outhdr = new_part1_hdr $
                 , missing=!values.f_nan

              cube_hastrom $
                 , data = pb_cube1 $
                 , hdr_in = pb1_hdr $
                 , target_hdr = target_hdr $
                 , outcube = new_pbcube1 $
                 , outhdr = new_pbcube1_hdr $
                 , missing=!values.f_nan
              pb_cube1 = new_pbcube1

              cube_hastrom $
                 , data = part2 $
                 , hdr_in = part2_hdr $
                 , target_hdr = target_hdr $
                 , outcube = new_part2 $
                 , outhdr = new_part2_hdr $
                 , missing=!values.f_nan

              cube_hastrom $
                 , data = pb_cube2 $
                 , hdr_in = pb2_hdr $
                 , target_hdr = target_hdr $
                 , outcube = new_pbcube2 $
                 , outhdr = new_pbcube2_hdr $
                 , missing=!values.f_nan
              pb_cube2 = new_pbcube2

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; MERGE THE TWO PARTS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

              cube_out = new_part1*!values.f_nan
              
              cov_part1 = finite(new_part1) and pb_cube1 gt 0.
              cov_part2 = finite(new_part2) and pb_cube2 gt 0.

              ind = where(cov_part1 eq 1 and cov_part2 eq 0, ct)
              if ct gt 0 then $
                 cube_out[ind] = new_part1[ind]
              
              ind = where(cov_part1 eq 0 and cov_part2 eq 1, ct)
              if ct gt 0 then $
                 cube_out[ind] = new_part2[ind]
              
              ind = where(cov_part1 eq 1 and cov_part2 eq 1, ct)
              if ct gt 0 then $
                 cube_out[ind] = (new_part1[ind]*pb_cube1[ind]^2 + $
                                  new_part2[ind]*pb_cube2[ind]^2) / $
                                 (pb_cube1[ind]^2 + pb_cube2[ind]^2)
              
              !p.multi=[0,3,1]
              disp, max(cube_out, dim=3, /nan), max=1
              disp, max(new_part1, dim=3, /nan), max=1
              disp, max(new_part2, dim=3, /nan), max=1
              
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; WRITE TO DISK
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

              writefits, out_dir+gals[ii]+'_co21'+ $
                         array+ext[jj]+'.fits', cube_out $
                         , target_hdr

           endfor

; MERGE THE TOTAL POWER - USUALLY JUST A COPYING STEP

           part1 = readfits(in_dir+galname[0]+'_tp_k.fits', part1_hdr)
           part2 = readfits(in_dir+galname[1]+'_tp_k.fits', part2_hdr)
           target_hdr = part1_hdr
           
           merge_ind = where(gals[ii] eq strcompress(merge_name,/rem) $
                             , merge_ct)
           if merge_ct eq 0 then begin
              message, 'Failed to find a match in the merge table. Stopping.', /info
              stop
           endif
           if merge_copy_tp[merge_ind] then begin
              writefits, out_dir+gals[ii]+'_tp_k.fits' $
                         , part1, part1_hdr
           endif else begin
              target_hdr = part1_hdr
              cdelt = abs(sxpar(target_hdr,'CDELT1'))

              npix_ra = abs(ceil(merge_dra[merge_ind]/3600. / cdelt))
              crpix_ra = npix_ra/2.0

              npix_dec = abs(ceil(merge_ddec[merge_ind]/3600. / cdelt))
              crpix_dec = npix_dec/2.0

              sxaddpar, target_hdr, 'CTYPE1', 'RA---SIN'
              sxdelpar, target_hdr, 'CRVAL1'
              sxaddpar, target_hdr, 'CRVAL1', double((merge_ra[merge_ind]*1.0)[0])
              sxdelpar, target_hdr, 'NAXIS1'              
              sxaddpar, target_hdr, 'NAXIS1', long(npix_ra[0]), after='NAXIS'
              sxdelpar, target_hdr, 'CRPIX1'
              sxaddpar, target_hdr, 'CRPIX1', crpix_ra[0]*1.0

              sxaddpar, target_hdr, 'CTYPE2', 'DEC--SIN'  
              sxdelpar, target_hdr, 'CRVAL2'
              sxaddpar, target_hdr, 'CRVAL2', double((merge_dec[merge_ind]*1.0)[0])
              sxdelpar, target_hdr, 'NAXIS2'
              sxaddpar, target_hdr, 'NAXIS2', long(npix_dec[0]), after='NAXIS1'
              sxdelpar, target_hdr, 'CRPIX2'
              sxaddpar, target_hdr, 'CRPIX2', crpix_dec[0]*1.0
              
              cube_hastrom $
                 , data = part1 $
                 , hdr_in = part1_hdr $
                 , target_hdr = target_hdr $
                 , outcube = new_part1 $
                 , outhdr = new_part1_hdr $
                 , missing=!values.f_nan
              
              cube_hastrom $
                 , data = part2 $
                 , hdr_in = part2_hdr $
                 , target_hdr = target_hdr $
                 , outcube = new_part2 $
                 , outhdr = new_part2_hdr $
                 , missing=!values.f_nan
              
              cube_out = new_part1*!values.f_nan
              
              cov_part1 = finite(new_part1)
              cov_part2 = finite(new_part2)

              ind = where(cov_part1 eq 1 and cov_part2 eq 0, ct)
              if ct gt 0 then $
                 cube_out[ind] = new_part1[ind]
              
              ind = where(cov_part1 eq 0 and cov_part2 eq 1, ct)
              if ct gt 0 then $
                 cube_out[ind] = new_part2[ind]
              
              ind = where(cov_part1 eq 1 and cov_part2 eq 1, ct)
              if ct gt 0 then $
                 cube_out[ind] = (new_part1[ind] + new_part2[ind])/2.0
              
              writefits, out_dir+gals[ii]+'_tp_k.fits' $
                         , cube_out, target_hdr

           endelse           

        endfor
        
     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; SANITIZE CUBES - TRIM EMPTY SPACE, ETC.
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_sanitize) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'SANITIZING CUBES FOR FINAL PRODUCT CONSTRUCTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_gals-1 do begin

        dir = release_dir+'process/'

        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        message, "Sanitizing the cubes for "+gals[ii], /info

        for kk = 0, 4 do begin
           
           if kk eq 0 then begin
              array = '_7m'
           endif

           if kk eq 1 then begin
              array = '_7m+tp'
           endif
           
           if kk eq 2 then begin
              array = '_12m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           if kk eq 3 then begin
              array = '_12m+7m'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           if kk eq 4 then begin
              array = '_12m+7m+tp'
              if keyword_set(only_7m) then $
                 continue
              if total(gals[ii] eq has_12m) eq 0B then $
                 continue
           endif

           ext_to_process = $
              ['_flat_round' $
               , '_pbcorr_round']
           n_ext = n_elements(ext_to_process)
           
           if total(gals[ii] eq two_part) eq 0 then begin                 
              galname = [gals[ii]]
           endif else begin
              galname = [gals[ii], gals[ii]+'north', gals[ii]+'south']
           endelse

           for zz = 0, n_elements(galname)-1 do begin

              for jj = 0, n_ext-1 do begin

                 in_cube = dir+galname[zz]+'_co21'+array+ext_to_process[jj]+".fits"
                 out_file = dir+galname[zz]+'_co21'+array+ext_to_process[jj]+"_trimmed.fits"
                 
                 message, '... cube '+in_cube, /info
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
                 
                 cube = readfits(out_file, hdr)           
                 jtok = calc_jtok(hdr=hdr)
                 cube *= jtok
                 sxaddpar, hdr, 'BUNIT', 'K'
                 writefits $
                    , dir+galname[zz]+'_co21'+array+ext_to_process[jj]+"_k.fits" $
                    , cube, hdr

              endfor

           endfor
           
        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MAKE CLEAN MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_clean_mask) then begin
     
     for ii = 0, n_gals-1 do begin
        
        message, "Making clean mask for "+gals[ii], /info

        dir = release_dir+'process/'
        
        if n_elements(just) gt 0 then $
           if total(just eq gals[ii]) eq 0 then continue

        cube = readfits(dir+gals[ii]+'_co21_7m+tp_flat_round_k.fits', hdr)
        
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
        disp, max(cube, dim=3, /nan), /xs, /ys
        contour, total(mask,3,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

        loadct, 33
        disp, max(cube, dim=2, /nan), /xs, /ys
        contour, total(mask,2,/nan) gt 0, /overplot, lev=[1], color=cgcolor('white')

        loadct, 33
        disp, max(cube, dim=1, /nan), /xs, /ys
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
        
        print, "Clean mask for "+gals[ii]+". Key to continue."
        test = get_kbrd(1)

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
     tol = 0.1

     s = gal_data(gals)

     for ii = 0, n_gals-1 do begin

        dir = release_dir+'process/'

        message, '', /info
        message, "Convolving the cubes for "+gals[ii], /info
        message, '', /info

        if n_elements(just) gt 0 then $
           if total(strlowcase(just) eq strlowcase(gals[ii])) eq 0 then continue

        for kk = 0, 1 do begin

           if kk eq 0 then begin
              array = '_7m+tp'
           endif
           
           if kk eq 1 then begin
              array = '_12m+7m+tp'
              if keyword_set(only_7m) then $
                 continue
              if total(strlowcase(gals[ii]) eq strlowcase(has_12m)) eq 0B then $
                 continue
           endif
           
           print, "ARRAY == ", array

           gal = gals[ii]
           
           ext_to_process = $
              ['_flat_round_k' $
               , '_pbcorr_round_k']
           n_ext = n_elements(ext_to_process)
           
           for jj = 0, n_ext-1 do begin                   
              
              cube = readfits(dir+strlowcase(gal)+ $
                              '_co21'+array+ext_to_process[jj]+'.fits', hdr)           
              sxaddpar, hdr, 'DIST', s[ii].dist_mpc, 'MPC / USED IN CONVOLUTION'           
              current_res_pc = s[ii].dist_mpc*!dtor*sxpar(hdr, 'BMAJ')*1d6
              
              for zz = 0, n_res -1 do begin
                 
                 res_str = strcompress(str(target_res[zz]),/rem)+'pc'
                 out_name = dir+strlowcase(gal)+ $
                            '_co21'+array+ext_to_process[jj]+'_'+res_str+'.fits'
                 target_res_as = target_res[zz]/(s[ii].dist_mpc*1d6)/!dtor*3600.d
                 
                 if current_res_pc gt (1.0+tol)*target_res[zz] then begin
                    print, strupcase(gal)+": Resolution too coarse. Skipping."
                    continue
                 endif
                 
                 if abs(current_res_pc - target_res[zz])/target_res[zz] lt tol then begin
                    print, strupcase(gal)+": I will call ", current_res_pc, " ", target_res[zz]
                    writefits, out_name, cube, hdr
                 endif else begin                 
                    print, strupcase(gal)+": I will convolve ", current_res_pc $
                           , " to ", target_res[zz]
                    conv_with_gauss $
                       , data=cube $
                       , hdr=hdr $
                       , target_beam=target_res_as*[1,1,0] $
                       , out_file=out_name
                 endelse
                 
              endfor

           endfor
           
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
           if total(strlowcase(just) eq strlowcase(gals[ii])) eq 0 then continue

        for kk = 0, 1 do begin

           if kk eq 0 then begin
              array = '_7m+tp'
           endif
           
           if kk eq 1 then begin
              array = '_12m+7m+tp'
              if keyword_set(only_7m) then $
                 continue
              if total(strlowcase(gals[ii]) eq strlowcase(has_12m)) eq 0B then $
                 continue
           endif
           
           print, "ARRAY == ", array

           gal = gals[ii]
           
           ext_to_process = $
              ['_flat_round_k' $
               , '_pbcorr_round_k']
           n_ext = n_elements(ext_to_process)
           
           for jj = 0, n_ext-1 do begin

              fname = file_search(dir+strlowcase(gal)+'_co21'+array $
                                  +ext_to_process[jj]+'.fits', count=fct)

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
              
              outfile = dir+strlowcase(gal)+'_co21'+array $
                        +ext_to_process[jj]+'_noise.fits'
              writefits, outfile, rms_cube, cube_hdr

              for zz = 0, n_res - 1 do begin

                 res_str = strcompress(str(target_res[zz]),/rem)+'pc'
                 print, "Resolution "+res_str

                 fname = file_search(dir+strlowcase(gal)+'_co21'+array $
                                     +ext_to_process[jj]+ $
                                     '_'+res_str+'.fits', count=fct)
                 
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

                 outfile = dir+strlowcase(gal)+'_co21'+array $
                           +ext_to_process[jj]+'_'+res_str+'_noise.fits'
                 writefits, outfile, rms_cube, cube_hdr

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

     dir = release_dir+'process/'
     
     for ii = 0, n_gals-1 do begin

        message, 'Building a mask for '+gals[ii], /info

        if n_elements(just) gt 0 then $
           if total(strlowcase(just) eq strlowcase(gals[ii])) eq 0 then continue

        for kk = 0, 1 do begin

           if kk eq 0 then begin
              array = '_7m+tp'
           endif
           
           if kk eq 1 then begin
              array = '_12m+7m+tp'
              if keyword_set(only_7m) then $
                 continue
              if total(strlowcase(gals[ii]) eq strlowcase(has_12m)) eq 0B then $
                 continue
           endif
           
           print, "ARRAY == ", array
           
           gal = gals[ii]
           
           ext_to_process = $
              ['_flat_round_k']
           n_ext = n_elements(ext_to_process)
           
           for jj = 0, n_ext-1 do begin

              for zz = 0, n_res do begin

                 if zz lt n_res then begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                    print, "Resolution "+res_str
                 endif else begin
                    res_str = ''
                    print, "Native resolution"
                 endelse
                 
                 fname = file_search(dir+strlowcase(gal)+'_co21'+array $
                                     +ext_to_process[jj] $
                                     +res_str+'.fits', count=fct)
                 
                 if fct eq 0 then continue

                 rms_name = $
                    file_search(dir+strlowcase(gal)+'_co21'+array $
                                +ext_to_process[jj] $
                                +res_str+'_noise.fits', count=noise_ct)

                 if noise_ct eq 0 then continue

                 cube = readfits(fname, cube_hdr)
                 
                 rms_cube = readfits(rms_name, noise_hdr)

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
                 writefits, dir+strlowcase(gal)+'_co21'+array $
                            +'_mask' $
                            +res_str+'.fits' $
                            , mask, mask_hdr
                 
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

     for ii = 0, n_gals-1 do begin

        message, 'Collapsing '+gals[ii], /info

        if n_elements(just) gt 0 then $
           if total(strlowcase(just) eq strlowcase(gals[ii])) eq 0 then continue

        for kk = 0, 1 do begin

           if kk eq 0 then begin
              array = '_7m+tp'
           endif
           
           if kk eq 1 then begin
              array = '_12m+7m+tp'
              if keyword_set(only_7m) then $
                 continue
              if total(strlowcase(gals[ii]) eq strlowcase(has_12m)) eq 0B then $
                 continue
           endif
           
           print, "ARRAY == ", array
           
           gal = gals[ii]
           
           ext_to_process = $
              ['_pbcorr_round_k']
           n_ext = n_elements(ext_to_process)
           
           for jj = 0, n_ext-1 do begin

              for zz = 0, n_res do begin

                 if zz lt n_res then begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                    print, "Resolution "+res_str
                 endif else begin
                    res_str = ''
                    print, "Native resolution"
                 endelse
                 
                 fname = file_search(dir+strlowcase(gal)+'_co21'+array $
                                     +ext_to_process[jj] $
                                     +res_str+'.fits', count=fct)
                 
                 if fct eq 0 then continue

                 rms_name = $
                    file_search(dir+strlowcase(gal)+'_co21'+array $
                                +ext_to_process[jj] $
                                +res_str+'_noise.fits', count=noise_ct)

                 if noise_ct eq 0 then continue

                 cube = readfits(fname, cube_hdr)
                 
                 rms_cube = readfits(rms_name, noise_hdr)

                 ppbeam = calc_pixperbeam(hdr=cube_hdr)         
                 
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

; Want to use a mask bigger than some size, deal with that here.

;                 if then begin

;                 endif else begin

;                 endelse
                 
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 mask = readfits(dir+strlowcase(gal)+'_co21'+array $
                                 +'_mask' $
                                 +res_str+'.fits' $
                                 , mask, mask_hdr)

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

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; PEAK TEMPERATURE IN ONE AND FIVE CHANNELS
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 tpeak_mask = mask
                 sz = size(mask)
                 for pp = 0, sz[3]-1 do $
                    tpeak_mask[*,*,pp] = total(mask[*,*,pp]) ge 1

                 tpeak = max(cube*tpeak_mask, dim=3, /nan)
;                 tpeak[blank_ind] = !values.f_nan
                 tpeak_hdr = twod_head(cube_hdr)
                 sxaddpar, tpeak_hdr, 'BUNIT', 'K'
                 
                 writefits, dir+strlowcase(gal)+'_co21'+array $
                            +'_tpeak'+res_str+'.fits' $
                            , tpeak, tpeak_hdr

                 tpeak_12p5 = max(smooth(cube,[1,1,5],/nan,/edge_wrap)*tpeak_mask, dim=3, /nan)
;                 tpeak_12p5[blank_ind] = !values.f_nan
                 tpeak_hdr = twod_head(cube_hdr)
                 sxaddpar, tpeak_hdr, 'BUNIT', 'K'

                 writefits, dir+strlowcase(gal)+'_co21'+array $
                            +'_tpeak12p5kms'+res_str+'.fits' $
                            , tpeak_12p5, tpeak_hdr

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; MOMENT 0
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 mom0_hdr = twod_head(cube_hdr)
                 sxaddpar, mom0_hdr, 'BUNIT', 'K*KM/S'
;                 mom0[blank_ind] = !values.f_nan

                 writefits, dir+strlowcase(gal)+'_co21'+array $
                            +'_mom0'+res_str+'.fits' $
                            , mom0, mom0_hdr

                 writefits, dir+strlowcase(gal)+'_co21'+array $
                            +'_emom0'+res_str+'.fits' $
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
;                 mom1[blank_ind] = !values.f_nan

                 writefits, dir+strlowcase(gal)+'_co21'+array $
                            +'_mom1'+res_str+'.fits' $
                            , mom1, mom1_hdr

                 writefits, dir+strlowcase(gal)+'_co21'+array $
                            +'_emom1'+res_str+'.fits' $
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

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COMPILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_compile) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COMPILE A RELEASE', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     in_dir = release_dir+'process/'
     out_dir = release_dir+'delivery/'

     ext_to_process = $
        ['_flat_round_k' $
         ,'_pbcorr_round_k' $
         ,'_mask' $
         ,'_mom0' $
         ,'_emom0' $
         ,'_mom1' $
         ,'_emom1' $
         ,'_tpeak' $
         ,'_tpeak12p5kms' $
        ]

     n_ext = n_elements(ext_to_process)

     for ii = 0, n_gals-1 do begin

        message, '-------------------------------------', /info
        message, 'Compiling files for '+gals[ii], /info
        message, '-------------------------------------', /info

        if n_elements(just) gt 0 then $
           if total(strlowcase(just) eq strlowcase(gals[ii])) eq 0 then continue

        for kk = 0, 1 do begin

           if kk eq 0 then begin
              array = '_7m+tp'
           endif
           
           if kk eq 1 then begin
              array = '_12m+7m+tp'
              if keyword_set(only_7m) then $
                 continue
              if total(strlowcase(gals[ii]) eq strlowcase(has_12m)) eq 0B then $
                 continue
           endif
           
           gal = gals[ii]
                      
           for jj = 0, n_ext-1 do begin

              for zz = 0, n_res do begin

                 if zz lt n_res then begin
                    res_str = '_'+strcompress(str(target_res[zz]),/rem)+'pc'
                 endif else begin
                    res_str = ''
                 endelse
                 
                 file_name = $
                    strlowcase(gal)+'_co21'+array $
                    +ext_to_process[jj] $
                    +res_str+'.fits'
                 fname = file_search(in_dir+file_name, count=fct)
                 
                 if fct eq 0 then continue

                 infile = fname
                 outfile = out_dir+file_name
                 
                 spawn, 'rm -rf '+outfile
                 spawn, 'cp '+infile+' '+outfile
                 print, '... ', outfile

;                LAZY BUT SHOULD WORK TO ALSO COPY NOISE FILES

                 file_name = $
                    strlowcase(gal)+'_co21'+array $
                    +ext_to_process[jj] $
                    +res_str+'_noise.fits'
                 fname = file_search(in_dir+file_name, count=fct)
                 
                 if fct eq 0 then continue

                 infile = fname
                 outfile = out_dir+file_name
                 
                 spawn, 'rm -rf '+outfile
                 spawn, 'cp '+infile+' '+outfile
                 print, '... ', outfile

              endfor

           endfor
           
        endfor
        
     endfor

  endif

end
