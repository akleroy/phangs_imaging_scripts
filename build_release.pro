pro build_release $
   , copy=copy $
   , stage=stage $
   , process=process $
   , convolve=do_convolve

;+
;
; Scripts to build the imaged data into a data release.
;
;-

; DIRECTORIES
  version = '0.4'
  vstring = 'v0p4'
  root_imaging_dir = '../'
  release_dir = root_imaging_dir+'release/'+vstring+'/'

; GALAXIES
  gals = $
     ['ngc0628' $
      , 'ngc1672' $
      , 'ngc3351' $
      , 'ngc3627' $
      , 'ngc4254' $
      , 'ngc4303' $
      , 'ngc4321' $
      , 'ngc4535' $
      , 'ngc5068' $
      , 'ngc6744' $
     ]
  n_gals = n_elements(gals)

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(copy) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COPYING DATA', /info     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     message, 'I am copying the raw data to build a release.', /info
     message, '... deleting and remaking the raw data directory', /info

     spawn, 'rm -rf '+release_dir+'raw/'
     spawn, 'mkdir '+release_dir+'raw/'

     ext_to_copy = $
        ['_co21_pb.fits' $
         , '_co21_residual.fits' $
         , '_co21_round.fits' $
         , '_co21_round_pbcor.fits']
     n_ext = n_elements(ext_to_copy)

     for ii = 0, n_gals-1 do begin
        message, '... copying data for '+gals[ii], /info        
        for jj = 0, n_ext-1 do begin
           spawn, 'cp '+root_imaging_dir+gal+'/'+ $
                  gals[ii]+ext_to_copy[jj]+'.fits '+ $
                  release_dir+'raw/.'
        endfor
     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED COPYING DATA', /info     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PROCESS DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(stage) then begin

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; Clean up NGC 6744
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=          

     ext = ['_cube' $
            , '_cube_pb' $
            , '_cube_pbcor' $
            , '_cube_residual' $
            , '_taper' $
            , '_taper_pb' $
            , '_taper_pbcor' $
            , '_taper_residual']
     out_ext = ['_cube_round' $
                , '_cube_pb' $
                , '_cube_round_pbcor' $
                , '_cube_residual' $
                , '_taper_round' $
                , '_taper_pb' $
                , '_taper_round_pbcor' $
                , '_taper_residual']
     n_ext = n_elements(ext)

     for ii = 0, n_ext-1 do begin

        for jj = 0, 2 do begin

           if jj eq 0 then $
              line = 'co21'
           if jj eq 1 then $
              line = 'c18o21'

           north = readfits(release_dir+ $
                            'raw/ngc6744north_'+line+ext[ii]+'_align.fits' $
                            , north_hdr)

           south = readfits(release_dir+ $
                            'raw/ngc6744south_'+line+ext[ii]+'_align.fits' $
                            , south_hdr)
           
           combo = north
           sz = size(north)
           if sz[2] gt 1300 then $
              combo[*,0:750,*] = south[*,0:750,*] $
           else $
              combo[*,0:600,*] = south[*,0:600,*]
           
           sxdelpar, north_hdr, 'BLANK'
           writefits, release_dir+'raw/ngc6744_'+line+out_ext[ii]+'.fits' $
                      , combo, north_hdr

           loadct, 33
           disp, max(combo, dim=3, /nan), /sq

        endfor

     endfor

     ext = ['' $
            , '_pb' $
            , '_pbcor' $
            , '_residual' $
            , '_taper' $
            , '_taper_pb' $
            , '_taper_pbcor' $
            , '_taper_residual']
     out_ext = ['_round' $
                , '_pb' $
                , '_round_pbcor' $
                , '_residual' $
                , '_taper_round' $
                , '_taper_pb' $
                , '_taper_round_pbcor' $
                , '_taper_residual']
     n_ext = n_elements(ext)

     for ii = 0, n_ext-1 do begin

        for jj = 0, 3 do begin

           if jj eq 0 then $
              line = 'chan0_co21'
           if jj eq 1 then $
              line = 'chan0_c18o21'
           if jj eq 2 then $
              line = 'cont'

           north = readfits(release_dir+ $
                            'raw/ngc6744north_'+line+ext[ii]+'_align.fits' $
                            , north_hdr)

           south = readfits(release_dir+ $
                            'raw/ngc6744south_'+line+ext[ii]+'_align.fits' $
                            , south_hdr)
           
           combo = north
           sz = size(north)
           if sz[2] gt 1300 then $
              combo[*,0:750] = south[*,0:750] $
           else $
              combo[*,0:600] = south[*,0:600]
           
           sxdelpar, north_hdr, 'BLANK'
           writefits, release_dir+'raw/ngc6744_'+line+out_ext[ii]+'.fits' $
                      , combo, north_hdr

           loadct, 33
           disp, combo, /sq

        endfor

     endfor

  endif

; TARGET RESOLUTIONS
  target_res = [60, 80, 100, 150]
  n_res = n_elements(target_res)
  res_tol = 0.1

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PROCESS CUBES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Loop over all of our targets and blank regions outside the main
; field of view. Convert to Kelvin. Note conversion in header.

  if keyword_set(process) then begin

     for ii = 0, n_gals-1 do begin

        gal = gals[ii]

        for jj = 0, 1 do begin
           
           if jj eq 0 then ext = '_cube'
           if jj eq 1 then ext = '_taper'
           
;--------------------------------------------------------------------------
; CO 2-1 AND C18O 2-1
;--------------------------------------------------------------------------

           for kk = 0, 2 do begin

              if kk eq 0 then $
                 line = 'co21' 
              if kk eq 1 then begin
;                NGC 0628 lacks C18O
                 if gal eq 'ngc0628' then $
                    continue
                 line = 'c18o21'
              endif
              
;             Read from disk

              pbcor_cube = readfits(release_dir+ $
                                    'raw/'+gal+'_'+line+ext+'_round_pbcor.fits', pbcor_hdr)
              cube = readfits(release_dir+ $
                              'raw/'+gal+'_'+line+ext+'_round.fits', hdr)
              pb = readfits(release_dir+ $
                            'raw/'+gal+'_'+line+ext+'_pb.fits', pb_hdr)

              loadct, 33
              disp, max(cube, dim=3, /nan), /sq
              contour, max(pb, dim=3, /nan), /overplot, lev=[0.2, 0.5]

;             Mask regions outside the main mosaic coverage

              ind = where(pb lt 0.5 or finite(pb) eq 0)
              pbcor_cube[ind] = !values.f_nan
              cube[ind] = !values.f_nan

;             Convert to kelvin

              jtok = calc_jtok(hdr=hdr)
              cube *= jtok
              pbcor_cube *= jtok

              sxaddpar, hdr, 'BUNIT', 'K'
              sxaddpar, pbcor_hdr, 'BUNIT', 'K'

;             Pare down to only a small set of NaN values outside the
;             main cube.

              mask_2d = total(finite(cube), 3) ge 1
              mask_x = total(mask_2d, 2) ge 1
              ind_x = where(mask_x)
              mask_y = total(mask_2d, 1) ge 1
              ind_y = where(mask_y)

              sz = size(mask_2d)
              xlo = (min(ind_x)-10) > 0
              xhi = (max(ind_x)+10) < (sz[1]-1)
              ylo = (min(ind_y)-10) > 0
              yhi = (max(ind_y)+10) < (sz[2]-1)

              new_pbcor_cube = $
                 cube_hextract(cube_in=pbcor_cube $
                               , hdr_in=pbcor_hdr $
                               , hdr_out=new_pbcor_hdr $
                               , x0=xlo, x1=xhi $
                               , y0=ylo, y1=yhi)
              new_cube = $
                 cube_hextract(cube_in=cube $
                               , hdr_in=hdr $
                               , hdr_out=new_hdr $
                               , x0=xlo, x1=xhi $
                               , y0=ylo, y1=yhi)
              new_pb = $
                 cube_hextract(cube_in=pb $
                               , hdr_in=pb_hdr $
                               , hdr_out=new_pb_hdr $
                               , x0=xlo, x1=xhi $
                               , y0=ylo, y1=yhi)

;             Write to disk

              writefits, release_dir+ $
                         'delivery/'+strupcase(gal)+'_'+line+ext+'.fits' $
                         , new_cube, new_hdr
              writefits, release_dir+ $
                         'delivery/'+strupcase(gal)+'_'+line+ext+'_pbcor.fits' $
                         , new_pbcor_cube, new_pbcor_hdr
              writefits, release_dir+ $
                         'delivery/'+strupcase(gal)+'_'+line+ext+'_pb.fits' $
                         , new_pb, new_pb_hdr

           endfor

        endfor

        for jj = 0, 1 do begin
           
           if jj eq 0 then ext = ''
           if jj eq 1 then ext = '_taper'

; ------------------------------------------------------------------------
; CONTINUUM
; -------------------------------------------------------------------------
           
           for kk = 0, 3 do begin

              if kk eq 0 then $
                 line = 'chan0_co21' 
              if kk eq 1 then begin
;                NGC 0628 lacks C18O
                 if gal eq 'ngc0628' then $
                    continue
                 line = 'chan0_c18o21'
              endif
              if kk eq 2 then begin
                 line = 'cont'
              endif

;             Read from disk

              pbcor_im = readfits(release_dir+ $
                                  'raw/'+gal+'_'+line+ext+'_round_pbcor.fits' $
                                  , pbcor_hdr)
              im = readfits(release_dir+ $
                            'raw/'+gal+'_'+line+ext+'_round.fits', hdr)
              pb = readfits(release_dir+ $
                            'raw/'+gal+'_'+line+ext+'_pb.fits', pb_hdr)

              loadct, 33
              disp, im, /sq
              contour, pb, /overplot, lev=[0.2, 0.5]

;             Mask regions outside the main mosaic coverage

              ind = where(pb lt 0.5 or finite(pb) eq 0)
              pbcor_im[ind] = !values.f_nan
              im[ind] = !values.f_nan

;             Leave these in Jy/beam. Later, we need to handle the
;             chan0 right to get K km/s. For now, the velocity channel
;             width has been dropped.

;             Pare down to only a small set of NaN values outside the
;             main cube.

              mask_2d = finite(im) ge 1
              mask_x = total(mask_2d, 2) ge 1
              ind_x = where(mask_x)
              mask_y = total(mask_2d, 1) ge 1
              ind_y = where(mask_y)

              sz = size(mask_2d)
              xlo = (min(ind_x)-10) > 0
              xhi = (max(ind_x)+10) < (sz[1]-1)
              ylo = (min(ind_y)-10) > 0
              yhi = (max(ind_y)+10) < (sz[2]-1)

              hextract, pbcor_im, pbcor_hdr, new_pbcor_im, new_pbcor_hdr $
                        , xlo, xhi, ylo, yhi

              hextract, im, hdr, new_im, new_hdr $
                        , xlo, xhi, ylo, yhi

              hextract, pb, pb_hdr, new_pb, new_pb_hdr $
                        , xlo, xhi, ylo, yhi

;             Write to disk

              writefits, release_dir+'delivery/'+strupcase(gal)+'_'+ $
                         line+ext+'.fits', new_im, new_hdr
              writefits, release_dir+'delivery/'+strupcase(gal)+'_'+ $
                         line+ext+'_pbcor.fits', new_pbcor_im, new_pbcor_hdr
              writefits, release_dir+'delivery/'+strupcase(gal)+'_'+line+ $
                         ext+'_pb.fits', new_pb, new_pb_hdr
              
           endfor

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; CONVOLVE CUBES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_convolve) then begin
     
     in_dir = release_dir+'delivery/'
     out_dir = release_dir+'delivery/'

     for ii = 0, n_gals-1 do begin

        gal = strupcase(gals[ii])

        for ww = 0, 1 do begin

           if ww eq 0 then begin
              ext = 'cube'
           endif

           if ww eq 1 then begin
              ext = 'cube_pbcor'
           endif

           cube = readfits(in_dir+gal+'_co21_'+ext+'.fits', cube_hdr)
           s = gal_data(gal)
           sxaddpar, cube_hdr, 'DIST', s.dist_mpc, 'MPC'        
           sxaddpar, cube_hdr, 'VERSION', version, 'SFNG Version'        
           current_res = s.dist_mpc*!dtor*sxpar(cube_hdr, 'BMAJ')*1d6

           for jj = 0, n_res -1 do begin

              res_str = strcompress(str(target_res[jj]),/rem)+'pc'
              out_name = out_dir+gal+'_co21_'+ext+'_'+res_str+'.fits'
              target_res_as = target_res[jj]/(s.dist_mpc*1d6)/!dtor*3600.d

              if target_res_as gt 3.0 then begin
                 print, "Switching to tapered cube."
                 if ww eq 0 then $
                    cube = readfits(in_dir+gal+'_co21_taper.fits', cube_hdr)
                 if ww eq 1 then $
                    cube = readfits(out_dir+gal+'_co21_taper_pbcor.fits', cube_hdr)
                 sxaddpar, cube_hdr, 'DIST', s.dist_mpc, 'MPC'
                 sxaddpar, cube_hdr, 'VERSION', version, 'SFNG Version'        
              endif

              if current_res gt (1.0+res_tol)*target_res[jj] then begin
                 print, "Resolution too coarse. Skipping."
                 continue
              endif

              if abs(current_res - target_res[jj])/target_res[jj] lt res_tol then begin
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
        
     endfor

  endif


end
