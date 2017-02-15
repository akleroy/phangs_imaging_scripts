pro build_release $
   , reset=reset $
   , copy=copy $
   , rescale=rescale $
   , process=process $
   , pbcorr=pbcorr $
   , roundbeam=roundbeam $
   , convolve=do_convolve $
   , smooth=do_smooth $
   , mask=do_masks $
   , collapse=do_collapse $
   , compile=do_compile

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
  dirs = $
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
      , 'ngc6744' $
     ]
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
; RESIDUAL RESCALING
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(rescale) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'RESIDUAL RESCALING', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     for ii = 0, n_gals-1 do begin
        message, '... rescaling data for '+gals[ii], /info        
        
        dirty_cube = readfits(release_dir+'raw/'+$
                              gals[ii]+'_co21_dirty.fits', dirty_hdr)

        resid_cube = readfits(release_dir+'raw/'+$
                              gals[ii]+'_co21_residual.fits', resid_hdr)

        clean_cube = readfits(release_dir+'raw/'+$
                              gals[ii]+'_co21.fits', clean_hdr)

        pb =  readfits(release_dir+'raw/'+$
                       gals[ii]+'_co21_pb.fits', pb_hdr)
        
        mask = pb gt 0.5

        rms = mad(clean_cube[where(mask)])
        mask = clean_cube gt 3.*rms
        mask = mask*(shift(mask,0,0,-1) or shift(mask,0,0,1))

        ind = where(mask)
        y = (clean_cube[ind]-resid_cube[ind])
        x = (dirty_cube[ind]-resid_cube[ind])
        
        print, "Rescaling factor would be: "
        sum_scale_factor = total(y,/nan) / total(x,/nan)
        print, "... sum: ", sum_scale_factor

        bins = bin_data(x, y, /nan, xmin=0., xmax=max(x, /nan), binsize=max(x,/nan)/20.)
        bin_scale_factor = median(bins.ymed/bins.xmed)

        ploterror, bins.xmid, bins.ymed, bins.ymad $
                   , ytitle="C-R", xtitle="D-R", ps=7
        equality, slop=sum_scale_factor, color=cgcolor('red')
        equality, lines=2, slope=bin_scale_factor, color=cgcolor('blue')
        print, "... bins: ", bin_scale_factor
        
        scale_factor = float(bin_scale_factor)

        rescale = float(clean_cube + resid_cube*(scale_factor-1.0))
        sxaddpar, clean_hdr, 'RESCALE', scale_factor, 'Residual emission rescaled by this factor.'
        sxdelpar, clean_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+gals[ii]+'_co21_rescale.fits', rescale, clean_hdr

        rescale = float(resid_cube*(scale_factor-1.0))
        sxaddpar, resid_hdr, 'RESCALE', scale_factor, 'Residual emission rescaled by this factor.'
        sxdelpar, resid_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+gals[ii]+'_co21_residual_rescale.fits', rescale, resid_hdr

     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED RESIDUAL RESCALING', /info     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; PRIMARY BEAM CORRECTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(pbcorr) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'PRIMARY BEAM CORRECTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     for ii = 0, n_gals-1 do begin

        message, '... primary beam correction for '+gals[ii], /info        
        
        clean_cube = readfits(release_dir+'raw/'+$
                              gals[ii]+'_co21.fits', clean_hdr)

        rescale_cube = readfits(release_dir+'process/'+$
                                gals[ii]+'_co21_rescale.fits', rescale_hdr)

        pb =  readfits(release_dir+'raw/'+$
                       gals[ii]+'_co21_pb.fits', pb_hdr)
        
        maxpb = max(pb,dim=3,/nan)
        disp, maxpb
        contour, maxpb, lev=[0.5], /overplot

        ind = where(pb lt 0.5)
        
        clean_cube[ind] = !values.f_nan
        rescale_cube[ind] = !values.f_nan
        
        sxdelpar, clean_hdr, 'BLANK'
        sxdelpar, clean_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_clip.fits', clean_cube, clean_hdr

        sxdelpar, rescale_hdr, 'BLANK'
        sxdelpar, rescale_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_rescale_clip.fits', rescale_cube, rescale_hdr

        clean_cube /= pb
        rescale_cube /= pb

        sxdelpar, clean_hdr, 'BLANK'
        sxdelpar, clean_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_pbcorr.fits', clean_cube, clean_hdr

        sxdelpar, rescale_hdr, 'BLANK'
        sxdelpar, rescale_hdr, 'CASAMBM'
        writefits, release_dir+'process/'+$
                   gals[ii]+'_co21_rescale_pbcorr.fits', rescale_cube, rescale_hdr

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
        
        beam_table = mrdfits(release_dir+'raw/'+gals[ii]+'_co21.fits',1)

        ext_to_convolve = $
           ['_co21_clip' $
            , '_co21_rescale_clip' $
            , '_co21_pbcorr' $
            , '_co21_rescale_pbcorr']
        n_ext = n_elements(ext_to_convolve)
        
        for jj = 0, n_ext-1 do begin

           cube = float(readfits(release_dir+'process/'+gals[ii]+ext_to_convolve[jj]+'.fits', hdr))
           pix = abs(sxpar(hdr,'CDELT1'))
           target_bmaj = sqrt(max(beam_table.bmaj,/nan)^2+(pix*1.5)^2)
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
; CONVOLVE TO 3" RESOLUTION
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_smooth) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLUTION TO LOWER RESOLUTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_gals-1 do begin
        
        target_bmaj = 3.
        
        ext_to_convolve = $
           ['_co21_clip_round' $
            , '_co21_rescale_clip_round' $
            , '_co21_pbcorr_round' $
            , '_co21_rescale_pbcorr_round']
        n_ext = n_elements(ext_to_convolve)
        
        for jj = 0, n_ext-1 do begin

           dir = '../release/v0p4/process/'
           in_cube = dir+gals[ii]+ext_to_convolve[jj]+".fits"
           out_cube = dir+gals[ii]+ext_to_convolve[jj]+"_smoothed.fits"
           print, "convolving "+in_cube

           conv_with_gauss $
              , data=in_cube $
              , target = [1.,1.,0.]*target_bmaj $
              , out_file=out_cube

        endfor

     endfor
     
  endif
  

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD A COUPLE OF MASKS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_masks) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'BUILD MASKS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     for ii = 0, n_gals-1 do begin
        
        gal = gals[ii]
        
        ext = $
           ['_co21_clip_round_smoothed' $
            , '_co21_clip_round']
        n_ext = n_elements(ext)
        
        for jj = 0, n_ext-1 do begin

           cube = readfits('../release/v0p4/process/'+gal+ext[jj]+'.fits', cube_hdr)

           rms = mad(cube)
           make_cprops_mask $
              , indata=cube $
              , inrms = rms $
              , lo_thresh = 2 $
              , hi_thresh = 4 $
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
     
     for ii = 0, n_gals-1 do begin
        
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
; COMPILE
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_compile) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COMPILE A RELEASE', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
 
     dir = '../release/v0p4/'

     ext_in = $
        ['_co21_correct.fits' $
;         ,'_co21_pbcorr_round_smoothed.fits' $
         ,'_co21_clip_round_smoothed_mask.fits' $
         ,'_co21_clip_round_mask.fits' $
         ,'_co21_mom0.fits' $
         ,'_co21_emom0.fits' $
         ,'_co21_mom1.fits' $
         ,'_co21_emom1.fits' $
         ,'_co21_tpeak.fits' $
         ,'_co21_tpeak_12p5kms.fits' $
        ]

     ext_out = $
        ['_co21.fits' $
;         , '_co21_smoothed.fits' $
         , '_co21_brightmask.fits' $
         , '_co21_smoothedmask.fits' $
         ,'_co21_mom0.fits' $
         ,'_co21_emom0.fits' $
         ,'_co21_mom1.fits' $
         ,'_co21_emom1.fits' $
         ,'_co21_tpeak.fits' $
         ,'_co21_tpeak_12p5kms.fits' $
        ]

     n_ext = n_elements(ext_in)

     for ii = 0, n_gals-1 do begin
        
        gal = gals[ii]

        print, 'Copying '+gal

        for jj = 0, n_ext-1 do begin

           spawn, 'rm -rf '+dir+'delivery/'+gal+ext_out[jj]

           spawn, 'cp '+dir+'process/'+gal+ext_in[jj]+' '+dir+'delivery/'+gal+ext_out[jj]

        endfor

     endfor

  endif

end
