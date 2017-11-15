pro build_cubes $
   , version=version $
   , only=only $
   , skip=skip $
   , just_array=just_array $
   , reset=do_reset $
   , copy=do_copy $
   , pbcorr=do_pbcorr $
   , roundbeam=do_roundbeam $
   , singledish=do_singledish $
   , stage=do_stage_feather $
   , feather=do_copy_feather $
   , merge=do_merge $
   , sanitize=do_sanitize $
   , convolve=do_conv_to_res $
   , target_res=target_res

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

; RESOLUTIONS

  if n_elements(target_res) eq 0 then begin
     target_res = [45, 60, 80, 100, 120, 500, 750, 1000]
  endif
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

        if n_elements(only) gt 0 then $
           if total(only eq gals[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq gals[ii]) gt 0 then continue

        message, '... copying data for '+gals[ii], /info        

        this_gal = gals[ii]
        this_dir = dir_for_gal[ii]

        for jj = 0, n_array - 1 do begin

           this_array = array_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue
           
           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_copy[ll]

                 in_file = root_imaging_dir+this_dir+'/'+this_gal+$
                           '_'+this_array+'_'+this_product+$
                           this_ext

                 out_file = release_dir+'raw/.'

                 test = file_search(in_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+in_file+' not found.', /info
                    continue
                 endif
                 
                 spawn, 'cp '+ in_file +' '+ out_file
                 
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
        
        if n_elements(only) gt 0 then $
           if total(only eq gals[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq gals[ii]) gt 0 then continue

        this_gal = gals[ii]
        this_dir = dir_for_gal[ii]

        message, "Applying primary beam correction for "+this_gal, /info

        for jj = 0, n_array - 1 do begin

           this_array = array_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue           

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]
              
              cube_file = $
                 release_dir+'raw/'+$
                 this_gal+'_'+this_array+$
                 '_'+this_product+'.fits'
              
              test = file_search(cube_file, count=found)
              if found eq 0 then begin
                 message, this_gal+" file "+cube_file+" not found.", /info
                 continue
              endif

              flat_cube = readfits(cube_file $
                                   , flat_hdr)

              pb_cube = readfits(release_dir+'raw/'+$
                                 this_gal+'_'+this_array+$
                                 '_'+this_product+'_pb.fits' $
                                 , pb_hdr)

;             Primary beam correct and blank outside the primary beam
;             limit. Determine the primary beam limit from the primary
;             beam cube.

              pb_limit = min(pb_cube, /nan)
              
              blank_ind = where(pb_cube lt pb_limit)
              flat_cube[blank_ind] = !values.f_nan
              
;             Write to the products directory.

              sxaddpar, flat_hdr, 'ARRAY', this_array
              writefits, release_dir+'process/'+$
                         this_gal+'_'+this_array+$
                         '_'+this_product+'_flat.fits' $
                         , flat_cube, flat_hdr
              
              pbcorr_cube = flat_cube/pb_cube
              pbcorr_hdr = flat_hdr

              writefits, release_dir+'process/'+$
                         this_gal+'_'+this_array+$
                         '_'+this_product+'_pbcorr.fits' $
                         , pbcorr_cube, pbcorr_hdr
              
           endfor
           
        endfor

     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED PRIMARY BEAM CORRECTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ROUND BEAM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
  if keyword_set(do_roundbeam) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'CONVOLUTION TO ROUND BEAM', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     ext_to_convolve = $
        ['flat' $
         , 'pbcorr']
     n_ext = n_elements(ext_to_convolve)
     
     for ii = 0, n_gals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq gals[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq gals[ii]) gt 0 then continue

        this_gal = gals[ii]
        this_dir = dir_for_gal[ii]

        message, "Applying primary beam correction for "+this_gal, /info

        for jj = 0, n_array - 1 do begin

           this_array = array_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue           

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              flat_infile = $
                 release_dir+'process/'+$
                 this_gal+'_'+this_array+'_'+$
                 this_product+'_flat.fits'
              
              test = file_search(flat_infile, count=found)
              if found eq 0 then begin
                 message, 'File '+flat_infile+' not found.', /info
                 continue
              endif

              hdr = headfits(flat_infile)
              
              bmaj = sxpar(hdr, 'BMAJ')*3600.
              pix = abs(sxpar(hdr,'CDELT1'))
              target_bmaj = sqrt(bmaj^2+(pix*2.*3600.)^2)

              print, "For "+flat_infile+" and related products, I target:"
              print, "BMAJ = ", target_bmaj

              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_convolve[ll]
                 
                 this_fname = $
                    release_dir+'process/'+$
                    this_gal+'_'+this_array+'_'+$
                    this_product+'_'+this_ext+'.fits'

                 test = file_search(this_fname, count=found)
                 if found eq 0 then begin
                    message, 'File '+this_fname+' not found.', /info
                    continue
                 endif

                 cube = float(readfits(this_fname, hdr))
                 print, "convolving "+this_fname

                 sz = size(cube)
                 conv_with_gauss $
                    , data=cube $
                    , hdr=hdr $
                    , target =[1,1,0.]*target_bmaj $
                    , out_data=out_cube $
                    , out_hdr=out_hdr $
                    , /perbeam

                 this_outfile = $
                    release_dir+'process/'+$
                    this_gal+'_'+this_array+'_'+$
                    this_product+'_'+this_ext+'_round.fits'

                 writefits, this_outfile, out_cube, out_hdr

              endfor

           endfor

        endfor
        
     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED CONVOLUTION TO ROUND BEAM', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY AND TAPER THE TOTAL POWER WITH THE PRIMARY BEAM
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_singledish) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COPYING AND PROCESSING TOTAL POWER DATA', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     

     readcol $
        , 'singledish_key.txt' $
        , format='A,A,A', comment='#' $
        , sd_gal, sd_fname, sd_product

     for ii = 0, n_gals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq gals[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq gals[ii]) gt 0 then continue

        this_gal = gals[ii]

        for kk = 0, n_product-1 do begin
           
           this_product = product_list[kk]

           sd_ind = where(sd_gal eq this_gal and $
                          sd_product eq this_product $
                          , sd_ct)

           if sd_ct eq 0 then begin
              message, 'I did not find a single dish key entry for '+ $
                       this_gal+' '+this_product, /info
              continue
           endif 

           message, "Found a single dish file for "+this_gal+' '+this_product, /info

                                ; Copy  to the process/ directory

           out_file = release_dir+'process/'+this_gal+'_tp_'+this_product+'.fits'
           spawn, 'rm -rf '+ out_file
           spawn, 'cp '+sd_fname[sd_ind]+' '+out_file

                                ; Make a couple weaks to the file

           cube = readfits(out_file, hdr)
           sxaddpar, hdr, 'ARRAY', 'TP'
           writefits, out_file, cube, hdr
           
                                ; Make a Kelvin version

           jtok = calc_jtok(hdr=hdr)
           hdr_k = hdr
           sxaddpar, hdr_k, 'JTOK', jtok
           cube_k = cube*jtok
           sxaddpar, hdr_k, 'BUNIT', 'K'

           writefits $
              , release_dir+'process/'+this_gal+'_tp_'+this_product+'_k.fits' $
              , cube_k, hdr_k

                                ; Now project a version in Jy/beam on to the astrometry of the
                                ; various interferometric cubes.

           for jj = 0, n_array - 1 do begin

              this_array = array_list[jj]
              if n_elements(just_array) gt 0 then $
                 if total(just_array eq this_array) eq 0 then continue           

              pb_file = $
                 release_dir+'raw/'+$
                 this_gal+'_'+this_array+ $
                 '_'+this_product+'_pb.fits'
              
              test = file_search(pb_file, count=found)
              if found eq 0 then begin
                 message, this_gal+" file "+pb_file+ $
                          " not found. Skipping TP alignment and tapering step.", /info
                 continue
              endif

              pb_cube =  readfits(pb_file, pb_hdr)
              
              cube_hastrom $
                 , data = cube $
                 , hdr_in = hdr $
                 , outcube = aligned_tp $
                 , outhdr = new_hdr_tp $
                 , target_hdr = pb_hdr
              
              writefits $
                 , release_dir+'process/'+ $
                 this_gal+'_tp_'+this_product+ $
                 '_aligned_'+this_array+'.fits' $
                 , aligned_tp, new_hdr_tp
              
              tapered_tp = pb_cube * aligned_tp

              writefits $
                 , release_dir+'process/'+ $
                 this_gal+'_tp_'+this_product+ $
                 '_tapered_'+this_array+'.fits' $
                 , tapered_tp, new_hdr_tp
              
           endfor
           
        endfor

     endfor
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED COPYING AND PROCESSING TOTAL POWER DATA', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     
     
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

     for jj = 0, n_array - 1 do begin
        
        this_array = array_list[jj]
        
        if n_elements(just_array) gt 0 then $
           if total(just_array eq this_array) eq 0 then continue           
        
        for kk = 0, n_product-1 do begin
           
           this_product = product_list[kk]

           gal_list = []

           for ii = 0, n_gals-1 do begin
              
              if n_elements(only) gt 0 then $
                 if total(only eq gals[ii]) eq 0 then continue
              
              if n_elements(skip) gt 0 then $
                 if total(skip eq gals[ii]) gt 0 then continue
              
              this_gal = gals[ii]
              this_dir = dir_for_gal[ii]
              
              interf_infile = $
                 release_dir+'process/'+this_gal+'_'+this_array+ $
                 '_'+this_product+'_flat_round.fits'
              
              test = file_search(interf_infile, count=found)
              if found eq 0 then begin
                 message, 'File '+interf_infile+' not found. Skipping feathering.', /info
                 continue
              endif

              tp_infile = $
                 release_dir+'process/'+this_gal+'_tp_'+ $
                 this_product+'_tapered_'+this_array+'.fits'

              test = file_search(tp_infile, count=found)
              if found eq 0 then begin
                 message, 'File '+tp_infile+' not found. Skipping feathering.', /info
                 continue
              endif

              gal_list = [gal_list, this_gal]

              spawn, 'cp -r '+interf_infile+' '+release_dir+'feather/.'
              spawn, 'cp -r '+tp_infile+' '+release_dir+'feather/.'

           endfor

           if n_elements(gal_list) eq 0 then $
              continue

;         Write the feathering script

           get_lun, lun
           openw, lun, release_dir+'feather/feather_script_'+this_array+'_'+ $
                  this_product+'.py'

           printf, lun, 'gal_list = ['
           for ii = 0, n_elements(gal_list)-1 do $
              printf, lun, '    "'+gal_list[ii]+'",'
           printf, lun, '    ]'

           printf, lun,'for gal in gal_list:'
           printf, lun,'    print "Importing "+gal'
           printf, lun,'    importfits(fitsimage=gal+"_'+this_array+'_'+this_product+'_flat_round.fits",'
           printf, lun,'               imagename=gal+"_'+this_array+'_'+this_product+'_flat_round.image",'
           printf, lun,'               zeroblanks=True, overwrite=True)'
           printf, lun,'    importfits(fitsimage=gal+"_tp_'+this_product+'_tapered_'+this_array+'.fits",'
           printf, lun,'               imagename=gal+"_tp_'+this_product+'_tapered_'+this_array+'.image",'
           printf, lun,'               zeroblanks=True,overwrite=True)'

           printf, lun,'for gal in gal_list:'
           printf, lun,'    print "Feathering "+gal    '
           printf, lun,'    feather(imagename=gal+"_'+this_array+'_'+this_product+'_feathered.image",'
           printf, lun,'            highres=gal+"_'+this_array+'_'+this_product+'_flat_round.image",'
           printf, lun,'            lowres=gal+"_tp_'+this_product+'_tapered_'+this_array+'.image")'

           printf, lun,'for gal in gal_list:'
           printf, lun,'    print "Exporting "+gal'
           printf, lun,'    exportfits(imagename=gal+"_'+this_array+'_'+this_product+'_feathered.image",'
           printf, lun,'               fitsimage=gal+"_'+this_array+'_'+this_product+'_feathered.fits",'
           printf, lun,'               velocity=True, dropdeg=True, dropstokes=True,'
           printf, lun,'               overwrite=True, bitpix=16)'

           close, lun        

        endfor

     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED STAGING FEATHER CALL - NOW RUN THE SCRIPT!', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FEATHERED AND SINGLE DISH DATA INTO THE DIRECTORY
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_copy_feather) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'COPYING FEATHERED DATA TO PROCESSED DIRECTORY', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

     in_dir = release_dir+'feather/'
     out_dir = release_dir+'process/'
     
     for ii = 0, n_gals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq gals[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq gals[ii]) gt 0 then continue

        this_gal = gals[ii]
        this_dir = dir_for_gal[ii]

        for jj = 0, n_array - 1 do begin

           this_array = array_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue           

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              feathered_infile = $
                 release_dir+'feather/'+$
                 this_gal+'_'+this_array+'_'+$
                 this_product+'_feathered.fits'
              
              test = file_search(feathered_infile, count=found)
              if found eq 0 then begin
                 message, 'File '+feathered_infile+' not found.', /info
                 continue
              endif
              
              cube = readfits(feathered_infile, hdr)
              
              interf_infile = $
                 release_dir+'process/'+this_gal+'_'+this_array+ $
                 '_'+this_product+'_flat_round.fits'            

              interf_cube = readfits(interf_infile, interf_hdr)
              
              blank_ind = where(finite(interf_cube) eq 0 or $
                                abs(cube - sxpar(hdr,'BLANK')) lt 1d-6, blank_ct)
              if blank_ct gt 0 then $
                 cube[blank_ind] = !values.f_nan

              sxaddpar, hdr, 'ARRAY', strupcase(this_array)+'+TP'

              writefits, $
                 out_dir+this_gal+'_'+this_array+'+tp_'+this_product+ $
                 '_flat_round.fits' $
                 , cube, hdr

              pb_cube = readfits(release_dir+'raw/'+$
                                 this_gal+'_'+this_array+$
                                 '_'+this_product+'_pb.fits' $
                                 , pb_hdr)
              cube = cube/pb_cube

              writefits, $
                 out_dir+this_gal+'_'+this_array+'+tp_'+this_product+ $
                 '_pbcorr_round.fits' $
                 , cube, hdr
              
           endfor
           
        endfor
        
     endfor

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED COPYING FEATHERED DATA.', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
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

        if n_elements(only) gt 0 then $
           if total(only eq gals[ii]) eq 0 then continue

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

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED MERGING MULTI-PART MOSAICS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

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

        if n_elements(only) gt 0 then $
           if total(only eq gals[ii]) eq 0 then continue

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

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'FINISHED SANITIZING CUBES.', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info

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

        if n_elements(only) gt 0 then $
           if total(strlowcase(only) eq strlowcase(gals[ii])) eq 0 then continue

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

end
