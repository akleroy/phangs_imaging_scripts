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

; This part of the script assumes that the directory name is the
; galaxy name (or at least the name of the synthesized product created
; by combining the fields).

  if keyword_set(do_merge) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'MERGING MULTI-PART MOSAICS', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     message, '', /info
     message, 'Merging two part cubes into a single cube.', /info
     message, '', /info
     
     dir = release_dir+'process/'
     
     readcol $
        , 'multipart_fields.txt', format='A,F,F,F,F,I' $
        , merge_name, merge_ra, merge_dec, merge_dra, merge_ddec, merge_copy_tp $
        , comment='#'

     in_dir = release_dir+'process/'
     out_dir = release_dir+'process/'

     n_merge = n_elements(merge_name)

     ext_to_merge = $
        ['flat_round' $
         , 'pbcorr_round' $
        ]
     n_ext = n_elements(ext_to_merge)           
     
     for ii = 0, n_merge-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq merge_name[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq merge_name[ii]) gt 0 then continue
        
        this_merge_gal = merge_name[ii]

        gal_ind = where(dir_for_gal eq this_merge_gal, n_part)
        if n_part lt 2 then begin
           message, 'Found less than 2 parts for '+this_merge_gal, /info
           continue
        endif
        
        for kk = 0, n_product-1 do begin
           
           this_product = product_list[kk]

           for jj = 0, n_fullarray - 1 do begin
              
              this_array = fullarray_list[jj]
              if n_elements(just_array) gt 0 then $
                 if total(just_array eq this_array) eq 0 then continue           
              
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; GET A LIST OF FILES TO MERGE AND READ THEM IN
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
              
              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_merge[ll]
                 
                 cube1_infile = $
                    release_dir+'process/'+ $
                    gals[gal_ind[0]]+'_'+ $
                    this_array+'_'+this_product+'_'+$
                    this_ext+'.fits'
                 
                 test = file_search(cube1_infile, count=found)
                 if found eq 0 then begin
                    message, 'File '+cube1_infile+' not found.', /info
                    continue
                 endif
                 cube1 = readfits(cube1_infile, cube1_hdr)

                 pb1_infile = $
                    release_dir+'raw/'+ $
                    gals[gal_ind[0]]+'_'+ $
                    this_array+'_'+this_product+'_'+$
                    'pb'+'.fits'

                 test = file_search(pb1_infile, count=found)
                 if found eq 0 then begin
                    message, 'File '+pb1_infile+' not found.', /info
                    continue
                 endif
                 pb_cube1 = readfits(pb1_infile, pb1_hdr)                 

                 cube2_infile = $
                    release_dir+'process/'+ $
                    gals[gal_ind[1]]+'_'+ $
                    this_array+'_'+this_product+'_'+$
                    this_ext+'.fits'

                 test = file_search(cube2_infile, count=found)
                 if found eq 0 then begin
                    message, 'File '+cube2_infile+' not found.', /info
                    continue
                 endif
                 cube2 = readfits(cube2_infile, cube2_hdr)

                 pb2_infile = $
                    release_dir+'raw/'+ $
                    gals[gal_ind[1]]+'_'+ $
                    this_array+'_'+this_product+'_'+$
                    'pb'+'.fits'

                 test = file_search(pb2_infile, count=found)
                 if found eq 0 then begin
                    message, 'File '+pb2_infile+' not found.', /info
                    continue
                 endif
                 pb_cube2 = readfits(pb2_infile, pb2_hdr)

                 if n_part eq 3 then begin
                    cube3_infile = $
                       release_dir+'process/'+ $
                       gals[gal_ind[2]]+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'.fits'
                    cube3 = readfits(cube3_infile, cube3_hdr)

                    test = file_search(cube3_infile, count=found)
                    if found eq 0 then begin
                       message, 'File '+cube3_infile+' not found.', /info
                       continue
                    endif
                    
                    pb3_infile = $
                       release_dir+'raw/'+ $
                       gals[gal_ind[2]]+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       'pb'+'.fits'

                    test = file_search(pb3_file, count=found)
                    if found eq 0 then begin
                       message, 'File '+pb3_file+' not found.', /info
                       continue
                    endif
                    pb_cube3 = readfits(pb3_infile, pb3_hdr)
                    
                 endif

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; FIGURE THE TARGET BEAM AND CONVOLVE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 beam1 = sxpar(cube1_hdr, 'BMAJ')
                 beam2 = sxpar(cube2_hdr, 'BMAJ')              
                 if n_part eq 2 then begin
                    target_bmaj = max([beam1, beam2],/nan)
                 endif else begin
                    beam3 = sxpar(cube3_hdr, 'BMAJ')
                    target_bmaj = max([beam1, beam2, beam3],/nan)
                 endelse
                 
                 conv_with_gauss $
                    , data=cube1 $
                    , hdr=cube1_hdr $
                    , target =[1,1,0.]*target_bmaj $
                    , out_data=smooth_cube1 $
                    , out_hdr=smooth_cube1_hdr $
                    , /perbeam $
                    , worked=worked
                 if worked eq 0 then begin
                    message, 'Problem with convolution in merging.', /info
                    stop
                 endif
                 
                 nan_ind = where(finite(cube1) eq 0, nan_ct)
                 if nan_ct gt 0 then smooth_cube1[nan_ind] = !values.f_nan
                 cube1 = smooth_cube1
                 cube1_hdr = smooth_cube1_hdr

                 conv_with_gauss $
                    , data=cube2 $
                    , hdr=cube2_hdr $
                    , target =[1,1,0.]*target_bmaj $
                    , out_data=smooth_cube2 $
                    , out_hdr=smooth_cube2_hdr $
                    , /perbeam $
                    , worked=worked
                 if worked eq 0 then begin
                    message, 'Problem with convolution in merging.', /info
                    stop
                 endif
                 
                 nan_ind = where(finite(cube2) eq 0, nan_ct)
                 if nan_ct gt 0 then smooth_cube2[nan_ind] = !values.f_nan
                 cube2 = smooth_cube2
                 cube2_hdr = smooth_cube2_hdr

                 if n_part eq 3 then begin

                    conv_with_gauss $
                       , data=cube3 $
                       , hdr=cube3_hdr $
                       , target =[1,1,0.]*target_bmaj $
                       , out_data=smooth_cube3 $
                       , out_hdr=smooth_cube3_hdr $
                       , /perbeam $
                       , worked=worked
                    if worked eq 0 then begin
                       message, 'Problem with convolution in merging.', /info
                       stop
                    endif
                    
                    nan_ind = where(finite(cube3) eq 0, nan_ct)
                    if nan_ct gt 0 then smooth_cube3[nan_ind] = !values.f_nan
                    cube3 = smooth_cube3
                    cube3_hdr = smooth_cube3_hdr

                 endif

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; BUILD A TARGET HEADER
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 target_hdr = cube1_hdr

                 cdelt = abs(sxpar(cube1_hdr,'CDELT1'))

                 npix_ra = abs(ceil(merge_dra[ii]/3600. / cdelt))
                 crpix_ra = npix_ra/2.0

                 npix_dec = abs(ceil(merge_ddec[ii]/3600. / cdelt))
                 crpix_dec = npix_dec/2.0

                 sxaddpar, target_hdr, 'CTYPE1', 'RA---SIN'
                 sxdelpar, target_hdr, 'CRVAL1'
                 sxaddpar, target_hdr, 'CRVAL1', double((merge_ra[ii]*1.0)[0])
                 sxdelpar, target_hdr, 'NAXIS1'              
                 sxaddpar, target_hdr, 'NAXIS1', long(npix_ra[0]), after='NAXIS'
                 sxdelpar, target_hdr, 'CRPIX1'
                 sxaddpar, target_hdr, 'CRPIX1', crpix_ra[0]*1.0

                 sxaddpar, target_hdr, 'CTYPE2', 'DEC--SIN'  
                 sxdelpar, target_hdr, 'CRVAL2'
                 sxaddpar, target_hdr, 'CRVAL2', double((merge_dec[ii]*1.0)[0])
                 sxdelpar, target_hdr, 'NAXIS2'
                 sxaddpar, target_hdr, 'NAXIS2', long(npix_dec[0]), after='NAXIS1'
                 sxdelpar, target_hdr, 'CRPIX2'
                 sxaddpar, target_hdr, 'CRPIX2', crpix_dec[0]*1.0

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; ALIGN TO THE TARGET HEADER
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 cube_hastrom $
                    , data = cube1 $
                    , hdr_in = cube1_hdr $
                    , target_hdr = target_hdr $
                    , outcube = aligned_cube1 $
                    , outhdr = aligned_cube1_hdr $
                    , missing=!values.f_nan
                 cube1 = aligned_cube1
                 cube1_hdr = aligned_cube1_hdr

                 cube_hastrom $
                    , data = pb_cube1 $
                    , hdr_in = pb1_hdr $
                    , target_hdr = target_hdr $
                    , outcube = aligned_pb_cube1 $
                    , outhdr = aligned_pb1_hdr $
                    , missing=!values.f_nan
                 pb_cube1 = aligned_pb_cube1
                 pb1_hdr = aligned_pb1_hdr
                 
                 cube_hastrom $
                    , data = cube2 $
                    , hdr_in = cube2_hdr $
                    , target_hdr = target_hdr $
                    , outcube = aligned_cube2 $
                    , outhdr = aligned_cube2_hdr $
                    , missing=!values.f_nan
                 cube2 = aligned_cube2
                 cube2_hdr = aligned_cube2_hdr

                 cube_hastrom $
                    , data = pb_cube2 $
                    , hdr_in = pb2_hdr $
                    , target_hdr = target_hdr $
                    , outcube = aligned_pb_cube2 $
                    , outhdr = aligned_pb2_hdr $
                    , missing=!values.f_nan
                 pb_cube2 = aligned_pb_cube2
                 pb2_hdr = aligned_pb2_hdr

                 if n_part eq 3 then begin

                    cube_hastrom $
                       , data = cube3 $
                       , hdr_in = cube3_hdr $
                       , target_hdr = target_hdr $
                       , outcube = aligned_cube3 $
                       , outhdr = aligned_cube3_hdr $
                       , missing=!values.f_nan
                    cube3 = aligned_cube3
                    cube3_hdr = aligned_cube3_hdr

                    cube_hastrom $
                       , data = pb_cube3 $
                       , hdr_in = pb3_hdr $
                       , target_hdr = target_hdr $
                       , outcube = aligned_pb_cube3 $
                       , outhdr = aligned_pb3_hdr $
                       , missing=!values.f_nan
                    pb_cube3 = aligned_pb_cube3
                    pb3_hdr = aligned_pb3_hdr

                 endif

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; COMBINE THE FILES INTO A SINGLE DATA PRODUCT
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 cube_out = cube1*!values.f_nan
                 
                 cov_part1 = finite(cube1) and pb_cube1 gt 0.
                 cov_part2 = finite(cube2) and pb_cube2 gt 0.
                 if n_part eq 3 then $
                    cov_part3 = finite(cube3) and pb_cube3 gt 0.
                 
                 if n_part eq 2 then begin
                    ind1 = where(cov_part1 eq 1 and cov_part2 eq 0, ct1)
                    if ct1 gt 0 then $
                       cube_out[ind1] = cube1[ind1]
                    
                    ind2 = where(cov_part1 eq 0 and cov_part2 eq 1, ct2)
                    if ct2 gt 0 then $
                       cube_out[ind2] = cube2[ind2]
                    
                    ind12 = where(cov_part1 eq 1 and cov_part2 eq 1, ct12)
                    if ct12 gt 0 then $
                       cube_out[ind12] = $
                       (cube1[ind12]*pb_cube1[ind12]^2 + $
                        cube2[ind12]*pb_cube2[ind12]^2) / $
                       (pb_cube1[ind12]^2 + pb_cube2[ind12]^2)
                 endif else begin

                    ind1 = where(cov_part1 eq 1 and $
                                 cov_part2 eq 0 and cov_part3 eq 0, ct1)
                    if ct1 gt 0 then $
                       cube_out[ind1] = cube1[ind1]
                    
                    ind2 = where(cov_part1 eq 0 and $
                                 cov_part2 eq 1 and cov_part3 eq 0, ct2)
                    if ct2 gt 0 then $
                       cube_out[ind2] = cube2[ind2]
                    
                    ind3 = where(cov_part1 eq 0 and $
                                 cov_part2 eq 0 and cov_part3 eq 1, ct3)
                    if ct3 gt 0 then $
                       cube_out[ind3] = cube3[ind3]
                    
                    ind12 = where(cov_part1 eq 1 and cov_part2 eq 1 and $
                                  cov_part3 eq 0, ct12)
                    if ct12 gt 0 then $
                       cube_out[ind12] = $
                       (cube1[ind12]*pb_cube1[ind12]^2 + $
                        cube2[ind12]*pb_cube2[ind12]^2) / $
                       (pb_cube1[ind12]^2 + pb_cube2[ind12]^2)

                    ind23 = where(cov_part1 eq 0 and cov_part2 eq 1 and $
                                  cov_part3 eq 1, ct23)
                    if ct23 gt 0 then $
                       cube_out[ind23] = $
                       (cube2[ind23]*pb_cube2[ind23]^2 + $
                        cube3[ind23]*pb_cube3[ind23]^2) / $
                       (pb_cube2[ind23]^2 + pb_cube3[ind23]^2)

                    ind13 = where(cov_part1 eq 1 and cov_part2 eq 0 and $
                                  cov_part3 eq 1, ct13)
                    if ct13 gt 0 then $
                       cube_out[ind13] = $
                       (cube1[ind13]*pb_cube1[ind13]^2 + $
                        cube3[ind13]*pb_cube3[ind13]^2) / $
                       (pb_cube1[ind13]^2 + pb_cube3[ind13]^2)

                    ind123 = where(cov_part1 eq 1 and cov_part2 eq 1 and $
                                   cov_part3 eq 1, ct123)
                    if ct123 gt 0 then $
                       cube_out[ind123] = $
                       (cube1[ind123]*pb_cube1[ind123]^2 + $
                        cube2[ind123]*pb_cube2[ind123]^2 + $
                        cube3[ind123]*pb_cube3[ind123]^2) / $
                       (pb_cube1[ind123]^2 + pb_cube2[ind123]^2 + pb_cube3[ind123]^2)
                    
                 endelse

                 !p.multi=[0,2,2]
                 disp, max(cube_out, dim=3, /nan), max=1
                 disp, max(cube1, dim=3, /nan), max=1
                 disp, max(cube2, dim=3, /nan), max=1
                 if n_part eq 3 then $
                    disp, max(cube3, dim=3, /nan), max=1

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; WRITE TO DISK
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 writefits $
                    , release_dir+'process/'+ $
                    merge_name[ii]+'_'+this_array+'_'+this_product+ $
                    '_'+this_ext+'.fits' $
                    , cube_out, target_hdr
                 
              endfor
              
           endfor

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; MERGE THE TOTAL POWER (USUALLY JUST A COPYING STEP)
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

           tp1_infile = $
              release_dir+'process/'+ $
              gals[gal_ind[0]]+'_tp_'+ $
              this_product+'_k.fits'

           test = file_search(tp1_infile, count=found)
           if found eq 0 then begin
              message, 'File '+tp1_infile+' not found.', /info
              continue
           endif

           tp1_cube = $
              readfits(tp1_infile, tp1_hdr)

           tp2_infile = $
              release_dir+'process/'+ $
              gals[gal_ind[1]]+'_tp_'+ $
              this_product+'_k.fits'

           test = file_search(tp2_infile, count=found)
           if found eq 0 then begin
              message, 'File '+tp2_infile+' not found.', /info
              continue
           endif

           tp2_cube = $
              readfits(tp2_infile, tp2_hdr)
           
           if n_part eq 3 then begin

              tp3_infile = $
                 release_dir+'process/'+ $
                 gals[gal_ind[2]]+'_tp_'+ $
                 this_product+'_k.fits'
              
              test = file_search(tp3_infile, count=found)
              if found eq 0 then begin
                 message, 'File '+tp3_infile+' not found.', /info
                 continue
              endif

              tp3_cube = $
                 readfits(tp3_infile, tp3_hdr)
              
           endif

;          In this case, just copy the TP from one case because it's
;          the same for all cases. This is the most common case.

           if merge_copy_tp[ii] then begin
              writefits, out_dir+merge_name[ii]+'_tp_'+this_product+'_k.fits' $
                         , tp1_cube, part1_hdr
           endif else begin

;          In this less common case align the different TP data sets
;          on to a new grid and average them where they overlap. This
;          generally shouldn't yield much overlap. Those cases appear
;          mostly as single cubes.
              
              target_hdr = tp1_hdr
              cdelt = abs(sxpar(target_hdr,'CDELT1'))

              npix_ra = abs(ceil(merge_dra[ii]/3600. / cdelt))
              crpix_ra = npix_ra/2.0
              
              npix_dec = abs(ceil(merge_ddec[ii]/3600. / cdelt))
              crpix_dec = npix_dec/2.0

              sxaddpar, target_hdr, 'CTYPE1', 'RA---SIN'
              sxdelpar, target_hdr, 'CRVAL1'
              sxaddpar, target_hdr, 'CRVAL1', double((merge_ra[ii]*1.0)[0])
              sxdelpar, target_hdr, 'NAXIS1'              
              sxaddpar, target_hdr, 'NAXIS1', long(npix_ra[0]), after='NAXIS'
              sxdelpar, target_hdr, 'CRPIX1'
              sxaddpar, target_hdr, 'CRPIX1', crpix_ra[0]*1.0

              sxaddpar, target_hdr, 'CTYPE2', 'DEC--SIN'  
              sxdelpar, target_hdr, 'CRVAL2'
              sxaddpar, target_hdr, 'CRVAL2', double((merge_dec[ii]*1.0)[0])
              sxdelpar, target_hdr, 'NAXIS2'
              sxaddpar, target_hdr, 'NAXIS2', long(npix_dec[0]), after='NAXIS1'
              sxdelpar, target_hdr, 'CRPIX2'
              sxaddpar, target_hdr, 'CRPIX2', crpix_dec[0]*1.0
              
              cube_hastrom $
                 , data = tp1_cube $
                 , hdr_in = tp1_hdr $
                 , target_hdr = target_hdr $
                 , outcube = new_tp1 $
                 , outhdr = new_tp1_hdr $
                 , missing=!values.f_nan
              
              cube_hastrom $
                 , data = tp2_cube $
                 , hdr_in = tp2_hdr $
                 , target_hdr = target_hdr $
                 , outcube = new_tp2 $
                 , outhdr = new_tp2_hdr $
                 , missing=!values.f_nan
              
              if n_part eq 3 then begin

                 cube_hastrom $
                    , data = tp3_cube $
                    , hdr_in = tp3_hdr $
                    , target_hdr = target_hdr $
                    , outcube = new_tp3 $
                    , outhdr = new_tp3_hdr $
                    , missing=!values.f_nan
                 
              endif

              cube_out = new_tp1*!values.f_nan
              
              cov_tp1 = finite(new_tp1)
              cov_tp2 = finite(new_tp2)
              if n_part eq 3 then cov_tp3 = finite(new_tp3)

              if n_part eq 2 then begin

                 ind1 = where(cov_tp1 eq 1 and cov_tp2 eq 0, ct1)
                 if ct gt 0 then $
                    cube_out[ind1] = new_tp1[ind1]
                 
                 ind2 = where(cov_tp1 eq 0 and cov_tp2 eq 1, ct2)
                 if ct2 gt 0 then $
                    cube_out[ind2] = new_tp2[ind2]
                 
                 ind12 = where(cov_tp1 eq 1 and cov_tp2 eq 1, ct12)
                 if ct12 gt 0 then $
                    cube_out[ind12] = (new_tp1[ind12] + new_tp2[ind12])/2.0

              endif else begin

                 ind1 = where(cov_tp1 eq 1 and cov_tp2 eq 0 and $
                              cov_tp3 eq 0, ct1)
                 if ct gt 0 then $
                    cube_out[ind1] = new_tp1[ind1]
                 
                 ind2 = where(cov_tp1 eq 0 and cov_tp2 eq 1 and $
                              cov_tp3 eq 0, ct2)
                 if ct2 gt 0 then $
                    cube_out[ind2] = new_tp2[ind2]

                 ind3 = where(cov_tp1 eq 0 and cov_tp2 eq 0 and $
                              cov_tp3 eq 1, ct3)
                 if ct3 gt 0 then $
                    cube_out[ind3] = new_tp2[ind3]
                 
                 ind12 = where(cov_tp1 eq 1 and cov_tp2 eq 1 and $
                               cov_tp3 eq 0, ct12)
                 if ct12 gt 0 then $
                    cube_out[ind12] = (new_tp1[ind12] + new_tp2[ind12])/2.0

                 ind23 = where(cov_tp1 eq 0 and cov_tp2 eq 1 and $
                               cov_tp3 eq 1, ct23)
                 if ct23 gt 0 then $
                    cube_out[ind23] = (new_tp2[ind23] + new_tp3[ind23])/2.0
                 
                 ind13 = where(cov_tp1 eq 1 and cov_tp2 eq 0 and $
                               cov_tp3 eq 1, ct13)
                 if ct13 gt 0 then $
                    cube_out[ind13] = (new_tp1[ind13] + new_tp3[ind13])/2.0

                 ind123 = where(cov_tp1 eq 1 and cov_tp2 eq 1 and $
                                cov_tp3 eq 1, ct123)
                 if ct123 gt 0 then $
                    cube_out[ind123] = (new_tp1[ind123] + new_tp2[ind123] + $
                                        new_tp3[ind123])/3.0
                 
              endelse
              
              writefits, release_dir+'process/'+ $
                         merge_name[ii]+'_tp_'+this_product+'_k.fits' $
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

     ext_to_process = $
        ['flat_round' $
         , 'pbcorr_round']
     n_ext = n_elements(ext_to_process)

     for ii = 0, n_fullgals-1 do begin
        
        if n_elements(only) gt 0 then $
           if total(only eq fullgal_list[ii]) eq 0 then continue

        if n_elements(skip) gt 0 then $
           if total(skip eq fullgal_list[ii]) gt 0 then continue

        this_gal = fullgal_list[ii]

        message, "Sanitizing cubes for "+this_gal, /info

        for jj = 0, n_fullarray - 1 do begin

           this_array = fullarray_list[jj]
           if n_elements(just_array) gt 0 then $
              if total(just_array eq this_array) eq 0 then continue

           for kk = 0, n_product-1 do begin
              
              this_product = product_list[kk]

              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_process[ll]

                 in_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+this_product+ $
                    '_'+this_ext+'.fits'

                 test = file_search(in_file, count=found)
                 if found eq 0 then begin
                    message, 'File '+in_file+' not found.', /info
                    continue
                 endif

                 out_file = $
                    release_dir+'process/'+ $
                    this_gal+'_'+this_array+'_'+this_product+ $
                    '_'+this_ext+'_trimmed.fits'
                 
                 message, '... cube '+in_file, /info
                 cube_trim $
                    , data = in_file $
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
                    , release_dir+'process/'+$
                    this_gal+'_'+this_array+'_'+this_product+ $
                    '_'+this_ext+'_k.fits' $
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
     message, 'CONVOLVING TO FIXED PHYSICAL RESOLUTION', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info     
     
     dir = release_dir+'process/'
     tol = 0.1

     s = gal_data(gal_for_fullgals)

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

                 sxaddpar, hdr, 'DIST', s[ii].dist_mpc, 'MPC / USED IN CONVOLUTION'
                 current_res_pc = s[ii].dist_mpc*!dtor*sxpar(hdr, 'BMAJ')*1d6
                 
                 for zz = 0, n_res -1 do begin
                    
                    res_str = strcompress(str(target_res[zz]),/rem)+'pc'
                    out_file = release_dir+'process/'+ $
                               this_gal+'_'+this_array+'_'+ $
                               this_product+'_'+this_ext+'_'+res_str+'.fits'
                    target_res_as = target_res[zz]/(s[ii].dist_mpc*1d6)/!dtor*3600.d
                    
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

end
