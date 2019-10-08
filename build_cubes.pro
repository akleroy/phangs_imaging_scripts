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
   , prepmerge=do_mergeprep $
   , stagemerge=do_stagemerge $
   , merge=do_merge $
   , sanitize=do_sanitize

;+
;
;Scripts to build the imaged data into data cubes suitable for
;analysis. Steps are: copy and primary beam correct, convolve to a
;round beam, align single dish to 
;
;-

; DIRECTORIES

  root_imaging_dir = '../'

  version_tag = '3.4'

  if n_elements(version) eq 0 then $
     version = '3'
  
  if version eq '3' then begin
     vstring = 'v3'
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
  valid_products = ['co21','c18o21','13co21']
  if n_elements(product_list) eq 0 then begin
     product_list = ['co21']
  endif
  n_product = n_elements(product_list)
  for ii = 0, n_product-1 do begin
     this_product = product_list[ii]
     if total(this_product eq valid_products) eq 0 then begin
        print, "Illegal product in product list: ", this_product
        print, "Stopping."
        stop
     endif
  endfor

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

     spawn, 'rm -rf '+release_dir+'products/'
     spawn, 'mkdir '+release_dir+'products/'

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
                 
                 spawn, '\cp '+ in_file +' '+ out_file
                 
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

        message, "Convolving to a round beam for "+this_gal, /info

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

        message, 'Copying feathered data for '+this_gal, /info

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
; HOMOGENIZE THE RESOLUTION OF MULTI-FIELD CUBES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
; Loop over all parts of each multi-part cube, read the pixel scale
; and beam size, and then convolve all parts to a common resolution
; for future linear mosaicking. Write the results to disk as an
; intermediate product.

  if keyword_set(do_mergeprep) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'PREPARING MULTI-PART MOSAICS FOR LINEAR MOSAICKING', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     message, '', /info
     message, 'Moving all parts of a galaxy to a common resolution.', /info
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

     tol_for_beam = 1.00

;    LOOP OVER MERGING GALAXIES
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
        
;       LOOP OVER DATA PRODUCTS
        for kk = 0, n_product-1 do begin
           
           this_product = product_list[kk]

;          LOOP OVER ARRAYS
           for jj = 0, n_fullarray - 1 do begin
              
              this_array = fullarray_list[jj]
              if n_elements(just_array) gt 0 then $
                 if total(just_array eq this_array) eq 0 then continue           
              
;             LOOP OVER EXTENSIONS 
              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_merge[ll]

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; LOOP OVER PARTS FOR THIS GALAXY TO GET A COMMON RESOLUTION
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 found = 0B

                 for mm = 0, n_part-1 do begin
                    
                    this_gal_part = gals[gal_ind[mm]]
                    
                    this_infile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'.fits'

                    if file_test(this_infile) eq 0 then begin
                       continue
                    endif
                    if found eq 0B then begin
                       first = 1B
                       found = 1B
                    endif else begin
                       first = 0B
                    endelse

                    this_hdr = headfits(this_infile)
                    this_bmaj = sxpar(this_hdr,'BMAJ')*3600.
                    this_pix = abs(sxpar(this_hdr,'CDELT1'))*3600.
                    print, this_infile, " has beam ", this_bmaj, " and pixel ", this_pix
                    if first then begin
                       bmaj_list = [this_bmaj] 
                       pix_list = [this_pix] 
                       hdr_list = [this_hdr]
                    endif else begin
                       bmaj_list = [bmaj_list, this_bmaj]
                       pix_list = [pix_list, this_pix]
                       hdr_list = [hdr_list, this_hdr]
                    endelse
                 endfor

                 if found eq 0B then begin
                    print, "No data found for this extension plus products. Continuing."
                    continue
                 endif

                 max_bmaj = max(bmaj_list, /nan)
                 max_pix = max(pix_list, /nan)

                 print, "Maximum beam size is: ", max_bmaj
                 target_bmaj = sqrt(max_bmaj^2+(2.0*max_pix)^2)
                 print, "I will target: ", target_bmaj

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; CONVOLVE THIS COMMON RESOLUTION
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 for mm = 0, n_part-1 do begin
                    
                    this_gal_part = gals[gal_ind[mm]]
                    
                    this_infile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'.fits'

                    if file_test(this_infile) eq 0 then begin
                       continue
                    endif

                    this_outfile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'_tomerge.fits'

                    this_hdr = headfits(this_infile)
                    this_bmaj = sxpar(this_hdr,'BMAJ')*3600.
                    if this_bmaj*tol_for_beam ge target_bmaj then begin
                       print, "... resolution within specified tolerance of "+str(tol_for_beam)+". I will copy the data."
                       spawn, 'rm -rf '+this_outfile
                       spawn, 'cp '+this_infile+' '+this_outfile
                    endif else begin                   

                       print, "... resolution not within specified tolerance. I will convolve the data."
                       conv_with_gauss $
                          , data=this_infile $
                          , target =[1,1,0.]*target_bmaj $
                          , out_file=this_outfile $
                          , /perbeam $
                          , worked=worked
                    endelse

                 endfor

              endfor

           endfor

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; BUILD A COMMON HEADER AND REPROJECT THE CUBES ONTO IT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Loop over all parts of each multi-part cube, build a target header
; with specifications set by a user input file, and then reproject
; both the cubes and the associated primary beam files onto this
; astrometry. The primary beam will be used as a weight.

  if keyword_set(do_stagemerge) then begin

     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'STAGING MULTI-PART MOSAICS FOR LINEAR MOSAICKING', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     
     message, '', /info
     message, 'Building a target header and reprojecting onto it.', /info
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

;    LOOP OVER MERGING GALAXIES
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
        
;       LOOP OVER DATA PRODUCTS
        for kk = 0, n_product-1 do begin
           
           this_product = product_list[kk]

;          LOOP OVER ARRAYS
           for jj = 0, n_fullarray - 1 do begin
              
              this_array = fullarray_list[jj]
              if n_elements(just_array) gt 0 then $
                 if total(just_array eq this_array) eq 0 then continue           

              array_for_pb = this_array
              if array_for_pb eq '7m+tp' then array_for_pb = '7m'
              if array_for_pb eq '12m+7m+tp' then array_for_pb = '12m+7m'

;             LOOP OVER EXTENSIONS 
              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_merge[ll]

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; LOOP OVER PARTS FOR THIS GALAXY TO BUILD A TARGET HEADER
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
                 
                 found = 0B

                 for mm = 0, n_part-1 do begin
                    
                    this_gal_part = gals[gal_ind[mm]]
                    
                    this_infile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'_tomerge.fits'

                    if file_test(this_infile) eq 0 then begin
                       continue
                    endif
                    if found eq 0B then begin
                       first = 1B
                       found = 1B
                    endif else begin
                       first = 0B
                    endelse

                    this_hdr = headfits(this_infile)
                    this_bmaj = sxpar(this_hdr,'BMAJ')*3600.
                    this_pix = abs(sxpar(this_hdr,'CDELT1'))*3600.
                    print, this_infile, " has beam ", this_bmaj, " and pixel ", this_pix
                    if first then begin
                       bmaj_list = [this_bmaj] 
                       pix_list = [this_pix] 
                       hdr_list = [this_hdr]
                       first_hdr = this_hdr
                    endif else begin
                       bmaj_list = [bmaj_list, this_bmaj]
                       pix_list = [pix_list, this_pix]
                       hdr_list = [hdr_list, this_hdr]
                    endelse
                 endfor

                 if found eq 0B then begin
                    print, "No data found for this extension plus products. Continuing."
                    continue
                 endif

                 print, "Beams (should all be the same): ", bmaj_list

                 target_hdr = first_hdr
                 
                 cdelt = abs(sxpar(first_hdr,'CDELT1'))
                 
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

; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; REPROJECT CUBE AND PRIMARY BEAM ONTO THIS HEADER
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 for mm = 0, n_part-1 do begin
                    
                    this_gal_part = gals[gal_ind[mm]]
                    
; REPROJECT THE RAW IMAGE FROM RAW DIRECTORY ONTO THE GRID

                    this_infile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'_tomerge.fits'

                    if file_test(this_infile) eq 0 then begin
                       continue
                    endif
                    
                    this_outfile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'_onmergegrid.fits'

                    cube_hastrom $
                       , data = this_infile $
                       , target_hdr = target_hdr $
                       , outfile = this_outfile $
                       , missing=!values.f_nan

; REPROJECT THE PRIMARY BEAM IMAGE FROM RAW DIRECTORY ONTO THE GRID

                    this_pbinfile = $
                       release_dir+'raw/'+ $
                       this_gal_part+'_'+ $
                       array_for_pb+'_'+this_product+'_'+$
                       'pb'+'.fits'

                    this_pboutfile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'_mergeweight.fits'

                    cube_hastrom $
                       , data = this_pbinfile $
                       , target_hdr = target_hdr $
                       , outfile = this_pboutfile $
                       , missing=!values.f_nan

                 endfor                 

              endfor

           endfor

        endfor

     endfor

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LINEARLY MOSAIC THE CUBES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
  
; Linearly mosaic all of the galaxy parts together for each
; project. Keep a running sum of weights and image values so that we
; can handle an arbitrary number of fields.

  if keyword_set(do_merge) then begin
     
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
     message, 'MERGING MULTI-PART MOSAICS VIA LINEAR MOSAICKING', /info
     message, '%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&', /info
    
     message, '', /info
     message, 'Linearly mosaicking together parts of multi-part cubes.', /info
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

;    LOOP OVER MERGING GALAXIES
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
        
;       LOOP OVER DATA PRODUCTS
        for kk = 0, n_product-1 do begin
           
           this_product = product_list[kk]

;          LOOP OVER ARRAYS
           for jj = 0, n_fullarray - 1 do begin
              
              this_array = fullarray_list[jj]
              if n_elements(just_array) gt 0 then $
                 if total(just_array eq this_array) eq 0 then continue           
              
              for ll = 0, n_ext-1 do begin

                 this_ext = ext_to_merge[ll]

                 found = 0B
                 first = 1B
                 
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; ACCUMULATE
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 for mm = 0, n_part-1 do begin
                    
                    this_gal_part = gals[gal_ind[mm]]
                    
                    this_infile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'_onmergegrid.fits'

                    this_pbinfile = $
                       release_dir+'process/'+ $
                       this_gal_part+'_'+ $
                       this_array+'_'+this_product+'_'+$
                       this_ext+'_mergeweight.fits'

                    if file_test(this_infile) eq 0 then begin
                       print, "... missing file "+this_infile
                       print, "... continuing."
                       continue
                    endif

                    if file_test(this_pbinfile) eq 0 then begin
                       print, "... missing file "+this_pbinfile
                       print, "... continuing."
                       continue
                    endif

                    this_cube = $
                       readfits(this_infile, this_hdr)
                    this_weight= $
                       readfits(this_pbinfile, this_weight_hdr)

                    if first then begin
                       first_cube_hdr = this_hdr
                       sum_cube = finite(this_cube)*0.0
                       weight_cube = finite(this_cube)*0.0
                       first = 0B
                    endif                       
                    
                    mosaic_ind = where(finite(this_cube), mosaic_ct)
                    if mosaic_ct gt 0 then begin
                       sum_cube[mosaic_ind] = $
                          sum_cube[mosaic_ind] + $
                          this_cube[mosaic_ind]*this_weight[mosaic_ind]^2
                       weight_cube[mosaic_ind] = $
                          weight_cube[mosaic_ind] + $
                          this_weight[mosaic_ind]^2
                    endif

                 endfor
                 
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
; DIVIDE, CLEAN UP, AND WRITE TO DISK
; -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                 combined_cube = sum_cube / weight_cube
                 covered_ind = where(weight_cube eq 0., covered_ct)
                 if covered_ct gt 0 then $
                    combined_cube[covered_ind] = !values.f_nan

                 outfile = $
                    release_dir+'process/'+ $
                    this_merge_gal+'_'+ $
                    this_array+'_'+this_product+'_'+$
                    this_ext+'.fits'
                 writefits $
                    , outfile $
                    , combined_cube, first_cube_hdr
 
                writefits $
                    , release_dir+'process/'+ $
                    this_merge_gal+'_'+ $
                    this_array+'_'+this_product+'_'+$
                    this_ext+'_coverage.fits' $
                    , weight_cube, first_cube_hdr

                loadct, 33
                peak_map = max(combined_cube, dim=3, /nan)
                disp, peak_map $
                      , /sq, title=outfile, min=0, max=10.*mad(combined_cube) $
                      , reserve=5, color=cgcolor('white',255)

              endfor
              
           endfor
           
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
                 sxaddpar, hdr, 'JTOK', jtok, 'Janskies to Kelvin conversion'

                 sxaddpar, hdr, 'DATAMAX', max(cube,/nan)
                 sxaddpar, hdr, 'DATAMIN', min(cube,/nan)
                 sxaddpar, hdr, 'OBJECT', strupcase(this_gal)
                 sxaddpar, hdr, 'ORIGIN', 'PHANGS-ALMA-'+version_tag

                 sxdelpar, hdr, 'BLANK'
                 sxdelpar, hdr, 'DATE-OBS'
                 sxdelpar, hdr, 'OBSERVER'

                 sxdelpar, hdr, 'O_BLANK'
                 sxdelpar, hdr, 'O_BSCALE'
                 sxdelpar, hdr, 'O_BZERO'

                 sxdelpar, hdr, 'OBSRA'
                 sxdelpar, hdr, 'OBSDEC'
                 sxdelpar, hdr, 'OBSGEO-X'
                 sxdelpar, hdr, 'OBSGEO-Y'
                 sxdelpar, hdr, 'OBSGEO-Z'

                 sxdelpar, hdr, 'DISTANCE'

                 sxdelpar, hdr, 'HISTORY'
                 sxaddpar, hdr, 'HISTORY', 'This cube was produced by the PHANGS-ALMA pipeline.'
                 sxaddpar, hdr, 'HISTORY', 'This is part of data release '+version_tag

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

end
