pro plot_integrated_spectrum $
   , pause=pause
  
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

; SKIP

  only = []
  skip = []

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; LOOP OVER GALAXIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  for ii = 0, n_gals-1 do begin

     if n_elements(skip) gt 0 then begin
        if total(gals[ii] eq skip) gt 0 then begin
           print, "Skipping "+gals[ii]
           continue
        endif
     endif

     if n_elements(only) gt 0 then begin
        if total(gals[ii] eq only) eq 0 then begin
           print, "Skipping "+gals[ii]
           continue
        endif
     endif

     message, "Plotting integrated spectrum for "+gals[ii], /info

     dir = release_dir+'process/'
     
     if n_elements(only) gt 0 then $
        if total(only eq gals[ii]) eq 0 then continue

     fname = dir+gals[ii]+'_7m_co21_flat_round_k.fits'
     if file_test(fname) eq 0 then begin
        print, "No file found for "+gals[ii]
        print, "Continuing."
        continue
     endif
        
     cube = readfits(fname, hdr)
     pk_map = max(cube, dim=3, /nan)
     dummy = max(pk_map, maxind, /nan)
     sz = size(pk_map)
     ind_to_xyv, maxind, x=x, y=y, sz=sz
     max_spec = cube[x, y, *]

     make_axes, hdr, vaxis=vaxis
     vaxis /= 1d3

     spec = total(total(cube,1,/nan),1,/nan) / $
            total(total(finite(cube)*1.0,1),1)

     plot, vaxis, spec, ps=10, title=gals[ii]
     oplot, vaxis, max_spec*max(spec,/nan)/max(max_spec,/nan) $
            , ps=10, color=cgcolor('red')
     
     if keyword_set(pause) then begin
        print, "Clean mask for "+gals[ii]+". Key to continue."
        test = get_kbrd(1)
     endif

  endfor


end
