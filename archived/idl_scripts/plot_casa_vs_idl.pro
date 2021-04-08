pro plot_casa_vs_idl $
   , pause=pause
;+
; 
; This is a diagnostic routine, probably not generally useful to
; anyone but me. The test could me migrated to be a general
; regression, but for now is modified by hand to allow me to benchmark
; the IDL vs CASA post-processing.
;
;-

  idl_dir = '../release/v3/process/'
  len_idl = strlen(idl_dir)
  casa_dir = '../release/v3_casa/process/'
  len_casa = strlen(casa_dir)

; Modify this search to get the list of files you want to compare.
  flist = file_search(casa_dir+'*_7m+tp_*_pbcorr_round_k.fits' $
                      , count=fct)
  
  for ii = 0, fct-1 do begin
     
     casa_file = flist[ii]

     idl_file = idl_dir+strmid(flist[ii],len_casa $
                               ,strlen(casa_file)-len_casa)

     if file_test(idl_file) eq 0 then begin
        if n_elements(no_match) eq 0 then $
           no_match = [casa_file] $
        else $
           no_match = [no_match, casa_file]
        continue
     endif

     print, 'Found: '
     print, '... casa file - ', casa_file
     print, '... idl file - ', idl_file
     
     casa_cube = readfits(casa_file, casa_hdr)
     casa_mom0 = total(casa_cube, 3, /nan)
     casa_peak = max(casa_cube, dim=3, /nan)
     casa_nan_ind = where(total(finite(casa_cube), 3) eq 0)
     casa_mom0[casa_nan_ind] = !values.f_nan
     casa_peak[casa_nan_ind] = !values.f_nan

     idl_cube = readfits(idl_file, idl_hdr)
     idl_mom0 = total(idl_cube, 3, /nan)
     idl_peak = max(idl_cube, dim=3, /nan)
     idl_nan_ind = where(total(finite(idl_cube), 3) eq 0)
     idl_mom0[idl_nan_ind] = !values.f_nan
     idl_peak[idl_nan_ind] = !values.f_nan

     hcopy1 = twod_head(casa_hdr)
     hastrom, casa_mom0, hcopy1, twod_head(idl_hdr) $
              , cubic=-0.5, interp=2, missing=!values.f_nan

     hcopy2 = twod_head(casa_hdr)
     hastrom, casa_peak, hcopy2, twod_head(idl_hdr) $
              , cubic=-0.5, interp=2, missing=!values.f_nan

     !p.multi=[0,2,1]

     loadct, 0
     plot, idl_mom0, casa_mom0, ps=3, title='Moment 0'
     equality, color=cgcolor('red')
     print, "Moment 0:"
     print, "... median ratio IDL/CASA: ", median(idl_mom0/casa_mom0)
     print, "... mad ratio IDL/CASA: ", mad(idl_mom0/casa_mom0)
     print, "... max difference IDL/CASA: ", max(abs(idl_mom0-casa_mom0),/nan)

     loadct, 0
     plot, idl_peak, casa_peak, ps=3, title='Peak Intensity'
     equality, color=cgcolor('red')
     print, "Peak:"
     print, "... median ratio IDL/CASA: ", median(idl_peak/casa_peak)
     print, "... mad ratio IDL/CASA: ", mad(idl_peak/casa_peak)
     print, "... max difference IDL/CASA: ", max(abs(idl_peak-casa_peak),/nan)

     !p.multi=0

     if n_elements(x) eq 0 then begin
        ind = where(finite(idl_mom0) and finite(casa_mom0), ct)
        if ct eq 0 then continue
        x = idl_mom0[ind]
        y = casa_mom0[ind]
        p = idl_peak[ind]
        q = casa_peak[ind]
     endif else begin
        ind = where(finite(idl_mom0) and finite(casa_mom0), ct)
        if ct eq 0 then continue
        x = [x, idl_mom0[ind]]
        y = [y, casa_mom0[ind]]
        p = [p, idl_peak[ind]]
        q = [q, casa_peak[ind]]
     endelse

     if keyword_set(pause) then begin
        print, "Hit a key to continue."
        ch = get_kbrd(1)
     endif

  endfor

; Make some general comparison plots.
  psfile='../plots/CASA_vs_IDL.eps'
  pnfile='../plots/CASA_vs_IDL.png'
  ps, /def, xs=8, ys=8, /encaps, /color, file=psfile, /ps
  plot, alog10(x), alog10(y), ps=3, xrange=[-1, 3], yrange=[-1,3] $
        , xtitle='log10 IDL Moment 0' $
        , ytitle='log10 CASA Moment 0'
  equality, color=cgcolor('red')
  ps, /xw
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

  psfile='../plots/CASA_vs_IDL_peak.eps'
  pnfile='../plots/CASA_vs_IDL_peak.png'
  ps, /def, xs=8, ys=8, /encaps, /color, file=psfile, /ps
  plot, alog10(p), alog10(q), ps=3, xrange=[-1, 3], yrange=[-1,3] $
        , xtitle='log10 IDL Peak' $
        , ytitle='log10 CASA Peak'
  equality, color=cgcolor('red')
  ps, /xw
  spawn, 'convert -density 300x300 '+psfile+' '+pnfile

  print, "Whole sample:"
  print, "... median CASA/IDL Mom0: ", median(y/x)
  print, "... mad CASA/IDL Mom0: ", mad(y/x)
  print, "... std CASA/IDL Mom0: ", stddev(y/x,/nan)
  print, "... madlog CASA/IDL Mom0: ", mad(alog10(y/x))
  print, "... stdlog CASA/IDL Mom0: ", stddev(alog10(y/x),/nan)

  print, "Whole sample:"
  print, "... median CASA/IDL Peak: ", median(q/p)
  print, "... mad CASA/IDL Peak: ", mad(q/p)
  print, "... std CASA/IDL Peak: ", stddev(q/p,/nan)
  print, "... madlog CASA/IDL Peak: ", mad(alog10(q/p))
  print, "... stdlog CASA/IDL Peak: ", stddev(alog10(q/p),/nan)

  stop

end
