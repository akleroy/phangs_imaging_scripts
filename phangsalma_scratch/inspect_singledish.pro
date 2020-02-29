pro inspect_singledish

  readcol, 'singledish_key.txt', format='A,A', comment='#', name, fname
  nf = n_elements(fname)

  for ii = 0, nf-1 do begin

     cube = readfits(fname[ii], hdr)
     mapmax = max(cube, dim=3, /nan)
     mapmin = min(cube, dim=3, /nan)
     loadct, 33
     !p.multi=[0,2,1]
     disp, mapmax
     disp, mapmin
     !p.multi=0
     print, name[ii], fname[ii]
     print, ' hit a key to continue.'
     ch = get_kbrd(1)

  endfor

end
