pro coadd_ngc6744 

  ext = ['_co21_round' $
         , '_co21_round_pbcor' $
         , '_co21_residual' $
         , '_co21_pb' $
        ]
  n_ext = n_elements(ext)

  for ii = 0, n_ext-1 do begin
     
     north = readfits('../ngc6744/ngc6744north'+ext[ii]+'_align.fits' $
                      , north_hdr)

     south = readfits('../ngc6744/ngc6744south'+ext[ii]+'_align.fits' $
                      , south_hdr)
    
     combo = north
     combo[*,0:1500,*] = south[*,0:1500,*]
     hdr = north_hdr
     
     sxdelpar, hdr, 'BLANK'
     writefits, '../ngc6744/ngc6744'+ext[ii]+'.fits' $
                , cube, hdr
     
     loadct, 33
     disp, max(combo, dim=3, /nan), /sq
     
  endfor

end
