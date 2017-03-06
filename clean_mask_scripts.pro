pro clean_mask_scripts
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; NGC 3351
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  print, "-=-=-=-=-=-=-=-=-=-=-=-=-=-="
  print, "NGC 3351"
  print, "-=-=-=-=-=-=-=-=-=-=-=-=-=-="

  mask1 = readfits('../clean_masks/ngc3351_co21_widemask.fits', mask_hdr)

  hera = readfits('../clean_masks/ngc3351_heracles_beta_at30.mask.fits', hera_hdr)
  reg = label_region(hera)
  sz = size(hera)
  hera = reg eq reg[sz[1]/2, sz[2]/2, sz[3]/2]

  cube_hastrom $
     , data = hera $
     , hdr_in = hera_hdr $
     , outcube = hera_on_alma $
     , target_hdr = mask_hdr
  
  combo_mask = grow_mask(mask1, constraint=(hera_on_alma ge 0.75))
  combo_mask = grow_mask(combo_mask, /xy, iters=5)
  writefits, '../clean_masks/ngc3351_clean_mask.fits', combo_mask, mask_hdr

  stop

end
