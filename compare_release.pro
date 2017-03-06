pro compare_release

  gals = $
     [ $
     'ngc0628' $
     , 'ngc1672' $
     , 'ngc3351' $
     , 'ngc3627' $
     , 'ngc4254' $
     , 'ngc4303' $
     , 'ngc4321' $
     , 'ngc4535' $
     , 'ngc5068' $
     ]
  n_gal = n_elements(gals)

  dir_1 = '../release/v0p2/delivery/'
  dir_2 = '../release/v0p4/delivery/'

  for ii = 0, n_gal-1 do begin
     
     loadct, 0

     m = readfits(dir_2+gals[ii]+'_co21_smoothedmask.fits', hm)

     v1 = readfits(dir_1+strupcase(gals[ii])+'_co21_cube.fits', h1)
     cube_hastrom, data=m, hdr=hm, target_hdr=h1, outcube=m1
     m1 = m1 gt 0.5
     v1 = total(total(v1*m1,1,/nan),1)*abs((sxpar(h1,'CDELT2'))^2*(sxpar(h1,'CDELT3')))
     make_axes, h1, vaxis=vaxis_1, /vonly

     v2 = readfits(dir_2+gals[ii]+'_co21.fits', h2)
     v2 = total(total(v2*m,1,/nan),1)*abs((sxpar(h2,'CDELT2'))^2*(sxpar(h2,'CDELT3')))
     make_axes, h2, vaxis=vaxis_2, /vonly

     v3 = readfits(dir_1+strupcase(gals[ii])+'_co21_taper.fits', h3)
     cube_hastrom, data=m, hdr=hm, target_hdr=h3, outcube=m3
     v3 = total(total(v3*m3,1,/nan),1)*abs((sxpar(h3,'CDELT2'))^2*(sxpar(h3,'CDELT3')))
     make_axes, h3, vaxis=vaxis_3, /vonly

     loadct, 0
     plot, vaxis_1, v1, title='!6'+gals[ii], yrange=[-0.1,1.1]*max([v1,v2,v3],/nan)
     oplot, vaxis_1, v1, color=cgcolor('blue')
     oplot, vaxis_2, v2, lines=2, color=cgcolor('green')
     oplot, vaxis_3, v3, color=cgcolor('red')

     stop
     continue

     v1 = readfits(dir_1+strupcase(gals[ii])+'_co21_cube_mom0.fits', h1)
     v2 = readfits(dir_2+gals[ii]+'_co21_mom0.fits', h2)

;     v1 = readfits(dir_1+strupcase(gals[ii])+'_co21_cube_tpeak_12p5kms.fits', h1)
;     v2 = readfits(dir_2+gals[ii]+'_co21_tpeak_12p5kms.fits', h2)

     hastrom, v2, h2, h1, interp=2, cubic=-0.5, missing=!values.f_nan
     
     plot, v1, v2, ps=3, /iso
     equality
     print, gals[ii]+' integrated new/old:' , total(v2,/nan)/total(v1,/nan)
     print, gals[ii]+' median new/old: ',median(v2/v1)

     stop

  endfor

end
