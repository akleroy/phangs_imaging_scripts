pro build_release_v2 $
   , wipe=do_wipe $
   , copy=do_copy

;+
;
;Assemble the cubes and maps made by the build_cubes and
;build_produts programs into a release for the team. Mostly just
;copying between directories.
;
;-

  out_dir = '../release/v2/delivery/'
  in_dir = '../release/v2/process/'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WIPE THE DELIVERY DIRECTORY AND RESET WITH NEW DIRECTORIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_wipe) then begin

     spawn, 'rm -rf '+out_dir+'*'

     spawn, 'mkdir '+out_dir+'strict_maps/'
     spawn, 'mkdir '+out_dir+'broad_maps/'
     spawn, 'mkdir '+out_dir+'cubes/'     

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FILES TO THE APPROPRIATE SUBDIRECTORIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_copy) then begin

     spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_strict_*.fits '+out_dir+'strict_maps/.'
     spawn, 'cp '+in_dir+'*_12m+7m_co21_strict_*.fits '+out_dir+'strict_maps/.'

     spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_broad_*.fits '+out_dir+'broad_maps/.'
     spawn, 'cp '+in_dir+'*_12m+7m_co21_broad_*.fits '+out_dir+'broad_maps/.'

     spawn, 'cp '+in_dir+'*_12m+7m_co21_pbcorr_round_k.fits '+out_dir+'cubes/.'
     spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_pbcorr_round_k.fits '+out_dir+'cubes/.'

     spawn, 'cp README_for_v2.txt '+out_dir+'.'

  endif

; Go tar things up.  

end
