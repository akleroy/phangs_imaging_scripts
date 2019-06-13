pro build_release_v3 $
   , wipe=do_wipe $
   , copy=do_copy $
   , array=just_array

;+
;
;Assemble the cubes and maps made by the build_cubes and
;build_produts programs into a release for the team. Mostly just
;copying between directories.
;
;-

  out_dir = '../release/v3/delivery/'
  in_dir = '../release/v3/process/'

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WIPE THE DELIVERY DIRECTORY AND RESET WITH NEW DIRECTORIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_wipe) then begin

     spawn, 'rm -rf '+out_dir+'*'

     spawn, 'mkdir '+out_dir+'strict_maps/'
     spawn, 'mkdir '+out_dir+'broad_maps/'
     spawn, 'mkdir '+out_dir+'cubes/'     
     spawn, 'mkdir '+out_dir+'support/'

  endif

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; COPY FILES TO THE APPROPRIATE SUBDIRECTORIES
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if keyword_set(do_copy) then begin

     if n_elements(just_array) gt 0 then begin
        
        if total(just_array eq '7m') gt 0 then begin
           print, "Copying maps."
           spawn, 'cp '+in_dir+'*_7m_co21_strict_*.fits '+out_dir+'strict_maps/.'
           spawn, 'cp '+in_dir+'*_7m_co21_broad_*.fits '+out_dir+'broad_maps/.'        
           print, "Copying cubes."
           spawn, 'cp '+in_dir+'*_7m_co21_pbcorr_round_k.fits '+out_dir+'cubes/.'
           print, "Copying noise cubes."
           spawn, 'cp '+in_dir+'*_7m_co21_noise_flat_round_k.fits '+out_dir+'cubes/.'
           print, "Copying signal masks."
           spawn, 'cp '+in_dir+'*_7m_co21_signalmask.fits '+out_dir+'cubes/.'
           print, "Copying hybrid masks."
           spawn, 'cp '+in_dir+'*_7m_co21_hybridmask.fits '+out_dir+'cubes/.'
           for ii = 1, 5 do $
              spawn, 'rm -rf '+out_dir+'*/*_'+str(ii)+'_7m_co21_*.fits'
        endif

        if total(just_array eq '7m+tp') gt 0 then begin
           spawn, 'cp '+in_dir+'*_7m+tp_co21_strict_*.fits '+out_dir+'strict_maps/.'
           spawn, 'cp '+in_dir+'*_7m+tp_co21_broad_*.fits '+out_dir+'broad_maps/.'        
           spawn, 'cp '+in_dir+'*_7m+tp_co21_pbcorr_round_k.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_7m+tp_co21_noise_flat_round_k.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_7m+tp_co21_signalmask.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_7m+tp_co21_hybridmask.fits '+out_dir+'cubes/.'
           for ii = 1, 5 do $
              spawn, 'rm -rf '+out_dir+'*/*_'+str(ii)+'_7m+tp_co21_*.fits'
        endif

        if total(just_array eq '12m+7m') gt 0 then begin
           spawn, 'cp '+in_dir+'*_12m+7m_co21_strict_*.fits '+out_dir+'strict_maps/.'
           spawn, 'cp '+in_dir+'*_12m+7m_co21_broad_*.fits '+out_dir+'broad_maps/.'        
           spawn, 'cp '+in_dir+'*_12m+7m_co21_pbcorr_round_k.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_12m+7m_co21_noise_pbcorr_round_k.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_12m+7m_co21_signalmask.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_12m+7m_co21_hybridmask.fits '+out_dir+'cubes/.'
           for ii = 1, 5 do $
              spawn, 'rm -rf '+out_dir+'*/*_'+str(ii)+'_12m+7m_co21_*.fits'
        endif

        if total(just_array eq '12m+7m+tp') gt 0 then begin
           spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_strict_*.fits '+out_dir+'strict_maps/.'
           spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_broad_*.fits '+out_dir+'broad_maps/.'
           spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_pbcorr_round_k.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_noise_pbcorr_round_k.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_signalmask.fits '+out_dir+'cubes/.'
           spawn, 'cp '+in_dir+'*_12m+7m+tp_co21_hybridmask.fits '+out_dir+'cubes/.'
           for ii = 1, 5 do $
              spawn, 'rm -rf '+out_dir+'*/*_'+str(ii)+'_12m+7m+tp_co21_*.fits'
        endif
        
     endif

     spawn, 'cp README_for_v3.txt '+out_dir+'.'

  endif

end
