# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# Build dirty maps that gives us the target astrometry
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

if False:

    #   CO 2-1

    os.system('rm -rf ngc6744_co21_cube*')
    tclean(vis='ngc6744north_co21.ms',
           imagename='ngc6744_co21_cube',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.3arcsec',
           imsize=[1000,1500],
           weighting='briggs',
           robust=0.5,
           specmode='cube',
           restfreq='230.53800GHz' ,
           outframe='lsrk',
           veltype='radio',
           niter=0,
           threshold='0.015Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   CO 2-1 TAPER

    os.system('rm -rf ngc6744_co21_taper*')
    tclean(vis='ngc6744north_co21.ms',
           imagename='ngc6744_co21_taper',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.5arcsec',
           imsize=[600,1200],
           weighting='briggs',
           robust=0.5,
           specmode='cube',
           restfreq='230.53800GHz' ,
           outframe='lsrk',
           veltype='radio',
           niter=0,
           threshold='0.015Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   C18O

    os.system('rm -rf ngc6744_c18o21_cube*')
    tclean(vis='ngc6744north_c18o21.ms',
           imagename='ngc6744_c18o21_cube',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.3arcsec',
           imsize=[1000,1500],
           weighting='briggs',
           robust=0.5,
           specmode='cube',
           restfreq='219.56035GHz',
           outframe='lsrk',
           veltype='radio',
           niter=0,
           threshold='0.015Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   C18O TAPER

    os.system('rm -rf ngc6744_c18o21_taper*')
    tclean(vis='ngc6744north_c18o21.ms',
           imagename='ngc6744_c18o21_taper',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.5arcsec',
           imsize=[600,1200],
           weighting='briggs',
           robust=0.5,
           specmode='cube',
           restfreq='219.56035GHz',
           outframe='lsrk',
           veltype='radio',
           niter=0,
           threshold='0.015Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   MM CONTINUUM

    os.system('rm -rf ngc6744_cont_map*')
    tclean(vis='ngc6744north_cont.ms',
           imagename='ngc6744_cont_map',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.3arcsec',
           imsize=[1000,1500],
           weighting='briggs',
           robust=0.5,
           specmode='mfs',
           niter=0,
           threshold='0.015Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   MM CONTINUUM TAPER

    os.system('rm -rf ngc6744_cont_taper*')
    tclean(vis='ngc6744north_cont.ms',
           imagename='ngc6744_cont_taper',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.5arcsec',
           imsize=[600,1200],
           weighting='briggs',
           robust=0.5,
           specmode='mfs',
           niter=0,
           threshold='0.0Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   CO 2-1 CHANNEL 0

    os.system('rm -rf ngc6744_chan0_co21_map*')
    tclean(vis='ngc6744north_chan0_co21.ms',
           imagename='ngc6744_chan0_co21_map',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.3arcsec',
           imsize=[1000,1500],
           weighting='briggs',
           robust=0.5,
           specmode='mfs',
           niter=0,
           threshold='0.015Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   CO 2-1 CHANNEL 0 TAPER

    os.system('rm -rf ngc6744_chan0_co21_taper*')
    tclean(vis='ngc6744north_chan0_co21.ms',
           imagename='ngc6744_chan0_co21_taper',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.5arcsec',
           imsize=[600,1200],
           weighting='briggs',
           robust=0.5,
           specmode='mfs',
           niter=0,
           threshold='0.0Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   C18O 2-1 CHANNEL 0 

    os.system('rm -rf ngc6744_chan0_c18o21_map*')
    tclean(vis='ngc6744north_chan0_c18o21.ms',
           imagename='ngc6744_chan0_c18o21_map',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.3arcsec',
           imsize=[1000,1500],
           weighting='briggs',
           robust=0.5,
           specmode='mfs',
           niter=0,
           threshold='0.015Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

    #   C18O 2-1 CHANNEL 0 TAPER

    os.system('rm -rf ngc6744_chan0_c18o21_taper*')
    tclean(vis='ngc6744north_cont.ms',
           imagename='ngc6744_chan0_c18o21_taper',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='clark',
           cell='0.5arcsec',
           imsize=[600,1200],
           weighting='briggs',
           robust=0.5,
           specmode='mfs',
           niter=0,
           threshold='0.0Jy/beam',
           interactive=False,
           usemask='pb', 
           pbmask=0.2,
           )

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# Align the other data to the template image
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

    

for field in ['north', 'south']:

#   --------------------------------------------------------------
#   CUBES
#   --------------------------------------------------------------

    print "-----------------------"
    print "LOOPING OVER DATA CUBES"
    print "-----------------------"

    for line in ['co21', 'c18o21']:
        
        #   --------------------------------------------------------------
        #   REGRID
        #   --------------------------------------------------------------

        print "Regridding native resolution images."

        os.system('rm -rf ngc6744'+field+'_'+line+'_align.image')
        imregrid(imagename='ngc6744'+field+'_'+line+'_cube_round.image'
                 , output='ngc6744'+field+'_'+line+'_cube_align.image'
                 , template='ngc6744_'+line+'_cube.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+line+'_align.pb')
        imregrid(imagename='ngc6744'+field+'_'+line+'_cube.pb'
                 , output='ngc6744'+field+'_'+line+'_cube_align.pb'
                 , template='ngc6744_'+line+'_cube.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+line+'_pbcor_align.image')
        imregrid(imagename='ngc6744'+field+'_'+line+'_cube_round_pbcor.image'
                 , output='ngc6744'+field+'_'+line+'_cube_pbcor_align.image'
                 , template='ngc6744_'+line+'_cube.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+line+'_align.residual')
        imregrid(imagename='ngc6744'+field+'_'+line+'_cube.residual'
                 , output='ngc6744'+field+'_'+line+'_cube_align.residual'
                 , template='ngc6744_'+line+'_cube.residual'
                 , overwrite=True)

        print "Regridding tapered resolution images."

        os.system('rm -rf ngc6744'+field+'_'+line+'_taper_align.image')
        imregrid(imagename='ngc6744'+field+'_'+line+'_taper_round.image'
                 , output='ngc6744'+field+'_'+line+'_taper_align.image'
                 , template='ngc6744_'+line+'_taper.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+line+'_taper_align.pb')
        imregrid(imagename='ngc6744'+field+'_'+line+'_taper.pb'
                 , output='ngc6744'+field+'_'+line+'_taper_align.pb'
                 , template='ngc6744_'+line+'_taper.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+line+'_taper_pbcor_align.image')
        imregrid(imagename='ngc6744'+field+'_'+line+'_taper_round_pbcor.image'
                 , output='ngc6744'+field+'_'+line+'_taper_pbcor_align.image'
                 , template='ngc6744_'+line+'_taper.residual'
                 , overwrite=True)

        os.system('rm -rf ngc6744'+field+'_'+line+'_taper_align.residual')
        imregrid(imagename='ngc6744'+field+'_'+line+'_taper.residual'
                 , output='ngc6744'+field+'_'+line+'_taper_align.residual'
                 , template='ngc6744_'+line+'_taper.residual'
                 , overwrite=True)

        #   --------------------------------------------------------------
        #   EXPORT TO FITS
        #   --------------------------------------------------------------

        print "Exporting native resolution images."
        
        exportfits(imagename='ngc6744'+field+'_'+line+'_cube_align.image',
                   fitsimage='ngc6744'+field+'_'+line+'_cube_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+line+'_cube_align.pb',
                   fitsimage='ngc6744'+field+'_'+line+'_cube_pb_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+line+'_cube_pbcor_align.image',
                   fitsimage='ngc6744'+field+'_'+line+'_cube_pbcor_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+line+'_cube_align.residual',
                   fitsimage='ngc6744'+field+'_'+line+'_cube_residual_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        print "Exporting tapered images."

        exportfits(imagename='ngc6744'+field+'_'+line+'_taper_align.image'
                   , fitsimage='ngc6744'+field+'_'+line+'_taper_align.fits'
                   , velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+line+'_taper_align.pb',
                   fitsimage='ngc6744'+field+'_'+line+'_taper_pb_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+line+'_taper_pbcor_align.image',
                   fitsimage='ngc6744'+field+'_'+line+'_taper_pbcor_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+line+'_taper_align.residual',
                   fitsimage='ngc6744'+field+'_'+line+'_taper_residual_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

#   --------------------------------------------------------------
#   MAPS
#   --------------------------------------------------------------

    print "-----------------------"
    print "LOOPING OVER IMAGES"
    print "-----------------------"

    for image in ['cont', 'chan0_co21', 'chan0_c18o21']:

#   --------------------------------------------------------------
#   REGRID
#   --------------------------------------------------------------

        print "Regridding native resolution images."

        os.system('rm -rf ngc6744'+field+'_'+image+'_align.image')
        imregrid(imagename='ngc6744'+field+'_'+image+'_round.image'
                 , output='ngc6744'+field+'_'+image+'_align.image'
                 , template='ngc6744_'+image+'_map.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+image+'_align.pb')
        imregrid(imagename='ngc6744'+field+'_'+image+'_map.pb'
                 , output='ngc6744'+field+'_'+image+'_align.pb'
                 , template='ngc6744_'+image+'_map.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+image+'_pbcor_align.image')
        imregrid(imagename='ngc6744'+field+'_'+image+'_round_pbcor.image'
                 , output='ngc6744'+field+'_'+image+'_pbcor_align.image'
                 , template='ngc6744_'+image+'_map.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+image+'_align.residual')
        imregrid(imagename='ngc6744'+field+'_'+image+'_map.residual'
                 , output='ngc6744'+field+'_'+image+'_align.residual'
                 , template='ngc6744_'+image+'_map.residual'
                 , overwrite=True)

        print "Regridding tapered resolution images."

        os.system('rm -rf ngc6744'+field+'_'+image+'_taper_align.image')
        imregrid(imagename='ngc6744'+field+'_'+image+'_taper_round.image'
                 , output='ngc6744'+field+'_'+image+'_taper_align.image'
                 , template='ngc6744_'+image+'_taper.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+image+'_taper_align.pb')
        imregrid(imagename='ngc6744'+field+'_'+image+'_taper.pb'
                 , output='ngc6744'+field+'_'+image+'_taper_align.pb'
                 , template='ngc6744_'+image+'_taper.residual'
                 , overwrite=True)
        
        os.system('rm -rf ngc6744'+field+'_'+image+'_taper_pbcor_align.image')
        imregrid(imagename='ngc6744'+field+'_'+image+'_taper_round_pbcor.image'
                 , output='ngc6744'+field+'_'+image+'_taper_pbcor_align.image'
                 , template='ngc6744_'+image+'_taper.residual'
                 , overwrite=True)

        os.system('rm -rf ngc6744'+field+'_'+image+'_taper_align.residual')
        imregrid(imagename='ngc6744'+field+'_'+image+'_taper.residual'
                 , output='ngc6744'+field+'_'+image+'_taper_align.residual'
                 , template='ngc6744_'+image+'_taper.residual'
                 , overwrite=True)

#   --------------------------------------------------------------
#   EXPORT TO FITS
#   --------------------------------------------------------------

        print "Exporting native resolution images."
        
        exportfits(imagename='ngc6744'+field+'_'+image+'_align.image',
                   fitsimage='ngc6744'+field+'_'+image+'_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+image+'_align.pb',
                   fitsimage='ngc6744'+field+'_'+image+'_pb_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+image+'_pbcor_align.image',
                   fitsimage='ngc6744'+field+'_'+image+'_pbcor_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+image+'_align.residual',
                   fitsimage='ngc6744'+field+'_'+image+'_residual_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        print "Exporting tapered images."

        exportfits(imagename='ngc6744'+field+'_'+image+'_taper_align.image'
                   , fitsimage='ngc6744'+field+'_'+image+'_taper_align.fits'
                   , velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+image+'_taper_align.pb',
                   fitsimage='ngc6744'+field+'_'+image+'_taper_pb_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+image+'_taper_pbcor_align.image',
                   fitsimage='ngc6744'+field+'_'+image+'_taper_pbcor_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
        exportfits(imagename='ngc6744'+field+'_'+image+'_taper_align.residual',
                   fitsimage='ngc6744'+field+'_'+image+'_taper_residual_align.fits',
                   velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
