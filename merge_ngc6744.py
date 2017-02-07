# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# Build dirty maps that gives us the target astrometry
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

# This is 100% an ad hoc script. It's just built to do what we need
# for this one case. We could build a general script to make a target
# astrometry

if False:

    #   CO 2-1
    os.system('rm -rf ngc6744_co21_cube*')
    tclean(vis='ngc6744north_956_co21.ms',
           imagename='ngc6744_co21_cube',
           phasecenter='J2000 19h09m46.1s -63d51m27',
           gridder='mosaic',
           deconvolver='hogbom',
           cell='0.15arcsec',
           imsize=[2000,3000],
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

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# Convolve the north and south cubes to have the same resolution
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

if False:

    # Match the southern image to the resolution of the northern one.

    imsmooth(imagename='ngc6744south_co21_round.image',
             outfile='ngc6744south_co21_round_matched.image',
             targetres=True,
             major='1.0arcsec', minor='1.0arcsec', pa='0deg',
             overwrite=True)

    imsmooth(imagename='ngc6744south_co21_round_pbcor.image',
             outfile='ngc6744south_co21_round_pbcor_matched.image',
             targetres=True,
             major='1.0arcsec', minor='1.0arcsec', pa='0deg',
             overwrite=True)

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# Align the other data to the template image
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

if False:
#if True: 
               
    print "Regridding to common grid."

    for ext in ['round', 'round_pbcor']:
        print ext

        os.system('rm -rf ngc6744north_co21_align.image')
        imregrid(imagename='ngc6744north_co21_'+ext+'.image'
                 , output='ngc6744north_co21_'+ext+'_align.image'
                 , template='ngc6744_co21_cube.residual'
                 , overwrite=True)

    for ext in ['round', 'round_pbcor']:
        print ext

        os.system('rm -rf ngc6744south_co21_align.image')
        imregrid(imagename='ngc6744south_co21_'+ext+'_matched.image'
                 , output='ngc6744south_co21_'+ext+'_align.image'
                 , template='ngc6744_co21_cube.residual'
                 , overwrite=True)

    os.system('rm -rf ngc6744north_co21_align.pb')
    imregrid(imagename='ngc6744north_co21.pb'
             , output='ngc6744north_co21_align.pb'
             , template='ngc6744_co21_cube.residual'
             , overwrite=True)

    os.system('rm -rf ngc6744south_co21_align.pb')
    imregrid(imagename='ngc6744south_co21.pb'
             , output='ngc6744south_co21_align.pb'
             , template='ngc6744_co21_cube.residual'
             , overwrite=True)

    os.system('rm -rf ngc6744north_co21_align.residual')
    imregrid(imagename='ngc6744north_co21.residual'
             , output='ngc6744north_co21_align.residual'
             , template='ngc6744_co21_cube.residual'
             , overwrite=True)

    os.system('rm -rf ngc6744south_co21_align.residual')
    imregrid(imagename='ngc6744south_co21.residual'
             , output='ngc6744south_co21_align.residual'
             , template='ngc6744_co21_cube.residual'
             , overwrite=True)

# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&
# EXPORT TO FITS
# %&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&

if False:

    print "Exporting to FITS."
        
    exportfits(imagename='ngc6744north_co21_round_align.image',
               fitsimage='ngc6744north_co21_round_align.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    exportfits(imagename='ngc6744south_co21_round_align.image',
               fitsimage='ngc6744south_co21_round_align.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    exportfits(imagename='ngc6744north_co21_round_pbcor_align.image',
               fitsimage='ngc6744north_co21_round_pbcor_align.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    exportfits(imagename='ngc6744south_co21_round_pbcor_align.image',
               fitsimage='ngc6744south_co21_round_pbcor_align.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    exportfits(imagename='ngc6744north_co21_align.residual',
               fitsimage='ngc6744north_co21_residual_align.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    exportfits(imagename='ngc6744south_co21_align.residual',
               fitsimage='ngc6744south_co21_residual_align.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    exportfits(imagename='ngc6744north_co21_align.pb',
               fitsimage='ngc6744north_co21_pb_align.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    exportfits(imagename='ngc6744south_co21_align.pb',
               fitsimage='ngc6744south_co21_pb_align.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

