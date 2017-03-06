# Script to process cubes produced by clean into FITS files useful for
# future analysis.

tested_versions = ['4.6.0','4.7.0']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
import numpy as np

start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_raw_fits
except NameError:
    do_raw_fits = True

try:
    do_rescale
except NameError:
    do_rescale = True

try:
    do_roundbeam
except NameError:    
    do_roundbeam = True

try:
    do_fits
except NameError:
    do_fits = True

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Inputs
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

abort = False
print "---------- postProcessCubes.py ----------"

# ......................................
# Input/output files
# ......................................

try:
    cube_root
except NameError:
    print "Please define a root output name for the cube via cube_root=XXX."
    abort = True

# ......................................
# Post processing options
# ......................................

try:
    target_beam
except NameError:
    print "I will try to calculate a target beam."    
    do_calc_beam = True

if abort:
    print "(Turning off other parts of the script)."
    do_process = False
    do_shrink = False
    do_fits = False

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXPORT RAW DATA TO FITS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_raw_fits:

    print ""
    print "........................................................"
    print "postProcessCubes: exporting the unprocessed data to FITS."
    print "........................................................"
    print ""

    exportfits(imagename=cube_root+'.image',
               fitsimage=cube_root+'.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)

    exportfits(imagename=cube_root+'_dirty.image',
               fitsimage=cube_root+'.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)

    exportfits(imagename=cube_root+'.model',
               fitsimage=cube_root+'_model.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)

    exportfits(imagename=cube_root+'.residual',
               fitsimage=cube_root+'_residual.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)
    
    exportfits(imagename=cube_root+'.pb',
               fitsimage=cube_root+'_pb.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)
    
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PICK A TARGET ROUND BEAM
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_calc_beam:

    print ""
    print "........................................................"
    print "postProcessCubes: picking a target resolution."
    print "........................................................"
    print ""

    header = imhead(cube_root+'.image')

    execfile('../scripts/extractBeam.py')
    
    print "postProcessCubes: found a beam of "+str(beam_arcsec)+" arcseconds"
    print "postProcessCubes: found "+str(pix_per_beam)+" pixels per beam."
    print "postProcessCubes: found a beam area of "+str(beam_area_pix)+" pixels."

    target_beam_arcsec = \
        sqrt(beam_arcsec**2+(pix_arcsec*2.0)**2)
    target_beam = str(target_beam_arcsec)+'arcsec'
    print "postProcessCubes: will target a round beam of "+str(target_beam_arcsec)+" arcseconds"

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PROCESS THE CUBES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
if do_process:

    print ""
    print "........................................................"
    print "postProcessCubes: processing cubes."
    print "........................................................"
    print ""

    # Smooth to have a round beam        
    imsmooth(imagename=cube_root+'.image',
             outfile=cube_root+'_round.image',
             targetres=True,
             major=target_beam, minor=target_beam, pa='0deg',
             overwrite=True)
        
    # Primary beam correct    
    os.system('rm -rf '+cube_root+'_pbcor.image')
    impbcor(imagename=cube_root+'.image',
            pbimage=cube_root+'.pb',
            outfile=cube_root+'_pbcor.image')

    # Smooth the primary beam corrected cube
    imsmooth(imagename=cube_root+'_pbcor.image',
             outfile=cube_root+'_round_pbcor.image',
             targetres=True,
             major=target_beam, minor=target_beam, pa='0deg',
             overwrite=True)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# SHRINK THE CUBES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_shrink:

    print ""
    print "........................................................"
    print "postProcessCubes: shrinking the cubes to save space."
    print "........................................................"
    print ""

    os.system('rm -rf '+cube_root+'round_pbcor_rebin.image')
    os.system('rm -rf '+cube_root+'_rebin.residual')
    os.system('rm -rf '+cube_root+'_round_rebin.image')
    os.system('rm -rf '+cube_root+'_rebin.pb')

    imrebin(imagename=cube_root+'_round_pbcor.image',
            outfile=cube_root+'_round_pbcor_rebin.image',
            factor=[2,2,1,1])
    
    imrebin(imagename=cube_root+'.residual',
            outfile=cube_root+'_rebin.residual',
            factor=[2,2,1,1])
    
    imrebin(imagename=cube_root+'_round.image',
            outfile=cube_root+'_round_rebin.image',
            factor=[2,2,1,1])
    
    imrebin(imagename=cube_root+'.pb',
            outfile=cube_root+'_rebin.pb',
            factor=[2,2,1,1])

if do_shrink == False:
    
    print ""
    print "........................................................"
    print "postProcessCubes: keeping original pixel scale."
    print "........................................................"
    print ""

    os.system('rm -rf '+cube_root+'round_pbcor_rebin.image')
    os.system('rm -rf '+cube_root+'_rebin.residual')
    os.system('rm -rf '+cube_root+'_round_rebin.image')
    os.system('rm -rf '+cube_root+'_rebin.pb')

    os.system('cp -r '+cube_root+'_round_pbcor.image '+ \
                  cube_root+'_round_pbcor_rebin.image')
    os.system('cp -r '+cube_root+'.residual '+ \
                  cube_root+'_rebin.residual')
    os.system('cp -r '+cube_root+'_round.image '+ \
                  cube_root+'_round_rebin.image')
    os.system('cp -r '+cube_root+'.pb '+ \
                  cube_root+'_rebin.pb')
    

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXPORT TO FITS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_fits:

    print ""
    print "........................................................"
    print "postProcessCubes: exporting to FITS."
    print "........................................................"
    print ""

    # Export to FITS        
    exportfits(imagename=cube_root+'_round_rebin.image',
               fitsimage=cube_root+'_round.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)

    exportfits(imagename=cube_root+'_round_pbcor_rebin.image',
               fitsimage=cube_root+'_round_pbcor.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)
    
    exportfits(imagename=cube_root+'_rebin.residual',
               fitsimage=cube_root+'_residual.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)
        
    exportfits(imagename=cube_root+'_rebin.pb',
               fitsimage=cube_root+'_pb.fits',
               velocity=True, overwrite=True, dropstokes=True, 
               dropdeg=True, bitpix=16)    

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

stop_time = time.time()

elapsed_time = (stop_time - start_time)/60.
print "This run took "+str(elapsed_time)+" minutes"
