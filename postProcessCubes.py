# Script to process cubes produced by clean into FITS files useful for
# future analysis.

tested_versions = ['4.6.0','4.7.0']
this_version = casa['build']['version']
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)
else:
    print "The script has been verified for this version of CASA."

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_process
except NameError:
    do_process = True

try:
    do_shrink
except NameError:    
    do_shrink = True

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
    print "Please specify a target round beam. Otherwise I will try to calculate it."    

if abort:
    print "(Turning off other parts of the script)."
    do_process = False
    do_shrink = False
    do_fits = False

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PROCESS THE CUBES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
if do_process:

    print "........................................................"
    print "postProcessCubes: processing cubes."
    print "........................................................"

    # Add logic to figure out round beam if not user supplied here.

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
             major=target_beam_co21, minor=target_beam_co21, pa='0deg',
             overwrite=True)

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# SHRINK THE CUBES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_shrink:

    print "........................................................"
    print "postProcessCubes: shrinking the cubes to save space."
    print "........................................................"
        
    os.system('rm -rf '+cube_root+'round_pbcor_rebin.image')
    imrebin(imagename=cube_root+'_round_pbcor.image',
            outfile=cube_root+'_round_pbcor_rebin.image',
            factor=[2,2,1,1])
    
    os.system('rm -rf '+cube_root+'_rebin.residual')
    imrebin(imagename=cube_root+'.residual',
            outfile=cube_root+'_rebin.residual',
            factor=[2,2,1,1])
    
    os.system('rm -rf '+cube_root+'_round_rebin.image')
    imrebin(imagename=cube_root+'_round.image',
            outfile=cube_root+'_round_rebin.image',
            factor=[2,2,1,1])
    
    os.system('rm -rf '+cube_root+'_rebin.pb')
    imrebin(imagename=cube_root+'.pb',
            outfile=cube_root+'_rebin.pb',
            factor=[2,2,1,1])

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# EXPORT TO FITS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_shrink:

    print "........................................................"
    print "postProcessCubes: exporting to FITS."
    print "........................................................"

    # Export to FITS        
    exportfits(imagename=cube_root+'_round_rebin.image',
               fitsimage=cube_root+'_round.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)

    exportfits(imagename=cube_root+'_round_pbcor_rebin.image',
               fitsimage=cube_root+'_round_pbcor.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
    
    exportfits(imagename=cube_root+'_rebin.residual',
               fitsimage=cube_root+'_residual.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)
        
    exportfits(imagename=cube_root+'_rebin.pb',
               fitsimage=cube_root+'_pb.fits',
               velocity=True, overwrite=True, dropstokes=True, dropdeg=True, bitpix=16)    

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

stop_time = time.time()

elapsed_time = (stop_time - start_time)/60.
print "This run took "+str(elapsed_time)+" minutes"
