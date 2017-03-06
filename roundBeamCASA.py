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
# OPTIONS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    cube_root
except NameError:
    print "Please define a root output name for the cube via cube_root=XXX."
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PROCESS THE CUBES
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
    
print ""
print "........................................................"
print "postProcessCubes: processing cubes."
print "........................................................"
print ""

for ext in ['_pbcorr','_rescale_pbcorr','_clip','_rescale_clip']:
    
    print "Calculating smooth beam for "+cube_root+ext+'.fits'
        
    imsmooth(imagename=cube_root+ext+'.fits',
             outfile=cube_root+ext+'_round.fits',
             targetres=True,
             major=target_beam, minor=target_beam, pa='0deg',
             overwrite=True)

stop_time = time.time()

elapsed_time = (stop_time - start_time)/60.
print "This run took "+str(elapsed_time)+" minutes"
