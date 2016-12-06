# Script to image one a line visibility into a data cube.

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
total_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_init
except NameError:
    do_init = True

try:
    do_mask
except NameError:    
    do_mask = True

try:
    do_clean
except NameError:
    do_clean = True

try:
    do_postprocess
except NameError:
    do_postprocess = False

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Inputs
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

abort = False
print "---------- imageLine.py ----------"

# ......................................
# Input/output files
# ......................................

try:
    input_vis
except NameError:
    print "Please define an input visibility via input_vis=XXX."
    abort = True

try:
    cube_root
except NameError:
    print "Please define a root output name for the cube via cube_root=XXX."
    abort = True

# ......................................
# Shape and center of the cube
# ......................................

try:
    phase_center
except NameError:
    print "Please specify a phase center for the map (phase_center string)."
    abort = True

try:
    im_size
except NameError:
    print "Please specify an image size (im_size [integer, integer]))."
    abort = True
    
try:
    cell_size
except NameError:
    print "Please specify a cell size (e.g., '0.1arcsec')."
    abort = True

try:
    restfreq
except NameError:
    print "Please specify a rest frequency (restfreq string)."
    abort = True

# ......................................
# Tuning parameters for CLEAN
# ......................................

try:
    deconvolver
except NameError:
    print "Defaulting deconvolver to MULTISCALE."
    deconvolver = 'multiscale'

try:
    threshold
except NameError:
    print "Defaulting to a threshold of 0.0mJy/beam"
    threshold = "0.0mJy/beam"

try:
    scales
except NameError:
    print "I will default to scales with [0,5,15]."
    scales = [0,5,15]

try:
    smallscalebias
except NameError:
    print "I will default to smallscalebias 0.6."
    smallscalebias = 0.6

try:
    briggs_weight
except NameError:
    print "I will default to briggs_weight 0.5."
    briggs_weight = 0.5

try:
    niter
except NameError:
    print "I will default to niter 10000."
    niter = 10000

try:
    cycle_niter
except NameError:
    print "I will default to cycle_niter 200."
    cycle_niter = 10000

try:
    uvtaper
except NameError:
    print "Defaulting to no uvtaper."
    uvtaper = False

# ......................................
# Post processing options
# ......................................

try:
    target_beam
except NameError:
    print "Please specify a target round beam. Otherwise I will calculate it."    

if abort:
    print "(Turning off other parts of the script)."
    do_init = False
    do_mask = False
    do_clean = False
    do_postprocess = False
 
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MAKE THE DIRTY CUBE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_makedirtycube:

    niter = 0
    do_reset = True
    do_callclean = True
    do_savecopy = True    
    execfile('../scripts/callClean.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MAKE A MASK
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_makemask:

    # smooth to large (~10-15") scales and make mask. Expand it even a
    # bit more. Also look for very bright things. 

    # Alternatively, regrid an external mask on to the parent cube.

    print "... ... aligning the mask to the map."        
    os.system('rm -rf '+gal+'_mask_for_native.image')
    imregrid(imagename=gal+'_co21_mask.image'
             , output=gal+'_mask_for_native.image'
             , template=gal+'_co21_cube.residual'
             , interpolation='nearest')

    pass

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CLEAN
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_cleancube:

    do_reset = False
    do_savecopy = False

    cube_hdr = (imhead(cube_root+'.residual'))
    nchan = -1
    if cube_hdr['axisnames'][2] == 'Frequency':
        nchan = cube_hdr['shape'][2]
    if cube_hdr['axisnames'][3] == 'Frequency':
        nchan = cube_hdr['shape'][3]        
    
    niter_per_channel = 100
    loop = 1
    max_loop = 5

    prev_flux = -1
    this_flux = 0
    while converged == False and loop < max_loop:
        do_callclean = True
        niter = niter_per_channel*(2**loop)*nchan
        execfile('../scripts/callClean.py')
        execfile('../scripts/statCleanCube.py')
        
        loop += 1

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# POST PROCESS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_postprocess:
    do_process = True
    do_shrink = True
    do_fits = True
    execfile('.../scripts/postProcessCubes.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

total_stop_time = time.time()

total_elapsed_time = (total_stop_time - total_start_time)/60.
print "This run took "+str(total_elapsed_time)+" minutes"
