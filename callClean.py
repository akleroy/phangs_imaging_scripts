# Script to call clean ONCE as part of producing a cube from
# visibility data.

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
start_time_clean = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_callclean
except NameError:
    do_savecopy = True

try:
    do_savecopy
except NameError:
    do_savecopy = False

try:
    do_reset
except NameError:
    do_reset = False

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# Inputs
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

abort = False
print "---------- callClean.py ----------"

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
# If we abort, turn off the script
# ......................................

if abort:
    print "(Turning off other parts of the script)."
    do_clean = False
    do_savecopy = False
 

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# PROCEED WITH THE CLEAN
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_callclean:

    print "........................................................"
    print "callClean: calling clean."
    print "........................................................"

    if do_reset:
        print "...wiping previous versions of the cube."
        os.system('rm -rf '+cube_root+'.*')
        
    try:
        logfile
    except NameError:
        print "... leaving the logfile unset."
    else:
        print "... setting the logfile to "+logfile
        casalog.setlogfile(logfile)
            
    try:
        mask_file
    except NameError:
        usemask = 'pb'
    else:
        if mask_file != '':
            usemask = 'user'
            mask = mask_file

    tclean(vis=input_vis,
           imagename=cube_root,
           # Spatial axes
           phasecenter=phase_center,
           cell=cell_size,
           imsize=im_size,
           gridder='mosaic',
           # Spectral axis
           specmode='cube',
           restfreq=restfreq,
           outframe='lsrk',
           veltype='radio',
           # Deconvolver
           deconvolver=deconvolver,
           scales=scales,
           smallscalebias=smallscalebias,
           # U-V plane gridding
           weighting='briggs',
           robust=briggs_weight,
           uvtaper=uvtaper,
           # Stopping criterion
           niter=niter,
           threshold=threshold,
           cycleniter=cycle_niter,
           cyclefactor=2.0,
           minpsffraction=0.1,
           # Mask
           usemask=usemask,
           pbmask=0.2,
           # UI
           interactive=False,
           )
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# POST PROCESS AND WRITE THE RESULTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_savecopy:

    print "........................................................"
    print "callClean: Saving a copy of the results."
    print "........................................................"

    try:
        bkup_root
    except NameError:
        print "Please define a root output name for the saved cube via bkup_root=XXX."
    
    os.system('rm -rf '+bkup_root+'.image')
    os.system('rm -rf '+bkup_root+'.model')
    os.system('rm -rf '+bkup_root+'.pb')
    os.system('rm -rf '+bkup_root+'.psf')
    os.system('rm -rf '+bkup_root+'.residual')

    os.system('cp -r '+cube_root+'.image '+bkup_root+'.image')
    os.system('cp -r '+cube_root+'.model '+bkup_root+'.model')
    os.system('cp -r '+cube_root+'.pb '+bkup_root+'.pb')
    os.system('cp -r '+cube_root+'.psf '+bkup_root+'.psf')
    os.system('cp -r '+cube_root+'.residual '+bkup_root+'.residual')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

stop_time_clean = time.time()

elapsed_time_clean = (stop_time_clean - start_time_clean)/60.
print "This CLEAN run took "+str(elapsed_time_clean)+" minutes"
