# Script to call clean ONCE as part of producing a cube from
# visibility data.

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
# Workflow
# ......................................

try:
    calcres
except NameError:
    calcres = True

# ......................................
# Shape and center of the cube
# ......................................

try:
    phase_center
except NameError:
    print "Please specify a phase center for the map (phase_center string)."
    abort = True

try:
    image_size
except NameError:
    print "Please specify an image size (image_size [integer, integer]))."
    abort = True
    
try:
    cell_size_string
except NameError:
    print "Please specify a cell size (e.g., '0.15arcsec')."
    abort = True

try:
    restfreq_ghz
except NameError:
    print "Please specify a rest frequency (restfreq_ghz)."
    print "Will default leaving this unset."
    restfreq_ghz = -1.0

try:
    specmode
except NameError:
    print "Spectral mode defaulting to a single continuum image."
    print "Should be 'cube' for a cube. I use 'mfs'."
    specmode = "mfs"

# ......................................
# Tuning parameters for CLEAN
# ......................................

try:
    deconvolver
except NameError:
    #print "Defaulting deconvolver to MULTISCALE."
    #deconvolver = 'multiscale'
    deconvolver = "hogbom"

try:
    threshold
except NameError:
    print "Defaulting to a threshold of 0.0mJy/beam"
    threshold = "0.0mJy/beam"

try:
    scales
except NameError:
    scales=[0]
    print "I will default to scales (in pixels) of ", scales
    #scales = [0,5,10,20,40]
    # Andreas - factors of a few 12m beam and a scale ~1 and ~2 times 7m beam
    # AKL - interaction with "taper" is weird. Need to think about this.

try:
    smallscalebias
except NameError:
    print "I will default to smallscalebias 0.9."
    smallscalebias = 0.9
    # Andreas ~0.8 or 0.9

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
    cycle_niter = 200

try:
    uv_taper_string
except NameError:
    print "Defaulting to no uvtaper."
    uv_taper_string = ''

try:
    minpsffraction
except NameError:
    print "Defaulting to a VERY HIGH minpsffraction."
    minpsffraction = 0.5

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
        oldlogfile = casalog.logfile()
        os.system('rm -rf '+logfile)
        casalog.setlogfile(logfile)
            
    #try:
    #    mask_file
    #except NameError:
    #    usemask = 'pb'
    #    mask = ''
    #else:
    #    if mask_file != '':
    #        usemask = 'user'
    #        mask = mask_file
    #    else:
    #        usemask = 'pb'
    #        mask = ''

    if restfreq_ghz < 0:
        restfreq_str = ''
    else:
        restfreq_str = str(restfreq_ghz)+'GHz'

    tclean(vis=input_vis,
           imagename=cube_root,
           # Spatial axes
           phasecenter=phase_center,
           cell=cell_size_string,
           imsize=image_size,
           gridder='mosaic',
           # Spectral axis
           specmode=specmode,
           restfreq=restfreq_str,
           outframe='lsrk',
           veltype='radio',
           # Workflow
           calcres=calcres,
           # Deconvolver
           deconvolver=deconvolver,
           scales=scales,
           smallscalebias=smallscalebias,
           # U-V plane gridding
           weighting='briggs',
           robust=briggs_weight,
           uvtaper=uv_taper_string,
           # Stopping criterion
           niter=niter,
           threshold=threshold,
           cycleniter=cycle_niter,
           cyclefactor=3.0,
           minpsffraction=minpsffraction,
           # Mask
           usemask=usemask,
           mask=mask,
           pbmask=0.2,
           # UI
           interactive=False,
           )

    try:
        logfile
    except NameError:
        pass
    else:
        print "... unsetting the logfile."
        casalog.setlogfile(oldlogfile)
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# POST PROCESS AND WRITE THE RESULTS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_savecopy:

    try:
        bkup_ext
    except NameError:
        print "do_savecopy is on but bkup_ext is not specified. No copy saved."
        bkup_ext = None

    if bkup_ext != None:

        print "........................................................"
        print "callClean: Saving a copy of the results."
        print "........................................................"
        
        os.system('rm -rf '+cube_root+'_'+bkup_ext+'.image')
        os.system('rm -rf '+cube_root+'_'+bkup_ext+'.model')
        os.system('rm -rf '+cube_root+'_'+bkup_ext+'.mask')
        os.system('rm -rf '+cube_root+'_'+bkup_ext+'.pb')
        os.system('rm -rf '+cube_root+'_'+bkup_ext+'.psf')
        os.system('rm -rf '+cube_root+'_'+bkup_ext+'.residual')
        
        os.system('cp -r '+cube_root+'.image '+cube_root+'_'+bkup_ext+'.image')
        os.system('cp -r '+cube_root+'.model '+cube_root+'_'+bkup_ext+'.model')
        os.system('cp -r '+cube_root+'.mask '+cube_root+'_'+bkup_ext+'.mask')
        os.system('cp -r '+cube_root+'.pb '+cube_root+'_'+bkup_ext+'.pb')
        os.system('cp -r '+cube_root+'.psf '+cube_root+'_'+bkup_ext+'.psf')
        os.system('cp -r '+cube_root+'.residual '+cube_root+'_'+bkup_ext+'.residual')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

stop_time_clean = time.time()

elapsed_time_clean = (stop_time_clean - start_time_clean)/60.
print "This CLEAN run took "+str(elapsed_time_clean)+" minutes"
