# Script to collapse a measurement set into a one-plane image.

tested_versions = ['4.6.0','4.7.0']
this_version = casa['build']['version']
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)
else:
    print "The script has been verified for this version of CASA."

# Hate this ... but directories
execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
total_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CONTROL FLOW
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

try:
    do_pickcellsize
except NameError:
    do_pickcellsize = True

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
# Other option
# ......................................

if abort:
    print "(Turning off other parts of the script)."
    do_init = False
    do_mask = False
    do_clean = False
    do_postprocess = False
 
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# MAKE THE DIRTY CUBE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

if do_pickcellsize:
    au_cellsize, au_imsize, au_centralField = \
        au.pickCellSize(input_vis, imsize=True)
    xextent = au_cellsize*au_imsize[0]
    yextent = au_cellsize*au_imsize[1]
    
    # Make the cell size a nice round number
    if au_cellsize < 0.1:
        cell_size = au_cellsize
    if au_cellsize >= 0.1 and au_cellsize < 0.5:
        cell_size = floor(au_cellsize/0.05)*0.05
    if au_cellsize >= 0.5 and au_cellsize < 1.0:
        cell_size = floor(au_cellsize/0.1)*0.1
    if au_cellsize >= 1.0 and au_cellsize < 2.0:
        cell_size = floor(au_cellsize/0.25)*0.25
    if au_cellsize >= 2.0 and au_cellsize < 5.0:
        cell_size = floor(au_cellsize/0.5)*0.5
    if au_cellsize >= 5.0:
        cell_size = floor(au_cellsize/1.0)*0.5

    # Now make the image size a nice round number
    need_cells_x = xextent / cell_size
    need_cells_y = yextent / cell_size

    cells_x = (floor(need_cells_x/10)+1*((need_cells_x % 10) != 0))*10
    cells_y = (floor(need_cells_y/10)+1*((need_cells_y % 10) != 0))*10

    image_size = [cells_x, cells_y]
    cell_size_string = str(cell_size)+'arcsec'

    print "I propose the following:"
    print "... cell size = "+cell_size_string
    print "... image size = ", image_size
    print "... central field = ", au_centralField
    
if do_makedirtycube:

    niter = 0
    do_reset = True
    do_callclean = True
    do_savecopy = True    

    # Single plane
    specmode = 'mfs'        

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
    
    nchan = 1
    base_niter_per_channel = 100
    base_cycle_niter = 100
    loop = 1
    max_loop = 5

    prev_flux = -1
    this_flux = 0
    while converged == False and loop <= max_loop:
        do_callclean = True
        niter = base_niter_per_channel*(2**loop)*nchan
        cycle_niter = base_cycle_niter*(2**loop)
        logfile = "clean_loop_"+str(loop)+".log"
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
