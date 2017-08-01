# Note that "IMPROVEMENT TBD" items are scattered throughout the
# program. Search on this to find areas where some clean up or a next
# logical step could be worked out.

tested_versions = ['4.6.0','4.7.0','4.7.1','4.7.2']
this_version = (casa['build']['version']).split('-')[0]
if this_version not in tested_versions:
    print "The script hasn't been verified for this version of CASA."
    print "This version of CASA is "+this_version
    print "Tested versions are "+str(tested_versions)

execfile('../scripts/auExtensions.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

import time
total_start_time = time.time()

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CALCULATE CELL SIZE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Factor by which to oversample the beam
try:
    oversamp
except NameError:
    oversamp = 5

# Cell size implied by baseline distribution
au_cellsize, au_imsize, au_centralField = \
    au.pickCellSize(input_vis, imsize=True, npix=oversamp)
xextent = au_cellsize*au_imsize[0]*1.2
yextent = au_cellsize*au_imsize[1]*1.2

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

# If we taper, then the cell size may be forced
print "I calculated cell size: ", cell_size
try:
    force_cell_size
except NameError:
    force_cell_size = None

if force_cell_size != None:
    cell_size = force_cell_size
    print "I will force cell size: ", cell_size

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# CALCULATE IMAGE SIZE
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

# Now make the image size a nice round number
need_cells_x = xextent / cell_size
need_cells_y = yextent / cell_size

# Valid image sizes are even and multiples of 3, 5, 7
valid_sizes = []
for ii in range(10):
    for kk in range(3):
        for jj in range(3):
            valid_sizes.append(2**(ii+1)*5**(jj)*3**(kk))
valid_sizes.sort()
valid_sizes = np.array(valid_sizes)

cells_x = np.min(valid_sizes[valid_sizes > need_cells_x])
cells_y = np.min(valid_sizes[valid_sizes > need_cells_y])

image_size = [int(cells_x), int(cells_y)]
cell_size_string = str(cell_size)+'arcsec'

print "I will use the following:"
print "... cell size = "+cell_size_string
print "... image size = ", image_size
    
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# ENSURE THAT WE HAVE A LIST OF TAPERS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

uvtaper_list = \
    ['','2.5','5','10','5','2.5','']

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# LOOP OVER TAPERS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

first = True

for uvind in range(len(uvtaper_list)):

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Making a dirty cube."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    uv_taper = uvtaper_list[uvind]
    print "UV TAPER :::", uv_taper
    if uv_taper == '':
        uv_taper_string = uv_taper
    else:
        uv_taper_string = [uv_taper+'arcsec',uv_taper+'arcsec','0deg']    

    if first:
        do_reset = True
    else:
        do_reset = False

    niter = 0
    do_callclean = True
    do_savecopy = False

    bkup_ext = "dirty"
    logfile = cube_root+"_dirty.log"

    usemask = 'pb'
    mask = ''
    pbmask=pb_limit
    calcres = True
    calcpsf = True    
    scales = [0]
    execfile('../scripts/callClean.py')
    reset = False

    if first:

        print ""
        print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
        print "Aligning the external clean mask to the dirty cube."
        print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
        print ""

        mask_in = clean_mask_file
        mask_root = cube_root
        execfile('../scripts/alignMask.py')

        ia.open(cube_root+".pb")
        pb = ia.getchunk()
        ia.close()
    
        ia.open(cube_root+".mask")
        mask = ia.getchunk()
        mask *= (pb > pb_limit)
        ia.putchunk(mask)
        ia.close()

        mask = ''
        pb = 0.0
        
    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Staging the clean call."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""
    
    execfile('../scripts/statCleanCube.py')    
    threshold_value = snr_thresh* \
        imstat_residual['medabsdevmed'][0]/0.6745
    threshold = str(threshold_value)+'Jy/beam'

    # Note the number of channels
    vm = au.ValueMapping(input_vis)
    nchan = vm.spwInfo[0]['numChannels']

    do_callclean = True
    do_savecopy = False
    this_niter = nchan*10000
    cycle_niter = 1000
    logfile = "dummy.log"
    
    # Clean and associated tuning parameters
    calcres = True
    calcpsf = True
    minpsffraction = 0.5
    usemask = 'user'
    mask = ''
    deconvolver = 'multiscale'
    restoringbeam = 'common'
    niter = this_niter
    if uv_taper == '':
        scales = [0]
    else:
        scales = [int(float(uv_taper)/cell_size)]

    do_savecopy = False
    bkup_ext = "dummy"

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Calling clean for scale "+uv_taper+" or "+str(scales[0])+" pixels"
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    execfile('../scripts/callClean.py')

    print ""
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print "Cleaning up loop then proceeding."
    print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
    print ""

    # Save a copy
    do_savecopy = False
    do_callclean = False
    bkup_ext = "uvtaper"+uv_taper
    execfile('../scripts/callClean.py')

    first = False

print ""
print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
print "Final imaging without a uv taper."
print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
print ""

uv_taper_string = ''
''
niter = 0
do_callclean = True
do_savecopy = False   

bkup_ext = "final"
logfile = cube_root+"_final.log"

usemask = 'pb'
mask = ''
pbmask=pb_limit
calcres = True
calcpsf = True    
scales = [0]
execfile('../scripts/callClean.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# POST PROCESS
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

print ""
print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
print "Exporting the data."
print "+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+"
print ""

execfile('../scripts/exportToFITS.py')

# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
# TIMER
# &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

total_stop_time = time.time()

total_elapsed_time = (total_stop_time - total_start_time)/60.
print "This run took "+str(total_elapsed_time)+" minutes"
